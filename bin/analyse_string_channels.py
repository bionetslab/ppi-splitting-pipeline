#!/usr/bin/env python3
"""
Standalone analysis: how well do STRING evidence channels predict from pair embeddings?

For each evidence channel, a Ridge regressor is fitted (with cross-validation) to
predict the channel score from the concatenated embeddings of the two proteins.
Spearman ρ between predicted and actual scores is reported.

Only POSITIVE pairs are used, because negatives have no STRING score — including
them at 0.0 would trivially confound the analysis with the label.

Evidence channels analysed (combined_score is excluded because it is derived from
the others and would just reflect their aggregate):
  experiments, experiments_transferred,
  database, database_transferred,
  textmining, textmining_transferred

Usage
-----
    python analyse_string_channels.py \\
        --ppis       ppis.csv \\
        --embeddings embeddings.npz \\
        [--cv 5] [--seed 42] [--out string_channel_analysis.tsv]
"""

import argparse
import csv
import os
import sys

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_predict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import load_embeddings

STRING_CHANNELS = [
    "experiments",
    "experiments_transferred",
    "database",
    "database_transferred",
    "textmining",
    "textmining_transferred",
]


def load_positives(ppis_path, channels):
    """Return (rows, available_channels) for positive pairs only.

    Each row is (protein1, protein2, {channel: score/1000}).
    Channels absent from the file are silently dropped from available_channels.
    """
    rows = []
    available = []
    with open(ppis_path) as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []
        available = [c for c in channels if c in fieldnames]
        if not available:
            sys.exit(f"None of the STRING channels found in {ppis_path}. "
                     "Make sure the file has STRING score columns.")
        for row in reader:
            try:
                label = int(row.get("label", row.get("interaction", "1")))
            except ValueError:
                label = 1
            if label != 1:
                continue
            p1 = row["protein1"].strip()
            p2 = row["protein2"].strip()
            scores = {}
            for c in available:
                try:
                    scores[c] = float(row[c]) / 1000.0
                except (ValueError, KeyError):
                    scores[c] = 0.0
            rows.append((p1, p2, scores))
    return rows, available


def build_matrices(rows, embeddings, available):
    """Build X (pair embeddings) and per-channel A arrays.

    Pairs where either protein has no embedding are skipped.
    Returns (X, {channel: A_array}).
    """
    X_list = []
    channel_lists = {c: [] for c in available}

    for p1, p2 in ((r[0], r[1]) for r in rows):
        if p1 not in embeddings or p2 not in embeddings:
            continue
        a, b = (p1, p2) if p1 <= p2 else (p2, p1)
        X_list.append(np.concatenate([embeddings[a], embeddings[b]]))

    # Re-iterate to keep alignment with X_list
    idx = 0
    for p1, p2, scores in rows:
        if p1 not in embeddings or p2 not in embeddings:
            continue
        for c in available:
            channel_lists[c].append(scores[c])
        idx += 1

    X = np.array(X_list, dtype=np.float32)
    A_by_channel = {c: np.array(v, dtype=np.float32) for c, v in channel_lists.items()}
    return X, A_by_channel


def main():
    ap = argparse.ArgumentParser(
        description="Predict STRING evidence channel scores from pair embeddings (positives only)."
    )
    ap.add_argument("--ppis",       required=True,
                    help="PPI CSV with 'protein1', 'protein2', 'label' and STRING score columns")
    ap.add_argument("--embeddings", required=True,
                    help="Pre-computed embeddings .npz (output of embed_sequences.py)")
    ap.add_argument("--cv",         type=int, default=5,
                    help="Cross-validation folds for Ridge regression (default: 5)")
    ap.add_argument("--seed",       type=int, default=42)
    ap.add_argument("--out",        default="string_channel_analysis.tsv",
                    help="Output TSV path (default: string_channel_analysis.tsv)")
    args = ap.parse_args()

    print("Loading embeddings ...", file=sys.stderr)
    embeddings = load_embeddings(args.embeddings)

    print(f"Loading positive pairs from {args.ppis} ...", file=sys.stderr)
    rows, available = load_positives(args.ppis, STRING_CHANNELS)
    print(f"  {len(rows)} positive pairs loaded", file=sys.stderr)
    print(f"  Available channels: {available}", file=sys.stderr)

    X, A_by_channel = build_matrices(rows, embeddings, available)
    print(f"  {X.shape[0]} pairs retained after embedding lookup", file=sys.stderr)

    if X.shape[0] < args.cv:
        sys.exit(f"Too few pairs ({X.shape[0]}) for {args.cv}-fold CV. "
                 "Reduce --cv or provide more data.")

    header = f"{'Channel':<30} {'Spearman ρ':>12} {'p-value':>14} {'Score mean':>12} {'Score std':>10}"
    print(f"\n{header}")
    print("-" * len(header))

    results = []
    for channel in available:
        A = A_by_channel[channel]
        if A.std() == 0:
            print(f"{channel:<30} {'(constant, skipped)':>12}")
            continue

        reg = Ridge(alpha=1.0, solver="sag", random_state=args.seed, max_iter=100)
        A_pred = cross_val_predict(reg, X, A, cv=args.cv)
        rho, pval = spearmanr(A, A_pred)

        print(f"{channel:<30} {rho:>12.4f} {pval:>14.4e} {A.mean():>12.4f} {A.std():>10.4f}")
        results.append((channel, rho, pval, float(A.mean()), float(A.std())))

    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["channel", "spearman_rho", "p_value", "score_mean", "score_std"])
        writer.writerows(results)
    print(f"\nResults written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
