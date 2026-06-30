#!/usr/bin/env python3
"""
Standalone analysis: how well do STRING evidence channels predict from pair embeddings?

A Ridge regressor is fitted on positive pairs from the training split and evaluated
on positive pairs from the test split. Spearman ρ between predicted and actual scores
is reported per channel.

Only POSITIVE pairs are used, because negatives have no STRING score — including
them at 0.0 would trivially confound the analysis with the label.

Evidence channels analysed (combined_score excluded — it is derived from the others):
  experiments, experiments_transferred,
  database, database_transferred,
  textmining, textmining_transferred

Usage
-----
    python analyse_string_channels.py \\
        --train      train.csv \\
        --test       test.csv \\
        --embeddings embeddings.npz \\
        [--seed 42] [--out string_channel_analysis.tsv]
"""

import argparse
import csv
import os
import sys

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge

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


def load_positives(path, channels):
    """Return (rows, available_channels) for positive pairs only.

    Each row is (protein1, protein2, {channel: score/1000}).
    Channels absent from the file are silently dropped from available_channels.
    """
    rows = []
    with open(path) as fh:
        reader = csv.DictReader(fh)
        available = [c for c in channels if c in (reader.fieldnames or [])]
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


def build_matrices(rows, embeddings, channels):
    """Build (X, {channel: A}) for pairs where both proteins have embeddings."""
    X_list = []
    channel_lists = {c: [] for c in channels}

    for p1, p2, scores in rows:
        if p1 not in embeddings or p2 not in embeddings:
            continue
        a, b = (p1, p2) if p1 <= p2 else (p2, p1)
        X_list.append(np.concatenate([embeddings[a], embeddings[b]]))
        for c in channels:
            channel_lists[c].append(scores[c])

    X = np.array(X_list, dtype=np.float32)
    A_by_channel = {c: np.array(v, dtype=np.float32) for c, v in channel_lists.items()}
    return X, A_by_channel


def main():
    ap = argparse.ArgumentParser(
        description="Predict STRING evidence channel scores from pair embeddings."
    )
    ap.add_argument("--train",      required=True,
                    help="Training CSV with 'protein1', 'protein2', 'label' and STRING columns")
    ap.add_argument("--test",       required=True,
                    help="Test CSV in the same format")
    ap.add_argument("--embeddings", required=True,
                    help="Pre-computed embeddings .npz (output of embed_sequences.py)")
    ap.add_argument("--seed",       type=int, default=42)
    ap.add_argument("--out",        default="string_channel_analysis.tsv",
                    help="Output TSV path (default: string_channel_analysis.tsv)")
    args = ap.parse_args()

    print("Loading embeddings ...", file=sys.stderr)
    embeddings = load_embeddings(args.embeddings)

    print(f"Loading train positives from {args.train} ...", file=sys.stderr)
    train_rows, available = load_positives(args.train, STRING_CHANNELS)
    if not available:
        sys.exit(f"None of the STRING channels found in {args.train}.")
    print(f"  {len(train_rows)} positive pairs, channels: {available}", file=sys.stderr)

    print(f"Loading test positives from {args.test} ...", file=sys.stderr)
    test_rows, _ = load_positives(args.test, available)
    print(f"  {len(test_rows)} positive pairs", file=sys.stderr)

    X_train, A_train = build_matrices(train_rows, embeddings, available)
    X_test,  A_test  = build_matrices(test_rows,  embeddings, available)
    print(f"  {X_train.shape[0]} train / {X_test.shape[0]} test pairs after embedding lookup",
          file=sys.stderr)

    header = (f"{'Channel':<30} {'Train ρ':>9} {'Test ρ':>9} "
              f"{'p-value':>12} {'Score mean (tr)':>16} {'Score std (tr)':>15}")
    print(f"\n{header}")
    print("-" * len(header))

    results = []
    for channel in available:
        A_tr = A_train[channel]
        A_te = A_test[channel]

        if A_tr.std() == 0:
            print(f"{channel:<30} {'(constant in train, skipped)'}")
            continue

        reg = Ridge(alpha=1.0, solver="sag", random_state=args.seed, max_iter=100)
        reg.fit(X_train, A_tr)

        train_rho, _     = spearmanr(A_tr, reg.predict(X_train))
        test_rho,  pval  = spearmanr(A_te, reg.predict(X_test))

        print(f"{channel:<30} {train_rho:>9.4f} {test_rho:>9.4f} "
              f"{pval:>12.4e} {A_tr.mean():>16.4f} {A_tr.std():>15.4f}")
        results.append((channel, train_rho, test_rho, pval, float(A_tr.mean()), float(A_tr.std())))

    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["channel", "train_spearman_rho", "test_spearman_rho",
                         "p_value", "score_mean_train", "score_std_train"])
        writer.writerows(results)
    print(f"\nResults written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
