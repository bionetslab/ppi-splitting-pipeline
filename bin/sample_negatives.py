#!/usr/bin/env python3
"""
Sample negative PPI pairs for a single split and write a combined labelled CSV.

Degree-weighted (default): both endpoints drawn from a stub pool (protein p
appears degree(p) times), preserving degree in expectation. Used for the
50/50 train/val/test_balanced splits.

Uniform (--uniform): endpoints drawn uniformly at random. Used for the
1:10 test_realistic split.

Output CSV has columns: protein1, protein2, label  (1 = positive, 0 = negative).

Per-split fan-out (train/val/test_balanced/test_realistic) is handled by
Nextflow channel logic, not by this script; each invocation samples exactly
one split.
"""

import argparse
import csv
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import mqc_sample, read_ppis


def _positive_set(rows):
    return {(min(r["protein1"], r["protein2"]), max(r["protein1"], r["protein2"])) for r in rows}


def _stub_pool(rows):
    degree = defaultdict(int)
    for r in rows:
        degree[r["protein1"]] += 1
        degree[r["protein2"]] += 1
    return [p for p, d in degree.items() for _ in range(d)]


def sample_negatives(rows, ratio=1, degree_weighted=True, seed=42):
    """Sample ratio * len(rows) negatives.

    degree_weighted=True  uses a stub pool (protein appears degree(p) times),
                          preserving positive degree in expectation.
    degree_weighted=False draws endpoints uniformly at random.
    """
    rng = random.Random(seed)
    positives = _positive_set(rows)
    if not rows:
        return []
    pool = _stub_pool(rows) if degree_weighted else sorted({p for r in rows for p in (r["protein1"], r["protein2"])})
    target = int(ratio * len(rows))
    negatives = set()
    for _ in range(target * 100):
        if len(negatives) >= target:
            break
        p1, p2 = rng.choice(pool), rng.choice(pool)
        key = (min(p1, p2), max(p1, p2))
        if key in positives or key in negatives:
            continue
        negatives.add(key)
    return sorted(negatives)


def write_combined(pos_rows, negatives, path):
    if pos_rows and "label" in pos_rows[0]:
        raise ValueError(
            "--positives already has a 'label' column -- refusing to relabel "
            "what looks like a previously-combined pos/neg file as all-positive"
        )
    extra_fields = [k for k in (pos_rows[0].keys() if pos_rows else [])
                    if k not in ("protein1", "protein2")]
    fieldnames = ["protein1", "protein2", "label"] + extra_fields
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in pos_rows:
            out = {f: row.get(f, "") for f in extra_fields}
            out["protein1"] = row["protein1"]
            out["protein2"] = row["protein2"]
            out["label"] = 1
            writer.writerow(out)
        for p1, p2 in negatives:
            out = {f: "" for f in extra_fields}
            out["protein1"] = p1
            out["protein2"] = p2
            out["label"] = 0
            writer.writerow(out)


def write_mqc(split_name, n_positives, n_negatives, gs_out, bar_out, id_):
    # Sharing the same 'id' across all datasets' files lets MultiQC merge
    # them into one combined general-stats table / bar plot; Sample is
    # qualified by dataset so rows/categories don't collide across datasets.
    sample = mqc_sample(id_, split_name)
    with open(gs_out, "w") as fh:
        fh.write(
            "# id: 'neg_generalstats'\n"
            "# plot_type: 'generalstats'\n"
            "# pconfig:\n"
            "#     - ID:\n"
            "#         title: 'ID'\n"
            "#         description: 'Dataset ID'\n"
            "#     - n_positives:\n"
            "#         title: 'Positives'\n"
            "#         description: 'Positive PPIs in the split'\n"
            "#         format: '{:,.0f}'\n"
            "#         scale: 'Blues'\n"
            "#     - n_negatives:\n"
            "#         title: 'Negatives'\n"
            "#         description: 'Sampled negatives for the split'\n"
            "#         format: '{:,.0f}'\n"
            "#         scale: 'Oranges'\n"
            "Sample\tID\tn_positives\tn_negatives\n"
            f"{sample}\t{id_}\t{n_positives}\t{n_negatives}\n"
        )

    with open(bar_out, "w") as fh:
        fh.write(
            "# id: 'pos_neg_bar'\n"
            "# section_name: 'Positive vs Negative Pairs'\n"
            "# description: 'Positive and sampled negative pairs per split, across every dataset.'\n"
            "# plot_type: 'bargraph'\n"
            "# pconfig:\n"
            "#     id: 'pos_neg_bar_plot'\n"
            "#     title: 'Positive vs Negative PPIs per Split'\n"
            "#     ylab: '# Pairs'\n"
            "Sample\tPositives\tNegatives\n"
            f"{sample}\t{n_positives}\t{n_negatives}\n"
        )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--positives",  required=True, help="Input CSV of positive PPIs for this split")
    ap.add_argument("--output",     required=True, help="Output combined CSV path")
    ap.add_argument("--split-name", required=True, help="Split label, e.g. train/val/test_balanced/test_realistic")
    ap.add_argument("--ratio",      type=float, default=1.0, help="negatives = ratio * positives")
    ap.add_argument("--uniform", action="store_true",
                     help="Draw negative endpoints uniformly at random instead of degree-weighted")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--id", required=True, help="Dataset ID, for MultiQC tagging")
    args = ap.parse_args()

    rows = read_ppis(args.positives)
    negs = sample_negatives(rows, ratio=args.ratio, degree_weighted=not args.uniform, seed=args.seed)
    write_combined(rows, negs, args.output)
    print(f"{args.split_name}: {len(rows)} positives → {len(negs)} negatives sampled", file=sys.stderr)

    write_mqc(
        args.split_name, len(rows), len(negs),
        gs_out=f"{args.split_name}_gs_mqc.tsv",
        bar_out=f"{args.split_name}_bar_mqc.tsv",
        id_=args.id,
    )


if __name__ == "__main__":
    main()