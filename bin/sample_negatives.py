#!/usr/bin/env python3
"""
Sample negative PPI pairs and write combined labelled CSVs.

train.csv / val.csv / test_balanced.csv — 50/50, degree-preserving in expectation:
  Both endpoints drawn from a stub pool (protein p appears degree(p) times).

test_realistic.csv — 1:10, uniform random:
  10 × len(positives) negatives, endpoints drawn uniformly at random.

All output files have columns: protein1, protein2, label  (1 = positive, 0 = negative).
"""

import argparse
import csv
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import read_ppis


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
    target = ratio * len(rows)
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


def write_mqc(split_results, n_random_test):
    with open("sample_negatives_gs_mqc.tsv", "w") as fh:
        fh.write(
            "# id: 'neg_generalstats'\n"
            "# plot_type: 'generalstats'\n"
            "# pconfig:\n"
            "#     - n_positives:\n"
            "#         title: 'Positives'\n"
            "#         description: 'Positive PPIs in the split'\n"
            "#         format: '{:,.0f}'\n"
            "#         scale: 'Blues'\n"
            "#     - n_negatives:\n"
            "#         title: 'Negatives'\n"
            "#         description: 'Degree-preserving (in expectation) negatives, 50/50 balanced'\n"
            "#         format: '{:,.0f}'\n"
            "#         scale: 'Oranges'\n"
            "#     - n_negatives_random:\n"
            "#         title: 'Negatives (random)'\n"
            "#         description: 'Randomly sampled test negatives, 1:10 ratio'\n"
            "#         format: '{:,.0f}'\n"
            "#         scale: 'Reds'\n"
            "Sample\tn_positives\tn_negatives\tn_negatives_random\n"
        )
        for r in split_results:
            rand = n_random_test if r["name"] == "test" else ""
            fh.write(f"{r['name']}\t{r['n_positives']}\t{r['n_negatives']}\t{rand}\n")

    with open("sample_negatives_bar_mqc.tsv", "w") as fh:
        fh.write(
            "# id: 'pos_neg_bar'\n"
            "# section_name: 'Positive vs Negative Pairs'\n"
            "# description: 'Positive and degree-preserving (in expectation) negative pairs per split (50/50 balanced).'\n"
            "# plot_type: 'bargraph'\n"
            "# pconfig:\n"
            "#     id: 'pos_neg_bar_plot'\n"
            "#     title: 'Positive vs Negative PPIs per Split'\n"
            "#     ylab: '# Pairs'\n"
            "Sample\tPositives\tNegatives\n"
        )
        for r in split_results:
            fh.write(f"{r['name']}\t{r['n_positives']}\t{r['n_negatives']}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True)
    ap.add_argument("--val",   required=True)
    ap.add_argument("--test",  required=True)
    ap.add_argument("--seed",  type=int, default=42)
    args = ap.parse_args()

    split_results = []
    for name, src, out in [
        ("train", args.train, "train.csv"),
        ("val",   args.val,   "val.csv"),
        ("test",  args.test,  "test_balanced.csv"),
    ]:
        rows = read_ppis(src)
        negs = sample_negatives(rows, seed=args.seed)
        write_combined(rows, negs, out)
        print(f"{name}: {len(rows)} positives → {len(negs)} negatives sampled", file=sys.stderr)
        split_results.append({
            "name": name,
            "n_positives": len(rows),
            "n_negatives": len(negs),
        })

    test_rows = read_ppis(args.test)
    random_negs = sample_negatives(test_rows, ratio=10, degree_weighted=False, seed=args.seed)
    write_combined(test_rows, random_negs, "test_realistic.csv")
    print(f"test (random 1:10): {len(test_rows)} positives → {len(random_negs)} negatives sampled", file=sys.stderr)

    write_mqc(split_results, n_random_test=len(random_negs))


if __name__ == "__main__":
    main()