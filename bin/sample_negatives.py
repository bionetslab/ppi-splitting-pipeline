#!/usr/bin/env python3
"""
Sample degree-preserving negative PPI pairs for train, val, and test sets.

For each set independently:
  - Only proteins that appear in positive pairs are eligible.
  - For each protein p with positive degree d_p, we aim to sample d_p negative
    partners (so that each protein has ~equal positive and negative annotations).
  - Partners are drawn from a stub pool where each protein appears proportional
    to its degree, making the expected negative degree of every protein equal to
    its positive degree (degree-preserving in expectation).
"""

import argparse
import csv
import random
import sys
from collections import defaultdict


def read_ppis(path):
    pairs = []
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            pairs.append((row["protein1"].strip(), row["protein2"].strip()))
    return pairs


def sample_negatives(ppis, seed=42):
    rng = random.Random(seed)
    positives = set()
    degree = defaultdict(int)

    for p1, p2 in ppis:
        key = (min(p1, p2), max(p1, p2))
        positives.add(key)
        degree[p1] += 1
        degree[p2] += 1

    if not degree:
        return []

    # Stub pool: protein appears degree[p] times so sampling is degree-weighted
    stub_pool = []
    for p, d in degree.items():
        stub_pool.extend([p] * d)

    negatives = set()

    for p1 in sorted(degree):
        target = degree[p1]
        found = 0
        max_tries = target * 100

        for _ in range(max_tries):
            if found >= target:
                break
            p2 = rng.choice(stub_pool)
            if p2 == p1:
                continue
            key = (min(p1, p2), max(p1, p2))
            if key in positives or key in negatives:
                continue
            negatives.add(key)
            found += 1

    return sorted(negatives)


def write_csv(pairs, path):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["protein1", "protein2"])
        writer.writerows(pairs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True)
    ap.add_argument("--val",   required=True)
    ap.add_argument("--test",  required=True)
    ap.add_argument("--seed",  type=int, default=42)
    args = ap.parse_args()

    for name, src, out in [
        ("train", args.train, "train_negatives.csv"),
        ("val",   args.val,   "val_negatives.csv"),
        ("test",  args.test,  "test_negatives.csv"),
    ]:
        ppis = read_ppis(src)
        negs = sample_negatives(ppis, seed=args.seed)
        write_csv(negs, out)
        print(
            f"{name}: {len(ppis)} positives → {len(negs)} negatives sampled",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
