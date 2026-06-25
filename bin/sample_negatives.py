#!/usr/bin/env python3
"""
Sample negative PPI pairs for train, val, test, and a second realistic test set.

Train / val / test  — balanced 50/50, degree-preserving in expectation:
  Both endpoints of each candidate negative are drawn independently from a
  stub pool where protein p appears degree(p) times.  This makes the expected
  negative degree of every protein equal to its positive degree without
  enforcing it exactly per protein.  Exactly len(positives) negatives are kept.

Second test set  — realistic 1:10, uniform random:
  10 * len(positives) negatives are sampled by drawing both endpoints
  uniformly at random from the set of proteins in the test positives.
  Written to test_negatives_random.csv.
"""

import argparse
import csv
import json
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


def _positive_set(ppis):
    return {(min(p1, p2), max(p1, p2)) for p1, p2 in ppis}


def _stub_pool(ppis):
    degree = defaultdict(int)
    for p1, p2 in ppis:
        degree[p1] += 1
        degree[p2] += 1
    return [p for p, d in degree.items() for _ in range(d)]


def sample_negatives(ppis, ratio=1, degree_weighted=True, seed=42):
    """Sample ratio * len(ppis) negatives.

    degree_weighted=True  uses a stub pool (protein appears degree(p) times),
                          preserving positive degree in expectation.
    degree_weighted=False draws endpoints uniformly at random.
    """
    rng = random.Random(seed)
    positives = _positive_set(ppis)
    if not ppis:
        return []
    pool = _stub_pool(ppis) if degree_weighted else sorted({p for pair in ppis for p in pair})
    target = ratio * len(ppis)
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


def per_protein_degrees(ppis, negs):
    """Return {protein: {"x": pos_degree, "y": neg_degree}}."""
    pos_deg = defaultdict(int)
    for p1, p2 in ppis:
        pos_deg[p1] += 1
        pos_deg[p2] += 1
    neg_deg = defaultdict(int)
    for p1, p2 in negs:
        neg_deg[p1] += 1
        neg_deg[p2] += 1
    return {p: {"x": pos_deg[p], "y": neg_deg.get(p, 0)} for p in pos_deg}


def write_csv(pairs, path):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["protein1", "protein2"])
        writer.writerows(pairs)


def write_mqc(split_results, n_random_test):
    gs_data = {r["name"]: {"n_negatives": r["n_negatives"]} for r in split_results}
    gs_data["test"]["n_negatives_random"] = n_random_test

    pos_neg_bar = {
        r["name"]: {"Positives": r["n_positives"], "Negatives": r["n_negatives"]}
        for r in split_results
    }
    sections = [
        {
            "id": "neg_generalstats",
            "plot_type": "generalstats",
            "pconfig": [
                {"n_negatives": {"title": "Negatives", "description": "Degree-preserving (in expectation) negatives, 50/50 balanced", "format": "{:,.0f}", "scale": "Oranges"}},
                {"n_negatives_random": {"title": "Negatives (random)", "description": "Randomly sampled test negatives, 1:10 ratio", "format": "{:,.0f}", "scale": "Reds"}},
            ],
            "data": gs_data,
        },
        {
            "id": "pos_neg_bar",
            "section_name": "Positive vs Negative Pairs",
            "description": "Positive and degree-preserving (in expectation) negative pairs per split (50/50 balanced).",
            "plot_type": "bargraph",
            "pconfig": {
                "id": "pos_neg_bar_plot",
                "title": "Positive vs Negative PPIs per Split",
                "ylab": "# Pairs",
            },
            "data": pos_neg_bar,
        },
    ]
    datasets_data, data_labels, extra_series = [], [], []
    for r in split_results:
        points = [
            {"x": v["x"], "y": v["y"], "name": prot}
            for prot, v in r["degree_scatter"].items()
        ]
        max_deg = max((pt[axis] for pt in r["degree_scatter"].values() for axis in ("x", "y")), default=1)
        identity_line = [
            {"x": v, "y": v, "name": "x=y", "color": "#aaaaaa", "marker_size": 3}
            for v in range(max_deg + 1)
        ]
        datasets_data.append({r["name"]: points})
        data_labels.append({"name": r["name"].capitalize(), "x_minrange": max_deg, "y_minrange": max_deg})
        extra_series.append(identity_line)

    sections.append({
        "id": "degree_scatter",
        "section_name": "Degree Distributions",
        "description": (
            "Each point is a protein. "
            "X = positive degree, Y = negative degree (50/50 balanced negatives). "
            "The grey diagonal is the identity line x=y."
        ),
        "plot_type": "scatter",
        "pconfig": {
            "id": "degree_scatter_plot",
            "title": "Positive vs negative degree per protein",
            "xlab": "Positive degree",
            "ylab": "Negative degree",
            "xmin": 0,
            "ymin": 0,
            "data_labels": data_labels,
            "extra_series": extra_series,
        },
        "data": datasets_data,
    })
    with open("sample_negatives_mqc.json", "w") as fh:
        json.dump(sections, fh, indent=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True)
    ap.add_argument("--val",   required=True)
    ap.add_argument("--test",  required=True)
    ap.add_argument("--seed",  type=int, default=42)
    args = ap.parse_args()

    split_results = []
    for name, src, out in [
        ("train", args.train, "train_negatives.csv"),
        ("val",   args.val,   "val_negatives.csv"),
        ("test",  args.test,  "test_negatives.csv"),
    ]:
        ppis = read_ppis(src)
        negs = sample_negatives(ppis, seed=args.seed)
        write_csv(negs, out)
        print(f"{name}: {len(ppis)} positives → {len(negs)} negatives sampled", file=sys.stderr)
        split_results.append({
            "name": name,
            "n_positives": len(ppis),
            "n_negatives": len(negs),
            "degree_scatter": per_protein_degrees(ppis, negs),
        })

    test_ppis = read_ppis(args.test)
    random_negs = sample_negatives(test_ppis, ratio=10, degree_weighted=False, seed=args.seed)
    write_csv(random_negs, "test_negatives_random.csv")
    print(f"test (random 1:10): {len(test_ppis)} positives → {len(random_negs)} negatives sampled", file=sys.stderr)

    write_mqc(split_results, n_random_test=len(random_negs))


if __name__ == "__main__":
    main()