#!/usr/bin/env python3
"""
Summarise classifier performance across the --robustness_test seed combinations.

Reads a manifest (file, split_seed, sampling_seed, training_seed) pointing at
TRAIN_CLASSIFIER's classifier_metrics_<test_set>_mqc.tsv files and writes:

- robustness_runs.tsv: one row per (seed combination, test set) with every metric.
- robustness_summary_mqc.tsv: per (test set, metric), the mean and standard
  deviation over all runs, decomposed into the three nested seed sources:

    Var_total = Var_split + Var_sampling + Var_training

  where Var_split is the variance of the per-split-seed means, Var_sampling the
  mean variance of the per-(split, sampling)-seed means within a split seed, and
  Var_training the mean variance over training seeds within a (split, sampling)
  pair. All variances are population variances (ddof=0), so for the full
  factorial design the three components add up exactly to Var_total.
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

METRICS = ["AUROC", "AUPRC", "F1", "MCC", "Precision", "Recall", "Accuracy"]
SOURCES = ["split", "sampling", "training"]


def read_metrics(path):
    """Return {metric: value} from a single-row classifier_metrics_*_mqc.tsv."""
    with open(path) as fh:
        rows = list(csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t"))
    if len(rows) != 1:
        sys.exit(f"{path}: expected one data row, found {len(rows)}")
    return {m: float(rows[0][m]) for m in METRICS}


def read_runs(manifest):
    """Return a list of dicts: split_seed, sampling_seed, training_seed, test_set, <metrics>."""
    runs = []
    with open(manifest) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            match = re.fullmatch(r"classifier_metrics_(.+)_mqc\.tsv", Path(row["file"]).name)
            if not match:
                sys.exit(f"Unexpected metrics file name: {row['file']}")
            runs.append(
                {
                    "split_seed": int(row["split_seed"]),
                    "sampling_seed": int(row["sampling_seed"]),
                    "training_seed": int(row["training_seed"]),
                    "test_set": match.group(1),
                    **read_metrics(row["file"]),
                }
            )
    return runs


def decompose(runs, metric):
    """Return (mean, var_total, {source: variance component}) for one metric over runs."""
    by_split_sampling = defaultdict(list)
    for r in runs:
        by_split_sampling[(r["split_seed"], r["sampling_seed"])].append(r[metric])

    var_training = np.mean([np.var(v) for v in by_split_sampling.values()])

    sampling_means_by_split = defaultdict(list)
    for (split_seed, _), v in by_split_sampling.items():
        sampling_means_by_split[split_seed].append(np.mean(v))
    var_sampling = np.mean([np.var(v) for v in sampling_means_by_split.values()])
    var_split = np.var([np.mean(v) for v in sampling_means_by_split.values()])

    values = [r[metric] for r in runs]
    components = {"split": var_split, "sampling": var_sampling, "training": var_training}
    return float(np.mean(values)), float(np.var(values)), components


def write_runs(runs, id_, path):
    cols = ["split_seed", "sampling_seed", "training_seed", "test_set"] + METRICS
    with open(path, "w") as fh:
        fh.write("ID\t" + "\t".join(cols) + "\n")
        for r in sorted(runs, key=lambda r: tuple(r[c] for c in cols[:4])):
            fh.write(id_ + "\t" + "\t".join(str(r[c]) for c in cols) + "\n")


def write_summary_mqc(runs, id_, path):
    seeds = {s: sorted({r[f"{s}_seed"] for r in runs}) for s in SOURCES}
    seed_text = ", ".join(f"{s} {seeds[s]}" for s in SOURCES)
    with open(path, "w") as fh:
        fh.write(
            "# id: 'robustness_summary'\n"
            "# section_name: 'Seed Robustness'\n"
            "# description: 'Baseline classifier performance over every combination of split, negative-sampling "
            f"and training seeds ({seed_text}). SD is the population standard deviation over all runs; "
            "SD split/sampling/training are the nested components (their squares add up to SD squared), "
            "and the % columns give the share of the variance per source.'\n"
            "# plot_type: 'table'\n"
            "# pconfig:\n"
            "#     id: 'robustness_summary_table'\n"
            "#     title: 'Seed robustness of the RF baseline'\n"
            "# headers:\n"
            "#     Runs:            {format: '{:,.0f}'}\n"
            "#     Mean:            {format: '{:.4f}'}\n"
            "#     SD:              {format: '{:.4f}'}\n"
            "#     SD_split:        {format: '{:.4f}', title: 'SD split'}\n"
            "#     SD_sampling:     {format: '{:.4f}', title: 'SD sampling'}\n"
            "#     SD_training:     {format: '{:.4f}', title: 'SD training'}\n"
            "#     Pct_split:       {format: '{:.0f}', suffix: '%', max: 100, title: '% split'}\n"
            "#     Pct_sampling:    {format: '{:.0f}', suffix: '%', max: 100, title: '% sampling'}\n"
            "#     Pct_training:    {format: '{:.0f}', suffix: '%', max: 100, title: '% training'}\n"
            "Sample\tRuns\tMean\tSD\tSD_split\tSD_sampling\tSD_training\tPct_split\tPct_sampling\tPct_training\n"
        )
        for test_set in sorted({r["test_set"] for r in runs}):
            subset = [r for r in runs if r["test_set"] == test_set]
            for metric in METRICS:
                mean, var_total, comps = decompose(subset, metric)
                pcts = [100 * comps[s] / var_total if var_total > 0 else 0.0 for s in SOURCES]
                fields = [len(subset), mean, np.sqrt(var_total)] + [np.sqrt(comps[s]) for s in SOURCES] + pcts
                fh.write(f"{id_} | {test_set} | {metric}\t" + "\t".join(f"{v:.6g}" for v in fields) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, help="TSV: file, split_seed, sampling_seed, training_seed")
    ap.add_argument("--id", required=True, help="Dataset ID")
    args = ap.parse_args()

    runs = read_runs(args.manifest)
    n_combinations = len({(r["split_seed"], r["sampling_seed"], r["training_seed"]) for r in runs})
    print(f"{args.id}: {n_combinations} seed combinations, {len(runs)} metric rows", file=sys.stderr)

    write_runs(runs, args.id, "robustness_runs.tsv")
    write_summary_mqc(runs, args.id, "robustness_summary_mqc.tsv")


if __name__ == "__main__":
    main()
