#!/usr/bin/env python3
"""
Assign PPIs to train/val/test by pure random shuffling -- no homology- or
topology-aware partitioning, no redundancy removal downstream.

This is a deliberately naive baseline matching how many PPI-splitting
publications split their data: since the same protein can (and typically
does) land in more than one split, a model can pick up a "topology
shortcut" -- a protein's positive-vs-negative degree ratio in the training
set alone becomes predictive of the label -- instead of learning real
interaction features. See bin/bias_analysis.py's "topology_shortcut"
attribute, which quantifies exactly this effect.

Unlike KaHIP partitioning, no PPI is ever discarded here.
"""

import argparse
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sort_ppis import write_mqc
from utils import expand_members, instances_by_family, read_fasta, read_instances, read_ppis, write_fasta, write_ppi_csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ppis", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument(
        "--instances",
        help="instances.tsv (DDI mode only): expands families to domain instances in the split FASTAs",
    )
    ap.add_argument("--train-split", type=float, default=0.8)
    ap.add_argument("--val-split", type=float, default=0.1)
    ap.add_argument("--test-split", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--id", required=True, help="Dataset ID, for MultiQC tagging")
    args = ap.parse_args()

    members = instances_by_family(read_instances(args.instances)) if args.instances else None
    node_label = "families" if members else "proteins"

    ppis = read_ppis(args.ppis)
    seqs = read_fasta(args.fasta)

    shuffled = ppis[:]
    random.Random(args.seed).shuffle(shuffled)

    n = len(shuffled)
    fracs = {"train": args.train_split, "val": args.val_split, "test": args.test_split}
    counts = {name: round(n * f) for name, f in fracs.items()}

    # Rounding drift goes to the largest non-zero split, so a 0-fraction split stays empty.
    active = [name for name, f in fracs.items() if f > 0]
    if active:
        largest = max(active, key=lambda name: fracs[name])
        counts[largest] += n - sum(counts.values())

    buckets, start = {}, 0
    for name in ("train", "val", "test"):
        buckets[name] = shuffled[start : start + counts[name]]
        start += counts[name]

    split_results = []
    for name, rows in buckets.items():
        proteins = {p for row in rows for p in (row["protein1"], row["protein2"])}
        write_ppi_csv(rows, f"{name}.csv")
        write_fasta(seqs, expand_members(proteins, members), f"{name}.fasta")
        print(f"{name}: {len(rows)} PPIs, {len(proteins)} {node_label}", file=sys.stderr)
        split_results.append({"name": name, "n_ppis": len(rows), "n_proteins": len(proteins)})

    # Random splitting never discards a PPI, by construction.
    write_mqc(split_results, args.id, n_ppis_discarded=0)


if __name__ == "__main__":
    main()
