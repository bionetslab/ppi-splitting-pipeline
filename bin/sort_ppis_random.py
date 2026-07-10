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
from utils import read_fasta, read_ppis, write_fasta, write_ppi_csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ppis",          required=True)
    ap.add_argument("--fasta",         required=True)
    ap.add_argument("--train-split",   type=float, default=0.8)
    ap.add_argument("--val-split",     type=float, default=0.1)
    ap.add_argument("--test-split",    type=float, default=0.1)
    ap.add_argument("--seed",          type=int, default=42)
    args = ap.parse_args()

    ppis = read_ppis(args.ppis)
    all_proteins = {p for row in ppis for p in (row["protein1"], row["protein2"])}
    seqs = read_fasta(args.fasta)

    shuffled = ppis[:]
    random.Random(args.seed).shuffle(shuffled)

    n = len(shuffled)
    n_train = round(n * args.train_split)
    n_val = round(n * args.val_split)
    # test absorbs whatever rounding drift is left over, rather than being
    # computed from --test-split directly.
    buckets = {
        "train": shuffled[:n_train],
        "val":   shuffled[n_train:n_train + n_val],
        "test":  shuffled[n_train + n_val:],
    }

    split_results = []
    for name, rows in buckets.items():
        proteins = {p for row in rows for p in (row["protein1"], row["protein2"])}
        write_ppi_csv(rows, f"{name}.csv")
        write_fasta(seqs, proteins, f"{name}.fasta")
        print(f"{name}: {len(rows)} PPIs, {len(proteins)} proteins", file=sys.stderr)
        split_results.append({"name": name, "n_ppis": len(rows), "n_proteins": len(proteins)})

    write_mqc(len(ppis), len(all_proteins), split_results)


if __name__ == "__main__":
    main()
