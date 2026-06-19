#!/usr/bin/env python3
"""
Remove cross-partition redundant sequences using CD-HIT-2D output.

cd-hit-2d writes to its output file the sequences from db2 (i2) that are NOT
similar to any sequence in db1 (i).  So the IDs present in each output file
are exactly the sequences to KEEP in the smaller partition.

Removal logic:
  val_nr   = val  ∩ {not similar to train}   (from sim_train_val.out)
  test_nr  = test ∩ {not similar to train}   (from sim_train_test.out)
  train_nr = train (unchanged; it is the reference)

PPIs are filtered so both partners must still be present.
"""

import argparse
import csv
import sys


def fasta_ids(path):
    """Return set of protein IDs present in a FASTA file."""
    ids = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].rstrip().split()[0])
    return ids


def read_fasta(path):
    seqs = {}
    acc = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(parts)
                acc = line[1:].split()[0]
                parts = []
            elif line:
                parts.append(line)
    if acc and parts:
        seqs[acc] = "".join(parts)
    return seqs


def read_ppis(path):
    pairs = []
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            pairs.append((row["protein1"].strip(), row["protein2"].strip()))
    return pairs


def filter_ppis(pairs, keep):
    return [(p1, p2) for p1, p2 in pairs if p1 in keep and p2 in keep]


def write_csv(pairs, path):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["protein1", "protein2"])
        writer.writerows(pairs)


def write_fasta(seqs, proteins, path):
    with open(path, "w") as fh:
        for p in sorted(proteins):
            if p in seqs:
                fh.write(f">{p}\n{seqs[p]}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train_ppis", required=True)
    ap.add_argument("--val_ppis", required=True)
    ap.add_argument("--test_ppis", required=True)
    ap.add_argument("--train_fasta", required=True)
    ap.add_argument("--val_fasta", required=True)
    ap.add_argument("--test_fasta", required=True)
    ap.add_argument("--sim_train_val", required=True)
    ap.add_argument("--sim_train_test", required=True)
    args = ap.parse_args()

    train_seqs = read_fasta(args.train_fasta)
    val_seqs   = read_fasta(args.val_fasta)
    test_seqs  = read_fasta(args.test_fasta)

    train_ppis = read_ppis(args.train_ppis)
    val_ppis   = read_ppis(args.val_ppis)
    test_ppis  = read_ppis(args.test_ppis)

    # IDs to keep (sequences in CD-HIT-2D output = NOT similar to the reference)
    val_keep  = fasta_ids(args.sim_train_val)
    test_keep = fasta_ids(args.sim_train_test)

    train_prot_nr = set(train_seqs)
    val_prot_nr   = set(val_seqs)  & val_keep
    test_prot_nr  = set(test_seqs) & test_keep

    train_ppis_nr = filter_ppis(train_ppis, train_prot_nr)
    val_ppis_nr   = filter_ppis(val_ppis,   val_prot_nr)
    test_ppis_nr  = filter_ppis(test_ppis,  test_prot_nr)

    write_csv(train_ppis_nr, "train_nr.csv")
    write_csv(val_ppis_nr,   "val_nr.csv")
    write_csv(test_ppis_nr,  "test_nr.csv")
    write_fasta(train_seqs, train_prot_nr, "train_nr.fasta")
    write_fasta(val_seqs,   val_prot_nr,   "val_nr.fasta")
    write_fasta(test_seqs,  test_prot_nr,  "test_nr.fasta")

    for name, ppis, prot in [
        ("train_nr", train_ppis_nr, train_prot_nr),
        ("val_nr",   val_ppis_nr,   val_prot_nr),
        ("test_nr",  test_ppis_nr,  test_prot_nr),
    ]:
        print(f"{name}: {len(ppis)} PPIs, {len(prot)} proteins", file=sys.stderr)


if __name__ == "__main__":
    main()
