#!/usr/bin/env python3
"""Turn a precomputed labelled split (e.g. a released train.csv) back into SPLIT_POSITIVES-style output.

Keeps the label == 1 rows without the label column (the negative samplers add their own and refuse an existing one)
and writes the FASTA of the split's proteins, so --precomputed_split can re-sample negatives on fixed positives.

Example:
    extract_positives.py train.csv sequences.fasta train_pos.csv train_pos.fasta
"""

import argparse
import csv
import sys

from utils import read_fasta, write_fasta, write_ppi_csv


def extract(labelled_csv):
    with open(labelled_csv) as fh:
        rows = [r for r in csv.DictReader(fh) if r["label"].strip() == "1"]
    for r in rows:
        del r["label"]
        r["protein1"], r["protein2"] = r["protein1"].strip(), r["protein2"].strip()
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("labelled_csv")
    ap.add_argument("sequences")
    ap.add_argument("out_csv")
    ap.add_argument("out_fasta")
    args = ap.parse_args()

    rows = extract(args.labelled_csv)
    if not rows:
        sys.exit(f"{args.labelled_csv}: no label == 1 rows")
    proteins = {p for r in rows for p in (r["protein1"], r["protein2"])}
    seqs = read_fasta(args.sequences)
    missing = proteins - seqs.keys()
    if missing:
        sys.exit(f"{args.labelled_csv}: {len(missing)} proteins have no sequence in {args.sequences}, "
                 f"e.g. {sorted(missing)[:3]}")
    write_ppi_csv(rows, args.out_csv)
    write_fasta(seqs, proteins, args.out_fasta)
    print(f"{args.labelled_csv}: {len(rows):,} positives, {len(proteins):,} proteins", file=sys.stderr)


if __name__ == "__main__":
    main()
