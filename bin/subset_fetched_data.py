#!/usr/bin/env python3
"""Filter a shared (multi-dataset) FETCH_DATA batch down to one dataset's own proteins.

Used after a single UniProt fetch covering the union of several PPI datasets'
proteins, so BLAST (whose background/E-value statistics depend on exactly
which proteins are in the search database) and per-dataset diagnostics like
the same_species bias check still see only this dataset's own protein set.
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import read_fasta, read_ppis, write_fasta


def get_protein_ids(ppis_path):
    proteins = set()
    for row in read_ppis(ppis_path):
        proteins.add(row["protein1"])
        proteins.add(row["protein2"])
    return proteins


def filter_tsv(in_path, out_path, keep_ids, id_col="protein_id"):
    with open(in_path) as fh_in, open(out_path, "w", newline="") as fh_out:
        reader = csv.DictReader(fh_in, delimiter="\t")
        writer = csv.DictWriter(fh_out, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            if row[id_col].strip() in keep_ids:
                writer.writerow(row)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ppis", required=True)
    ap.add_argument("--sequences", required=True)
    ap.add_argument("--go_annotations", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--lengths", required=True)
    ap.add_argument("--out_sequences", required=True)
    ap.add_argument("--out_go_annotations", required=True)
    ap.add_argument("--out_species", required=True)
    ap.add_argument("--out_lengths", required=True)
    args = ap.parse_args()

    proteins = get_protein_ids(args.ppis)
    print(
        f"Subsetting shared fetch batch to {len(proteins)} proteins for this dataset...",
        file=sys.stderr,
    )

    seqs = read_fasta(args.sequences)
    write_fasta(seqs, proteins, args.out_sequences)

    filter_tsv(args.go_annotations, args.out_go_annotations, proteins)
    filter_tsv(args.species, args.out_species, proteins)
    filter_tsv(args.lengths, args.out_lengths, proteins)

    missing = proteins - set(seqs)
    if missing:
        print(
            f"Warning: {len(missing)} proteins from this dataset were not found in "
            f"the shared fetch batch: {sorted(missing)}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
