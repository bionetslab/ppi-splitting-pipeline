#!/usr/bin/env python3
"""Filter a shared (multi-dataset) FETCH_DATA batch down to one dataset's own proteins.

Used after a single UniProt fetch covering the union of several PPI datasets'
proteins, so BLAST (whose background/E-value statistics depend on exactly
which proteins are in the search database) and per-dataset diagnostics like
the same_species bias check still see only this dataset's own protein set.

With --instances it does the same for DDI mode, mapping the CSV's families to
their domain instances through instances.tsv.
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import INSTANCE_COLUMNS, read_fasta, read_instances, read_ppis, write_fasta


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


def write_instances(rows, path):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(INSTANCE_COLUMNS)
        for row in rows:
            writer.writerow([row[c] for c in INSTANCE_COLUMNS])


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
    ap.add_argument(
        "--instances",
        help="DDI mode: shared instances.tsv (family -> instance ids)",
    )
    ap.add_argument("--out_instances", help="DDI mode: this dataset's instances.tsv")
    args = ap.parse_args()

    ids = get_protein_ids(args.ppis)

    if args.instances:
        rows = [row for row in read_instances(args.instances) if row["family"] in ids]
        # species.tsv has a row per instance and per family; the ILP sampler reads the family rows.
        seq_ids = {row["instance_id"] for row in rows}
        species_ids = seq_ids | ids
        print(
            f"Subsetting shared fetch batch to {len(seq_ids)} instances "
            f"across {len(ids)} families for this dataset...",
            file=sys.stderr,
        )
        write_instances(rows, args.out_instances)

        no_instance = ids - {row["family"] for row in rows}
        if no_instance:
            print(
                f"Warning: {len(no_instance)} families from this dataset have no sampled "
                f"instance in the shared batch: {sorted(no_instance)}",
                file=sys.stderr,
            )
    else:
        seq_ids = species_ids = ids
        print(
            f"Subsetting shared fetch batch to {len(ids)} proteins for this dataset...",
            file=sys.stderr,
        )

    seqs = read_fasta(args.sequences)
    write_fasta(seqs, seq_ids, args.out_sequences)

    filter_tsv(args.go_annotations, args.out_go_annotations, seq_ids)
    filter_tsv(args.species, args.out_species, species_ids)
    filter_tsv(args.lengths, args.out_lengths, seq_ids)

    missing = seq_ids - set(seqs)
    if missing:
        print(
            f"Warning: {len(missing)} proteins from this dataset were not found in "
            f"the shared fetch batch: {sorted(missing)}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
