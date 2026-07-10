#!/usr/bin/env python3
"""
Assign PPIs to train / val / test using KaHIP partition output.

Rules:
  - PPIs whose two proteins belong to different partitions are discarded.
  - Partitions are ranked by their intra-partition PPI count.
  - Largest → train, second → val, smallest → test.
"""

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import mqc_category, read_fasta, read_node_mapping, read_partition, read_ppis, write_fasta, write_ppi_csv


def write_mqc(split_results, id_, n_ppis_discarded):
    """Write the PPI-Partitioning-bar MultiQC contribution for a dataset whose
    split assignment is NOT followed by CD-HIT redundancy removal -- i.e.
    only ever called from sort_ppis_random.py (split_method=random), where
    n_ppis_discarded is trivially 0 (random never discards a PPI). For
    kahip/ilp, REMOVE_REDUNDANT is the sole contributor to this same
    per-dataset "split_bar_{id_}" id instead, since it alone knows the
    post-CD-HIT counts needed for the "Discarded (CD-HIT-2D)" series."""
    with open("sort_ppis_bar_mqc.tsv", "w") as fh:
        fh.write(
            f"# id: 'split_bar_{id_}'\n"
            f"# section_name: 'PPI Partitioning: {id_}'\n"
            "# description: 'PPI counts per split. The discarded bar is coloured by why "
            "a PPI never made it into a split: cross-partition (KaHIP/ILP) or removed by "
            "CD-HIT-2D redundancy filtering.'\n"
            "# plot_type: 'bargraph'\n"
            "# pconfig:\n"
            f"#     id: 'split_bar_plot_{id_}'\n"
            f"#     title: 'PPI Partitioning: edges per split ({id_})'\n"
            "#     ylab: '# PPIs'\n"
            "Sample\tKept\tDiscarded (KaHIP/ILP)\tDiscarded (CD-HIT-2D)\n"
        )
        for r in split_results:
            fh.write(f"{mqc_category(r['name'])}\t{r['n_ppis']}\t0\t0\n")
        fh.write(f"{mqc_category('discarded')}\t0\t{n_ppis_discarded}\t0\n")



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ppis", required=True)
    ap.add_argument("--partition", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--node_mapping", required=True)
    args = ap.parse_args()

    node_to_prot = read_node_mapping(args.node_mapping)
    partition_list = read_partition(args.partition)

    # node_id is 1-indexed; partition_list is 0-indexed
    prot_to_part = {
        node_to_prot[nid]: partition_list[nid - 1]
        for nid in node_to_prot
        if nid - 1 < len(partition_list)
    }

    ppis = read_ppis(args.ppis)
    seqs = read_fasta(args.fasta)

    # Bucket intra-partition PPIs
    part_ppis = defaultdict(list)
    for row in ppis:
        p1, p2 = row["protein1"], row["protein2"]
        part1 = prot_to_part.get(p1)
        part2 = prot_to_part.get(p2)
        if part1 is None or part2 is None or part1 != part2:
            continue
        part_ppis[part1].append(row)

    # Rank partitions by PPI count (largest → train)
    ranked = sorted(part_ppis.keys(), key=lambda p: len(part_ppis[p]), reverse=True)
    if len(ranked) < 3:
        print(
            f"Warning: only {len(ranked)} non-empty partition(s) found; "
            "remaining splits will be empty",
            file=sys.stderr,
        )

    split_names = ["train", "val", "test"]
    split_results = []
    for name, part in zip(split_names, ranked + [None] * (3 - len(ranked))):
        rows = part_ppis[part] if part is not None else []
        proteins = {p for row in rows for p in (row["protein1"], row["protein2"])}
        write_ppi_csv(rows, f"{name}.csv")
        write_fasta(seqs, proteins, f"{name}.fasta")
        print(f"{name}: {len(rows)} PPIs, {len(proteins)} proteins", file=sys.stderr)
        split_results.append({"name": name, "n_ppis": len(rows), "n_proteins": len(proteins)})

    n_ppis_assigned = sum(r["n_ppis"] for r in split_results)
    print(
        f"{len(ppis) - n_ppis_assigned} of {len(ppis)} PPIs discarded (cross-partition); "
        "PPI Partitioning chart is written by REMOVE_REDUNDANT, which runs next.",
        file=sys.stderr,
    )



if __name__ == "__main__":
    main()
