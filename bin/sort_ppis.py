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
from utils import read_fasta, read_node_mapping, read_partition, read_ppis, write_fasta, write_ppi_csv


def write_mqc(n_input, n_proteins_input, split_results):
    n_ppis_assigned     = sum(r["n_ppis"]     for r in split_results)
    n_proteins_assigned = sum(r["n_proteins"] for r in split_results)
    n_ppis_discarded     = n_input         - n_ppis_assigned
    n_proteins_discarded = n_proteins_input - n_proteins_assigned

    with open("sort_ppis_gs_mqc.tsv", "w") as fh:
        fh.write(
            "# id: 'split_generalstats'\n"
            "# plot_type: 'generalstats'\n"
            "# pconfig:\n"
            "#     - n_ppis_discarded_kahip:\n"
            "#         title: 'PPIs discarded (KaHIP)'\n"
            "#         description: 'Cross-partition PPIs discarded during KaHIP partitioning (total, same for all splits)'\n"
            "#         format: '{:,.0f}'\n"
            "#         scale: 'Greys'\n"
            "#     - n_proteins_discarded_kahip:\n"
            "#         title: 'Proteins discarded (KaHIP)'\n"
            "#         description: 'Proteins whose every PPI was cross-partition and thus discarded by KaHIP (total, same for all splits)'\n"
            "#         format: '{:,.0f}'\n"
            "#         scale: 'Greys'\n"
            "Sample\tn_ppis_discarded_kahip\tn_proteins_discarded_kahip\n"
        )
        for r in split_results:
            fh.write(f"{r['name']}\t{n_ppis_discarded}\t{n_proteins_discarded}\n")

    with open("sort_ppis_bar_mqc.tsv", "w") as fh:
        fh.write(
            "# id: 'split_bar'\n"
            "# section_name: 'PPI Partitioning'\n"
            f"# description: 'Of {n_input:,} input PPIs ({n_proteins_input:,} proteins), {n_ppis_assigned:,} PPIs ({n_proteins_assigned:,} proteins) were assigned to a split and {n_ppis_discarded:,} PPIs ({n_proteins_discarded:,} proteins) were discarded because their interactions were cross-partition.'\n"
            "# plot_type: 'bargraph'\n"
            "# pconfig:\n"
            "#     id: 'split_bar_plot'\n"
            "#     title: 'PPI Partitioning: edges per split'\n"
            "#     ylab: '# PPIs'\n"
            "Sample\tPPIs\n"
        )
        for r in split_results:
            fh.write(f"{r['name']}\t{r['n_ppis']}\n")
        fh.write(f"discarded\t{n_ppis_discarded}\n")



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
    all_proteins = {p for row in ppis for p in (row["protein1"], row["protein2"])}
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

    write_mqc(len(ppis), len(all_proteins), split_results)



if __name__ == "__main__":
    main()
