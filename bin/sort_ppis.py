#!/usr/bin/env python3
"""
Assign PPIs to train / val / test using KaHIP partition output.

Rules:
  - PPIs whose two proteins belong to different partitions are discarded.
  - Partitions are ranked by their intra-partition PPI count.
  - Largest → train, second → val, smallest → test.
"""

import argparse
import csv
import sys
from collections import defaultdict


def read_node_mapping(path):
    """Return {node_id (int, 1-indexed): protein_id}."""
    mapping = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mapping[int(row["node_id"])] = row["protein_id"]
    return mapping


def read_partition(path):
    """Return list where element i is the partition of node (i+1)."""
    partitions = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                partitions.append(int(line))
    return partitions


def read_ppis(path):
    pairs = []
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            pairs.append((row["protein1"].strip(), row["protein2"].strip()))
    return pairs


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


def write_csv(pairs, path):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["protein1", "protein2"])
        writer.writerows(pairs)


def write_mqc(n_input, split_results):
    n_assigned = sum(r["n_ppis"] for r in split_results)
    n_discarded = n_input - n_assigned

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
            "Sample\tn_ppis_discarded_kahip\n"
        )
        for r in split_results:
            fh.write(f"{r['name']}\t{n_discarded}\n")

    with open("sort_ppis_bar_mqc.tsv", "w") as fh:
        fh.write(
            "# id: 'split_bar'\n"
            "# section_name: 'PPI Partitioning'\n"
            f"# description: 'Of {n_input:,} input PPIs, {n_assigned:,} were assigned to a split and {n_discarded:,} were discarded because their two proteins landed in different KaHIP partitions.'\n"
            "# plot_type: 'bargraph'\n"
            "# pconfig:\n"
            "#     id: 'split_bar_plot'\n"
            "#     title: 'PPI Partitioning: edges per split'\n"
            "#     ylab: '# PPIs'\n"
            "Sample\tPPIs\n"
        )
        for r in split_results:
            fh.write(f"{r['name']}\t{r['n_ppis']}\n")
        fh.write(f"discarded\t{n_discarded}\n")


def write_fasta(seqs, proteins, path):
    with open(path, "w") as fh:
        for p in sorted(proteins):
            if p in seqs:
                fh.write(f">{p}\n{seqs[p]}\n")


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
    for p1, p2 in ppis:
        part1 = prot_to_part.get(p1)
        part2 = prot_to_part.get(p2)
        if part1 is None or part2 is None or part1 != part2:
            continue
        part_ppis[part1].append((p1, p2))

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
        pairs = part_ppis[part] if part is not None else []
        proteins = {p for pair in pairs for p in pair}
        write_csv(pairs, f"{name}.csv")
        write_fasta(seqs, proteins, f"{name}.fasta")
        print(f"{name}: {len(pairs)} PPIs, {len(proteins)} proteins", file=sys.stderr)
        split_results.append({"name": name, "n_ppis": len(pairs), "n_proteins": len(proteins)})

    write_mqc(len(ppis), split_results)


if __name__ == "__main__":
    main()
