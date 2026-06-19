#!/usr/bin/env python3
"""
Build a METIS graph file from BLAST all-vs-all output and sequence lengths.

Edge weights are integers (METIS requirement).  Float scores are multiplied by
WEIGHT_SCALE and rounded so that 3 significant decimal places are preserved.

Outputs:
  similarity.graph  – METIS format (fmt=1: edge weights present)
  node_mapping.tsv  – node_id (1-indexed) <-> protein_id
"""

import argparse
import csv
import sys
from collections import defaultdict

WEIGHT_SCALE = 1000  # float → int conversion factor


def parse_lengths(path):
    lengths = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            lengths[row["protein_id"]] = int(row["length"])
    return lengths


def parse_blast(path, lengths, edge_weight):
    """Return {(prot_a, prot_b): int_weight} with prot_a < prot_b lexicographically."""
    edges = {}
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip().split("\t")
            if len(cols) < 4:
                continue
            q, s, _, bitscore_str = cols[0], cols[1], cols[2], cols[3]
            if q == s:
                continue
            if q not in lengths or s not in lengths:
                continue
            bitscore = float(bitscore_str)
            if edge_weight == "normalized_bitscore":
                weight = bitscore / min(lengths[q], lengths[s])
            else:
                weight = bitscore
            key = (min(q, s), max(q, s))
            if weight > edges.get(key, -1):
                edges[key] = weight
    return edges


def write_metis(proteins, edges, graph_path, mapping_path):
    prot_to_node = {p: i + 1 for i, p in enumerate(proteins)}

    # Build per-node adjacency: node_id -> {neighbor_id: int_weight}
    adj = defaultdict(dict)
    for (pa, pb), weight in edges.items():
        if pa not in prot_to_node or pb not in prot_to_node:
            continue
        na, nb = prot_to_node[pa], prot_to_node[pb]
        iw = max(1, round(weight * WEIGHT_SCALE))
        adj[na][nb] = iw
        adj[nb][na] = iw

    n = len(proteins)
    m = len(edges)

    with open(graph_path, "w") as fh:
        fh.write(f"{n} {m} 1\n")
        for i, _ in enumerate(proteins):
            node = i + 1
            neighbors = adj.get(node, {})
            line = " ".join(f"{nb} {w}" for nb, w in sorted(neighbors.items()))
            fh.write(line + "\n")

    with open(mapping_path, "w") as fh:
        fh.write("node_id\tprotein_id\n")
        for i, p in enumerate(proteins):
            fh.write(f"{i + 1}\t{p}\n")

    print(f"METIS graph: {n} nodes, {m} edges", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("blast_results")
    ap.add_argument("lengths")
    ap.add_argument("output_graph")
    ap.add_argument("node_mapping")
    ap.add_argument(
        "--edge_weight",
        choices=["bitscore", "normalized_bitscore"],
        default="normalized_bitscore",
    )
    args = ap.parse_args()

    lengths = parse_lengths(args.lengths)
    proteins = sorted(lengths.keys())
    edges = parse_blast(args.blast_results, lengths, args.edge_weight)
    write_metis(proteins, edges, args.output_graph, args.node_mapping)


if __name__ == "__main__":
    main()
