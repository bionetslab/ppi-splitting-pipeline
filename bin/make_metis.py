#!/usr/bin/env python3
"""
Build a METIS graph file from BLAST all-vs-all output and sequence lengths.

Edge weights are integers (METIS requirement).  Float scores are multiplied by
WEIGHT_SCALE and rounded so that 3 significant decimal places are preserved.

With --instances (DDI mode) one node is one Pfam clan: hits are scored per
instance pair and a clan pair keeps its strongest instance pair.

Outputs:
  similarity.graph  – METIS format (fmt=1: edge weights present)
  node_mapping.tsv  – node_id (1-indexed) <-> protein_id (clan id in DDI mode)
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from utils import read_instances

WEIGHT_SCALE = 1000  # float → int conversion factor


def parse_lengths(path):
    lengths = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            lengths[row["protein_id"]] = int(row["length"])
    return lengths


def parse_blast(path, lengths, edge_weight, groups=None):
    """Return {(node_a, node_b): float_weight} with node_a < node_b lexicographically.

    `groups` maps a BLAST id to its node (instance -> clan); None makes every id its own node.
    """
    edges = {}
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip().split("\t")
            if len(cols) < 4:
                continue
            q, s, bitscore_str = cols[0], cols[1], cols[3]
            if q == s:
                continue
            if q not in lengths or s not in lengths:
                continue
            if groups is None:
                node_q, node_s = q, s
            else:
                node_q, node_s = groups[q], groups[s]
                # Two instances of one clan would otherwise write a (clan, clan) edge.
                if node_q == node_s:
                    continue
            bitscore = float(bitscore_str)
            if edge_weight == "normalized_bitscore":
                weight = bitscore / min(lengths[q], lengths[s])
            else:
                weight = bitscore
            key = (min(node_q, node_s), max(node_q, node_s))
            if weight > edges.get(key, -1):
                edges[key] = weight
    return edges


def write_metis(nodes, edges, graph_path, mapping_path):
    node_of = {p: i + 1 for i, p in enumerate(nodes)}

    # Build per-node adjacency: node_id -> {neighbor_id: int_weight}
    adj = defaultdict(dict)
    for (pa, pb), weight in edges.items():
        if pa not in node_of or pb not in node_of:
            continue
        na, nb = node_of[pa], node_of[pb]
        iw = max(1, round(weight * WEIGHT_SCALE))
        adj[na][nb] = iw
        adj[nb][na] = iw

    n = len(nodes)
    m = len(edges)

    with open(graph_path, "w") as fh:
        fh.write(f"{n} {m} 1\n")
        for i, _ in enumerate(nodes):
            node = i + 1
            neighbors = adj.get(node, {})
            line = " ".join(f"{nb} {w}" for nb, w in sorted(neighbors.items()))
            fh.write(line + "\n")

    with open(mapping_path, "w") as fh:
        fh.write("node_id\tprotein_id\n")
        for i, p in enumerate(nodes):
            fh.write(f"{i + 1}\t{p}\n")

    print(f"METIS graph: {n} nodes, {m} edges", file=sys.stderr)
    print(f"isolated nodes: {n - len(adj)}/{n}", file=sys.stderr)


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
    ap.add_argument(
        "--instances",
        help="instances.tsv (DDI mode only): one node per Pfam clan instead of per protein",
    )
    args = ap.parse_args()

    lengths = parse_lengths(args.lengths)
    if args.instances:
        groups = {r["instance_id"]: r["clan"] for r in read_instances(args.instances)}
        # Every clan in the table is a node, isolated ones included.
        nodes = sorted(set(groups.values()))
        print(f"clan contraction: {len(groups)} instances -> {len(nodes)} clans", file=sys.stderr)
    else:
        groups = None
        nodes = sorted(lengths.keys())
    edges = parse_blast(args.blast_results, lengths, args.edge_weight, groups)
    write_metis(nodes, edges, args.output_graph, args.node_mapping)


if __name__ == "__main__":
    main()
