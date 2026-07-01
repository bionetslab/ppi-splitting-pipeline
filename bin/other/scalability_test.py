#!/usr/bin/env python3
"""
Scalability test: line-graph KaHiP partitioning at increasing PPI dataset sizes.

Line graph construction
-----------------------
Nodes  : unique PPIs (unordered protein pairs)
Edges  : between every pair of PPI-nodes with weight > 0, where weight is:
    1.0   if the two PPIs share a protein
    max(normalised bitscore across all four protein–protein combinations)
          if a BLAST hit exists between any protein in PPI_i and any in PPI_j
    (pairs with neither a shared protein nor any BLAST hit are omitted)

Normalisation: raw bitscore / max bitscore across the full BLAST file → [0, 1].
METIS requires integer weights; floats are multiplied by 1000 and clamped to ≥1.

kaffpa is called with k=10. Results (n_nodes, n_edges, build time, kaffpa time)
are printed and saved to scalability_results.csv in --outdir.
"""

import argparse
import csv
import random
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

WEIGHT_SCALE = 1000


def load_ppis(path):
    """Return a list of unique canonical (min, max) protein pairs."""
    seen = set()
    ppis = []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            p1, p2 = row["protein1"].strip(), row["protein2"].strip()
            if p1 == p2:
                continue
            key = (min(p1, p2), max(p1, p2))
            if key not in seen:
                seen.add(key)
                ppis.append(key)
    return ppis


def load_blast(path):
    """
    Return ({(pa, pb): max_bitscore}, global_max_bitscore)
    where pa < pb lexicographically.
    """
    sim = {}
    global_max = 0.0
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip().split("\t")
            if len(cols) < 4:
                continue
            q = cols[0]
            s = cols[1].split("|")[1] if "|" in cols[1] else cols[1]
            if q == s:
                continue
            try:
                bs = float(cols[3])
            except ValueError:
                continue
            key = (min(q, s), max(q, s))
            if bs > sim.get(key, 0.0):
                sim[key] = bs
            if bs > global_max:
                global_max = bs
    return sim, global_max


def _flat_idx(i, j, n):
    """Index into an upper-triangle flat array for pair (i, j) with i < j."""
    return i * n - i * (i + 1) // 2 + j - i - 1


def build_line_graph(ppis, blast_sim, global_max_bitscore):
    """
    Return (n_nodes, weights) where weights is a float32 numpy array of length
    n*(n-1)/2 storing the upper-triangle of the fully connected line graph.

    All pairs start at 0.0 and are updated to:
      1.0          – PPIs share a protein
      bs/max_bs    – max normalised BLAST bitscore across the four protein pairs
    Pairs with no shared protein and no BLAST hit stay at 0.0.
    """
    n = len(ppis)
    weights = np.zeros(n * (n - 1) // 2, dtype=np.float32)

    protein_to_ppis = defaultdict(list)
    for i, (p1, p2) in enumerate(ppis):
        protein_to_ppis[p1].append(i)
        protein_to_ppis[p2].append(i)

    # Shared-protein edges: weight = 1.0 (maximum possible)
    for ppi_list in protein_to_ppis.values():
        for a in range(len(ppi_list)):
            for b in range(a + 1, len(ppi_list)):
                i, j = ppi_list[a], ppi_list[b]
                if i > j:
                    i, j = j, i
                weights[_flat_idx(i, j, n)] = 1.0

    # BLAST-similarity edges
    if global_max_bitscore > 0:
        proteins_in_sample = set(protein_to_ppis)
        for (pa, pb), bs in blast_sim.items():
            if pa not in proteins_in_sample or pb not in proteins_in_sample:
                continue
            w = np.float32(bs / global_max_bitscore)
            for i in protein_to_ppis[pa]:
                for j in protein_to_ppis[pb]:
                    if i == j:
                        continue
                    if i > j:
                        i, j = j, i
                    k = _flat_idx(i, j, n)
                    if w > weights[k]:
                        weights[k] = w

    return n, weights


def write_metis(n, weights, graph_path):
    """Write a fully connected weighted graph in METIS format.

    Edge weights are floats in [0, 1]; they are scaled by WEIGHT_SCALE and
    clamped to ≥1 so that zero-similarity pairs still carry a positive integer
    weight (METIS requirement).
    """
    m = n * (n - 1) // 2
    with open(graph_path, "w") as fh:
        fh.write(f"{n} {m} 1\n")
        for i in range(n):
            parts = []
            for j in range(n):
                if i == j:
                    continue
                a, b = (i, j) if i < j else (j, i)
                iw = max(1, round(float(weights[_flat_idx(a, b, n)]) * WEIGHT_SCALE))
                parts.append(f"{j + 1} {iw}")
            fh.write(" ".join(parts) + "\n")


def run_kaffpa(graph_path, k, seed, preconfiguration, partition_path):
    cmd = [
        "kaffpa", str(graph_path),
        f"--k={k}",
        f"--seed={seed}",
        f"--output_filename={str(partition_path)}",
        f"--preconfiguration={preconfiguration}",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"kaffpa failed (exit {result.returncode}):\n{result.stderr}")
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ppis",              required=True, help="PPI CSV (protein1,protein2)")
    ap.add_argument("--blast",             required=True, help="BLAST all-vs-all TSV (qseqid sseqid evalue bitscore …)")
    ap.add_argument("--sizes",             nargs="+", type=int,
                    default=[1_000, 5_000, 10_000, 50_000, 100_000],
                    help="PPI subsample sizes to benchmark")
    ap.add_argument("--k",                 type=int, default=10,
                    help="Number of partitions for kaffpa (default: 10)")
    ap.add_argument("--seed",              type=int, default=42)
    ap.add_argument("--preconfiguration",  default="strong",
                    choices=["fast", "eco", "strong", "ultrasetting"],
                    help="kaffpa preconfiguration (default: strong)")
    ap.add_argument("--outdir",            default="scalability_results")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    print("Loading PPIs …", file=sys.stderr)
    all_ppis = load_ppis(args.ppis)
    print(f"  {len(all_ppis):,} unique PPIs loaded", file=sys.stderr)

    print("Loading BLAST similarities …", file=sys.stderr)
    blast_sim, global_max_bitscore = load_blast(args.blast)
    print(f"  {len(blast_sim):,} BLAST edges, max bitscore = {global_max_bitscore:.1f}", file=sys.stderr)

    rng = random.Random(args.seed)

    col_w = [8, 8, 10, 10, 10, 10]
    header = ["Size", "Nodes", "Edges", "Build(s)", "KaHiP(s)", "Total(s)"]
    print("\n" + "  ".join(h.rjust(w) for h, w in zip(header, col_w)))
    print("-" * 62)

    records = []
    for size in sorted(args.sizes):
        if size > len(all_ppis):
            print(f"  [skip] size {size:,} exceeds available PPIs ({len(all_ppis):,})", file=sys.stderr)
            continue

        sample = rng.sample(all_ppis, size)
        n_edges = size * (size - 1) // 2

        try:
            t0 = time.perf_counter()
            n_nodes, weights = build_line_graph(sample, blast_sim, global_max_bitscore)
            graph_path     = outdir / f"line_graph_{size}.metis"
            partition_path = outdir / f"partition_{size}.txt"
            write_metis(n_nodes, weights, graph_path)
            t_build = time.perf_counter() - t0

            t1 = time.perf_counter()
            run_kaffpa(graph_path, args.k, args.seed, args.preconfiguration, partition_path)
            t_kaffpa = time.perf_counter() - t1

            t_total = t_build + t_kaffpa
            records.append(dict(size=size, n_nodes=size, n_edges=n_edges,
                                build_s=t_build, kaffpa_s=t_kaffpa, total_s=t_total))

            row = [str(size), str(size), str(n_edges),
                   f"{t_build:.2f}", f"{t_kaffpa:.2f}", f"{t_total:.2f}"]
            print("  ".join(v.rjust(w) for v, w in zip(row, col_w)))

        except MemoryError:
            print(f"  [OOM] size {size:,}: n*(n-1)/2 = {n_edges:,} edges exceeded available memory")
            records.append(dict(size=size, n_nodes=size, n_edges=n_edges,
                                build_s="OOM", kaffpa_s="OOM", total_s="OOM"))

    csv_path = outdir / "scalability_results.csv"
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["size", "n_nodes", "n_edges", "build_s", "kaffpa_s", "total_s"],
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    print(f"\nResults written to {csv_path}")


if __name__ == "__main__":
    main()