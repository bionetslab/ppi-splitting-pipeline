#!/usr/bin/env python3
"""
Assign CD-HIT clusters to train/val/test splits by solving an ILP.

Reads a CD-HIT .clstr file plus a PPI CSV and protein FASTA.
Clusters are assigned to splits to minimise discarded cross-cluster PPIs
while keeping each split within epsilon of its target protein fraction.
"""

import argparse
import re
import sys
from collections import defaultdict

import cvxpy as cp
import numpy as np

from utils import read_fasta, read_ppis, write_fasta, write_ppi_csv


def parse_clstr(path):
    """Parse a CD-HIT .clstr file into {protein_id: cluster_id}."""
    assignment = {}
    cluster_id = -1
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">Cluster"):
                cluster_id += 1
            else:
                m = re.search(r">(\S+?)\.\.\.", line)
                if m:
                    assignment[m.group(1)] = cluster_id
    return assignment


def build_matrices(clusters_list, protein_to_cluster, ppi_rows):
    n = len(clusters_list)
    cluster_to_idx = {c: i for i, c in enumerate(clusters_list)}

    weights = np.zeros(n, dtype=np.float64)
    for p, c in protein_to_cluster.items():
        if c in cluster_to_idx:
            weights[cluster_to_idx[c]] += 1

    cross_ppi = np.zeros((n, n), dtype=np.float64)
    for row in ppi_rows:
        p1, p2 = row["protein1"], row["protein2"]
        c1 = protein_to_cluster.get(p1)
        c2 = protein_to_cluster.get(p2)
        if c1 is None or c2 is None or c1 == c2:
            continue
        i, j = cluster_to_idx.get(c1), cluster_to_idx.get(c2)
        if i is None or j is None:
            continue
        cross_ppi[i, j] += 1
        cross_ppi[j, i] += 1

    cross_ppi /= 2
    return weights, cross_ppi


def solve_ilp(clusters_list, n_ppis, cross_ppi, splits, names, epsilon, max_sec, solver):
    """
    Variables: x[s, c] ∈ {0,1}  — cluster c assigned to split s.

    Constraints:
        (1) Σ_{s=1}^S x[s,c]  = 1
        ∀c  (each cluster in exactly one split)

        (2) Σ_{i=1}^n x[s,c_i] * k(c_i,c_i) + Σ_{i=1}^{n-1}Σ_{j=i+1}^{n} x[s,c_i] * x[s,c_j] * k(c_i,c_j)
        ≤ (1-ε) * f_s * Σ_{s=1}^S Σ_{i=1}^{n} Σ_{j=i}^{n} x[s,c_i] * x[s,c_j] * k(c_i,c_j)
        ∀s  (minimum split size)

        with k(c_i,c_j) := number of PPIs between clusters c_i and c_j and f_s := fraction of PPIs in split s.

    Objective (minimize): the data loss, i.e., PPIs between clusters assigned to different splits.
        min_X Σ_{i=1}^{n-1}Σ_{j=i+1}^{n} k(c_i,c_j) * (1 - Σ_{s=1}^S x[s,c_i] * x[s,c_j])
    """
    n_splits   = len(splits)
    n_clusters = len(clusters_list)

    # Matrix variable with dimensionality number of splits (3) * number of clusters (n_clusters)
    # showing which clusters are in which splits.
    x = cp.Variable((n_splits, n_clusters), boolean=True)

    # Constraint 1: every cluster is in exactly one split.
    constraints = [cp.sum(x, axis=0) == np.ones(n_clusters)]

    # Constraint 2: each split receives ≥ (f_s − ε) × n_ppis PPIs.
    # LHS = Σ_i k(c_i,c_i)·x[s,i] + Σ_{i<j} k(c_i,c_j)·x[s,i]·x[s,j]
    # The quadratic term x[s,i]·x[s,j] is linearised with z[s,k] ∈ {0,1}
    # for each loss pair k=(i,j) where k(c_i,c_j) > 0:
    #   z[s,k] ≤ x[s,i],  z[s,k] ≤ x[s,j],  z[s,k] ≥ x[s,i] + x[s,j] − 1
    # k(c_i,c_j) = 2·cross_ppi[i,j]  (cross_ppi stores halved counts due to /=2).
    # k(c_i,c_i) = intra-cluster PPIs per cluster, approximated by distributing
    # total_intra = n_ppis − total_cross proportionally to cross-cluster involvement.
    loss_pairs = [
        (e1, e2)
        for e1 in range(n_clusters)
        for e2 in range(e1 + 1, n_clusters)
        if cross_ppi[e1, e2] > 0
    ]

    total_cross  = float(np.sum(cross_ppi))
    total_intra  = max(0.0, float(n_ppis) - total_cross)
    cross_row_sums = np.sum(cross_ppi, axis=1)
    if total_cross > 0:
        intra_ppi = total_intra * (cross_row_sums / total_cross)
    else:
        intra_ppi = np.full(n_clusters, total_intra / max(n_clusters, 1))

    if loss_pairs:
        z = cp.Variable((n_splits, len(loss_pairs)), boolean=True)
        for k, (e1, e2) in enumerate(loss_pairs):
            for s in range(n_splits):
                constraints += [
                    z[s, k] <= x[s, e1],
                    z[s, k] <= x[s, e2],
                    z[s, k] >= x[s, e1] + x[s, e2] - 1,
                ]
        cross_counts = np.array([2.0 * cross_ppi[e1, e2] for e1, e2 in loss_pairs])
    else:
        z            = None
        cross_counts = np.array([])

    for s, frac in enumerate(splits):
        lo = max(0.0, frac - epsilon) * float(n_ppis)
        ppi_in_s = cp.sum(cp.multiply(intra_ppi, x[s]))
        if z is not None:
            ppi_in_s = ppi_in_s + cross_counts @ z[s]
        constraints.append(lo <= ppi_in_s)

    # Objective: minimize discarded cross-cluster PPIs.
    # Since each cluster is in exactly one split (constraint 1),
    # cp.max(x[s,e1] − x[s,e2]) over s = 1 iff e1,e2 are in different splits, 0 otherwise.
    # This is equivalent to k(c_i,c_j) · (1 − Σ_s x[s,c_i]·x[s,c_j]) from the docstring.
    if loss_pairs:
        dl_terms = [
            2.0 * cross_ppi[e1, e2] * cp.max(cp.vstack([x[s, e1] - x[s, e2]
                                                          for s in range(n_splits)]))
            for (e1, e2) in loss_pairs
        ]
        objective = cp.Minimize(cp.sum(dl_terms))
    else:
        objective = cp.Minimize(cp.sum(x) * 0.0)

    problem = cp.Problem(objective, constraints)

    kwargs = dict(time_limit=max_sec)
    if solver:
        problem.solve(solver=solver, **kwargs)
    else:
        problem.solve(**kwargs)

    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        print(f"Solver status: {problem.status}", file=sys.stderr)
        return None

    return {
        clusters_list[c]: names[s]
        for s in range(n_splits)
        for c in range(n_clusters)
        if x[s, c].value is not None and x[s, c].value > 0.5
    }


def write_mqc(n_input, n_proteins_input, split_results, total_cross):
    n_ppis_assigned  = sum(r["n_ppis"] for r in split_results)
    n_ppis_discarded = n_input - n_ppis_assigned

    with open("sort_ppis_gs_mqc.tsv", "w") as fh:
        fh.write(
            "# id: 'split_generalstats'\n"
            "# plot_type: 'generalstats'\n"
            "# pconfig:\n"
            "#     - n_ppis_discarded_ilp:\n"
            "#         title: 'PPIs discarded (ILP)'\n"
            "#         description: 'Cross-cluster PPIs discarded because their proteins were assigned to different splits'\n"
            "#         format: '{:,.0f}'\n"
            "#         scale: 'Greys'\n"
            "Sample\tn_ppis_discarded_ilp\n"
        )
        for r in split_results:
            fh.write(f"{r['name']}\t{n_ppis_discarded}\n")

    with open("sort_ppis_bar_mqc.tsv", "w") as fh:
        fh.write(
            "# id: 'split_bar'\n"
            "# section_name: 'PPI Partitioning'\n"
            f"# description: 'Of {n_input:,} input PPIs ({n_proteins_input:,} proteins), "
            f"{n_ppis_assigned:,} were assigned to a split and {n_ppis_discarded:,} were "
            f"discarded (cross-cluster PPIs). {total_cross:,} cross-cluster PPI pairs penalised in ILP.'\n"
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
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ppis",        required=True, help="PPI CSV (protein1,protein2)")
    ap.add_argument("--fasta",       required=True, help="Protein FASTA")
    ap.add_argument("--clusters",    required=True, help="CD-HIT .clstr file")
    ap.add_argument("--train-split", type=float, default=0.8)
    ap.add_argument("--val-split",   type=float, default=0.1)
    ap.add_argument("--test-split",  type=float, default=0.1)
    ap.add_argument("--epsilon",     type=float, default=0.05,
                    help="Allowed fractional deviation from target split size (default 0.05)")
    ap.add_argument("--max-sec",     type=int,   default=300,
                    help="ILP solver time limit in seconds (default 300)")
    ap.add_argument("--solver",      default=None,
                    help="CVXPY solver name, e.g. SCIP, GLPK_MI (default: auto)")
    args = ap.parse_args()

    splits = [args.train_split, args.val_split, args.test_split]
    names  = ["train", "val", "test"]
    assert abs(sum(splits) - 1.0) < 1e-6, "Split fractions must sum to 1"

    print("Loading PPIs …", file=sys.stderr)
    all_rows = read_ppis(args.ppis)
    seen, ppi_rows = set(), []
    for row in all_rows:
        p1, p2 = row["protein1"], row["protein2"]
        if p1 == p2:
            continue
        key = (min(p1, p2), max(p1, p2))
        if key not in seen:
            seen.add(key)
            ppi_rows.append(row)
    print(f"  {len(ppi_rows):,} unique PPIs", file=sys.stderr)

    print("Reading FASTA …", file=sys.stderr)
    seqs = read_fasta(args.fasta)
    all_proteins = sorted(
        {p for row in ppi_rows for p in (row["protein1"], row["protein2"])} & set(seqs)
    )
    print(f"  {len(all_proteins):,} proteins with sequences", file=sys.stderr)

    print("Parsing CD-HIT clusters …", file=sys.stderr)
    protein_to_cluster = parse_clstr(args.clusters)
    protein_to_cluster = {p: protein_to_cluster[p] for p in all_proteins if p in protein_to_cluster}

    clusters_list = sorted(set(protein_to_cluster.values()))
    n_clusters = len(clusters_list)
    cluster_counts = defaultdict(int)
    for v in protein_to_cluster.values():
        cluster_counts[v] += 1
    sizes = sorted(cluster_counts.values(), reverse=True)
    print(f"  {n_clusters:,} clusters; largest has {sizes[0]:,} proteins, "
          f"median {sizes[len(sizes)//2]:,}", file=sys.stderr)

    if n_clusters > 1000:
        print(f"  Warning: {n_clusters:,} clusters may make the ILP slow.", file=sys.stderr)

    print("Building problem matrices …", file=sys.stderr)
    weights, cross_ppi = build_matrices(clusters_list, protein_to_cluster, ppi_rows)
    n_loss_pairs = int(np.sum(cross_ppi > 0)) // 2
    total_cross  = int(np.sum(cross_ppi))
    print(f"  {n_loss_pairs:,} cluster pairs with cross-cluster PPIs "
          f"({total_cross:,} PPIs at risk)", file=sys.stderr)

    print("Solving ILP …", file=sys.stderr)
    assignment = solve_ilp(
        clusters_list, len(ppi_rows), cross_ppi,
        splits, names, args.epsilon, args.max_sec, args.solver,
    )
    if assignment is None:
        print("ILP did not find a feasible solution.", file=sys.stderr)
        sys.exit(1)

    protein_to_split = {
        p: assignment[c]
        for p, c in protein_to_cluster.items()
        if c in assignment
    }

    split_rows = defaultdict(list)
    for row in ppi_rows:
        s1 = protein_to_split.get(row["protein1"])
        s2 = protein_to_split.get(row["protein2"])
        if s1 is not None and s1 == s2:
            split_rows[s1].append(row)

    split_results = []
    for name in names:
        rows     = split_rows[name]
        proteins = {p for row in rows for p in (row["protein1"], row["protein2"])}
        write_ppi_csv(rows, f"{name}.csv")
        write_fasta(seqs, proteins, f"{name}.fasta")
        print(f"  {name}: {len(rows):,} PPIs, {len(proteins):,} proteins", file=sys.stderr)
        split_results.append({"name": name, "n_ppis": len(rows)})

    write_mqc(len(ppi_rows), len(all_proteins), split_results, total_cross)


if __name__ == "__main__":
    main()