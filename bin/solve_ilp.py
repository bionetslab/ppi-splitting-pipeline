#!/usr/bin/env python3
"""
Assign KaHIP partitions to train/val/test splits by solving an ILP.

Reads a KaHIP partition file plus its node_mapping.tsv, a PPI CSV and a
protein FASTA. Partitions ("clusters") are assigned to splits to minimise
discarded cross-cluster PPIs while keeping each split within epsilon of its
target protein fraction.
"""

import argparse
import sys
from collections import defaultdict

import cvxpy as cp
import numpy as np

from utils import (
    read_fasta,
    read_node_mapping,
    read_partition,
    read_ppis,
    write_fasta,
    write_ppi_csv,
)


def parse_kahip_partition(partition_path, node_mapping_path):
    """Return {protein_id: cluster_id} from a KaHIP partition + node mapping."""
    node_to_prot = read_node_mapping(node_mapping_path)
    partition_list = read_partition(partition_path)
    return {
        node_to_prot[nid]: partition_list[nid - 1]
        for nid in node_to_prot
        if nid - 1 < len(partition_list)
    }


def build_matrices(clusters_list, protein_to_cluster, ppi_rows):
    n = len(clusters_list)
    cluster_to_idx = {c: i for i, c in enumerate(clusters_list)}

    weights = np.zeros(n, dtype=np.float64)
    for p, c in protein_to_cluster.items():
        if c in cluster_to_idx:
            weights[cluster_to_idx[c]] += 1

    intra_ppi = np.zeros(n, dtype=np.float64)
    cross_ppi = np.zeros((n, n), dtype=np.float64)  # upper triangle; cross_ppi[i,j] = count for i < j
    for row in ppi_rows:
        p1, p2 = row["protein1"], row["protein2"]
        c1 = protein_to_cluster.get(p1)
        c2 = protein_to_cluster.get(p2)
        if c1 is None or c2 is None:
            continue
        i, j = cluster_to_idx.get(c1), cluster_to_idx.get(c2)
        if i is None or j is None:
            continue
        if i == j:
            intra_ppi[i] += 1
        else:
            cross_ppi[min(i, j), max(i, j)] += 1

    return weights, cross_ppi, intra_ppi


def solve_ilp(clusters_list, intra_ppi, cross_ppi, splits, names, epsilon, max_sec, solver):
    """
    Variables: x[s, c] ∈ {0,1}  — cluster c assigned to split s.

    Constraints:
        (1) Σ_{s=1}^S x[s,c]  = 1
        ∀c  (each cluster in exactly one split)

        (2) Σ_{i=1}^n x[s,c_i] * k(c_i,c_i) + Σ_{i=1}^{n-1}Σ_{j=i+1}^{n} x[s,c_i] * x[s,c_j] * k(c_i,c_j)
        ≥ (1-ε) * f_s * Σ_{s=1}^S Σ_{i=1}^{n} Σ_{j=i}^{n} x[s,c_i] * x[s,c_j] * k(c_i,c_j)
        ∀s  (minimum split size)

        with k(c_i,c_j) := number of PPIs between clusters c_i and c_j and f_s := fraction of PPIs in split s.
        intra_ppi[i] = k(c_i,c_i); cross_ppi[i,j] = k(c_i,c_j) for i < j (upper triangle, actual counts).

    Objective (minimize): the data loss, i.e., PPIs between clusters assigned to different splits.
        min_X Σ_{i=1}^{n-1}Σ_{j=i+1}^{n} k(c_i,c_j) * (1 - Σ_{s=1}^S x[s,c_i] * x[s,c_j])
    """
    n_splits   = len(splits)
    n_clusters = len(clusters_list)

    # Matrix variable: x[s, c] = 1 iff cluster c is assigned to split s.
    x = cp.Variable((n_splits, n_clusters), boolean=True)

    # Constraint 1: every cluster is in exactly one split.
    constraints = [cp.sum(x, axis=0) == np.ones(n_clusters)]

    # Constraint 2: each split receives ≥ (1-ε)·f_s of all selected PPIs.
    #
    # PPIs in split s = intra-cluster PPIs of clusters in s
    #                 + cross-cluster PPIs where BOTH clusters are in s
    #
    # The product x[s,i]·x[s,j] (both clusters in same split) is linearised:
    #   introduce z[s,k] ∈ {0,1} for each pair k=(i,j) with k(c_i,c_j) > 0
    #   z[s,k] ≤ x[s,i],  z[s,k] ≤ x[s,j],  z[s,k] ≥ x[s,i] + x[s,j] − 1
    loss_pairs = [
        (i, j)
        for i in range(n_clusters)
        for j in range(i + 1, n_clusters)
        if cross_ppi[i, j] > 0
    ]
    cross_counts = np.array([cross_ppi[i, j] for i, j in loss_pairs])  # actual PPI counts

    if loss_pairs:
        z = cp.Variable((n_splits, len(loss_pairs)), boolean=True)
        for k, (i, j) in enumerate(loss_pairs):
            for s in range(n_splits):
                constraints += [
                    z[s, k] <= x[s, i],
                    z[s, k] <= x[s, j],
                    z[s, k] >= x[s, i] + x[s, j] - 1,
                ]
        # total_assigned: intra PPIs (always kept) + co-assigned cross-cluster PPIs
        z_sum          = cp.sum(z, axis=0)  # z_sum[k] = 1 iff pair k ends up in the same split
        total_assigned = float(np.sum(intra_ppi)) + cross_counts @ z_sum
    else:
        z              = None
        total_assigned = float(np.sum(intra_ppi))

    for s, frac in enumerate(splits):
        ppi_in_s = cp.sum(cp.multiply(intra_ppi, x[s]))
        if z is not None:
            ppi_in_s = ppi_in_s + cross_counts @ z[s]
        constraints.append((1.0 - epsilon) * frac * total_assigned <= ppi_in_s)

    # Objective: minimise discarded cross-cluster PPIs.
    # Since each cluster is in exactly one split (constraint 1),
    # cp.max(x[s,i] − x[s,j]) over s = 1 iff i and j are in different splits, 0 otherwise,
    # which equals (1 − Σ_s x[s,i]·x[s,j]) from the docstring.
    if loss_pairs:
        dl_terms = [
            cross_ppi[i, j] * cp.max(cp.vstack([x[s, i] - x[s, j] for s in range(n_splits)]))
            for (i, j) in loss_pairs
        ]
        objective = cp.Minimize(cp.sum(dl_terms))
    else:
        objective = cp.Minimize(cp.sum(x) * 0.0)

    problem = cp.Problem(objective, constraints)

    kwargs = dict(time_limit=max_sec, verbose=True)
    if solver:
        problem.solve(solver=solver, **kwargs)
    else:
        problem.solve(**kwargs)

    if problem.status not in cp.settings.SOLUTION_PRESENT:
        print(f"Solver status: {problem.status}", file=sys.stderr)
        return None

    if problem.status == cp.settings.USER_LIMIT:
        print(
            "Solver hit the time limit before proving optimality; "
            "using the best incumbent found (suboptimal).",
            file=sys.stderr,
        )

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
    ap.add_argument("--partition",    required=True, help="KaHIP partition file")
    ap.add_argument("--node_mapping", required=True, help="KaHIP node_mapping.tsv (node_id -> protein_id)")
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
    ppi_rows = read_ppis(args.ppis)

    print("Reading FASTA …", file=sys.stderr)
    seqs = read_fasta(args.fasta)
    all_proteins = sorted(
        {p for row in ppi_rows for p in (row["protein1"], row["protein2"])} & set(seqs)
    )
    print(f"  {len(all_proteins):,} proteins with sequences", file=sys.stderr)

    print("Parsing KaHIP partition …", file=sys.stderr)
    protein_to_cluster = parse_kahip_partition(args.partition, args.node_mapping)
    protein_to_cluster = {p: protein_to_cluster[p] for p in all_proteins if p in protein_to_cluster}

    clusters_list = sorted(set(protein_to_cluster.values()))
    n_clusters = len(clusters_list)
    cluster_counts = defaultdict(int)
    for v in protein_to_cluster.values():
        cluster_counts[v] += 1
    sizes = sorted(cluster_counts.values(), reverse=True)
    print(f"  {n_clusters:,} clusters; largest has {sizes[0]:,} proteins, "
          f"median {sizes[len(sizes)//2]:,}", file=sys.stderr)

    print("Building problem matrices …", file=sys.stderr)
    weights, cross_ppi, intra_ppi = build_matrices(clusters_list, protein_to_cluster, ppi_rows)
    n_loss_pairs = int(np.sum(cross_ppi > 0))
    total_cross  = int(np.sum(cross_ppi))
    print(f"  {n_loss_pairs:,} cluster pairs with cross-cluster PPIs "
          f"({total_cross:,} PPIs at risk)", file=sys.stderr)

    print("Solving ILP …", file=sys.stderr)
    assignment = solve_ilp(
        clusters_list, intra_ppi, cross_ppi,
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