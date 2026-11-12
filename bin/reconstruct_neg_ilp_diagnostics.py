#!/usr/bin/env python3
"""Reconstruct the ILP negative sampler's MultiQC diagnostics table for robustness runs, which don't export it.

The bias terms only depend on the positives, the chosen negatives and the candidate pool the ILP picked from. The
pool isn't saved but is deterministic (sorted protein index, --seed = sampling seed, --max-candidates = neg_ilp_pool_factor (default 4) * #positives,
see processes/negative_sampling.nf), so it is rebuilt with sample_negatives_ilp.py's own functions and every term is
evaluated at the saved negatives instead of re-solving.

Writes one row per (ILP-sampled samplesheet row, split seed, sampling seed, split). test_realistic negatives are
sampled uniformly, so it has no row (as in the pipeline).

Example:
    reconstruct_neg_ilp_diagnostics.py --samplesheet samplesheets/samplesheet.csv --results results_robustness \\
        --out results_robustness/multiqc/multiqc_report_data/multiqc_neg_ilp_diagnostics_table.txt
"""

import argparse
import csv
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sample_negatives_ilp import (  # noqa: E402
    _TERM_KEY,
    DIAG_COLUMNS,
    SamplingConfig,
    _build_go_membership,
    assemble_active_biases,
    build_candidate_set,
    build_context,
    build_protein_index,
    load_candidate_network,
    load_go_bp,
    load_species,
    pos_pairs_from_rows,
)
from utils import mqc_sample  # noqa: E402

SPLITS = ["train", "val", "test_balanced"]  # ILP-sampled, each with neg_ratio 1
LAMBDA_DEFAULT = 1.0  # nextflow.config neg_ilp_lambda_* when the samplesheet leaves it empty


def read_labelled(path):
    pos, neg = [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            pair = {"protein1": row["protein1"].strip(), "protein2": row["protein2"].strip()}
            (pos if row["label"] == "1" else neg).append(pair)
    return pos, neg


def reconstruct(split, csv_path, lambdas, seed, species_path, go_path, candidate_network=None, pool_factor=4):
    """Return one diagnostics row, rebuilding the context exactly like sample_negatives_ilp.sample_negatives_ilp()."""
    pos_rows, neg_rows = read_labelled(csv_path)
    cfg = SamplingConfig(**lambdas, seed=seed, max_candidates=pool_factor * len(pos_rows))
    protein_to_idx, idx_to_protein = build_protein_index(pos_rows)
    pos_pairs = pos_pairs_from_rows(pos_rows, protein_to_idx)

    taxonomy_codes = n_taxa = None
    if cfg.lambda_taxon_pair > 0:
        uniq, taxonomy_codes = np.unique(load_species(species_path, protein_to_idx), return_inverse=True)
        taxonomy_codes, n_taxa = taxonomy_codes.astype(np.int64), int(len(uniq))
    go_bp = go_membership = go_sizes = None
    if cfg.lambda_jaccard > 0:
        go_bp = load_go_bp(go_path, protein_to_idx)
        go_membership, go_sizes = _build_go_membership(go_bp)

    if candidate_network is not None:  # same branch as sample_negatives_ilp.sample_negatives_ilp()
        pos_pairs_set = {tuple(p) for p in pos_pairs.tolist()}
        candidates = load_candidate_network(candidate_network, protein_to_idx, pos_pairs_set, cfg)
    else:
        candidates = build_candidate_set(len(protein_to_idx), pos_pairs, cfg, max_candidates=cfg.max_candidates,
                                         seed=cfg.seed, taxon_codes=taxonomy_codes, go_membership=go_membership,
                                         go_sizes=go_sizes)
    ctx = build_context(pos_pairs, protein_to_idx, idx_to_protein, candidates, 1.0, species_path=species_path,
                        go_annotations_path=go_path, taxonomy_codes=taxonomy_codes, n_taxa=n_taxa, go_bp=go_bp)

    # x = indicator of the saved negatives within the rebuilt pool; a negative outside it means the pool differs
    n = ctx.n_proteins
    neg_keys = np.array([min(protein_to_idx[r["protein1"]], protein_to_idx[r["protein2"]]) * n
                         + max(protein_to_idx[r["protein1"]], protein_to_idx[r["protein2"]]) for r in neg_rows])
    cand_keys = candidates[:, 0] * n + candidates[:, 1]
    order = np.argsort(cand_keys)
    loc = np.searchsorted(cand_keys, neg_keys, sorter=order)
    found = (loc < len(cand_keys)) & (cand_keys[order[np.minimum(loc, len(cand_keys) - 1)]] == neg_keys)
    if not found.all() or len(neg_rows) != ctx.n_neg:
        raise RuntimeError(f"{csv_path}: {np.sum(~found)} of {len(neg_rows)} negatives not in the rebuilt candidate "
                           f"pool (expected {ctx.n_neg} negatives); the pool does not match the original run.")
    x = np.zeros(len(candidates))
    x[order[loc]] = 1.0

    terms = {key: 0.0 for key in _TERM_KEY.values()}
    status = "reconstructed"
    if ctx.n_neg == len(candidates):
        status = "optimal (all candidates forced)"  # sampler skips the ILP; all terms reported as 0
    else:
        # Same deviation as each bias's build() objective, evaluated at x
        for b in assemble_active_biases(cfg):
            b.precompute(ctx)
            if not b.is_active():
                continue
            if b.name == "self":
                dev = abs(x[b.self_idx].sum() - b.target)
            elif b.name == "jaccard":
                dev = abs(b.J_cand @ x / ctx.n_neg - b.j_bar_pos)
            else:  # degree / taxon_pair: weighted absolute residual per protein / taxon pair
                weight = b.coef if b.name == "degree" else b.gamma
                dev = float(np.sum(weight * np.abs(b.M @ x - b.target)))
            terms[_TERM_KEY[b.name]] += b.lambda_weight * dev / b.U

    return {"split": split, "n_pos": ctx.n_pos, "n_neg": ctx.n_neg, "r": 1.0, "n_candidates": len(candidates),
            "obj_value": sum(terms.values()), **terms, "solver": "reconstructed", "wall_time_s": "",
            "mip_gap": "", "status": status}


def lambdas_from(row):
    return {f"lambda_{name}": float(row[f"neg_ilp_lambda_{name}"] or LAMBDA_DEFAULT)
            for name in ("degree", "taxon_pair", "self_loop", "jaccard")}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", required=True, help="Samplesheet the robustness run was started with")
    ap.add_argument("--results", required=True, type=Path, help="Robustness --outdir")
    ap.add_argument("--data-dir", type=Path, default=None,
                    help="Where the samplesheet's species/go_annotations files are, if not at the paths it gives "
                         "(default: the fetched <results>/<id>/data/ files)")
    ap.add_argument("--out", required=True, help="Output TSV")
    args = ap.parse_args()

    rows = []
    with open(args.samplesheet) as fh:
        for ss in csv.DictReader(fh):
            if ss["negative_sampling_method"] != "ilp":
                continue
            id_, run_dir = ss["id"], args.results / ss["id"]
            if ss.get("go_annotations"):
                species, go = (Path(ss[c]) if args.data_dir is None else args.data_dir / Path(ss[c]).name
                               for c in ("species", "go_annotations"))
            else:
                species, go = run_dir / "data" / "species.tsv", run_dir / "data" / "go_annotations.tsv"
            network = ss.get("candidate_network") or None
            if network and args.data_dir is not None:
                network = args.data_dir / Path(network).name
            for seed_dir in sorted(run_dir.glob("robustness/split*/sampling*")):
                split_seed = int(seed_dir.parent.name.removeprefix("split"))
                sampling_seed = int(seed_dir.name.removeprefix("sampling"))
                for split in SPLITS:
                    print(f"{id_} {split_seed}/{sampling_seed} {split} ...", file=sys.stderr, flush=True)
                    row = reconstruct(split, seed_dir / f"{split}.csv", lambdas_from(ss), sampling_seed, species, go,
                                      network, int(ss.get("neg_ilp_pool_factor") or 4))
                    rows.append({"Sample": f"{mqc_sample(id_, split)}_split{split_seed}_sampling{sampling_seed}",
                                 "ID": id_, "split_seed": split_seed, "sampling_seed": sampling_seed, **row})

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["Sample", "ID", "split_seed", "sampling_seed"] + DIAG_COLUMNS,
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
