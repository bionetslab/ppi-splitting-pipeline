#!/usr/bin/env python3
"""ILP-based bias-aware negative sampling for PPI splits.

Alternative to sample_negatives.py: solves a MILP that matches per-protein degree, per-taxon-pair interaction
counts, self-interaction count, and mean GO-BP Jaccard similarity between positive and negative sets.
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import cvxpy as cp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import mqc_sample, read_ppis  # noqa: E402

# Hard per-protein cap on negative degree (see build_problem): mx_p <= neg_ratio * (1 + MAX_DEGREE_SLACK) * d_plus_p
# so that individual proteins don't absorb all the --lambda-degree penalty
MAX_DEGREE_SLACK = 5


# ============================================================
# 1. Config & CLI
# ============================================================


@dataclass
class SamplingConfig:
    lambda_degree: float = 0.6
    lambda_taxon_pair: float = 0.0
    lambda_self_loop: float = 0.1
    lambda_jaccard: float = 0.3
    solver: str = "auto"  # "auto" | "gurobi" | "scip" | "highs"
    time_limit: float = 3600
    mip_gap: float = 0.01
    threads: int = 1
    seed: int = 42
    max_candidates: int = 50_000_000
    verbose: bool = False


def parse_args(argv=None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    ap.add_argument("--positives", required=True, help="Positive PPI CSV for this split")
    ap.add_argument("--output", required=True, help="Output labelled CSV for this split")
    ap.add_argument(
        "--split-name",
        default=None,
        help="Label for this split in diagnostics output " "(default: derived from --output filename)",
    )
    ap.add_argument(
        "--neg-ratio",
        type=float,
        default=1.0,
        help="|NEG| / |POS| for this split (default 1.0)",
    )

    # shared inputs
    ap.add_argument("--species", default=None)
    ap.add_argument("--go-annotations", default=None)
    ap.add_argument("--candidate-network", default=None)
    ap.add_argument("--gurobi-license", default=None)

    # weights
    ap.add_argument("--lambda-degree", type=float, default=0.0)
    ap.add_argument("--lambda-taxon-pair", type=float, default=0.0)
    ap.add_argument("--lambda-self-loop", type=float, default=0.0)
    ap.add_argument("--lambda-jaccard", type=float, default=0.0)

    # solver
    ap.add_argument("--solver", type=str.lower, choices=["auto", "gurobi", "scip", "highs"], default="auto")
    ap.add_argument("--time-limit", type=float, default=200)
    ap.add_argument("--mip-gap", type=float, default=0.01)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--max-candidates", type=int, default=50_000_000)

    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--diagnostics-out", default="neg_sampling_ilp_mqc.tsv")
    ap.add_argument(
        "--residuals-out",
        default="neg_sampling_ilp_residuals_mqc.tsv",
        help="Per-protein degree residual TSV, written only with --verbose",
    )
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--id", required=True, help="Dataset ID, for MultiQC tagging")

    return ap.parse_args(argv)


def _validate_config(cfg: SamplingConfig) -> None:
    for flag, val in [
        ("lambda-degree", cfg.lambda_degree),
        ("lambda-taxon-pair", cfg.lambda_taxon_pair),
        ("lambda-self-loop", cfg.lambda_self_loop),
        ("lambda-jaccard", cfg.lambda_jaccard),
    ]:
        if val < 0:
            raise ValueError(f"--{flag} must be >= 0 (got {val})")


def config_from_args(args: argparse.Namespace) -> SamplingConfig:
    """Build a SamplingConfig from CLI args (argparse already carries the defaults)."""
    cfg = SamplingConfig(
        lambda_degree=args.lambda_degree,
        lambda_taxon_pair=args.lambda_taxon_pair,
        lambda_self_loop=args.lambda_self_loop,
        lambda_jaccard=args.lambda_jaccard,
        solver=args.solver,
        time_limit=args.time_limit,
        mip_gap=args.mip_gap,
        threads=args.threads,
        seed=args.seed,
        max_candidates=args.max_candidates,
        verbose=args.verbose,
    )
    _validate_config(cfg)
    return cfg


# ============================================================
# 2. Data loading
# ============================================================


def build_protein_index(rows):
    """Return (protein_to_idx, idx_to_protein) covering every protein in `rows`."""
    proteins = sorted({p for r in rows for p in (r["protein1"], r["protein2"])})
    return {p: i for i, p in enumerate(proteins)}, proteins


def pos_pairs_from_rows(rows, protein_to_idx) -> np.ndarray:
    """Return (n_pos, 2) int64 array of (i, j) with i <= j."""
    if not rows:
        return np.zeros((0, 2), dtype=np.int64)
    pairs = [(protein_to_idx[r["protein1"]], protein_to_idx[r["protein2"]]) for r in rows]
    return np.array([(min(i, j), max(i, j)) for i, j in pairs], dtype=np.int64)


def load_species(path, protein_to_idx) -> np.ndarray:
    """Return an object array of taxon-id strings, one per protein index. Proteins absent from the file get "" """
    taxon_map = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            taxon_map[row["protein_id"].strip()] = row["taxon_id"].strip()
    taxonomy = np.empty(len(protein_to_idx), dtype=object)
    for p, idx in protein_to_idx.items():
        taxonomy[idx] = taxon_map.get(p, "")
    return taxonomy


def load_go_bp(path, protein_to_idx) -> list:
    """Return a list of frozensets of GO-BP term IDs, one per protein index.

    Reads the go_bp column of the go_annotations.tsv produced by fetch_data.py
    (columns: protein_id, go_bp, go_mf, go_cc; ';'-separated term lists).
    """
    go_map = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            terms = frozenset(t.strip() for t in row.get("go_bp", "").split(";") if t.strip())
            go_map[row["protein_id"].strip()] = terms
    result = [frozenset()] * len(protein_to_idx)
    for p, idx in protein_to_idx.items():
        result[idx] = go_map.get(p, frozenset())
    return result


def load_candidate_network(path, protein_to_idx, pos_pairs_set, cfg: SamplingConfig):
    """Read a pre-supplied candidate network CSV (protein1,protein2).

    Returns candidates restricted to the protein universe, excluding positives. Oversized files are subsampled
    via _subsample_candidate_pairs (same degree-weighted logic, and forced-in self-pairs, as the auto-generated
    pool)."""
    pairs = set()
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            p1, p2 = row["protein1"].strip(), row["protein2"].strip()
            if p1 not in protein_to_idx or p2 not in protein_to_idx:
                continue
            i, j = protein_to_idx[p1], protein_to_idx[p2]
            i, j = (min(i, j), max(i, j))
            if (i, j) in pos_pairs_set:
                continue
            pairs.add((i, j))

    network_pairs = np.array(sorted(pairs), dtype=np.int64) if pairs else np.zeros((0, 2), dtype=np.int64)

    if len(network_pairs) > cfg.max_candidates:
        logging.warning(
            "--candidate-network has %s pairs, exceeding --max-candidates=%s; "
            "subsampling %s pairs from it (weighted toward each protein's own "
            "negative-degree cap, not uniformly).",
            f"{len(network_pairs):,}",
            f"{cfg.max_candidates:,}",
            f"{cfg.max_candidates:,}",
        )

    n_proteins = len(protein_to_idx)
    pos_pairs_arr = (
        np.array(sorted(pos_pairs_set), dtype=np.int64) if pos_pairs_set else np.zeros((0, 2), dtype=np.int64)
    )
    return _subsample_candidate_pairs(
        n_proteins,
        pos_pairs_arr,
        cfg.max_candidates,
        cfg,
        seed=cfg.seed,
        given_pairs=network_pairs,
    )


# ============================================================
# 3. Candidate enumeration
# ============================================================


def build_candidate_set(
    n_proteins,
    pos_pairs,
    cfg: SamplingConfig,
    max_candidates=50_000_000,
    seed=42,
    taxon_codes=None,
    go_membership=None,
    go_sizes=None,
) -> np.ndarray:
    """Return (n_cand, 2) sorted (i<=j) upper-triangle array, excluding positives. Vectorized; delegates to
    _subsample_candidate_pairs instead of enumerating the full complement if that would exceed max_candidates.
    """
    print(f"Number of unique PPIs in the positive set: {len(pos_pairs)}")
    n_pairs_full = n_proteins * (n_proteins + 1) // 2
    n_est = n_pairs_full - len(pos_pairs)
    if n_est > max_candidates:
        logging.warning(
            "Default candidate set would have ~%s pairs, exceeding --max-candidates=%s; "
            "subsampling %s random candidate pairs instead of the full complement. "
            "Supply --candidate-network to restrict the pool deliberately, or raise "
            "--max-candidates if you have the memory.",
            f"{n_est:,}",
            f"{max_candidates:,}",
            f"{max_candidates:,}",
        )
        return _subsample_candidate_pairs(
            n_proteins,
            pos_pairs,
            max_candidates,
            cfg,
            seed=seed,
            taxon_codes=taxon_codes,
            go_membership=go_membership,
            go_sizes=go_sizes,
        )
    i_idx, j_idx = np.triu_indices(n_proteins)
    if len(pos_pairs):
        keys_all = i_idx.astype(np.int64) * n_proteins + j_idx.astype(np.int64)
        pos_keys = np.sort(pos_pairs[:, 0].astype(np.int64) * n_proteins + pos_pairs[:, 1].astype(np.int64))
        mask = ~np.isin(keys_all, pos_keys, assume_unique=True)
        i_idx, j_idx = i_idx[mask], j_idx[mask]
    return np.stack([i_idx, j_idx], axis=1).astype(np.int64)


def _missing_self_pairs(n_proteins, pos_pairs) -> np.ndarray:
    """(k, 2) array of every (i, i) that is not already a positive self-interaction."""
    all_i = np.arange(n_proteins, dtype=np.int64)
    if len(pos_pairs):
        pi, pj = pos_pairs[:, 0], pos_pairs[:, 1]
        pos_self = np.unique(pi[pi == pj])
        all_i = all_i[~np.isin(all_i, pos_self, assume_unique=True)]
    return np.stack([all_i, all_i], axis=1)


def _positive_same_species_fraction(pos_pairs, taxon_codes) -> float:
    """Fraction of positive pairs that are same-species (a self-interaction counts as same-species), edge-level
    (not averaged per protein)."""
    if len(pos_pairs) == 0:
        return 0.0
    i_arr, j_arr = pos_pairs[:, 0], pos_pairs[:, 1]
    return float(np.mean(taxon_codes[i_arr] == taxon_codes[j_arr]))


def _filter_same_species(keys, bi, bj, taxon_codes, want_same):
    """Restrict (keys, bi, bj) to same-species (want_same=True) or cross-species (want_same=False) pairs
    per taxon_codes; shared by _fill_stratum and _fill_stratum_from_pool."""
    same_mask = taxon_codes[bi] == taxon_codes[bj]
    keep_mask = same_mask if want_same else ~same_mask
    return keys[keep_mask], bi[keep_mask], bj[keep_mask]


def _fill_stratum_from_pool(
    pool_keys,
    n_proteins,
    exclude_keys,
    quota,
    rng,
    taxon_codes=None,
    want_same=None,
    go_membership=None,
    go_sizes=None,
    weights=None,
) -> np.ndarray:
    """Like _fill_stratum, but samples from an already-known `pool_keys` (e.g. --candidate-network) via direct weighted
    sampling, instead of _fill_stratum's generate-and-reject loop (too inefficient for a pool sparse relative to
    n_proteins**2)."""
    keys = pool_keys
    if len(exclude_keys):
        keys = keys[~np.isin(keys, exclude_keys, assume_unique=True)]
    if len(keys) == 0:
        return np.empty(0, dtype=np.int64)

    bi, bj = keys // n_proteins, keys % n_proteins
    if want_same is not None:
        keys, bi, bj = _filter_same_species(keys, bi, bj, taxon_codes, want_same)
    if len(keys) == 0:
        return np.empty(0, dtype=np.int64)

    def _weighted_sample(idx, n):
        """idx: positions into keys/bi/bj to draw up to n (<= len(idx)) from."""
        if len(idx) <= n:
            return idx
        if weights is None:
            return rng.choice(idx, size=n, replace=False)
        w = weights[bi[idx]] + weights[bj[idx]]
        return rng.choice(idx, size=n, replace=False, p=w / w.sum())

    all_idx = np.arange(len(keys))
    if go_membership is not None:
        jac = _pairwise_jaccard(np.stack([bi, bj], axis=1), go_membership, go_sizes)
        priority_idx, filler_idx = all_idx[jac > 0], all_idx[jac == 0]
    else:
        priority_idx, filler_idx = np.empty(0, dtype=np.int64), all_idx

    if len(priority_idx) >= quota:
        chosen = _weighted_sample(priority_idx, quota)
    else:
        chosen = np.concatenate([priority_idx, _weighted_sample(filler_idx, quota - len(priority_idx))])

    return np.sort(keys[chosen])


def _fill_stratum(
    rng,
    n_proteins,
    exclude_keys,
    quota,
    taxon_codes=None,
    want_same=None,
    go_membership=None,
    go_sizes=None,
    weights=None,
    pool_keys=None,
    max_rounds=200,
) -> np.ndarray:
    """Draw up to `quota` unique (i,j) keys (encoded i*n_proteins+j) via repeated random batches, excluding
    `exclude_keys` (e.g. positives).

    :param rng: random number generator used to sample from the full n_proteins**2 space.
    :param n_proteins: number of proteins in the dataset
    :param exclude_keys: (optional) keys to exclude from the sampling process (e.g. positives)
    :param quota: number of unique (i,j) pairs to draw.
    :param taxon_codes: taxon code for each protein, used to stratify the same-species/cross-species
    :param want_same: restrict to same-/cross-species pairs (stratifies the rare same-species population against the
    much larger cross-species one).
    :param go_membership: and
    :param go_sizes: keep nonzero GO-BP Jaccard pairs first (rare under a uniform draw).
    :param weights: sample proteins proportional to their negative-degree cap instead of uniformly, so low-cap
    proteins aren't wasted past their cap.
    :param pool_keys: sample from a fixed known pool (--candidate-network) via _fill_stratum_from_pool
    :param max_rounds: maximum number of sampling rounds to try before giving up. 200 is a good default, but

    :returns: (quota,) array of unique (i,j) keys (encoded i*n_proteins+j) drawn from the full n_proteins**2 space.
    """
    if quota <= 0:
        return np.empty(0, dtype=np.int64)

    if pool_keys is not None:
        return _fill_stratum_from_pool(
            pool_keys,
            n_proteins,
            exclude_keys,
            quota,
            rng,
            taxon_codes=taxon_codes,
            want_same=want_same,
            go_membership=go_membership,
            go_sizes=go_sizes,
            weights=weights,
        )

    priority_keys = np.empty(0, dtype=np.int64)
    filler_keys = np.empty(0, dtype=np.int64)

    for _ in range(max_rounds):
        have = len(priority_keys) + len(filler_keys)
        if have >= quota:
            break
        batch_size = max(10_000, (quota - have) * 5)
        if weights is None:
            i = rng.integers(0, n_proteins, size=batch_size, dtype=np.int64)
            j = rng.integers(0, n_proteins, size=batch_size, dtype=np.int64)
        else:
            i = rng.choice(n_proteins, size=batch_size, p=weights)
            j = rng.choice(n_proteins, size=batch_size, p=weights)
        lo, hi = np.minimum(i, j), np.maximum(i, j)
        batch_keys = np.unique(lo * n_proteins + hi)
        if len(exclude_keys):
            batch_keys = batch_keys[~np.isin(batch_keys, exclude_keys, assume_unique=True)]
        already = np.union1d(priority_keys, filler_keys)
        if len(already):
            batch_keys = batch_keys[~np.isin(batch_keys, already, assume_unique=True)]
        if len(batch_keys) == 0:
            continue

        bi, bj = batch_keys // n_proteins, batch_keys % n_proteins
        if want_same is not None:
            batch_keys, bi, bj = _filter_same_species(batch_keys, bi, bj, taxon_codes, want_same)
        if len(batch_keys) == 0:
            continue

        if go_membership is not None:
            jac = _pairwise_jaccard(np.stack([bi, bj], axis=1), go_membership, go_sizes)
            priority_keys = np.union1d(priority_keys, batch_keys[jac > 0])
            filler_keys = np.union1d(filler_keys, batch_keys[jac == 0])
        else:
            filler_keys = np.union1d(filler_keys, batch_keys)
    # Not a for-else: that also fired when the last round filled the quota.
    if len(priority_keys) + len(filler_keys) < quota:
        logging.warning(
            "Could not fill a candidate stratum (quota=%d) after %d sampling "
            "rounds; its true population is likely smaller than requested. "
            "Using the %d candidates found.",
            quota,
            max_rounds,
            len(priority_keys) + len(filler_keys),
        )

    if len(priority_keys) >= quota:
        return priority_keys[:quota]
    return np.union1d(priority_keys, filler_keys[: quota - len(priority_keys)])


def _degree_weights(pos_pairs, n_proteins) -> np.ndarray:
    """Per-protein sampling probability, proportional to positive degree (i.e. to its negative-degree cap,
    a fixed multiple of that degree). Every protein has degree > 0, so no zero-weight case to guard against."""
    d_plus = _degree_array(pos_pairs, n_proteins)
    return d_plus / d_plus.sum()


def _subsample_candidate_pairs(
    n_proteins,
    pos_pairs,
    n_target,
    cfg: SamplingConfig,
    seed=42,
    taxon_codes=None,
    go_membership=None,
    go_sizes=None,
    given_pairs=None,
) -> np.ndarray:
    """Randomly draw n_target unique (i,j) pairs (i<=j), excluding positives, without materializing the full
    upper-triangle. Used when the full complement is too large to enumerate.

    given_pairs: sample from this fixed pool (--candidate-network) instead of the whole n_proteins**2 space
    * Self-pairs are force-kept if --lambda-self-loop > 0.
    * The same-/cross-species budget is stratified to the positive set's ratio if taxon_codes is given,
    * Nonzero-Jaccard pairs are prioritized if go_membership/go_sizes is given.
    * Throughout, proteins are drawn proportional to their negative-degree cap (_degree_weights), not uniformly,
    so low-cap proteins aren't wasted on unusable candidates."""
    rng = np.random.default_rng(seed)
    weights = _degree_weights(pos_pairs, n_proteins)
    pos_keys = (
        np.sort(pos_pairs[:, 0].astype(np.int64) * n_proteins + pos_pairs[:, 1].astype(np.int64))
        if len(pos_pairs)
        else np.empty(0, dtype=np.int64)
    )

    pool_keys = None
    if given_pairs is not None:
        pool_keys = (
            np.sort(given_pairs[:, 0].astype(np.int64) * n_proteins + given_pairs[:, 1].astype(np.int64))
            if len(given_pairs)
            else np.empty(0, dtype=np.int64)
        )

    self_keys = np.empty(0, dtype=np.int64)
    if cfg.lambda_self_loop > 0:
        self_pairs = _missing_self_pairs(n_proteins, pos_pairs)
        self_keys = np.sort(self_pairs[:, 0] * n_proteins + self_pairs[:, 1])
        if pool_keys is not None:
            # ilp_candidates: only self-pairs the supplied network contains.
            self_keys = np.intersect1d(self_keys, pool_keys)
    remaining = max(n_target - len(self_keys), 0)

    if taxon_codes is not None:
        same_ratio = _positive_same_species_fraction(pos_pairs, taxon_codes)
        n_same = int(round(same_ratio * remaining))
        n_cross = remaining - n_same
        same_keys = _fill_stratum(
            rng,
            n_proteins,
            pos_keys,
            n_same,
            taxon_codes=taxon_codes,
            want_same=True,
            go_membership=go_membership,
            go_sizes=go_sizes,
            weights=weights,
            pool_keys=pool_keys,
        )
        cross_keys = _fill_stratum(
            rng,
            n_proteins,
            pos_keys,
            n_cross,
            taxon_codes=taxon_codes,
            want_same=False,
            go_membership=go_membership,
            go_sizes=go_sizes,
            weights=weights,
            pool_keys=pool_keys,
        )
        drawn_keys = np.union1d(same_keys, cross_keys)
        shortfall = remaining - len(drawn_keys)
        if shortfall > 0:
            exclude = np.union1d(pos_keys, np.union1d(self_keys, drawn_keys))
            drawn_keys = np.union1d(
                drawn_keys,
                _fill_stratum(
                    rng,
                    n_proteins,
                    exclude,
                    shortfall,
                    weights=weights,
                    pool_keys=pool_keys,
                ),
            )
    else:
        drawn_keys = _fill_stratum(
            rng,
            n_proteins,
            pos_keys,
            remaining,
            go_membership=go_membership,
            go_sizes=go_sizes,
            weights=weights,
            pool_keys=pool_keys,
        )

    keys = np.union1d(self_keys, drawn_keys)

    if len(keys) > n_target:
        is_self = np.isin(keys, self_keys, assume_unique=True)
        non_self_keys = keys[~is_self]
        n_non_self_target = max(n_target - int(is_self.sum()), 0)
        keys = np.sort(np.concatenate([keys[is_self], non_self_keys[:n_non_self_target]]))

    i_idx, j_idx = keys // n_proteins, keys % n_proteins
    return np.stack([i_idx, j_idx], axis=1).astype(np.int64)


def _build_incidence(n_proteins, candidates) -> sp.csr_matrix:
    """(n_proteins, n_cand) 0/1 matrix; self-loops contribute 1 (not 2)."""
    i_arr, j_arr = candidates[:, 0], candidates[:, 1]
    self_mask = i_arr == j_arr
    n_cand = len(candidates)
    rows = np.concatenate([i_arr, j_arr[~self_mask]])
    cols = np.concatenate([np.arange(n_cand), np.arange(n_cand)[~self_mask]])
    data = np.ones(len(rows), dtype=np.float64)
    return sp.csr_matrix((data, (rows, cols)), shape=(n_proteins, n_cand))


def _pairwise_jaccard(pairs, membership, sizes) -> np.ndarray:
    i_arr, j_arr = pairs[:, 0], pairs[:, 1]
    inter = np.asarray(membership[i_arr].multiply(membership[j_arr]).sum(axis=1)).ravel()
    union = sizes[i_arr] + sizes[j_arr] - inter
    jac = np.zeros(len(i_arr), dtype=np.float64)
    nz = union > 0
    jac[nz] = inter[nz] / union[nz]
    return jac


def _build_go_membership(go_bp) -> tuple[sp.csr_matrix, np.ndarray]:
    """(n_proteins, n_terms) 0/1 membership matrix and per-protein term counts."""
    terms = sorted({t for s in go_bp for t in s})
    term_to_col = {t: k for k, t in enumerate(terms)}
    rows, cols = [], []
    for p, s in enumerate(go_bp):
        for t in s:
            rows.append(p)
            cols.append(term_to_col[t])
    membership = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(go_bp), len(terms)))
    sizes = np.asarray(membership.sum(axis=1)).ravel()
    return membership, sizes


# ============================================================
# 3b. Descriptive dataset stats (pos vs. neg)
# ============================================================


def _unique_proteins(pairs) -> np.ndarray:
    if len(pairs) == 0:
        return np.zeros(0, dtype=np.int64)
    return np.unique(pairs.ravel())


def _degree_array(pairs, n_proteins) -> np.ndarray:
    """Per-protein interaction degree; self-loops contribute 1 (not 2)"""
    deg = np.zeros(n_proteins, dtype=np.float64)
    if len(pairs):
        i_arr, j_arr = pairs[:, 0], pairs[:, 1]
        self_mask = i_arr == j_arr
        np.add.at(deg, i_arr, 1.0)
        np.add.at(deg, j_arr[~self_mask], 1.0)
    return deg


def _same_species_ratio(pairs, n_proteins, taxon_codes) -> np.ndarray:
    """Per-protein fraction of interactions that are same-species. NaN for proteins with none."""
    same = np.zeros(n_proteins, dtype=np.float64)
    total = np.zeros(n_proteins, dtype=np.float64)
    if len(pairs):
        i_arr, j_arr = pairs[:, 0], pairs[:, 1]
        self_mask = i_arr == j_arr
        is_same = (taxon_codes[i_arr] == taxon_codes[j_arr]).astype(np.float64)
        rows = np.concatenate([i_arr, j_arr[~self_mask]])
        vals = np.concatenate([is_same, is_same[~self_mask]])
        np.add.at(total, rows, 1.0)
        np.add.at(same, rows, vals)
    ratio = np.full(n_proteins, np.nan)
    nz = total > 0
    ratio[nz] = same[nz] / total[nz]
    return ratio


def _fmt_degree_stats(deg: np.ndarray) -> str:
    deg = deg[deg > 0]
    if len(deg) == 0:
        return "n/a"
    return f"median={np.median(deg):.2f} mean={np.mean(deg):.2f} max={np.max(deg):.0f}"


def print_objective_breakdown(name, diag: dict) -> None:
    """Print each term's share of the objective, so it's obvious when one term (e.g. an unsatisfiable self-loop
    target) is silently dominating and crowding out the others, even when their lambdas are set equal."""
    obj = diag["obj_value"]
    terms = [
        ("degree", diag["bias_deg_term"]),
        ("taxon_pair", diag["bias_tax_term"]),
        ("self_loop", diag["bias_self_term"]),
        ("jaccard", diag["bias_jac_term"]),
    ]
    print(f"=== Objective breakdown: {name} (total={obj:.4f}) ===")
    for label, val in terms:
        pct = (val / obj * 100) if obj else 0.0
        print(f"  {label:<10} {val:.4f}  ({pct:5.1f}%)")
    print()


def print_dataset_stats(name, pos_pairs, neg_pairs, ctx: "BuildContext", cfg: SamplingConfig) -> None:
    """Print pos-vs-neg stats for one split. Protein counts always shown; the rest gated on which --lambda-* biases
    are active, since that data is only loaded when needed."""
    print(f"\n=== Dataset stats: {name} ===")

    print(f"PPIs                  -- positive: {len(pos_pairs)}   negative: {len(neg_pairs)}")

    pos_proteins = _unique_proteins(pos_pairs)
    neg_proteins = _unique_proteins(neg_pairs)
    print(f"Unique proteins       -- positive: {len(pos_proteins)}   negative: {len(neg_proteins)}")

    if cfg.lambda_self_loop > 0:
        n_self_pos = int(np.sum(pos_pairs[:, 0] == pos_pairs[:, 1])) if len(pos_pairs) else 0
        n_self_neg = int(np.sum(neg_pairs[:, 0] == neg_pairs[:, 1])) if len(neg_pairs) else 0
        print(f"Self-interactions     -- positive: {n_self_pos}   negative: {n_self_neg}")

    if cfg.lambda_degree > 0:
        deg_pos = _degree_array(pos_pairs, ctx.n_proteins)
        deg_neg = _degree_array(neg_pairs, ctx.n_proteins)
        print(f"Degree (positive)     -- {_fmt_degree_stats(deg_pos)}")
        print(f"Degree (negative)     -- {_fmt_degree_stats(deg_neg)}")

    species_used = cfg.lambda_taxon_pair > 0
    if species_used and ctx.species_path is not None:
        taxon_codes, _ = ctx.ensure_taxonomy()
        same_pos = _same_species_ratio(pos_pairs, ctx.n_proteins, taxon_codes)
        same_neg = _same_species_ratio(neg_pairs, ctx.n_proteins, taxon_codes)
        mean_same_pos, mean_same_neg = np.nanmean(same_pos), np.nanmean(same_neg)
        print(f"Same-species ratio    -- positive: {mean_same_pos:.3f}   negative: {mean_same_neg:.3f}")
        print(f"Cross-species ratio   -- positive: {1 - mean_same_pos:.3f}   negative: {1 - mean_same_neg:.3f}")

    if cfg.lambda_jaccard > 0 and ctx.go_annotations_path is not None:
        membership, sizes = _build_go_membership(ctx.ensure_go_bp())
        jac_pos = _pairwise_jaccard(pos_pairs, membership, sizes) if len(pos_pairs) else np.zeros(0)
        jac_neg = _pairwise_jaccard(neg_pairs, membership, sizes) if len(neg_pairs) else np.zeros(0)
        mean_jac_pos = np.mean(jac_pos) if len(jac_pos) else float("nan")
        mean_jac_neg = np.mean(jac_neg) if len(jac_neg) else float("nan")
        print(f"Mean GO-BP Jaccard    -- positive: {mean_jac_pos:.3f}   negative: {mean_jac_neg:.3f}")
    print()


# ============================================================
# 4. BuildContext + BiasTerm interface
# ============================================================


@dataclass
class BuildContext:
    """Everything the bias terms may need. Expensive derived fields (taxonomy codes, GO membership) are populated
    lazily via the ensure_* methods, only when a bias term actually requests them."""

    n_proteins: int
    candidates: np.ndarray
    pos_pairs: np.ndarray
    n_pos: int
    n_neg: int
    r: float
    incidence: sp.csr_matrix
    protein_to_idx: dict
    idx_to_protein: list
    species_path: object = None
    go_annotations_path: object = None

    taxonomy: np.ndarray | None = field(default=None, init=False, repr=False)
    taxonomy_codes: np.ndarray | None = field(default=None, init=False, repr=False)
    n_taxa: int | None = field(default=None, init=False, repr=False)
    go_bp: list | None = field(default=None, init=False, repr=False)

    def ensure_taxonomy(self):
        if self.taxonomy_codes is None:
            if self.species_path is None:
                raise ValueError("--species is required for this bias term")
            self.taxonomy = load_species(self.species_path, self.protein_to_idx)
            uniq, codes = np.unique(self.taxonomy, return_inverse=True)
            self.taxonomy_codes = codes.astype(np.int64)
            self.n_taxa = int(len(uniq))
        return self.taxonomy_codes, self.n_taxa

    def ensure_go_bp(self):
        if self.go_bp is None:
            if self.go_annotations_path is None:
                raise ValueError("--go-annotations is required when --lambda-jaccard > 0")
            self.go_bp = load_go_bp(self.go_annotations_path, self.protein_to_idx)
        return self.go_bp


def build_context(
    pos_pairs,
    protein_to_idx,
    idx_to_protein,
    candidates,
    neg_ratio,
    species_path=None,
    go_annotations_path=None,
    taxonomy_codes=None,
    n_taxa=None,
    go_bp=None,
) -> BuildContext:
    """taxonomy_codes/n_taxa/go_bp let a caller that already loaded them pass them straight through, instead of
    ensure_taxonomy()/ensure_go_bp() re-reading the same files."""
    n_proteins = len(protein_to_idx)
    n_pos = len(pos_pairs)
    n_neg = int(round(neg_ratio * n_pos))
    incidence = _build_incidence(n_proteins, candidates)
    ctx = BuildContext(
        n_proteins=n_proteins,
        candidates=candidates,
        pos_pairs=pos_pairs,
        n_pos=n_pos,
        n_neg=n_neg,
        r=neg_ratio,
        incidence=incidence,
        protein_to_idx=protein_to_idx,
        idx_to_protein=idx_to_protein,
        species_path=species_path,
        go_annotations_path=go_annotations_path,
    )
    if taxonomy_codes is not None:
        ctx.taxonomy_codes = taxonomy_codes
        ctx.n_taxa = n_taxa
    if go_bp is not None:
        ctx.go_bp = go_bp
    return ctx


class BiasTerm:
    name = "base"

    def __init__(self, lambda_weight: float):
        self.lambda_weight = float(lambda_weight)
        self._active = False

    def is_active(self) -> bool:
        return self._active and self.lambda_weight > 0

    def precompute(self, ctx: BuildContext) -> None:
        raise NotImplementedError

    def build(self, x: cp.Variable, ctx: BuildContext):
        """Return (aux_vars, constraints, objective_expr) already scaled by lambda_weight / U."""
        raise NotImplementedError

    def debug_rows(self, x_value, ctx: BuildContext) -> list:
        return []


class SelfLoopBias(BiasTerm):
    name = "self"

    def precompute(self, ctx: BuildContext) -> None:
        i_arr, j_arr = ctx.candidates[:, 0], ctx.candidates[:, 1]
        self.self_idx = np.flatnonzero(i_arr == j_arr)
        d_size = len(self.self_idx)
        if len(ctx.pos_pairs):
            s_plus = int(np.sum(ctx.pos_pairs[:, 0] == ctx.pos_pairs[:, 1]))
        else:
            s_plus = 0
        self.target = ctx.r * s_plus
        self.U = max(self.target, d_size - self.target)
        self._active = self.U > 0

    def build(self, x, ctx):
        tau = cp.Variable(nonneg=True)
        d_sum = cp.sum(x[self.self_idx]) if len(self.self_idx) else 0.0
        constraints = [tau >= d_sum - self.target, tau >= self.target - d_sum]
        obj = self.lambda_weight * tau / self.U
        return [tau], constraints, obj


class JaccardMeanBias(BiasTerm):
    name = "jaccard"

    def precompute(self, ctx: BuildContext) -> None:
        go_bp = ctx.ensure_go_bp()
        membership, sizes = _build_go_membership(go_bp)

        self.J_cand = _pairwise_jaccard(ctx.candidates, membership, sizes)
        if len(ctx.pos_pairs):
            j_pos = _pairwise_jaccard(ctx.pos_pairs, membership, sizes)
            self.j_bar_pos = float(np.mean(j_pos))
        else:
            self.j_bar_pos = 0.0

        self.U = max(self.j_bar_pos, 1.0 - self.j_bar_pos)
        self._active = self.U > 0

    def build(self, x, ctx):
        z = cp.Variable(nonneg=True)
        coef = self.J_cand / ctx.n_neg
        term = coef @ x
        constraints = [z >= term - self.j_bar_pos, z >= self.j_bar_pos - term]
        obj = self.lambda_weight * z / self.U
        return [z], constraints, obj


class DegreeBias(BiasTerm):
    """Per-protein aggregate degree (no taxon)."""

    name = "degree"

    def precompute(self, ctx: BuildContext) -> None:
        n = ctx.n_proteins
        dplus = np.zeros(n, dtype=np.float64)
        pos = ctx.pos_pairs
        if len(pos):
            pi, pj = pos[:, 0], pos[:, 1]
            pself = pi == pj
            np.add.at(dplus, pi, 1.0)
            np.add.at(dplus, pj[~pself], 1.0)
        n_cand_per_p = np.asarray(ctx.incidence.sum(axis=1)).ravel()
        active_mask = (dplus > 0) | (n_cand_per_p > 0)
        if not np.any(active_mask):
            self._active = False
            return

        coef = np.zeros(n, dtype=np.float64)
        coef[active_mask] = 1.0 / np.log1p(dplus[active_mask])

        active_idx = np.flatnonzero(active_mask)
        target = ctx.r * dplus[active_idx]
        n_cand_active = n_cand_per_p[active_idx]
        coef_active = coef[active_idx]
        U = float(np.sum(coef_active * np.maximum(target, n_cand_active - target)))

        self.active_idx = active_idx
        self.M = ctx.incidence[active_idx, :]
        self.coef = coef_active
        self._dplus = dplus[active_idx]
        self.n_cand_per_p = n_cand_active
        self.target = target
        self.U = U
        self.n_groups = len(active_idx)
        self._active = U > 0

    def build(self, x, ctx):
        u = cp.Variable(self.n_groups, nonneg=True)
        mx = self.M @ x
        constraints = [u >= mx - self.target, u >= self.target - mx]
        obj = self.lambda_weight * cp.sum(cp.multiply(self.coef, u)) / self.U
        return [u], constraints, obj

    def debug_rows(self, x_value, ctx):
        mx = np.asarray(self.M @ np.round(np.asarray(x_value))).ravel()
        rows = []
        for k, p in enumerate(self.active_idx):
            rows.append(
                {
                    "protein_id": ctx.idx_to_protein[int(p)],
                    "taxon": "",
                    "d_plus": float(self._dplus[k]),
                    "d_minus": float(mx[k]),
                    "residual": float(mx[k] - self.target[k]),
                }
            )
        return rows


class TaxonPairBias(BiasTerm):
    """Global taxon-pair counts."""

    name = "taxon_pair"

    def precompute(self, ctx: BuildContext) -> None:
        taxon, T = ctx.ensure_taxonomy()
        cand = ctx.candidates
        ti, tj = taxon[cand[:, 0]], taxon[cand[:, 1]]
        t_lo, t_hi = np.minimum(ti, tj), np.maximum(ti, tj)
        key_c = t_lo.astype(np.int64) * T + t_hi.astype(np.int64)

        pos = ctx.pos_pairs
        if len(pos):
            pti, ptj = taxon[pos[:, 0]], taxon[pos[:, 1]]
            pt_lo, pt_hi = np.minimum(pti, ptj), np.maximum(pti, ptj)
            key_p = pt_lo.astype(np.int64) * T + pt_hi.astype(np.int64)
        else:
            key_p = np.zeros(0, dtype=np.int64)

        uniq_key_c = np.unique(key_c)
        if len(key_p):
            uniq_key_p, pos_counts = np.unique(key_p, return_counts=True)
        else:
            uniq_key_p, pos_counts = np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
        active_keys = np.union1d(uniq_key_c, uniq_key_p)
        n_groups = len(active_keys)
        if n_groups == 0:
            self._active = False
            return

        group_c = np.searchsorted(active_keys, key_c)
        group_p_idx = np.searchsorted(active_keys, uniq_key_p)

        m_plus = np.zeros(n_groups, dtype=np.float64)
        m_plus[group_p_idx] = pos_counts
        n_cand_per_group = np.bincount(group_c, minlength=n_groups).astype(np.float64)
        with np.errstate(divide="ignore"):
            gamma = np.where(m_plus > 0, 1.0 / np.log1p(m_plus), 1.0 / math.log(2.0))
        target = ctx.r * m_plus
        U = float(np.sum(gamma * np.maximum(target, n_cand_per_group - target)))

        self.M = sp.csr_matrix(
            (np.ones(len(group_c)), (group_c, np.arange(len(cand)))),
            shape=(n_groups, len(cand)),
        )
        self.gamma = gamma
        self.target = target
        self.U = U
        self.n_groups = n_groups
        self._active = U > 0

    def build(self, x, ctx):
        mu = cp.Variable(self.n_groups, nonneg=True)
        mx = self.M @ x
        constraints = [mu >= mx - self.target, mu >= self.target - mx]
        obj = self.lambda_weight * cp.sum(cp.multiply(self.gamma, mu)) / self.U
        return [mu], constraints, obj


# ============================================================
# 5. Model assembly and solve
# ============================================================


def assemble_active_biases(cfg: SamplingConfig):
    """Return [bias terms with lambda > 0]. Whether each ends up active (nonzero U, data available) is decided
    later by precompute()/is_active()."""
    biases = []
    if cfg.lambda_degree > 0:
        biases.append(DegreeBias(cfg.lambda_degree))
    if cfg.lambda_taxon_pair > 0:
        biases.append(TaxonPairBias(cfg.lambda_taxon_pair))
    if cfg.lambda_self_loop > 0:
        biases.append(SelfLoopBias(cfg.lambda_self_loop))
    if cfg.lambda_jaccard > 0:
        biases.append(JaccardMeanBias(cfg.lambda_jaccard))
    return biases


_TERM_KEY = {
    "self": "bias_self_term",
    "jaccard": "bias_jac_term",
    "degree": "bias_deg_term",
    "taxon_pair": "bias_tax_term",
}


def _max_degree_cap(ctx: BuildContext) -> np.ndarray:
    """Per-protein hard cap on total negative degree: neg_ratio * (1 + MAX_DEGREE_SLACK) * d_plus. See
    MAX_DEGREE_SLACK at the top of the file."""
    dplus = _degree_array(ctx.pos_pairs, ctx.n_proteins)
    return ctx.r * (1.0 + MAX_DEGREE_SLACK) * dplus


def build_problem(ctx: BuildContext, active_biases):
    x = cp.Variable(len(ctx.candidates), boolean=True)
    constraints = [cp.sum(x) == ctx.n_neg, ctx.incidence @ x <= _max_degree_cap(ctx)]

    objective_terms = []
    term_exprs = []  # list of (bias, scaled_expr)
    for b in active_biases:
        _, cons, raw_expr = b.build(x, ctx)
        constraints += cons
        term_exprs.append((b, raw_expr))
        objective_terms.append(raw_expr)

    objective = cp.Minimize(sum(objective_terms))
    problem = cp.Problem(objective, constraints)
    return problem, x, term_exprs


# ============================================================
# 5b. Solver selection
# ============================================================


def _solver_options(solver_name, cfg: SamplingConfig) -> dict:
    if solver_name == cp.GUROBI:
        return {
            "TimeLimit": cfg.time_limit,
            "MIPGap": cfg.mip_gap,
            "Threads": cfg.threads,
            "Seed": cfg.seed,
        }
    if solver_name == cp.HIGHS:
        return {
            "time_limit": cfg.time_limit,
            "mip_rel_gap": cfg.mip_gap,
            "threads": cfg.threads,
            "random_seed": cfg.seed,
        }
    if solver_name == cp.SCIP:
        return {
            "scip_params": {
                "limits/time": cfg.time_limit,
                "limits/gap": cfg.mip_gap,
                "randomization/randomseedshift": cfg.seed,
            }
        }
    return {}


def select_solver(cfg: SamplingConfig, gurobi_license, verbose: bool):
    if gurobi_license:
        os.environ["GRB_LICENSE_FILE"] = str(Path(gurobi_license).resolve())

    prefer = cfg.solver.lower()
    if prefer == "gurobi":
        return cp.GUROBI, _solver_options(cp.GUROBI, cfg)
    if prefer == "scip":
        return cp.SCIP, _solver_options(cp.SCIP, cfg)
    if prefer == "highs":
        return cp.HIGHS, _solver_options(cp.HIGHS, cfg)
    if prefer != "auto":
        raise ValueError(f"Unknown --solver {cfg.solver!r}; expected auto/gurobi/scip/highs")

    try:
        import gurobipy

        gurobipy.Model()  # triggers a license check
        return cp.GUROBI, _solver_options(cp.GUROBI, cfg)
    except Exception as exc:
        if verbose:
            logging.info(
                "Gurobi unavailable (%s); falling back to an open-source solver.",
                type(exc).__name__,
            )

    installed = cp.installed_solvers()
    for cand in (cp.SCIP, cp.HIGHS, cp.CBC, cp.GLPK_MI):
        if cand in installed:
            return cand, _solver_options(cand, cfg)
    raise RuntimeError("No MIP solver installed (need one of GUROBI, SCIP, HIGHS, CBC, GLPK_MI).")


def solve(problem: cp.Problem, solver, options: dict, verbose: bool = False) -> dict:
    t0 = time.time()
    problem.solve(solver=solver, verbose=verbose, **options)
    wall = time.time() - t0
    if problem.status not in cp.settings.SOLUTION_PRESENT:
        raise RuntimeError(f"Solver {solver} failed to find a solution (status={problem.status})")
    return {"status": problem.status, "wall_time_s": wall, "obj_value": problem.value}


# ============================================================
# 6. Output
# ============================================================


def extract_negatives(x_value, ctx: BuildContext) -> np.ndarray:
    x_val = np.asarray(x_value).ravel()
    if np.any(np.abs(np.round(x_val) - x_val) > 1e-4):
        raise RuntimeError("Solver returned a non-integral solution for x.")
    x_round = np.round(x_val).astype(np.int64)
    n_selected = int(x_round.sum())
    if n_selected != ctx.n_neg:
        raise RuntimeError(f"Rounded selection has {n_selected} pairs, expected |NEG|={ctx.n_neg}.")
    return ctx.candidates[x_round == 1]


def write_split_csv(pos_rows, negative_pairs, idx_to_protein, out_path) -> None:
    extra_fields = [k for k in (pos_rows[0].keys() if pos_rows else []) if k not in ("protein1", "protein2")]
    fieldnames = ["protein1", "protein2", "label"] + extra_fields
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in pos_rows:
            out = {f: row.get(f, "") for f in extra_fields}
            out["protein1"] = row["protein1"]
            out["protein2"] = row["protein2"]
            out["label"] = 1
            writer.writerow(out)
        for i, j in negative_pairs:
            out = {f: "" for f in extra_fields}
            out["protein1"] = idx_to_protein[i]
            out["protein2"] = idx_to_protein[j]
            out["label"] = 0
            writer.writerow(out)


DIAG_COLUMNS = [
    "split",
    "n_pos",
    "n_neg",
    "r",
    "n_candidates",
    "obj_value",
    "bias_deg_term",
    "bias_tax_term",
    "bias_self_term",
    "bias_jac_term",
    "solver",
    "wall_time_s",
    "mip_gap",
    "status",
]


def write_diagnostics(rows, out_path, id_) -> None:
    """Sample/ID qualify each row so this merges cleanly across datasets in the combined MultiQC report,
    like Classifier Performance."""
    fieldnames = ["Sample", "ID"] + DIAG_COLUMNS
    with open(out_path, "w", newline="") as fh:
        fh.write(
            "# id: 'neg_ilp_diagnostics'\n"
            "# section_name: 'ILP Negative Sampler Diagnostics'\n"
            "# plot_type: 'table'\n"
            "# pconfig:\n"
            "#     id: 'neg_ilp_diagnostics_table'\n"
            "#     title: 'ILP Negative Sampler Diagnostics'\n"
        )
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            out = {k: row.get(k, "") for k in DIAG_COLUMNS}
            out["Sample"] = mqc_sample(id_, row.get("split", ""))
            out["ID"] = id_
            writer.writerow(out)


def write_neg_generalstats(diag, out_path, id_) -> None:
    """Writes to the same shared 'neg_generalstats' id as sample_negatives.py, so ILP-sampled datasets appear in the
    combined General Statistics table too."""
    sample = mqc_sample(id_, diag.get("split", ""))
    with open(out_path, "w") as fh:
        fh.write(
            "# id: 'neg_generalstats'\n"
            "# plot_type: 'generalstats'\n"
            "# pconfig:\n"
            "#     - ID:\n"
            "#         title: 'ID'\n"
            "#         description: 'Dataset ID'\n"
            "#     - n_positives:\n"
            "#         title: 'Positives'\n"
            "#         description: 'Positive PPIs in the split'\n"
            "#         format: '{:,.0f}'\n"
            "#         scale: 'Blues'\n"
            "#     - n_negatives:\n"
            "#         title: 'Negatives'\n"
            "#         description: 'Sampled negatives for the split'\n"
            "#         format: '{:,.0f}'\n"
            "#         scale: 'Oranges'\n"
            "Sample\tID\tn_positives\tn_negatives\n"
            f"{sample}\t{id_}\t{diag.get('n_pos', '')}\t{diag.get('n_neg', '')}\n"
        )


RESIDUAL_COLUMNS = ["split", "protein_id", "taxon", "d_plus", "d_minus", "residual"]


def write_residuals(rows, out_path, id_) -> None:
    """Same Sample/ID qualification as write_diagnostics; residuals are per-protein so Sample also includes
    protein_id to stay unique."""
    fieldnames = ["Sample", "ID"] + RESIDUAL_COLUMNS
    with open(out_path, "w", newline="") as fh:
        fh.write(
            "# id: 'neg_ilp_residuals'\n"
            "# section_name: 'ILP Negative Sampler Residuals'\n"
            "# plot_type: 'table'\n"
            "# pconfig:\n"
            "#     id: 'neg_ilp_residuals_table'\n"
            "#     title: 'ILP Negative Sampler: per-protein degree residuals'\n"
        )
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            out = {k: row.get(k, "") for k in RESIDUAL_COLUMNS}
            out["Sample"] = f"{mqc_sample(id_, row.get('split', ''))}_{row.get('protein_id', '')}"
            out["ID"] = id_
            writer.writerow(out)


# ============================================================
# 7. Split driver
# ============================================================


def sample_negatives_ilp(
    name,
    pos_ppis,
    output_path,
    cfg: SamplingConfig,
    neg_ratio,
    species_path=None,
    go_annotations_path=None,
    candidate_network_path=None,
    gurobi_license_path=None,
    protein_to_idx=None,
    idx_to_protein=None,
    verbose_rows_out=None,
):
    """Sample negatives for one split. Returns (diagnostics_row, ctx)."""
    pos_pairs = pos_pairs_from_rows(pos_ppis, protein_to_idx)
    pos_pairs_set = {tuple(p) for p in pos_pairs.tolist()}

    # Pre-load what active biases need, so an over-budget subsample can use it (not a uniform draw); also passed to
    # build_context to avoid re-reads.
    taxonomy_relevant = cfg.lambda_taxon_pair > 0
    taxonomy_codes = n_taxa = None
    if taxonomy_relevant and species_path is not None:
        taxonomy = load_species(species_path, protein_to_idx)
        uniq, taxonomy_codes = np.unique(taxonomy, return_inverse=True)
        taxonomy_codes = taxonomy_codes.astype(np.int64)
        n_taxa = int(len(uniq))

    go_bp = go_membership = go_sizes = None
    if cfg.lambda_jaccard > 0 and go_annotations_path is not None:
        go_bp = load_go_bp(go_annotations_path, protein_to_idx)
        go_membership, go_sizes = _build_go_membership(go_bp)

    if candidate_network_path is not None:
        candidates = load_candidate_network(candidate_network_path, protein_to_idx, pos_pairs_set, cfg)
    else:
        candidates = build_candidate_set(
            len(protein_to_idx),
            pos_pairs,
            cfg,
            max_candidates=cfg.max_candidates,
            seed=cfg.seed,
            taxon_codes=taxonomy_codes,
            go_membership=go_membership,
            go_sizes=go_sizes,
        )
    n_pos = int(len(pos_pairs))
    n_neg = int(round(neg_ratio * n_pos))
    if n_pos > 0 and len(candidates) == 0:
        raise RuntimeError(
            f"{name}: need {n_neg} negatives but no candidate pairs are available. "
            f"Supply a larger --candidate-network or lower the negative ratio."
        )
    print(f"Candidate pool size: {len(candidates)}", file=sys.stderr)
    print(
        f"Number of unique proteins in the candidate set: {len(np.unique(candidates))}",
        file=sys.stderr,
    )
    print(
        f"Number of positive self-interactions: {sum(pos_pairs[:,0] == pos_pairs[:, 1])}",
        file=sys.stderr,
    )
    print(
        f"Number of possible self-interactions in the candidates: {sum(candidates[:,0] == candidates[:, 1])}",
        file=sys.stderr,
    )

    ctx = build_context(
        pos_pairs,
        protein_to_idx,
        idx_to_protein,
        candidates,
        neg_ratio,
        species_path=species_path,
        go_annotations_path=go_annotations_path,
        taxonomy_codes=taxonomy_codes,
        n_taxa=n_taxa,
        go_bp=go_bp,
    )

    if ctx.n_neg > len(ctx.candidates):
        raise RuntimeError(
            f"{name}: need {ctx.n_neg} negatives but only {len(ctx.candidates)} "
            f"candidate pairs are available. Supply a larger --candidate-network "
            f"or lower the negative ratio."
        )

    base_diag = {
        "split": name,
        "n_pos": ctx.n_pos,
        "n_neg": ctx.n_neg,
        "r": neg_ratio,
        "n_candidates": len(ctx.candidates),
    }

    if ctx.n_pos == 0:
        raise ValueError(f"{name}: no positive pairs found in the input.")
    if ctx.n_neg == 0:
        raise ValueError(f"{name}: no negative pairs found in the input.")

    if ctx.n_neg == len(ctx.candidates):
        logging.info(
            "%s: |NEG| == |C| (%d); selecting all candidates without solving.",
            name,
            ctx.n_neg,
        )
        write_split_csv(pos_ppis, ctx.candidates.tolist(), ctx.idx_to_protein, output_path)
        diag = {
            **base_diag,
            "obj_value": 0.0,
            "bias_deg_term": 0.0,
            "bias_tax_term": 0.0,
            "bias_self_term": 0.0,
            "bias_jac_term": 0.0,
            "solver": "trivial",
            "wall_time_s": 0.0,
            "mip_gap": 0.0,
            "status": "optimal (all candidates forced)",
        }
        print_objective_breakdown(name, diag)
        print_dataset_stats(name, ctx.pos_pairs, ctx.candidates, ctx, cfg)
        return diag, ctx

    biases = assemble_active_biases(cfg)
    for b in biases:
        b.precompute(ctx)
    active = [b for b in biases if b.is_active()]

    problem, x, term_exprs = build_problem(ctx, active)
    solver, options = select_solver(cfg, gurobi_license_path, cfg.verbose)
    try:
        result = solve(problem, solver, options, verbose=cfg.verbose)
    except RuntimeError as exc:
        raise RuntimeError(
            f"{exc} Note: every split enforces a hard per-protein negative-degree "
            f"cap of neg_ratio * {1 + MAX_DEGREE_SLACK} * d_plus; infeasibility "
            f"often means the positive set's structure (e.g. too many "
            f"self-interactions relative to available negative self-pairs) makes "
            f"exact ratio-matched negatives impossible. Consider lowering --neg-ratio."
        ) from exc
    negatives_idx = extract_negatives(x.value, ctx)

    term_values = {
        "bias_deg_term": 0.0,
        "bias_tax_term": 0.0,
        "bias_self_term": 0.0,
        "bias_jac_term": 0.0,
    }
    for b, raw_expr in term_exprs:
        term_values[_TERM_KEY[b.name]] += float(raw_expr.value)

    diag = {
        **base_diag,
        "obj_value": result["obj_value"],
        **term_values,
        "solver": str(solver),
        "wall_time_s": result["wall_time_s"],
        "mip_gap": cfg.mip_gap,
        "status": str(result["status"]),
    }

    if verbose_rows_out is not None:
        for b, _ in term_exprs:
            for row in b.debug_rows(x.value, ctx):
                verbose_rows_out.append({"split": name, **row})

    write_split_csv(pos_ppis, negatives_idx.tolist(), ctx.idx_to_protein, output_path)
    print_objective_breakdown(name, diag)
    print_dataset_stats(name, ctx.pos_pairs, negatives_idx, ctx, cfg)
    return diag, ctx


if __name__ == "__main__":
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="[%(asctime)s] %(levelname)s %(message)s",
    )
    cfg = config_from_args(args)
    split_name = args.split_name or Path(args.output).stem

    pos_ppis = read_ppis(args.positives)
    protein_to_idx, idx_to_protein = build_protein_index(pos_ppis)

    residual_rows = [] if cfg.verbose else None
    diag, _ = sample_negatives_ilp(
        split_name,
        pos_ppis,
        args.output,
        cfg,
        args.neg_ratio,
        species_path=args.species,
        go_annotations_path=args.go_annotations,
        candidate_network_path=args.candidate_network,
        gurobi_license_path=args.gurobi_license,
        protein_to_idx=protein_to_idx,
        idx_to_protein=idx_to_protein,
        verbose_rows_out=residual_rows,
    )
    write_diagnostics([diag], args.diagnostics_out, args.id)
    write_neg_generalstats(diag, f"{split_name}_gs_mqc.tsv", args.id)
    if cfg.verbose and residual_rows:
        write_residuals(residual_rows, args.residuals_out, args.id)
