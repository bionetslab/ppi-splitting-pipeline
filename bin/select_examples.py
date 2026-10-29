#!/usr/bin/env python3
"""
DDI mode: turn family-pair DDIs into domain-instance-pair examples.

Reads the three split CSVs (Pfam family pairs), the split FASTAs (domain instances),
instances.tsv and optionally a candidate_network of negative family pairs. Picks up to
N instance pairs per DDI and writes {split}_sel.csv, {split}_examples.csv,
{split}_candidate_examples.csv, {split}_universe.txt, {split}_reserve.txt and
unclaimed.txt, plus MultiQC tables.

Key constraint: a parent protein is used by at most one split. Selection is solved as a
lexicographic ILP per connected component of contested parents (keep every DDI, fill
DDIs to N, maximise examples less a parent-reuse penalty, then candidate negatives);
DDIs touching no contested parent are picked greedily, and components over the
candidate cap or time budget fall back to a deterministic greedy.
"""

import argparse
import math
import os
import random
import statistics
import sys
import time
from collections import defaultdict
from dataclasses import dataclass, field

import cvxpy as cp
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import (  # noqa: E402
    diverse_pick,
    instances_by_family,
    mqc_dataset,
    mqc_sample,
    pair_candidates,
    read_fasta,
    read_instances,
    read_ppis,
    write_ppi_csv,
)
from solve_ilp import _solver_options  # noqa: E402

SPLITS = ["train", "val", "test"]
EXAMPLE_COLUMNS = ["protein1", "protein2", "family1", "family2"]


@dataclass
class Unit:
    """One DDI (positive) or candidate_network pair, with its example candidates."""

    split: str
    kind: str  # "pos" | "cand"
    row: dict  # the source CSV row; protein1/protein2 hold the two families
    fam1: str
    fam2: str
    cands: list  # [(instance_a, instance_b)] after shortlisting
    n_raw: int  # candidates before the shortlist trim
    picked: list = field(default_factory=list)


# ---------------------------------------------------------------------------
# Candidate enumeration
# ---------------------------------------------------------------------------


def build_units(kind, rows_by_split, members, split_ids, parent_of, shortlist_k, seed):
    """Wrap each interaction row as a Unit, enumerating and shortlisting candidates."""
    units = []
    for split in SPLITS:
        for row in rows_by_split[split]:
            fam1, fam2 = row["protein1"], row["protein2"]
            cands = pair_candidates(fam1, fam2, members, split_ids[split])
            n_raw = len(cands)
            if n_raw > shortlist_k:
                # Seeded per unit, so the shortlist does not depend on processing order.
                rng = random.Random(f"{seed}:{kind}:{split}:{fam1}:{fam2}")
                cands = diverse_pick(cands, shortlist_k, parent_of, rng)
            units.append(Unit(split, kind, row, fam1, fam2, cands, n_raw))
    return units


# ---------------------------------------------------------------------------
# The selection ILP
# ---------------------------------------------------------------------------


def _sum(var, idx):
    """cp.sum over an index list, tolerating the empty list."""
    return cp.sum(var[idx]) if idx else cp.Constant(0.0)


_WARNED = set()


def _warn_once(msg, key=None):
    """Print a warning only the first time its `key` (default: the message) is seen."""
    key = msg if key is None else key
    if key not in _WARNED:
        _WARNED.add(key)
        print(msg, file=sys.stderr)


@dataclass
class SolverEnv:
    """One solver environment reused across component solves, so the Gurobi licence is checked once.

    Dropped after the first solve that rejects the `env` argument.
    """

    env: object = None

    def kwargs(self):
        return {"env": self.env} if self.env is not None else {}

    def disable(self, exc):
        print(
            f"Note: the solver did not accept a shared environment ({exc!r}); building one per solve.",
            file=sys.stderr,
        )
        self.env = None


def make_solver_env(solver):
    """Return a shared gurobipy.Env when the solver is Gurobi, an inert holder otherwise."""
    if (solver or "").upper() != cp.GUROBI:
        return SolverEnv(None)
    try:
        import gurobipy

        return SolverEnv(gurobipy.Env())
    except Exception as exc:  # gurobipy absent, or the licence check failed here
        print(f"Note: no shared Gurobi environment ({exc!r}); each solve builds its own.", file=sys.stderr)
        return SolverEnv(None)


@dataclass
class SolveContext:
    """Solver settings and run-level counters, threaded through every stage solve."""

    solver: str = None
    seed: int = 42
    verbose: bool = False
    env: SolverEnv = field(default_factory=SolverEnv)
    stages: int = 0  # stage solves attempted
    timeouts: int = 0  # of those, the ones that returned an incumbent at the time limit


def _run_stage(problem, label, secs, ctx):
    """Solve one lexicographic stage; return its objective value, or None if no usable solution."""
    # Solver-specific names for time limit and seed; SCIP rejects the generic time_limit kwarg.
    kwargs = {"verbose": ctx.verbose}
    if ctx.solver:
        solver_opts = _solver_options(ctx.solver, ctx.seed, secs)
        if not solver_opts:
            _warn_once(
                f"Warning: no solver-specific options known for {ctx.solver}; time limit passed generically "
                f"and the selection is unseeded."
            )
            solver_opts = {"time_limit": secs}
        kwargs.update(solver_opts)
        kwargs["solver"] = ctx.solver
    else:
        _warn_once("Warning: no --solver given, so the selection is unseeded; pass --solver for reproducibility.")
        kwargs["time_limit"] = secs

    ctx.stages += 1
    try:
        try:
            problem.solve(**kwargs, **ctx.env.kwargs())
        except Exception as exc:
            if not ctx.env.kwargs():
                raise
            ctx.env.disable(exc)
            problem.solve(**kwargs)
    except cp.error.SolverError as exc:
        # e.g. a time limit hit before the first incumbent: fall back to the greedy.
        _warn_once(
            f"Warning: stage '{label}' raised {type(exc).__name__} ({exc}) at a {secs:.0f}s limit; that "
            f"component falls back to the greedy.",
            key=f"solvererror:{label}",
        )
        return None

    if problem.status not in cp.settings.SOLUTION_PRESENT:
        _warn_once(
            f"Warning: stage '{label}' returned status {problem.status} with no usable solution at a "
            f"{secs:.0f}s limit; that component falls back to the greedy.",
            key=f"nosolution:{label}",
        )
        return None
    if problem.status == cp.settings.USER_LIMIT:
        ctx.timeouts += 1
        _warn_once(
            f"Warning: stage '{label}' hit its {secs:.0f}s limit; using the best incumbent found.",
            key=f"userlimit:{label}",
        )
    if ctx.verbose:
        print(f"  stage '{label}': objective {problem.value:,.4g} ({problem.status})", file=sys.stderr)
    return float(problem.value)


def solve_selection(units, parent_of, contested, n, lam, max_sec, ctx):
    """Choose <= n examples per unit, one split per parent protein. Fills Unit.picked.

    Variables
        y[e]   in {0,1}  candidate example e is selected
        nz[d]  in {0,1}  unit d ends up with >= 1 example
        r[d]   in {0,1}  unit d reaches the full n
        c[p,s] in {0,1}  contested parent p is claimed by split s
        o[d,p] >= 0      how often unit d reuses parent p beyond the first time

    Constraints
        (1) sum_{e in d} y[e] <= min(n, |C_d|)          per-unit cap (n is a cap, not a quota)
        (2) sum_{e in d} y[e] >= nz[d]
        (3) sum_{e in d} y[e] >= n * r[d]
        (4) sum_s c[p,s] <= 1                           one split per parent protein
        (5) y[e] <= c[parent(e), split(d(e))]           both parents, contested only
        (6) o[d,p] >= (uses of p by d's selected examples) - 1

    Returns False (units untouched) when a stage had no usable solution.
    """
    cand_index, unit_rows, unit_cols, caps = [], [], [], []
    for ui, u in enumerate(units):
        for pair in u.cands:
            unit_rows.append(ui)
            unit_cols.append(len(cand_index))
            cand_index.append((ui, pair))
        caps.append(min(n, len(u.cands)))
    n_units, n_cand = len(units), len(cand_index)

    # Claim variables only for contested parents; for the others constraint 5 is slack.
    claim_index = {}
    for ui, u in enumerate(units):
        for a, b in u.cands:
            for p in (parent_of[a], parent_of[b]):
                if p in contested:
                    claim_index.setdefault((p, u.split), len(claim_index))

    y = cp.Variable(n_cand, boolean=True)
    nz = cp.Variable(n_units, boolean=True)
    r = cp.Variable(n_units, boolean=True)

    per_unit = sp.coo_matrix((np.ones(n_cand), (unit_rows, unit_cols)), shape=(n_units, n_cand)).tocsr()
    unit_examples = per_unit @ y
    cons = [unit_examples <= np.array(caps, dtype=float), unit_examples >= nz, unit_examples >= n * r]

    c = cp.Variable(len(claim_index), boolean=True)
    by_prot = defaultdict(list)
    for (p, _), idx in claim_index.items():
        by_prot[p].append(idx)
    rows, cols = [], []
    for pi, (_p, idxs) in enumerate(sorted(by_prot.items())):
        for idx in idxs:
            rows.append(pi)
            cols.append(idx)
    one_split_per_parent = sp.coo_matrix(
        (np.ones(len(rows)), (rows, cols)), shape=(len(by_prot), len(claim_index))
    ).tocsr()
    cons.append(one_split_per_parent @ c <= np.ones(len(by_prot)))

    e_idx, cl_idx = [], []
    for k, (ui, (a, b)) in enumerate(cand_index):
        for p in {parent_of[a], parent_of[b]}:
            if p in contested:
                e_idx.append(k)
                cl_idx.append(claim_index[(p, units[ui].split)])
    cons.append(y[e_idx] <= c[cl_idx])

    # Diversity: one overflow variable per (unit, parent); coo_matrix sums duplicates, so a
    # same-parent pair counts twice. Keys with a single candidate are always slack and dropped.
    key_count = defaultdict(int)
    for ui, (a, b) in cand_index:
        for p in (parent_of[a], parent_of[b]):
            key_count[(ui, p)] += 1
    over_index, o_rows, o_cols = {}, [], []
    for k, (ui, (a, b)) in enumerate(cand_index):
        for p in (parent_of[a], parent_of[b]):
            key = (ui, p)
            if key_count[key] < 2:
                continue
            if key not in over_index:
                over_index[key] = len(over_index)
            o_rows.append(over_index[key])
            o_cols.append(k)
    o = None
    if over_index:
        o = cp.Variable(len(over_index), nonneg=True)
        per_parent = sp.coo_matrix((np.ones(len(o_rows)), (o_rows, o_cols)), shape=(len(over_index), n_cand)).tocsr()
        cons.append(per_parent @ y - 1.0 <= o)

    pos_u = [i for i, u in enumerate(units) if u.kind == "pos"]
    cand_u = [i for i, u in enumerate(units) if u.kind == "cand"]
    pos_y = [k for k, (ui, _) in enumerate(cand_index) if units[ui].kind == "pos"]
    cand_y = [k for k, (ui, _) in enumerate(cand_index) if units[ui].kind == "cand"]
    pos_o = [idx for (ui, _), idx in over_index.items() if units[ui].kind == "pos"]
    cand_o = [idx for (ui, _), idx in over_index.items() if units[ui].kind == "cand"]

    # Lexicographic objective, one stage per level, so every coefficient stays 1 or lam.
    stages = []
    if pos_u:
        stages.append(("keep every positive DDI", _sum(nz, pos_u), 0.15))
        stages.append(("fill positives to N", _sum(r, pos_u), 0.35))
        stages.append(("maximise positive examples", _sum(y, pos_y) - lam * _sum(o, pos_o), 0.25))
    if cand_u:
        stages.append(("expand candidate negatives", _sum(y, cand_y) - lam * _sum(o, cand_o), 0.25))

    total_share = sum(s for _, _, s in stages)
    for label, expr, share in stages:
        secs = max(1.0, max_sec * share / total_share)
        value = _run_stage(cp.Problem(cp.Maximize(expr), cons), label, secs, ctx)
        if value is None:
            return False
        # Freeze this level before optimising the next, with a relative tolerance.
        cons = cons + [expr >= value - max(1e-4, 1e-6 * abs(value))]

    y_val = np.asarray(y.value).ravel()
    for k, (ui, pair) in enumerate(cand_index):
        if y_val[k] > 0.5:
            units[ui].picked.append(pair)
    for ui, u in enumerate(units):
        # Trim in case the solver returns 0.5+eps on more candidates than the cap.
        u.picked = sorted(u.picked)[: caps[ui]]
    return True


# ---------------------------------------------------------------------------
# Decomposition, fallback and the driver over both
# ---------------------------------------------------------------------------


def component_partition(units, parent_of, contested):
    """Cut the ILP units into connected components over the contested parents they share.

    Components share no variable or constraint, so solving them separately is exact.
    Returns lists of Units, largest candidate count first.
    """
    uf = list(range(len(units)))

    def find(x):
        while uf[x] != x:
            uf[x] = uf[uf[x]]  # path halving
            x = uf[x]
        return x

    first_seen = {}
    for ui, u in enumerate(units):
        for a, b in u.cands:
            for p in (parent_of[a], parent_of[b]):
                if p not in contested:
                    continue
                if p in first_seen:
                    ra, rb = find(ui), find(first_seen[p])
                    if ra != rb:
                        uf[ra] = rb
                else:
                    first_seen[p] = ui

    groups = defaultdict(list)
    for ui in range(len(units)):
        groups[find(ui)].append(ui)
    # Deterministic order, since the time budget decides where the greedy takes over.
    keyed = sorted((-sum(len(units[i].cands) for i in idxs), idxs[0], idxs) for idxs in groups.values())
    return [[units[i] for i in idxs] for _, _, idxs in keyed]


def greedy_component(units, parent_of, contested, n, seed, tag):
    """Deterministic greedy fallback for one component; keeps one split per parent, loses optimality.

    Positives first, seeded shuffle within each kind; each unit picks from candidates whose
    parents no other split has claimed yet.
    """
    rng = random.Random(f"{seed}:greedy:{tag}")
    order = []
    for kind in ("pos", "cand"):
        group = [u for u in units if u.kind == kind]
        rng.shuffle(group)
        order += group

    claimed = {}
    for u in order:
        usable = [
            (a, b)
            for a, b in u.cands
            if claimed.get(parent_of[a], u.split) == u.split and claimed.get(parent_of[b], u.split) == u.split
        ]
        pick_rng = random.Random(f"{seed}:greedy:{u.kind}:{u.split}:{u.fam1}:{u.fam2}")
        u.picked = sorted(diverse_pick(usable, min(n, len(usable)), parent_of, pick_rng))
        for a, b in u.picked:
            for p in (parent_of[a], parent_of[b]):
                if p in contested:
                    claimed[p] = u.split


def run_selection(units, parent_of, contested, n, lam, max_sec, max_cand, ctx):
    """Solve every component, falling back to the greedy over the candidate cap or time budget.

    The budget is shared in proportion to candidate count, recomputed after each component.
    Returns the numbers for the DDI Example Selection ILP MultiQC table.
    """
    comps = component_partition(units, parent_of, contested)
    sizes = [sum(len(u.cands) for u in c) for c in comps]
    n_units = [len(c) for c in comps]
    print(
        f"  {len(comps):,} independent component(s); units min/median/max "
        f"{min(n_units):,}/{statistics.median(n_units):,.0f}/{max(n_units):,}, candidates "
        f"{min(sizes):,}/{statistics.median(sizes):,.0f}/{max(sizes):,}. Largest: "
        f"{sizes[0]:,} candidates over {n_units[0]:,} units.",
        file=sys.stderr,
    )

    started = time.monotonic()
    fallback, fallback_units, remaining_cand = 0, 0, sum(sizes)
    for ci, comp in enumerate(comps):
        elapsed = time.monotonic() - started
        kind, reason = None, None
        if sizes[ci] > max_cand:
            kind = "candidate cap"
            reason = f"{sizes[ci]:,} candidates is over --max-ilp-candidates ({max_cand:,})"
        elif elapsed >= max_sec:
            kind = "time budget"
            reason = f"the {max_sec}s ILP budget ran out after {ci:,} of {len(comps):,} components"
        if kind is None:
            share = (max_sec - elapsed) * sizes[ci] / remaining_cand
            if not solve_selection(comp, parent_of, contested, n, lam, share, ctx):
                kind = "no solution"
                reason = "the solver returned no usable solution"
        remaining_cand -= sizes[ci]
        if kind is not None:
            fallback += 1
            fallback_units += len(comp)
            _warn_once(
                f"WARNING: at least one component fell back to the greedy selection ({reason}); "
                f"its selection is no longer optimal.",
                key=f"fallback:{kind}",
            )
            greedy_component(comp, parent_of, contested, n, ctx.seed, str(ci))

    info = {
        "units": len(units),
        "components": len(comps),
        "largest_units": n_units[0],
        "largest_cands": sizes[0],
        "fallback": fallback,
        "fallback_units": fallback_units,
        "timeouts": ctx.timeouts,
        "seconds": time.monotonic() - started,
    }
    print(
        f"  ILP done in {info['seconds']:,.1f}s over {ctx.stages:,} stage solves "
        f"({info['timeouts']:,} hit their time limit); {fallback:,} component(s) "
        f"covering {fallback_units:,} units fell back to the greedy.",
        file=sys.stderr,
    )
    return info


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def example_rows(units):
    """Return instance-pair rows (protein1/protein2 instances, family1/family2, extra input columns)."""
    rows = []
    for u in units:
        extra = {k: v for k, v in u.row.items() if k not in ("protein1", "protein2")}
        for a, b in u.picked:
            rows.append({"protein1": a, "protein2": b, "family1": u.fam1, "family2": u.fam2, **extra})
    return rows


def write_examples(units, path):
    """write_ppi_csv, but keep the DDI header for an empty split so it is not read as PPI mode."""
    rows = example_rows(units)
    if rows:
        write_ppi_csv(rows, path)
        return len(rows)
    with open(path, "w", newline="") as fh:
        fh.write(",".join(EXAMPLE_COLUMNS) + "\n")
    return 0


def write_ids(ids, path):
    with open(path, "w") as fh:
        for i in sorted(ids):
            fh.write(f"{i}\n")


def partition_reserve(unclaimed, weights, seed, allow_shared):
    """Hand each never-in-play parent protein to one split, weighted by the DDIs each split kept.

    Under --allow-shared-parents every split gets the whole pool.
    """
    if allow_shared:
        return {s: set(unclaimed) for s in SPLITS}
    out = {s: set() for s in SPLITS}
    active = [s for s in SPLITS if weights.get(s, 0) > 0]
    if not active:
        return out
    total = float(sum(weights[s] for s in active))
    rng = random.Random(f"{seed}:reserve")
    for p in sorted(unclaimed):
        x = rng.random() * total
        acc = 0.0
        for s in active:
            acc += weights[s]
            if x < acc:
                out[s].add(p)
                break
    return out


def write_mqc(stats, ilp, id_):
    """Write the per-split outcome bar, per-split stats table and per-dataset ILP table.

    Shared section ids: MultiQC merges every dataset into one section each."""
    with open("select_examples_bar_mqc.tsv", "w") as fh:
        fh.write(
            "# id: 'ddi_examples_bar'\n"
            "# section_name: 'DDI Example Selection'\n"
            "# description: 'Positive DDIs per split, coloured by how many domain-instance "
            "examples each one ended up with. A DDI that reached zero examples -- every "
            "candidate blocked because another split already had its parent proteins, or the "
            "family had no usable instance -- is dropped from the split.'\n"
            "# plot_type: 'bargraph'\n"
            "# pconfig:\n"
            "#     id: 'ddi_examples_bar_plot'\n"
            "#     title: 'DDI Example Selection: DDIs per outcome'\n"
            "#     ylab: '# DDIs'\n"
            "Sample\tFull N examples\tPartial\tDropped (0 examples)\n"
        )
        for split in SPLITS:
            st = stats[split]
            fh.write(f"{mqc_sample(id_, split)}\t{st['full']}\t{st['partial']}\t{st['dropped']}\n")

    with open("select_examples_stats_mqc.tsv", "w") as fh:
        fh.write(
            "# id: 'ddi_examples_stats'\n"
            "# section_name: 'DDI Example Selection Stats'\n"
            "# description: 'Per split: the DDIs that survived selection, the instance-pair "
            "examples written for them, the parent proteins claimed for this split (its protein "
            "universe, which the negative sampler draws from), and how many DDIs had their "
            "candidate pool shortlisted before the ILP saw it.'\n"
            "# plot_type: 'table'\n"
            "# pconfig:\n"
            "#     id: 'ddi_examples_stats_table'\n"
            "#     title: 'DDI Example Selection Stats'\n"
            "Sample\tID\tDDIs in\tDDIs kept\tExamples\tExamples per DDI\tProtein universe\t"
            "Contested parents\tShortlisted DDIs\tCandidate pairs\tCandidate examples\n"
        )
        for split in SPLITS:
            st = stats[split]
            per = st["examples"] / st["kept"] if st["kept"] else 0.0
            fh.write(
                f"{mqc_sample(id_, split)}\t{id_}\t{st['ddis_in']}\t{st['kept']}\t{st['examples']}\t"
                f"{per:.2f}\t{st['universe']}\t{st['contested']}\t{st['shortlisted']}\t"
                f"{st['cand_pairs']}\t{st['cand_examples']}\n"
            )

    with open("select_examples_ilp_mqc.tsv", "w") as fh:
        fh.write(
            "# id: 'ddi_examples_ilp'\n"
            "# section_name: 'DDI Example Selection ILP'\n"
            "# description: 'Shape of the selection ILP, one row per dataset. Units are the "
            "interactions left after the reduction that settles anything touching no contested "
            "parent locally; those are cut into independent components over the contested parents "
            "they share, and each component is solved on its own. A component that exceeds "
            "--max-ilp-candidates, or that is reached after the time budget is spent, falls back to "
            "a deterministic greedy: one split per parent protein still holds there, but that "
            "component is no longer optimal, so a nonzero fallback count is worth acting on.'\n"
            "# plot_type: 'table'\n"
            "# pconfig:\n"
            "#     id: 'ddi_examples_ilp_table'\n"
            "#     title: 'DDI Example Selection ILP'\n"
            "Sample\tID\tUnits in ILP\tComponents\tLargest component (units)\tLargest component (candidates)\t"
            "Greedy fallback components\tGreedy fallback units\tStages at time limit\tILP seconds\n"
        )
        fh.write(
            f"{mqc_dataset(id_)}\t{id_}\t{ilp['units']}\t{ilp['components']}\t{ilp['largest_units']}\t"
            f"{ilp['largest_cands']}\t{ilp['fallback']}\t{ilp['fallback_units']}\t{ilp['timeouts']}\t"
            f"{ilp['seconds']:.1f}\n"
        )


# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--train_ppis", required=True, help="train split CSV (Pfam family pairs)")
    ap.add_argument("--val_ppis", required=True)
    ap.add_argument("--test_ppis", required=True)
    ap.add_argument("--train_fasta", required=True, help="train split FASTA (domain instances)")
    ap.add_argument("--val_fasta", required=True)
    ap.add_argument("--test_fasta", required=True)
    ap.add_argument("--instances", required=True, help="instances.tsv: family/clan/instance/parent table")
    ap.add_argument(
        "--candidate-network",
        default=None,
        help="high-confidence negative family pairs, selected in the same ILP as the positives",
    )
    ap.add_argument("--examples-target", type=int, default=5, help="N: cap on examples per DDI (default 5)")
    ap.add_argument(
        "--shortlist-factor",
        type=int,
        default=4,
        help="cap a DDI's candidate pool at this multiple of N before the ILP (default 4)",
    )
    ap.add_argument(
        "--candidate-factor",
        type=float,
        default=4.0,
        help="cap candidate_network pairs per split at this multiple of the split's DDI count "
        "(default 4), mirroring SAMPLE_NEGATIVES_ILP's own candidate cap",
    )
    ap.add_argument(
        "--lambda-diversity",
        type=float,
        default=0.1,
        help="weight on parent reuse; must be < 0.5 so diversity never costs an example (default 0.1)",
    )
    ap.add_argument(
        "--allow-shared-parents",
        action="store_true",
        help="let a parent protein be used by several splits (for split_method=random)",
    )
    ap.add_argument(
        "--max-sec",
        type=int,
        default=300,
        help="total ILP time budget in seconds, shared across components (default 300)",
    )
    ap.add_argument(
        "--max-ilp-candidates",
        type=int,
        default=200000,
        help="use the greedy for any component with more candidate pairs than this (default 200000)",
    )
    ap.add_argument("--solver", default=None, help="CVXPY solver name, e.g. GUROBI, SCIP (default: auto)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--id", required=True, help="Dataset ID, for MultiQC tagging")
    ap.add_argument("--verbose", action="store_true", help="let the solver print its own log")
    args = ap.parse_args()

    n = args.examples_target
    if n < 1:
        sys.exit("--examples-target must be at least 1")
    if args.lambda_diversity >= 0.5:
        sys.exit(f"--lambda-diversity must be < 0.5 (got {args.lambda_diversity})")

    inst_rows = read_instances(args.instances)
    members = instances_by_family(inst_rows)
    parent_of = {r["instance_id"]: r["protein_id"] for r in inst_rows}
    all_parents = set(parent_of.values())
    print(
        f"{len(inst_rows):,} instances over {len(members):,} families and {len(all_parents):,} parent proteins",
        file=sys.stderr,
    )

    ppi_paths = {"train": args.train_ppis, "val": args.val_ppis, "test": args.test_ppis}
    fasta_paths = {"train": args.train_fasta, "val": args.val_fasta, "test": args.test_fasta}
    rows_by_split = {s: read_ppis(p) for s, p in ppi_paths.items()}
    split_ids = {s: set(read_fasta(p)) for s, p in fasta_paths.items()}

    shortlist_k = max(n, args.shortlist_factor * n)
    units = build_units("pos", rows_by_split, members, split_ids, parent_of, shortlist_k, args.seed)

    cand_by_split = {s: [] for s in SPLITS}
    if args.candidate_network:
        # A set per family: under split_method=random a family can be in several splits.
        fam_splits = defaultdict(set)
        for s in SPLITS:
            for row in rows_by_split[s]:
                fam_splits[row["protein1"]].add(s)
                fam_splits[row["protein2"]].add(s)
        n_unplaced = 0
        for row in read_ppis(args.candidate_network):
            shared = fam_splits.get(row["protein1"], set()) & fam_splits.get(row["protein2"], set())
            for s in shared:
                cand_by_split[s].append(row)
            if not shared:
                n_unplaced += 1
        for s in SPLITS:
            limit = int(math.ceil(args.candidate_factor * len(rows_by_split[s])))
            if len(cand_by_split[s]) > limit:
                rng = random.Random(f"{args.seed}:candidate_network:{s}")
                keep = rng.sample(range(len(cand_by_split[s])), limit)
                cand_by_split[s] = [cand_by_split[s][i] for i in sorted(keep)]
                print(
                    f"  {s}: candidate_network capped to {limit:,} pairs "
                    f"(--candidate-factor {args.candidate_factor})",
                    file=sys.stderr,
                )
        print(
            f"  candidate_network: {sum(len(v) for v in cand_by_split.values()):,} pairs placed in a "
            f"split, {n_unplaced:,} skipped (families in different splits or in none)",
            file=sys.stderr,
        )
        units += build_units("cand", cand_by_split, members, split_ids, parent_of, shortlist_k, args.seed)

    # A parent is contested when it carries candidates in more than one split.
    splits_of = defaultdict(set)
    for u in units:
        for a, b in u.cands:
            splits_of[parent_of[a]].add(u.split)
            splits_of[parent_of[b]].add(u.split)
    contested = set() if args.allow_shared_parents else {p for p, ss in splits_of.items() if len(ss) > 1}
    if args.allow_shared_parents:
        print(
            "--allow-shared-parents: parent proteins may be shared across splits; the ILP is skipped.",
            file=sys.stderr,
        )
    print(
        f"{len(units):,} units ({sum(1 for u in units if u.kind == 'pos'):,} positive); "
        f"{len(splits_of):,} parent proteins in play, {len(contested):,} contested",
        file=sys.stderr,
    )

    free = [
        u for u in units if all(parent_of[a] not in contested and parent_of[b] not in contested for a, b in u.cands)
    ]
    # By identity: two Units of the same DDI in different splits may compare equal.
    free_ids = {id(u) for u in free}
    ilp_units = [u for u in units if id(u) not in free_ids]
    print(f"  {len(free):,} decoupled (local greedy), {len(ilp_units):,} in the ILP", file=sys.stderr)

    for u in free:
        rng = random.Random(f"{args.seed}:free:{u.kind}:{u.split}:{u.fam1}:{u.fam2}")
        u.picked = sorted(diverse_pick(u.cands, min(n, len(u.cands)), parent_of, rng))

    ilp = dict.fromkeys(
        ("units", "components", "largest_units", "largest_cands", "fallback", "fallback_units", "timeouts"), 0
    )
    ilp["seconds"] = 0.0
    if ilp_units:
        print("Solving the selection ILP …", file=sys.stderr)
        ctx = SolveContext(solver=args.solver, seed=args.seed, verbose=args.verbose, env=make_solver_env(args.solver))
        ilp = run_selection(
            ilp_units,
            parent_of,
            contested,
            n,
            args.lambda_diversity,
            args.max_sec,
            args.max_ilp_candidates,
            ctx,
        )
    else:
        why = "--allow-shared-parents" if args.allow_shared_parents else "no contested parent protein"
        print(f"ILP skipped ({why}): every unit's examples are decided locally.", file=sys.stderr)

    # Per-split protein universe: uncontested parents in play for the split plus the contested
    # parents its selected examples used (not c[p,s], which may claim unused parents).
    universe = {s: set() for s in SPLITS}
    for p, ss in splits_of.items():
        if p in contested:
            continue
        for s in ss:
            universe[s].add(p)
    for u in units:
        for a, b in u.picked:
            for p in (parent_of[a], parent_of[b]):
                if p in contested:
                    universe[u.split].add(p)

    stats = {}
    for split in SPLITS:
        pos = [u for u in units if u.kind == "pos" and u.split == split]
        cand = [u for u in units if u.kind == "cand" and u.split == split]
        kept = [u for u in pos if u.picked]
        write_ppi_csv([u.row for u in kept], f"{split}_sel.csv")
        n_examples = write_examples(kept, f"{split}_examples.csv")
        n_cand_examples = write_examples([u for u in cand if u.picked], f"{split}_candidate_examples.csv")
        write_ids(universe[split], f"{split}_universe.txt")
        stats[split] = {
            "ddis_in": len(pos),
            "kept": len(kept),
            "dropped": len(pos) - len(kept),
            "full": sum(1 for u in pos if len(u.picked) >= n),
            "partial": sum(1 for u in pos if 0 < len(u.picked) < n),
            "examples": n_examples,
            "universe": len(universe[split]),
            "contested": len(universe[split] & contested),
            "shortlisted": sum(1 for u in pos if u.n_raw > len(u.cands)),
            "cand_pairs": len(cand),
            "cand_examples": n_cand_examples,
        }
        st = stats[split]
        print(
            f"  {split}: {st['kept']:,}/{st['ddis_in']:,} DDIs kept "
            f"({st['full']:,} at full N, {st['partial']:,} partial, {st['dropped']:,} dropped), "
            f"{st['examples']:,} examples, {st['universe']:,} proteins claimed",
            file=sys.stderr,
        )

    # Parents no candidate ever reached. A contested parent the ILP left unused is in neither,
    # since handing it to a split now would reintroduce the leak.
    unclaimed = all_parents - set(splits_of)
    write_ids(unclaimed, "unclaimed.txt")
    reserve = partition_reserve(unclaimed, {s: stats[s]["kept"] for s in SPLITS}, args.seed, args.allow_shared_parents)
    for split in SPLITS:
        write_ids(reserve[split], f"{split}_reserve.txt")
    print(
        f"{len(unclaimed):,} parent proteins never in play, split "
        + ", ".join(f"{len(reserve[s]):,} {s}" for s in SPLITS)
        + " (free for the negative sampler to extend a split's universe with)",
        file=sys.stderr,
    )

    write_mqc(stats, ilp, args.id)


if __name__ == "__main__":
    main()
