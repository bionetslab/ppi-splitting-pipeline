#!/usr/bin/env python3
"""
DDI mode: turn a family-level labelled CSV into domain-instance pairs.

Positives reuse the instance pairs SELECT_EXAMPLES chose. Negatives reuse
SELECT_EXAMPLES' candidate_network pairs where available; otherwise up to N
instance pairs are drawn from this split's protein universe, and, if that falls
short, from its reserve of unclaimed proteins. Writes the instance-level
labelled CSV and MultiQC tables; the example-level class ratio is reported, not
rebalanced.
"""

import argparse
import csv
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import (  # noqa: E402
    diverse_pick,
    instances_by_family,
    mqc_sample,
    pair_candidates,
    read_fasta,
    read_instances,
    read_ppis,
)

SPLITS = ["train", "val", "test"]
BASE_COLUMNS = ["protein1", "protein2", "label", "family1", "family2"]


def unordered(a, b):
    """Order-independent key for a pair, so (A, B) and (B, A) match."""
    return (min(a, b), max(a, b))


def universe_split(label):
    """Return the SELECT_EXAMPLES split a negative-sampling label draws from (test_* -> test)."""
    return "test" if label.startswith("test") else label


def read_ids(path):
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def read_header(path):
    """Column names of a CSV, so an empty split still knows its extra columns."""
    with open(path) as fh:
        return next(csv.reader(fh))


def examples_by_family_pair(path):
    """{(fam1, fam2): [(instance_a, instance_b)]} from a SELECT_EXAMPLES table."""
    out = defaultdict(list)
    for row in read_ppis(path):
        out[unordered(row["family1"], row["family2"])].append((row["protein1"], row["protein2"]))
    return out


def pick_pairs(fam1, fam2, members, available, parent_of, n, label, seed):
    """Up to n diverse instance pairs for one family pair, from `available`."""
    cands = pair_candidates(fam1, fam2, members, available)
    rng = random.Random(f"{seed}:neg:{label}:{fam1}:{fam2}")
    return sorted(diverse_pick(cands, min(n, len(cands)), parent_of, rng))


def write_expanded(records, extras, path):
    """Write instance-pair rows with label, family pair and extras; an empty split keeps the full header."""
    fieldnames = BASE_COLUMNS + extras
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for a, b, label, fam1, fam2, extra in records:
            out = {f: extra.get(f, "") for f in extras}
            out.update({"protein1": a, "protein2": b, "label": label, "family1": fam1, "family2": fam2})
            writer.writerow(out)


def write_mqc(label, st, id_):
    """One table row and one bar row per task; MultiQC merges them run-wide."""
    sample = mqc_sample(id_, label)
    with open(f"{label}_expand_table_mqc.tsv", "w") as fh:
        fh.write(
            "# id: 'ddi_expand_stats'\n"
            "# section_name: 'DDI Instance Expansion'\n"
            "# description: 'Family-level pairs from the negative sampler turned into domain-instance "
            "pairs. Positives reuse SELECT_EXAMPLES choices; negatives draw carriers from the split "
            "protein universe, so a negative pair can reach fewer examples than a positive one -- the "
            "two ratio columns are what shows whether the split is still balanced at example level.'\n"
            "# plot_type: 'table'\n"
            "# pconfig:\n"
            "#     id: 'ddi_expand_stats_table'\n"
            "#     title: 'DDI Instance Expansion'\n"
            "Sample\tID\tPositive pairs\tPositive examples\tNegative pairs\tNegative examples\t"
            "Negative pairs dropped\tFrom candidate_network\tExtended with unclaimed\t"
            "Pair ratio (family)\tExample ratio (instance)\n"
            f"{sample}\t{id_}\t{st['pos_pairs']}\t{st['pos_examples']}\t{st['neg_pairs']}\t"
            f"{st['neg_examples']}\t{st['neg_dropped']}\t{st['from_cand']}\t{st['extended']}\t"
            f"{st['target_ratio']:.3f}\t{st['achieved_ratio']:.3f}\n"
        )

    with open(f"{label}_expand_bar_mqc.tsv", "w") as fh:
        fh.write(
            "# id: 'ddi_expand_bar'\n"
            "# section_name: 'DDI Instance Expansion: Negative Pairs'\n"
            "# description: 'Sampled negative family pairs per split, by how many domain-instance "
            "examples each reached. A pair at zero had no instance of either family whose parent "
            "protein belongs to this split -- another split already had those parents -- and is "
            "dropped from the labelled CSV.'\n"
            "# plot_type: 'bargraph'\n"
            "# pconfig:\n"
            "#     id: 'ddi_expand_bar_plot'\n"
            "#     title: 'DDI Instance Expansion: negative pairs per outcome'\n"
            "#     ylab: '# Negative pairs'\n"
            "Sample\tFull N examples\tPartial\tDropped (0 examples)\n"
            f"{sample}\t{st['neg_full']}\t{st['neg_partial']}\t{st['neg_dropped']}\n"
        )


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--labelled", required=True, help="family-level labelled CSV from SAMPLE_NEGATIVES")
    ap.add_argument("--examples", required=True, help="SELECT_EXAMPLES {split}_examples.csv (positives)")
    ap.add_argument(
        "--candidate-examples",
        required=True,
        help="SELECT_EXAMPLES {split}_candidate_examples.csv -- pairs already claimed inside its ILP",
    )
    ap.add_argument("--universe", required=True, help="SELECT_EXAMPLES {split}_universe.txt")
    ap.add_argument(
        "--reserve",
        required=True,
        help="SELECT_EXAMPLES {split}_reserve.txt: this split's share of the unclaimed parent proteins",
    )
    ap.add_argument("--instances", required=True, help="instances.tsv")
    ap.add_argument("--fasta", required=True, help="this split's instance FASTA")
    ap.add_argument("--output", required=True, help="output instance-level labelled CSV")
    ap.add_argument("--split-name", required=True, help="train | val | test_balanced | test_realistic")
    ap.add_argument("--examples-target", type=int, default=5, help="N: cap on examples per pair (default 5)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--id", required=True, help="Dataset ID, for MultiQC tagging")
    args = ap.parse_args()

    n = args.examples_target
    if n < 1:
        sys.exit("--examples-target must be at least 1")
    label = args.split_name
    split = universe_split(label)

    inst_rows = read_instances(args.instances)
    members = instances_by_family(inst_rows)
    parent_of = {r["instance_id"]: r["protein_id"] for r in inst_rows}

    rows = read_ppis(args.labelled)
    header = read_header(args.labelled)
    extras = [c for c in header if c not in BASE_COLUMNS]
    pos_rows = [r for r in rows if int(r["label"]) == 1]
    neg_rows = [r for r in rows if int(r["label"]) == 0]

    pos_examples = examples_by_family_pair(args.examples)
    cand_examples = examples_by_family_pair(args.candidate_examples)
    universe = read_ids(args.universe)
    share = read_ids(args.reserve)

    # Instances usable for a negative: in the split's FASTA and parented by an allowed protein.
    in_split = set(read_fasta(args.fasta))
    from_universe = {i for i in in_split if parent_of[i] in universe}
    from_extended = from_universe | {i for i in in_split if parent_of[i] in share}
    has_reserve = from_extended != from_universe
    print(
        f"{label}: {len(pos_rows):,} positive and {len(neg_rows):,} negative family pairs; "
        f"{len(universe):,} proteins in the {split} universe -> {len(from_universe):,} usable "
        f"instances ({len(share):,} unclaimed proteins in reserve)",
        file=sys.stderr,
    )

    records = []
    for row in pos_rows:
        extra = {k: v for k, v in row.items() if k in extras}
        for a, b in pos_examples[unordered(row["protein1"], row["protein2"])]:
            records.append((a, b, 1, row["protein1"], row["protein2"], extra))
    n_pos_examples = len(records)

    from_cand, extended, per_neg = 0, 0, []
    for row in neg_rows:
        fam1, fam2 = row["protein1"], row["protein2"]
        extra = {k: v for k, v in row.items() if k in extras}
        reused = cand_examples.get(unordered(fam1, fam2))
        if reused:
            picked = sorted(reused)
            from_cand += 1
        else:
            picked = pick_pairs(fam1, fam2, members, from_universe, parent_of, n, label, args.seed)
            if len(picked) < n and has_reserve:
                wider = pick_pairs(fam1, fam2, members, from_extended, parent_of, n, label, args.seed)
                if len(wider) > len(picked):
                    picked = wider
                    extended += 1
        per_neg.append(len(picked))
        for a, b in picked:
            records.append((a, b, 0, fam1, fam2, extra))
    n_neg_examples = len(records) - n_pos_examples

    write_expanded(records, extras, args.output)

    target_ratio = len(neg_rows) / len(pos_rows) if pos_rows else 0.0
    achieved_ratio = n_neg_examples / n_pos_examples if n_pos_examples else 0.0
    st = {
        "pos_pairs": len(pos_rows),
        "pos_examples": n_pos_examples,
        "neg_pairs": len(neg_rows),
        "neg_examples": n_neg_examples,
        "neg_full": sum(1 for k in per_neg if k >= n),
        "neg_partial": sum(1 for k in per_neg if 0 < k < n),
        "neg_dropped": sum(1 for k in per_neg if k == 0),
        "from_cand": from_cand,
        "extended": extended,
        "target_ratio": target_ratio,
        "achieved_ratio": achieved_ratio,
    }
    print(
        f"  {n_pos_examples:,} positive and {n_neg_examples:,} negative examples "
        f"({st['neg_full']:,} negative pairs at full N, {st['neg_partial']:,} partial, "
        f"{st['neg_dropped']:,} dropped; {from_cand:,} reused from candidate_network, "
        f"{extended:,} needed the unclaimed reserve)",
        file=sys.stderr,
    )
    if target_ratio and abs(achieved_ratio - target_ratio) > 0.1 * target_ratio:
        print(
            f"Warning: {label} is balanced {target_ratio:.3f} negatives per positive at family level "
            f"but {achieved_ratio:.3f} at example level; raise ddi_examples_pool_factor to widen the "
            f"split's protein universe.",
            file=sys.stderr,
        )

    write_mqc(label, st, args.id)


if __name__ == "__main__":
    main()
