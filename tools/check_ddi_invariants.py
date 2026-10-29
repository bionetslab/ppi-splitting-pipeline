#!/usr/bin/env python3
"""Check the invariants a DDI-mode run must satisfy, over a published results tree.

    python tools/check_ddi_invariants.py --results results/results_test_ddi

--results is a single dataset directory or the whole outdir (`_shared` is skipped).
Per dataset and negative set it checks that no family, clan or parent protein is in
two splits, every instance id matches its row's family and its split's universe or
reserve, no (family pair, label) has more than N examples, no instance pair carries
both labels, and the universes are disjoint; across negative sets, that the positives
are identical. For --leaky-ids datasets (split_method=random) the leak checks must fail.
Exits 1 on any failure.
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin"))

from utils import read_instances

LABELS = ["train", "val", "test_balanced", "test_realistic"]
SUFFIX = "_instances.csv"


def split_of(label):
    """Return the SELECT_EXAMPLES split a sampling label draws its universe from."""
    return "test" if label.startswith("test") else label


def discover_negsets(path):
    """Return {negset: {label: csv path}} for one dataset directory.

    A single-negative-set run comes back under "". The unsuffixed test_realistic is
    attached to every set. Labels are matched longest-first.
    """
    by_negset, shared_realistic = defaultdict(dict), None
    longest_first = sorted(LABELS, key=len, reverse=True)
    for name in sorted(os.listdir(path)):
        if not name.endswith(SUFFIX):
            continue
        stem = name[: -len(SUFFIX)]
        label = next((lbl for lbl in longest_first if stem == lbl or stem.startswith(lbl + "_")), None)
        if label is None:
            continue
        negset = "" if stem == label else stem[len(label) + 1 :]
        if label == "test_realistic" and not negset:
            shared_realistic = os.path.join(path, name)
        by_negset[negset][label] = os.path.join(path, name)
    if shared_realistic is not None:
        for files in by_negset.values():
            files.setdefault("test_realistic", shared_realistic)
        # With named sets present, "" only holds the shared test_realistic.
        if len(by_negset) > 1 and set(by_negset.get("", {})) == {"test_realistic"}:
            del by_negset[""]
    return dict(by_negset)


def find_datasets(root):
    """Return [(dataset_id, path)], accepting either one dataset dir or a whole outdir."""
    if discover_negsets(root):
        return [(os.path.basename(os.path.abspath(root)), root)]
    found = []
    for name in sorted(os.listdir(root)):
        path = os.path.join(root, name)
        if name == "_shared" or not os.path.isdir(path):
            continue
        if discover_negsets(path):
            found.append((name, path))
    return found


def positive_pairs(rows):
    """Return the positive instance pairs of one labelled table as a set."""
    return {
        (row["family1"], row["family2"], tuple(sorted((row["protein1"], row["protein2"]))))
        for row in rows
        if row.get("label") == "1"
    }


def find_instances(dataset_dir, explicit):
    if explicit:
        return explicit
    candidates = [
        os.path.join(dataset_dir, "data", "instances.tsv"),
        os.path.join(os.path.dirname(os.path.abspath(dataset_dir)), "_shared", "data", "instances.tsv"),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    return None


def read_rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def read_ids(path):
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


class Report:
    """Collects check outcomes so one dataset reports all its failures, not just the first."""

    def __init__(self, dataset):
        self.dataset = dataset
        self.failures = 0

    def check(self, ok, name, detail=""):
        status = "PASS" if ok else "FAIL"
        if not ok:
            self.failures += 1
        print(f"  [{status}] {name}{(' -- ' + detail) if detail else ''}")

    def note(self, text):
        print(f"  [note] {text}")


def check_dataset(dataset, path, instances_path, target, leaky, negset, files):
    """Check one negative set of one dataset; returns (failures, {label: rows})."""
    shown = f"{dataset} [{negset}]" if negset else dataset
    print(f"\n== {shown}{'  (leaky baseline: leak-shaped checks are inverted)' if leaky else ''}")
    rep = Report(dataset)

    rows = read_instances(instances_path)
    inst_family = {r["instance_id"]: r["family"] for r in rows}
    inst_parent = {r["instance_id"]: r["protein_id"] for r in rows}
    fam_clan = {r["family"]: r["clan"] for r in rows}

    labelled, missing = {}, []
    for lbl in LABELS:
        f = files.get(lbl)
        if f and os.path.exists(f):
            labelled[lbl] = read_rows(f)
        else:
            missing.append(lbl)
    if missing:
        rep.check(False, "all four labelled instance tables present", f"missing {', '.join(missing)}")
        return rep.failures, labelled
    rep.check(True, "all four labelled instance tables present")

    unclaimed = read_ids(os.path.join(path, "examples", "unclaimed.txt"))
    # Each split may only draw from its own share of the reserve.
    universes, reserves = {}, {}
    for split in ("train", "val", "test"):
        universes[split] = read_ids(os.path.join(path, "examples", f"{split}_universe.txt"))
        reserves[split] = read_ids(os.path.join(path, "examples", f"{split}_reserve.txt"))

    bad_family, bad_parent, over, both = [], [], [], []
    fam_by_split, clan_by_split, parents_by_split = defaultdict(set), defaultdict(set), defaultdict(set)
    all_pairs = set()

    for lbl, table in labelled.items():
        split = split_of(lbl)
        # Per table: test_balanced and test_realistic share the test positives.
        per_pair, pairs_by_label = defaultdict(int), defaultdict(set)
        for row in table:
            if not {"protein1", "protein2", "label", "family1", "family2"} <= set(row):
                rep.check(False, f"{lbl}: instance table carries the DDI columns", "expected family1/family2/label")
                return rep.failures, labelled
            for col, fam_col in (("protein1", "family1"), ("protein2", "family2")):
                iid = row[col]
                if inst_family.get(iid) != row[fam_col]:
                    bad_family.append(f"{lbl}:{iid}")
                    continue
                parent = inst_parent[iid]
                parents_by_split[split].add(parent)
                if parent not in universes[split] and parent not in reserves[split]:
                    bad_parent.append(f"{lbl}:{parent}")
                fam_by_split[split].add(row[fam_col])
                clan_by_split[split].add(fam_clan[row[fam_col]])
            per_pair[(row["family1"], row["family2"], row["label"])] += 1
            pairs_by_label[row["label"]].add(tuple(sorted((row["protein1"], row["protein2"]))))
        over += [f"{lbl}:{f1}-{f2}(label {y})={n}" for (f1, f2, y), n in per_pair.items() if n > target]
        both += [f"{lbl}:{p}" for p in pairs_by_label.get("1", set()) & pairs_by_label.get("0", set())]
        all_pairs |= {(f1, f2) for f1, f2, _ in per_pair}

    rep.check(not bad_family, "every id is an instance of the family its row declares", f"{len(bad_family)} bad")
    rep.check(
        not bad_parent,
        "every parent is in its own split's universe or its own share of the reserve",
        f"{len(bad_parent)} outside",
    )
    rep.check(not over, f"no (family pair, label) carries more than N={target} examples", ", ".join(over[:4]))
    rep.check(not both, "no instance pair carries both labels", f"{len(both)} do")

    counts = ", ".join(f"{lbl}={len(labelled[lbl])}" for lbl in LABELS)
    rep.note(f"rows: {counts}; distinct family pairs: {len(all_pairs)}")

    # Leak-shaped checks; inverted for a leaky (split_method=random) dataset.
    shared = {}
    for a, b in (("train", "val"), ("train", "test"), ("val", "test")):
        shared[("family", a, b)] = fam_by_split[a] & fam_by_split[b]
        shared[("clan", a, b)] = clan_by_split[a] & clan_by_split[b]
        shared[("parent", a, b)] = parents_by_split[a] & parents_by_split[b]
        shared[("universe", a, b)] = universes[a] & universes[b]

    for kind in ("family", "clan", "parent", "universe"):
        overlaps = {k: v for k, v in shared.items() if k[0] == kind}
        total = sum(len(v) for v in overlaps.values())
        detail = ", ".join(f"{a}/{b}: {len(v)}" for (_, a, b), v in sorted(overlaps.items()))
        name = {
            "family": "no Pfam family in two splits",
            "clan": "no clan in two splits",
            "parent": "no parent protein in two splits",
            "universe": "protein universes pairwise disjoint",
        }[kind]
        if leaky:
            rep.check(total > 0, f"leaky baseline still leaks ({kind})", detail)
        else:
            rep.check(total == 0, name, detail)

    if not leaky:
        bleed = {s: len(u & unclaimed) for s, u in universes.items()}
        rep.check(
            not any(bleed.values()),
            "unclaimed.txt disjoint from every universe",
            ", ".join(f"{s}: {n}" for s, n in sorted(bleed.items())),
        )
        rep.note(f"{len(unclaimed)} proteins never claimed by any split")

    return rep.failures, labelled


def check_shared_positives(dataset, tables):
    """Check that every negative set carries the same positive rows."""
    if len(tables) < 2:
        return 0
    print(f"\n== {dataset}  (across {len(tables)} negative sets)")
    rep = Report(dataset)
    ref_negset = sorted(tables)[0]
    for lbl in LABELS:
        ref = positive_pairs(tables[ref_negset].get(lbl, []))
        differing = []
        for ns, rows in sorted(tables.items()):
            if ns == ref_negset:
                continue
            got = positive_pairs(rows.get(lbl, []))
            if got != ref:
                differing.append(f"{ns}: {len(got ^ ref)} differ")
        rep.check(
            not differing, f"{lbl}: same positives in every negative set ({len(ref)} pairs)", "; ".join(differing)
        )
    return rep.failures


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, help="a dataset directory, or the whole outdir")
    ap.add_argument("--instances", help="instances.tsv; found under the results tree when omitted")
    ap.add_argument("--examples-target", type=int, default=5, help="N, the per-DDI example cap (default 5)")
    ap.add_argument(
        "--leaky-ids",
        default=r"(?i)leak|random",
        help="regex for dataset ids run with split_method=random, whose leak-shaped invariants are inverted",
    )
    args = ap.parse_args()

    datasets = find_datasets(args.results)
    if not datasets:
        sys.exit(f"No dataset with *_instances.csv found under {args.results} -- is this a DDI-mode run?")

    leaky_re = re.compile(args.leaky_ids)
    failures = 0
    for dataset, path in datasets:
        instances_path = find_instances(path, args.instances)
        if not instances_path:
            print(f"\n== {dataset}\n  [FAIL] instances.tsv not found -- pass --instances")
            failures += 1
            continue
        leaky = bool(leaky_re.search(dataset))
        negsets = discover_negsets(path)
        tables = {}
        for negset in sorted(negsets):
            n_failed, labelled = check_dataset(
                dataset, path, instances_path, args.examples_target, leaky, negset, negsets[negset]
            )
            failures += n_failed
            tables[negset] = labelled
        failures += check_shared_positives(dataset, tables)

    print(f"\n{len(datasets)} dataset(s) checked, {failures} failed check(s)")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
