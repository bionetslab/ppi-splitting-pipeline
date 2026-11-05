#!/usr/bin/env python3
"""Shared I/O utilities for PPI pipeline scripts."""

import csv
import heapq
import sys
from collections import defaultdict

import numpy as np

# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------


def read_fasta(path):
    """Return {protein_id: sequence} from a FASTA file."""
    seqs = {}
    acc, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(parts)
                acc, parts = line[1:].split()[0], []
            elif line:
                parts.append(line)
    if acc and parts:
        seqs[acc] = "".join(parts)
    return seqs


def write_fasta(seqs, proteins, path):
    """Write proteins from seqs to a FASTA file, sorted by ID."""
    with open(path, "w") as fh:
        for p in sorted(proteins):
            if p in seqs:
                fh.write(f">{p}\n{seqs[p]}\n")


# ---------------------------------------------------------------------------
# PPI CSV I/O
# ---------------------------------------------------------------------------


def read_ppis(path):
    """Return list of row dicts from a PPI CSV, stripping protein ID whitespace."""
    seen, rows = set(), []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            row["protein1"] = row["protein1"].strip()
            row["protein2"] = row["protein2"].strip()
            p1, p2 = row["protein1"], row["protein2"]
            key = (min(p1, p2), max(p1, p2))
            if key not in seen:
                seen.add(key)
                rows.append(row)
    print(f"  {len(rows):,} unique PPIs", file=sys.stderr)
    return rows


def write_ppi_csv(rows, path):
    """Write PPI row dicts to CSV, preserving all columns from the input."""
    fieldnames = list(rows[0].keys()) if rows else ["protein1", "protein2"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# Domain instance table I/O (DDI mode)
# ---------------------------------------------------------------------------

INSTANCE_COLUMNS = ["instance_id", "family", "clan", "protein_id", "start", "end", "taxon_id", "source_db"]


def read_instances(path):
    """Return list of row dicts from an instances.tsv."""
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def instances_by_family(rows):
    """Invert read_instances() output to {family: {instance_id}}."""
    members = {}
    for row in rows:
        members.setdefault(row["family"], set()).add(row["instance_id"])
    return members


def instance_parent(instance_id):
    """Return the parent protein of an instance id ({family}_{protein}_{start}_{end})."""
    return instance_id.rsplit("_", 3)[1]


EMPTY = frozenset()


def pair_candidates(fam1, fam2, members, available):
    """Return every instance pair of (fam1, fam2) within `available`; self-interactions include (i, i)."""
    # A family can be missing from instances.tsv (e.g. a dead Pfam accession).
    la = sorted(members.get(fam1, EMPTY) & available)
    if fam1 == fam2:
        return [(a, b) for i, a in enumerate(la) for b in la[i:]]
    lb = sorted(members.get(fam2, EMPTY) & available)
    return [(a, b) for a in la for b in lb]


def diverse_pick(pairs, k, parent_of, rng):
    """Pick up to k pairs, preferring pairs whose parent proteins are least used so far.

    Each pick takes the lowest (used[parent_a] + used[parent_b], shuffled index), kept in a lazy heap.
    """
    order = sorted(pairs)
    rng.shuffle(order)

    touching = defaultdict(list)
    for i, (a, b) in enumerate(order):
        pa, pb = parent_of[a], parent_of[b]
        touching[pa].append(i)
        if pb != pa:
            touching[pb].append(i)

    score = [0] * len(order)
    taken = [False] * len(order)
    heap = [(0, i) for i in range(len(order))]  # already in heap order
    used, chosen = defaultdict(int), []
    while len(chosen) < k and heap:
        s, i = heapq.heappop(heap)
        if taken[i] or s != score[i]:
            continue
        taken[i] = True
        a, b = order[i]
        chosen.append(order[i])
        # A pair whose instances share a parent uses that parent twice.
        used[parent_of[a]] += 1
        used[parent_of[b]] += 1
        for p in {parent_of[a], parent_of[b]}:
            for j in touching[p]:
                if not taken[j]:
                    score[j] = used[parent_of[order[j][0]]] + used[parent_of[order[j][1]]]
                    heapq.heappush(heap, (score[j], j))
    return chosen


def expand_members(groups, members):
    """Return the union of the members of `groups`; identity when `members` is falsy (PPI mode)."""
    if not members:
        return set(groups)
    return {m for g in groups for m in members.get(g, ())}


def strict_survivors(present, keep, members):
    """Return the groups whose every member in `present` is also in `keep`.

    Plain `present & keep` when `members` is falsy (PPI mode).
    """
    if not members:
        return present & keep
    lost = present - keep
    return {g for g, ms in members.items() if (ms & present) and not (ms & lost)}


# ---------------------------------------------------------------------------
# KaHIP partition I/O
# ---------------------------------------------------------------------------


def read_node_mapping(path):
    """Return {node_id (int, 1-indexed): protein_id}."""
    mapping = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mapping[int(row["node_id"])] = row["protein_id"]
    return mapping


def read_partition(path):
    """Return list where element i is the partition of node (i+1)."""
    partitions = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                partitions.append(int(line))
    return partitions


def read_labelled_csv(path):
    """Return (pairs, labels) from a CSV with protein1, protein2, label columns."""
    pairs, labels = [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            pairs.append((row["protein1"].strip(), row["protein2"].strip()))
            labels.append(int(row["label"]))
    return pairs, np.array(labels)


def read_family_pairs(path):
    """Return the sorted (family1, family2) pair per row, or None if the CSV has no family columns.

    Rows are in file order, so they line up with read_labelled_csv().
    """
    with open(path) as fh:
        reader = csv.DictReader(fh)
        cols = reader.fieldnames or []
        if "family1" not in cols or "family2" not in cols:
            return None
        pairs = []
        for row in reader:
            f1, f2 = row["family1"].strip(), row["family2"].strip()
            pairs.append((min(f1, f2), max(f1, f2)))
    return pairs


# ---------------------------------------------------------------------------
# Embeddings
# ---------------------------------------------------------------------------


def load_embeddings(path):
    """Return {protein_id: embedding_array} from an NPZ file."""
    raw = np.load(path, allow_pickle=False)
    return {k: raw[k] for k in raw.files}


# ---------------------------------------------------------------------------
# MultiQC sample naming
# ---------------------------------------------------------------------------

# MultiQC sorts rows alphabetically by Sample name, so a plain f"{id}_{split}"
# would show test_balanced, test_realistic, train, val. This numeric prefix
# forces train/val/test/test_balanced/test_realistic order per dataset.
SPLIT_ORDER = ["train", "val", "test", "test_balanced", "test_realistic", "discarded"]


def mqc_sample(label, split):
    """Build a MultiQC Sample name for one dataset+split that sorts correctly."""
    return f"{label}_{SPLIT_ORDER.index(split) + 1}_{split}"


def mqc_dataset(label):
    """Build a MultiQC Sample name for a row that covers a whole dataset."""
    return label


def mqc_category(name):
    """Same sort-key prefix as mqc_sample, for bar-chart categories that are
    already scoped to one dataset's own chart (no id_ qualification needed)."""
    return f"{SPLIT_ORDER.index(name) + 1}_{name}"
