#!/usr/bin/env python3
"""Shared I/O utilities for PPI pipeline scripts."""

import csv

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
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            row["protein1"] = row["protein1"].strip()
            row["protein2"] = row["protein2"].strip()
            rows.append(row)
    return rows


def write_ppi_csv(rows, path):
    """Write PPI row dicts to CSV, preserving all columns from the input."""
    fieldnames = list(rows[0].keys()) if rows else ["protein1", "protein2"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_labelled_csv(path):
    """Return (pairs, labels) from a CSV with protein1, protein2, label columns."""
    pairs, labels = [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            pairs.append((row["protein1"].strip(), row["protein2"].strip()))
            labels.append(int(row["label"]))
    return pairs, np.array(labels)


# ---------------------------------------------------------------------------
# Embeddings
# ---------------------------------------------------------------------------

def load_embeddings(path):
    """Return {protein_id: embedding_array} from an NPZ file."""
    raw = np.load(path, allow_pickle=False)
    return {k: raw[k] for k in raw.files}
