#!/usr/bin/env python3
"""
Bias analysis: utility and detectability for one PPI dataset attribute.

Run once per attribute (--attribute flag); Nextflow calls this in parallel for
all attributes.

Attributes
----------
sequence_similarity           – BLAST pident between the pair, normalised to [0,1]
embedding_similarity          – cosine similarity of the two protein embeddings
functional_relatedness_BP/MF/CC – Jaccard similarity of GO term sets
self_interactions             – 1 if both proteins are identical, 0 otherwise

Utility       = NMI(A; Y) = MI / sqrt(H(A)·H(Y)), MI estimated by sklearn kNN, H(A) by histogram
Detectability = Spearman ρ of a Ridge regressor predicting A from the pair embedding X

Output
------
{attribute}_bias_mqc.tsv  – MultiQC table with splits as rows
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.feature_selection import mutual_info_classif

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import load_embeddings, read_labelled_csv




def build_pair_X(pairs, embeddings):
    """Return (X, mask) where mask[i]=True when both proteins have embeddings."""
    rows, mask = [], []
    for p1, p2 in pairs:
        a, b = (p1, p2) if p1 <= p2 else (p2, p1)
        if a in embeddings and b in embeddings:
            rows.append(np.concatenate([embeddings[a], embeddings[b]]))
            mask.append(True)
        else:
            mask.append(False)
    X = np.array(rows, dtype=np.float32) if rows else np.empty((0, 0), dtype=np.float32)
    return X, np.array(mask, dtype=bool)


def parse_blast_pident(path):
    """Return {prot: {other_prot: max_pident}} from BLAST outfmt-6 with pident at col 4."""
    sim = defaultdict(dict)
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip().split("\t")
            if len(cols) < 5:
                continue
            q = cols[0]
            s = cols[1].split("|")[1] if "|" in cols[1] else cols[1]
            if q == s:
                continue
            try:
                pident = float(cols[4])
            except ValueError:
                continue
            if pident > sim[q].get(s, -1):
                sim[q][s] = pident
            if pident > sim[s].get(q, -1):
                sim[s][q] = pident
    return sim


def seq_sim_within_pair(pairs, blast_sim):
    """BLAST pident between the two proteins in each pair, normalised to [0, 1]."""
    A = []
    for p1, p2 in pairs:
        pident = blast_sim.get(p1, {}).get(p2, 0.0)
        A.append(pident / 100.0)
    return np.array(A, dtype=np.float32)


def emb_sim_within_pair(pairs, embeddings):
    """Cosine similarity between the two individual protein embeddings in each pair."""
    sims = []
    for p1, p2 in pairs:
        e1, e2 = embeddings[p1], embeddings[p2]
        n1, n2 = np.linalg.norm(e1), np.linalg.norm(e2)
        sims.append(float(np.dot(e1, e2) / (n1 * n2)) if n1 > 0 and n2 > 0 else 0.0)
    return np.array(sims, dtype=np.float32)


def load_go_annotations(path):
    """Return {protein_id: {BP/MF/CC: frozenset[go_term]}} from the four-column TSV."""
    result = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            acc = row["protein_id"].strip()
            result[acc] = {
                cat: frozenset(t.strip() for t in row[col].strip().split(";") if t.strip())
                for cat, col in [("BP", "go_bp"), ("MF", "go_mf"), ("CC", "go_cc")]
            }
    return result


def self_interactions(pairs):
    """1 if both proteins in the pair are identical, 0 otherwise."""
    return np.array([1.0 if p1 == p2 else 0.0 for p1, p2 in pairs], dtype=np.float32)


def _prepare_split(pairs, y, embeddings):
    """Return (X, filtered_pairs, filtered_y) keeping only pairs with embeddings."""
    X, mask = build_pair_X(pairs, embeddings)
    return X, [p for p, ok in zip(pairs, mask) if ok], y[mask]


def func_relatedness(pairs, go_anns, category):
    """Jaccard similarity of GO term sets for one category per pair; 0.0 if union is empty."""
    empty = frozenset()
    sims = []
    for p1, p2 in pairs:
        s1 = go_anns.get(p1, {}).get(category, empty)
        s2 = go_anns.get(p2, {}).get(category, empty)
        union_size = len(s1 | s2)
        sims.append(len(s1 & s2) / union_size if union_size > 0 else 0.0)
    return np.array(sims, dtype=np.float32)



def _entropy(counts):
    """Shannon entropy in nats; ignores zero-count bins."""
    c = counts[counts > 0].astype(float)
    p = c / c.sum()
    return float(-np.sum(p * np.log(p)))


def analyse(A, X, y, name, seed=42, n_bins=10):
    """Return dict with nmi, related, detectability (train-set Spearman ρ).

    NMI = MI(A; Y) / sqrt(H(A) · H(Y)) — symmetric, bounded in [0, 1].
    H(A) is estimated by binning A into n_bins equal-width histogram bins
    (or using unique values directly when A is already discrete).
    """
    discrete_features = name == "self_interactions"
    mi = float(
        mutual_info_classif(A.reshape(-1, 1), y, discrete_features=discrete_features, random_state=seed)[0]
    )

    _, y_counts = np.unique(y, return_counts=True)
    h_y = _entropy(y_counts)

    if discrete_features:
        _, a_counts = np.unique(A, return_counts=True)
    else:
        a_counts, _ = np.histogram(A, bins=n_bins)
    h_a = _entropy(a_counts)

    denom = np.sqrt(h_a * h_y)
    nmi = float(np.clip(mi / denom, 0.0, 1.0)) if denom > 0 else 0.0
    print(f"    NMI(A;Y) = {nmi:.4f}  (related? {'Yes' if nmi > 0 else 'No'})", file=sys.stderr)

    if X.shape[0] < 10:
        detectability = 0.0
    else:
        reg = Ridge(alpha=1.0, solver="sag", random_state=seed, max_iter=1000)
        reg.fit(X, A)
        detectability = float(spearmanr(A, reg.predict(X))[0])

    return {"nmi": nmi, "related": nmi > 0, "detectability": detectability}


ATTRIBUTES = {
    "sequence_similarity":       lambda pairs, blast_sim, emb, go: seq_sim_within_pair(pairs, blast_sim),
    "embedding_similarity":      lambda pairs, blast_sim, emb, go: emb_sim_within_pair(pairs, emb),
    "functional_relatedness_BP": lambda pairs, blast_sim, emb, go: func_relatedness(pairs, go, "BP"),
    "functional_relatedness_MF": lambda pairs, blast_sim, emb, go: func_relatedness(pairs, go, "MF"),
    "functional_relatedness_CC": lambda pairs, blast_sim, emb, go: func_relatedness(pairs, go, "CC"),
    "self_interactions":         lambda pairs, blast_sim, emb, go: self_interactions(pairs),
}


def write_mqc(attribute, results):
    fname = f"{attribute}_bias_mqc.tsv"
    with open(fname, "w") as fh:
        fh.write(
            f"# id: 'bias_{attribute}'\n"
            f"# section_name: 'Bias: {attribute}'\n"
            f"# description: 'Utility NMI(A;Y) = MI/sqrt(H(A)·H(Y)) and detectability (Ridge Spearman ρ) for {attribute}.'\n"
            "# plot_type: 'table'\n"
            "# pconfig:\n"
            f"#     id: 'bias_{attribute}_table'\n"
            f"#     title: 'Bias: {attribute}'\n"
            "Split\tNMI(A;Y)\tRelated?\tDetectability (Spearman ρ)\n"
        )
        for split, r in results:
            fh.write(f"{split}\t{r['nmi']:.4f}\t{'Yes' if r['related'] else 'No'}\t{r['detectability']:.4f}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--attribute",      required=True,
                    choices=list(ATTRIBUTES),
                    help="attribute to analyse")
    ap.add_argument("--train",           required=True)
    ap.add_argument("--val",             required=True)
    ap.add_argument("--test_balanced",   required=True)
    ap.add_argument("--test_realistic",  required=True)
    ap.add_argument("--blast",           required=True)
    ap.add_argument("--embeddings",      required=True)
    ap.add_argument("--go_annotations",  required=True)
    ap.add_argument("--seed",            type=int, default=42)
    args = ap.parse_args()

    print(f"=== {args.attribute} ===", file=sys.stderr)

    compute = ATTRIBUTES[args.attribute]

    # Load data ---------------------------------------------------------------
    print("Loading embeddings ...", file=sys.stderr)
    embeddings = load_embeddings(args.embeddings)

    blast_sim = None
    if args.attribute == "sequence_similarity":
        print("Parsing BLAST similarities ...", file=sys.stderr)
        blast_sim = parse_blast_pident(args.blast)

    go_anns = None
    if args.attribute.startswith("functional_relatedness"):
        go_anns = load_go_annotations(args.go_annotations)

    split_paths = [
        ("train",          args.train),
        ("val",            args.val),
        ("test_balanced",  args.test_balanced),
        ("test_realistic", args.test_realistic),
    ]

    results = []
    for split, path in split_paths:
        pairs, y = read_labelled_csv(path)
        X, pairs_f, y_f = _prepare_split(pairs, y, embeddings)
        print(f"  {X.shape[0]} {split} pairs retained", file=sys.stderr)
        A = compute(pairs_f, blast_sim, embeddings, go_anns)
        r = analyse(A, X, y_f, args.attribute, seed=args.seed)
        results.append((split, r))
        print(f"  [{split}] NMI={r['nmi']:.4f}  ρ={r['detectability']:.4f}", file=sys.stderr)

    write_mqc(args.attribute, results)


if __name__ == "__main__":
    main()
