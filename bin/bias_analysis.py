#!/usr/bin/env python3
"""
Bias analysis: utility and detectability for PPI dataset attributes.

Attributes
----------
sequence_similarity      – BLAST pident between the two proteins in the pair,
                           normalised to [0, 1]
embedding_similarity     – cosine similarity between the two individual protein
                           embeddings in the pair
functional_relatedness_BP/MF/CC – Jaccard similarity of GO Biological Process,
                                  Molecular Function, and Cellular Component term
                                  sets between the two proteins in a pair

Utility       = MI(A; Y) via sklearn kNN estimator (continuous A)
Detectability = 5-fold CV Spearman ρ of Random Forest Regressor predicting A
                from the pair embedding X

Outputs
-------
bias_table_{split}_mqc.tsv    – MultiQC table per split (MI, relationship flag, detectability)
bias_scatter_{split}_mqc.json – MultiQC scatter per split (x = detectability, y = utility)
"""

import argparse
import csv
import json
import sys
from collections import defaultdict

import numpy as np
import plotly.graph_objects as go
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import mutual_info_classif
from sklearn.metrics import make_scorer
from sklearn.model_selection import cross_val_score


def read_labelled_csv(path):
    pairs, labels = [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            pairs.append((row["protein1"].strip(), row["protein2"].strip()))
            labels.append(int(row["label"]))
    return pairs, np.array(labels)


def load_embeddings(path):
    raw = np.load(path, allow_pickle=False)
    return {k: raw[k] for k in raw.files}


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


_spearman_scorer = make_scorer(lambda y_true, y_pred: spearmanr(y_true, y_pred)[0])

# ColorBrewer Set2 qualitative palette (8 colours)
_SET2 = ["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
         "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"]


def analyse(A, X, y, name, seed=42):
    """Return dict with mi, related, detectability (5-fold CV Spearman ρ)."""
    if name == "self_interactions":
        discrete_features = True
    else:
        discrete_features = False
    mi = float(
        mutual_info_classif(A.reshape(-1, 1), y, discrete_features=discrete_features, random_state=seed)[0]
    )

    if X.shape[0] < 10:
        detectability = 0.0
    else:
        rf = RandomForestRegressor(n_estimators=100, random_state=seed, n_jobs=-1)
        A_jit = A + np.random.default_rng(seed).random(len(A)).astype(np.float32) * 1e-6
        detectability = float(cross_val_score(rf, X, A_jit, cv=5, scoring=_spearman_scorer).mean())

    return {"mi": mi, "related": mi > 0, "detectability": detectability}


def _write_tables_mqc(all_results):
    """MultiQC JSON table with one dataset tab per split.

    Passing data as a list + data_labels in pconfig makes MultiQC render
    dataset-switch tabs, the same mechanism used for bar/line plots.
    The native table engine keeps sorting and per-column colour scales.
    """
    splits = list(all_results.keys())

    data_list = [
        {name: {"MI(A;Y)": r["mi"],
                "Related?": "Yes" if r["related"] else "No",
                "Detectability (Spearman ρ)": r["detectability"]}
         for name, r in all_results[split]}
        for split in splits
    ]

    doc = {
        "id": "bias_tables",
        "section_name": "Bias Analysis Tables",
        "description": "Utility MI(A;Y) and detectability (5-fold RF Spearman ρ) per split.",
        "plot_type": "table",
        "pconfig": {
            "id": "bias_tables_plot",
            "title": "Attribute Utility and Detectability",
            "data_labels": [{"name": s} for s in splits],
        },
        "headers": {
            "MI(A;Y)": {"title": "MI(A;Y)", "format": "{:.4f}", "scale": "Greens"},
            "Related?": {"title": "Related?", "scale": False},
            "Detectability (Spearman ρ)": {
                "title": "Detectability (Spearman ρ)",
                "format": "{:.4f}",
                "scale": "Blues",
            },
        },
        "data": data_list,
    }

    with open("bias_tables_mqc.json", "w") as fh:
        json.dump(doc, fh, indent=2)


_SPLIT_SYMBOLS = {
    "train":          "circle",
    "val":            "square",
    "test_balanced":  "diamond",
    "test_realistic": "triangle-up",
}


def _write_scatter_mqc(all_results):
    """Interactive Plotly scatter: colour = attribute, marker shape = split.

    One trace per attribute with one point per split; marker.symbol is a list
    so each point gets the correct split shape.  Standard Plotly legend
    behaviour (click to hide, double-click to isolate) works on attributes.
    Split dummy traces provide a shape key but are purely informational.
    """
    first = next(iter(all_results.values()))
    attr_names = [name for name, _ in first]
    splits = list(all_results.keys())
    color_map = {name: _SET2[i % len(_SET2)] for i, name in enumerate(attr_names)}

    fig = go.Figure()

    for attr in attr_names:
        xs, ys, symbols, hover_splits = [], [], [], []
        for split in splits:
            r = dict(all_results[split])[attr]
            xs.append(r["detectability"])
            ys.append(r["mi"])
            symbols.append(_SPLIT_SYMBOLS[split])
            hover_splits.append(split)
        fig.add_trace(go.Scatter(
            x=xs, y=ys,
            mode="markers",
            name=attr,
            customdata=hover_splits,
            marker=dict(
                color=color_map[attr],
                symbol=symbols,
                size=10,
                line=dict(color="white", width=1),
            ),
            hovertemplate=(
                f"<b>{attr}</b><br>"
                "Split: %{customdata}<br>"
                "Detectability: %{x:.3f}<br>"
                "MI: %{y:.4f}<extra></extra>"
            ),
        ))

    # Informational split shape key (dummies; clicking only hides the key marker).
    for j, (split, symbol) in enumerate(_SPLIT_SYMBOLS.items()):
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode="markers",
            name=split,
            legendgroup="_splits",
            legendgrouptitle=dict(text="Split") if j == 0 else None,
            showlegend=True,
            marker=dict(color="gray", symbol=symbol, size=10),
        ))

    fig.update_layout(
        xaxis_title="Detectability (RF Spearman ρ)",
        yaxis_title="Utility MI(A;Y)",
        legend=dict(title=dict(text="Attribute"), groupclick="toggleitem"),
        margin=dict(l=60, r=20, t=30, b=50),
        height=420,
        autosize=True,
    )

    div = fig.to_html(full_html=False, include_plotlyjs="cdn", config={"responsive": True})

    with open("bias_scatter_mqc.html", "w") as fh:
        fh.write(
            "<!--\n"
            "id: 'bias_scatter'\n"
            "section_name: 'Bias Analysis – Utility vs. Detectability'\n"
            "description: 'Colour = attribute (click/double-click to filter), "
            "shape = split. "
            "x = detectability (RF Spearman ρ), y = utility MI(A;Y).'\n"
            "-->\n"
        )
        fh.write(div)


def write_mqc(all_results):
    _write_tables_mqc(all_results)
    _write_scatter_mqc(all_results)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train",           required=True, help="train.csv with labels")
    ap.add_argument("--val",             required=True, help="val.csv with labels")
    ap.add_argument("--test_balanced",   required=True, help="test_balanced.csv with labels")
    ap.add_argument("--test_realistic",  required=True, help="test_realistic.csv with labels")
    ap.add_argument("--blast",           required=True, help="all_vs_all.tsv with pident at col 4")
    ap.add_argument("--embeddings",     required=True, help="embeddings.npz")
    ap.add_argument("--go_annotations", required=True, help="go_annotations.tsv with GO term sets")
    ap.add_argument("--seed",           type=int, default=42)
    args = ap.parse_args()

    print("Loading embeddings ...", file=sys.stderr)
    embeddings = load_embeddings(args.embeddings)
    print("Parsing BLAST similarities ...", file=sys.stderr)
    blast_sim = parse_blast_pident(args.blast)
    go_anns = load_go_annotations(args.go_annotations)

    split_paths = [
        ("train",          args.train),
        ("val",            args.val),
        ("test_balanced",  args.test_balanced),
        ("test_realistic", args.test_realistic),
    ]

    all_results = {}
    for split, path in split_paths:
        pairs, y = read_labelled_csv(path)
        X, pairs_f, y_f = _prepare_split(pairs, y, embeddings)
        print(f"  {X.shape[0]} {split} pairs retained for analysis", file=sys.stderr)

        attrs = [
            ("sequence_similarity",        seq_sim_within_pair(pairs_f, blast_sim)),
            ("embedding_similarity",       emb_sim_within_pair(pairs_f, embeddings)),
            ("functional_relatedness_BP",  func_relatedness(pairs_f, go_anns, "BP")),
            ("functional_relatedness_MF",  func_relatedness(pairs_f, go_anns, "MF")),
            ("functional_relatedness_CC",  func_relatedness(pairs_f, go_anns, "CC")),
            ("self_interactions",          self_interactions(pairs_f)),
        ]

        results = []
        for name, A in attrs:
            r = analyse(A, X, y_f, name, seed=args.seed)
            results.append((name, r))
            print(f"  [{split}] {name}: MI={r['mi']:.4f}  related={r['related']}  ρ={r['detectability']:.4f}", file=sys.stderr)
        all_results[split] = results

    write_mqc(all_results)


if __name__ == "__main__":
    main()
