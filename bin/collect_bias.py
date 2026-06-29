#!/usr/bin/env python3
"""
Collect per-attribute bias TSVs and produce a combined Plotly scatter.

Reads all *_bias_mqc.tsv files written by bias_analysis.py and emits
bias_scatter_mqc.html for MultiQC.

Each TSV has the header:
    Split  MI(A;Y)  Related?  Detectability (Spearman ρ)
The attribute name is derived from the filename: {attribute}_bias_mqc.tsv
"""

import argparse
import csv
import sys
from pathlib import Path

import plotly.graph_objects as go

_SET2 = ["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
         "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"]

_SPLIT_SYMBOLS = {
    "train":          "circle",
    "val":            "square",
    "test_balanced":  "diamond",
    "test_realistic": "triangle-up",
}


def parse_tsv(path):
    """Return list of (split, mi, detectability) from a bias TSV."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "Split":
                continue
            split, mi, _, det = parts[0], float(parts[1]), parts[2], float(parts[3])
            rows.append((split, mi, det))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsvs", nargs="+", help="*_bias_mqc.tsv files from bias_analysis.py")
    args = ap.parse_args()

    # Parse all files; derive attribute name from filename
    data = {}
    for p in sorted(args.tsvs):
        stem = Path(p).name  # e.g. sequence_similarity_bias_mqc.tsv
        attribute = stem.replace("_bias_mqc.tsv", "")
        rows = parse_tsv(p)
        if rows:
            data[attribute] = rows
        else:
            print(f"  [warn] no data in {p}", file=sys.stderr)

    if not data:
        print("No data found; skipping scatter.", file=sys.stderr)
        sys.exit(0)

    color_map = {attr: _SET2[i % len(_SET2)] for i, attr in enumerate(sorted(data))}

    fig = go.Figure()

    for attr, rows in sorted(data.items()):
        xs, ys, symbols, hover_splits = [], [], [], []
        for split, mi, det in rows:
            xs.append(det)
            ys.append(mi)
            symbols.append(_SPLIT_SYMBOLS.get(split, "circle"))
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
        xaxis_title="Detectability (Ridge Spearman ρ)",
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
            "shape = split. x = detectability (Ridge Spearman ρ), y = utility MI(A;Y).'\n"
            "-->\n"
        )
        fh.write(div)

    print("Wrote bias_scatter_mqc.html", file=sys.stderr)


if __name__ == "__main__":
    main()
