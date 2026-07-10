#!/usr/bin/env python3
"""
Clustered bitscore heatmap with train / val / test colour annotations.

Samples up to --max_per_split proteins from each split, builds a symmetric
bitscore matrix from the all-vs-all BLAST file, runs average-linkage
hierarchical clustering, and writes a MultiQC-compatible HTML file.
"""

import argparse
import sys

import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist


def read_fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def build_bitscore_matrix(blast_path, protein_ids):
    id_set = set(protein_ids)
    id_idx = {pid: i for i, pid in enumerate(protein_ids)}
    n = len(protein_ids)
    mat = np.zeros((n, n), dtype=np.float32)

    with open(blast_path) as fh:
        for line in fh:
            cols = line.rstrip().split("\t")
            if len(cols) < 4:
                continue
            q = cols[0]
            s = cols[1].split("|")[1] if "|" in cols[1] else cols[1]
            if q == s or q not in id_set or s not in id_set:
                continue
            try:
                bitscore = float(cols[3])
            except ValueError:
                continue
            i, j = id_idx[q], id_idx[s]
            if bitscore > mat[i, j]:
                mat[i, j] = bitscore
                mat[j, i] = bitscore

    return mat


def cluster_order(mat):
    """Return leaf ordering from average-linkage clustering on bitscore profiles."""
    Z = linkage(pdist(mat, metric="euclidean"), method="average", optimal_ordering=True)
    # scipy's dendrogram is recursive; with n leaves the stack can reach ~n frames.
    needed = max(sys.getrecursionlimit(), len(mat) * 10)
    sys.setrecursionlimit(needed)
    return dendrogram(Z, no_plot=True)["leaves"]


def make_html(mat, split_labels, output, id_):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    palette = {"train": "#4C72B0", "val": "#DD8452", "test": "#55A868"}
    label_to_num = {"train": 0, "val": 1, "test": 2}

    # Discrete 3-colour scale: 0→train, 1→val, 2→test
    ann_colorscale = [
        [0.000, palette["train"]], [0.333, palette["train"]],
        [0.334, palette["val"]],   [0.666, palette["val"]],
        [0.667, palette["test"]],  [1.000, palette["test"]],
    ]

    col_nums = [[label_to_num[l] for l in split_labels]]
    row_nums = [[label_to_num[l]] for l in split_labels]

    shared_ann = dict(
        colorscale=ann_colorscale,
        zmin=0, zmax=2,
        showscale=False,
        hoverinfo="skip",
    )

    fig = make_subplots(
        rows=2, cols=2,
        column_widths=[0.025, 0.975],
        row_heights=[0.025, 0.975],
        horizontal_spacing=0.003,
        vertical_spacing=0.003,
    )

    # Top annotation strip (column labels)
    fig.add_trace(go.Heatmap(z=col_nums, **shared_ann), row=1, col=2)

    # Left annotation strip (row labels)
    fig.add_trace(go.Heatmap(z=row_nums, **shared_ann), row=2, col=1)

    # Clip colorscale at 95th percentile of non-zero values so rare high-scoring
    # pairs don't compress all contrast into the bottom of the scale.
    nonzero = mat[mat > 0]
    zmax = float(np.percentile(nonzero, 95)) if nonzero.size > 0 else 1.0

    # Main heatmap
    fig.add_trace(
        go.Heatmap(
            z=mat.tolist(),
            colorscale="Purples",
            zmin=0,
            zmax=zmax,
            colorbar=dict(
                title=dict(text="Bitscore", side="right"),
                x=1.02,
                ticksuffix="",
                outlinewidth=0,
            ),
            hovertemplate="bitscore: %{z:.1f}<extra></extra>",
        ),
        row=2, col=2,
    )

    # Legend for annotation colours (dummy scatter traces)
    for split, colour in palette.items():
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None],
                mode="markers",
                name=split,
                marker=dict(color=colour, size=10, symbol="square"),
                showlegend=True,
            )
        )

    fig.update_layout(
        title=dict(text="Sequence similarity heatmap (bitscore)", x=0.5),
        height=700,
        autosize=True,
        legend=dict(title=dict(text="Split"), x=1.05, y=0.5),
        margin=dict(l=10, r=120, t=50, b=10),
    )

    # Hide axes ticks/labels on annotation strips and main heatmap
    for axis in ("xaxis", "yaxis", "xaxis2", "yaxis2", "xaxis3", "yaxis3"):
        fig.update_layout(**{axis: dict(showticklabels=False, ticks="")})

    div = fig.to_html(full_html=False, include_plotlyjs="cdn", config={"responsive": True})

    with open(output, "w") as fh:
        fh.write(
            "<!--\n"
            f"id: 'similarity_heatmap_{id_}'\n"
            f"section_name: 'Sequence Similarity Heatmap: {id_}'\n"
            "description: 'Clustered bitscore matrix for a sample of proteins, coloured by "
            "train / val / test split. Proteins from the same split should cluster together "
            "with high within-split similarity and reduced cross-split similarity.'\n"
            "-->\n"
        )
        fh.write(div)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train_fasta",    required=True)
    ap.add_argument("--val_fasta",      required=True)
    ap.add_argument("--test_fasta",     required=True)
    ap.add_argument("--blast",          required=True)
    ap.add_argument("--max_per_split",  type=int, default=200)
    ap.add_argument("--seed",           type=int, default=42)
    ap.add_argument("--output",         default="similarity_heatmap_mqc.html")
    ap.add_argument("--id", required=True, help="Dataset ID, for MultiQC tagging")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    split_order = [("train", args.train_fasta), ("val", args.val_fasta), ("test", args.test_fasta)]
    all_ids, split_labels = [], []
    for split, path in split_order:
        ids = read_fasta_ids(path)
        if len(ids) > args.max_per_split:
            ids = list(rng.choice(ids, size=args.max_per_split, replace=False))
        all_ids.extend(ids)
        split_labels.extend([split] * len(ids))
        print(f"  {split}: {len(ids)} proteins sampled", file=sys.stderr)

    print(f"Building {len(all_ids)}×{len(all_ids)} bitscore matrix …", file=sys.stderr)
    mat = build_bitscore_matrix(args.blast, all_ids)

    print("Clustering …", file=sys.stderr)
    order = cluster_order(mat)
    mat_clustered = mat[np.ix_(order, order)]
    labels_clustered = [split_labels[i] for i in order]

    print(f"Writing {args.output} …", file=sys.stderr)
    make_html(mat_clustered, labels_clustered, args.output, args.id)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()