#!/usr/bin/env python3
"""
Train a PPI classifier on protein embeddings and evaluate on both test sets.

Feature construction: concatenation of the two sorted protein embeddings.
Sorting by protein ID ensures the feature vector is the same regardless of
the order in which the pair appears in the CSV.

Models (--model, definitions in bin/models/):
  RFC (default)   RandomForest. 3 pre-defined configs (max_depth 5/10/30,
                  max_samples 0.2) evaluated by AUROC on val set; the best
                  config is retrained on train+val.
  Crossattention  cross-attention network (PPI_prediction_study); needs
                  per-residue embeddings, i.e. an (L, d) array per protein in the NPZ
  Richoux-like    fully connected two-tower network (Richoux et al.)
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from models import MODELS, get_model
from utils import load_embeddings, load_embeddings_lazy, read_labelled_csv
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)


def build_X(pairs, labels, embeddings):
    """Return (X, y), silently dropping pairs where either protein is missing."""
    rows, y = [], []
    skipped = 0
    for (p1, p2), label in zip(pairs, labels):
        a, b = (p1, p2) if p1 <= p2 else (p2, p1)
        if a in embeddings and b in embeddings:
            rows.append(np.concatenate([embeddings[a], embeddings[b]]))
            y.append(label)
        else:
            skipped += 1
    if skipped:
        print(f"  skipped {skipped} pairs with missing embeddings", file=sys.stderr)
    return np.array(rows, dtype=np.float32), np.array(y)


def build_pairs(pairs, labels, embeddings, min_len=1, max_len=None):
    """Like build_X, but keep the two per-residue embeddings (L, d) of every pair
    separate (they have different lengths). Pairs with a missing protein, or a protein
    shorter than min_len / longer than max_len residues, are dropped."""
    rows, y = [], []
    skipped = 0
    for (p1, p2), label in zip(pairs, labels):
        a, b = (p1, p2) if p1 <= p2 else (p2, p1)
        if a in embeddings and b in embeddings and all(
            len(embeddings[p]) >= min_len and (max_len is None or len(embeddings[p]) <= max_len) for p in (a, b)
        ):
            rows.append((embeddings[a], embeddings[b]))
            y.append(label)
        else:
            skipped += 1
    if skipped:
        print(f"  skipped {skipped} pairs with missing embeddings or protein length outside [{min_len}, {max_len}]", file=sys.stderr)
    return rows, np.array(y)


def compute_metrics(y_true, y_prob):
    y_pred = (y_prob >= 0.5).astype(int)
    return {
        "auroc": roc_auc_score(y_true, y_prob),
        "auprc": average_precision_score(y_true, y_prob),
        "f1": f1_score(y_true, y_pred, zero_division=0),
        "mcc": matthews_corrcoef(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall": recall_score(y_true, y_pred, zero_division=0),
        "accuracy": accuracy_score(y_true, y_pred),
    }


def write_mqc(results, id_, model_name, description):
    """One table per test split (test_balanced / test_realistic), each
    merging across datasets -- keeps either table from getting cluttered
    with the other's rows, matching the split evaluation of the two test
    sets."""
    cols = ["auroc", "auprc", "f1", "mcc", "precision", "recall", "accuracy"]
    for name, metrics in results:
        table_id = f"classifier_metrics_{name}"
        with open(f"{table_id}_mqc.tsv", "w") as fh:
            fh.write(
                f"# id: '{table_id}'\n"
                f"# section_name: 'Classifier Performance ({name})'\n"
                f"# description: '{description}'\n"
                "# plot_type: 'table'\n"
                "# pconfig:\n"
                f"#     id: '{table_id}_table'\n"
                f"#     title: '{model_name} Classifier - {name} Test Performance'\n"
                "# headers:\n"
                "#     ID:         {description: 'Dataset ID'}\n"
                "#     AUROC:      {format: '{:.4f}'}\n"
                "#     AUPRC:      {format: '{:.4f}'}\n"
                "#     F1:         {format: '{:.4f}'}\n"
                "#     MCC:        {format: '{:.4f}'}\n"
                "#     Precision:  {format: '{:.4f}'}\n"
                "#     Recall:     {format: '{:.4f}'}\n"
                "#     Accuracy:   {format: '{:.4f}'}\n"
                "Sample\tID\tAUROC\tAUPRC\tF1\tMCC\tPrecision\tRecall\tAccuracy\n"
            )
            row = "\t".join(f"{metrics[c]:.4f}" for c in cols)
            fh.write(f"{id_}\t{id_}\t{row}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True)
    ap.add_argument("--val", required=True)
    ap.add_argument("--test_balanced", required=True)
    ap.add_argument("--test_realistic", required=True)
    ap.add_argument("--embeddings", required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--id", required=True, help="Dataset ID, for MultiQC tagging")
    ap.add_argument("--model", default="RFC", choices=MODELS, help="Classifier to train (default: RFC)")
    args = ap.parse_args()

    model = get_model(args.model, seed=args.seed)
    print("Loading embeddings ...", file=sys.stderr)
    # per-residue embeddings are read from disk on demand (they do not fit in memory)
    embeddings = load_embeddings_lazy(args.embeddings) if model.per_residue else load_embeddings(args.embeddings)
    print(f"  {len(embeddings)} proteins", file=sys.stderr)

    if model.per_residue:
        bad = [k for k, v in embeddings.items() if v.ndim != 2]
        if bad:
            sys.exit(f"{args.model} needs per-residue embeddings (L, d) per protein; got 1-D vectors, e.g. {bad[0]}")

        def build(pairs, labels):
            return build_pairs(pairs, labels, embeddings, model.min_len, model.max_len)
    else:

        def build(pairs, labels):
            return build_X(pairs, labels, embeddings)

    train_pairs, y_train = read_labelled_csv(args.train)
    val_pairs, y_val = read_labelled_csv(args.val)
    X_train, y_train = build(train_pairs, y_train)
    X_val, y_val = build(val_pairs, y_val)

    print(f"Training {args.model} ...", file=sys.stderr)
    model.fit(X_train, y_train, X_val, y_val)

    # Evaluate on both test sets
    results = []
    for name, path in [("test_balanced", args.test_balanced), ("test_realistic", args.test_realistic)]:
        pairs, y_test_raw = read_labelled_csv(path)
        X_test, y_test = build(pairs, y_test_raw)
        y_prob = model.predict_proba(X_test)
        metrics = compute_metrics(y_test, y_prob)
        results.append((name, metrics))
        print(
            f"{name}: " + "  ".join(f"{k}={v:.4f}" for k, v in metrics.items()),
            file=sys.stderr,
        )

    #write_mqc(results, args.id, args.model, model.description)
    print(results)


if __name__ == "__main__":
    main()
