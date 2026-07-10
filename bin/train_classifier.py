#!/usr/bin/env python3
"""
Train a RandomForest PPI classifier with limited hyperparameter tuning.

Feature construction: concatenation of the two sorted protein embeddings.
Sorting by protein ID ensures the feature vector is the same regardless of
the order in which the pair appears in the CSV.

Hyperparameter search: 3 pre-defined configs (max_depth 5/10/30, max_samples 0.2)
evaluated by AUROC on val set. The best config is retrained on train+val, then
evaluated on both test sets.
"""

import argparse
import os
import sys

import numpy as np
from sklearn.ensemble import RandomForestClassifier

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import load_embeddings, read_labelled_csv
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)

HP_CONFIGS = [
    {"n_estimators": 100, "max_depth": 5, "max_samples": 0.2},
    {"n_estimators": 100, "max_depth": 10, "max_samples": 0.2},
    {"n_estimators": 100, "max_depth": 30, "max_samples": 0.2},
]



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
    return np.array(rows), np.array(y)


def compute_metrics(y_true, y_prob):
    y_pred = (y_prob >= 0.5).astype(int)
    return {
        "auroc":     roc_auc_score(y_true, y_prob),
        "auprc":     average_precision_score(y_true, y_prob),
        "f1":        f1_score(y_true, y_pred, zero_division=0),
        "mcc":       matthews_corrcoef(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall":    recall_score(y_true, y_pred, zero_division=0),
        "accuracy":  accuracy_score(y_true, y_pred),
    }


def write_mqc(results, id_):
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
                "# description: 'RandomForest PPI classifier. Hyperparameters tuned on val AUROC (3 configs: max_depth 5/10/30, max_samples 0.2), then retrained on train+val.'\n"
                "# plot_type: 'table'\n"
                "# pconfig:\n"
                f"#     id: '{table_id}_table'\n"
                f"#     title: 'RF Classifier - {name} Test Performance'\n"
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
    ap.add_argument("--train",          required=True)
    ap.add_argument("--val",            required=True)
    ap.add_argument("--test_balanced",  required=True)
    ap.add_argument("--test_realistic", required=True)
    ap.add_argument("--embeddings",     required=True)
    ap.add_argument("--seed",           type=int, default=42)
    ap.add_argument("--id", required=True, help="Dataset ID, for MultiQC tagging")
    args = ap.parse_args()

    print("Loading embeddings ...", file=sys.stderr)
    embeddings = load_embeddings(args.embeddings)
    print(f"  {len(embeddings)} proteins", file=sys.stderr)

    train_pairs, y_train = read_labelled_csv(args.train)
    val_pairs,   y_val   = read_labelled_csv(args.val)

    X_train, y_train = build_X(train_pairs, y_train, embeddings)
    X_val,   y_val   = build_X(val_pairs,   y_val,   embeddings)

    # Hyperparameter search on val AUROC
    print("Tuning hyperparameters ...", file=sys.stderr)
    best_auroc, best_cfg = -1.0, None
    for cfg in HP_CONFIGS:
        clf = RandomForestClassifier(**cfg, random_state=args.seed, n_jobs=-1)
        clf.fit(X_train, y_train)
        auroc = roc_auc_score(y_val, clf.predict_proba(X_val)[:, 1])
        print(f"  {cfg}  →  val AUROC {auroc:.4f}", file=sys.stderr)
        if auroc > best_auroc:
            best_auroc, best_cfg = auroc, cfg

    print(f"Best: {best_cfg}  (val AUROC {best_auroc:.4f})", file=sys.stderr)

    # Retrain on train + val combined
    X_all     = np.concatenate([X_train, X_val])
    y_all     = np.concatenate([y_train, y_val])
    final_clf = RandomForestClassifier(**best_cfg, random_state=args.seed, n_jobs=-1)
    final_clf.fit(X_all, y_all)

    # Evaluate on both test sets
    results = []
    for name, path in [("test_balanced", args.test_balanced), ("test_realistic", args.test_realistic)]:
        pairs, y_test_raw = read_labelled_csv(path)
        X_test, y_test    = build_X(pairs, y_test_raw, embeddings)
        y_prob            = final_clf.predict_proba(X_test)[:, 1]
        metrics = compute_metrics(y_test, y_prob)
        results.append((name, metrics))
        print(
            f"{name}: " + "  ".join(f"{k}={v:.4f}" for k, v in metrics.items()),
            file=sys.stderr,
        )

    write_mqc(results, args.id)


if __name__ == "__main__":
    main()
