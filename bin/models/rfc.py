"""RandomForest baseline with limited hyperparameter tuning.

3 pre-defined configs (max_depth 5/10/30, max_samples 0.2) are evaluated by
AUROC on the val set. The best config is retrained on train+val.
"""

import sys

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

HP_CONFIGS = [
    {"n_estimators": 100, "max_depth": 5, "max_samples": 0.2},
    {"n_estimators": 100, "max_depth": 10, "max_samples": 0.2},
    {"n_estimators": 100, "max_depth": 30, "max_samples": 0.2},
]


class RFCClassifier:
    description = (
        "RandomForest PPI classifier. Hyperparameters tuned on val AUROC "
        "(3 configs: max_depth 5/10/30, max_samples 0.2), then retrained on train+val."
    )

    per_residue = False  # works on pooled per-protein embeddings
    min_len = 1
    max_len = None

    def __init__(self, seed=42):
        self.seed = seed
        self.clf = None

    def fit(self, X_train, y_train, X_val, y_val):
        print("Tuning hyperparameters ...", file=sys.stderr)
        best_auroc, best_cfg = -1.0, None
        for cfg in HP_CONFIGS:
            clf = RandomForestClassifier(**cfg, random_state=self.seed, n_jobs=-1)
            clf.fit(X_train, y_train)
            auroc = roc_auc_score(y_val, clf.predict_proba(X_val)[:, 1])
            print(f"  {cfg}  →  val AUROC {auroc:.4f}", file=sys.stderr)
            if auroc > best_auroc:
                best_auroc, best_cfg = auroc, cfg

        print(f"Best: {best_cfg}  (val AUROC {best_auroc:.4f})", file=sys.stderr)

        # Retrain on train + val combined
        self.clf = RandomForestClassifier(**best_cfg, random_state=self.seed, n_jobs=-1)
        self.clf.fit(np.concatenate([X_train, X_val]), np.concatenate([y_train, y_val]))
        return self

    def predict_proba(self, X):
        return self.clf.predict_proba(X)[:, 1]
