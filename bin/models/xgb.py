"""XGBoost baseline with limited hyperparameter tuning.

3 pre-defined configs (max_depth 3/6/10) are evaluated by AUROC on the val set; each
fit uses early stopping on val AUROC to pick the number of trees (up to 500). The best
config is retrained on train+val with that number of trees.
"""

import sys

import numpy as np
from sklearn.metrics import roc_auc_score
from xgboost import XGBClassifier

HP_CONFIGS = [
    {"max_depth": 3, "learning_rate": 0.1, "subsample": 0.8, "colsample_bytree": 0.5},
    {"max_depth": 6, "learning_rate": 0.1, "subsample": 0.8, "colsample_bytree": 0.5},
    {"max_depth": 10, "learning_rate": 0.1, "subsample": 0.8, "colsample_bytree": 0.5},
]
MAX_TREES = 500
EARLY_STOPPING_ROUNDS = 20


class XGBoostClassifier:
    description = (
        "XGBoost PPI classifier. Hyperparameters tuned on val AUROC (3 configs: max_depth 3/6/10, "
        "early stopping for the number of trees), then retrained on train+val."
    )

    per_residue = False  # works on pooled per-protein embeddings
    min_len = 1
    max_len = None

    def __init__(self, seed=42):
        self.seed = seed
        self.clf = None

    def _make(self, cfg, n_estimators, **kwargs):
        return XGBClassifier(
            **cfg, n_estimators=n_estimators, random_state=self.seed, n_jobs=-1, tree_method="hist", **kwargs
        )

    def fit(self, X_train, y_train, X_val, y_val):
        print("Tuning hyperparameters ...", file=sys.stderr)
        best_auroc, best_cfg, best_trees = -1.0, None, None
        for cfg in HP_CONFIGS:
            clf = self._make(cfg, MAX_TREES, eval_metric="auc", early_stopping_rounds=EARLY_STOPPING_ROUNDS)
            clf.fit(X_train, y_train, eval_set=[(X_val, y_val)], verbose=False)
            auroc = roc_auc_score(y_val, clf.predict_proba(X_val)[:, 1])
            print(f"  {cfg}  trees={clf.best_iteration + 1}  →  val AUROC {auroc:.4f}", file=sys.stderr)
            if auroc > best_auroc:
                best_auroc, best_cfg, best_trees = auroc, cfg, clf.best_iteration + 1

        print(f"Best: {best_cfg} with {best_trees} trees  (val AUROC {best_auroc:.4f})", file=sys.stderr)

        # Retrain on train + val combined
        self.clf = self._make(best_cfg, best_trees)
        self.clf.fit(np.concatenate([X_train, X_val]), np.concatenate([y_train, y_val]))
        return self

    def predict_proba(self, X):
        return self.clf.predict_proba(X)[:, 1]
