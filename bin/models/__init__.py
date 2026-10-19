"""Classifier definitions for train_classifier.py.

Every model exposes the same interface:

    model.fit(X_train, y_train, X_val, y_val)
    model.predict_proba(X)  ->  1-D array with P(interaction)
    model.description       ->  one-line description for the MultiQC table

X is the feature matrix built by train_classifier.build_X: the concatenation
of the two (ID-sorted) protein embeddings, i.e. shape (n_pairs, 2 * embed_dim).

Models are imported lazily so that RFC runs do not need to import torch.
"""

MODELS = ("RFC", "XGBoost", "TUnA", "Richoux-like")


def get_model(name, seed=42):
    if name == "RFC":
        from .rfc import RFCClassifier as cls
    elif name == "XGBoost":
        from .xgb import XGBoostClassifier as cls
    elif name == "TUnA":
        from .tuna import TUnAClassifier as cls
    elif name == "Richoux-like":
        from .richoux import RichouxClassifier as cls
    else:
        raise ValueError(f"Unknown model {name!r}; choose from {list(MODELS)}")
    return cls(seed=seed)
