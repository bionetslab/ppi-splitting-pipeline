"""Shared training loop for the PyTorch baselines (Adam + BCE, early stopping on val AUROC)."""

import copy
import sys

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import roc_auc_score


def get_device():
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


class TorchClassifier:
    """Subclasses set the hyperparameters below and implement build_net(embed_dim).

    By default the network receives a tensor of shape (batch, 2, embed_dim)
    (protein 1 and protein 2) and returns probabilities of shape (batch,).
    Subclasses with other input formats override _prepare / _batch.
    """

    per_residue = False  # True: model needs per-residue embeddings (see train_classifier.build_pairs)
    min_len = 1          # shortest per-residue sequence the model can handle
    eval_batch_size = 4096

    description = ""
    lr = 1e-3
    batch_size = 256
    num_epochs = 25
    patience = 8  # epochs without val-AUROC improvement before stopping

    def __init__(self, seed=42):
        self.seed = seed
        self.device = get_device()
        self.net = None

    def build_net(self, embed_dim):
        raise NotImplementedError

    def _embed_dim(self, Xp):
        return Xp.shape[-1]

    def _prepare(self, X):
        X = torch.as_tensor(np.asarray(X), dtype=torch.float32)
        return X.view(X.shape[0], 2, -1)

    def _batch(self, Xp, idx):
        return Xp[idx].to(self.device)

    def fit(self, X_train, y_train, X_val, y_val):
        torch.manual_seed(self.seed)
        np.random.seed(self.seed)
        Xt, Xv = self._prepare(X_train), self._prepare(X_val)
        yt = torch.as_tensor(np.asarray(y_train), dtype=torch.float32)
        self.net = self.build_net(self._embed_dim(Xt)).to(self.device)
        opt = torch.optim.Adam(self.net.parameters(), lr=self.lr)
        loss_fn = nn.BCELoss()

        best_auroc, best_state, stale = -1.0, None, 0
        for epoch in range(self.num_epochs):
            self.net.train()
            perm = torch.randperm(len(Xt))
            total = 0.0
            for i in range(0, len(Xt), self.batch_size):
                idx = perm[i : i + self.batch_size]
                if len(idx) < 2:  # BatchNorm cannot handle a batch of one
                    continue
                opt.zero_grad()
                out = self.net(self._batch(Xt, idx))
                loss = loss_fn(out, yt[idx].to(self.device))
                loss.backward()
                opt.step()
                total += loss.item() * len(idx)
            val_auroc = roc_auc_score(y_val, self._predict(Xv))
            print(
                f"  epoch {epoch + 1}/{self.num_epochs}  loss {total / len(Xt):.4f}  val AUROC {val_auroc:.4f}",
                file=sys.stderr,
            )
            if val_auroc > best_auroc:
                best_auroc, stale = val_auroc, 0
                best_state = copy.deepcopy(self.net.state_dict())
            else:
                stale += 1
                if stale >= self.patience:
                    print(f"  early stopping (best val AUROC {best_auroc:.4f})", file=sys.stderr)
                    break
        self.net.load_state_dict(best_state)
        return self

    @torch.no_grad()
    def _predict(self, X):
        self.net.eval()
        bs = self.eval_batch_size
        out = [
            self.net(self._batch(X, torch.arange(i, min(i + bs, len(X))))).cpu()
            for i in range(0, len(X), bs)
        ]
        return torch.cat(out).numpy()

    def predict_proba(self, X):
        return self._predict(self._prepare(X))
