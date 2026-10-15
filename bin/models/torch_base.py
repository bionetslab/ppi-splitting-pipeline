"""Shared training loop for the PyTorch baselines (Adam + BCE, early stopping on val AUROC)."""

import copy
import sys
import time

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import roc_auc_score


def _fmt(seconds):
    seconds = int(seconds)
    return f"{seconds // 3600}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}"


def get_device():
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():  # Apple silicon GPU
        return torch.device("mps")
    return torch.device("cpu")


class TorchClassifier:
    """Subclasses set the hyperparameters below and implement build_net(embed_dim).

    By default the network receives a tensor of shape (batch, 2, embed_dim)
    (protein 1 and protein 2) and returns probabilities of shape (batch,).
    Subclasses with other input formats override _prepare / _batch.
    """

    per_residue = False  # True: model needs per-residue embeddings (see train_classifier.build_pairs)
    min_len = 1          # shortest per-residue sequence the model can handle
    max_len = None       # longest per-residue sequence; pairs with longer proteins are dropped
    eval_batch_size = 4096

    description = ""
    lr = 1e-3
    batch_size = 256
    num_epochs = 25
    patience = 8  # epochs without val-AUROC improvement before stopping
    # fraction of the training pairs drawn per epoch. 1.0 = every pair once per epoch (random order);
    # < 1 = random draw with replacement of int(n * subset_size) pairs, as in PPI_prediction_study
    subset_size = 1.0
    log_every = None  # if set: print progress (batches, loss, speed, ETA) every this many seconds within an epoch

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
        print(f"Using device: {self.device}", file=sys.stderr)
        self.net = self.build_net(self._embed_dim(Xt)).to(self.device)
        opt = torch.optim.Adam(self.net.parameters(), lr=self.lr)
        loss_fn = nn.BCELoss()

        best_auroc, best_state, stale = -1.0, None, 0
        for epoch in range(self.num_epochs):
            self.net.train()
            n_epoch = int(len(Xt) * self.subset_size)
            if self.subset_size < 1.0:
                perm = torch.randint(len(Xt), (n_epoch,))
            else:
                perm = torch.randperm(len(Xt))
            total, seen = 0.0, 0
            t_start = t_last = time.time()
            for i in range(0, n_epoch, self.batch_size):
                idx = perm[i : i + self.batch_size]
                if len(idx) < 2:  # BatchNorm cannot handle a batch of one
                    continue
                opt.zero_grad()
                out = self.net(self._batch(Xt, idx))
                loss = loss_fn(out, yt[idx].to(self.device))
                loss.backward()
                opt.step()
                total += loss.item() * len(idx)
                seen += len(idx)
                if self.log_every and time.time() - t_last >= self.log_every:
                    t_last = time.time()
                    rate = seen / (t_last - t_start)
                    print(
                        f"    epoch {epoch + 1}: {seen}/{n_epoch} pairs ({seen / n_epoch:.0%})  "
                        f"loss {total / seen:.4f}  {rate:.1f} pairs/s  "
                        f"elapsed {_fmt(t_last - t_start)}  ETA {_fmt((n_epoch - seen) / rate)}",
                        file=sys.stderr,
                    )
            if self.log_every:
                print(f"    epoch {epoch + 1} training took {_fmt(time.time() - t_start)}; validating ...", file=sys.stderr)
            val_auroc = roc_auc_score(y_val, self._predict(Xv))
            print(
                f"  epoch {epoch + 1}/{self.num_epochs}  loss {total / n_epoch:.4f}  val AUROC {val_auroc:.4f}",
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
        out, t_last = [], time.time()
        for i in range(0, len(X), bs):
            out.append(self.net(self._batch(X, torch.arange(i, min(i + bs, len(X))))).cpu())
            if self.log_every and time.time() - t_last >= self.log_every:
                t_last = time.time()
                print(f"    predicted {min(i + bs, len(X))}/{len(X)} pairs", file=sys.stderr)
        return torch.cat(out).numpy()

    def predict_proba(self, X):
        return self._predict(self._prepare(X))
