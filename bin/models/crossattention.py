"""Cross-attention model (original architecture).

Taken from https://github.com/daisybio/PPI_prediction_study (models/attention.py,
CrossAttInteraction). It needs PER-RESIDUE embeddings: every protein is an
(L, embed_dim) array with its own length L. Each residue embedding is compressed
by 3 FC layers to `h3` dims, passed through a cross-attention encoder layer
attending to the other protein, combined in a pairwise (L1 x L2) interaction
map, convolved, pooled, and reduced by a max to one probability.

As in the original, every pair is processed on its own (lengths differ), and a
mini-batch is the stack of its per-pair predictions. Hyperparameters follow the
repo's best "crossattention" config (lr 1e-4, batch size 16, dropout 0.1,
8 heads, h3 64, ff_dim 256, avg pooling, kernel 3, max_seq_len 1000, subset_size 0.5,
early stopping patience 8, up to 25 epochs). As in the original, pairs where either
protein has more than max_len (1000) residues are dropped (not truncated).
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.utils.parametrizations import spectral_norm

from .torch_base import TorchClassifier


class Attention(nn.Module):
    def __init__(self, hid_dim, n_heads, dropout):
        super().__init__()
        assert hid_dim % n_heads == 0, "hid_dim must be divisible by n_heads"
        self.hid_dim, self.n_heads = hid_dim, n_heads
        self.w_q = spectral_norm(nn.Linear(hid_dim, hid_dim))
        self.w_k = spectral_norm(nn.Linear(hid_dim, hid_dim))
        self.w_v = spectral_norm(nn.Linear(hid_dim, hid_dim))
        self.fc = spectral_norm(nn.Linear(hid_dim, hid_dim))
        self.do = nn.Dropout(dropout)
        self.register_buffer("scale", torch.sqrt(torch.tensor(float(hid_dim // n_heads))))

    def forward(self, query, key, value, mask=None):
        bsz = query.shape[0]
        d_head = self.hid_dim // self.n_heads
        Q = self.w_q(query).view(bsz, -1, self.n_heads, d_head).permute(0, 2, 1, 3)
        K = self.w_k(key).view(bsz, -1, self.n_heads, d_head).permute(0, 2, 1, 3)
        V = self.w_v(value).view(bsz, -1, self.n_heads, d_head).permute(0, 2, 1, 3)

        energy = torch.matmul(Q, K.permute(0, 1, 3, 2)) / self.scale
        if mask is not None:
            energy = energy.masked_fill(mask == 0, -1e10)
        attention = self.do(F.softmax(energy, dim=-1))

        x = torch.matmul(attention, V).permute(0, 2, 1, 3).contiguous()
        return self.fc(x.view(bsz, -1, self.hid_dim))


class Feedforward(nn.Module):
    def __init__(self, hid_dim, ff_dim, dropout):
        super().__init__()
        self.fc_1 = spectral_norm(nn.Linear(hid_dim, ff_dim))
        self.fc_2 = spectral_norm(nn.Linear(ff_dim, hid_dim))
        self.do = nn.Dropout(dropout)
        self.activation = nn.SiLU()  # "swish"

    def forward(self, x):
        return self.fc_2(self.do(self.activation(self.fc_1(x))))


class CrossEncoderLayer(nn.Module):
    def __init__(self, hid_dim, n_heads, ff_dim, dropout):
        super().__init__()
        self.ln1, self.ln2 = nn.LayerNorm(hid_dim), nn.LayerNorm(hid_dim)
        self.do1, self.do2 = nn.Dropout(dropout), nn.Dropout(dropout)
        self.sa = Attention(hid_dim, n_heads, dropout)
        self.ff = Feedforward(hid_dim, ff_dim, dropout)

    def forward(self, trg, cross, mask=None):
        trg = self.ln1(trg + self.do1(self.sa(trg, cross, cross, mask)))
        return self.ln2(trg + self.do2(self.ff(trg)))


class CrossAttInteraction(nn.Module):
    def __init__(self, embed_dim, num_heads=8, h3=64, dropout=0.1, ff_dim=256, pooling="avg", kernel_size=3):
        super().__init__()
        h = int(embed_dim // 4)
        h2 = int(h // 4)

        self.cross_encoder = CrossEncoderLayer(h3, num_heads, ff_dim, dropout)
        self.conv = nn.Conv2d(h3, 1, kernel_size=kernel_size, padding="same")
        pool = nn.MaxPool2d if pooling == "max" else nn.AvgPool2d
        self.pool = pool(kernel_size=kernel_size)

        self.fc1 = nn.Linear(embed_dim, h)
        self.fc2 = nn.Linear(h, h2)
        self.fc3 = nn.Linear(h2, h3)

    def _compress(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        return F.relu(self.fc3(x))

    def forward(self, protein1, protein2, mask1=None, mask2=None):
        """protein1: (L1, embed_dim), protein2: (L2, embed_dim) -> (prob, interaction map)"""
        x1 = self._compress(protein1.to(torch.float32).unsqueeze(0))
        x2 = self._compress(protein2.to(torch.float32).unsqueeze(0))

        x1 = self.cross_encoder(x1, x2, mask1)
        x2 = self.cross_encoder(x2, x1, mask2)

        mat = torch.einsum("bik,bjk->bijk", x1, x2).permute(0, 3, 1, 2)
        mat = self.conv(mat)
        x = self.pool(mat)

        return torch.sigmoid(torch.max(x))[None], mat


class _PairNet(nn.Module):
    """Runs the model pair by pair on a list of (protein1, protein2) tensors -> (batch,) probabilities."""

    def __init__(self, net):
        super().__init__()
        self.net = net

    def forward(self, pairs):
        return torch.cat([self.net(p1, p2)[0] for p1, p2 in pairs])


class CrossAttentionClassifier(TorchClassifier):
    description = (
        "Cross-attention PPI classifier on per-residue embeddings (PPI_prediction_study). "
        "Adam lr 1e-4, early stopping on val AUROC; trained on train only."
    )
    lr = 1e-4
    batch_size = 16
    log_every = 30  # seconds between progress lines (an epoch takes very long)
    subset_size = 0.1  # as in the original: each epoch draws half of the training pairs (with replacement)
    eval_batch_size = 16
    per_residue = True
    min_len = 3  # the 3x3 average pooling needs an interaction map of at least 3x3
    max_len = 1000  # original max_seq_len: longer pairs are dropped

    def _embed_dim(self, Xp):
        return Xp[0][0].shape[-1]

    def _prepare(self, X):
        return X  # list of (array(L1, d), array(L2, d))

    def _batch(self, Xp, idx):
        def to_t(a):
            return torch.as_tensor(np.asarray(a), dtype=torch.float32, device=self.device)

        return [(to_t(Xp[i][0]), to_t(Xp[i][1])) for i in idx.tolist()]

    def build_net(self, embed_dim):
        return _PairNet(CrossAttInteraction(embed_dim))
