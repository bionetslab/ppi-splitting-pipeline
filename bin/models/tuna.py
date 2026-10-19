"""TUnA model (batched, with padding masks).

Re-implementation of TUnA (Ko et al.) as in https://github.com/daisybio/PPI_prediction_study
(models/attention.py, class TUnA, non-cross variant, which is the repo's best "TUnA" config).
It needs PER-RESIDUE embeddings: every protein is an (L, embed_dim) array. Proteins are
truncated at max_len (1000) residues and zero-padded to the longest protein of the batch;
explicit validity masks keep the padding out of attention and pooling.

Architecture: spectral-norm linear projection to hid_dim (64) -> intra-protein encoder layer
(self-attention, per protein) -> both proteins concatenated in both orders (x12, x21) ->
inter-protein encoder layer -> masked mean pooling -> element-wise max of both orders ->
random-Fourier-feature head -> sigmoid.
As in the original, the attention mask of the inter encoder is block-diagonal, so the two
proteins only meet in the pooling step and the output head. Only the RFF head's logits are used
(no precision/covariance update, no variance), as in the original forward pass.

Hyperparameters follow the repo's best "TUnA" config (lr 5e-4, batch size 32 -> 16 here to fit
the 2L x 2L attention matrices in GPU memory, dropout 0.15, 8 heads, hid_dim 64, ff_dim 256,
1028 RFFs, subset_size 0.5, early stopping patience 8, up to 25 epochs).
"""

import math

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
        """mask: bool, broadcastable to (batch, heads, len_q, len_k); False = not attended."""
        bsz = query.shape[0]
        d_head = self.hid_dim // self.n_heads
        Q = self.w_q(query).view(bsz, -1, self.n_heads, d_head).permute(0, 2, 1, 3)
        K = self.w_k(key).view(bsz, -1, self.n_heads, d_head).permute(0, 2, 1, 3)
        V = self.w_v(value).view(bsz, -1, self.n_heads, d_head).permute(0, 2, 1, 3)

        energy = torch.matmul(Q, K.permute(0, 1, 3, 2)) / self.scale
        if mask is not None:
            energy = energy.masked_fill(~mask, -1e10)
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


class EncoderLayer(nn.Module):
    def __init__(self, hid_dim, n_heads, ff_dim, dropout):
        super().__init__()
        self.ln1, self.ln2 = nn.LayerNorm(hid_dim), nn.LayerNorm(hid_dim)
        self.do1, self.do2 = nn.Dropout(dropout), nn.Dropout(dropout)
        self.sa = Attention(hid_dim, n_heads, dropout)
        self.ff = Feedforward(hid_dim, ff_dim, dropout)

    def forward(self, x, mask=None):
        x = self.ln1(x + self.do1(self.sa(x, x, x, mask)))
        return self.ln2(x + self.do2(self.ff(x)))


class RFFHead(nn.Module):
    """Random-Fourier-feature output layer (logits of the original VanillaRFFLayer): fixed random
    projection (seeded, not trained) -> cos/sin features -> trained linear layer -> 1 logit."""

    def __init__(self, in_features, rffs=1028, amplitude=1.0, random_seed=123):
        super().__init__()
        assert rffs % 2 == 0 and rffs > 1, "rffs must be an even number greater than 1"
        num_freqs = rffs // 2
        gen = torch.Generator().manual_seed(random_seed)
        self.register_buffer("weight_mat", torch.randn(in_features, num_freqs, generator=gen))
        self.output_weights = nn.Parameter(torch.randn(rffs, 1, generator=gen))
        self.amplitude = amplitude
        self.feature_scale = math.sqrt(2.0 / num_freqs)

    def forward(self, x):
        proj = self.amplitude * x @ self.weight_mat
        feats = self.feature_scale * torch.cat([torch.cos(proj), torch.sin(proj)], dim=1)
        return feats @ self.output_weights


class TUnA(nn.Module):
    def __init__(self, embed_dim, num_heads=8, hid_dim=64, dropout=0.15, ff_dim=256, rffs=1028):
        super().__init__()
        self.lin1 = spectral_norm(nn.Linear(embed_dim, hid_dim))
        self.intra = EncoderLayer(hid_dim, num_heads, ff_dim, dropout)
        self.inter = EncoderLayer(hid_dim, num_heads, ff_dim, dropout)
        self.head = RFFHead(hid_dim, rffs)

    @staticmethod
    def _block_diag_mask(v_a, v_b):
        """(B, 1, 2L, 2L) bool: valid-residue square mask of protein a (top-left) and b (bottom-right)."""
        B, L = v_a.shape
        mask = torch.zeros(B, 1, 2 * L, 2 * L, dtype=torch.bool, device=v_a.device)
        mask[:, 0, :L, :L] = v_a[:, :, None] & v_a[:, None, :]
        mask[:, 0, L:, L:] = v_b[:, :, None] & v_b[:, None, :]
        return mask

    def forward(self, batch):
        """batch = (x1 (B, L, d), x2 (B, L, d), valid1 (B, L) bool, valid2 (B, L) bool) -> (B,) probabilities"""
        x1, x2, v1, v2 = batch
        # intra-protein encoder with square masks over the real residues of each protein
        x1 = self.intra(self.lin1(x1), (v1[:, :, None] & v1[:, None, :])[:, None])
        x2 = self.intra(self.lin1(x2), (v2[:, :, None] & v2[:, None, :])[:, None])

        # both orders, and the inter-protein encoder
        x12, x21 = torch.cat((x1, x2), dim=1), torch.cat((x2, x1), dim=1)
        v12, v21 = torch.cat((v1, v2), dim=1), torch.cat((v2, v1), dim=1)
        x12 = self.inter(x12, self._block_diag_mask(v1, v2))
        x21 = self.inter(x21, self._block_diag_mask(v2, v1))

        # mean over the real residues, then max over the two orders
        pooled = []
        for x, v in ((x12, v12), (x21, v21)):
            w = v.unsqueeze(-1).to(x.dtype)
            pooled.append((x * w).sum(dim=1) / w.sum(dim=1))
        feat = torch.maximum(pooled[0], pooled[1])
        return torch.sigmoid(self.head(feat)).squeeze(-1)


class TUnAClassifier(TorchClassifier):
    description = (
        "TUnA (spectral-norm transformer encoders + RFF head) on per-residue embeddings "
        "(PPI_prediction_study). Adam lr 5e-4, proteins truncated/padded to 1000 residues, "
        "early stopping on val AUROC; trained on train only."
    )
    lr = 5e-4
    batch_size = 16
    log_every = 30  # seconds between progress lines
    subset_size = 0.5  # as in the original: each epoch draws half of the training pairs (with replacement)
    eval_batch_size = 32
    per_residue = True
    max_len = None  # nothing is dropped: proteins are truncated in _batch instead
    trunc_len = 1000  # original max_seq_len

    def _embed_dim(self, Xp):
        return Xp[0][0].shape[-1]

    def _prepare(self, X):
        return X  # list of (array(L1, d), array(L2, d))

    def _batch(self, Xp, idx):
        pairs = [Xp[i] for i in idx.tolist()]
        out = []
        for k in (0, 1):
            seqs = [np.asarray(p[k])[: self.trunc_len] for p in pairs]
            L = max(len(s) for s in seqs)
            x = np.zeros((len(seqs), L, seqs[0].shape[-1]), dtype=np.float32)
            valid = np.zeros((len(seqs), L), dtype=bool)
            for j, s in enumerate(seqs):
                x[j, : len(s)] = s
                valid[j, : len(s)] = True
            out.append((x, valid))
        # both proteins padded to the same length (needed to concatenate them in the inter encoder)
        L = max(out[0][0].shape[1], out[1][0].shape[1])
        tensors = []
        for x, valid in out:
            pad = L - x.shape[1]
            if pad:
                x = np.pad(x, ((0, 0), (0, pad), (0, 0)))
                valid = np.pad(valid, ((0, 0), (0, pad)))
            tensors.append((torch.from_numpy(x).to(self.device), torch.from_numpy(valid).to(self.device)))
        (x1, v1), (x2, v2) = tensors
        return x1, x2, v1, v2

    def build_net(self, embed_dim):
        return TUnA(embed_dim)
