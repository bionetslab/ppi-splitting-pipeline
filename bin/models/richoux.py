"""Richoux-like model: two parallel 2-layer FC towers (one per protein) -> concat -> FC -> 1.

Re-implementation of the fully connected model of Richoux et al. as used in
https://github.com/daisybio/PPI_prediction_study (models/fc2_20_2_dense.py).
Hyperparameters follow that repo's best "richoux" config (lr 1e-3, 20-unit
layers, spectral norm, early stopping patience 8, up to 25 epochs).
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.utils.parametrizations import spectral_norm

from .torch_base import TorchClassifier


class FC2_20_2Dense(nn.Module):
    def __init__(self, embed_dim, ff_dim1=20, ff_dim2=20, ff_dim3=20, spec_norm=True):
        super().__init__()
        sn = spectral_norm if spec_norm else (lambda layer: layer)

        # tower for protein 1
        self.fc1 = sn(nn.Linear(embed_dim, ff_dim1))
        self.fc2 = sn(nn.Linear(ff_dim1, ff_dim2))
        self.bn1 = nn.BatchNorm1d(ff_dim1)
        self.bn2 = nn.BatchNorm1d(ff_dim2)

        # tower for protein 2
        self.fc3 = sn(nn.Linear(embed_dim, ff_dim1))
        self.fc4 = sn(nn.Linear(ff_dim1, ff_dim2))
        self.bn3 = nn.BatchNorm1d(ff_dim1)
        self.bn4 = nn.BatchNorm1d(ff_dim2)

        # joint layers
        self.fc5 = sn(nn.Linear(2 * ff_dim2, ff_dim3))
        self.bn5 = nn.BatchNorm1d(ff_dim3)
        self.fc6 = sn(nn.Linear(ff_dim3, 1))

    def forward(self, x):
        """x: (batch, 2, embed_dim) -> probabilities (batch,)"""
        x1, x2 = x[:, 0], x[:, 1]

        x1 = self.bn1(F.relu(self.fc1(x1)))
        x1 = self.bn2(F.relu(self.fc2(x1)))

        x2 = self.bn3(F.relu(self.fc3(x2)))
        x2 = self.bn4(F.relu(self.fc4(x2)))

        x = torch.cat((x1, x2), dim=1)
        x = self.bn5(F.relu(self.fc5(x)))
        return torch.sigmoid(self.fc6(x)).squeeze(1)


class RichouxClassifier(TorchClassifier):
    description = (
        "Richoux-like PPI classifier (two FC towers, 20 units, spectral norm). "
        "Adam lr 1e-3, early stopping on val AUROC; trained on train only."
    )
    lr = 1e-3
    batch_size = 1024

    def build_net(self, embed_dim):
        return FC2_20_2Dense(embed_dim)
