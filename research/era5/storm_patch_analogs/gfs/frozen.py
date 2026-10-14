"""Frozen ERA5 pattern features: scaler, patch mean and PCA fitted once on ERA5 seasons 2014-15..2020-21 (see fit_frozen.py)."""
import os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
NP = 13
_a = (np.arange(NP) - NP // 2) * 154.3333
PR = np.hypot(*np.meshgrid(_a, _a, indexing="xy"))


class Frozen:
    def __init__(self, path):
        z = np.load(path)
        self.sd, self.mask, self.mu, self.comp, self.pm = z["sd"], z["mask"], z["mu"], z["comp"], z["pm"]

    def flat(self, X):
        Z = np.nan_to_num(X[:, :4].astype(np.float32) / self.sd[:, :4], nan=0.0)
        return Z[:, :, self.mask].reshape(len(X), -1)

    def pcs(self, X):
        return (self.flat(X) - self.mu - self.pm) @ self.comp.T
