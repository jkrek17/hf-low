"""Sector masks, model specifications and region-restricted EOFs for the upstream/local ablation.
Every definition here is the one in PREREGISTRATION.md; nothing is chosen from an outcome."""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hemispheric"))
import hemlib as H

LONC = ((5.625 * np.arange(64) + 180.0) % 360.0) - 180.0
LONC[LONC == -180.0] = 180.0                                  # centre at 180E belongs to (a, 180]


def cols(a, b):
    """Boolean over 64 longitudes: centre in the half-open arc (a, b], going east, wrapping the dateline."""
    return (LONC > a) | (LONC <= b) if a >= b else (LONC > a) & (LONC <= b)


def cells(a, b):
    m = np.zeros((12, 64), bool)
    m[:, cols(a, b)] = True
    return m


# sectors per target basin: (west edge, east edge) of the half-open arc
SECTORS = {
    "atl": {"UP": (140, -100), "LOC": (-100, 20), "DN": (20, 140)},
    "pac": {"UP": (0, 120), "LOC": (120, -120), "DN": (-120, 0)},
}
V1 = {"UP1": (140, -70), "LOC1": (-70, 20)}                      # Atlantic only: trough core upstream
ALLV = tuple(H.VARS)


def specs(basin):
    """List of dicts: name, basin, cells (12x64 bool), vars, base ('b1' or 'b2'), family."""
    S = SECTORS[basin]
    out = []

    def add(name, m, vars_=ALLV, base="b1", fam="primary"):
        out.append(dict(name=name, basin=basin, cells=m, vars=tuple(vars_), base=base, fam=fam))

    full = np.ones((12, 64), bool)
    sec = {k: cells(*v) for k, v in S.items()}
    add("FULL", full)
    for k in ("UP", "LOC", "DN"):
        add("only_" + k, sec[k])
    add("no_UP", sec["LOC"] | sec["DN"])
    add("no_LOC", sec["UP"] | sec["DN"])
    add("no_DN", sec["UP"] | sec["LOC"])
    if basin == "atl":
        add("only_UP1", cells(*V1["UP1"]), fam="V1")
        add("only_LOC1", cells(*V1["LOC1"]), fam="V1")
    for k in ("UP", "LOC"):
        add(f"only_{k}_z500", sec[k], ("z500",), fam="V2")
        add(f"only_{k}_u250", sec[k], ("u250",), fam="V2")
    for j in range(16):
        c = -180.0 + 22.5 * j
        a, b = c - 60, c + 60
        wrap = lambda x: ((x + 180) % 360) - 180 if ((x + 180) % 360) != 0 else 180.0
        add(f"win_{c:+06.1f}", cells(wrap(a), wrap(b)), fam="V3")
    for nm, m in (("FULL", full), ("only_UP", sec["UP"]), ("only_LOC", sec["LOC"])):
        add(nm + "_b2", m, base="b2", fam="imprint")
    return out


class REOFs:
    """H.EOFs restricted to a cell mask and a variable list (K leading EOFs per variable)."""

    def __init__(self, maps, ocean, region, vars_, k=H.K):
        self.vars = tuple(vars_)
        self.w = np.sqrt(np.cos(np.radians(H.LAT)))[:, None] * np.ones((1, 64))
        self.cell = {v: (region & ocean) if v == "sst" else region for v in self.vars}
        self.mean, self.vec, self.sd = {}, {}, {}
        for v in self.vars:
            X = self._flat(maps[v], v)
            self.mean[v] = X.mean(0)
            U, S, Vt = np.linalg.svd(X - self.mean[v], full_matrices=False)
            self.vec[v] = Vt[:k]
            self.sd[v] = ((X - self.mean[v]) @ Vt[:k].T).std(0, ddof=1)

    def _flat(self, M, v):
        return (M * self.w)[:, self.cell[v]]

    def scores(self, maps):
        return np.hstack([((self._flat(maps[v], v) - self.mean[v]) @ self.vec[v].T) / self.sd[v] for v in self.vars])
