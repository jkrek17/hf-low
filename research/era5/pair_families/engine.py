"""Pair engine for RA-12: candidate (parent P, daughter D) pairs among HF tracks of one basin, optionally with whole tracks
shifted in time (permutation null). Definitions are those in PREREGISTRATION.md."""
import numpy as np
from common import hav, bearing

MAXD, WIN, RMIN, RMAX = 1500.0, 48.0, 120.0, 240.0


def prep(G, F):
    """G: tracks (one basin) DataFrame; F: dict track->fixes. Returns dict of arrays."""
    fx = []
    for tr in G.track.values:
        f = F[tr]
        fx.append((f.t.values, f.lat.values, np.degrees(np.unwrap(np.radians(f.lon.values)))))
    return dict(t0=G.t0.values, t1=G.t1.values, tpk=G.tpeak.values, lat0=G.lat0.values, lon0=G.lon0.values, fx=fx,
                gen=G.gen_obs.fillna(False).astype(bool).values, lag0=G.lag0.values, n=len(G))


def pairs(P, shift=None):
    """All candidate pairs. Returns arrays i (parent), j (daughter), dist, rel (bearing of D from P minus P heading, 0-360), dt."""
    n = P["n"]
    s = np.zeros(n) if shift is None else shift
    t0, t1, tpk = P["t0"] + s, P["t1"] + s, P["tpk"] + s
    M = (t0[:, None] < t0[None, :]) & (t1[:, None] >= t0[None, :]) & (np.abs(t0[None, :] - tpk[:, None]) <= WIN)
    np.fill_diagonal(M, False)
    ii, jj = np.nonzero(M)
    out = []
    for i, j in zip(ii, jj):
        t, la, lo = P["fx"][i]
        tq = t0[j] - s[i]                     # D's first-fix time in P's own (unshifted) clock
        if tq < t[0] or tq > t[-1]:
            continue
        pla, plo = np.interp(tq, t, la), np.interp(tq, t, lo)
        k = min(max(np.searchsorted(t, tq), 1), len(t) - 1) if len(t) > 1 else 0
        if len(t) < 2:
            continue
        hd = bearing(la[k - 1], lo[k - 1] % 360, la[k], lo[k] % 360)
        d = hav(pla, plo % 360, P["lat0"][j], P["lon0"][j])
        if d <= MAXD:
            b = bearing(pla, plo % 360, P["lat0"][j], P["lon0"][j])
            out.append((i, j, d, (b - hd) % 360, t0[j] - tpk[i]))
    if not out:
        return np.zeros((0, 5))
    return np.array(out)


def strata(G):
    ts = (np.datetime64("2000-01-01") + (G.tpeak.values * 3600).astype("timedelta64[s]")).astype("datetime64[M]")
    return np.array([f"{G.season.values[k]}-{ts[k]}" for k in range(len(G))])


def perm_shift(P, st, rng):
    s = np.zeros(P["n"])
    for u in np.unique(st):
        idx = np.nonzero(st == u)[0]
        new = rng.permutation(P["tpk"][idx])
        s[idx] = new - P["tpk"][idx]
    return s


def summarise(pr, P, rear=True):
    """Counts: all candidates; daughters with observed genesis; rear-sector (RMIN..RMAX); each over all/obs."""
    if len(pr) == 0:
        return dict(all=0, obs=0, rear=0, rear_obs=0)
    j = pr[:, 1].astype(int)
    obs = P["gen"][j] & (P["lag0"][j] <= 12)
    r = (pr[:, 3] >= RMIN) & (pr[:, 3] <= RMAX)
    return dict(all=len(pr), obs=int(obs.sum()), rear=int(r.sum()), rear_obs=int((r & obs).sum()))
