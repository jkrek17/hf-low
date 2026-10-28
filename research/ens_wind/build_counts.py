"""Turn pulled member maps + truth into compact count tables: for each (model, truth, basin, threshold, run date, lead) the number of
cases n[k] and events y[k] among cells where k members reached the threshold. Everything downstream (reliability, slope, Brier,
isotonic calibration, bootstrap) works from these counts. Truth 'ascat': 1 degree ocean cell, ASCAT max >= T (cells with >= 50 % of
their ocean sub-cells observed within +-3 h). Truth 'opc': neighbourhood (600 km ~ 324 nm around every 3rd cell), forecast = k members
whose neighbourhood max >= T, event = an OPC archive low fix of the category set within 600 km in valid time +-6 h."""
import os, sys, datetime as dt
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ens_fetch as F

OPC_DIR = os.path.join(HERE, "..", "..", "data", "hf_lows")
CATS = {48: {"S", "DS", "DHF", "HF"}, 64: {"HF"}}
STRIDE = 3
RAD = 600.0
RE = 6371.0


def opc_fixes():
    out = []
    for b, f in (("natl", "HF_Data_-_Atl.csv"), ("npac", "HF_Data_-_Pac.csv")):
        d = pd.read_csv(os.path.join(OPC_DIR, f), dtype={"date": str})
        d = d[d.date.str.len() == 10].copy()          # drops the malformed 11-character Pacific date
        d["t"] = pd.to_datetime(d.date, format="%Y%m%d%H", errors="coerce")
        d = d.dropna(subset=["t"]); d["basin"] = b
        d["lon"] = d.Longitude % 360
        out.append(d[["basin", "t", "Latitude", "lon", "Category"]].rename(columns={"Latitude": "lat"}))
    return pd.concat(out)


def gc(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * RE * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def cell_centres(b):
    ny = b["n"] - b["s"]; nx = b["e"] - b["w"]
    la = b["n"] - np.arange(ny) - 0.5; lo = (b["w"] + np.arange(nx) + 0.5) % 360
    return np.meshgrid(la, lo, indexing="ij")


_NB = {}


def neighbours(basin):
    """centres (every STRIDE-th cell) and, for each, flat indices of cells within RAD."""
    if basin in _NB: return _NB[basin]
    b = F.BOX[basin]; LA, LO = cell_centres(b); ny, nx = LA.shape
    cr = [(r, c) for r in range(0, ny, STRIDE) for c in range(0, nx, STRIDE)]
    idx = []
    for r, c in cr:
        dla = int(RAD / 111.0) + 2; dlo = int(RAD / (111.0 * max(np.cos(np.radians(LA[r, c])), 0.2))) + 2
        r0, r1 = max(0, r - dla), min(ny, r + dla + 1); c0, c1 = max(0, c - dlo), min(nx, c + dlo + 1)
        rr, cc = np.meshgrid(np.arange(r0, r1), np.arange(c0, c1), indexing="ij")
        ok = gc(LA[r, c], LO[r, c], LA[rr, cc], LO[rr, cc]) <= RAD
        idx.append((rr[ok] * nx + cc[ok]).astype(np.int64))
    _NB[basin] = (np.array([cr_[0] for cr_ in cr]), np.array([cr_[1] for cr_ in cr]), idx, LA, LO)
    return _NB[basin]


def counts_from(kvec, yvec, N):
    n = np.bincount(kvec, minlength=N + 1); y = np.bincount(kvec, weights=yvec.astype(float), minlength=N + 1)
    return n.astype(np.int64), y.astype(np.int64)


def build(model, seasons, leads, fixes):
    rows = []; keys = []
    for season in seasons:
        for ymd in F.run_dates(season):
            run = dt.datetime.strptime(ymd, "%Y%m%d")
            for lead in leads:
                shift = dt.timedelta(days=int(os.environ.get("PLACEBO_DAYS", "0")))      # placebo: truth taken from another date; 0 in all real runs
                p = os.path.join(F.WORK, model, "%s_%03d.npz" % (ymd, lead))
                if not os.path.exists(p): continue
                z = np.load(p); V = run + dt.timedelta(hours=lead)
                tz = os.path.join(F.WORK, "ascat", (V + shift).strftime("%Y%m%d") + ".npz")
                a = np.load(tz) if os.path.exists(tz) else None
                fx = fixes[(fixes.t >= V + shift - dt.timedelta(hours=6)) & (fixes.t <= V + shift + dt.timedelta(hours=6))]
                for basin in F.BOX:
                    m = z[basin]; N = m.shape[0]
                    # truth 1
                    if a is not None:
                        wnd, cov, oc = a[basin + "_wind"], a[basin + "_cov"], a[basin + "_oc"]
                        ok = (oc >= 8) & (cov >= 0.5 * oc)
                        for T in (48, 64):
                            k = (m >= 2 * T).sum(0)[ok]; y = (wnd >= 2 * T)[ok]
                            n_, y_ = counts_from(k, y, N); rows.append((n_, y_)); keys.append((model, "ascat", basin, T, ymd, lead, N))
                    # truth 2
                    cr, cc, idx, LA, LO = neighbours(basin)
                    ny, nx = LA.shape; flat = m.reshape(N, -1)
                    nbmax = np.stack([flat[:, ix].max(1) for ix in idx], 1)       # (N, centres)
                    clat, clon = LA[cr, cc], LO[cr, cc]
                    oc = (a[basin + "_oc"][cr, cc] >= 8) if a is not None else np.ones(len(cr), bool)
                    # fixes of this basin only
                    fb = fx[fx.basin == basin]
                    for T in (48, 64):
                        kk = (nbmax >= 2 * T).sum(0)
                        f = fb[fb.Category.isin(CATS[T])]
                        if len(f):
                            d = gc(clat[:, None], clon[:, None], f.lat.values[None, :], f.lon.values[None, :])
                            y = (d <= RAD).any(1)
                        else: y = np.zeros(len(cr), bool)
                        n_, y_ = counts_from(kk[oc], y[oc], N); rows.append((n_, y_)); keys.append((model, "opc", basin, T, ymd, lead, N))
    return keys, rows


if __name__ == "__main__":
    model = sys.argv[1]; seasons = [int(x) for x in sys.argv[2].split(",")]; leads = [int(x) for x in sys.argv[3].split(",")]
    out = sys.argv[4]
    keys, rows = build(model, seasons, leads, opc_fixes())
    N = max(k[-1] for k in keys)
    n = np.zeros((len(keys), N + 1), np.int64); y = np.zeros_like(n)
    for i, (a, b) in enumerate(rows): n[i, :len(a)] = a; y[i, :len(b)] = b
    np.savez_compressed(out, keys=np.array([list(map(str, k)) for k in keys]), n=n, y=y)
    print(len(keys), "case groups ->", out)
