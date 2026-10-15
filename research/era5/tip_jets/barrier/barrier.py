"""Barrier add-on (ERA5 PROXY): Denmark Strait pressure difference vs hurricane force in the 73 northerly-class tip jet episodes.
Plan: PREREGISTRATION.md (committed before the pull). usage: barrier.py pull | barrier.py analyse
Pull: MSLP only, one ARCO chunk per episode time (about 2.5 MB each), cached to $ERA5_WORK/tip_jets/barrier/msl_<YYYYMMDDHH>.npy.
"""
import os, sys, time, urllib.request
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EP = os.path.join(HERE, "..", "results", "episodes.csv")
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "tip_jets", "barrier")
ARCO = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
T0 = np.datetime64("1900-01-01T00")
ICE, GRN = (64.5, -18.0), (65.6, -37.6)       # Iceland point, Greenland east-coast point (lat, lon east)


def fetch(ti):
    import numcodecs
    codec = numcodecs.Blosc()
    url = f"{ARCO}/mean_sea_level_pressure/{ti}.0.0"
    for k in range(7):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                raw = r.read()
            return np.frombuffer(codec.decode(raw), "<f4").reshape(721, 1440)
        except Exception:
            if k == 6:
                raise
            time.sleep(2 ** k)


def cellmean(a, lat, lon):
    r = int(round((90 - lat) / 0.25)); c = int(round((lon % 360) / 0.25))
    return float(a[r - 1:r + 2, c - 1:c + 2].mean()) / 100.0


def episodes():
    e = pd.read_csv(EP, parse_dates=["time"])
    return e[e.cls == "E"].reset_index(drop=True)


def pull():
    os.makedirs(WORK, exist_ok=True)
    tot = 0
    for t in episodes().time:
        f = os.path.join(WORK, t.strftime("%Y%m%d%H") + ".npy")
        if os.path.exists(f):
            continue
        ti = int((np.datetime64(t.strftime("%Y-%m-%dT%H")) - T0) / np.timedelta64(1, "h"))
        a = fetch(ti)
        np.save(f, np.array([cellmean(a, *ICE), cellmean(a, *GRN)]))
        tot += 1
    print("pulled", tot, "times")


def analyse():
    from scipy.optimize import minimize
    from scipy.stats import spearmanr
    e = episodes()
    vals = np.array([np.load(os.path.join(WORK, t.strftime("%Y%m%d%H") + ".npy")) for t in e.time])
    e["DP"] = vals[:, 0] - vals[:, 1]
    e[["time", "season", "hf", "box_gust_kt", "LOW_DEPTH", "GH", "DP"]].to_csv(os.path.join(HERE, "barrier_values.csv"), index=False)
    z = lambda x: (x - x.mean()) / x.std()
    def fit(X, y):
        X1 = np.column_stack([np.ones(len(y)), X])
        nll = lambda b: np.sum(np.logaddexp(0, X1 @ b) - y * (X1 @ b)) + 0.5 * 1e-3 * np.sum(b[1:] ** 2)
        return minimize(nll, np.zeros(X1.shape[1]), method="BFGS").x
    y = e.hf.astype(float).values
    rng = np.random.default_rng(248)
    seasons = e.season.unique()
    out = []
    def run(cols, label):
        X = np.column_stack([z(e[c]).values for c in cols])
        b = fit(X, y)[1:]
        bs = []
        for _ in range(2000):
            pick = rng.choice(seasons, len(seasons))
            idx = np.concatenate([np.where(e.season.values == s)[0] for s in pick])
            if y[idx].sum() < 3 or y[idx].sum() > len(idx) - 3:
                continue
            Xb = np.column_stack([z(e.iloc[idx][c]).values for c in cols])
            bs.append(fit(Xb, y[idx])[1:])
        bs = np.array(bs); lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
        for i, c in enumerate(cols):
            out.append(f"{label:34s} {c:10s} OR/SD {np.exp(b[i]):5.2f}  [{np.exp(lo[i]):5.2f}, {np.exp(hi[i]):5.2f}]  (n_boot {len(bs)})")
    run(["DP", "LOW_DEPTH"], "PRIMARY hf ~ DP + depth")
    run(["DP", "LOW_DEPTH", "GH"], "secondary: + ice-sheet pressure")
    r, p = spearmanr(e.DP, e.box_gust_kt)
    out.append(f"SECONDARY spearman(DP, peak box gust) = {r:.2f} (p {p:.3f}, episodes not independent)")
    out.append(f"spearman(DP, GH) = {spearmanr(e.DP, e.GH)[0]:.2f}; spearman(DP, depth) = {spearmanr(e.DP, e.LOW_DEPTH)[0]:.2f}")
    out.append(f"DP mean/SD {e.DP.mean():.1f}/{e.DP.std():.1f} hPa; hf mean DP {e.DP[e.hf].mean():.1f}, storm-force mean DP {e.DP[~e.hf].mean():.1f}; n {len(e)}, hf {int(e.hf.sum())}")
    txt = "\n".join(out); print(txt)
    open(os.path.join(HERE, "barrier_result.txt"), "w").write("ERA5 PROXY, post hoc; plan: PREREGISTRATION.md\n" + txt + "\n")


if __name__ == "__main__":
    {"pull": pull, "analyse": analyse}[sys.argv[1]]()
