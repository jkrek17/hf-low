"""RA-20 Arm A extraction: the Tier 2 groups of PR 56's plan (stability, trough, omega700) at the stratum fixes at the times in
results/times_armB.csv that WeatherBench2 covers (to 2023-01-09). Same grid, discs and weighting as intensity_extra/extract.py;
only temperature, geopotential and vertical velocity are read (26.96 GB). Ocean = land-sea mask < 0.5 (no SST/sea-ice term,
because SST is not read). ERA5 proxy, 1.5 degree. One CSV per WB2 chunk under $ERA5_WORK/explosive_onset/armA/, resumable.

Run: python3 research/era5/explosive_onset/extract_armA.py <repo_root> [n_proc]
"""
import os, sys, time, numpy as np, pandas as pd, numcodecs
from multiprocessing import Pool
root = sys.argv[1]
HERE = os.path.join(root, "research/era5/explosive_onset")
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "armA")
os.makedirs(WORK, exist_ok=True)
sys.path.insert(0, os.path.join(root, "research/era5/intensity"))
import env, hart
from env import LAT2, LON2, R0, R1, WB2, WB2_T0, WB2_END, WB2_LEVELS, cat, wb2meta, wmean
VARS = {"temperature": (925, 700), "geopotential": (500,), "vertical_velocity": (700,)}
_lsm = None


def lsm():
    global _lsm
    if _lsm is None:
        m = wb2meta()["land_sea_mask/.zarray"]
        a = np.frombuffer(numcodecs.get_codec(m["compressor"]).decode(cat(f"{WB2}/land_sea_mask/0.0")), m["dtype"]).reshape(m["chunks"])
        _lsm = a.T[R0:R1].astype(np.float64)
    return _lsm


def read_chunk(chunk):
    m = wb2meta(); out = {}; nb = 0
    for var, levs in VARS.items():
        za = m[f"{var}/.zarray"]; codec = numcodecs.get_codec(za["compressor"])
        raw = cat(f"{WB2}/{var}/{chunk}.0.0.0"); nb += len(raw)
        a = np.frombuffer(codec.decode(raw), za["dtype"]).reshape(za["chunks"])
        a = a[:min(za["chunks"][0], za["shape"][0] - chunk * za["chunks"][0])]
        for L in levs:
            out[f"{var}{L}"] = a[:, WB2_LEVELS.index(L)].transpose(0, 2, 1)[:, R0:R1]
    return out, nb


def features(F, la, lo):
    d = hart.great_circle_km(LAT2, LON2, la, lo)
    m5, m10 = d <= 500, d <= 1000
    ocean = lsm() < 0.5
    th925 = F["temperature925"] * (1000 / 925) ** 0.2857
    th700 = F["temperature700"] * (1000 / 700) ** 0.2857
    z = F["geopotential500"] / env.G
    ann = (d > 1500) & (d <= 2500)
    return dict(stab=wmean((th700 - th925) / 225.0, m5 & ocean),
                trough=wmean(z, ann) - float(np.nanmin(np.where(m10, z, np.nan))),
                omega700=float(np.nanmin(np.where(m5, F["vertical_velocity700"], np.nan))))


def run_block(job):
    key, fixes = job
    path = os.path.join(WORK, f"wb2_{key:05d}.csv")
    if os.path.exists(path):
        return key, 0
    Fall, nb = read_chunk(key)
    t0 = WB2_T0 + np.timedelta64(key * 8 * 6, "h"); rows = []
    for i in range(Fall["temperature925"].shape[0]):
        t = t0 + np.timedelta64(6 * i, "h"); stamp = int(str(t).replace("-", "").replace("T", "")[:10])
        sub = fixes[fixes.time == stamp]
        if not len(sub): continue
        F = {k: v[i].astype(np.float64) for k, v in Fall.items()}
        for r in sub.itertuples():
            f = features(F, r.lat, r.lon); f.update(track=int(r.track), time=stamp); rows.append(f)
    tmp = path + ".tmp"; pd.DataFrame(rows).to_csv(tmp, index=False); os.replace(tmp, path)
    with open(os.path.join(WORK, "bytes.log"), "a") as f: f.write(f"{key},{nb}\n")
    return key, nb


if __name__ == "__main__":
    F = pd.read_csv(os.path.join(root, "research/era5/intensity/results/fixes_2004.csv.gz"))
    F = F[(F.g800 < 55) & (~F.hf_now)]
    T = pd.read_csv(os.path.join(HERE, "results/times_armB.csv"))
    F = F[F.time.isin(T.time)]
    t = pd.to_datetime(F.time.astype(str), format="%Y%m%d%H").values.astype("datetime64[h]")
    F = F[t < WB2_END.astype("datetime64[h]")].copy(); t = t[t < WB2_END.astype("datetime64[h]")]
    F["key"] = ((t - WB2_T0) / np.timedelta64(6, "h")).astype(int) // 8
    J = [(int(k), g.drop(columns="key")) for k, g in F.groupby("key")]
    print(len(F), "fixes", len(J), "chunks", flush=True)
    if "--dry" in sys.argv: sys.exit()
    n = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 6
    t0, tot = time.time(), 0
    with Pool(n) as p:
        for i, (k, nb) in enumerate(p.imap_unordered(run_block, J)):
            tot += nb
            if i % 50 == 0: print(f"{i + 1}/{len(J)} {tot / 1e9:.2f} GB {time.time() - t0:.0f}s", flush=True)
    print(f"done {len(J)} chunks, {tot / 1e9:.2f} GB this run")
