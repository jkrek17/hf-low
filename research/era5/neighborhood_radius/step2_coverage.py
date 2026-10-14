"""Step 2 (PREREGISTRATION.md): coverage of shading by radius, position error only. ERA5 PROXY gust (cache), Oct-Apr, 22 seasons.
Reads /mnt/project-files/era5-cache/gust_nh_12h/*.npz (Blosc+bitshuffle+zstd; decoded here with node's zlib zstd, no extra libraries).
Usage: python3 step2_coverage.py [first_file last_file]   -> results/season_sums_<tag>.csv"""
import os, sys, subprocess, tempfile, numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
CACHE = "/mnt/project-files/era5-cache"; BLOSC_JS = os.path.join(HERE, "blosc_decode.js")
RADII = [30, 60, 120, 270, 432]                       # nautical miles
CLASSES = {"gale": 144 * 38 / 72, "storm": 108, "hf": 144}   # stored value (0.5 kt steps): 76, 108, 144 (HF = 72.0 kt)
CLASSES["gale"] = 76
LEADS = [0, 24, 48, 72, 96]
NDRAW = {"gale": 1, "storm": 1, "hf": 2}              # deviation logged in the PREREG: one draw for gale and storm
CELL = 0.25; NM_PER_DEG = 60.0
ROW0, ROW1 = 60, 281                                   # 75N .. 20N (row i = 90N - 0.25 i)
DOM = {"atl": np.r_[1040:1440, 0:41], "pac": np.arange(480, 961)}
BAND_ROWS = 20; KROWS = 29                             # 5 degree bands; 432 nm = 7.2 deg = 29 rows
lsm = np.load(os.path.join(CACHE, "lsm_nh_0p25.npy"))
OCEAN = {b: (lsm[ROW0:ROW1][:, cols] < 0.5) for b, cols in DOM.items()}
VEC = pd.read_csv(os.path.join(RES, "displacement_vectors_primary.csv"))

def decode(path):
    z = np.load(path); shape = tuple(int(x) for x in z["shape"])
    with tempfile.TemporaryDirectory() as td:
        a, b = os.path.join(td, "f.blosc"), os.path.join(td, "f.raw")
        open(a, "wb").write(z["data"].tobytes())
        subprocess.run(["node", BLOSC_JS, a, b], check=True, stdout=subprocess.DEVNULL)
        raw = np.fromfile(b, np.uint8)
    out = np.empty_like(raw); bs = 262144
    for s in range(0, len(raw), bs):
        blk = raw[s:s + bs]; k = (len(blk) // 8) * 8
        if k:
            bits = np.unpackbits(blk[:k].reshape(8, k // 8), axis=1, bitorder="little")
            v = np.zeros(k, np.uint8)
            for j in range(8): v |= (bits[j] << j).astype(np.uint8)
            out[s:s + k] = v
        out[s + k:s + len(blk)] = blk[k:]
    return z["hours_since_1900"], out.reshape(shape)

_KF = {}
def kernel_fft(width, latc):
    key = (width, round(latc))
    if key not in _KF:
        dx = 15.0 * np.cos(np.radians(latc)); dy = 15.0
        H = BAND_ROWS + 2 * KROWS; ks = []
        for R in RADII:
            ky = int(R // dy); kx = int(R // dx)
            yy = np.arange(-ky, ky + 1)[:, None] * dy; xx = np.arange(-kx, kx + 1)[None, :] * dx
            d = ((yy ** 2 + xx ** 2) <= R ** 2).astype(np.float32)
            K = np.zeros((H, width + 2 * 112), np.float32)   # centred kernel placed at origin with wrap
            ky2, kx2 = d.shape[0] // 2, d.shape[1] // 2
            K[:d.shape[0], :d.shape[1]] = d
            K = np.roll(K, (-ky2, -kx2), axis=(0, 1))
            ks.append(np.fft.rfft2(K))
        _KF[key] = ks
    return _KF[key]

_TPL = {}
def templates(nr):
    """Per source row: nm distance from the source cell to each cell in a (2*KROWS+1) x (2*kc+1) window."""
    if nr not in _TPL:
        t = []
        for r in range(nr):
            lat = 75 - CELL * r; dx = 15.0 * np.cos(np.radians(lat)); kc = int(432 // dx) + 1
            yy = np.arange(-KROWS, KROWS + 1)[:, None] * 15.0; xx = np.arange(-kc, kc + 1)[None, :] * dx
            t.append((kc, np.sqrt(yy ** 2 + xx ** 2).astype(np.float32)))
        _TPL[nr] = t
    return _TPL[nr]

def dilations_sparse(B):
    nr, nc = B.shape; md = np.full((nr, nc), 1e9, np.float32); tp = templates(nr)
    for r, c in zip(*np.nonzero(B)):
        kc, D = tp[r]
        r0, r1 = max(0, r - KROWS), min(nr, r + KROWS + 1); c0, c1 = max(0, c - kc), min(nc, c + kc + 1)
        w = md[r0:r1, c0:c1]
        np.minimum(w, D[r0 - r + KROWS:r1 - r + KROWS, c0 - c + kc:c1 - c + kc], out=w)
    return [md <= R for R in RADII]

def dilations(B):
    return dilations_sparse(B)
    return dilations_fft(B)

def dilations_fft(B):
    """B: bool (rows, cols). Returns list per radius of bool arrays: cells within R nm of any True cell."""
    nr, nc = B.shape; out = [np.zeros_like(B) for _ in RADII]
    if not B.any(): return out
    for r0 in range(0, nr, BAND_ROWS):
        r1 = min(r0 + BAND_ROWS, nr); lo = max(0, r0 - KROWS); hi = min(nr, r1 + KROWS)
        if not B[lo:hi].any(): continue
        latc = 75 - CELL * (r0 + r1) / 2
        ks = kernel_fft(nc, latc)
        H = BAND_ROWS + 2 * KROWS; W = nc + 224
        pad = np.zeros((H, W), np.float32)
        pad[(lo - (r0 - KROWS)):(lo - (r0 - KROWS)) + (hi - lo), 112:112 + nc] = B[lo:hi]
        Fp = np.fft.rfft2(pad)
        for i, K in enumerate(ks):
            c = np.fft.irfft2(Fp * K, s=(H, W))
            out[i][r0:r1] = c[KROWS:KROWS + (r1 - r0), 112:112 + nc] > 0.5
    return out

def shift(B, dr, dc):
    out = np.zeros_like(B); nr, nc = B.shape
    if abs(dr) >= nr or abs(dc) >= nc: return out
    rs, re = max(0, dr), min(nr, nr + dr); cs, ce = max(0, dc), min(nc, nc + dc)
    out[rs:re, cs:ce] = B[rs - dr:re - dr, cs - dc:ce - dc]
    return out

def vec_for(basin, lead):
    v = VEC[(VEC.lead == lead) & (VEC.basin == basin)]
    if len(v) < 100: v = VEC[VEC.lead == lead]
    return v[["dlat", "dlon"]].to_numpy()

def do_file(args):
    fi, path = args
    hrs, data = decode(path)
    t = pd.to_datetime("1900-01-01") + pd.to_timedelta(hrs.astype(np.int64), unit="h")
    rng = np.random.default_rng(1000 + fi)
    acc = {}
    for ti, tt in enumerate(t):
        if tt.month not in (10, 11, 12, 1, 2, 3, 4): continue
        season = tt.year if tt.month >= 10 else tt.year - 1
        for basin, cols in DOM.items():
            sub = data[ti, ROW0:ROW1][:, cols]; oc = OCEAN[basin]
            for cname, thr in CLASSES.items():
                T = (sub >= thr) & oc
                if not T.any(): continue
                dT = dilations(T)
                for L in LEADS:
                    V = vec_for(basin, L)
                    for _ in range(NDRAW[cname]):
                        dlat, dlon = V[rng.integers(len(V))]
                        F = shift(T, -int(round(dlat / CELL)), int(round(dlon / CELL))) & oc
                        nF = int(F.sum()); nT = int(T.sum())
                        dF = dilations(F) if nF else [np.zeros_like(F) for _ in RADII]
                        for i, R in enumerate(RADII):
                            a = acc.setdefault((season, basin, cname, L, R), np.zeros(6))
                            a += [nT, (T & dF[i]).sum(), (dF[i] & oc).sum(), nF, (F & dT[i]).sum(), 1]
    return acc

if __name__ == "__main__":
    files = sorted(f for f in os.listdir(os.path.join(CACHE, "gust_nh_12h")) if f.endswith(".npz"))
    a, b = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (0, len(files))
    tag = f"{a}_{b}"
    jobs = [(i, os.path.join(CACHE, "gust_nh_12h", files[i])) for i in range(a, b)]
    tot = {}
    with Pool(int(os.environ.get("NPROC", 4))) as p:
        for k, acc in enumerate(p.imap_unordered(do_file, jobs)):
            for key, v in acc.items(): tot.setdefault(key, np.zeros(6)); tot[key] += v
            print("done", k + 1, "of", len(jobs), flush=True)
    rows = [dict(season=k[0], basin=k[1], cls=k[2], lead=k[3], R_nm=k[4], T=v[0], T_covered=v[1], shaded=v[2], F=v[3], F_near_T=v[4], draws=v[5]) for k, v in sorted(tot.items())]
    pd.DataFrame(rows).to_csv(os.path.join(RES, f"season_sums_{tag}.csv"), index=False)
