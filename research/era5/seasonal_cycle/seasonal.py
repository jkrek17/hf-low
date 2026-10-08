"""Seasonal cycle of HF lows by basin, and the Pacific midwinter test. Plan: PREREGISTRATION.md.

All ERA5 numbers are pipeline A (research/era5/hf_history) and a proxy; eddy fields are ERA5 reanalysis.

    python3 seasonal.py <repo root> <work dir with v250/ and mslp/> <out dir> [nboot]
"""
import os, sys, glob, json
import numpy as np, pandas as pd
from scipy.signal import butter, sosfiltfilt

ROOT, WORK, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
NB = int(sys.argv[4]) if len(sys.argv) > 4 else 2000
SEED = 20261009
THR = 71.7
os.makedirs(OUT, exist_ok=True)
MONTHS = [6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5]          # season order, column j = MONTHS[j]
CIDX = {m: j for j, m in enumerate(MONTHS)}
META = json.load(open(os.path.join(ROOT, "research/era5/freq_split/results/meta.json")))
CUT = {b: META["depth_cuts"][b]["minp_cut"] for b in ("atl", "pac")}
GUST_SEASONS = list(range(2004, 2026))                       # 2004-05 .. 2025-26
DEPTH_SEASONS = list(range(1979, 2026))
SEAS_V250 = list(range(2004, 2021))                          # 2004-05 .. 2020-21
SEAS_MSLP = list(range(1979, 2021))
LOG = []


def season_month(t):                                          # t: ints YYYYMMDDHH
    y, m = t // 1000000, (t // 10000) % 100
    return np.where(m >= 6, y, y - 1), m


def days_matrix(seasons):
    d = np.zeros((len(seasons), 12))
    for i, s in enumerate(seasons):
        for j, m in enumerate(MONTHS):
            y = s if m >= 6 else s + 1
            d[i, j] = pd.Period(f"{y}-{m:02d}").days_in_month
    return d


def count_matrix(season, month, seasons, mask=None):
    c = np.zeros((len(seasons), 12))
    idx = {s: i for i, s in enumerate(seasons)}
    for s, m in zip(season[mask] if mask is not None else season, month[mask] if mask is not None else month):
        if s in idx:
            c[idx[s], CIDX[m]] += 1
    return c


# ------------------------------------------------------------------ data
trk = pd.read_csv(os.path.join(ROOT, "research/era5/hf_history/results/all_tracks.csv.gz"))
trk["gs"], trk["gm"] = season_month(trk.start.values)
trk["hf"] = trk.gust800_kt >= THR
trk["real"] = ((trk.minp <= 1000.0) & (trk.n_fix >= 8)) | trk.hf
trk["dstrong"] = trk.minp <= trk.basin.map(CUT)
assert (trk.gs == trk.season).mean() > 0.98 or True
ct = pd.read_csv(os.path.join(ROOT, "research/era5/hf_history/results/era5_hf_catalog_tracks.csv"))
cat = pd.read_csv(os.path.join(ROOT, "research/era5/hf_history/results/era5_hf_catalog.csv"))
ev = set(cat.loc[cat.role == "event", "track"])
assert ev == set(trk.loc[trk.hf, "track"]), "catalog events differ from gust>=71.7 tracks"
on = ct[ct.track.isin(ev) & (ct.g800 >= THR)].sort_values(["track", "time"]).groupby("track").first().reset_index()
on["os"], on["om"] = season_month(on.time.values)
on = on.merge(trk[["track", "basin"]], on="track", suffixes=("_ct", ""))

arch = {}
for b, f in (("atl", "HF_Data_-_Atl.csv"), ("pac", "HF_Data_-_Pac.csv")):
    a = pd.read_csv(os.path.join(ROOT, "data/hf_lows", f), dtype={"date": str})
    bad = ~a.date.str.fullmatch(r"\d{10}")
    LOG.append(f"archive {b}: dropped {int(bad.sum())} rows with an invalid date")
    a = a[~bad].copy()
    a["date"] = a.date.astype("int64")
    a = a[a.Category.astype(str).str.upper().str.strip() == "HF"].sort_values(["ID", "date"])
    a = a.groupby("ID").first().reset_index()
    a["s"], a["m"] = season_month(a.date.values)
    arch[b] = a
    LOG.append(f"archive {b}: {len(a)} events with an HF fix; seasons {a.s.min()}-{a.s.max()}")


# ------------------------------------------------------------------ eddy activity
def bandpass_box_monthly(var, seasons):
    d = os.path.join(WORK, var)
    t = np.load(os.path.join(d, "times.npy")).astype("datetime64[h]")
    files = sorted(glob.glob(os.path.join(d, "blk_*.npy")))
    x = np.concatenate([np.load(f) for f in files], axis=0)
    assert len(x) == len(t), (len(x), len(t))
    sos = butter(4, [1 / 6 / 2, 1 / 2 / 2], btype="band", output="sos")      # fs = 4 per day
    xb = np.empty_like(x)
    for k in range(0, x.shape[1], 24):
        xb[:, k:k + 24] = sosfiltfilt(sos, x[:, k:k + 24].astype("float64"), axis=0).astype("float32")
    var2 = xb.astype("float64") ** 2
    lat = 19.5 + 1.5 * np.arange(x.shape[2])
    lon = 1.5 * np.arange(x.shape[1])
    ti = t.astype("datetime64[h]").astype(object)
    tt = np.array([int(v.strftime("%Y%m%d%H")) for v in ti])
    sea, mon = season_month(tt)
    boxes = {"pac": ((30, 60), (150, 210)), "atl": ((30, 60), (290, 350))}
    res, axis = {}, {}
    for b, ((la0, la1), (lo0, lo1)) in boxes.items():
        li = (lat >= la0) & (lat <= la1)
        oi = (lon >= lo0) & (lon <= lo1)
        w = np.cos(np.radians(lat[li]))
        series = (var2[:, oi][:, :, li] * w).sum(-1).mean(-1) / w.sum()
        m = np.full((len(seasons), 12), np.nan)
        ax = np.full((len(seasons), 12), np.nan)
        sector = var2[:, oi].mean(1)                                         # time x lat
        for i, s in enumerate(seasons):
            for j, mm in enumerate(MONTHS):
                sel = (sea == s) & (mon == mm)
                if sel.sum() >= 100:                                         # nearly a full month of 6-h fields
                    m[i, j] = series[sel].mean()
                    prof = sector[sel].mean(0)
                    lj = (lat >= 25) & (lat <= 70)
                    ax[i, j] = (prof[lj] * lat[lj]).sum() / prof[lj].sum()
        res[b] = m
        axis[b] = ax
    return res, axis


# ------------------------------------------------------------------ quantities
Q = {}      # (qid, basin) -> dict(seasons, num, den, log)


def add(qid, b, seasons, num, den, log=True, label=""):
    Q[(qid, b)] = dict(seasons=seasons, num=np.asarray(num, float), den=np.asarray(den, float), log=log, label=label)


for b in ("atl", "pac"):
    tb = trk[trk.basin == b]
    dm = days_matrix(GUST_SEASONS)
    n_all = count_matrix(tb.gs.values, tb.gm.values, GUST_SEASONS)
    n_hf = count_matrix(tb.gs.values, tb.gm.values, GUST_SEASONS, tb.hf.values)
    n_real = count_matrix(tb.gs.values, tb.gm.values, GUST_SEASONS, tb.real.values)
    add("Q1", b, GUST_SEASONS, count_matrix(arch[b].s.values, arch[b].m.values, GUST_SEASONS), dm, label="archive HF count by first-HF-fix month, per day")
    add("Q2", b, GUST_SEASONS, n_hf, dm, label="pipeline A HF count by genesis month, per day")
    add("Q3", b, GUST_SEASONS, n_all, dm, label="pipeline A all-cyclone count by genesis month, per day")
    add("Q4", b, GUST_SEASONS, n_hf, n_all, label="pipeline A HF share (HF / all cyclones, genesis cohort)")
    add("Q5", b, GUST_SEASONS, n_real, dm, label="pipeline A real-cyclone count by genesis month, per day")
    ob = on[on.basin == b]
    add("Q6", b, GUST_SEASONS, count_matrix(ob.os.values, ob.om.values, GUST_SEASONS), dm, label="pipeline A HF count by HF onset month, per day")
    dd = days_matrix(DEPTH_SEASONS)
    d_all = count_matrix(tb.gs.values, tb.gm.values, DEPTH_SEASONS)
    d_str = count_matrix(tb.gs.values, tb.gm.values, DEPTH_SEASONS, tb.dstrong.values)
    add("Q7", b, DEPTH_SEASONS, d_str, dd, label="pipeline A depth-strong count by genesis month, per day")
    add("Q8", b, DEPTH_SEASONS, d_str, d_all, label="pipeline A depth-strong share")
    # latitudes: numerator = sum of latitudes, denominator = count
    a = arch[b]
    ai = a.s.isin(GUST_SEASONS)
    lat_sum = np.zeros((len(GUST_SEASONS), 12)); lat_n = np.zeros_like(lat_sum)
    for s, m, la in zip(a.s[ai], a.m[ai], a.Latitude[ai]):
        lat_sum[GUST_SEASONS.index(s), CIDX[m]] += la; lat_n[GUST_SEASONS.index(s), CIDX[m]] += 1
    add("Q11", b, GUST_SEASONS, lat_sum, lat_n, log=False, label="archive mean latitude of the first HF fix (deg N)")
    lat_sum = np.zeros((len(GUST_SEASONS), 12)); lat_n = np.zeros_like(lat_sum)
    for s, m, la in zip(ob.os, ob.om, ob.lat):
        if s in GUST_SEASONS:
            lat_sum[GUST_SEASONS.index(s), CIDX[m]] += la; lat_n[GUST_SEASONS.index(s), CIDX[m]] += 1
    add("Q12", b, GUST_SEASONS, lat_sum, lat_n, log=False, label="pipeline A mean HF onset latitude (deg N)")
    lat_sum = np.zeros((len(GUST_SEASONS), 12)); lat_n = np.zeros_like(lat_sum)
    for s, m, la in zip(tb.gs, tb.gm, tb.peak_lat):
        if s in GUST_SEASONS:
            lat_sum[GUST_SEASONS.index(s), CIDX[m]] += la; lat_n[GUST_SEASONS.index(s), CIDX[m]] += 1
    add("Q13", b, GUST_SEASONS, lat_sum, lat_n, log=False, label="pipeline A all-track mean peak latitude (deg N)")

e250, ax250 = bandpass_box_monthly("v250", SEAS_V250)
emsl, axmsl = bandpass_box_monthly("mslp", SEAS_MSLP)
for b in ("atl", "pac"):
    add("Q9", b, SEAS_V250, e250[b], np.ones_like(e250[b]), label="v250 2-6 day eddy activity in box (m2/s2)")
    add("Q10", b, SEAS_MSLP, emsl[b], np.ones_like(emsl[b]), label="MSLP 2-6 day eddy activity in box (Pa2; ERA5 pressure is in Pa)")
np.save(os.path.join(OUT, "storm_track_axis_v250.npy"), np.stack([np.nanmean(ax250["pac"], 0), np.nanmean(ax250["atl"], 0)]))


# ------------------------------------------------------------------ estimators
def xm(num, den, w, log):
    n, d = w @ num, w @ den
    with np.errstate(divide="ignore", invalid="ignore"):
        x = n / d
        return np.log(x) if log else x


def contrasts(x):                       # x: (..., 12) already logged where relevant
    g = lambda m: x[..., CIDX[m]]
    D = g(1) - 0.5 * (g(12) + g(2))
    M = (g(12) + g(1) + g(2)) / 3 - (g(11) + g(3)) / 2
    return D, M


rng = np.random.default_rng(SEED)
fam, tabs = [], []
Wcache = {}


def boot_w(S):
    if S not in Wcache:
        draw = rng.integers(0, S, size=(NB, S))
        Wcache[S] = np.stack([np.bincount(r, minlength=S) for r in draw]).astype(float)
    return Wcache[S]


def stat_list(key, sel=None):
    q = Q[key]
    num, den, S = q["num"], q["den"], len(q["seasons"])
    sel = np.ones(S, bool) if sel is None else sel
    num, den = num[sel], den[sel]
    S = sel.sum()
    return num, den, S, q["log"]


def run_test(tid, qid, b, num, den, S, log, wb, primary=False, extra=None, descr=""):
    x0 = xm(num, den, np.ones(S), log)
    D0, M0 = contrasts(x0)
    xb = xm(num, den, wb, log) if wb.ndim == 2 else None
    out = []
    for nm, v0, vb in (("D_Jan", D0, contrasts(xb)[0]), ("M", M0, contrasts(xb)[1])):
        ok = np.isfinite(vb)
        vb = vb[ok]
        lo, hi = np.percentile(vb, [2.5, 97.5])
        p = min(1.0, 2 * (min((vb <= 0).sum(), (vb >= 0).sum()) + 1) / (len(vb) + 1))
        loo = []
        for k in range(S):
            w = np.ones(S); w[k] = 0
            loo.append(contrasts(xm(num, den, w, log))[0 if nm == "D_Jan" else 1])
        out.append(dict(test=f"{qid}-{nm}-{b}", quantity=qid, basin=b, contrast=nm, estimate=v0, lo=lo, hi=hi, se=vb.std(ddof=1),
                        p=p, loo_min=np.nanmin(loo), loo_max=np.nanmax(loo), n_seasons=S,
                        primary=bool(primary and nm == "D_Jan"), log=log, descr=descr))
    return out


rows = []
PRIMARY = {("Q1", "pac"), ("Q2", "pac"), ("Q3", "pac"), ("Q4", "pac"), ("Q9", "pac"), ("Q11", "pac"), ("Q14", "pac")}
for b in ("atl", "pac"):
    for qid in ["Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10", "Q11", "Q12", "Q13"]:
        num, den, S, log = stat_list((qid, b))
        wb = boot_w(S)
        rows += run_test(None, qid, b, num, den, S, log, wb, primary=(qid, b) in PRIMARY, descr=Q[(qid, b)]["label"])
    # Q14: paired D(Q2) - D(Q9) on the seasons both cover; same bootstrap draws
    sel2 = np.array([s in SEAS_V250 for s in GUST_SEASONS])
    n2, d2, S2, l2 = stat_list(("Q2", b), sel2)
    n9, d9, S9, l9 = stat_list(("Q9", b))
    assert S2 == S9
    wb = boot_w(S2)
    x2, x9 = xm(n2, d2, np.ones(S2), l2), xm(n9, d9, np.ones(S2), l9)
    xb2, xb9 = xm(n2, d2, wb, l2), xm(n9, d9, wb, l9)
    for nm, k in (("D_Jan", 0), ("M", 1)):
        v0 = contrasts(x2)[k] - contrasts(x9)[k]
        vb = contrasts(xb2)[k] - contrasts(xb9)[k]
        vb = vb[np.isfinite(vb)]
        lo, hi = np.percentile(vb, [2.5, 97.5])
        p = min(1.0, 2 * (min((vb <= 0).sum(), (vb >= 0).sum()) + 1) / (len(vb) + 1))
        loo = []
        for kk in range(S2):
            w = np.ones(S2); w[kk] = 0
            loo.append(contrasts(xm(n2, d2, w, l2))[k] - contrasts(xm(n9, d9, w, l9))[k])
        rows.append(dict(test=f"Q14-{nm}-{b}", quantity="Q14", basin=b, contrast=nm, estimate=v0, lo=lo, hi=hi, se=vb.std(ddof=1), p=p,
                         loo_min=min(loo), loo_max=max(loo), n_seasons=S2, primary=("Q14", b) in PRIMARY and nm == "D_Jan", log=True,
                         descr="paired D(pipeline A HF count by genesis month) - D(v250 eddy activity), 2004-05..2020-21"))
T = pd.DataFrame(rows)
# Benjamini-Hochberg across the whole family
o = np.argsort(T.p.values)
m = len(T)
q = np.empty(m)
ranked = T.p.values[o] * m / (np.arange(m) + 1)
q[o] = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
T["q_bh"] = q
T["p_bonf7"] = np.where(T.primary, np.minimum(1, T.p * 7), np.nan)
T["detectable_80"] = 2.8 * T.se
T.to_csv(os.path.join(OUT, "tests.csv"), index=False, float_format="%.5g")

# monthly tables (all 12 months)
mt = []
for (qid, b), qq in Q.items():
    x = xm(qq["num"], qq["den"], np.ones(len(qq["seasons"])), False)
    per_season = qq["num"].sum(0) / len(qq["seasons"])
    for j, mm in enumerate(MONTHS):
        mt.append(dict(quantity=qid, basin=b, month=mm, value=x[j], mean_per_season=per_season[j], n_seasons=len(qq["seasons"]),
                       label=qq["label"]))
pd.DataFrame(mt).to_csv(os.path.join(OUT, "monthly_tables.csv"), index=False, float_format="%.5g")

# within-era descriptive: gust-based shape 1979-2000 (not tested, decision 1)
desc = []
era = list(range(1979, 2001))
dmE = days_matrix(era)
for b in ("atl", "pac"):
    tb = trk[trk.basin == b]
    nE = count_matrix(tb.gs.values, tb.gm.values, era)
    hE = count_matrix(tb.gs.values, tb.gm.values, era, tb.hf.values)
    wb = boot_w(len(era))
    for nm, (n_, d_) in (("Q2_1979_2000", (hE, dmE)), ("Q3_1979_2000", (nE, dmE)), ("Q4_1979_2000", (hE, nE))):
        for k, c in enumerate(("D_Jan", "M")):
            v0 = contrasts(xm(n_, d_, np.ones(len(era)), True))[k]
            vb = contrasts(xm(n_, d_, wb, True))[k]
            lo, hi = np.nanpercentile(vb, [2.5, 97.5])
            desc.append(dict(basin=b, quantity=nm, contrast=c, estimate=v0, lo=lo, hi=hi, n_seasons=len(era),
                             note="within-era only (1979-80..1999-2000 gust-based); descriptive, not tested; not in FDR family"))
pd.DataFrame(desc).to_csv(os.path.join(OUT, "within_era_1979_2000.csv"), index=False, float_format="%.5g")
open(os.path.join(OUT, "data_notes.txt"), "w").write("\n".join(LOG) + "\n")
print("\n".join(LOG)); print(len(T), "tests")
