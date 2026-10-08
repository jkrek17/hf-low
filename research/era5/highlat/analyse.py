"""High-latitude Atlantic HF lows: the tests planned in PLAN.md.

ERA5 numbers are pipeline A (research/era5/hf_history) and are a proxy.
Seasons 2004-05 on (June-May, labelled by the June year). Inference unit is
the season: every interval is a season block bootstrap (2,000 draws, seed 1).

Inputs
  docs/data/hf-lows.json                      archive payload (built by tools/build_hf_lows.py)
  ../hf_history/results/era5_hf_catalog*.csv  pipeline A events and their track points
  work/gustloc/*.csv                          from gustloc.py (ignored; about 3,100 small files)
Outputs
  gustloc_fixes.csv   one row per Atlantic HF-strength fix (small, committed)
  results.txt         every number quoted from this study
"""
import os, sys, glob, json
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RES = os.path.join(HERE, "..", "hf_history", "results")
WORK = os.environ.get("ERA5_WORK", os.path.join(HERE, "work"))
THR = 71.7
S0 = 2004
NB = 2000
RE = 6371.0
rng = np.random.default_rng(1)
out = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)


def gc(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * RE * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def region(lat, lon):
    """lon in -180..180. Regions fixed in PLAN.md before results were seen."""
    if lat < 60:
        return "S<60N"
    if lon < -50:
        return "Labrador/Davis"
    if lon < -30:
        return "CapeFarewell/Irminger"
    if lon < -10:
        return "DenmarkStr/Iceland"
    return "Norwegian"


def boot(seasons, stat, nb=NB):
    """Season block bootstrap of stat(list_of_seasons) -> (2.5, 97.5) percentiles."""
    u = np.array(sorted(set(seasons)))
    vals = []
    for _ in range(nb):
        v = stat(rng.choice(u, len(u), replace=True))
        if v == v:
            vals.append(v)
    return np.percentile(vals, [2.5, 97.5])


def fmt(ci, d=3):
    return f"[{ci[0]:.{d}f}, {ci[1]:.{d}f}]"


# ---------------------------------------------------------------- archive
d = json.load(open(os.path.join(REPO, "docs", "data", "hf-lows.json")))
LOWS = [dict(zip(d["lowFields"], r)) for r in d["lows"]]
ARC = [l for l in LOWS if l["season"] >= S0]

say("HIGH-LATITUDE ATLANTIC HF LOWS  (archive 2004-05..2025-26; ERA5 pipeline A proxy)")
say("=" * 78)
say("\n1. Share of HF lows with an HF fix at or north of 60N (archive)")
say("   The archive has no Atlantic fix north of %.1fN: OPC's high-seas area ends at 67N." %
    max(f[1] for l in LOWS if l["basin"] == "atl" for f in l["fixes"]))
for b in ("atl", "pac"):
    E = [l for l in ARC if l["basin"] == b and any(f[3] == "HF" for f in l["fixes"])]
    hi = [any(f[3] == "HF" and f[1] >= 60 for f in l["fixes"]) for l in E]
    s = np.array([l["season"] for l in E]); h = np.array(hi)
    ci = boot(s, lambda ss: np.mean(np.concatenate([h[s == x] for x in ss])))
    say(f"   {b}: {sum(hi)} of {len(E)} = {100 * np.mean(hi):.1f}%  95% CI {fmt(100 * ci, 1)}")
EA = [l for l in ARC if l["basin"] == "atl" and any(f[3] == "HF" for f in l["fixes"])]
for c in ("low", "tipjet", "nocentre"):
    sub = [l for l in EA if l["cls"] == c]
    say(f"   atl class {c:9s}: {len(sub):4d} events, {sum(any(f[3] == 'HF' and f[1] >= 60 for f in l['fixes']) for l in sub):4d} with an HF fix >= 60N")

# 3. archive pressure by region, one value per event per region
say("\n2. Archive central pressure at HF fixes by region (one median per event per region)")
rows = []
for l in EA:
    if l["cls"] != "low":
        continue
    per = {}
    for t, la, lo, cat, p in l["fixes"]:
        if cat == "HF" and p is not None:
            per.setdefault(region(la, lo), []).append(p)
    for r, ps in per.items():
        rows.append((l["season"], r, float(np.median(ps)), l["id"]))
P = pd.DataFrame(rows, columns=["season", "region", "pres", "id"])
south = P[P.region == "S<60N"]
for r in ["S<60N", "Labrador/Davis", "CapeFarewell/Irminger", "DenmarkStr/Iceland", "Norwegian"]:
    x = P[P.region == r]
    if len(x) == 0:
        say(f"   {r:22s} n=0")
        continue
    line = f"   {r:22s} n={len(x):4d} events  median {x.pres.median():6.1f} hPa  IQR {x.pres.quantile(.25):.0f}-{x.pres.quantile(.75):.0f}"
    if r != "S<60N" and len(x) >= 10:
        def st(ss, x=x):
            a = pd.concat([x[x.season == q] for q in ss]); b = pd.concat([south[south.season == q] for q in ss])
            return a.pres.median() - b.pres.median() if len(a) and len(b) else np.nan
        line += f"  minus south {x.pres.median() - south.pres.median():+.1f}  95% CI {fmt(boot(P.season, st), 1)}"
    say(line)
say("   (Positive = weaker storms at HF there. Background: the Iceland Low makes high-latitude")
say("    pressures lower for the same storm, so a positive difference is conservative.)")

# ---------------------------------------------------------------- pipeline A vs archive
say("\n3. Pipeline A against the archive by region (2004-05 on; Atlantic; any archive low)")
cat = pd.read_csv(os.path.join(RES, "era5_hf_catalog.csv"))
trk = pd.read_csv(os.path.join(RES, "era5_hf_catalog_tracks.csv"))
ev = cat[(cat.role == "event") & (cat.basin == "atl") & (cat.season >= S0)].copy()
tp = trk[trk.track.isin(ev.track)].copy()
tp["lon180"] = np.where(tp.lon > 180, tp.lon - 360, tp.lon)
by_time = {t: g for t, g in tp.groupby("time")}
AA = [l for l in ARC if l["basin"] == "atl"]
hit_tracks = {}
arc_rows = []
for l in AA:
    rad = 400 if l["cls"] == "low" else 800
    m = set()
    for t, la, lo, cat_, p in l["fixes"]:
        g = by_time.get(t)
        if g is None:
            continue
        dd = gc(la, lo, g.lat.values, g.lon180.values)
        m.update(g.track.values[dd <= rad])
    for k in m:
        hit_tracks.setdefault(k, []).append(l["id"])
    hi = any(f[1] >= 60 for f in l["fixes"])
    hihf = any(f[1] >= 60 and f[3] == "HF" for f in l["fixes"])
    arc_rows.append((l["season"], l["id"], l["cls"], hi, hihf, len(m) > 0))
AR = pd.DataFrame(arc_rows, columns=["season", "id", "cls", "hi", "hihf", "hit"])
ev["hit"] = ev.track.isin(hit_tracks)
ev["hi"] = ev.peak_lat >= 60
tp_in = tp[tp.basin == "atl"]
ev["anyhi"] = ev.track.isin(tp_in[(tp_in.lat >= 60) & (tp_in.g800 >= THR)].track)


def podfar(sub_a, sub_e):
    pod = sub_a.hit.mean(); far = 1 - sub_e.hit.mean()
    return pod, far


say(f"   all:  archive lows {len(AR)}, POD {AR.hit.mean():.3f};  ERA5 events {len(ev)}, FAR {1 - ev.hit.mean():.3f};"
    f"  bias (ERA5/archive) {len(ev) / len(AR):.2f}")
for lab, a, e in (("archive fix >= 60N / ERA5 peak >= 60N", AR[AR.hi], ev[ev.hi]),
                  ("archive fix < 60N only / ERA5 peak < 60N", AR[~AR.hi], ev[~ev.hi]),
                  ("ERA5 has an HF-strength fix >= 60N", None, ev[ev.anyhi]),
                  ("ERA5 has none", None, ev[~ev.anyhi])):
    s = f"   {lab:42s}"
    if a is not None:
        s += f" POD {a.hit.mean():.3f} (n={len(a)})"
    s += f"  FAR {1 - e.hit.mean():.3f} (n={len(e)})"
    if a is not None:
        s += f"  ERA5/archive {len(e) / len(a):.2f}"
    say(s)


def fardiff(ss):
    a = pd.concat([ev[(ev.season == q) & ev.hi] for q in ss]); b = pd.concat([ev[(ev.season == q) & ~ev.hi] for q in ss])
    return (1 - a.hit.mean()) - (1 - b.hit.mean())


def poddiff(ss):
    a = pd.concat([AR[(AR.season == q) & AR.hi] for q in ss]); b = pd.concat([AR[(AR.season == q) & ~AR.hi] for q in ss])
    return a.hit.mean() - b.hit.mean()


fd = (1 - ev[ev.hi].hit.mean()) - (1 - ev[~ev.hi].hit.mean())
pdd = AR[AR.hi].hit.mean() - AR[~AR.hi].hit.mean()
say(f"   FAR(>=60N) - FAR(<60N) = {fd:+.3f}  95% CI {fmt(boot(ev.season, fardiff))}")
say(f"   POD(>=60N) - POD(<60N) = {pdd:+.3f}  95% CI {fmt(boot(AR.season, poddiff))}")
say("   Matching as the calibration: a fix within 400 km at the same time (800 km for pressure-less")
say("   archive events). Peak = pipeline A's peak_lat (time of the track's maximum index).")
say(f"   ERA5 share of events peaking >= 60N: {ev.hi.mean() * 100:.1f}% (2004-05 on); "
    f"all seasons {(cat[(cat.role == 'event') & (cat.basin == 'atl')].peak_lat >= 60).mean() * 100:.1f}%")

# ---------------------------------------------------------------- gust location
files = sorted(glob.glob(os.path.join(WORK, "gustloc", "*.csv")))
if not files:
    say("\n(no gustloc output yet)")
    open(os.path.join(HERE, "results.txt"), "w").write("\n".join(out) + "\n")
    sys.exit()
G = pd.concat([pd.read_csv(f) for f in files if os.path.getsize(f) > 0], ignore_index=True)
G.to_csv(os.path.join(HERE, "gustloc_fixes.csv"), index=False)
G["season"] = G.track.map(ev.set_index("track").season)
G["hi"] = G.fix_lat >= 60
G["lon180"] = np.where(G.fix_lon > 180, G.fix_lon - 360, G.fix_lon)
G["reg"] = [region(a, b) for a, b in zip(G.fix_lat, G.lon180)]
say(f"\n4. Where the gust maximum sits (pipeline A, {len(G)} Atlantic HF-strength fixes, {G.time.nunique()} times)")
say(f"   reproduction: |g800 recomputed - catalog| max {np.abs(G.g800 - G.fix_g800).max():.2f} kt, "
    f"centre match max {G.cen_dkm.max():.1f} km")
G["remote"] = G.max_dkm > 400
G["nearGL"] = G.max_dgl <= 300
G["coast100"] = G.max_dgl <= 100
G["terrain"] = G.remote & G.nearGL
for lab, x in (("fix >= 60N", G[G.hi]), ("fix < 60N", G[~G.hi])):
    say(f"   {lab:11s} n={len(x):5d}  median dist centre->max {x.max_dkm.median():5.0f} km  "
        f"max >400 km out {100 * x.remote.mean():5.1f}%  max <=300 km of Greenland {100 * x.nearGL.mean():5.1f}%  "
        f"both {100 * x.terrain.mean():5.1f}%  <=100 km of Greenland {100 * x.coast100.mean():5.1f}%  "
        f"<=50 km of any land {100 * (x.max_dland <= 50).mean():5.1f}%")
for r in ["S<60N", "Labrador/Davis", "CapeFarewell/Irminger", "DenmarkStr/Iceland", "Norwegian"]:
    x = G[G.reg == r]
    if len(x):
        say(f"   {r:22s} n={len(x):5d}  centre->max median {x.max_dkm.median():5.0f} km  terrain {100 * x.terrain.mean():5.1f}%  "
            f"near Greenland {100 * x.nearGL.mean():5.1f}%  centre MSLP median {x.msl.median():6.1f}")


def tdiff(ss):
    a = pd.concat([G[(G.season == q) & G.hi] for q in ss]); b = pd.concat([G[(G.season == q) & ~G.hi] for q in ss])
    return a.terrain.mean() - b.terrain.mean()


say(f"   terrain share (max >400 km from centre and <=300 km from Greenland), >=60N minus <60N: "
    f"{G[G.hi].terrain.mean() - G[~G.hi].terrain.mean():+.3f}  95% CI {fmt(boot(G.season, tdiff))}")
hiw = G[G.hi & G.max_ice.notna()]
say(f"   at >=60N fixes: gust max over sea ice > 0.15 in {100 * (hiw.max_ice > 0.15).mean():.1f}% (n={len(hiw)})")
tj = G[G.terrain & G.hi & G.max_wdir.notna()]
if len(tj):
    wd = tj.max_wdir
    say(f"   terrain-type maxima at >=60N (n={len(tj)}): wind from W (225-315) {100 * ((wd >= 225) & (wd < 315)).mean():.1f}%, "
        f"from N-NE (315-90) {100 * ((wd >= 315) | (wd < 90)).mean():.1f}%, from E-S (90-225) {100 * ((wd >= 90) & (wd < 225)).mean():.1f}%")
    cf = tj[(tj.max_lat < 61.5) & (tj.max_lon > 310) & (tj.max_lon < 325)]
    say(f"   of those, maximum within 59-61.5N 35-50W (Cape Farewell): {100 * len(cf) / len(tj):.1f}%")

# 5. masks on the track index
say("\n5. Track index with terrain, ice, or remote points removed (events 2004-05 on)")
mx = G.groupby("track")[["g800", "g800_gl100", "g800_gl300", "g800_land50", "g800_ice", "g800_r400"]].max()
E5 = ev.set_index("track").join(mx, how="left")
# events with no HF-strength in-domain fix should not exist; report if any
say(f"   events {len(E5)}, with gustloc rows {E5.g800.notna().sum()}")
E5 = E5[E5.g800.notna()]
for col, lab in (("g800_gl100", "ocean within 100 km of Greenland removed"),
                 ("g800_gl300", "ocean within 300 km of Greenland removed"),
                 ("g800_land50", "ocean within 50 km of any land removed"),
                 ("g800_ice", "points under sea ice > 0.15 removed (fixes >= 55N only)"),
                 ("g800_r400", "only points within 400 km of the centre")):
    v = E5[col]
    if col == "g800_ice":
        v = v.fillna(E5.g800)          # ice only fetched >= 55N; south of that no ice in this domain
    lost = v < THR
    hi = E5.peak_lat >= 60
    say(f"   {lab:58s} lose HF: all {100 * lost.mean():5.1f}% ({lost.sum()})  peak>=60N {100 * lost[hi].mean():5.1f}%  "
        f"peak<60N {100 * lost[~hi].mean():5.1f}%")
    E5["lost_" + col] = lost
for col in ("g800_gl100", "g800_ice"):
    def fr(ss, col=col):
        return np.mean(np.concatenate([E5[E5.season == q]["lost_" + col].values for q in ss]))
    say(f"   share lost under {col}: 95% CI {fmt(100 * boot(E5.season, fr), 1)} %")
newhi = (E5.peak_lat >= 60) & ~E5.lost_g800_gl300
say(f"   peak>=60N share if Greenland 300 km band removed: {100 * newhi.sum() / (~E5.lost_g800_gl300).sum():.1f}% "
    f"(now {100 * (E5.peak_lat >= 60).mean():.1f}%)")

# 6. does the terrain dependence line up with archive misses?
say("\n6. Terrain-dependent ERA5 events against the archive")
dep = E5.lost_g800_gl300
for lab, m in (("index needs the Greenland 300 km band", dep), ("does not", ~dep)):
    x = E5[m]
    say(f"   {lab:40s} n={len(x):4d}  matched to an archive event {100 * x.hit.mean():5.1f}%  FAR {1 - x.hit.mean():.3f}")
dep1 = E5.lost_g800_gl100
x = E5[dep1]
say(f"   index needs the Greenland 100 km band       n={len(x):4d}  matched {100 * x.hit.mean():5.1f}%")

open(os.path.join(HERE, "results.txt"), "w").write("\n".join(out) + "\n")
