"""POST HOC: GH, GRAD and NAO on all 4,559 candidate times (no control sampling).

GH, GRAD and NAO are known at every candidate time (stage 1), so the matched and the adjusted estimates can be
checked without the 3:1 control sample the plan used for the SST-based ingredient. Not pre-registered; the
pre-registered results stay as the primary ones. ERA5 proxy. Needs analyse.py to have run. Writes results/posthoc2.txt.
"""
import os, sys
import numpy as np, pandas as pd
from statsmodels.discrete.conditional_models import ConditionalLogit
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "greenland_jets")
RES = os.path.join(HERE, "results")
sys.path.insert(0, HERE)
import analyse as an
out = []


def say(s):
    print(s); out.append(s)


d = pd.read_csv(f"{RES}/stage1_times.csv")
d["w"] = 1.0
d["STAB"] = np.nan
d["GBI"] = np.nan
nao = pd.read_csv(f"{an.CPC}/norm.daily.nao.index.b500101.current.ascii", sep=r"\s+|(?<=\d)(?=-99)", header=None,
                  names=["y", "m", "d", "NAO"], engine="python")
nao["date"] = pd.to_datetime(dict(year=nao.y, month=nao.m, day=nao.d))
t = pd.to_datetime(d.time.astype(str), format="%Y%m%d%H")
d["month"] = t.dt.month
d["date"] = t.dt.normalize()
d = d.merge(nao[["date", "NAO"]], on="date", how="left")
d["x"] = (d.ref_lon - an.CF[1]) * np.cos(np.radians(d.ref_lat)) * an.KM
d["y"] = (d.ref_lat - an.CF[0]) * an.KM
d["Dcf"] = np.hypot(d.x, d.y)
h = np.radians(d.ref_heading)
d["u"], d["v"] = d.ref_speed * np.sin(h), d.ref_speed * np.cos(h)
d["mmiss"] = d.ref_speed.isna().astype(float)
st = pd.read_pickle(f"{WORK}/zstats.pkl").to_dict()
z = an.prep(d, {k: tuple(v) for k, v in st.items()}).reset_index(drop=True)
say(f"POST HOC on all candidate times: {len(z)} times, {int(z.case.sum())} cases (ERA5 proxy, outcome Y_T)")
cols = an.BASE + ["GH", "GRAD", "NAO", "u", "v"]
y = z.case.values.astype(float)
X = an.mat(z, cols)
b = an.logit(X, y)
sid = an.seasons_idx(z)
rng = np.random.default_rng(3)
bb = np.array([an.logit(X[ix], y[ix]) for ix in (an.boot_idx(sid, rng) for _ in range(500))])
say("adjusted logistic on all times (pre-registered baseline), OR per SD [season-block 95% CI]:")
for c in ("GH", "GRAD", "NAO"):
    j = cols.index(c) + 1
    say(f"  {c:5s} {np.exp(b[j]):.2f} [{np.exp(np.percentile(bb[:, j], 2.5)):.2f}, {np.exp(np.percentile(bb[:, j], 97.5)):.2f}]")

# matching on all non-case times
def match(z, rng, k, cal):
    sets, used, dist = [], set(), []
    cs = z[z.case == 1].sample(frac=1.0, random_state=int(rng.integers(1 << 31)))
    ctrl = z[z.case == 0]
    sdv = z[["ref_msl", "ref_lat", "ref_lon"]].std().values
    for ci, c in cs.iterrows():
        pool = ctrl[(ctrl.month == c.month) & (~ctrl.index.isin(used))]
        dd = np.sqrt((((pool[["ref_msl", "ref_lat", "ref_lon"]].values - c[["ref_msl", "ref_lat", "ref_lon"]].values.astype(float)) / sdv) ** 2).sum(1))
        o = np.argsort(dd)[:k]
        o = o[dd[o] <= cal]
        if len(o) == 0:
            continue
        pick = pool.index[o]
        used.update(pick)
        sets.append((ci, list(pick)))
        dist += list(dd[o])
    return sets, np.array(dist)


def fit(z, sets, cc):
    rows, g, yy = [], [], []
    for n, (ci, ps) in enumerate(sets):
        for j, i in enumerate([ci] + ps):
            rows.append(z.loc[i, cc].values.astype(float)); g.append(n); yy.append(1 if j == 0 else 0)
    return ConditionalLogit(np.array(yy), np.array(rows), groups=np.array(g)).fit(disp=0).params


cc = ["GH", "GRAD", "NAO", "u", "v"]
for k, cal in ((5, 99), (5, 0.5), (2, 0.25)):
    rng = np.random.default_rng(an.SEED + 70)
    sets, dist = match(z, rng, k, cal)
    ms = np.mean([len(p) for _, p in sets])
    dm = np.mean([abs(z.loc[ci, "ref_msl"] - z.loc[p, "ref_msl"].mean()) for ci, p in sets]) * 1  # z units
    b0 = fit(z, sets, cc)
    sidx = {s: [n for n, (ci, _) in enumerate(sets) if z.loc[ci, "season"] == s] for s in sorted(z.season.unique())}
    bs = []
    for _ in range(300):
        ss = rng.choice(list(sidx), len(sidx))
        try:
            bs.append(fit(z, [sets[n] for s in ss for n in sidx[s]], cc))
        except Exception:
            pass
    bs = np.array(bs)
    msd = np.sqrt(((z.loc[[ci for ci, _ in sets], ["ref_lat", "ref_lon"]].values - np.array([z.loc[p, ["ref_lat", "ref_lon"]].values.mean(0) for _, p in sets])) ** 2).mean(0))
    say(f"\nmatched on MSLP, lat, lon within month: up to {k} nearest controls, caliper {cal} SD-distance: "
        f"{len(sets)} sets, {ms:.2f} controls per set, mean match distance {dist.mean():.2f}; "
        f"case-minus-control reference-low MSLP {np.mean([z.loc[ci,'ref_msl'] - z.loc[p,'ref_msl'].mean() for ci,p in sets]) * st['ref_msl'][1]:+.1f} hPa, "
        f"lat {np.mean([z.loc[ci,'ref_lat'] - z.loc[p,'ref_lat'].mean() for ci,p in sets]):+.2f}, lon {np.mean([z.loc[ci,'ref_lon'] - z.loc[p,'ref_lon'].mean() for ci,p in sets]):+.2f}")
    for i, c in enumerate(cc[:3]):
        say(f"  {c:5s} OR per SD {np.exp(b0[i]):.2f} [{np.exp(np.percentile(bs[:, i], 2.5)):.2f}, {np.exp(np.percentile(bs[:, i], 97.5)):.2f}]")
open(f"{RES}/posthoc2.txt", "w").write("\n".join(out) + "\n")
