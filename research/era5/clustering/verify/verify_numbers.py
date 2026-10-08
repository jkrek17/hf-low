import json, numpy as np, pandas as pd, statsmodels.api as sm, datetime as dt
rng = np.random.default_rng(12345)
R = '/home/claude/hf-low/'
NSIM_PHI, NSIM_GAP, NPERM = 500, 2000, 4000
def season_hours(t):
    t = int(t); y, m, d, h = t//1000000, t//10000 % 100, t//100 % 100, t % 100
    s = y if m >= 10 else y-1
    return s, (dt.datetime(y, m, d, h) - dt.datetime(s, 10, 1)).total_seconds()/3600
def umean(lat, lon):
    la, lo = np.radians(np.asarray(lat, float)), np.radians(np.asarray(lon, float) % 360)
    v = np.array([(np.cos(la)*np.cos(lo)).mean(), (np.cos(la)*np.sin(lo)).mean(), np.sin(la).mean()])
    return np.degrees(np.arctan2(v[2], np.hypot(v[0], v[1]))), np.degrees(np.arctan2(v[1], v[0])) % 360
# ---- events
J = json.load(open(R+'docs/data/hf-lows.json')); lf = J['lowFields']; ff = J['fixFields']
rows = []
for l in J['lows']:
    r = dict(zip(lf, l))
    if r['cls'] != 'low': continue
    hf = [dict(zip(ff, f)) for f in r['fixes'] if f[ff.index('cat')] == 'HF']
    if not hf: continue
    hf.sort(key=lambda f: f['date'])
    s, h = season_hours(hf[0]['date']); la, lo = umean([f['lat'] for f in hf], [f['lon'] for f in hf])
    rows.append((r['basin'], s, h, la, lo))
ARCH = pd.DataFrame(rows, columns=['basin', 'season', 'hours', 'lat', 'lon'])
T = pd.read_csv(R+'research/era5/hf_history/results/all_tracks.csv.gz')
T = T[T.basin.isin(['atl', 'pac']) & (T.gust800_kt >= 71.7)]
C = pd.read_csv(R+'research/era5/hf_history/results/era5_hf_catalog_tracks.csv')
C = C[C.g800 >= 71.7]
pos = {k: umean(g.lat, g.lon) for k, g in C.groupby('track')}
rows = []
for r in T.itertuples():
    s, h = season_hours(r.peak_time)
    la, lo = pos.get(r.track, (np.nan, np.nan))
    rows.append((r.basin, s, h, la, lo))
PIPE = pd.DataFrame(rows, columns=['basin', 'season', 'hours', 'lat', 'lon'])
print('pipeA events with no position:', int(PIPE.lat.isna().sum()))
W = 5040
def sub(df, b, s0, s1, inwin=True):
    d = df[(df.basin == b) & df.season.between(s0, s1)]
    return d[(d.hours >= 0) & (d.hours < W)] if inwin else d[~((d.hours >= 0) & (d.hours < W))]
# ---- A1
print('\n=== A1 ===')
for nm, df, s0, s1 in [('ARCH 2004-25', ARCH, 2004, 2025), ('PIPE 2004-25', PIPE, 2004, 2025), ('PIPE 1979-2003', PIPE, 1979, 2003)]:
    for b in ['atl', 'pac']:
        print(nm, b, 'in', len(sub(df, b, s0, s1)), 'out', len(sub(df, b, s0, s1, False)))
# ---- N1
mon = np.array([(dt.date(2001, 10, 1)+dt.timedelta(d)).month for d in range(210)])
mons = sorted(set(mon))
def counts(d, s0, s1):
    ss = list(range(s0, s1+1)); Y = np.zeros((len(ss), 210))
    for s, h in zip(d.season, d.hours): Y[s-s0, int(h//24)] += 1
    return ss, Y
def design(ss, extra=None):
    n = len(ss); cols = [np.tile((mon == m).astype(float), n) for m in mons]
    tr = np.repeat(np.array(ss, float) - np.mean(ss), 210); X = np.column_stack(cols+[tr])
    if extra is not None: X = np.column_stack([X, extra.ravel()])
    return X
def fit(y, X): return sm.GLM(y, X, family=sm.families.Poisson()).fit()
def phi(Y, M, L, k=8):
    nb = 210//L; yb = Y[:, :nb*L].reshape(len(Y), nb, L).sum(2); mb = M[:, :nb*L].reshape(len(Y), nb, L).sum(2)
    return ((yb-mb)**2/mb).sum()/(yb.size-k)
data = {}
for nm, df, s0, s1 in [('ARCH', ARCH, 2004, 2025), ('PIPE', PIPE, 2004, 2025)]:
    for b in ['atl', 'pac']:
        d = sub(df, b, s0, s1); ss, Y = counts(d, s0, s1); X = design(ss)
        M = fit(Y.ravel(), X).fittedvalues.reshape(Y.shape); data[(nm, b)] = (d, ss, Y, X, M)
print('\n=== A2 (archive) ===')
rep = {('atl', 7): -0.056, ('atl', 30): 0.113, ('pac', 7): -0.028, ('pac', 30): 0.077}
for b in ['atl', 'pac']:
    d, ss, Y, X, M = data[('ARCH', b)]
    sims = {7: [], 30: []}
    for _ in range(NSIM_PHI):
        Ys = rng.poisson(M); Ms = fit(Ys.ravel(), X).fittedvalues.reshape(Y.shape)
        for L in (7, 30): sims[L].append(phi(Ys, Ms, L))
    for L in (7, 30):
        p = phi(Y, M, L); m = np.mean(sims[L])
        print(f'{b} L={L}: phi={p:.4f} simmean={m:.4f} E={p-m:.4f} (reported {rep[(b,L)]})')
# ---- A3
def gaps(d, ss):
    n = 0
    for s in ss:
        h = np.sort(d[d.season == s].hours.values); n += int((np.diff(h) <= 48).sum())
    return n
def gap_exp(M, ss, nsim):
    tot = []
    for _ in range(nsim):
        n = 0
        for i in range(len(ss)):
            c = rng.poisson(M[i]); days = np.repeat(np.arange(210), c)
            h = np.sort(days*24 + 6*rng.integers(0, 4, len(days))); n += int((np.diff(h) <= 48).sum())
        tot.append(n)
    return np.mean(tot)
print('\n=== A3 ===')
for nm in ['ARCH', 'PIPE']:
    for b in ['atl', 'pac']:
        d, ss, Y, X, M = data[(nm, b)]
        print(nm, b, 'observed', gaps(d, ss), 'expected', round(gap_exp(M, ss, NSIM_GAP), 1))
# ---- A4 Knox
def hav(la1, lo1, la2, lo2):
    a = np.radians; p1, p2 = a(la1), a(la2); dl = a(lo2)-a(lo1)
    h = np.sin((p2-p1)/2)**2+np.cos(p1)*np.cos(p2)*np.sin(dl/2)**2
    return 2*6371*np.arcsin(np.sqrt(np.minimum(h, 1)))
def knox(d, nperm, bins=None, dmax=72):
    d = d.dropna(subset=['lat']).reset_index(drop=True)
    h = d.hours.values.copy(); s = d.season.values; la = d.lat.values; lo = d.lon.values
    mo = mon[(h//24).astype(int)]
    iu = [np.triu_indices(len(d), 1)][0]; I, Jx = iu
    ok = s[I] == s[Jx]; I, Jx = I[ok], Jx[ok]
    near = hav(la[I], lo[I], la[Jx], lo[Jx]) <= 1000; I, Jx = I[near], Jx[near]
    def stat(hh):
        dtt = np.abs(hh[I]-hh[Jx])
        if bins is None: return int((dtt <= dmax).sum())
        return np.array([((dtt > a) & (dtt <= b)).sum() for a, b in bins])
    obs = stat(h); groups = [np.where((s == a) & (mo == m))[0] for a in np.unique(s) for m in mons]
    sims = []
    for _ in range(nperm):
        hh = h.copy()
        for g in groups:
            if len(g) > 1: hh[g] = h[rng.permutation(g)]
        sims.append(stat(hh))
    sims = np.array(sims, float); return obs, sims
print('\n=== A4 ===')
for nm, df, s0, s1 in [('ARCH 2004-25', ARCH, 2004, 2025), ('PIPE 2004-25', PIPE, 2004, 2025), ('PIPE 1979-2003', PIPE, 1979, 2003)]:
    for b in ['atl', 'pac']:
        d = sub(df, b, s0, s1); o, sm_ = knox(d, NPERM)
        print(nm, b, 'obs', o, 'exp', round(sm_.mean(), 1), 'ratio', round(o/sm_.mean(), 2), 'p1', (1+(sm_ >= o).sum())/(1+len(sm_)))
# ---- A5
def cpc(f):
    out = {}
    for line in open(f):
        y, m, d = int(line[0:4]), int(line[4:8]), int(line[8:12]) if False else None, None
    return out
def readcpc(f):
    s = {}
    for line in open(f):
        t = line.split()
        if len(t) < 4: continue
        y, m, d, v = int(t[0]), int(t[1]), int(t[2]), float(t[3])
        s[dt.date(y, m, d)] = np.nan if v <= -99 else v
    return pd.Series(s)
def readcpc_fw(f):
    s = {}
    for line in open(f):
        if len(line.strip()) < 10: continue
        y, m, d = int(line[0:4]), int(line[4:7]), int(line[7:10])
        v = float(line[10:]); s[dt.date(y, m, d)] = np.nan if v <= -99 else v
    return pd.Series(s)
idx = {}
for k, n in [('atl', 'nao'), ('pac', 'pna')]:
    idx[k] = readcpc_fw(f'/mnt/project-files/teleconnection-test/cpc_indices/norm.daily.{n}.index.b500101.current.ascii')
print('\n=== A5 ===  (CPC rows', {k: len(v) for k, v in idx.items()}, 'NaN', {k: int(v.isna().sum()) for k, v in idx.items()}, ')')
for nm in ['ARCH', 'PIPE']:
    for b in ['atl', 'pac']:
        d, ss, Y, X, M = data[(nm, b)]; s = idx[b]; lag = np.zeros(Y.shape)
        for i, sea in enumerate(ss):
            for dd in range(210):
                day = dt.date(sea, 10, 1)+dt.timedelta(dd)
                lag[i, dd] = np.nanmean([s.get(day-dt.timedelta(k), np.nan) for k in range(4, 11)])
        z = (lag-lag.mean())/lag.std(ddof=1)
        for ddof in (1,):
            r = fit(Y.ravel(), np.column_stack([X, z.ravel()]))
            print(nm, b, 'exp(coef)=%.4f' % np.exp(r.params[-1]), 'p=%.4f' % r.pvalues[-1], '(sd ddof=1)', ' with ddof0: %.4f' % np.exp(fit(Y.ravel(), np.column_stack([X, ((lag-lag.mean())/lag.std()).ravel()])).params[-1]))
# ---- A6
print('\n=== A6 ===')
bins = [(0, 24), (24, 48), (48, 72), (72, 96), (96, 168)]
d = sub(ARCH, 'atl', 2004, 2025); o, sm_ = knox(d, NPERM, bins=bins)
print('obs', o, 'exp', sm_.mean(0).round(1), 'ratio', (o/sm_.mean(0)).round(2))
