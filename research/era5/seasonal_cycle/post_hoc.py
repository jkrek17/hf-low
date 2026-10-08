"""Post hoc checks, run after the pre-registered family had been seen. Every row here is labelled post hoc in the README
and is not in the 56-test family or its FDR.   python3 post_hoc.py <repo root> <work dir> <out dir>"""
import sys, os, runpy
import numpy as np, pandas as pd

ROOT, WORK, OUT = sys.argv[1:4]
sys.argv = ["seasonal.py", ROOT, WORK, OUT, "2000"]
g = runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "seasonal.py"), run_name="seasonal_import")
Q, xm, contrasts, boot_w, CIDX, GUST_SEASONS, SEAS_V250 = (g[k] for k in ("Q", "xm", "contrasts", "boot_w", "CIDX", "GUST_SEASONS", "SEAS_V250"))
rows = []


def ci(vb):
    vb = vb[np.isfinite(vb)]
    p = min(1.0, 2 * (min((vb <= 0).sum(), (vb >= 0).sum()) + 1) / (len(vb) + 1))
    return np.percentile(vb, 2.5), np.percentile(vb, 97.5), p


for b in ("pac", "atl"):
    sel = np.array([s in SEAS_V250 for s in GUST_SEASONS])
    n1, d1 = Q[("Q1", b)]["num"][sel], Q[("Q1", b)]["den"][sel]
    n9, d9 = Q[("Q9", b)]["num"], Q[("Q9", b)]["den"]
    S = sel.sum()
    wb = boot_w(S)
    for k, c in ((0, "D_Jan"), (1, "M")):
        v0 = contrasts(xm(n1, d1, np.ones(S), True))[k] - contrasts(xm(n9, d9, np.ones(S), True))[k]
        vb = contrasts(xm(n1, d1, wb, True))[k] - contrasts(xm(n9, d9, wb, True))[k]
        lo, hi, p = ci(vb)
        rows.append(dict(check="PH1 archive HF count minus v250 eddy (paired, 17 seasons)", basin=b, contrast=c, estimate=v0, lo=lo, hi=hi, p=p))
    # PH2: fraction of the pipeline A HF midwinter excess carried by the share, M(share)/M(HF)
    n2, d2 = Q[("Q2", b)]["num"], Q[("Q2", b)]["den"]
    n4, d4 = Q[("Q4", b)]["num"], Q[("Q4", b)]["den"]
    S = len(GUST_SEASONS)
    wb = boot_w(S)
    for k, c in ((0, "D_Jan"), (1, "M")):
        v0 = contrasts(xm(n4, d4, np.ones(S), True))[k] / contrasts(xm(n2, d2, np.ones(S), True))[k]
        vb = contrasts(xm(n4, d4, wb, True))[k] / contrasts(xm(n2, d2, wb, True))[k]
        rows.append(dict(check="PH2 fraction of pipeline A HF contrast carried by the HF share, f", basin=b, contrast=c, estimate=v0,
                         lo=np.nanpercentile(vb, 2.5), hi=np.nanpercentile(vb, 97.5), p=np.nan))
    # PH4: Jan vs Dec, Feb vs Jan for the HF counts
    for qid in ("Q1", "Q2"):
        n_, d_ = Q[(qid, b)]["num"], Q[(qid, b)]["den"]
        for nm, (m1, m2) in (("Jan-Dec", (1, 12)), ("Feb-Jan", (2, 1))):
            f = lambda w: (lambda x: x[..., CIDX[m1]] - x[..., CIDX[m2]])(xm(n_, d_, w, True))
            lo, hi, p = ci(f(wb))
            rows.append(dict(check=f"PH4 {qid} ln rate difference {nm}", basin=b, contrast=nm, estimate=f(np.ones(S)), lo=lo, hi=hi, p=p))
# PH3: archive first-HF latitude, Pacific, per season Jan minus mean(Dec, Feb): sign count
a = g["arch"]["pac"]
a = a[a.s.isin(GUST_SEASONS)]
per = []
for s, d in a.groupby("s"):
    j, dc, f = (d[d.m == k].Latitude.mean() for k in (1, 12, 2))
    per.append(j - (dc + f) / 2)
per = np.array(per)
rows.append(dict(check="PH3 archive first-HF latitude, Pacific: seasons with Jan below mean(Dec,Feb) / seasons with all three months",
                 basin="pac", contrast="D_Jan by season", estimate=np.nanmedian(per), lo=int((per < 0).sum()), hi=int(np.isfinite(per).sum()), p=np.nan))
# PH5: which month holds the maximum (bootstrap share of draws), after Jason's reframing to the peak month
MN = [6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5]
for b in ("pac", "atl"):
    for qid in ("Q1", "Q2", "Q3", "Q4", "Q7"):
        q = Q[(qid, b)]
        S = len(q["seasons"])
        wb = boot_w(S)
        x0 = xm(q["num"], q["den"], np.ones(S), True)
        xb = xm(q["num"], q["den"], wb, True)
        am = np.nanargmax(xb, axis=1)
        for mon in (12, 1, 2):
            rows.append(dict(check=f"PH5 {qid}: share of bootstrap draws in which the maximum month is this one", basin=b, contrast=f"month {mon}",
                             estimate=(am == CIDX[mon]).mean(), lo=MN[int(np.nanargmax(x0))], hi=np.nan, p=np.nan))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "post_hoc.csv"), index=False, float_format="%.4g")
print(pd.DataFrame(rows).round(3).to_string())
