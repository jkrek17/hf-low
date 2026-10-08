"""Winter ENSO flavor indices from the ERA5 SST box means (pre-registered definitions).

Anomaly = monthly box mean minus its 1991-2020 calendar-month mean. DJF of winter s
= mean of Dec(s), Jan(s+1), Feb(s+1). Also DJF PNA from the CPC daily index and a
check against the official ONI in docs/data/teleconnections.json.

usage: indices.py SST_RAW_CSV CPC_DIR REPO_ROOT OUT_CSV
"""
import sys, os, json, re, glob
import numpy as np, pandas as pd

def main(raw, cpc, repo, out):
    d = pd.read_csv(raw)
    for k in ["nino3", "nino4", "nino34", "emi_a", "emi_b", "emi_c"]:
        clim = d[(d.year >= 1991) & (d.year <= 2020)].groupby("month")[k].mean()
        d[k + "_a"] = d[k] - d.month.map(clim)
    d["emi"] = d.emi_a_a - 0.5 * d.emi_b_a - 0.5 * d.emi_c_a
    d["n4mn3"] = d.nino4_a - d.nino3_a
    rows = []
    for s in range(1979, 2026):
        sel = d[((d.year == s) & (d.month == 12)) | ((d.year == s + 1) & (d.month.isin([1, 2])))]
        assert len(sel) == 3
        rows.append(dict(winter=s, emi=sel.emi.mean(), n34=sel.nino34_a.mean(), n3=sel.nino3_a.mean(),
                         n4=sel.nino4_a.mean(), n4mn3=sel.n4mn3.mean()))
    W = pd.DataFrame(rows)
    # PNA DJF
    rec = []
    for line in open(os.path.join(cpc, "norm.daily.pna.index.b500101.current.ascii")):
        m = re.match(r"\s*(\d{4})\s+(\d{1,2})\s+(\d{1,2})\s*(-?\d+\.\d+)\s*$", line)
        if m and float(m.group(4)) > -90:      # CPC glues a -99.000 missing value onto the day
            rec.append((int(m.group(1)), int(m.group(2)), int(m.group(3)), float(m.group(4))))
    p = pd.DataFrame(rec, columns=["y", "m", "d", "v"])
    p["date"] = pd.to_datetime(dict(year=p.y, month=p.m, day=p.d))
    W["pna"] = [p[(p.date >= f"{s}-12-01") & (p.date < f"{s+1}-03-01")].v.mean() for s in W.winter]
    # official ONI DJF
    t = json.load(open(os.path.join(repo, "docs/data/teleconnections.json")))["oni"]
    y0, m0 = int(t["start"][:4]), int(t["start"][5:7])
    oni = {}
    for i, v in enumerate(t["values"]):
        oni[(y0 + (m0 - 1 + i) // 12, (m0 - 1 + i) % 12 + 1)] = v
    SEAS = ["DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND", "NDJ"]
    for fn in glob.glob(os.path.join(cpc, "oni*.txt")):
        for line in open(fn):
            q = line.split()
            if len(q) == 4 and q[0] in SEAS:
                oni.setdefault((int(q[1]), SEAS.index(q[0]) + 1), float(q[3]))
    # ONI is a 3-month running mean labelled by the centre month; DJF = centred on Jan
    W["oni_official"] = [oni.get((s + 1, 1), np.nan) for s in W.winter]
    W.to_csv(out, index=False)
    print(W.round(2).to_string())
    print(W[["emi", "n34", "n3", "n4", "n4mn3", "pna", "oni_official"]].corr().round(2).to_string())

if __name__ == "__main__":
    main(*sys.argv[1:5])
