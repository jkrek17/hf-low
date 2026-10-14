#!/usr/bin/env python3
"""Write tools/gfs_gust_rule.html: the "GFS gust rule" section of docs/bomb-environment.html.

    python3 -I tools/gfs_gust_rule.py --fixes raw_gust_fixes.csv.gz --raw-json raw_gust.json --tc tc_candidates.csv

Inputs: the compact product of research/era5/gfs_raw_gust (hf-low research PR 182; shared folder
/mnt/project-files/gfs-raw-gust/) and research/era5/tc_candidates.csv (IBTrACS points within 400 km of a pipeline A
candidate). One row per GFS track fix (pipeline A detector on GFS forecasts, 5 seasons 2021-22 to 2025-26, October to April,
00/12 UTC cycles, leads 0-48 h). G = highest ocean gust within 800 km (kt). Y0 = OPC archive hurricane-force fix within +-3 h
and 600 km of the valid time.

Threshold rules, reproduced from the stored scores in raw_gust.json: "matched" = G strictly above the leave-one-season-out
count-matched threshold T_heldout_Y0[basin][season]; "64 kt" = G strictly above 64.
Tropical flag (POST HOC, added after the band table was seen): the fix is within 400 km of a tc_candidates point at the same
valid time and basin. It is a first-look proximity match, not a tightened tropical-cyclone definition.
"""
import argparse, csv, gzip, json, math, os, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BANDS = [(0, 60), (60, 64), (64, 68), (68, 72), (72, 76), (76, 80), (80, 85), (85, 90), (90, 999)]


def km(a, b, c, d):
    a, b, c, d = map(math.radians, (a, b, c, d))
    return 6371 * math.acos(min(1, math.sin(a) * math.sin(c) + math.cos(a) * math.cos(c) * math.cos(b - d)))


def label(lo, hi):
    return "under 60" if lo == 0 else ("90 and over" if hi == 999 else f"{lo} to {hi}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixes", required=True)
    ap.add_argument("--raw-json", required=True)
    ap.add_argument("--tc", required=True)
    ap.add_argument("--out", default=os.path.join(ROOT, "tools", "gfs_gust_rule.html"))
    a = ap.parse_args()
    T = json.load(open(a.raw_json))["T_heldout_Y0"]
    tc = collections.defaultdict(list)
    for r in csv.DictReader(open(a.tc)):
        t = r["time_utc"]
        tc[(t[:4] + t[5:7] + t[8:10] + t[11:13], r["basin"])].append((float(r["lat"]), float(r["lon"]) % 360))
    rows = list(csv.DictReader(gzip.open(a.fixes, "rt")))
    for r in rows:
        r["G"] = float(r["G"]); r["Y0"] = int(r["Y0"]); r["lead"] = int(r["lead"])
        la, lo = float(r["lat"]), float(r["lon"]) % 360
        r["tc"] = any(km(la, lo, p, q) <= 400 for p, q in tc.get((r["valid_time"], r["basin"]), []))

    def cell(sel):
        n, y = len(sel), sum(r["Y0"] for r in sel)
        return f"{100 * y / n:.0f}% <span class=\"muted\">({y}/{n})</span>" if n else "&ndash;"

    band_rows = []
    for lo, hi in BANDS:
        tds = []
        for lead in (0, 24):
            sel = [r for r in rows if r["lead"] == lead and lo <= r["G"] < hi]
            tds += [cell(sel), cell([r for r in sel if not r["tc"]])]
        band_rows.append(f"<tr><td>{label(lo, hi)}</td>" + "".join(f"<td class=\"num\">{t}</td>" for t in tds) + "</tr>")

    def score(sel, hit):
        H = sum(1 for r in sel if hit(r) and r["Y0"]); FA = sum(1 for r in sel if hit(r) and not r["Y0"])
        M = sum(1 for r in sel if not hit(r) and r["Y0"]); CN = len(sel) - H - FA - M
        N = H + FA + M + CN
        e = ((H + FA) * (H + M) + (M + CN) * (FA + CN)) / N
        return H / (H + M), H / (H + FA), ((H + CN) - e) / (N - e)

    rate_rows = []
    for basin, name in (("atl", "Atlantic"), ("pac", "Pacific")):
        for lead in (0, 24, 48):
            sel = [r for r in rows if r["basin"] == basin and r["lead"] == lead]
            m = score(sel, lambda r: r["G"] > T[basin][r["season"]])
            f = score(sel, lambda r: r["G"] > 64)
            rate_rows.append(f"<tr><td>{name}</td><td class=\"num\">{lead} h</td>"
                             + "".join(f"<td class=\"num\">{100 * v:.0f}%</td>" for v in (m[0], m[1])) + f"<td class=\"num\">{m[2]:.2f}</td>"
                             + "".join(f"<td class=\"num\">{100 * v:.0f}%</td>" for v in (f[0], f[1])) + f"<td class=\"num\">{f[2]:.2f}</td></tr>")
    thr = sorted({round(v, 1) for b in T.values() for v in b.values()})
    out = f"""<section id="gfs-rule">
  <div class="tag">GFS gust rule</div>
  <h2>When the GFS forecast gust reaches about 72 kt, how often is OPC calling hurricane force?</h2>
  <p>This is the raw GFS forecast gust, with no model on top: for each forecast low, the highest ocean gust within 800 km, set against whether OPC&rsquo;s archive had a hurricane-force fix near the same place and time. Five seasons (2021-22 to 2025-26), October to April, forecasts from 00 and 12 UTC, 0 to 48 hours ahead; {len(rows):,} forecast fixes. The cut-off that matches OPC&rsquo;s hurricane-force counts is about 72 kt ({thr[0]} to {thr[-1]} kt depending on the season held out).</p>
  <p><strong>Read it as odds, not a switch.</strong> A bigger gust does not make OPC&rsquo;s call certain: the share peaks near two in three and does not reach 100%. No gust band is a sure thing.</p>
  <div class="tablewrap"><table id="gfsBands">
    <caption class="muted">Share of GFS forecast lows, by their peak gust, where OPC had a hurricane-force fix (hits/lows in brackets)</caption>
    <thead><tr><th>GFS peak gust (kt)</th><th class="num">Forecast hour 0, all lows</th><th class="num">Hour 0, tropical-linked removed</th><th class="num">Hour 24, all lows</th><th class="num">Hour 24, tropical-linked removed</th></tr></thead>
    <tbody>{"".join(band_rows)}</tbody>
  </table></div>
  <p class="muted">&ldquo;Tropical-linked removed&rdquo; drops forecast fixes within 400 km of a hurricane or tropical-storm track point (IBTrACS) at the same time. That cut was added after seeing the first table, so treat it as a post hoc check; it changes little. The matching is a first look, not a tightened tropical definition. The two highest bands hold few lows, so their shares are rough.</p>
  <div class="tablewrap"><table id="gfsRates">
    <caption class="muted">Two cut-offs compared: share of OPC hurricane-force lows the forecast caught, share of its calls where OPC also had hurricane force, and a skill score (Heidke, 0 = no skill, 1 = perfect)</caption>
    <thead><tr><th>Basin</th><th class="num">Forecast hour</th><th class="num">72 kt: caught</th><th class="num">72 kt: OPC agreed</th><th class="num">72 kt: skill</th><th class="num">64 kt: caught</th><th class="num">64 kt: OPC agreed</th><th class="num">64 kt: skill</th></tr></thead>
    <tbody>{"".join(rate_rows)}</tbody>
  </table></div>
  <p class="limit"><strong>OPC is not a measurement.</strong> Its call is an analyst&rsquo;s judgement of all the data available, and can lean towards caution for customers, so it can sit above or below the wind the model shows. A GFS gust with no OPC call is not necessarily a GFS error, and a missed OPC call may be a GFS under-forecast or an OPC over-call. This table cannot tell which: that needs observed winds (scatterometer, ships, buoys). &ldquo;OPC agreed&rdquo; means OPC also called it, not that the forecast was correct.</p>
  <p class="limit">At about 72 kt the forecast catches roughly half of OPC&rsquo;s hurricane-force lows and OPC agrees with about half of its calls. Lowering the cut to 64 kt catches about nine in ten but OPC agrees with only about three in ten of the calls. ERA5 proxy lows are not used here: these are GFS forecast fixes found by the same tracker, scored against OPC&rsquo;s archive. Five seasons; fixes within one storm are not independent, so no intervals are shown. The &ldquo;72 kt&rdquo; cut-off was matched to OPC counts on these same seasons, one season held out at a time. Why the share stops rising above about 85 kt is not known.</p>
</section>
"""
    open(a.out, "w").write(out)
    print("wrote", a.out, len(rows), "fixes")


if __name__ == "__main__":
    main()
