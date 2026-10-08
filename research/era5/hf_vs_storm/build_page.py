"""Build the Artifact HTML page from results/ (figures embedded as data URIs). usage: build_page.py REPO_ROOT OUT_HTML VERIF_TEXT_FILE"""
import sys, os, base64, html, numpy as np, pandas as pd
root, out, vf = sys.argv[1:4]
R = os.path.join(root, "research/era5/hf_vs_storm/results")
def img(name, alt):
    b = base64.b64encode(open(os.path.join(R, "figs", name), "rb").read()).decode()
    return f'<img alt="{html.escape(alt)}" src="data:image/png;base64,{b}">'
sc = pd.read_csv(os.path.join(R, "storm/scalars.csv")); s2 = pd.read_csv(os.path.join(R, "storm/scalars_secondary.csv"))
bx = pd.read_csv(os.path.join(R, "large/boxtests.csv"))
BAS = {"atl": "North Atlantic", "pac": "North Pacific"}
def f(x, d=1): return f"{x:.{d}f}"
def qf(q): return "&lt;0.01" if q < 0.01 else f"{q:.2f}"
def scal_table(df):
    rows = []
    for b in ("atl", "pac"):
        for r in df[(df.basin == b)].itertuples():
            d = 2 if abs(r.hf) < 10 else 1
            rows.append(f"<tr><td>{BAS[b]}</td><td>{html.escape(r.label)}</td><td class=n>{f(r.hf,d)}</td><td class=n>{f(r.sf,d)}</td>"
                        f"<td class=n>{f(r.diff,d)} <span class=ci>[{f(r.lo,d)}, {f(r.hi,d)}]</span></td><td class=n>{r.d_std:+.2f}</td><td class=n>{r.mde_std:.2f}</td>"
                        f"<td class=n>{qf(r.q)}</td></tr>")
    return "<table><thead><tr><th>Basin</th><th>Measure</th><th>HF</th><th>Storm-force</th><th>HF minus SF [95%]</th><th>Effect (SD)</th><th>Detectable (SD)</th><th>q</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>"
c1 = sc[sc.cmp == "C1"]
s1 = s2[(s2.family == "S1_depth_matched") & (s2.cmp == "C2")]
def box_table():
    nm = {"z500": "Z500 (m)", "mslp": "MSLP (hPa)", "u250": "250 hPa zonal wind (m/s)", "sst": "SST (K)", "tcwv": "Column water vapour (kg/m2)"}
    rows = []
    p = bx[(bx.cmp == "C1") & bx.cfg.isin(["L", "k0"])]
    for b in ("atl", "pac"):
        for v in nm:
            cells = []
            for cfg in ("L", "k0"):
                r = p[(p.basin == b) & (p["var"] == v) & (p.cfg == cfg)].iloc[0]
                d = 3 if v == "sst" else 2 if v in ("mslp", "tcwv", "u250") else 1
                cells.append(f"<td class=n>{f(r['diff'],d)} <span class=ci>[{f(r.lo,d)}, {f(r.hi,d)}]</span></td><td class=n>{qf(r.q_primary)}</td><td class=n>{f(r.mde80,d)}</td>")
            rows.append(f"<tr><td>{BAS[b]}</td><td>{nm[v]}</td>{''.join(cells)}</tr>")
    return "<table><thead><tr><th>Basin</th><th>Box mean anomaly, HF onset minus SF peak</th><th>Days -10 to -4 [95%]</th><th>q</th><th>MDE</th><th>Onset day [95%]</th><th>q</th><th>MDE</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>"
verif = open(vf).read()
page = f'''<title>HF vs Storm-Force Composites</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,500;6..72,600&family=Public+Sans:wght@400;500;600&family=JetBrains+Mono:wght@400;500&display=swap">
<style>
/* Layout: single reading column; figures bleed to the wide column; HF in warm red, storm-force in slate blue throughout. */
:root {{ --bg:#f5f6f7; --panel:#ffffff; --ink:#1b2430; --muted:#586374; --rule:#d6dbe1; --hf:#b23a25; --sf:#33597f; --note:#fff4e5; --noteink:#5a3a07;
  --serif:'Newsreader',Georgia,serif; --sans:'Public Sans',system-ui,sans-serif; --mono:'JetBrains Mono',ui-monospace,monospace; }}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{ --bg:#12171d; --panel:#1b222b; --ink:#e6eaef; --muted:#9ba6b5; --rule:#2c3541; --hf:#ee8a73; --sf:#8db4dc; --note:#2b2316; --noteink:#e9cf9d; color-scheme:dark }} }}
:root[data-theme="dark"] {{ --bg:#12171d; --panel:#1b222b; --ink:#e6eaef; --muted:#9ba6b5; --rule:#2c3541; --hf:#ee8a73; --sf:#8db4dc; --note:#2b2316; --noteink:#e9cf9d; color-scheme:dark }}
body {{ background:var(--bg); color:var(--ink); font-family:var(--sans); font-size:16px; line-height:1.55; padding-inline:16px; padding-block:32px 64px }}
main {{ max-width:1100px; margin-inline:auto; display:flex; flex-direction:column; gap:40px }}
h1,h2,h3 {{ font-family:var(--serif); font-weight:600; text-wrap:balance; line-height:1.15; margin:0 }}
h1 {{ font-size:clamp(1.9rem,4vw,2.7rem) }} h2 {{ font-size:1.6rem }} h3 {{ font-size:1.15rem }}
p {{ margin:0 0 .8em; max-width:70ch }} .lede {{ font-size:1.12rem; max-width:62ch }}
.tag {{ font-family:var(--mono); font-size:.78rem; letter-spacing:.06em; text-transform:uppercase; color:var(--muted) }}
.hf {{ color:var(--hf); font-weight:600 }} .sf {{ color:var(--sf); font-weight:600 }}
section {{ display:flex; flex-direction:column; gap:16px; min-width:0 }}
figure {{ margin:0; background:var(--panel); border:1px solid var(--rule); border-radius:6px; padding:12px; display:flex; flex-direction:column; gap:8px; min-width:0 }}
figure img {{ width:100%; height:auto; background:#fff; border-radius:3px }}
figcaption {{ font-size:.92rem; color:var(--muted); max-width:90ch }}
.answers {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(min(100%,300px),1fr)); gap:16px }}
.answers > div {{ background:var(--panel); border:1px solid var(--rule); border-radius:6px; padding:16px; min-width:0 }}
.answers ul {{ margin:.5em 0 0; padding-left:1.1em }} .answers li {{ margin-bottom:.4em }}
.scroll {{ overflow-x:auto; background:var(--panel); border:1px solid var(--rule); border-radius:6px }}
table {{ border-collapse:collapse; width:100%; font-size:.86rem }} th,td {{ padding:6px 10px; text-align:left; border-bottom:1px solid var(--rule); vertical-align:top }}
th {{ font-weight:600; color:var(--muted); white-space:nowrap }} td.n {{ font-family:var(--mono); font-variant-numeric:tabular-nums; white-space:nowrap }}
.ci {{ color:var(--muted) }}
.note {{ background:var(--note); color:var(--noteink); border-radius:6px; padding:12px 16px; max-width:90ch }}
.note p {{ margin:0 0 .5em }} .note p:last-child {{ margin:0 }}
code {{ font-family:var(--mono); font-size:.85em }}
pre.v {{ white-space:pre-wrap; font-family:var(--mono); font-size:.8rem; background:var(--panel); border:1px solid var(--rule); border-radius:6px; padding:12px; margin:0 }}
</style>
<main>
<header style="display:flex;flex-direction:column;gap:12px">
<span class="tag">HF Low Project &middot; ERA5 proxy, pipeline A &middot; seasons 2004-05 to 2025-26</span>
<h1>What a hurricane-force low looks like, against a storm-force low</h1>
<p class="lede"><span class="hf">HF lows</span> (gust index at or above 71.7 kt, 2,000 storms) are compared with <span class="sf">storm-force-only lows</span> (peak gust index 54 to 71.7 kt, 8,136 storms), Atlantic and Pacific separately, matched by basin and month. All of it is the ERA5 reanalysis used as a proxy for the archive, and ERA5 reads low in the strongest storms.</p>
</header>

<section>
<h2>Plain answer</h2>
<div class="answers">
<div><h3>Before formation</h3><ul>
<li>HF lows form under a slightly deeper trough and a stronger 250 hPa jet than storm-force lows do. In the Atlantic, Z500 over the box around Greenland and Iceland is 7 m lower and the jet 1.1 m/s stronger over days -10 to -4, and 22 m lower and 1.9 m/s stronger on the onset day.</li>
<li>The difference is small next to the weather itself. Over days -10 to -4 it is real in both basins (q&lt;0.01 for Z500 and jet) but only a few metres of Z500, and it grows steadily from about day -7 to the onset day.</li>
<li>Sea-surface temperature does not differ. Column water vapour does not differ days ahead; at onset it is lower in the Atlantic box (-0.28 kg/m2) and not different in the Pacific.</li></ul></div>
<div><h3>At the storm</h3><ul>
<li>An HF low is deeper and tighter. Central MSLP is 971 against 984 hPa in the Atlantic and 974 against 985 hPa in the Pacific, the pressure gradient within 500 km is about 60% steeper, and the strongest gust sits closer in (302 against 407 km Atlantic, 233 against 333 km Pacific).</li>
<li>The wind area is wider at onset: the 48-kt gust area reaches about 155 km farther out. Matched on central pressure, the Atlantic difference is +39 km at onset (q 0.041) and +42 km at peak (q 0.024); the Pacific one disappears (-0.4 and +4 km).</li>
<li>Gusts are a little peakier relative to sustained wind (gust factor 1.51 against 1.47 Atlantic, 1.52 against 1.48 Pacific) and the air within 500 km is moister.</li>
<li>No difference was detected in where the wind sits relative to the track: about three quarters of the 48-kt area is right of motion in both groups (detectable only above about 0.2 SD).</li></ul></div>
<div><h3>What this cannot say</h3><ul>
<li>Differences in gust and 10 m wind between the groups are true by construction, since HF is defined by the gust index.</li>
<li>The groups are at different points of life at the anchor (HF onset about 34 to 37 h after the track starts, storm-force peak 29 to 33 h) and the HF onsets sit 2.6 degrees (Atlantic) and 5.7 degrees (Pacific) farther south than the storm-force peaks, and the dewpoint is the raw value, not an anomaly, so the moisture difference is confounded with latitude. A larger share of the HF anchors' 1,200 km circle is open ocean (owned-ocean fraction 0.13 to 0.15 higher), which affects the gust-area measures.</li>
<li>These are composites, not forecasts or causes. 2 m temperature was not pulled.</li></ul></div>
</div>
</section>

<section>
<h2>The pattern before formation</h2>
<figure>{img("large_typical_atl.png", "North Atlantic mean Z500, jet and MSLP for HF and storm-force lows")}
<figcaption>North Atlantic, typical pattern. Top: before an HF low forms. Bottom: before a storm-force low peaks. White lines are 500 hPa height, colour is 250 hPa zonal wind, right-hand panels are sea-level pressure. Fields are 5.6 degree daily means, so they show the large pattern only. Both groups look alike: a strong zonal jet over eastern North America entering the ocean basin and the Icelandic low. The HF lows have a stronger jet and a deeper trough on the onset day.</figcaption></figure>
<figure>{img("large_diff_atl_C1_L.png", "North Atlantic difference maps days -10 to -4")}
<figcaption>North Atlantic, HF minus storm-force, days -10 to -4 before the anchor. Dots mark cells where the difference passes the multiplicity test (season-block bootstrap, q&lt;0.05). Heights are lower from Greenland to Scandinavia and higher to the south, and the jet is stronger along 40 to 55N. The green outline is the test box fixed before the run.</figcaption></figure>
<figure>{img("large_diff_atl_C1_k0.png", "North Atlantic difference maps onset day")}
<figcaption>North Atlantic, onset day. The same dipole is three times stronger.</figcaption></figure>
<figure>{img("large_diff_pac_C1_L.png", "North Pacific difference maps days -10 to -4")}
<figcaption>North Pacific, days -10 to -4. A weaker version of the same signal: lower heights to the north and a stronger jet.</figcaption></figure>
<figure>{img("large_diff_pac_C1_k0.png", "North Pacific difference maps onset day")}
<figcaption>North Pacific, onset day.</figcaption></figure>
<figure>{img("large_series_C1.png", "Lead-up of box-mean differences")}
<figcaption>How the difference builds up. Box-mean HF minus storm-force anomaly by day before the anchor, with 95% season-block bootstrap intervals. The shaded band is the lead window days -10 to -4. In the Atlantic Z500 becomes distinguishable around day -7, which is the lead a forecaster could use, and in the Pacific it holds near -9 to -11 m from day -4 on.</figcaption></figure>
<div class="scroll">{box_table()}</div>
<p style="font-size:.9rem;color:var(--muted)">20 pre-registered tests (2 basins x 5 fields x 2 times), Benjamini-Hochberg across the 20; q is shown. MDE is the difference detectable at 80% power (2.8 x bootstrap standard error). A difference below the MDE with q above 0.05 is "can't tell". Anomalies are against the 2004-2025 calendar-day climatology.</p>
</section>

<section>
<h2>The storm itself</h2>
<figure>{img("storm_comp_atl.png", "North Atlantic storm-relative composites")}
<figcaption>North Atlantic, composites rotated so the storm moves to the right (+x), 0.25 degree. Columns: HF at onset, HF at peak, storm-force at peak, then the two differences with dots where the difference passes the pixelwise test. The HF low has its strongest wind to the right of the track, a deeper and more compact low, and moister air within 500 km. ERA5 reads low in the strongest storms.</figcaption></figure>
<figure>{img("storm_comp_pac.png", "North Pacific storm-relative composites")}
<figcaption>North Pacific, same layout. The structure is the same as in the Atlantic, with a larger moisture difference (HF onsets are farther south).</figcaption></figure>
<h3>Structure measures, HF onset against storm-force peak</h3>
<div class="scroll">{scal_table(c1)}</div>
<p style="font-size:.9rem;color:var(--muted)">Pre-registered family of 28 tests (7 measures x 2 basins x 2 comparisons), BH over all 28; 23 pass q&lt;0.05. Sample 400 HF and 400 storm-force lows per basin, matched by month. Differences in maximum gust and 10 m wind are left out of the family because they exist by construction. The intervals use 2,000 season-block bootstraps, so q cannot go below about 0.001.</p>
<h3>At the same depth (HF peak against storm-force peak, matched on central MSLP)</h3>
<div class="scroll">{scal_table(s1)}</div>
<p style="font-size:.9rem;color:var(--muted)">Secondary, own family. The tighter gradient, smaller radius of maximum gust, higher gust factor and moister air remain when central pressure is matched. The wider 48-kt area remains in the Atlantic (+42 km) and not in the Pacific (+4 km). The same match at HF onset gives +39 km (q 0.041) and -0.4 km.</p>
</section>

<section>
<h2>Design, tests and checks</h2>
<div class="note"><p><b>Pre-registered</b> before any field was read (commit 1d11f8c): groups, the 54 kt lower bound, anchors, the 20 environment tests and 28 structure tests, and the decision rules. <b>Deviation, logged after the first run:</b> the SST box mean was missing because land cells are empty; it is now averaged over ocean cells only, which changes no other result.</p>
<p><b>Counts:</b> 12 of the 20 environment tests and 23 of the 28 structure tests pass q&lt;0.05. Secondary: 25 of 40 earlier-lag tests, 18 of 24 depth-matched, 23 of 28 with the stronger storm-force group (60 kt and above), 12 of 14 without Greenland and Iceland.</p>
<p><b>Not accounted for:</b> pixelwise maps are tested pixel by pixel and are strongly spatially correlated, so the dotted areas are generous. Storms that occur within the same week share an environment; the season-block bootstrap reduces but does not remove that. Short-lived storms have no heading and are missing from the rotated composites (about 10% of HF and 12% of storm-force storms). The lower bound of 54 kt is a choice.</p></div>
<p class="tag">Verification</p>
<pre class="v">{html.escape(verif)}</pre>
<p style="font-size:.9rem;color:var(--muted)">Pull: 16.7 GB large-scale (5.625 degree daily, WeatherBench2 and ARCO-ERA5) plus 27.2 GB storm-scale (ARCO-ERA5 0.25 degree, 1,894 times), 43.9 GB in total, under the 50 GB gate. Code and tables: <code>research/era5/hf_vs_storm/</code> on the research branch.</p>
</section>
</main>
'''
open(out, "w").write(page)
print(len(page) / 1e6, "MB")
