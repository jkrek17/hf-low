# Storm review template (v2, preview)

One standalone page per storm, `storm-review-<track>.html`, linked from the Storm reviews index. Everything is an ERA5 proxy on pipeline A tracks. Text is generated from the numbers, never typed per storm.

| # | Section | Answers | Filled from |
|---|---|---|---|
| 0 | Header tiles | How big was it? | lowest pressure, fastest 24 h fall, peak gust index, hours at HF strength, peak ocean area at 71.7 kt+ (nm2) |
| 1 | Storm in plain words | What happened, in one paragraph | pipeline A fixes (existing text) |
| 2 | Track and wind field | Where did it go, where was the wind | track coloured by pressure, 48 h before to 24 h after gust swath (ARCO-ERA5 0.25 deg) |
| 3 | How it developed | When did it deepen, and what was its structure | six-panel series: pressure, 24 h fall vs bomb yardstick, gust index and 71.7 kt area, 250 hPa jet, 300 hPa divergence and 500 hPa vorticity advection, Hart B and lower-level core with structure class strip (1.5 deg) |
| 4 | Upper-level support | What was driving it aloft | SLP / 500 hPa height / 250 hPa jet maps every 12 h; mean jet, divergence, vorticity advection, Eady growth and water vapour over the 24 h of fastest fall |
| 5 | Low-level structure | Where were the fronts and the low-level jet | 850 hPa temperature, SLP, 925 hPa 50/64 kt contours; ground above 1,000 m masked |
| 6 | Where the strongest wind was | Which part of the storm had the gusts | Hart phase-space path; storm-relative map of 64 and 71.7 kt+ ocean cells; split by front-right, rear-right, rear-left, front-left; peak area and strongest gust location |
| 7 | Set-up compared with other HF-strength lows | Was the set-up unusual | the storm's median jet, divergence, vorticity advection, Eady rate, water vapour, SST gradient at its HF-strength fixes, ranked against all HF-strength fixes of the basin since 2004-05 (`storm_review_ref.csv`) |
| 8 | What OPC had (archive era only) | Did OPC carry it, and when | OPC archive event matched by time and position within 200 nm: fixes, HF hours, lowest pressure, overlap with the proxy |
| 9 | Read with care | What not to over-read | fixed wording per storm era |

Wording rules: "OPC also carried it", never "verified". Nautical miles. Pre-2004 storms say "not direct observation". The structure class is Hart-style at 1.5 deg with 900/850/700/600/500/400/300 hPa levels, so it is a guide.

Code: `storm_review_extract.py` (pull, about 1.1 GB a storm), `storm_review_env.py` and `hart_cps.py` (1.5 deg diagnostics), `storm_review_v2.py` (figures and page).
