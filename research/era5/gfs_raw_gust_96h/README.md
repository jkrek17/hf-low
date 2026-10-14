# GFS raw gust against OPC hurricane-force calls, 12-hourly leads to 96 h

Extends `/mnt/project-files/gfs-raw-gust/` (research PR 182, leads 0 to 48 h) to 72 and 96 h. **PROXY**: the pipeline A detector run on GFS PRMSL and GUST forecasts; truth is the OPC archive Category HF fix within +-3 h and 600 km of the valid time (Y0). Seasons 2021-22 to 2025-26, 00 and 12 UTC cycles, 1 Oct to 30 Apr (2,122 cycles), 5 seasons.

- `raw_gust_fixes_96h.csv.gz`: one row per GFS track fix (180,784), columns cycle, lead (0..96 step 12), valid_time, track, lat, lon, msl, g500, g800, g1200, basin, season, G (= g800, kt), Y0.
- `raw_gust_96h.json`: held-out thresholds, scores by basin and lead (matched / 64 kt / 72 kt: hits, fa, miss, cn, pod, agree = hits/(hits+fa), hss), gust-band counts by lead.
- `a_fetch.py`, `a_tracks.py`, `track.py`, `score96.py`: the code used (GFS_WORK=<dir> a_fetch.py pulls; score96.py LOWS_DIR ARCHIVE_DIR OUT_PREFIX scores). Raw fields are not kept (31.2 GB streamed from the NOAA GFS open-data bucket by byte range, 2026-10-14).

## Deviations from PR 182 (made after the card choice of 2026-10-14; not pre-registered)
Leads every 12 h (PR 182: every 6 h to 48 h); linker MAXD 1800 km per 12 h and 2400 km gate (900 and 1200); MIN_LEN 3 fixes instead of 4. Check: leads 0 to 48 reproduce PR 182 within about 0.02 (Atlantic lead 0 matched POD 0.52 vs 0.54, HSS 0.50 vs 0.52).

## Limits
Five seasons, one model with version changes inside the window; no intervals; fixes within one storm are not independent. The tracker must find the low in each forecast, so by 72 to 96 h the fix set is a different sample of lows. OPC's call is an analyst judgement, not a measurement: "agree" is not "correct". June to September not covered.

## Code layout
`a_fetch.py`, `a_tracks.py`, `track.py` are copies of `research/era5/gfs_test/a_fetch.py`, `a_tracks.py` and `research/era5/hf_history/track.py` with three edits: `a_fetch.LEADS = range(0, 97, 12)`, `track.MAXD = 1800`, `track` gate 2400 km. They expect the sibling layout `gfs_test/` and `hf_history/` (see the top of `score96.py`); to rerun, place them there or adjust the two `sys.path` lines. Written 2026-10-14; the raw GFS fields are not kept, so reproducing the pull needs the NOAA open-data bucket (31.2 GB).

## Verification
A fresh Sonnet agent recomputed the fix count, gust-band shares and the Atlantic scores from `results/raw_gust_fixes_96h.csv.gz` without reading `score96.py`: they match. Pacific agreement and skill cells differ by up to 0.004 (tie convention at the matched threshold: this code uses the (k+1)-th largest as in PR 182, the verifier the k-th). Not independently checked: the pull itself, the tracker changes.
