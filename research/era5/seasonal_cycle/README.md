# Seasonal cycle of hurricane-force lows, and what sets the timing of each basin's peak

Plan, committed before any outcome was computed: [PREREGISTRATION.md](PREREGISTRATION.md) (commit `d5c9f44`).
Code: [extract_eddy.py](extract_eddy.py) (ERA5 pull), [seasonal.py](seasonal.py) (all 56 tests),
[make_report.py](make_report.py) (tables, figure), [post_hoc.py](post_hoc.py) (checks run after seeing results).
Results: [results/](results/): `summary.txt` (readable), `tests.csv` (every test), `monthly_tables.csv`, `post_hoc.csv`,
`within_era_1979_2000.csv`, `seasonal_cycle.png`.

Reproduce (about 2 min of ERA5 pull, then about 3 min; the pull is about 34 GB for v250 and 4 GB for MSLP, under the 50 GB gate):

    python3 research/era5/seasonal_cycle/extract_eddy.py v250 2004-06-01 2022-01-01
    python3 research/era5/seasonal_cycle/extract_eddy.py mslp 1979-06-01 2022-01-01
    python3 -I research/era5/seasonal_cycle/seasonal.py . research/era5/seasonal_cycle/work research/era5/seasonal_cycle/results 2000
    python3 -I research/era5/seasonal_cycle/make_report.py research/era5/seasonal_cycle/results
    python3 -I research/era5/seasonal_cycle/post_hoc.py . research/era5/seasonal_cycle/work research/era5/seasonal_cycle/results

Counts are **pipeline A** (800 km ocean gust index, HF-equivalent at 71.7 kt) and a **proxy**, or the OPC archive where
labelled. Eddy activity is ERA5 reanalysis (1.5 degree, 6-hourly), a different use of ERA5 from the gust index. The pipeline A
tracks and the archive are Oct-Apr dominated; "season" is 1 June to 31 May and the sample is seasons, not storms.

## Answer

**Framing (changed after the results were computed; see Departures).** Jason's steer arrived after the first run: start from the
observed fact that HF lows peak in January in the Pacific and February in the Atlantic, describe the cycle and what sets its
timing (cyclone count, HF share, latitude, jet and baroclinicity), and treat Nakamura's band-pass eddy-variance suppression
as one neutral secondary comparison. The numbers below are the same ones the pre-registered plan produced; only the lead
and the reading changed, and that change was made with the results in view.

**Timing of the peaks (post hoc, `post_hoc.csv` PH5; bootstrap share of draws in which each month holds the per-day maximum).**
- Pacific, archive: January 81%, December 19%, February 0%. Pipeline A (genesis month): December 63%, January 37%.
  So January is the archive's peak but the proxy cannot separate December from January: a December-January plateau.
- Atlantic, archive: February 58%, January 37%, December 5% (per day; by raw monthly count January 9.36 is above February 8.77
  because February is short). Pipeline A: February 54%, January 34%. The Atlantic February maximum is the point estimate, not
  established against January.
- What carries the timing: the cyclone count is nearly flat across the winter in both basins (it does not set the peak); the HF
  share rises from about 7-8% in October-November to 11.4% (Pacific, December and January) and 15-16% (Atlantic, December to
  February, February highest at 15.7%) and falls in March. In the Pacific the share peak is a December-January plateau
  (December 47%, January 48%, February 5% of draws); in the Atlantic it is February (50%).
- Latitude and jet/baroclinicity: latitude below. Jet and baroclinicity were not examined (no u-wind or Eady fields pulled;
  not in scope of this run).

**Secondary, neutral comparison with eddy variance.** Does 2-6 day eddy variance dip while HF counts do not? In the Pacific
v250 box, yes: variance is lowest in January while HF counts and share are highest in December-January. In the Atlantic it
does not dip. Details in items 2, 3 and 5.

1. **The Pacific January peak is a December-January plateau, and it is firm only in the archive.** Per day, January is above the
   mean of December and February by 22% in the archive (95% interval +4% to +41%, 22 seasons) but by 10% in the ERA5 proxy
   (-4% to +26%): "archive only" under the plan's rule. The excess comes from February being low (archive Feb/Jan: -27%, interval
   -42% to -9%; January against December alone: +7%, -8% to +24%, post hoc).
2. **Our ERA5 reproduces the Pacific lull in eddy activity** (2-6 day v' variance at 250 hPa, 30-60N, 150E-150W, 17 seasons):
   January is 12% below the mean of December and February (-22% to -2%) and December-February is 16% below the November and March
   mean (-25% to -5%). The surface pressure variance shows less of it (January -7%, -13% to 0%; December-February vs shoulders
   +4%, -1% to +10%, not distinguishable from zero). The Atlantic shows no lull, as the literature says (v250: January -4%, -17%
   to +9%; December-February vs shoulders +8%, -1% to +19%).
3. **HF counts are in antiphase with that cycle.** Pacific HF counts are 64% (archive, +37% to +103%) to 72% (proxy, +47% to
   +107%) higher in December-February than in November and March; eddy activity is 16% lower. On the same 17 seasons the HF
   cycle departs from the eddy cycle by a factor of 2.0 in the December-February contrast (1.6 to 2.6) and by +23% for January
   (+2% to +44%).
4. **What reconciles them, by the plan's own rules:**
   - **R1, fewer cyclones but a higher HF share: not supported on the pre-registered January contrast** (share +6%, -7% to
     +22%) and the cyclone count does not dip (+4%, -2% to +9%). Read on the December-February contrast (pre-registered as a test,
     but the plan's verdict rule was written for January, so this reading is post hoc) the share does carry the excess: HF share
     is 5.8% in March, 8.3% in November and 11.4% in both December and January, and the share term is 84% of the log excess
     (75% to 92%). Cyclone counts rise only 9% (5% to 13%).
   - **R2, storms farther south: partial.** Archive first-HF fixes are 2.3 degrees farther south in January than in the mean of
     December and February (1.3 to 3.3 degrees; January was south of the December-February mean in 17 of 21 seasons that have all
     three months, post hoc). The proxy's HF onset latitude moves only 0.5 degrees (-1.8 to +0.6), not distinguishable from zero.
     The v250 storm-track axis (descriptive) moves from 45.1N in November to 41.9N in January.
   - **R3, eddy activity and extremes differ: supported by the rule** (paired contrast above, interval above zero), and the
     archive version, run post hoc, agrees (+37% for January, 11% to 64%; factor 1.9 for December-February, 1.4 to 2.6).
5. **It is not Pacific-specific.** The Atlantic HF cycle is also far steeper than its eddy cycle (December-February vs shoulders:
   +74% in the proxy, against +8% for v250 variance; paired contrast 0.51 in logs). What is special in the Pacific is that the
   mean eddy activity actually dips while HF counts do not.

**What this does not show.** It does not show why a 16% lower eddy variance can sit with a 60-70% higher HF count. The
data say the HF share rises sharply in midwinter; no jet, Eady growth or moisture field was examined (see the plan's scope).
A tail-sensitivity argument (the HF threshold sits far in the tail of the intensity distribution, so a modest rise in the
supply of baroclinic energy can raise the share a lot) is plausible and untested here.

## Monthly tables (mean per season, Oct-Apr shown; all twelve months in `results/summary.txt`)

Archive and pipeline A: 2004-05 to 2025-26 (22 seasons). "Cyclones" in every table and test means all pipeline A tracks (a low below 1010 hPa for at least
24 h with at least two fixes in the basin domain), counted in its genesis month; the share is HF / cyclones in that cohort.

| Pacific | Oct | Nov | Dec | Jan | Feb | Mar | Apr |
|---|---|---|---|---|---|---|---|
| Archive HF lows (first HF fix month) | 3.95 | 5.14 | 8.00 | 8.59 | 5.68 | 4.00 | 1.50 |
| Pipeline A HF lows (genesis month) | 4.41 | 5.23 | 8.45 | 8.27 | 6.14 | 3.77 | 1.27 |
| Pipeline A cyclones | 60.14 | 62.86 | 74.27 | 72.64 | 60.23 | 65.23 | 65.96 |
| Pipeline A HF share, % | 7.33 | 8.32 | 11.38 | 11.39 | 10.19 | 5.78 | 1.93 |
| v250 eddy activity, m2/s2 (17 seasons) | 121.05 | 128.38 | 111.53 | 95.52 | 106.67 | 119.68 | 133.05 |

| Atlantic | Oct | Nov | Dec | Jan | Feb | Mar | Apr |
|---|---|---|---|---|---|---|---|
| Archive HF lows (first HF fix month) | 3.82 | 4.50 | 8.45 | 9.36 | 8.77 | 5.95 | 1.77 |
| Pipeline A HF lows (genesis month) | 4.05 | 4.50 | 9.32 | 9.86 | 9.27 | 6.82 | 1.95 |
| Pipeline A cyclones | 58.73 | 62.32 | 62.23 | 64.77 | 59.00 | 67.27 | 70.18 |
| Pipeline A HF share, % | 6.89 | 7.22 | 14.97 | 15.23 | 15.72 | 10.14 | 2.78 |
| v250 eddy activity, m2/s2 (17 seasons) | 106.73 | 133.15 | 132.25 | 128.79 | 136.34 | 112.78 | 109.72 |

The cyclone count column is the whole population of weak and strong lows, so it describes how many lows there are, not how
intense. Depth-based (pressure) counts, the real-cyclone count, onset-month counts and latitude tables are in `summary.txt`.
Figure: [results/seasonal_cycle.png](results/seasonal_cycle.png).

## The literature's prediction, and what was found

| | Predicted | Found (Pacific) | Found (Atlantic) |
|---|---|---|---|
| v250 eddy, January vs Dec and Feb | below | -12% (-22% to -2%), q = 0.054 | -4% (-17% to +9%) |
| v250 eddy, Dec-Feb vs Nov and Mar | below | -16% (-25% to -5%), q = 0.014 | +8% (-1% to +19%) |
| MSLP eddy, January vs Dec and Feb (secondary) | same sign, weaker | -7% (-13% to 0%), q = 0.105 | +3% (-4% to +12%) |

## Tests, multiplicity, power

Family of 56 (14 quantities x 2 basins x 2 contrasts); Benjamini-Hochberg q is across all 56; bootstrap p is two-sided from
2,000 season-block draws (smallest possible 0.001). The pre-registered primary set (Pacific January contrast, 7 tests) with
Bonferroni x7:

| test | estimate (log) | 95% interval | p | BH q | Bonferroni x7 |
|---|---|---|---|---|---|
| Q1 archive HF count | +0.196 | 0.038 to 0.346 | 0.014 | 0.036 | 0.098 |
| Q2 pipeline A HF count | +0.092 | -0.040 to 0.227 | 0.184 | 0.289 | 1.0 |
| Q3 all cyclones | +0.036 | -0.022 to 0.088 | 0.221 | 0.326 | 1.0 |
| Q4 HF share | +0.056 | -0.076 to 0.198 | 0.424 | 0.505 | 1.0 |
| Q9 v250 eddy activity | -0.133 | -0.242 to -0.020 | 0.023 | 0.054 | 0.161 |
| Q11 archive first-HF latitude (degrees) | -2.277 | -3.342 to -1.277 | 0.001 | 0.003 | 0.007 |
| Q14 paired HF minus eddy | +0.204 | 0.018 to 0.362 | 0.038 | 0.085 | 0.266 |

So the archive peak, the eddy lull and the HF-minus-eddy contrast are each nominally significant and none survives
Bonferroni; only the latitude shift does. On the December-February contrast (M) the same quantities have BH q of 0.003 (HF
counts, share), 0.014 (v250), 0.003 (paired). Of the 56 tests, the whole list with leave-one-season-out ranges is in
`tests.csv`.

**Power.** The January contrast on a count has bootstrap SE 0.07-0.08, so a 80%-power detectable |D| of about 0.19-0.22 (about
a 21-25% departure). The proxy's +10% January excess is therefore neither established nor excluded: effects from a 4% dip to
a 26% rise fit. The v250 eddy lull is detectable at about 0.16 (17%) and measured at 0.13. Effective n is the number of
seasons (22, 17, 42, 47), not the storm count.

## Limits

- **The ERA5 proxy is not the archive.** The proxy's January excess is smaller than the archive's and it is weaker in latitude.
  The agreement in the December-February contrast (+64% archive, +72% proxy) is the safe comparison.
- **Eddy activity ends in 2020-21** (the 1.5 degree dataset stops on 2021-12-31) and v250 covers 17 seasons, not the 22 of the
  HF counts; the v250 pull was restricted to 2004-05 on to stay under 50 GB. MSLP covers 42 seasons.
- **Fixed box.** The Pacific box (30-60N) is not the storm-track axis and the axis moves south in midwinter; a latitude-following
  measure was only described (axis latitude), not tested.
- **The genesis-month cohort** differs from the archive's first-HF-fix month for storms born late in a month; both are in
  `summary.txt` (Q2 and Q6 agree: +0.092 and +0.102).
- **Archive detection and practice** could depend on the month and cannot be tested here.
- **Domain edges** (20-75N, fixed) are not handled, as in `freq_split`; the Pacific storm track moves south in midwinter, so
  some of the Pacific cyclone count change could be edge crossings.
- **Gust before 2004.** The 1979-80 to 1999-2000 gust-based Pacific shape (`within_era_1979_2000.csv`, descriptive,
  within-era only, not in the family) has HF January contrast -0.036 (-0.197 to +0.125) and December-February
  contrast +0.345 (0.186 to 0.515): no January peak, and a weaker midwinter excess than in 2004-05 on. Not tested against the
  later era.

## Departures from the plan

None changes a pre-registered estimate.

0. **Reframing after the data were seen.** The plan (`d5c9f44`) was built around reconciling a January peak with a midwinter lull. Jason's relayed steer (2026-10-08 09:12Z, "there is no storm track lull; Pacific peaks in January and Atlantic in February") reframes the work as a description of the cycle and its timing, with eddy variance as a neutral secondary comparison. It reached this thread after the pre-registered family had been run once, so it is **not** a change made before seeing data. The test family and verdict rules are unchanged; the peak-month analysis (PH5) and the lead text are post hoc. The v250 eddy lull is a statement about ERA5 band-pass variance only, not about HF counts, where Jason's reading (no lull) is what the counts show.

1. The H0-peak rule said "if the peak is not established, the reconciliation questions are not interpreted." It was established in
   the archive and not in the proxy ("archive only"). The reconciliation hypotheses were therefore read as written, and the
   proxy-based ones (R1, R3) are reported with that caveat; the archive version of R3 (PH1) is post hoc.
2. The R1 verdict rule used the January contrast. The December-February reading and the share fraction f are post hoc.
3. MSLP eddy activity is in Pa^2 (the dataset's units), not hPa^2; a label only.
4. Run once with 2,000 draws; no debug run was looked at. The post hoc script re-runs the whole primary script (the results
   file is byte-identical, `tests.csv` md5 `a71410fe...`), so the post hoc rows use the same bootstrap draws.
5. The archive latitude (Q11) depends on a tie rule the plan did not state: 15 Pacific HF rows share an ID and time with different positions; the committed run keeps file order among ties (stable sort). Other tie rules move Q11 D_Jan from -2.277 to -2.267 or -2.303 (verifier).
6. `data/hf_lows/HF_Data_-_Pac.csv` has one invalid date (`20241101018`, 11 digits); the row was dropped. It belongs on the
   sheet-fix list.

## Verification

A fresh Sonnet agent recomputed from the raw inputs, without reading the code or results: all monthly counts and shares quoted, archive totals (856 Pacific, 1005 Atlantic), the Pacific D_Jan and M point estimates for every quantity quoted (to 3 decimals), the Atlantic eddy and archive contrasts, the paired contrast (+0.204, +0.709), the v250 monthly means, and the 95% bootstrap intervals (own seed; within 0.013, latitude within 0.01 degrees). It also re-read three raw v250 blocks from the zarr: identical. All matched. **Not independently checked:** the BH q-values, Bonferroni values, leave-one-season-out ranges, the Atlantic paired contrast, MSLP monthly means, Q5 to Q8, Q10, Q12, Q13 tables, the 1979-2000 within-era rows, and every post hoc row (PH1 to PH5). The verifier flagged that "all cyclones" means all pipeline A tracks (relabelled) and the tie rule above.

## Post hoc checks (`results/post_hoc.csv`; not in the 56-test family or its FDR)

- PH1 archive HF count minus v250 eddy, 17 seasons, Pacific: D_Jan +0.314 (0.100 to 0.492), M +0.622 (0.336 to 0.938).
  Atlantic: D_Jan +0.076 (-0.144 to +0.323), M +0.491 (0.352 to 0.635).
- PH2 share fraction f = contrast(share) / contrast(HF), pipeline A: Pacific M 0.842 (0.747 to 0.919); Pacific D_Jan 0.611
  (-2.8 to 4.4, uninformative); Atlantic M 1.053 (0.968 to 1.130).
- PH3 archive first-HF latitude, Pacific, per season: January south of the mean of December and February in 17 of 21
  seasons; median difference -2.07 degrees.
- PH5 which month holds the maximum: see the timing section above.
- PH4 Pacific HF count log-rate differences: January minus December +0.071 (-0.078 to 0.217) archive and -0.022 (-0.163 to
  0.115) proxy; February minus January -0.320 (-0.546 to -0.091) archive and -0.205 (-0.409 to -0.029) proxy.
