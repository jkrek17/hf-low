# Near-storm upper-level environment: 24 h gale / storm / HF outlook and bomb-environment diagrams

ERA5 proxy, pipeline A tracks (`../intensity`). Written 2026-10-13 before any 4-class gust outcome or any new
model was scored. Looked at before writing: the committed `intensity/results/skill.txt` (old classes and HF24),
column names and row counts of `fixes_2004.csv.gz` and `env_2004.csv.gz`. Not looked at: the gust-class outcome
below, any fit of the models below.

## Question (Jason, 2026-10-13)
"A near storm environment check, like the Hart phase space, that looks at the upper-level parameters and tells if the
storm is in an environment for bomb, or otherwise indicate the growth rate." Then: "does it best fit a gale vs storm
vs HF for growth potential over the next 24 hours?" Answer to give: yes / no (well powered) / can't tell, with the gain
over what the storm's own state already tells.

## Outcomes
- **Gust outlook (4 ordered classes).** Peak gust index (`g800`, max ocean gust within 800 km, kt) over the next
  24 h, sampled at the next two 00/12 UTC fixes (t+12, t+24; the 6-hourly fixes are not in the committed table, so
  this differs from `hf24` in `intensity`): below gale (<38), gale (38 to <54), storm (54 to <71.7), HF (>=71.7).
  Cutoffs are nominal proxy classes on the gust index. Only 71.7 kt is tied to the archive; 38 and 54 are not
  calibrated on observed events (the archive lists HF lows only). 54 = 0.75 x 71.7 as in `hf_vs_storm`; 38 as in
  `era5-cache`. A fix whose track ends within 24 h is scored on the fixes that exist (no fix = below gale), as in
  `hf24`; sensitivity: drop fixes with `left24`.
- **Rapid deepening (binary)**, class from `intensity` (NDR >= 1 Bergeron in 24 h), as scored in PR 12.

## Sample, unit, n
Seasons 2004-05 to 2025-26 only (decision 1), 159,430 fixes, 36,598 tracks. Effective n = 22 seasons. Uncertainty:
resampling whole seasons (2,000 draws). Tropical cyclones in. Atlantic fixes >60N are lower confidence.

## Models (multinomial logistic, L2, standardised, median-imputed, as `intensity/model.py`; leave-one-season-out)
- clim: basin x month frequencies
- state: storm state (msl, dp12, young, lat, speed, g800, logage, pac, doy_c, doy_s)
- A "jet + trough": state + jet250 + vadv500
- B "jet + baroclinic": state + jet250 + eady
- full: state + Hart (B, VTL, VTU) + all environment terms = the Storm model's predictors, refit to the new target
The gauge on the site uses `full`. The two diagrams use binned observed rates (6 x 6 training-season quantile
bins of each axis pair, per basin pooled; cells under 30 fixes are blank); their skill is that of A and B.

## Tests (gain in score over state; score = RPSS vs clim for the outlook, BSS vs clim for rapid deepening)
Family of 6 (BH within the family, q < 0.05): A-state, B-state, full-state, each for the outlook and for rapid
deepening. Secondary (not in the family, reported): full-A, full-B; HF-class-only BSS; gale-only; each basin;
sensitivity without `left24`; sensitivity without `g800` in state (as `full_nogust`).
Reliability: observed vs forecast frequency by decile for each class, `full` model, pooled seasons; report the
slope. Rapid-deepening bins of the diagrams: shown with n and Wilson intervals.

## Decision rules
- Yes: BH q < 0.05 and the 90 % season-bootstrap interval of the gain excludes 0.
- No: the interval upper bound is below 0.005 RPSS/BSS (the smallest gain worth a forecaster's attention; stated
  now, not tuned) -- then report as a well-powered null with that bound as the MDE.
- Otherwise can't tell.
Prior expectation (from PRs 12, 90, 149, 154): small gains from upper-level terms (RPSS about +0.01 to +0.03 for
full over state); jet and trough axes alone likely smaller. The display is a communication tool even if the gain is small.

## Power
22 seasons, thousands of fixes per season; the season-bootstrap interval width reported in `skill.txt` is the
detectable gain. No planted-effect simulation (intervals are from resampling the data themselves; an upper bound on power
because fixes within a storm are not independent).

## Pull plan
None. Uses `intensity/results/{fixes,env}_2004.csv.gz` only. The era5-cache gust table is not needed.

## Known conflicts
PR 90: no gain going aloft. PR 149: jet quadrant and upstream trough small or null. The `g800` term in state already
carries current intensity, so gains are over a strong persistence baseline (full-nogust sensitivity covers this).

## Deviations (post hoc)
(none yet)
