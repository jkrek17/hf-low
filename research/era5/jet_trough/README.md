# Jet speed and upstream trough depth against bomb and hurricane-force onset (ERA5 proxy, pipeline A)

Question (Jason, 2026-10-08): is there a threshold jet speed and upper-level trough depth that tilts a
cyclone toward being a bomb or reaching hurricane force? Did we come up with the probability of HF yet?

**This is a proxy.** "HF" is pipeline A's gust index reaching 71.7 kt in ERA5, not the archive.
Pipeline A (`../hf_history`), the near-storm intensity framework's fixes (`../intensity`),
159,430 low fixes at 00/12 UTC, seasons 2004-05 to 2025-26. Plan committed before any outcome was joined:
[PREREGISTRATION.md](PREREGISTRATION.md). Fit seasons 2004-05..2014-15, held out 2015-16..2025-26, per basin.

## Answer in plain words

**Jet speed: yes, it matters a lot, but there is no sharp threshold. The odds rise along a smooth ramp.**

- A bomb (24 hPa in 24 h scaled to 60 N) needs a strong jet. In the held-out seasons, bombs are rare below
  about 100 kt of 250 hPa wind within 1000 km of the low (0.3% of lows in the Atlantic, 0.1% in the Pacific),
  8% at 120 kt and above, 14% at 140 kt and above, 20% at 160 kt and above, in both basins.
- At fixed storm state (pressure, pressure tendency, latitude, season, and current gust for HF) each extra
  standard deviation of jet (Atlantic 27 kt, Pacific 34 kt) multiplies the odds of a bomb by 2.9 (Atlantic) and 2.7
  (Pacific): odds double about every 18 kt in the Atlantic and 24 kt in the Pacific. Tests of this against
  nothing: q = 0.002 in each basin.
- A one-knot (threshold) curve does not fit better than a straight line on the logit scale on held-out
  seasons for bombs (q = 0.97 Atlantic, 0.17 Pacific), and a smooth spline beats it (p = 0.005 and 0.001). The
  data show a ramp, not a break. Where a one-knot fit puts its knot (Atlantic 129 kt, 90% interval 126-131;
  Pacific 153 kt, 133-158), the slope flattens above it, which is saturation, not a threshold.
- **HF onset**: the jet matters too, less after the storm's current state is accounted for. Odds per SD: 1.29
  Atlantic, 1.48 Pacific (q = 0.002 each). In the Atlantic the curve is nonlinear (one knot beats a line,
  q = 0.002; knot 105 kt, interval 79-120) and a kink cannot be told apart from a sharp smooth curve
  (see "Power" below). In the Pacific it is a smooth ramp.

**Trough depth: not detectable once storm state is accounted for.** Upstream 500 hPa trough depth, measured
as height below the monthly mean 500-2500 km to the west, has adjusted odds per SD of 1.23 and 1.24 for bombs
(Atlantic, Pacific) and 1.15 and 0.90 for HF onset. None passes FDR (q = 0.14, 0.34, 0.67, 0.19). Raw
frequencies do rise with trough depth (Atlantic HF within 24 h: 2.5%, 3.1%, 7.4% across trough terciles), but
the deep-trough storms are already deeper and nearer HF, so the rise is mostly storm state. This is **not a
well-powered null**: a per-SD odds ratio of 1.15 is detected 21-41% of the time (see Power). The trough's point estimates (at most 1.24 per SD) are far below the jet's for bombs
(2.7-2.9 per SD), but effects around 1.15-1.25 per SD are neither found nor excluded.

**Jet and trough together: no interaction detected.** Odds ratio of the product term 0.96, 0.95 (Atlantic bomb,
HF), 1.00, 0.99 (Pacific), q = 0.38, 0.33, 0.23, 0.21. Effects add on the logit scale as far as these data
resolve; the interaction test was not power-simulated.

**Joint probability, held-out observed frequency.** Terciles from the fit seasons. Bomb (class-eligible fixes) / HF
within 24 h (all fixes, the framework's target). Full 3 x 3 table with four-way splits in `results/grid_terciles.csv`.

| Basin | Jet tercile | Bomb | HF within 24 h | HF onset (not already HF) | Framework mean P(HF 24 h) |
|---|---|---|---|---|---|
| Atlantic | below 95 kt | 0.20% | 1.5% | 0.77% | 1.3% |
| Atlantic | 95-119 kt | 1.4% | 2.4% | 1.4% | 2.3% |
| Atlantic | above 119 kt | 7.8% | 8.6% | 6.0% | 8.1% |
| Pacific | below 92 kt | 0.09% | 0.50% | 0.19% | 0.53% |
| Pacific | 92-124 kt | 1.1% | 1.1% | 0.62% | 1.4% |
| Pacific | above 124 kt | 9.3% | 7.2% | 4.9% | 7.0% |

Within the top jet tercile the trough tercile moves the Atlantic HF-within-24 h frequency from 5.4% and 6.2% to
12.5% and bomb not at all (8.3%, 7.0%, 8.0%); in the Pacific HF goes 6.5%, 5.9%, 8.9% and bomb 10.9%, 8.8%, 8.5%.
The Atlantic jump at the deepest trough is raw, not adjusted: it does not survive storm state (q = 0.67).

## Did we come up with the probability of HF?

Yes: the near-storm framework (hf-low PR 12, `../intensity`) gives P(HF within 24 h) from storm state, Hart's phase
space and the environment, jet speed included. Trained on 2004-05..2014-15 and scored on 2015-16..2025-26
(binary version of its `full` model, pooled, this PR's `results/added_skill.csv`): Brier skill score against
basin-month climatology 0.426 (90% interval 0.410-0.442), POD 0.60, FAR 0.38, HSS 0.59, bias 0.97. For bombs:
BSS 0.292, POD 0.49, FAR 0.52, HSS 0.47. Its forecast matches the observed frequency in the jet terciles above
(last column) to within about 1 point. Perfect prognosis on ERA5's own fields: it has not been tested on GFS.

**Does trough depth add to it?** Slightly. Adding upstream trough depth to the framework's `full` predictors
raises held-out BSS from 0.4264 to 0.4305 (HF) and 0.2924 to 0.2960 (bomb), with q = 0.02 and 0.002 (HF,
Atlantic and Pacific) and 0.04 and 0.002 (bomb); HSS moves 0.593 to 0.594 and 0.470 to 0.474. A gain of
four thousandths of Brier skill is detectable over 11 seasons and has no practical weight. Adding the local trough
and the upstream jet as well gives 0.4309 and 0.3013.

## Power: what these data can and cannot separate

Simulated on the real fixes (`results/power.txt`, 200 simulations per row; p < 0.01 per test as a stand-in for FDR).
Truth = the fitted line plus a hinge with odds ratio 2.0 across knot +/- 0.5 SD.

- A hinge is separated from a straight line (T2) 37-90% of the time for jet (Atlantic bomb 55%, Atlantic HF 90%, Pacific bomb 47%, Pacific HF 37%), 99% for trough depth.
- A hinge is separated from a smooth 4-df spline (T3) only 4-12% of the time. The false-positive rate under a
  smooth ramp is 4-7%, so the 11-season bootstrap p is somewhat liberal.
- **A kink cannot be told apart from a sharp smooth curve with this sample.** The literal pre-registered rule can
  therefore almost never declare a threshold, and "smooth ramp" is a weaker finding than it sounds. The verdict
  wording was amended before the held-out run (PREREGISTRATION.md, deviations 2). Both labels are in
  `results/verdicts.csv`. In practice: for bombs and jet speed the one-knot form is *worse* than the spline
  and no better than a line, which is positive evidence of no sharp break. For HF onset it is not.
- A straight-line effect of 1.15 per SD on the odds is detected 21-41% of the time. "No relation detected" for trough
  depth means "no effect large enough to see; effects of 1.15-1.25 per SD are neither found nor excluded".

## All tests

One BH family of 32 primary tests (`results/fdr_table.csv`): 24 threshold tests, 4 interaction tests, 4 added-skill
tests. Passing at q < 0.05: the four jet T1 tests (q = 0.002), the Atlantic HF jet T2 (0.002), the Atlantic bomb
trough T3 (0.016: the hinge beats the spline there, but T1 and T2 fail, so it is not a threshold),
and the four added-skill tests (0.02, 0.002, 0.04, 0.002). 24 further tests on the two secondary variables (upstream jet
and the local trough, their own family): the upstream jet is detected as strongly as the jet maximum for bombs
in both basins and for Pacific HF (q = 0.003); the local trough only for Pacific bombs.

## Sensitivities (not in the family; `results/thresholds.csv`, tags)

- **Lead 24 h** (predictors at t, outcome 24-48 h later): the jet effect persists, odds per SD 1.8-2.4 (T1 p = 0.0005);
  trough depth none (odds per SD 1.03 Atlantic bomb, 1.20 Atlantic HF, 0.79 Pacific bomb, 0.91 Pacific HF).
- **October-April only**: same picture as all months (jet odds per SD 2.5-2.9 for bombs, 1.4-1.6 for HF).
- **Unadjusted** (no covariates): jet odds per SD 2.6-4.3, larger than adjusted, as it is partly the storm's own
  strength; again no hinge better than the line for bombs.
- **Eddy-height trough** (Z500 minus its zonal mean, instead of the departure from the fit-season monthly mean): a trough effect
  becomes detectable, odds per SD 1.78 (Atlantic bomb), 1.37 (Atlantic HF), 1.68 (Pacific bomb), 1.17 (Pacific HF), T1
  p = 0.0005, 0.0005, 0.0005, 0.10. This definition keeps the stationary-wave mean, so it also tells where on the globe the
  storm is (Icelandic and Aleutian lows) and the models carry latitude and season but not longitude. It is a sensitivity,
  reported as run; the primary definition did not detect an effect. The two are not reconciled here.

## Caveats

- `jet250` is the maximum 250 hPa wind within 1000 km at the fix, which includes the storm's own outflow wind, so
  part of the jet effect is storm intensity. The upstream jet (`jet_up`, secondary) gives the same T1 result in
  three of four cases.
- The trough anomaly baseline is the 2004-05..2014-15 monthly mean of Z500; any warming trend over the later seasons
  appears as a small positive anomaly shift (order 10 m against anomalies of 100 m and more).
- Predictors at t, outcomes in (t, t+24 h]: lead 0-24 h, mean about 12 h, not exactly the 12-24 h in the brief.
- The 5.625 degree grid is coarse for trough depth: a minimum over about 20 cells in a 500-2500 km sector.
- Intervals resample 11 held-out seasons, 2000 times. Fixes within a storm are not independent; the season blocks absorb that.
- Seasons before 2004-05 are not used (decision 1). HF figures before 2001 would carry the gust drift caveat.
- Bomb frequency conditions on the track surviving 24 h (the framework's class-eligible fixes).

## Independent check

A fresh Sonnet agent, given only the three input tables and the definitions, recomputed the tercile and threshold
frequencies, jet SDs, event counts and adjusted odds ratios. Everything matched within rounding except three small
differences from its own covariate conventions: Atlantic bomb odds per SD of jet 2.90 against 2.87 here, Pacific 2.70
against 2.67, Pacific HF-onset 1.50 against 1.48; the figures here sit inside its season-bootstrap intervals. Not
independently recomputed: the hinge and spline held-out gains and FDR q values, the knot intervals, the power simulation,
the interaction tests, the added-skill Brier scores, and the sensitivities.

## Analysis history

Pre-registration committed first (`7cd1aac`); fields extracted and predictors computed without any outcome (`8fedea0`);
power simulation committed with the verdict-wording amendment before the held-out run (`aa42b3c`). The code was exercised
in a mechanics mode (`JT_DEV=1`) whose "held-out" seasons were 2010-11..2014-15, inside the fit period. The held-out seasons
were scored in one full run (`results/looks.log`); no choice was changed after it. The extraction streamed 31.6 GB
(WeatherBench2 64 x 32 about 2 GB, ARCO-ERA5 about 29.5 GB), under the 50 GB gate.

## Files

    PREREGISTRATION.md      plan, definitions, decision rules, deviations
    extract_fields.py       Z500, U250, V250 at 00/12 UTC, 5.625 deg (resumable; work/ is ignored)
    features.py             trough_up, trough_loc, jet_up, trough_up_eddy per fix (outcome-free)
    jtlib.py                logistic fits, designs, season bootstrap
    power.py                simulation of the threshold rule
    analyse.py              the held-out analysis; figures.py, marginals.py
    results/features_2004.csv.gz  predictors per fix
    results/jet_trough-result.txt, thresholds.csv, verdicts.csv, fdr_table.csv, interaction.csv,
            added_skill.csv, grid_terciles.csv, marginals.csv, curves.csv, raw_deciles.csv, power.*, looks.log,
            curves_jet250.png, curves_trough_up.png

## Reproducing

    python3 extract_fields.py work/fields        # about 32 GB streamed
    python3 features.py ../intensity/results/fixes_2004.csv.gz work/fields results/features_2004.csv.gz
    python3 power.py 200 && python3 analyse.py --sens && python3 figures.py results && python3 marginals.py
