# Near-storm intensity framework: Hart phase space plus environment (ERA5 proxy, pipeline A tracks)

For each ERA5 low at 00 and 12 UTC, this framework takes Hart's (2003)
cyclone phase space position and the near-storm environment. From those it
gives two sets of probabilities:

- **Intensity class over the next 24 h:** rapid decay, decay, steady,
  deepening or rapid deepening, on the normalised deepening rate in
  Bergerons. One Bergeron is 24 hPa in 24 h scaled by sin(latitude)/sin(60°).
- **Hurricane-force-equivalent gust:** the chance that the low reaches
  pipeline A's gust index (71.7 kt) within 24 h and within 48 h.

Research code: nothing here feeds `docs/` or the site.

**This is a proxy.** "HF" here means pipeline A's gust index, a cyclone whose
ERA5 gust field looks like the ones OPC warned for as hurricane force in
2021-26. It is not the archive. The intensity classes come from ERA5 central
pressure on pipeline A's tracks, and the archive does not record them.

## Population and labels

- **Fixes.** Every fix from pipeline A's tracker (`../hf_history/extract.py`,
  `track.py`) at 00 or 12 UTC, inside pipeline A's basin domains, in seasons
  2004-05 to 2025-26. That is 159,430 fixes on 36,598 tracks. Every low below
  1010 hPa is included, not only the catalog's events and null cases, so the
  base rates are those of all lows.
- **Re-extraction.** The full pipeline A re-run reproduces the catalog's
  event counts exactly: 2,254 Atlantic and 1,903 Pacific for 1979-2025.
- **Classes.** The cut points are -1, -0.3, +0.3 and +1 Bergeron.
  - Frequencies: rapid decay 0.7%, decay 18.6%, steady 60.7%, deepening
    16.5%, rapid deepening 3.4%.
  - A fix gets no class when its track ends within 24 h (59,400 fixes). Five
    fixes with |NDR| > 3 are dropped as tracker relinks.
- **HF.** The gust index reaches 71.7 kt at any fix in (t, t+24 h] or
  (t, t+48 h]. A track that ends inside the window counts as not reaching it.
  - Base rates: 3.4% within 24 h, 4.3% within 48 h. 1.8% of fixes are
    already HF at t.

## Predictors

Everything is measured relative to the fix, from ERA5 on a 1.5° grid, the
resolution class Hart used (`env.py`).

| Group | Predictor | Radius |
|---|---|---|
| Storm state | central pressure, 12 h pressure change, track age, translation speed (previous 6 h), current gust index, latitude, basin, day of year | at the fix |
| Hart | B (900-600 hPa thickness asymmetry across motion) | 500 km |
| Hart | -V_T lower (900-600 hPa) and -V_T upper (600-300 hPa) | 500 km |
| Upper | 250 hPa maximum wind | 1000 km |
| Upper | 300 hPa divergence, 500 hPa absolute vorticity advection | 500 km mean |
| Baroclinic | 850-500 hPa Eady growth rate | 500 km mean |
| Surface | SST, maximum SST gradient, SST minus 500 hPa temperature, upward sensible plus latent heat flux | 500 km |
| Moisture | total column water vapour | 500 km mean |

Two departures from Hart's own setup:

- **Levels.** The thermal wind uses 900 (interpolated), 850, 700 and
  600 hPa, and 600, 500, 400 and 300 hPa. Hart used 50 hPa spacing, so the
  magnitudes are not FSU's. The signs and zero crossings carry over.
- **Motion.** The direction of motion for B comes from the previous 6 h, not
  a centred difference, so no future position is used.

`hart.py` is an unchanged copy of awips-tools' `cyclone_phase_space/cps/hart.py`.

**Sources.** Up to 2023-01-09 the fields come from WeatherBench2's
conservative 1.5° regrid of ERA5. After that, ARCO-ERA5 0.25° is coarsened
here onto the same cells. At an overlap time the two agree to 0.08-0.16 m in
height and 0.01-0.03 m/s in wind (RMS).

## Models and validation

Logistic regression with L2 penalty (C = 1) on standardised predictors.
Classes use a multinomial fit; HF is binary. There are four predictor sets:

- `state`: storm state only
- `hart`: state plus B, -V_T lower and -V_T upper
- `full`: state plus Hart plus environment
- `full_nogust`: full minus the current gust index

Each set is cross-validated leave-one-season-out over the 22 seasons. Skill
is measured against basin-by-month climatology (RPSS for the ordered classes,
BSS for HF), with 90% intervals from resampling whole seasons. Yes/no scores
use the probability cut at which the training seasons' forecast count equals
their observed count, applied unchanged to the held-out season.

## Results (`results/skill.txt`, seasons 2004-05 to 2025-26, n = 22 seasons)

| | RPSS / BSS | AUC | POD | FAR | CSI | HSS | bias | season HSS range |
|---|---|---|---|---|---|---|---|---|
| Classes, state | 0.242 | | | | | | | |
| Classes, hart | 0.290 | | | | | | | |
| Classes, full | 0.305 | | | | | | | |
| Rapid deepening, state | | 0.918 | 0.36 | 0.64 | 0.22 | 0.34 | 1.00 | 0.27-0.40 |
| Rapid deepening, hart | | 0.945 | 0.47 | 0.53 | 0.30 | 0.45 | 1.00 | 0.37-0.52 |
| Rapid deepening, full | | 0.950 | 0.49 | 0.51 | 0.32 | 0.47 | 1.00 | 0.41-0.53 |
| Rapid decay, full | | 0.948 | 0.18 | 0.82 | 0.10 | 0.18 | 1.01 | 0.05-0.30 |
| HF within 24 h, state | 0.370 | 0.963 | 0.56 | 0.44 | 0.39 | 0.55 | 1.00 | 0.49-0.61 |
| HF within 24 h, full | 0.423 | 0.972 | 0.60 | 0.40 | 0.43 | 0.59 | 1.00 | 0.53-0.65 |
| HF within 24 h, onset only, state | 0.221 | 0.950 | | | | | | |
| HF within 24 h, onset only, full | 0.291 | 0.963 | | | | | | |
| HF within 48 h, full | 0.374 | 0.957 | 0.57 | 0.43 | 0.40 | 0.55 | 1.00 | 0.48-0.61 |

"Onset only" means the fix is not HF at t (base rate 2.2%). Reliability
tables are in `skill.txt`. Rapid deepening and HF are reliable to within
about 0.03 across the range. Rapid decay over-forecasts above 0.2, where
there are few cases.

### What the phase space adds

- **Hart's terms carry most of what the environment adds for rapid
  deepening.**
  - HSS goes from 0.34 to 0.45 when B, -V_T lower and -V_T upper are added.
    The rest of the environment adds only 0.02 more.
  - Rapid deepening concentrates in the asymmetric, cold-core corner (large B
    with strongly negative -V_T lower and upper), the frontal-wave stage of
    Hart's diagram (`results/phase_rapid_deepening.png`).
- **The environment matters more for hurricane force than Hart does.**
  - For onset within 24 h, BSS rises from 0.221 with storm state to 0.246
    with Hart and 0.291 with everything.
  - The largest environmental terms are column water vapour and the 850-500
    hPa Eady growth rate (`results/coefficients.csv`).
- **Rapid decay is mostly a warm-core signal.**
  - It sits where -V_T lower is positive and B is small
    (`results/phase_rapid_decay.png`).
  - Tropical cyclones and warm seclusions dominate it, and the Bergeron's
    1/sin(latitude) factor inflates their filling rates at low latitude.
  - Its yes/no skill is weak (HSS 0.18).
- **The class model barely uses the gust index.** Dropping the gust predictor
  costs 0.003 RPSS. The HF model does lean on it (0.042 BSS).

Individual coefficients are partial effects among correlated predictors. SST
and water vapour, for example, enter with opposite signs. Read them for
which groups matter, not as physical signs.

## Checks

- Leave-one-season-out, with every season held out once.
- **Era drift.** Trained on 2004-2014 and tested on 2015-2025, the full model
  scores RPSS 0.309 and BSS (HF 24 h) 0.426. Leave-one-season-out skill by
  decade is 0.299, 0.302, 0.314 (classes) and 0.403, 0.421, 0.444 (HF 24 h),
  so skill does not fall over time.
- **Basin.** A Pacific term is in every model. Skill by basin is not broken
  out yet.
- **Transitioning tropical cyclones are in.** They lie in the warm-core half
  of the phase diagrams.

## Gates and limits

- **Fit and test window.** Everything is fitted and scored on seasons 2004-05
  onward only, per the project rule.
  - Pipeline A tracks and environment fields exist for 1979-2003 and the
    fitted model can be applied there.
  - Gust index values before 2001 drift upward at fixed depth (`c9dc994`,
    `c55e74d`). Any HF number before 2001 carries that caveat, and none is
    reported here.
- **Every outcome is ERA5's own.** The model forecasts what ERA5 does 24-48 h
  later from what ERA5 shows now. That is a perfect-prognosis setup: applied
  to a forecast model's analysis, its skill depends on how close that
  analysis is to ERA5. It has not been tested on GFS or any other operational
  input.
- **Classes condition on survival.** A track that ends within 24 h (filled,
  merged, or left the 20-75N band) gets no class. Rapid decay is therefore
  under-counted.
- **Skill is dominated by persistence.** The headline AUCs (0.95-0.97) are
  mostly the storm-state terms. The useful figure is the gain over `state`,
  which is reported with its interval.
- **The intervals are narrow because each season holds thousands of fixes.**
  They describe season-to-season sampling of the skill, not uncertainty in
  the predictor definitions. The per-season HSS range is the more honest
  spread.

## Files

    fixes.py   pipeline A track points -> forecast fixes, storm state, outcomes
    env.py     Hart terms and environment at each fix (streams ERA5; see below)
    model.py   fits, leave-one-season-out validation, model.json, predict()
    figures.py phase diagrams shaded by observed outcome
    hart.py    copy of awips-tools' Hart module
    results/fixes_2004.csv.gz  fixes and outcomes, seasons 2004-2025 (159,430)
    results/env_2004.csv.gz    their predictors
    results/skill.txt, coefficients.csv, model.json, phase_*.png

`model.json` holds the fitted coefficients and the scaling.
`model.predict(model["cls"], rows)` gives probabilities for new rows built
the way `model.load` builds them.

## Reproducing

    # pipeline A lows and tracks (about 370 GB streamed; needs go-ahead)
    ERA5_WORK=work python3 ../hf_history/extract.py 1979 2025
    ERA5_WORK=work python3 ../hf_history/track.py 197906 202605 work/track_points.csv
    python3 fixes.py work/track_points.csv work/fixes.csv
    # predictors (about 250 GB from WeatherBench2 plus 0.8 TB from ARCO for 2023-26; needs go-ahead)
    python3 env.py work/fixes.csv work/env 8
    # models and figures, from the committed tables (about 25 min)
    python3 model.py results/fixes_2004.csv.gz results/env_2004.csv.gz results
    python3 figures.py results/fixes_2004.csv.gz results/env_2004.csv.gz results
