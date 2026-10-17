# Pre-registration: Gulf Stream effects on marine wind and waves

Written and committed before any buoy wind, temperature or wave value, and before any ERA5 pressure field,
was read for this question. Observations are NDBC buoy reports; ERA5 is used only for a pressure-gradient
control. This is **not** a hurricane-force-low study: it uses all wind states Oct-Apr, so the HF definition,
Decision 1 and the two ERA5 gust pipelines do not enter. No ERA5 gust is used.

## What has been seen

- Station coverage only: which NDBC buoys exist in the western Atlantic and which years they report (file
  names, no values). Chunk sizes of the ARCO-ERA5 pressure field (sizing only).
- Not seen: any buoy wind, air or sea temperature, wave height, or any pressure field for this question.
- The earlier Kuroshio test (null) is a different current and is not reused.

## Question (Jason, 2026-10-17)

"The Gulf Stream is a major factor in the western Atlantic. Can we study the effects and quantify them? We
typically bump winds up in instability over the Gulf Stream and increase wave heights, especially in wind-
against-current situations."

Plain-words answers to be given: **how many knots to add** for unstable air (water warmer than air) at a given
pressure-gradient wind, and **how many feet to add** for wind against the current, each as yes, no
(well-powered) or can't tell.

## Data

- **Observations.** NDBC historical standard meteorological files (hourly: wind speed `WSPD`, direction `WDIR`,
  `GST`, wave height `WVHT`, dominant period `DPD`, pressure `PRES`, air temperature `ATMP`, water temperature
  `WTMP`) for 2004-05 to 2024-25 winters, October to April (21 winters; the 2025-26 winter is not in the
  historical files). Fit seasons 2004-05 on, as everywhere in this repository. Hours 00 and 12 UTC only, so the
  samples are about 12 h apart.
- **Stations.** *Stream-influenced:* 41002 (34.0N 75.2W), 41001 (34.7N 72.5W), 41025 (35.0N 75.4W).
  *Shelf-water contrast (same weather, no strong current):* 41013 (33.4N 77.8W), 44014 (36.6N 74.8W), 41004
  (32.5N 79.1W), 44009 (38.5N 74.7W). Coastal fixed stations and nearshore buoys are not used (land effects).
  Stations are fixed now; none is added or dropped after seeing results. A station with under 5 winters of data
  is dropped and the drop is reported.
- **Wind speed** is the buoy's reported 8-minute mean, brought to 10 m with a neutral log law from the
  anemometer height in NDBC's station metadata (typically 5 m). Raw wind is the sensitivity.
- **Pressure-gradient wind** `Vg`: ERA5 mean sea level pressure (ARCO-ERA5, 0.25 degree), central differences
  over +-0.5 degree around each buoy, `Vg = |grad p| / (rho f)` with rho = 1.25 kg/m^3. It is a **control for
  the weather, not an outcome**. ERA5 assimilates buoy pressures and winds, so Vg is not independent of the
  buoys, but wind speed and temperature are not entered into it.
- **Instability** `dT = WTMP - ATMP` from the buoy (deg C). Clipped to [-3, 12].

## Pull plan

NDBC: 7 stations x 21 years x about 0.3 MB, under 0.1 GB. ERA5 MSLP: one global 2.2 MB field per time (measured
2.18 MB), 00 and 12 UTC, 1 Oct to 30 Apr, 2004-05 to 2024-25 = 8,988 times = about 19.6 GB; below the 50 GB
gate. Stop and ask if the measured total passes 25 GB. Existing caches hold gust only, so none covers this.

## Wind: how much to add for instability

**Model.** For each station-time with `Vg >= 15 kt` and observed 10 m wind present:
`log(V10) = a_station + b * log(Vg) + g(dT) + e`, where `g` is five bins of `dT` (stable < 0; 0-3; 3-6; 6-9;
>= 9 deg C) with 0-3 as the reference. The effect is reported in knots at `Vg` = 30 and 45 kt and as per cent.

**Primary tests (family W, 5 tests).**
- W1: pooled, mean of the 6-9 and >= 9 bins minus the 0-3 bin (log scale), all 7 stations.
- W2: same contrast, Stream-influenced stations only.
- W3: same contrast, shelf-water stations only.
- W4: W2 minus W3 (does the same air-sea temperature difference add more wind over the Stream than over shelf
  water?).
- W5: the >= 9 bin minus the 6-9 bin, pooled (does the effect keep growing, or level off?).

**Secondary (not in the family, labelled as such).**
- S1: effect separately for flow from 270-360 deg (cold air blowing off the continent over warm water) and
  other directions.
- S2: Dec-Mar only (drops October and April tropical and transition cases; hurricanes in October are otherwise in).
- S3: effect on buoy gust (`GST`), reported beside mean wind to show how much of a gust bump is mean wind.
- S4: unadjusted raw wind instead of 10 m-adjusted.

**Inference.** Effective n is winters (21) times about two independent weather states per week, not hours.
Resample **winters** with replacement (10,000 draws, 95% interval, percentile). Hours within a winter are never
resampled as independent.

## Waves: how many feet to add for wind against the current

The ERA5 wave model, as far as we know, does not include ocean currents, so ERA5 is not used for waves.

**Current state.** No current analysis is reachable from the work container. A climatological current
direction is used, fixed now: current flows toward 045 deg true at 41002, 065 at 41001, 030 at 41025. A buoy is
counted as **in the Stream** when `WTMP >= 22 deg C`; otherwise out. This misclassifies meanders and
warm-core-ring water, so the Stream effect is probably diluted (bias toward null). No shelf station has a
current direction; the shelf stations get the same direction label as 41002 for the sectors below, as a placebo.

**Sectors** (wind direction `WDIR`, the direction the wind blows from): *opposed* if `WDIR` is within 45 deg of
the direction the current flows toward; *following* if `WDIR` is within 45 deg of the direction the current comes
from; the remaining sectors are *cross* and are not used in the primary tests.

**Sample.** Steady wind only: `V10 >= 20 kt`, and `V10` and `WDIR` at the previous three hourly reports within
30% and 30 deg of the current report, so the sea has had time to grow; `WVHT` present; `DPD` present. Wind bins
(kt): 20-30, 30-40, >= 40.

**Model.** `WVHT = a_station + c(Vbin) + d * opposed + e` on Stream-state hours at the Stream stations, with
`DPD` as a swell-contamination covariate. Reported in feet.

**Primary tests (family V, 3 tests).**
- V1: opposed minus following wave height (ft) at Stream-state hours, Stream stations.
- V2: difference-in-differences: V1 minus the same opposed-minus-following difference at the same stations when
  `WTMP < 22` (not in the Stream).
- V3: placebo: opposed-minus-following at the shelf stations, using 41002's direction labels. Expected near zero.

**Inference.** Winter-block bootstrap as above. Because the Stream-state hours at three buoys may be few,
**report the number of qualifying winters and hours first**, and the minimum detectable difference
(see Power). If fewer than 8 winters or 150 hours qualify for V1, V1 is reported as "can't tell" without a
point estimate and no threshold is lowered afterwards.

## Multiplicity

Families W (5) and V (3) are each corrected with Benjamini-Hochberg; q is also given across all 8. Secondary tests
are listed with raw p and labelled secondary.

## Decision rules

- **Yes:** q < 0.05 within the family, and the sign is the forecaster's expectation (more wind for larger `dT`;
  higher waves when opposed).
- **No (well-powered):** only if the 95% interval lies within +-0.5 kt (wind, at Vg 30 kt) or +-0.5 ft (waves),
  which the power section must show was detectable.
- **Otherwise: can't tell**, with the interval reported.

## Power

Before any coefficient is read: a planted-effect simulation. Take the actual station-time pattern, plant `+3%`,
`+6%` and `+10%` of wind in the `dT >= 6` bins and `+0.5`, `+1`, `+2 ft` for opposed wave states in noise drawn
from the buoys' own residual spread (residuals taken after fitting without the effect terms, not the effect
estimates), resampled by winter; report the share of draws passing the rule above. If the simulation ignores
within-winter weather-regime dependence, say so and call the result an upper bound.

## What this cannot show

- Whether a given model (GFS, GFE blend) already includes the effect. This measures the real atmosphere against
  pressure-gradient wind. A forecaster's bump depends on what the guidance already shows, which is a different
  test (see the GFE blend project).
- Gulf Stream position meanders and ring water from a buoy temperature alone.
- Three buoys sit in the strongest current, so the wave result can be coarse. Wave-current coupling also depends
  on swell, which the steady-wind filter and `DPD` only partly remove.
- Stability from `WTMP - ATMP` at a buoy only; no marine-boundary-layer profile.
- Buoy mean wind is a sustained measure. The ERA5 gust formula inflates gusts in unstable air, so ERA5 gust is
  not used here, and buoy `GST` is secondary S3 only.

## Deviations (post hoc)

Logged in order. All tagged *before any outcome* (no buoy coefficient had been read; only counts and SEs).

1. *Before any outcome.* Wind sample also requires `V10 >= 1 kt`. A zero-speed report in a 15-kt gradient is a
   sensor or gap value, and the model uses log wind.
2. *Before any outcome.* Planted-effect residuals are taken after the **full** fit (with the effect terms), not
   after the null fit as written above. Residuals from the null fit would contain any real effect and the
   planted effect would sit on top of it. Only SEs and detection rates were read.
3. *Before any outcome.* Coverage seen when the dataset was built: 44,838 station-times; 25,160 pass the wind
   filter; 587 opposed-or-following hours in Stream water at the three Stream buoys (12+ winters), so V1 clears
   the pre-set floor of 8 winters and 150 hours.
