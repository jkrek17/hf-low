# Do combined teleconnection states move where HF lows peak, or how many peak in a region? (pre-registration)

Written 2026-10-08, before any interaction or main effect in this directory was fitted. ERA5 **proxy**,
**pipeline A** (`research/era5/hf_history`: 800 km ocean gust index, HF-equivalent at 71.7 kt). Not observations.
Catalog and rationale: `CATALOG.md` (commit `643deac`).

## What was looked at before this was written

- Column headers and marginal counts of `all_tracks.csv.gz` (Oct-Apr counts, mean and SD of peak latitude and
  longitude by basin; table in `CATALOG.md`). No count or position was related to any index.
- Predictor-to-predictor correlations over Oct 2004-Apr 2026 (design check): ONI-MJOWP 0.02, ONI-MJODL -0.01,
  NAO-PNA 0.02, NAO-ONI 0.01, PNA-ONI 0.07, PNA-MJOWP -0.14, MJOWP-MJODL 0.33. The raw MJO series correlate
  0.4-0.5 with ONI monthly (`teleconnections.json`); the causal 90-day high-pass used here removes that.
- What was already known from earlier threads, which could bias choices: the basin-count additive test (no
  interaction in totals; ONI x MJO archive hit gone in ERA5), Q1 (NAO shifts Atlantic storm position about
  +1.5 deg lat and +1.9 deg lon per SD), Q2 (`freq_split`). The sibling threads' plans
  (`research/era5/enso_pna/PREREGISTRATION.md` on its own branch) were read for conventions only.

## Questions

For two teleconnection states A and B, with both main effects in the model:

- **L (location).** Does the A x B interaction shift the position (longitude, latitude) of the 800 km gust-index
  peak of HF-equivalent lows in the basin? *Constructive* means the joint-state shift is larger than the sum of
  the two separate shifts; *destructive* means smaller or opposite.
- **C (count).** Does the A x B interaction change the number of HF-equivalent lows in the basin, beyond the
  product of the separate effects?

## Primary tests (ten: five pairs x {L, C})

| Test | Basin | A | B |
|---|---|---|---|
| T1 | Pacific | ONI | MJOWP (120E+140E, enhanced convection positive) |
| T2 | Pacific | PNA | MJODL (160E+120W, enhanced convection positive) |
| T3 | Atlantic | NAO | PNA |
| T4 | Atlantic | NAO | ONI |
| T5 | Atlantic | NAO | SPV (10 hPa, 60N zonal-mean zonal wind anomaly) |

**T5 is contingent.** It runs only if the derived SPV series passes its own check (below). If it does not, T5L and
T5C are reported as "not run" and the family has eight tests. Nothing about T5 depends on outcomes.

Hypotheses (two-sided; no sign is claimed for any interaction):
- H0 for every test: gamma = 0 (the A x B coefficient).
- H1 for T1: jet-extension physics predicts a non-zero longitude and latitude interaction; T2: wave-train
  physics predicts a shift of P-E (Gulf of Alaska) peaks; T3: a +PNA ridge changes Atlantic cyclogenesis
  longitude; T4: ENSO changes the NAO-related track; T5: polar vortex state changes the NAO-related poleward shift.
  These are the mechanisms in `CATALOG.md`; they motivate the test and set no sign.
- Against: the basin-count additive test found nothing, the interaction is a second-order effect, and 22
  seasons is a small sample. A null or inconclusive result is likely and will be reported as such.

## Data and variables (fixed)

- Tracks: `hf_history/results/all_tracks.csv.gz`. Each track sits on its genesis day (`start`, first fix).
  HF-equivalent: `gust800_kt >= 71.7`. Outcome position = `peak_lon`, `peak_lat` (where the gust index peaks;
  **not genesis position**, which is not committed before 2004). Atlantic longitude unwrapped (+360 if < 180).
- Window and seasons (primary): Oct 1 + 210 days; seasons 2004-05..2025-26 (22 seasons). Gust-based quantities
  fitted and tested only from 2004-05 (decision 1).
- Predictors, all standardised over the analysis days:
  - NAO, PNA: mean of the CPC daily index over days -10..-4 before the genesis day.
  - ONI: monthly ONI of the calendar month containing day -7.
  - MJOWP, MJODL: the CPC pentad longitude indices expanded to daily by nearest pentad centre, sign flipped so
    positive = enhanced convection, each longitude scaled to unit SD over 1979-2025, averaged as above, then
    high-passed *causally* (value minus the trailing 90-day mean), then the value of day -7.
  - SPV (T5): daily 00 UTC zonal-mean zonal wind at 10 hPa and 60N, ERA5 (WeatherBench2 1.5 deg to 2023-01-09,
    ARCO 0.25 deg after; ERA5 is a proxy), turned into an anomaly by removing the day-of-year mean and SD over
    1979-2025, then the mean of days -30..-11.
- **SPV check (rule).** Pass only if (i) on the one-year overlap of the two sources (2022-01-01..2023-01-09) the
  daily correlation is >= 0.98, and (ii) the 10 hPa 60N wind is < 0 on at least four of the five known winter
  reversals within +/-5 days of 2009-01-24, 2013-01-06, 2018-02-12, 2019-01-02, 2021-01-05. A failed check
  withdraws T5; the SPV definition does not change.

## Models (fixed)

Month fixed effects and a linear season trend are in every model. y is always for one basin.

- **L**: OLS, y = (lon, lat) of HF-equivalent peaks, on A + B + A.B (+ month FE + trend); one row per event.
  The test of gamma is **joint over (lon, lat), 2 df**. The two coordinates are reported separately.
- **C**: Poisson (log link), daily HF-equivalent genesis counts in the basin, same right-hand side.
- Inference for gamma: a sign-flip score test (restricted wild cluster bootstrap with Rademacher weights,
  seasons as clusters, 50,000 draws, efficient score with the other terms projected out). This is the
  *permutation* p. Intervals: season-pairs bootstrap, 2,000 draws, 95% percentile. CR1 clustered SE reported
  alongside. Randomisation seed 20261008.
- Descriptive four-corner table from the fitted model: predicted outcome at (A, B) = (+1,+1), (+1,-1), (-1,+1),
  (-1,-1) SD, each beside the **additive** prediction (main effects only). Reported for every test; untested.

## Decision rules (fixed now)

Smallest effects of interest (SESOI): **L** 0.10 of the outcome's SD per SD x SD for each coordinate (this is
roughly half the Q1 NAO main effect on position; the SD is the pooled SD of the analysis events); **C** RR 1.05 per SD x SD.

- **Detected:** BH q < 0.05 over the family of primary tests AND the bootstrap interval of gamma excludes 0
  (for L: the joint test, and at least one coordinate interval excludes 0).
- **Well-powered null:** the whole 95% bootstrap interval of gamma lies inside +/-SESOI (L: both coordinates).
- **Inconclusive:** anything else. Said in those words. A test whose detectable effect (2.8 x SE) exceeds the
  SESOI cannot return a well-powered null.
- "Constructive" or "destructive" is only claimed for a **detected** interaction, and only in the sense defined
  above, read from the four-corner table.
- Power is the detectable effect at 80% power from `power.py` (season-block permutation spread of gamma, 2.8 x SD),
  reported with every result. Effective n = seasons x independent predictor values per season.

## Multiplicity

- Primary family: the ten tests (eight if T5 is withdrawn); BH q on the permutation p (L: joint p).
- Whole-directory family: every gamma fitted in this directory, primary and secondary, with BH q. Stated as such.
- Sibling threads ("El Niño and PNA together"; "El Niño flavor and Kuroshio genesis") are members of the family.
  If their primary p-values are published before this closes, a combined BH is added beside both.
  Neither is re-run here.

## Secondary and exploratory analyses (fixed now; labelled so)

1. **S1 same-time predictors** (NAO/PNA days -3..+3; ONI of genesis month; MJO of day 0; SPV days -10..+10) for all
   primary tests. Circularity is the exposure; the contrast with the lagged result is the point.
2. **S2 box counts.** Poisson interaction for HF-equivalent peaks in each of the three boxes of the basin
   (`BOXES` in `core.py`), days as above. Boxes fixed in `CATALOG.md`.
3. **S3 longer record, HF-equivalent peaks, location and count**: seasons 1979-80..2025-26 with an era term
   (season <= 2000), and 2001-02..2025-26 on its own. Within-era variation only (decision 1).
4. **S4 all cyclones** (lows below 1010 hPa, >= 24 h): storm-track position (peak location) and daily count, same two
   record choices as S3. The peak position is gust-defined, so the era term and the 2001+ run apply here too.
   The count of lows is pressure-based and needs no era term but the same layout is kept.
5. **S5 HF-equivalent minus storm track:** does the HF-equivalent position shift differ from the all-cyclone shift?
   Difference of the gamma estimates from S4 (all cyclones) and primary (HF), 2004-05 on, season-pairs bootstrap.
6. **S6 Jun-May window**, primary seasons. Not run for the Pacific MJO tests where the MJO impact is winter-weighted;
   reported for completeness otherwise.

Nothing else is added after the first fit. Any later change goes in `README.md` under "Departures from this plan",
labelled post hoc, with both versions reported.

## Not in scope

New pulls beyond the SPV series (under 10 GB, well below the 50 GB gate); genesis position before 2004; QBO and
three-way terms (not powerable); ONI x PNA and ENSO flavour (sibling threads); intensity given HF; awips-tools.
