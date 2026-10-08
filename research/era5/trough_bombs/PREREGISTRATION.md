# Pre-registration: upstream trough depth and bombs, 1979-2003 seasons (RA-24)

Written and committed before any ERA5 field was read for this question and before any 1979-2003 track,
fix or outcome existed. ERA5 proxy, **pipeline A** (`research/era5/hf_history`) re-run on pressure only for
1979-06-01 to 2004-05-31. "Bomb" is a pressure-defined outcome (the intensity framework's rapid-deepening
class), not the archive and not a gust quantity, so Decision 1 does not bar the pre-2004 seasons.

## What has been seen

- Exact chunk sizes of the pull (sizing only, no field values).
- PR 78's reported results (`jet_trough/results/jet_trough-result.txt`): held-out 2015-26 adjusted trough
  odds ratios for BOMB, Atlantic 1.231, Pacific 1.238 per SD, none passing FDR; 21-41% power at 1.15; and a
  zonal-eddy sensitivity that could not be reconciled with the primary. The 2004-2014 fit seasons were used by
  PR 78 to fit models and the trough climatology; their trough odds ratios were not reported there.
- Not seen: any Z500 or pressure field for 1979-2003, any 1979-2003 track, bomb frequency or trough value.

## Question (Jason, via RA-24)

Is there a trough-depth effect on bombs once the test has power, on seasons not yet used?
Plain-words answer to be given: yes, no (well-powered), or can't tell, and by how much.

## Data and design

**Fix population.** Pipeline A's rules run on 0.25 degree ERA5 surface pressure only (ARCO-ERA5, 6-hourly,
1979-06-01 00 UTC to 2004-05-31 18 UTC): smoothed MSLP minima below 1010 hPa, 20-75 N, 400 km de-duplication
(`hf_history/extract.py` logic with the gust step removed), linked with `hf_history/track.py` (<= 900 km per
6 h, tracks of >= 24 h), then `intensity/fixes.py` for the 00/12 UTC in-domain fixes (Atlantic 30-67 N,
98 W-10 E; Pacific 27-67 N, 135 E-120 W). Gust is not used. Transitioning tropical cyclones are in.

**Tracker check (outcome-free, before any use).** (a) Re-extract Jan 2008, Dec 2012, Jan 2020 and compare with
the committed 2004+ fixes; (b) compare the 1979-2003 catalog-track points (`era5_hf_catalog_tracks.csv`,
event and null tracks) with the re-tracked lows. Required: at least 99% of committed fixes recovered with the
same time, position within 0.25 degree and pressure within 0.2 hPa; every shortfall explained in the README.
If it fails, stop and report; no bomb statistic is computed.

**Outcome BOMB.** `ndr24` >= +1.0 Bergeron (24 h pressure fall scaled by sin(lat)/sin 60 deg) on the framework's
class-eligible fixes (track survives 24 h, |ndr24| <= 3). Same as PR 78. Pressure only.

**Seasons.** *New*: 1979-80 to 2003-04 (25 seasons). *Reference* (baseline and standardisation only): 2004-05 to
2014-15, PR 78's fit seasons. 2015-16 to 2025-26 are **not scored here**: no look at that block is spent.

**Fields.** WeatherBench2 64 x 32 (5.625 degree) 6-hourly ERA5, Z500 and U250, 00 and 12 UTC, June 1979 to
May 2015 (PR 78's grid; the brief's "1.5 degree" is not what PR 78 used for the trough). Reproduction check:
`trough_up` and `trough_up_eddy` recomputed for the 2004-2014 fixes must equal `jet_trough/results/features_2004.csv.gz`
to 0.01 m.

## The two trough definitions (both primary, equal status, fixed now)

Sector: 5.625 degree cells 500-2500 km from the fix at bearings 225-315 degrees (upstream in westerlies);
cells south of 14 N ignored; a fix with no cell in the sector is dropped (count reported). Trough depth is minus the
lowest value in the sector, so larger is deeper.

- **D1, climatological anomaly.** Z500 minus the calendar-month mean at that cell over 00 and 12 UTC of the
  reference seasons 2004-05 to 2014-15 (PR 78's primary, unchanged, applied to 1979-2003 too).
- **D2, zonal eddy.** Z500 minus its zonal mean at the same latitude and time (PR 78's sensitivity 4). It needs
  no climatology, so it carries no Z500 trend offset.

Neither is chosen after outcomes. Their results are reported side by side, and the verdict rule below uses both.

**Standardisation.** Per basin, the trough variable is clipped at the 1st and 99th percentile and standardised
with mean and SD of the reference-season eligible fixes of that basin (so an odds ratio is per reference-season
SD, as in PR 78).

## Model

Unpenalised logistic regression on the new-season eligible fixes, both basins pooled:

`BOMB ~ msl + dp12 + young + lat + cos(doy) + sin(doy) + pac + yr + z(trough)`

Covariates standardised with a `Prep` fitted on the reference-season pooled fixes; `pac` is a basin indicator;
`yr` is the season start year centred at 1991 in decades (absorbs slow ERA5 pressure-record changes, which are
not an effect to attribute to the trough). No gust, no jet, no interactions in the primary.

**Inference.** The odds ratio per SD with a 95% interval from a Wald test with season-clustered (CR1) variance
and t reference on 24 df. Direction is pre-specified positive; the p-value is two-sided. The 25 seasons are the
clusters, because fixes within a season share synoptic regimes. A season bootstrap (2,000 resamples) is
reported beside it and is not used in the rule.

## Tests and multiplicity

- **Primary family (2 tests):** D1 pooled OR, D2 pooled OR, new seasons.
- **Secondary family (6 tests):** the same by basin (Atlantic, Pacific; 4 tests), and pooled in 1979-80..2000-01
  only (2 tests; the stretch where ERA5's gust-index drift is documented and open-ocean pressure is least
  constrained).
- BH-FDR within each family and across all 8; every test is reported, passing or not.

**Sensitivities (reported, not in a family):** (1) the trough effect with jet_up (upstream max U250, 5.625 degree)
added; (2) without the `yr` term; (3) basin x covariate interactions; (4) Oct-Apr fixes only; (5) lead 24 (outcome
window (t+24, t+48 h]); (6) 2004-2014 alone and 1979-2014 pooled with the same estimator (the 2004-2014 seasons
are not fresh); (7) a nonlinear check, spline (4 df) vs linear on the new seasons, season-bootstrap; (8) bomb
frequency by basin and decade, as a record-quality diagnostic that selects nothing.

## Decision rule (fixed now)

- **Yes:** D1 and D2 both have positive pooled OR with q < 0.05 in the primary family. Report the OR and interval.
- **No (well-powered):** both definitions have an upper 95% bound below 1.10, and the simulated power for a planted
  OR of 1.15 is at least 80%. Report the smallest effect that was detectable.
- **Definitions disagree:** one passes and the other does not, or signs differ. The answer is "can't tell", with
  both results shown; the thread does not pick the one that fits.
- **Can't tell:** otherwise. Report the interval and the minimum detectable effect.

The agenda's prediction was a pooled OR of 1.10-1.25 per SD with an interval excluding 1 under the primary
definition (D1). It is compared with the result; it does not alter the rule.

## Power (run before the trough coefficient is read)

On the actual new-season design, outcomes are redrawn from the covariate-only fitted model plus a planted
trough effect (D1 and D2 separately) at OR 1.05, 1.10, 1.15, 1.20, 1.25 per SD, 400 draws each. A draw "detects"
when the primary inference gives a positive OR with p < 0.05. Minimum detectable effect = smallest planted OR with
at least 80% power (linear interpolation). Limitation: draws are independent given covariates, so they omit the
season-level clustering that the real inference pays for; the simulated power is an upper bound, and the real
interval width is reported next to it.

## Looks

Scores 1979-80..2003-04 once (a look at the pre-2001 block for 1979-80..2000-01 and the 2001-2003 seasons).
Written to `research/era5/looks/ra24-trough-bombs.log` with a UTC stamp, not to the shared log. 2015-26 not scored.

## Deviations (post hoc)

None yet.
