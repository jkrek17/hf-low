# More cyclones, or a larger share reaching hurricane force? (pre-registration)

Written 2026-10-08, before any model in this directory was fitted. The only data
looked at so far are the column headers and first rows of
`research/era5/hf_history/results/all_tracks.csv.gz` (commit `8b109b5`, hf-low PR 12)
and of the CPC index files. Anything below that changes after the first fit is
logged in `README.md` under "Departures from this plan", with the reason.

Question Q2 in the project's teleconnection-intensity list. Everything here is an
**ERA5 proxy** (pipeline A, `research/era5/hf_history`: 800 km ocean gust index,
HF-equivalent at 71.7 kt), not observations.

## Question

The warning archive (2004-05 to 2025-26, Oct-Apr) shows Atlantic HF-low frequency
rising 12.8% per SD of NAO and Pacific frequency rising 21.7% per SD of PNA
(`/mnt/project-files/teleconnection-test/REPORT.md`). Pipeline A's events show the
same main effects (RR 1.16-1.24 per SD). Split each effect into

    RR(HF lows) = RR(all cyclones) x RR(share of cyclones reaching HF)

exactly, on the log scale, and say which factor carries it.

## Hypotheses, stated before results

- **H-Atl.** NAO+ acts mostly through the share. A positive NAO strengthens and
  shifts the Atlantic jet poleward and raises baroclinicity along the storm track,
  which should make more of the cyclones that form become intense, while the
  total number of lows in a 30-67N domain changes less (more in the north, fewer
  to the south, both inside the box). Prediction: the share factor carries more
  than half of log RR(HF).
- **H-Pac.** PNA+ acts mostly through the share, for the same reason: a deeper
  Aleutian low and a stronger, eastward-extended jet. Same prediction.
- Against both: the Q1 thread's early read is that NAO and PNA shift where storms
  peak and their apparent depth, but not their deepening rate. If intensification
  is unchanged, more of the effect could be in how many cyclones form or how
  long they stay in the domain. These hypotheses could therefore fail, and a
  count-dominated or mixed result will be reported in full.

Decision rule, fixed now: let f = log RR(share) / log RR(HF), with a 95%
season-block bootstrap interval.
- "Mostly share" if the whole interval is above 0.5.
- "Mostly count" if it is below 0.5.
- Otherwise "both contribute / not resolved".
If RR(HF) itself is not distinguishable from 1 under the lagged index, f is not
interpreted and the result is reported as such.

## Data

- Cyclones: every pipeline A track in `all_tracks.csv.gz` (lows below 1010 hPa,
  linked 6-hourly, at least 24 h long, at least two fixes in a basin domain).
  Basin is the track's majority in-domain basin. HF = `gust800_kt >= 71.7`.
- Genesis: the track's first fix (`start`). Each track is placed on its genesis
  day; season and window are taken from the genesis date.
- Indices: CPC daily NAO and PNA (`norm.daily.{nao,pna}.index.b500101.current.ascii`,
  as used by the additive test and by Q1), and monthly ONI (exploratory only).

## Design

- **Lag (primary).** The index is the mean of days -10 to -4 before the genesis
  day, the same window Q1 uses, so the storm cannot feed its own index.
  Same-time index (days -3 to +3 around genesis) is reported as secondary, only to
  connect with the frequency effect as published.
- **Count model.** Daily genesis counts per basin, Poisson GLM, log link, with the
  lagged index, calendar-month fixed effects and a linear season trend (as in the
  additive test). Fitted separately to (i) all cyclones and (ii) HF cyclones,
  with the same design matrix. Then log RR(share) = b(HF) - b(all) exactly.
- **Share model.** Track-level logistic regression of HF on the same lagged
  index, month fixed effects and season trend, giving the HF-share odds ratio.
  The exact decomposition uses the risk ratio above; the odds ratio is reported
  alongside because it was asked for. They differ little when the share is small.
- **Standardisation.** Each index is standardised over the analysis days of the
  window (mean 0, SD 1 of the lagged 7-day mean), so effects are per SD.
- **Primary scope.** Atlantic with NAO and Pacific with PNA, one index per model.
  Window Oct 1 to Apr 28, seasons 2004-05 to 2025-26 (22 seasons), matching the
  published frequency effect and the rule that gust-based quantities before 2004
  are not used to fit or test.
- **Inference.**
  - Season-clustered (sandwich) SEs.
  - 2,000 season-block permutations: each season's index series is assigned to a
    different season's counts, keeping day-of-season. This gives p for each factor
    and for the share term.
  - A 2,000-draw season-block bootstrap, with seasons resampled jointly for all
    three quantities, gives intervals for RR(all), RR(share) and f.
  - Leave-one-season-out ranges for every primary number.
  - Effective n = seasons x independent index values per season: about 22 x 5.3
    = 116 (NAO) and 22 x 5.7 = 126 (PNA). Not the cyclone count.
- **Multiplicity.** Two primary decompositions (Atl-NAO, Pac-PNA). Each factor's
  p-value is reported raw and Bonferroni x2. Everything else is exploratory.

## Secondary and exploratory analyses (fixed now, labelled so in results)

1. Same-time index (days -3..+3), primary scope.
2. June-May window, 2004-05 to 2025-26.
3. Three-stage chain, primary scope: all cyclones -> cyclones reaching 1000 hPa
   or deeper (`minp <= 1000`, and at least 48 h long: the catalog's null-case
   population) -> HF. Shows whether the share effect enters at "becomes a real
   cyclone" or at "a real cyclone becomes HF".
4. Depth version, 1979-80 to 2025-26. HF replaced by a depth-equivalent: per
   basin, the `minp` cut that gives the same number of cyclones as the 71.7 kt
   gust cut over 2004-05 to 2025-26 Oct-Apr. Pressure depth did not drift before
   2001 (drift_check), so all 47 seasons are used; also reported for 1979-2000
   alone as an out-of-sample check. Caveat: NAO and PNA are defined by pressure
   at the same centres, so a depth outcome is more exposed to circularity than
   the gust one; the lag is the guard.
5. The cross-basin index (PNA in the Atlantic, NAO in the Pacific) and ONI, one
   at a time, primary window. ONI has about one independent value per season
   (effective n about 22), so a null is expected.
6. NAO and PNA fitted together, primary window.

## Not in scope

New ERA5 pulls; any change to pipeline A's definitions; intensity given HF
(that is Q1); interactions between indices (the additive test found none).
