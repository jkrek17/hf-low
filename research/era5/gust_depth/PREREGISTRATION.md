# Why do gust-based and depth-based Pacific HF counts disagree? (pre-registration)

Written 2026-10-08, before any gust-versus-depth skill, disagreement table or teleconnection rerun was computed.
Item 17 of "Future tasks" in `STATUS.md`. Everything is an **ERA5 proxy, pipeline A**
(`research/era5/hf_history`: 800 km ocean gust index; tracks in `results/all_tracks.csv.gz`). Not observations.

**What was looked at before this file was committed.** Column headers and first rows of `all_tracks.csv.gz`, the
catalog and the archive; archive event counts per season and class (archive only, no ERA5); and one pipeline check: how
well the track-level matching rule below recovers the 1,410 archive-to-track pairs already in the catalog
(`archive_events` column). That check was run for several (distance, penalty) settings and the setting below was fixed
after seeing it (recovery 0.88 of the catalog pairs at 800 km / 25 km per 6 h; 0.85 at 500; 0.90 at 1200). It could only
use gust-matched pairs, so it cannot tell us whether the rule is neutral between the definitions: a sensitivity run at 500
and 1200 km is planned for that reason. **Prior knowledge that shapes H3:** I know from the finished work that the
depth version of the Pacific PNA share is about +11% per SD against +4% for gust (PR 14), and that ONI's total effect is
1.066 per SD for depth against 1.02 for gust (PR 34). H3 is therefore a hypothesis about an already-seen disagreement,
not a blind prediction; H1 and H2 are blind to ERA5-versus-archive skill and to the disagreement geography.

## Question

In the ERA5 proxy there are two ways to say a cyclone is "hurricane force": **G**, the gust index (800 km ocean gust,
>= 71.7 kt) and **D**, depth (track minimum ERA5 MSLP at or below a cut that gives the same count). They are not the
same storms, and the Pacific frequency work gives different answers under each. (1) How often and where do they
disagree? (2) Against the archive, which is the better stand-in, and for what? (3) Does the disagreement explain the two
contested Pacific results? (4) Do other finished results depend on the choice?

## Definitions, all fixed here

- **Population.** Every row of `all_tracks.csv.gz` (lows below 1010 hPa, >= 24 h, >= 2 fixes in a basin domain), basin as
  given. Season = June-May labelled by the June year, from the track's first fix, as in `freq_split`.
- **G (published).** `gust800_kt >= 71.7`.
- **D (published-style).** `minp <= cut`, where the cut is the one for which the count of D tracks equals the count of G
  tracks, per basin, over Oct-Apr genesis in seasons 2004-2025 (the `freq_split` S4 construction; it was 965.0 hPa in
  the Pacific there).
- **D' (location-adjusted depth, diagnostic).** `minp - clim(cell) <= cut'`. `clim` = mean `minp` of all tracks of the
  same basin in the same calendar month of `peak_time`, 10 degrees of latitude band and 30 degrees of longitude band of
  the peak position, over seasons 2004-2025 (cells with fewer than 25 tracks fall back to basin by month). `cut'` is
  count-matched to G exactly as D.
- **Archive reference.** Events in `docs/data/hf-lows.json` of class `low` (1,890; tip-jet and centreless events are
  treated separately, below), seasons 2004-2025. Every archive event is a positive regardless of peak category, as in the
  calibration (sensitivity: peak category `HF` only).
- **Matching archive event to track** (`match.py`). Same basin; track [start, end] within 12 h of the archive [start, end];
  distance = great-circle km from the track's peak position to the archive position (interpolated) at the track's peak
  time, plus 25 km per 6 h by which the peak time falls outside the archive's fixes; the nearest candidate within 800 km
  is the match. The gust value and pressure do not enter. Archive events with no match are misses for both definitions.
  Sensitivity: 500 km and 1200 km.
- **Fit seasons** 2021-22 to 2025-26 (the pipeline A calibration seasons). **Test seasons** 2004-05 to 2020-21
  (17 seasons, held out of every threshold). Per-basin thresholds for the skill comparison, G* and D*, are the values giving
  forecast count = archive count of matched-eligible events in the fit seasons (June-May), per basin. (G* is therefore not
  exactly 71.7; the published G is also scored.)

## Hypotheses

### H1. Which definition tracks the archive (blind)

The archive is a wind-warning record, so the gust index should agree with it better than depth, in both basins.
Prediction: on the test seasons, AUC(G) - AUC(D) >= +0.02, and HSS(G*) > HSS(D*), in the Pacific and the Atlantic; depth's
loss is mostly extra false alarms (higher FAR) rather than extra misses. Wrong if the paired difference interval includes 0
or is negative.

Tests (each basin, paired season-block bootstrap, 2,000 draws, 17 seasons):
- T1a dAUC on all tracks (positives = tracks matched to archive `low` events; AUC by Mann-Whitney over all tracks);
- T1b dAUC on "real cyclones" (`minp <= 1000` and `n_fix >= 8`), the `freq_split` S3 population;
- T1c dHSS at the fit-season count-matched thresholds, contingency as in `hf_history/calib.py` (unmatched archive events
  are misses).
Also reported, untested: POD, FAR, CSI, bias for G, G*, D*; leave-one-season-out range; recall of archive tip-jet and
centreless events by G and by D (D cannot see a tip jet by construction).
Primary comparison: **Pacific T1a**.

### H2. Where and when they disagree (blind)

Gust at a fixed depth is not constant. Prediction: at fixed `minp`, the gust index falls with latitude (the same pressure
drop spread over a larger Coriolis parameter and a lower background pressure gives weaker wind), so depth over-flags
high-latitude lows and gust over-flags lower-latitude, colder-season ones. No sign is predicted for longitude sector or
season; they are reported two-sided.

Tests:
- T2a-d, per basin. OLS of `gust800_kt` on `minp` and `minp^2`, `peak_lat`, west-sector indicator (Pacific: longitude
  below 180 E; Atlantic: below 315 E), cold-season indicator (peak month Nov-Feb), `log(n_fix)`, season trend; real
  cyclones (`minp <= 1000`, `n_fix >= 8`) in test + fit seasons 2004-2025; season-clustered SEs. The coefficients on
  latitude (predicted negative), sector, cold season and `log(n_fix)` are the four tests.
- T2e-i, per basin. Among tracks flagged by exactly one definition at the published thresholds (G-only versus D-only),
  difference in mean `peak_lat`, share in the west sector, share in the cold season, median `n_fix`, and the share matched
  to an archive event (the archive's verdict on each kind of disagreement). Season-block bootstrap intervals and p-values.
Also reported, untested: counts of both / G-only / D-only / neither by season, their latitude and longitude histograms, and
the archive-match rate of each group.

### H3. The disagreement in the two contested results is carried by location and background pressure (informed, see above)

Q1 found that NAO and PNA shift where storms peak and the background pressure. A fixed pressure cut converts a shift in
background pressure into a change in count; a wind cut does not. Prediction, Pacific only: (i) the lagged-PNA share term
under D' (RR) is nearer to G's than D's is; (ii) the ONI total effect on the HF count under D' is nearer to G's than D's is.
"Nearer" = absolute difference in log RR from G is smaller for D' than for D. H3 is supported if (i) and (ii) both hold, partly
supported if one holds, and rejected otherwise. If rejected, the alternative reading is that the depth-only storms are a
different population (larger or slower lows) and the responses differ for that reason; the group decomposition below
then carries the answer.

Design: exactly the `freq_split` / `enso_pna` machinery (daily genesis counts, Poisson with calendar-month fixed effects and
a linear season trend, season-clustered SEs, season-block bootstrap, index = mean of days -10 to -4 before genesis,
ONI of the month of day -7, standardised over the analysis days), Oct-Apr, seasons 2004-05 to 2025-26. Outcomes: all cyclones, G, D,
D', and the archive's own events. The archive event is placed on the genesis day of its matched track; an unmatched archive
event is placed at its archive start date minus the median (start date minus track genesis) of the matched pairs.
- T3a-b, log RR(G), RR(D), RR(D'), RR(archive) per SD of PNA, and of ONI (8 estimates).
- T3c-h, paired differences in log RR (season-block bootstrap, 2,000 draws): D - G, D' - G, archive - G, archive - D,
  archive - D' for PNA and for ONI (10 differences), plus the share terms RR(share) = RR(outcome)/RR(all) for G, D, D' under PNA
  (3 estimates; the share of the archive is not defined).
- Group decomposition (descriptive): RR per SD of PNA and ONI for tracks that are both, G-only, D-only.
- A 1979-80 to 2025-26 version with an era term (before 2001-02 / after) for D and D' only, as secondary, reproducing the
  `enso_pna` 47-season depth contrast; G is not used across the 2001 break (decision 1).
Replication check: the G, D and all-cyclone numbers must reproduce `freq_split` and `enso_pna` (PNA share 1.043 for G; Pacific D
share 1.111; ONI total 1.018 for G, 1.066 for D over 47 seasons) to the third decimal, or every difference is explained.

### H4. Other finished results (read-through, no new fit)

For every finished result in `STATUS.md` (PRs 11, 12, 14, 16, 17, 18, 27, 29, 34, 38, 40, 41) state whether it used G, D or
both and whether a switch would plausibly change its conclusion; list those where both were reported and agree, disagree, or
that never had a depth version. A fresh agent checks the table against the READMEs.

## Multiplicity

All tests T1a-c, T2a-i, T3 estimates and differences (those marked as tests above) go into one Benjamini-Hochberg family
with q computed over the whole set and also within H1, H2 and H3. Atlantic rows are in the same family. Nothing outside the
list above is a test. Anything added after this file is committed is post hoc, labelled so in the README, and reported
next to the pre-registered version.

## Decision rules for "which should we use"

Written now so the verdict cannot be tuned to the results.
- G is "more trustworthy against the archive" where T1a/T1c favour it with the paired interval excluding 0 for the
  Pacific; D is where it favours D; "no preference" if the interval includes 0 in both T1a and T1c.
- Where the two responses to a teleconnection differ (T3 differences with q < 0.10, or intervals excluding 0), the
  archive's own response arbitrates: the definition whose RR is nearer to the archive's, with the archive-minus-definition
  interval excluding 0 for the other, is preferred for that question. If the archive interval is too wide to separate them, the
  answer is "can't tell", with the width stated.
- A combined or new definition (e.g. a logistic of the archive on gust, depth and latitude) is built only if the held-out
  AUC of the combination exceeds the better of G and D by 0.01 or more in the Pacific; it is then exploratory, reported
  alongside, and not used for any claim in this study.

## Power and effective n

For the skill comparisons the unit is the season (17 test seasons, about 36-40 Pacific archive events each). For the
teleconnection comparisons the effective sample is 22 seasons x about 5.7 independent PNA values per season, and about
one ONI value per season; `freq_split` found that a Pacific effect of about 1.09-1.10 per SD is the smallest detectable at
80% power, so differences between definitions below that are expected to be inconclusive and will be reported as such.

## Out of scope

New ERA5 pulls (none needed); changing the published thresholds; pre-2001 gust levels; the Atlantic NAO
(Atlantic rows are reported as a check, not a goal).
