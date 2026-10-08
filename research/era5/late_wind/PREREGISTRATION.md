# Why does the HF wind come after the lowest pressure in about 29% of storms? Pre-registration (RA-7)

Written and committed **before any predictor was joined to the late-onset outcome in this thread**.

What had been looked at: the published marginal shares of PR 16 (onset before / at the minimum-pressure fix: Atlantic 66% / 21%,
Pacific 77% / 16%, both 71% / 18%; so "at or after" is 29%), the column names, row counts and missing-value counts of
`hf_structure/results/fixes.csv`, `hf_history/results/lifecycle_events.csv` and `intensity/results/env_2004.csv.gz`, and the
**prevalence of each binary predictor** (listed below, computed without the outcome). Nothing else: no outcome has been tabulated
against any predictor, and no model has been fitted.

All results are an ERA5 **proxy**. Pipeline A = `research/era5/hf_history` (800 km ocean gust index, threshold 71.7 kt). The archive
is not used: "onset" and "minimum pressure" are ERA5 quantities on pipeline A tracks. Transitioning tropical cyclones are **out**
(the 114 TC-linked events of `lifecycle.py`). Atlantic fixes north of 60N are lower confidence (no terrain mask); they are
the subject of T2 and are removed in a sensitivity run.

## Question (agenda RA-7), and the answer I will give

In about 29% of storms the first HF fix comes at or after the fix of lowest pressure. What sets those storms apart: a gust maximum far
from the centre (barrier flow, post-occlusion wind), the high-latitude Atlantic near Greenland and Iceland, or a threshold artefact
or ordinary phase variant (marginal storms, whose crossing time is set by noise in the wind-pressure phase)?

Answer given in plain words: which of those is supported, which are not (well powered), and **by how much**: the late share inside and
outside each factor, the share of all late storms that carry it, and how well the factors together separate late from early storms
(leave-one-season-out AUC and Brier skill).

Prediction (agenda, written before data): the late share is higher when the gust maximum is more than 400 km from the centre (T1) and
higher in the Atlantic north of 60N (T2). If both fail, the late-onset storms are an ordinary phase variant. My own added
expectation, no evidence for it yet: marginal storms (T3) are late more often, because a storm that barely crosses the threshold crosses it
only near its gust peak, which is scattered about the pressure minimum, while a strong storm is above the threshold for a long window
that usually contains the minimum. A null on all three is plausible.

## Data (all committed; no pull)

- Events: `hf_history/results/lifecycle_events.csv` (one row per pipeline A event, from `lifecycle.py`), rows with `era == "2004+"`
  and `tc == False`: **1,886 events** (Atlantic 1,039, Pacific 847), 22 seasons 2004-05 to 2025-26 (decision 1: gust-based, 2004-05
  on; no pre-2004 season is used, because "onset" is a gust threshold crossing and pre-2004 gust-based onset times cannot be
  defended. A pressure-only analogue is not defined for onset, so there is no pre-2001 replication).
- Storm structure at the onset fix: `hf_structure/results/fixes.csv`, the first row (earliest time) of each track, which is the
  first in-domain HF fix, the onset of `lifecycle.py`. Check (not an outcome check): onset time in `fixes.csv` equals `t_on` of the
  events table for every event; any that differ are dropped and counted.
- Peak event gust: `gust800_kt` of `era5_hf_catalog.csv`.
- Environment (S family only): `intensity/results/env_2004.csv.gz` at the last 00 or 12 UTC fix at or before onset (at most 6 h earlier).

## Outcome

`h_on_minp` of `lifecycle_events.csv`: onset time minus time of the lowest MSLP fix of the whole track (hours, 6 h resolution).

- **LATE (primary)** = `h_on_minp >= 0`, the 29%.
- LATE6 = `h_on_minp >= 6` (onset strictly after the minimum; about 11% by PR 16) and LATE12 = `h_on_minp >= 12`, as the
  agenda's "windows (at least 6 h and at least 12 h)" for the 6-hourly quantisation. Sensitivity only.

## Predictors (all at the onset fix unless stated; thresholds fixed here, not tuned)

Prevalences below are of the predictor alone.

| id | predictor | definition | expected sign on LATE | prevalence |
|---|---|---|---|---|
| T1 | far gust maximum | `gmax_r > 400` km (distance of the gust maximum from the centre) | + (agenda) | Atl 28.2%, Pac 17.1%; 438 events |
| T2 | high-latitude Atlantic | Atlantic events only; `gmax_lat >= 60` N | + (agenda) | 27.2% of Atlantic; 283 events |
| T3 | marginal storm | event peak `gust800_kt < 76.7` (within 5 kt of 71.7) | + (my addition) | Atl 44.0%, Pac 46.6% |
| T4 | gust maximum in the rear half | `90 <= gmax_rel < 270` (relative to motion) | + (post-occlusion wind) | Atl 66.8%, Pac 75.7% |
| T5 | explosive deepening | track max 24 h deepening `maxdeep >= 1` Bergeron (38 events without it are dropped from T5) | - (bombs reach HF while deepening) | Atl 64.3%, Pac 77.4% |
| T6 | translation speed | `speed_kt` at onset, per SD (SD 13.0 kt, mean 27.4) | two-sided, no prior | continuous |
| T7 | basin | Pacific against Atlantic | - (PR 16 descriptive: Pacific 23% against Atlantic 34% late-or-at) | 847 / 1,039 |

Secondary exploratory family **S1-S9**: the nine predictors `B`, `VTL`, `VTU`, `jet250`, `div300`, `vadv500`, `eady`, `flux`, `tcwv` of
`env_2004.csv.gz`, each standardised within basin, per SD, two-sided, no prior. (`sst` and `sstgrad` have 11-20% missing and are not
used.) These are descriptive of the state at onset and several are partly set by the storm itself.

## Model and test

For each of T1-T7 one logistic regression of LATE on the predictor, with **basin** and **onset-month group** (Sep-Nov, Dec-Feb, Mar-May)
effects (T2: Atlantic only, no basin term; T7: no basin term). Cluster-robust (CR1) sandwich standard errors with **season as the
cluster** (22 clusters; effective n is seasons, not storms), Wald test with t on 21 degrees of freedom, two-sided p. Report the odds
ratio with 95% interval, the late share inside and outside the factor (basin-month standardised to the pooled sample), and the share of
all late storms that carry the factor. A season-block bootstrap (2,000 draws) interval is reported alongside as a check on the
sandwich interval; the sandwich p is the test.

Multiplicity: Benjamini-Hochberg q within the primary family T1-T7, within S1-S9, and across all 16. I report the count passing q < 0.05
out of the total in each.

### Decision rules (thresholds, not judgement)

- A factor is **supported** if its test has q < 0.05 (within its family) in the expected sign (T6 and S: either sign, labelled a lead
  only).
- A factor is a **well-powered no** if the 95% interval of its odds ratio lies within [0.67, 1.5] and the planted-effect power
  at OR = 1.5 (below) is at least 80%. Otherwise it is **can't tell**, and I state the minimum detectable OR.
- Mechanism verdict: T1 supported means "displaced gust maximum"; T2 supported means "high-latitude Atlantic / orographic". If neither T1
  nor T2 is supported, the agenda rule reads "ordinary phase variant"; I add T3 to that reading (threshold artefact) only if T3 is
  supported.
- **How much is explained**, from a joint logistic model of T1, T3, T4, T5, T6 and the basin and month terms (T2 added in a second,
  Atlantic-only fit), scored leave-one-season-out: AUC at least 0.70 means the late storms are largely identified by these factors; 0.60 to
  0.70 partly; below 0.60 they are an ordinary phase variant whatever the single tests say. Brier skill against the basin-month base rate is
  reported next to it.

## Power

Planted-effect simulation (2,000 draws): outcomes drawn from the fitted null (basin-month rates) with a planted log-odds shift on each
predictor, season-clustered Wald test as above. Report power at OR 1.5 and 2.0 and the minimum detectable OR at 80% for every
test; a null without it is "can't tell". Detectable effects are about what the group sizes imply (smallest group: T2, 283 events).

## Descriptive block (no tests, no multiplicity)

- Decomposition of onset-minus-minimum: peak-gust-time minus minimum-pressure time against onset-time minus peak-gust-time; the share
  of late storms that have a single HF fix; the share of late storms that would be early with a 6 h tolerance.
- Late share by basin, month group, number of HF fixes, and stage composition; the distribution of peak-gust-minus-minimum-pressure
  time (to see whether late storms are the tail of one continuous distribution or a separate class).
- Heterogeneity: each T1-T5 estimate within each basin, and in 2004-14 against 2015-25 (sign agreement only; a stability check inside the
  same look, not a replication).

## Sensitivity (reported in full, not counted in the families)

LATE6 and LATE12 as outcomes; events with at least 2 HF fixes; Atlantic events with `gmax_lat >= 60` removed; Atlantic events
with the terrain flag (gust maximum within 100 km of Greenland or Iceland) removed.

## Held-out look

No model is trained on one set of seasons and scored on another, and no pattern is selected. All 22 seasons are used once, so I
count **one look at 2015-25 (look 11 of this log by the coordinator's count)** and **none at pre-2001 seasons**. Logged in
`hemispheric/results/heldout_looks.log` with a UTC timestamp.

## What this cannot show

- Whether the late storms differ in the real atmosphere: onset and minimum pressure are ERA5 quantities at 6 h, from a 31 km
  reanalysis that under-resolves the gust peak (README of `hf_structure`: ERA5 sustained wind reaches 64 kt at about 1% of these fixes).
- Causes. The tests are of association at the onset fix; "far gust maximum" and "late life-cycle stage" cannot be separated, because the
  gust maximum moves out as the storm fills (RMG 175 km deepening, 337 km filling in the Atlantic, PR 27).
- The Atlantic north of 60N has no terrain mask; T2 and the Greenland sensitivity are the only handling.
- Nothing here transfers to the archive's own HF fixes, which are OPC's judgement.

## Deviations (post hoc)

None yet.
