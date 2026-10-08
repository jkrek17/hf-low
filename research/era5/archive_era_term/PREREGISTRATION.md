# RA-29: do finished results that score archive HF by fixes per event survive an era term? (pre-registration)

Written 2026-10-08 before any era-adjusted outcome was computed. Everything ERA5 is a **proxy**; the pipeline is named where used
(A = `research/era5/hf_history`). No ERA5 pull. Deviations after the first outcome go under "Deviations (post hoc)" at the end.

## Question and the answer I will give

RA-11 (PR 105, `fixes_per_event/`) found archive HF fixes per event fell by about 0.35 in one step at the end of QuikSCAT
(2009-11-23); event counts and total fixes are flat; the ERA5 proxy's fixes per event rose. For each finished result that scores
the archive by fixes per event, hours at HF, or "at least 2 or 3 fixes", does the estimate keep its sign and size when a
before/after-2009-11-23 term is added? Plain-words answer per result: unchanged / changed (by how much) / no longer resolved.

## Audit: what was read to build the list (before outcomes; reading code and READMEs only)

Scripts that read archive HF fixes (`grep HF_Data|hf-lows.json|events_archive|linked_ids`) were each checked for a per-event fix
measure.

**Re-fitted (archive per-event fix measure):**
- T-A. PR 79 / RA-22 (`sustained_hf/`), archive design (PR 41 weekly archive counts, hemispheric pattern index): HF_k with k = 1, 2, 3,
  brief (exactly 1 fix), middle (exactly 2), and the paired contrasts. Certain.
- T-B. Count reconciliation (`count_reconcile/`): archive events per year at >= 1, >= 2, >= 3 HF fixes (the "duration" step of the 46.0 -> 33.8 -> 22.4 chain) by era.
- T-C. Climatology atlas (`climatology_atlas/`): archive hours at HF per event (the 12 h / 18 h medians against the proxy's 12 h), by era.
- T-D. HF skill limits (`hf_skill_limits/`): the archive-label ceiling (ERA5 hf24 label scored against the archive's "any HF fix in the next 24 h"), by era. The label is built from archive fixes, so fewer fixes per event changes its positives.
- T-E. HF probability tracks (`hf_probability_tracks/`): the archive-listed share of proxy events (a match to an archive HF fix at the peak position) by era.

**Excluded, with the reason (listed so nobody reruns them for this):**
- Counts and total fixes, flat by RA-11: PR 41 headline (k = 1), weekly counts, basin counts (RA-25), pair counts (RA-12), clustering (PR 37), seasonal cycle (PR 45), `front_hours` S8 (archive fixes in a box per season; a total, already has a 2013-14 dummy).
- Proxy-side fix counts (the proxy's fixes per event rose, RA-11 S7): the proxy half of PR 79 / RA-22 (PR 14, 63, 64 at k = 2, 3), `lifecycle.py` durations, tele_intensity HF duration, late_wind (PR 86), threshold axis, gust scaling.
- Labels that never use archive fixes: PR 12/56/68 (ERA5 labels; T-D covers the one archive-label table), RA-20, RA-27, RA-9, P(HF) tracks other than T-E.
- Archive late-onset check: already done in RA-11 addendum A.
- Atlas distance while HF and speed: derived from the same fix count as T-C; reported as a note, not tested separately.

Not claimed exhaustive: a result I did not find reading `research/era5/` could still use an archive fix count.

## Era definition

Step date 2009-11-23 (QuikSCAT end, verified in RA-11). Weekly designs: era = 1 for weeks whose start date is on or after 2009-11-23
(first post week 2009-11-26). Event designs: era = 1 if the event's first HF fix is on or after 2009-11-23 (proxy: first HF time).
Sample: seasons 2004-05..2025-26 (Decision 1; nothing earlier). Pre-step is about 5.5 seasons, post about 16.5.

## T-A design (primary)

Exactly the model of `sustained_hf/contrast.py archive` (Poisson, month dummies, log1p previous-week count, lagged leave-one-season-out
pattern index, 660 weeks per basin), plus one era column. Control: with the era column dropped, the code must reproduce
`sustained_hf/results/contrast_archive_withm.csv` (`log_est`) to 4 decimals before the era-adjusted run is read.
Outcomes h1, h2, h3, b, m; pairs h2-h1, h3-h1, h3-b. 2000 season-block bootstrap resamples (percentile 95%), 2000 season-block
permutations of the index (p), season-clustered SE for the MDE (2.8 x SE). Also logged for every class: the era coefficient.
Families (BH): F1 headline = RR(h1), RR(h2), RR(h3) x 2 basins (6); F2 = sustained minus brief (h3-b) x 2 (2); F3 = other pairs and b, m
(descriptive, 10). Reported also across all.

**Prediction (committed).** Era coefficients: h1 about 0 (counts), h2 and h3 negative, b positive. Each RR and pair keeps its sign and moves
by less than its interval: the era-adjusted point lies inside the original 95% bootstrap interval. The index is not collinear with the
era (|corr| small), so intervals widen by less than 1.3 times.

**Decision rules.** A cell is
- *unchanged* if the sign is kept, the adjusted point is inside the original interval, and its pass/fail at q < 0.05 (F1, F2) is unchanged;
- *changed* if the sign flips or the adjusted point leaves the original interval (the ledger sentence changes by the amount);
- *no longer resolved* if the sign and size are kept but q passes originally and fails adjusted (power, not direction).
- The overall verdict for PR 79's archive side is "holds with an era term" only if all four F1 cells at k = 2, 3 are unchanged or no worse than "no longer resolved" with the sign kept.
Power: report the adjusted SE and MDE next to the original; the era term is identified by about 5.5 pre seasons, so say so.

Sensitivity (labelled secondary, all reported, no re-picking): S1 step at ASCAT-B (2012-09-01) instead; S2 a season-linear trend instead of the step; S3 era and index
interaction is NOT fitted (underpowered).

## T-B to T-E (descriptive era splits, same step date; each with a season-block bootstrap, 95%)

- T-B: archive events per year at HF fixes >= 1, >= 2, >= 3, by basin and era (rate = count / years of exposure in the era within 1 June 2004..31 May 2026), plus the post/pre ratio. Prediction: ratio at >= 1 within 0.9 to 1.1; at >= 3 below 1. Also the season counts that the published comparison uses (Von Ahn et al. and Jelenak et al. cover QuikSCAT years): the pre-step >= 3 rate against the post-step rate, to say which era the "about 20" published counts should be compared with.
- T-C: archive mean and median hours at HF per event by era and basin, and the archive minus proxy difference (proxy events from `events_proxy.csv`, same eras). Prediction: archive mean falls by about 2 h (0.35 fixes x 6 h) post; proxy does not fall; the gap narrows but archive stays above the proxy median 12 h.
- T-D: HSS of ERA5 hf24 label against the archive label, by era, 2004-05..2025-26 fix-level table, with the PR 68 control (overall 0.66 reproduced before splitting). Prediction: HSS lower post (fewer later fixes), within the overall interval.
- T-E: share of proxy events archive-listed (2004+ events), by era. Prediction: not lower post; if it is, the matching rule depends on fixes per event.
Tests in T-B..T-E: each pre/post contrast is one test; BH over the 4 studies' contrasts together (listed in the README, with counts).

## Looks

The archive-side re-fits score no held-out block: no model is fitted on one block and scored on another, no entry in any looks log. T-D re-reads the 2004-05..2025-26 fix table it was originally scored on with fixed labels and no refit. If the plan changes so that something is scored on 2015-25, it will be said first and logged in `research/era5/looks/archive_era_term.log`.

## What this cannot show

Why the step happened (RA-11 left that open); anything about the proxy side; whether the era term absorbs other changes in archive practice that coincide with 2009-11. The pre-step sample is 5.5 seasons, so an era-adjusted interval is wider and "no longer resolved" is a power statement.

## Deviations (post hoc)

None yet.
