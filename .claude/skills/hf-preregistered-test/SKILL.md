---
name: hf-preregistered-test
description: Use when planning, running or reporting any hypothesis test on HF-low counts, position or intensity against teleconnections or other predictors, so the test is pre-registered, lagged, FDR-corrected and honest about power.
---

# Pre-registered test protocol (hf-low)

Load `hf-harness` first. This is the statistics layer for a question of the form "does X change HF lows?". Worked examples on the research branch: `research/era5/{freq_split,enso_pna,enso_kuroshio,tele_combos,hemispheric,clustering,seasonal_cycle,jet_trough,trough_bombs}/PREREGISTRATION.md`. Copy the nearest.

Jason's principle: a null is plausible, and a well-powered null is a full result. Nothing is picked after seeing outcomes.

## 1. Write the plan before any outcome is joined to a predictor

Commit `research/era5/<topic>/PREREGISTRATION.md` first and cite its hash in the README. State what had and had not been looked at (predictor extraction and overlap checks are fine; counts and fits are not). Sections:

- Question in Jason's words, then the answer you will give: yes, no (well-powered) or can't tell, and by how much.
- Outcome (archive primary; pipeline A proxy secondary), sample, seasons, basin split, unit (week, track, season).
- Predictors and lag. Default: index averaged over days -10..-4 before the track's first fix; weekly designs use the 7 days before the window starts (second lag 14..8). A same-time run is a labelled contrast only, because storms feed the indices. If a predictor must be read at the fix time (trough depth), call it a same-time association and register a lead test (outcome 24-48 h on) as a secondary. (RA-24: 1.21 same-time, 0.92 at lead 24, filed only as a sensitivity.)
- Model, covariates (month effects, previous-week count where storms could imprint on fields), comparison models.
- Primary and secondary tests and the multiplicity families.
- Decision rules as thresholds (q < 0.05 for yes; no only if 80% power at a stated effect; otherwise can't tell).
- Power plan and what the study cannot show.
- **Pull plan.** GB per variable from a measured chunk size, counting every component the sensitivities need (a jet speed needs U and V), dry runs and checks. State the approved GB; stop and ask if the total will pass it by over 10%. (RA-24: 105 GB streamed against 85 approved; a 9.6 GB run was redone for V250.)
- **Known conflicts from earlier work.** Read the caveats of the result you build on. If it showed two definitions disagreeing, or a stronger correlated predictor, register an agreement rule on size (not sign alone), a place control (position-cell fixed effects) and the overlapping predictor as secondary tests, and how the verdict reads if the effect shrinks with the overlap added. (RA-24 met all three post hoc though PR 78 had flagged them: D1 1.21 fell to 1.05 with the jet; D2 was 1.75.)
- An empty "Deviations (post hoc)" section; log in order, tagged before any outcome or after the primary result.

## 2. Fixed project rules

- **Seasons.** Fit and test on 2004-05 onward (decision 1). Pre-2004 gust-based pipeline A counts only for variation within that era, with an era term; levels and trends start 2001-02; no 1979-2025 HF trend. Pressure-depth counts may go back to 1979 (cuts 966.2 Atlantic, 965.0 Pacific hPa, fixed in `freq_split`).
- **Effective n** = seasons x independent predictor values per season, never cyclone or fix count. ONI is about one value per season. Resample or permute by season block (cluster by season when fixes are pooled); resample storms when fixes within a storm are used.
- **Distinct regions only** when combining teleconnections (NAO, PNA, MJO as a longitude not RMM phases, ENSO via monthly ONI).
- **Splits.** If a hold-out is used, freeze the models and commit them before the single look, log every look (`heldout_looks.log` or the topic's file under `looks/`) and state the count. A swap-split replicate is not independent.
- **Name the pipeline and say proxy.** Pipeline A = `research/era5/hf_history`, B = `event_fields/criterion/series`. Say whether transitioning tropical cyclones are in. Atlantic fixes north of 60N are a lower-confidence group.
- **Re-running a pipeline in a new era.** Pre-register an outcome-free recovery check against committed points (RA-24: 99% within 0.25 degree and 0.2 hPa; got 100% of 71,472); stop if it fails. Write tables at full precision and check the table the analysis uses (RA-24's first fix table, rounded to 4 digits, was re-linked).

## 3. Report every test

- Every pre-registered test appears, with p, Benjamini-Hochberg q within its family, and q across all tests. Count tests passing FDR out of the total (0 of 55).
- Effect sizes with intervals, per SD of the index, on the scale a forecaster reads. Beside a "yes", show the estimate adjusted for the overlapping predictor.
- **Power and minimum detectable effect** by planted-effect simulation (`trough_bombs/power.py`, `clustering/power.py`, `hemispheric/`, `freq_split/fdr_power.py`), run before any coefficient is read. If draws ignore season clustering, call power an upper bound and show the real interval width. A null without a stated MDE is "inconclusive", not "no effect".
- Anything not in the plan is "post hoc", labelled in the table and the prose, with both versions shown. A post hoc lead is a lead, not a finding; the pre-registered verdict stands.

## 4. Before you call it done

Hand off to `hf-result-closeout` (fresh verifier, ledger PR, Q&A forward). A README that quotes numbers lists which were independently recomputed and which were not; put the power run and bootstrap intervals in the verifier's scope (RA-24 left both unchecked).
