---
name: hf-preregistered-test
description: Use when planning, running or reporting any hypothesis test on HF-low counts, position or intensity against teleconnections or other predictors, so the test is pre-registered, lagged, FDR-corrected and honest about power.
---

# Pre-registered test protocol (hf-low)

Load `hf-harness` first (orientation, ledger, ERA5 proxy rules). This skill is the statistics layer for a question of the form "does X change HF lows?". Worked examples on the research branch: `research/era5/{freq_split,enso_pna,enso_kuroshio,tele_combos,hemispheric,clustering,seasonal_cycle}/PREREGISTRATION.md`. Copy the structure of the nearest one.

Jason's principle: a null is plausible, and a well-powered null is a full result. Nothing is picked after seeing outcomes.

## 1. Write the plan before any outcome is joined to a predictor

Commit `research/era5/<topic>/PREREGISTRATION.md` first and cite its commit hash in the README. State in it what had and had not been looked at (predictor extraction and overlap checks are fine; counts and fits are not). Sections used so far:

- Question in Jason's words, then the answer you will give in plain words: yes, no (well-powered), or can't tell, and by how much.
- Outcome (archive primary; pipeline A proxy secondary), sample, seasons, basin split, unit (week, track, season).
- Predictors and the lag. Default: index averaged over days -10..-4 before the track's first fix; weekly designs use the 7 days before the window starts (second lag 14..8). A same-time run is a labelled contrast only, because storms feed the indices.
- Model, covariates (month effects, previous-week count where storms could imprint on fields), comparison models.
- Primary tests, secondary tests, and the multiplicity families.
- Decision rules, written as thresholds (q < 0.05 for yes; no only if 80% power at a stated effect; otherwise can't tell).
- Power plan and what the study cannot show.
- An empty "Deviations (post hoc)" section.

## 2. Fixed project rules

- **Seasons.** Fit and test on 2004-05 onward (decision 1). Pre-2004 gust-based pipeline A counts only for variation within that era, with an era term; levels and trends start 2001-02; no 1979-2025 HF trend. Pressure-depth counts may go back to 1979 (cuts 966.2 hPa Atlantic, 965.0 hPa Pacific, fixed in `freq_split`, not re-chosen).
- **Effective n** = seasons x independent predictor values per season, never cyclone or fix count. ONI is about one value per season. Resample or permute by season block; resample storms when fixes within a storm are used.
- **Distinct regions only** when combining teleconnections (NAO, PNA, MJO as a longitude not RMM phases, ENSO via monthly ONI).
- **Splits.** If a hold-out is used, freeze the models and commit them before the single look, log every look (`heldout_looks.log`) and state the count. A swap-split replicate is a replication, not independent.
- **Name the pipeline and say proxy.** Pipeline A = `research/era5/hf_history`, B = `event_fields/criterion/series`. Say whether transitioning tropical cyclones are in. Atlantic fixes north of 60N are a lower-confidence group.

## 3. Report every test

- Every pre-registered test appears, with its p, Benjamini-Hochberg q within its family, and q across all tests. Count tests passing FDR out of the total (for example 0 of 55).
- Effect sizes with intervals, per SD of the index, on the scale a forecaster reads.
- **Power and minimum detectable effect** by planted-effect simulation (`hemispheric/`, `clustering/power.py`, `freq_split/fdr_power.py`). A null without a stated MDE is "inconclusive", not "no effect".
- Anything not in the plan is "post hoc", labelled in the table and the prose, with both versions shown. A post hoc lead is a lead, not a finding.

## 4. Before you call it done

Hand off to `hf-result-closeout` (fresh verifier, ledger PR, Q&A forward). A README that quotes numbers lists which were independently recomputed and which were not.
