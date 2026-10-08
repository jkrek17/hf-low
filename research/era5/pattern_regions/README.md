# Upstream or local: where does the lagged hemispheric pattern's skill sit? (ERA5 proxy fields, archive outcome)

Plan, committed before any region model was fitted: [PREREGISTRATION.md](PREREGISTRATION.md) (`06ff6af`). Builds on
`research/era5/hemispheric/` (hf-low PR 41). Code: `regions.py` (sectors, region-restricted EOFs), `ablate.py` (nested
leave-one-season-out scoring), `stats.py` (skill, shares, Shapley, bootstrap, permutation, BH), `r1.py` (replication).
Results: [results/](results/) (`summary.json`, `tests.csv`, `window_curve.csv/.png`, `r1.json`, `looks.log`, and the per-model
scores in `results/loso/`, 1.3 MB, enough to recompute everything but R1 without the fields).

Outcome: the archive's weekly counts of hurricane-force lows per basin (Oct-Apr, 22 seasons 2004-05..2025-26). Predictors: ERA5
fields, **a proxy** for the atmosphere (Z500, U250, SST at 5.6 degrees, 7-day mean over days -7..-1 before the week).
Pipeline A (`research/era5/hf_history`, a proxy) appears only as the replication outcome (R1: depth counts 1979-2000) and as a
covariate in the imprint check. No gust-based pipeline A count before 2004 is used.

## Answer

**Atlantic: upstream and local cannot be separated. The skill is spread redundantly across the North Pacific, North America and
the Atlantic itself; none of the three regions is needed on its own, and Eurasia (downstream) carries none. The agenda's
prediction (upstream keeps at least half, local less than half) is half met and not supported as a statement of remote
control. Pacific: local.** (Pipeline-free statement of what was cut: fields only; the outcome is the archive.)

### Atlantic (full pattern: +5.74% leave-one-season-out deviance skill, p = 0.0001; 19 of 22 seasons beat the baseline in PR 41)

| fields kept (all of Z500, U250, SST) | skill | share of full (90% interval) | Shapley share (90%) |
|---|---|---|---|
| upstream, 140E to 100W (N Pacific, Alaska, western N America) | +4.18% | **73%** (45-105%) | 41% (11-69%) |
| local, 100W to 20E (the basin; the 90W trough core and the whole jet) | +5.37% | **94%** (48-145%) | 57% (27-87%) |
| downstream, 20E to 140E (Europe, Siberia) | +0.09% (p = 0.08) | 2% (-17 to 23%) | 2% (-13 to 20%) |
| drop upstream (local + downstream) | +5.39% | loses 0.36 points, p = 0.39 | |
| drop local (upstream + downstream) | +4.84% | loses 0.91 points, p = 0.18 | |
| drop downstream | +5.80% | gains 0.05 | |

- By the pre-registered rule this is **SHARED** (both keep-only shares at or above one half), and **unresolved**: both 90%
  intervals straddle 0.5. Upstream minus local keep-only skill is -1.2 points (95% interval includes 0, p = 0.56; the smallest
  difference this record resolves at 80% power is 5.7 points). Both beat downstream (4.1 and 5.3 points, p = 0.010 and 0.016).
- Dropping either region costs under a point and neither loss is detectable. That is the signature of a wave train seen in
  more than one place, not of one region being the source.
- **Trough upstream (V1, pre-registered variant):** put the North American trough core (to 70W) in the upstream sector and
  upstream-only keeps **84%**, while the rest of the Atlantic (70W to 20E) alone keeps **26%**. The skill sits in the North Pacific to
  North America half of the pattern, not in the basin east of the coast.
- **Window curve** (120-degree windows, `results/window_curve.png`): skill peaks for windows centred 135W to 112W (89-94% of full),
  is still 81% for one centred on 45W, and falls to 0-15% for windows centred from 0 to 112E. Smooth, one broad maximum over the
  Pacific-America-Atlantic arc.
- **Variables:** Z500 alone keeps 79% upstream and 110% local; U250 alone keeps 76% upstream and 66% local (local jet).
- **Storm imprint (post-registered check S-imprint, baseline adds the previous week's all-cyclone count):** full +6.16%; upstream
  keeps 68%, local 84% of that. Not a clustering artefact of the local fields.
- **Replication on 1979-80..2000-01 (R1, pipeline A depth counts, a proxy outcome; within-era variation only):** the refitted
  full pattern has rate ratio 1.25 per SD (p = 0.0002). Upstream-only keeps **30%** (1-48%), local-only **75%** (37-131%), downstream
  15%. Local is the larger share, as in the primary result, but upstream-only is below one half here, so the primary's
  "upstream keeps at least half" does not replicate. By the pre-registered rule (same larger of the two) R1 agrees on direction;
  by category it reads LOCAL where the primary reads SHARED. Read together: local fields carry at least as much as upstream
  ones, and an upstream-only account is not supported.

### Pacific (full pattern: +3.40%, p = 0.0001; 95% interval -0.4 to +6.5%, so ratios are wide)

| fields kept | skill | share of full (90%) | Shapley share (90%) |
|---|---|---|---|
| upstream, Eurasia 0-120E | 0.00% (no component survived the ridge penalty in any season) | 0% | -17% (-69 to -1%) |
| local, 120E to 120W | +3.88% | **114%** (73-214%) | 99% (75-155%) |
| downstream, 120W to 0 (N America, Atlantic) | +0.88% | 26% (-28 to 75%) | 18% (-13 to 52%) |
| drop local | +0.34% | loses 3.06 points, p = 0.039 (q 0.13) | |
| drop upstream | +4.33% | gains 0.93 | |

- **LOCAL, resolved:** local keeps at least half with its whole interval above one half, upstream keeps none. Local-only is
  better than the full pattern (extra fields add noise at this skill level). The dateline trough and the extended jet are
  local; the central North American ridge in PR 41's map is downstream and adds little alone.
- Replication on 1979-2000: local-only rate ratio 1.17 per SD, 122% of the full pattern's deviance gain (full 1.15); downstream
  only 28%; Eurasia has no pattern to test. Agrees.
- The Pacific full-pattern interval includes 0 and the keep-only differences are mostly short of resolution (smallest resolvable
  difference 4.4-5.0 points): the direction is clear, the sizes are not.

## Tests and multiplicity

74 pre-registered tests in the leave-one-season-out set (38 Atlantic, 36 Pacific; `results/tests.csv`), Benjamini-Hochberg across all
74: **49 pass q < 0.05** (Atlantic 31 of 38 at basin level, Pacific 17 of 36). Skill tests (SS > 0) that fail are mostly sliding windows
over Eurasia and the Atlantic-Pacific gap, which is the result. Family 2 (the pairwise and drop-one comparisons, 6 per
basin): Atlantic upstream-versus-downstream and local-versus-downstream pass (q = 0.047 each), nothing else passes; Pacific none
pass (smallest q = 0.13). R1: 13 tests (the Pacific upstream model has no pattern to test), all p <= 0.015, all pass BH.

## What was and was not held out

- **No fresh archive split exists.** Every archive season 2004-05..2025-26 was in a PR 41 look; this thread scored all 22 again,
  once, under nested leave-one-season-out (logged in `results/looks.log`). The comparison is between models that are each
  out of sample on every season and none was chosen from outcomes, but the absolute skill is not new evidence for the pattern.
- **FULL reproduced PR 41's leave-one-season-out numbers exactly** (+5.74% Atlantic, +3.40% Pacific; the pre-condition).
- The 1979-2000 depth counts were a second look at that outcome for the full pattern (PR 41's S5 was the first) and a first look
  for the sector models. The 2026-27 season is the only fresh archive data (1 of 30 weeks at the time of writing); re-running
  `ablate.py` and `stats.py` after it adds a seventh of a season of power per month. It will not settle the Atlantic split.

## Limits

- Not causal, and one week is long enough for a wave packet to circle much of the hemisphere, so "upstream" and "downstream" are
  labels for where the fields are, not proof of direction of propagation. Downstream (Eurasia, for the Atlantic) carrying nothing is
  the cleanest part of the result.
- Sector edges are a choice. V1 (trough upstream) and the window curve show the Atlantic split moves with them: local keeps 94% or
  26% depending on whether the 90W trough is counted as local.
- Skill is modest (a weekly count is mostly noise); intervals for shares are wide because 22 seasons carry 3-6 points of skill.
- Atlantic fixes near Greenland and the archive's recording process are unchanged confounds from PR 41.

## Deviations (post hoc)

None to estimates or tests. Implementation note: the inner leave-one-season-out penalty search is `fast_loso` in `ablate.py`
(design matrix built once); it equals `hemlib.loso_lambda` to 1e-13 on a test fold. Pacific V1 was not run, as planned.

## Reproduce

    python3 research/era5/hemispheric/fields.py FIELDS          # 16.7 GB, predictors only
    OMP_NUM_THREADS=1 python3 research/era5/pattern_regions/ablate.py FIELDS OUT 4    # about 15 min, 4 cores
    python3 research/era5/pattern_regions/stats.py OUT research/era5/pattern_regions/results
    python3 research/era5/pattern_regions/r1.py FIELDS research/era5/pattern_regions/results
    # or recompute from results/loso/ with stats.py (copy the npz files to OUT)

## Verification

A fresh Sonnet agent that had not seen the code or this README recomputed from `results/loso/` with its own implementation
(own bootstrap seed): **matched** FULL skill (both basins), keep-only and drop-one skills (all 12), retained shares, Shapley shares
(both basins), the FULL and f_UP / f_LOC Atlantic intervals (within 0.02), seven sliding-window skills, seasons beating the baseline
(19 of 22 Atlantic, 14 of 22 Pacific), and the sector tiling (Atlantic 22/21/21 longitudes, Pacific 21/21/22, no gap or overlap).

**Not independently checked:** the R1 numbers (rate ratios, deviance gains, shares and p; they need the fields), the V1, V2 and
imprint values, the permutation p and BH q values, the Family 2 p values and minimum resolvable differences, the Pacific share
intervals, the Shapley intervals, the split-half table in `summary.json`, and the plain-language reading.
