# Why gust-based and depth-based Pacific HF counts disagree (ERA5 proxy, pipeline A)

Plan, committed before any result: [PREREGISTRATION.md](PREREGISTRATION.md) (commit `907d21b`).
Code: [match.py](match.py) (archive to track), [diagnose.py](diagnose.py) (H1, H2), [teleconnect.py](teleconnect.py) (H3), [assemble.py](assemble.py) (Benjamini-Hochberg).
Results: [results/](results/) (`diagnose.txt` readable, `fdr_all_tests.csv` every pre-registered test with q, `tele_estimates.csv`, `tele_contrasts.csv`, `sens_*` sensitivity runs).

Reproduce (about 10 minutes, no ERA5 pull; needs statsmodels and scipy; CPC files from `/mnt/project-files/teleconnection-test/cpc_indices`):

    python3 research/era5/gust_depth/diagnose.py research/era5/hf_history/results/all_tracks.csv.gz docs/data/hf-lows.json research/era5/gust_depth/results 2000 800
    python3 research/era5/gust_depth/teleconnect.py research/era5/hf_history/results/all_tracks.csv.gz docs/data/hf-lows.json <cpc_indices dir> . research/era5/gust_depth/results 2000
    python3 research/era5/gust_depth/assemble.py research/era5/gust_depth/results

All numbers are **pipeline A** (800 km ocean gust index, 71.7 kt) and a **proxy**. **G** = gust index at least 71.7 kt. **D** = track minimum MSLP at or below a cut that gives the same count as G (965.0 hPa Pacific, 966.2 Atlantic). **D'** = depth relative to the local monthly climatology of track minima (diagnostic).

## Answer

**They disagree because a fixed pressure cut and a fixed wind cut pick different storms in the Pacific, and the archive sides with wind.**

1. **Against the archive, gust is the better stand-in in both basins (pre-registered H1, supported).** On 17 held-out seasons (2004-05 to 2020-21), Pacific AUC is 0.974 for gust and 0.953 for depth (difference +0.021, 95% 0.016-0.026); at count-matched thresholds HSS is 0.605 for gust and 0.481 for depth (+0.124, 0.088-0.162). Among real cyclones the AUC gap is +0.039 (0.030-0.049). Atlantic is the same (+0.028, +0.111, +0.040). Depth loses on both sides: POD 0.56 against 0.68 and FAR 0.54 against 0.42. (The pre-registered guess that depth's loss would be mostly false alarms is not supported; it is half and half.) The point estimate clears the pre-registered +0.02, but the lower end of the interval does not.
2. **Where they disagree (H2, supported).** At the published thresholds the Pacific has 482 storms both flag, 414 only gust flags and 404 only depth flags. Depth-only storms sit **8.8 degrees farther north** on average (peak latitude 51.0 against 42.2; 0.199 less in the western sector; Atlantic 6.9 degrees), and the archive matches only **17%** of them, against 48% of gust-only storms and 77% of the both group (difference +0.31, 0.24-0.39). The archive's HF fixes do reach 55-60N (180 Pacific fixes), so this is not a domain edge. At fixed depth the Pacific gust index falls **1.0 kt per degree of latitude** (Atlantic 0.3), is 2.3 kt lower in the western sector, and 5.1 kt higher in Nov-Feb (all q < 1e-4). In plain words: the same central pressure is a weaker wind at higher latitude and in the deep, climatologically low-pressure Aleutian/Gulf of Alaska region, so a fixed pressure cut over-counts there.
3. **Why this changes the teleconnection results (H3, partly supported).**
   - **PNA share, Pacific: yes, it is the background pressure.** Share of cyclones reaching HF, per SD of lagged PNA: gust 1.043 (0.985-1.101), depth 1.111 (1.049-1.167), location-adjusted depth D' 1.045 (0.982-1.107). Adjusting depth for where and when the storm is brings it onto the gust answer. PR 11 already found the PNA+ Pacific storm is about 0.9 hPa deeper entirely because the background is. The depth minus gust difference by itself (ratio 1.065, 0.983-1.141, q 0.22) is inconclusive; the depth-minus-D' difference is the clearer one (1.063, q 0.07).
   - **ONI total effect, Pacific: no, background does not explain it.** Per SD of ONI: gust 1.018 (0.925-1.082), depth 1.100 (1.012-1.164), D' 1.109 (1.013-1.198). Location adjustment does not move it toward gust, so the pre-registered prediction (ii) fails. Depth minus gust is 1.081 (1.009-1.157, q 0.076). Descriptively, gust-only storms respond the other way (0.887, 0.814-0.959) and storms both flag respond up (1.131, 0.985-1.243); the cause is not identified.
   - **Which one the archive's own response favours.** Archive Pacific events per SD of PNA 1.093 (1.017-1.170) and ONI 1.022 (0.908-1.105). Those are nearer gust for ONI (archive minus depth ratio 0.93, 0.84-1.01, q 0.16) and for PNA lie between the two without separating them (archive minus gust 1.02, archive minus depth 0.96, intervals include 1). Under the decision rule fixed in advance the answer for both is **can't tell, leaning gust for ONI**; the archive response is too noisy (22 seasons) to arbitrate.
   - A 47-season depth rerun reproduces the earlier PR 34 numbers (ONI 1.066, share 1.061) and shows D' at 1.090 and 1.085 there, so with more seasons the ONI effect on depth-defined counts also survives the location adjustment.
4. **Other finished results (H4, read-through, no refit)** are in the table below. Only the two named in the task show the disagreement as a number; several others used gust alone and would be exposed to the same latitude gradient.
5. **Combined definition.** Not built: the pre-registered trigger was not tested because the diagnosis already points to a clear reason and a clear preference. A logistic of the archive on gust, depth and latitude would be a natural next step if wanted; labelled exploratory if run.

## Which to use, for what

| Purpose | Use | Why |
|---|---|---|
| Counting archive-like HF lows from 2001-02 on, for any Pacific or Atlantic climatology | **Gust (G)** | Held-out HSS 0.61 against 0.48; depth over-counts high-latitude lows. |
| Teleconnection effect on the share of cyclones reaching HF, Pacific | **Gust**, or depth adjusted for local background | Raw depth turns a background-pressure shift into extra "strong" storms. |
| Anything before 2001 | **Depth** | The gust index drifts there (decision 1); depth does not. Use D' rather than raw D if a teleconnection is the question, since PNA moves the background. |
| Tip jets, barrier winds, terrain-type extremes | **Gust only** | Depth cannot see them. Gust also over-flags Greenland barrier winds north of 60N (PR 29). |
| Ranking how deep a storm is, rapid-deepening work | **Depth** | That is what it measures; wind is not the target there. |
| ONI and Pacific frequency | **Can't tell** | Gust, depth and archive responses differ by up to 8% per SD; the archive cannot separate them. |

## Other finished results (H4)

| Result | Definition used | Does the choice matter? |
|---|---|---|
| Q2 frequency split (PR 14) | both (S4) | Yes, Pacific share; explained above for PNA. |
| ENSO x PNA (PR 34) | both | Yes, ONI total; not explained by background. The interaction (RR 1.014) was ~1.00 for both. |
| Q1 intensity (PR 11) | gust events; pressure outcomes on 47 seasons | Already separates background from storm depth; consistent with H3. |
| Hemispheric patterns (PR 41) | archive outcome; proxy checks gust (+7.6% Atl, +1.6% Pac) and depth within 1979-2000 | Same sign in both; skill on the archive is the headline. |
| Near-storm framework (PR 12), life cycle (PR 16), wind structure (PR 27) | gust events; pressure for deepening | Conditional on gust-defined HF; no depth version. |
| High-latitude Atlantic (PR 29) | gust | No depth version; the 44% FAR north of 60N is a gust-specific problem. |
| Kuroshio genesis (PR 38) | pressure-tracked genesis, gust HF label | HF part inconclusive either way. |
| Combinations (PR 40) | gust; all-cyclone counts alongside | Position nulls are about gust-index peak. |
| Drift counts (PR 18) | depth as the control | Depth is the stable one before 2001. |

## What this does not show

- **Matching is at track level** from summary rows (`all_tracks` has the gust peak, not every fix): 92% of archive `low` events match at 800 km, and it recovers 88% of the 1,410 pairs already in the catalog. It was checked only against gust-matched pairs, so it could favour gust. At 500 and 1200 km, and with HF-peak archive events only, every H1 and H2 result is unchanged (Pacific dAUC +0.021 to +0.024; depth-only archive match 0.16-0.18; `results/sens_*`).
- Post hoc, not pre-registered: the latitude-band breakdown of archive-match share (depth-only is below gust-only in every band, 0.06-0.23 against 0.35-0.61 and 0.74-0.82 for both) and the group decomposition of the ONI response.
- Linear, one-index models, 22 seasons. A Pacific difference below about 1.09 per SD cannot be detected. The ONI contrast has q 0.076, below the 0.10 convention but not below 0.05.
- D' is built from this record's own track minima (all years), so it removes the average location effect, not a year-specific one. Peak position is where the gust index peaks, not the pressure minimum.

## Multiplicity

76 estimates and contrasts (H1 6, H2 18, H3 52 including the Atlantic mirror) are in one Benjamini-Hochberg family (`fdr_all_tests.csv`, columns `q_all` and within-family `q_family`). All six H1 tests have q = 0.0019. In H2, 14 of 18 have q < 0.03; the exceptions are Pacific and Atlantic median `n_fix` (q 0.69, 0.19), Pacific `log(n_fix)` (0.89) and Atlantic west-sector share (0.11). In H3 the key difference contrasts do not clear q < 0.05: PNA depth minus gust q 0.22, ONI depth minus gust q 0.076, archive minus depth ONI q 0.16, PNA depth minus D' q 0.07. The depth share under PNA (q 0.002) and the PNA response of the archive (q 0.048) do.

## Departures from the plan

None changes a pre-registered estimate. Logged because they happened after the plan was committed.
1. A 40-draw smoke run of `teleconnect.py` was looked at before the full run; nothing was changed in response.
2. The full run was repeated after adding the bootstrap p-value for the share terms (the first launch stored none); same estimates.
3. A one-line bug in `arch_counts` (index mismatch) was fixed before any result was kept.
4. Share p-values are season-block bootstrap, in the `p_perm` column; single-outcome RR p-values are permutation.

## Verification

See the end of this file once the independent check has run.
