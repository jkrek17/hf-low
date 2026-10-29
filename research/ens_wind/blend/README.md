# GEFS + ECMWF ENS blend of 10 m wind chances (winter 2025-26)

Asked by Jason 2026-10-29 as the follow-up to `../` (PR 285). Plan: `PREREGISTRATION.md`, committed before the 2025-26 maps were re-pulled (the PR 285 maps were kept in scratch). Sustained 10 m model wind, no ERA5. Plain-language version first, method after.

## What to take from it
- **The blend beats GEFS alone on weeks it was not fitted on, but only for 48 kt.** Pre-registered test (mean over 16 cells of the Brier skill difference, leave-one-week-out, 37 run dates): blend minus GEFS +0.033, interval +0.016 to +0.050, so the rule is met.
- **48 kt, checked against ASCAT (cells ASCAT saw):** blend minus GEFS is positive in all eight cells, +0.02 to +0.10; six of eight have intervals above zero and survive BH at q 0.05 (all Atlantic; Pacific days 2 and 3); Pacific days 1 and 4 cross zero. A plain 50/50 average does as well as the fitted weights (+0.034). Fitted weight on GEFS is about 0.3 to 0.6 for most cells, so ENS is carrying more than half in most of them (Pacific day 3: 0.1).
- **64 kt, checked against OPC hurricane-force lows within 324 nm: no gain.** Differences run -0.02 to +0.03 and no interval is above zero (Atlantic day 1 only just touches it, -0.007 to 0.000); on days 1 and 2 the fitted weight goes to GEFS (0.8 to 1.0). With 100 to 340 events per cell in one winter this is "cannot tell", not "no benefit".
- **Reading the blended chance (held-out):** 48 kt: under 5 % happened 0.1 %; 5 to 20 % happened 11 %; 20 to 50 % happened 37 %; 50 to 80 % happened 73 %; over 80 % happened 94 %. 64 kt: under 5 % 0.4 %; 5 to 20 % 12 %; 20 to 50 % 45 %; 50 to 80 % 65 %; over 80 % 89 %.
- **Rough weights, one winter.** ENS has only 2025-26. Weights and the ENS calibration are rough; refit when winter 2026-27 is in (planned spring 2027).

## Caveats that matter
- The GEFS baseline table (PR 285) was fitted on all six winters including this one, which favours GEFS, so the gain is if anything understated. Refitting GEFS on the training weeks only gives the same GEFS skill (`gefs-recal` in the results file), and the blend built on it gains +0.039 [+0.024, +0.058].
- Placebo (truth taken 23 days off): the like-for-like version (blend on recalibrated GEFS minus recalibrated GEFS) is -0.001 [-0.007, +0.001]. The primary metric is not meaningful under the placebo, because the committed GEFS table is, by construction, miscalibrated against unrelated truth. 48 kt ASCAT cells at lead 48 had no ASCAT file on the shifted days and were dropped (`results/placebo_23d.txt`).
- 48 kt truth is ASCAT stress-equivalent wind (no COARE), saturating near 50 kt, seeing about half of any area. 64 kt truth is an analyst neighbourhood event, "OPC agreed".
- 37 run dates, every 5th day; intervals resample ISO weeks; neighbouring run dates share storms. Fitting and testing inside one winter is why leave-one-week-out is used.
- Marine wind project (GFE blend) numbers are reported there, not rechecked here.

## Checked
A fresh Sonnet agent, with its own code and the written spec only, recomputed cases, events, GEFS-alone skill and the fixed 50/50 blend minus GEFS for four cells (Atlantic 48 kt day 2, Pacific 48 kt day 3, Atlantic 64 kt day 3, Pacific 64 kt day 2). Cases, events and GEFS skill matched exactly; the 50/50 differences matched to within 0.002 (+0.073 vs +0.071, +0.085 vs +0.085, +0.031 vs +0.025, -0.024 vs -0.024; the small gaps come from how the ENS k = 0 rate is smoothed). **Not rechecked:** the fitted weights, bootstrap intervals, BH, the primary mean, reliability bands, the placebo. Run dates are 37 (the last is 13 Apr 2026, every 5th day from 15 Oct).

## Files
`build_joint.py` (joint GEFS x ENS count tables from the pulled cell maxima and truth), `blend.py` (leave-one-week-out scorer and bootstrap), `test_blend.py` (planted-effect test, no network), `results/blend_2025-26.{txt,json}`, `results/placebo_23d.txt`. Pulled cell maxima (61 MB) are not committed; they are in `/mnt/project-files/ens-wind-blend/` (`joint.npz` is the table `blend.py` reads). Reproduce: `ENS_WORK=<dir> python3 ../ens_fetch.py` for 2025-26 dates, `../ascat_truth.py 2025`, then `build_joint.py` and `blend.py JOINT.npz ../results/gefs_primary.json results/blend_2025-26`.
