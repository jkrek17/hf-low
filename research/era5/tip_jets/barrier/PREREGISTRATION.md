# Barrier-flow add-on to the tip jet study (written before the pull; ERA5 PROXY)

Asked by Jason 2026-10-15 ("is there any parameter that would point to stronger jet flow" in the barrier case), after hf-low issue 248.
Post hoc in the sense of `../PREREGISTRATION.md`: the 73 northerly-class episodes and their hurricane-force labels (18) already exist and
their depth / ice-sheet-pressure results were seen. The new predictor below and its direction were fixed here before any of its values exist.

**Population:** the 73 episodes with `cls == "E"` in `../results/episodes.csv` (northerly flow, gust maximum east of Cape Farewell), 20 winters 2004-05 to 2025-26.
Time per episode: the 00/12 UTC time of the episode's highest box gust (`time` of the row; one row per episode in this class).
**Outcomes:** (1) `hf` = episode reaches 71.7 kt (18 of 73), (2) peak box gust in kt (continuous).
**New predictor (the only one):** Denmark Strait pressure difference DP = MSLP(Iceland point 64.5N 18W) - MSLP(Greenland east-coast point 65.6N 37.6W), each the mean of the
3x3 nearest 0.25 deg ERA5 cells, at the episode time. Reason: the cross-strait difference is the usual geostrophic measure of how hard the flow along the coast is pushed.
Expected sign: larger DP (higher on the Iceland side) goes with stronger flow. Sea-level pressure over the ice margin is a reduction, so DP is a proxy.
**Primary test (one):** logistic of `hf` on standardised DP and standardised low depth (`LOW_DEPTH`), odds ratio per SD for DP, 95% interval by resampling whole winters (2,000 draws).
**Secondary (labelled secondary):** rank correlation of DP with peak box gust within the class, and the same logistic with ice-sheet pressure (`GH`) added, to see whether DP only repeats it.
**Decision:** "a lead" if the interval excludes 1 in the expected direction; otherwise "cannot tell". Nothing here is a forecast rule. Two tests in all; no multiplicity correction beyond reporting both.
**Power:** with 18 hurricane-force episodes only odds of about 2 per SD or more can be seen; a smaller effect will read as "cannot tell", not as no effect.
**Not tested (stays untested):** low-level stability (needs pressure-level temperature), upper-air patterns, ASCAT.
**Pull:** MSLP only, 73 times, ARCO ERA5 0.25 deg full-globe chunk per time (about 2.5 MB each, under 0.2 GB in all). No existing cache holds MSLP at these times
(`/mnt/project-files/tip-jets/` has only the composite figure; the 33.6 GB field tar is not on this mount).
