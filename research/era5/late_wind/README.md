# Why HF wind comes at or after the lowest pressure (RA-7)

ERA5 **proxy**, **pipeline A** (`research/era5/hf_history`, 800 km ocean gust index, threshold 71.7 kt). Seasons 2004-05 to 2025-26,
1,886 events (Atlantic 1,039, Pacific 847); events linked to tropical cyclones are out. Plan committed first: `PREREGISTRATION.md`
(`33039f0`). Atlantic fixes north of 60N are lower confidence (no terrain mask). Nothing was pulled.

## Answer in plain words

**Most of the "29%" is not wind arriving after the low.** In the proxy, 29.3% of events reach HF at or after the fix of lowest pressure
(Atlantic 34.0%, Pacific 23.5%). But 18.4% reach HF on the *same* 6-hourly fix as the minimum (63% of the late storms), and 78% of the
late storms have onset pressure within 2 hPa of their minimum. Strictly after by at least 6 h is 10.9%, by at least 12 h 4.7%.

**What sets the late ones apart** (odds ratios for being late, with 95% intervals clustered by season; all in `results/tests.csv`):

| factor | late share with / without | odds ratio | share of all late storms with it |
|---|---|---|---|
| Marginal storm (peak gust within 5 kt of the threshold), T3 | 41.0% / 19.5% | 2.96 [2.43, 3.62] | 64% |
| Atlantic, gust maximum at or north of 60N, T2 | 52.7% / 26.9% | 3.10 [2.24, 4.31] | 43% of Atlantic late storms |
| Slow translation, T6 (per SD, 13 kt) | | 0.49 [0.44, 0.55] | |
| Not an explosive deepener (< 1 Bergeron), T5 | 44.2% / 22.2% | 0.35 [0.28, 0.44] for bombs | 45% |
| Gust maximum more than 400 km from the centre, T1 | 34.7% / 27.6% | 1.41 [1.16, 1.71] | 28% |
| Gust maximum in the rear half, T4 | 24.2% / 41.2% | 0.45 [0.36, 0.56] (**opposite** to what I expected) | |
| Pacific rather than Atlantic, T7 | 23.5% / 34.0% | 0.59 [0.48, 0.72] | |

Reading it: late storms are mostly **short, marginal, slow-moving, non-explosive events** that cross the threshold once, at or just
after the lowest pressure (61% have a single HF fix, against 24% of early storms; 80% of late storms have their onset on the peak-gust
fix). The agenda's two predictions both hold: T1 (far gust maximum) and T2 (north of 60N in the Atlantic) are positive and pass FDR.
But they are not equal. T2 is large and survives every adjustment (joint OR 1.82 [1.09, 3.03], Atlantic only; 2.85 with Greenland- and
Iceland-terrain-flagged fixes removed). T1 is modest, falls to OR 0.81 [0.64, 1.04] once marginal, rear-half, deepening and speed are in the model, and is
not significant in events with two or more HF fixes (OR 1.30, p 0.15). Rear-half gust maxima, which I expected to mark post-occlusion wind, are
*less* often late: the rear-right maximum belongs to deepening storms (PR 27). So there is no sign that post-occlusion wind is what makes a storm late.

**How much is explained.** The five storm factors with basin and month separate late from early storms with leave-one-season-out
AUC **0.76** (basin and month alone 0.57; Brier skill 0.14). By the plan's rule (AUC at least 0.70) that is "largely identified by these
factors", not an unexplained class. In the Atlantic, adding the 60N factor gives AUC 0.77.

**Two kinds, with the evidence only partly separating them.** (1) A threshold effect: weak or short-lived storms whose gust only just
reaches 71.7 kt do so near their pressure minimum, with a tie or one step either way. (2) A high-latitude Atlantic effect, where wind persists
or arrives at a low's mature stage near Greenland and Iceland. I did not test its physics (barrier flow, a Greenland high, orography) here; PR 29 and PR 58 are
the places for that.

## All tests

Primary family T1-T7: 7 of 7 pass BH q < 0.05, and 6 are in the expected direction (T4 is significant in the opposite direction; T6 had no
expected sign). Exploratory family S1-S9 (environment at the last 00/12 UTC fix before onset, per SD): 9 of 9 pass q < 0.05 (16 of 16 across both
families). The S results are mostly a restatement of life-cycle stage: more thermal wind aloft (VTU OR 2.50, VTL 2.28), less jet and ascent
(div300 0.48, vadv500 0.46), less Eady growth (0.55), smaller storm asymmetry B (0.45) all mean "older, filling storm", which late storms are by
definition. They are not independent evidence for a mechanism.

**Power.** For T1-T7 the power at OR 1.5 is 0.90, 0.75 (T2), 0.96, 0.92, 0.90, 1.00 and 0.97; the minimum detectable OR at 80% is 1.4 to 1.6
(1.2 for T6). T2 is the least powered, and still significant. No test here is a null, so the well-powered-no rule was not needed.

## Post hoc (not in the plan; both versions shown)

- `posthoc.py`: the joint model with intervals (above), and the pressure gap at onset. Not pre-registered, labelled as such.
- Logged deviations are in `PREREGISTRATION.md`: onset time taken from the catalog tracks, Jun-Aug folded into Sep-Nov, T4's sign.
- Heterogeneity: every sign agrees between Atlantic and Pacific and between 2004-14 and 2015-25, except that T1 is weaker in the Pacific
  (1.23 [0.96, 1.56]). T2 is Atlantic-only by design.

## Caveats

- ERA5 at 31 km under-resolves the gust peak and the sampling is 6-hourly, so "same fix" is a tie the data cannot break. Where pressure is
  flat, which is the case for non-deepening storms, the time of the minimum is poorly defined, which inflates late counts for exactly the storms
  T3, T5 and T6 flag. I could not remove this with the data on hand.
- Atlantic events north of 60N are lower confidence (no terrain mask). T2 is a statement about those events; the sensitivity run with the
  flagged fixes removed leaves it at OR 2.85.
- Onset and minimum pressure are ERA5 quantities. Nothing here says the archive's own HF fixes (OPC's judgement) follow the same timing.
- No pre-2001 check: onset is a gust threshold crossing and decision 1 rules out gust-based use before 2004-05.

## Held-out looks

One look at 2015-25 (all 22 seasons used once, no model trained on one set and scored on another; 11th by the coordinator's count) and none at pre-2001 seasons.
Logged in `../hemispheric/results/heldout_looks.log`.

## Verification

See `verify/VERIFICATION.md` (a fresh Sonnet agent recomputed the numbers it lists from the committed inputs without my code).

## Reproduce

`python3 late_wind.py` (about 8 minutes, mostly the power simulation), `python3 posthoc.py`, `python3 figure.py`. Needs numpy, pandas, scipy,
matplotlib; reads committed files only. Figure: `results/late_wind.png`.
