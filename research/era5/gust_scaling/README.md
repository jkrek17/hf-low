# Does the gust in hurricane-force fixes follow gradient-wind scaling? (RA-9 first pass; ERA5 proxy, pipeline A)

Pipeline A = `research/era5/hf_history`. Every number is an ERA5 **proxy** (31 km gust and MSLP, pipeline A's 800 km
gust index g800), 5,934 HF-strength fixes on 1,988 storms, seasons 2004-05 to 2025-26. Transitioning tropical cyclones
are in. No pull: committed tables only. Plan committed first (`3acbd3e`, `PREREGISTRATION.md`); deviations are logged
there. **No held-out look was spent** (an estimation on all 22 seasons, nothing frozen or scored).

## Answer in plain words

**Partly.** Among storms that reach hurricane-force gust, the gust rises with depth, falls with latitude and falls with
the radius of the maximum gust, all with the sign gradient-wind balance predicts, in both basins. But two of the three
slopes are much weaker than the idealised balance gives, and the size result depends on how size is measured.
Indices do not matter at fixed depth, size and latitude, and ENSO does not move the size of these storms.

Slopes are elasticities of the underlying (before selection) gust: the change in ln gust per unit ln of the
predictor. 95% intervals from a season-block bootstrap (22 seasons, 2,000 draws).

| | Atlantic | Pacific | Predicted by gradient-wind balance |
|---|---|---|---|
| Depth (ring-mean pressure minus centre) | **+0.34** [0.28, 0.40] | **+0.26** [0.21, 0.31] | +0.55 to +0.72 |
| Radius of maximum gust | **-0.053** [-0.071, -0.037] | **-0.043** [-0.064, -0.022] | -0.11 to -0.45 |
| Latitude (as sin lat, the Coriolis term) | **-0.16** [-0.25, -0.07] | **-0.48** [-0.62, -0.38] | -0.11 to -0.45 |

- **Depth: follows in sign, only about half as strong.** A storm 10% deeper (about 3.6 hPa) has a gust about 3.3%
  (Atlantic) or 2.5% (Pacific) higher, roughly 2.5 and 1.9 kt at 77 kt. The interval excludes the predicted range.
  Plausible reasons, not separated here: Dp is approximate (ring mean at the storm's deepest fix, applied to all its
  fixes, which attenuates a slope), ERA5 gust is not gradient wind, and the same storm's gust is capped by ERA5's
  resolution.
- **Size: negative as predicted, but small and not robust to the definition.** Radius of maximum gust (the
  pre-registered size): -0.05 per ln unit, so 10% more radius is about 0.5% less gust. The sign **flips** (+0.085
  Atlantic, +0.141 Pacific) when size is the equivalent radius of the area with gust of 48 kt or more, which is built
  from the same gust field and grows with the maximum gust by construction. So the data cannot say that bigger storms
  have weaker peak gust at fixed depth; they say the radius of the peak does.
- **Latitude: follows, inside the predicted range, in both basins.** Moving from 45 N to 55 N lowers the
  gust by about 2.4% (Atlantic) and 7.1% (Pacific), 1.9 and 5.5 kt. Using raw depth (1013 hPa minus centre) instead of
  ring-relative depth makes the Atlantic latitude slope 2.4 times steeper (-0.38) and the Pacific one almost unchanged
  (-0.53), so PR 52's "1.0 kt less gust per degree at fixed depth" is mostly background pressure in the Atlantic and
  not in the Pacific.
- **ONI does not change the size of these storms.** Per SD of ONI the radius of maximum gust changes by +0.5%
  [-2.7, +3.2] (Atlantic) and +1.9% [-1.7, +4.7] (Pacific). The smallest effect the test could detect with 80% power
  (2.8 x SE, an approximation) is 4.0% and 4.5%, below the pre-registered 5%: **a well-powered null for a 5% ENSO
  effect on radius**, so the unresolved ONI gust-versus-depth gap of PR 52 is not about the radius of the peak.
  Smaller changes are not excluded, and size is observed only in storms that reached HF.
- **ONI and PNA do not add gust at fixed depth, size and latitude (cannot tell below about 2%).** Direct slopes
  -0.7%, +0.2% (Atlantic) and -0.6%, +0.6% (Pacific) per SD, none near significance; the minimum detectable effect is
  1.8% to 2.5% per SD (about 1.4 to 2 kt), above the pre-registered 1% needed to call it "no". Adding radius to the model
  changes the ONI slope by 0.0006 to 0.0011 (intervals include 0), so size is not a pathway.
- **One lead, not a finding: a +PNA Pacific has a larger radius of maximum gust** (+4.4% per SD [0.5, 8.1], p 0.026,
  q 0.052 across the 14 primary tests, so it just misses). With the area-based size it is +4.0% (p 0.007, secondary,
  q 0.014 across all 193 secondary and primary rows, which are not independent). A +PNA Pacific has a deeper
  background (Aleutian low; `tele_intensity`, PR 11), which may broaden the storm. Needs new seasons to confirm.

## What this sample can and cannot show

Every fix has g800 >= 71.7 kt by construction (selection on the outcome), so a plain regression is biased towards zero.
OLS of ln g800 on the same covariates gives depth +0.083 / +0.070 (Atlantic / Pacific), a quarter of the truncated
estimates above: **the naive regression would have said the scaling is four times weaker.** The estimator is a
truncated-normal MLE of ln g800 (known truncation point), which recovers population slopes only if the tail of ln g800
is normal. It is: the probability-integral-transform KS distance is 0.010 (Atlantic) and 0.016 (Pacific) against a 5%
critical value of 0.024 and 0.026, while a generalised Pareto tail fails (0.027, 0.033; post hoc check, `diag.py`).
Treat the slopes as model-based: the profile likelihood over sigma (`results/profile_sigma.csv`) shows the fitted slope
rises with sigma along a ridge, with a clear maximum (log-likelihood falls by 14 units at sigma 0.16 and 43 at 0.20 in the
Atlantic).

It cannot say how storms that never reach HF scale, whether the scaling holds below 71.7 kt, or whether an index
changes the *probability* of reaching HF (PR 85: it does, through conversion of deepening storms). Size exists only here
(agenda RA-9: all-cyclone size needs a pull, size it first). RMG is read off the same gust field as the outcome.

## Tests and multiplicity

Primary family: 14 tests (per basin five gust slopes and two size slopes). **6 of 14 pass q < 0.05** (depth, radius,
latitude in both basins, all q <= 0.0008); the other 8 are the ONI and PNA slopes, the
smallest q being 0.052 (Pacific PNA to radius). Secondary S1-S10 are reported in `results/results.csv` with their own
q across all 193 rows; they do not change any primary sign. Sensitivities (Atlantic / Pacific, depth, radius, latitude
slopes): one deepest fix per storm (+0.42, -0.058, -0.14 / +0.34, -0.031, -0.49); drop terrain, coast and Atlantic
fixes north of 60 N (+0.33, -0.062, -0.21 / +0.27, -0.042, -0.48); climatological depth (+0.35 / +0.28); month effects
(+0.32 / +0.28); by stage (all signs kept; Pacific filling storms steeper); same-time PNA instead of lagged (no change);
pooled with a basin term (+0.30, -0.050, -0.30). One isolated nominal result: Atlantic ONI on gust in the one-fix-per-storm
sample (-1.4% per SD, q 0.039 across all rows) is not in the primary family and the primary Atlantic ONI slope is -0.7%
(p 0.34); read it as noise.

## Files

| File | What |
|---|---|
| `PREREGISTRATION.md` | Plan, committed first (`3acbd3e`) and the logged deviations |
| `scaling.py` | All fits and the season-block bootstrap. `python3 scaling.py 2000` (about 12 min, 4 cores; needs scipy) |
| `diag.py` | Tail check, six-start check and profile over sigma |
| `results/results.csv` | Every test: estimate, SE, interval, p, q (primary family and all rows), minimum detectable effect |
| `results/descriptive.csv`, `meta.json`, `tail_diagnostic.csv`, `trunc_normal_multistart.csv`, `profile_sigma.csv`, `run.log` | Supporting |

## Verification

See `VERIFICATION.md`: every headline slope, sigma, count, correlation and KS value was recomputed to 4 decimals by a fresh agent; the minimum detectable effects, secondary models, conversions and the 2,000-draw intervals were not.

## Deviations from the plan

See the log in `PREREGISTRATION.md`: an optimiser fix (estimator unchanged), the added GPD and tail diagnostics, and
implementation notes. The GPD (S10) slopes are within-HF conditional effects, comparable to the OLS contrast, not to the
primary slopes; the plan wrongly described them as a second estimate of the same quantity.
