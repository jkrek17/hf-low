# Upper-level bomb-environment diagrams and 24 h gale / storm / HF outlook (ERA5 proxy, pipeline A)

**This is a proxy.** Everything here is ERA5 on pipeline A's tracks (`../hf_history`), 2004-05 to 2025-26, tropical
cyclones in. Outcomes are ERA5's own gust index, not the archive. Plan committed first: `PREREGISTRATION.md`
(`7c0669a`). Predictors are the existing `../intensity` environment table; no new pull.

Files: `score.py` (leave-one-season-out scoring), `diagrams.py` (binned observed rates), `results/skill.txt`,
`results/family.csv` (the pre-registered family), `results/reliability.csv`, `results/diagram_cells.csv`,
`results/diagrams.txt`, `results/oof_outlook_full.csv.gz` (held-out-season outlook probabilities per fix).

## What the outlook is
For each 00/12 UTC fix, the chance that the peak gust index (max ocean gust within 800 km) over the next 24 h falls
in: below gale (<38 kt), gale (38 to <54), storm (54 to <71.7), hurricane force (>=71.7). Cutoffs are nominal
proxy classes: only 71.7 kt is tied to the archive, which lists HF lows only. The outcome uses the t+12 and t+24 h
fixes (the 6-hourly fixes are not in the committed table), so its HF class is not the same label as `hf24` in `../intensity`.
A track that ends inside 24 h is scored on the fixes it has.

## Results (leave-one-season-out, 22 seasons, season-bootstrap 90 % intervals)

| Score vs basin-month climatology | state only | + jet and trough (A) | + jet and baroclinicity (B) | full (Storm model predictors) |
|---|---|---|---|---|
| Outlook RPSS | +0.368 | +0.380 | +0.376 | +0.399 |
| Rapid-deepening BSS | +0.161 | +0.212 | +0.234 | +0.285 |

Pre-registered family, gain over state (BH q within the 6, all q = 0.0005 = the bootstrap floor):

| Target | A | B | full |
|---|---|---|---|
| Outlook RPSS | +0.0115 [+0.0106, +0.0123] | +0.0084 [+0.0079, +0.0089] | +0.0306 [+0.0295, +0.0318] |
| Rapid-deepening BSS | +0.0504 [+0.0434, +0.0571] | +0.0733 [+0.0667, +0.0798] | +0.1234 [+0.1154, +0.1316] |

All six: yes. Gains are largest for the rarest class: hurricane-force Brier-skill gain over state +0.067, storm
+0.036, gale +0.030, below gale +0.023. The full model beats A by 0.019 and B by 0.022 RPSS, so the two upper-level
axes capture part of the environment's value but not most of it. Atlantic and Pacific gains are alike (+0.031,
+0.031). The reliability slope is 1.00 to 1.13 for the four classes: probabilities are close to calibrated.
The current gust term matters: dropping it lowers the outlook RPSS by 0.044. Most of the skill (0.37 of 0.40) is
the storm's own state, as in PR 12; the environment adds the +0.03.

## Reading it
- The fraction of fixes where the most likely class is correct is shown in `skill.txt`; hurricane-force is
  under-called as storm (2,351 of 3,880 HF-outcome fixes get "storm" as the most likely class), which is expected
  for a 2.4 % class.
- The diagrams show the observed rapid-deepening share in each cell of 6 x 6 quantile bins (all fixes pooled, not
  conditioned on the storm's state). Highest cells, strong jet with positive vorticity advection or high Eady growth,
  reach 21 % and 25 % against 3.4 % overall; low cells are near 0. A cell is a marginal rate: it mixes storms that were already deep with
  storms that were not, so the model, not the cell, is the skill claim.

## Limits
- Gale and storm cutoffs are not calibrated on observations. The 12-hourly outcome samples miss peaks between samples.
- Passing FDR at this n says the gains are not zero, not that they are large; RPSS gain over state is +0.01 to +0.03.
- Interval widths treat seasons as independent; fixes within a storm are not independent (the interval is the season-level one).
- Perfect-prognosis setup (ERA5 predictors and outcomes); not tested on GFS.
- Looks at these data: one scoring run for the pre-registered family. The old `../intensity` scores on the same fixes are PRs 7 and 12 (not independent of this).

## Verification
A fresh Sonnet agent, writing its own scoring code from the committed tables, recomputed: class frequencies
(0.5157 / 0.3120 / 0.1479 / 0.0243), outlook RPSS state +0.3680 and full +0.3987, rapid-deepening n 100,025,
base rate 0.0343 and BSS state +0.1612, full +0.2846. All match.
Not independently checked: the A and B rows, the bootstrap intervals, p and q values, per-class and per-basin
gains, reliability slopes, the confusion table, and the diagram cells.

## Four regions of the vorticity-advection x Eady plane (added after the first run, not pre-registered)
`diagrams.py` also writes `results/regions.csv`: all fixes split at the upper third of 500 hPa vorticity advection
(0.372) and of Eady growth (0.760 per day) into quiet, baroclinic only, forcing only and coupled. Coupled is 19 % of
fixes, deepens rapidly in 10.6 % (base 3.4 %) and holds 72 % of the rapid-deepening fixes; quiet is 52 % of fixes at 0.4 %;
baroclinic only 14 % at 4.6 %; forcing only 14 % at 1.2 %. P(HF next 24 h): 8.5 / 2.1 / 1.3 / 0.6 %. The cuts are
empirical, not physical thresholds; fixes are not independent (a storm changes region within a day). Descriptive, in-sample.
Third diagram `trough_baroclinic` added to `results/diagram_cells.csv`. Not independently checked until the verifier note below is added.

## Pairing the regions with the models (post hoc)
`pair_regions.py` (output `results/pair_regions.txt`): adding the three region dummies changes the full model's outlook RPSS by 0.0000 [-0.0002, +0.0001] and rapid-deepening BSS by +0.0011 [+0.0004, +0.0018]; state+Hart by +0.0047 and +0.0012; state alone by +0.0098 and +0.0366. The regions are a coarse copy of variables the full model already uses. A fresh agent recomputed the model scores and gains (not the intervals).
