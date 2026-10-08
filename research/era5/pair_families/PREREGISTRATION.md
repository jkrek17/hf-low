# RA-12 plan: are close HF pairs parent-daughter, or a shared environment?

Written before the test statistic below was computed. Pipeline A (`research/era5/hf_history`), an ERA5 **proxy**; transitioning tropical cyclones are in (not separated). Follows hf-low PR 37 (`research/era5/clustering`).

## What had been looked at before this plan (disclosure)

- `count_pairs.py` printed the number of candidate pairs (definition below): Atlantic 160, Pacific 60 (2004-05 to 2025-26, 1,104 and 896 HF tracks). Of the daughters with observed genesis (first fix at 1000 hPa or above, north of 21N, within 12 h of the tracker start): 15 Atlantic, 9 Pacific.
- `power.py 40` (40 permutations) then printed, **unintentionally including the rear-sector counts**: observed rear-sector pairs Atlantic 70 of 160, Pacific 18 of 60 (observed-genesis daughters: 10 of 15, 5 of 9); permutation-null means of the all-pairs count 197.8 (Atlantic) and 98.6 (Pacific), of the rear fraction about 0.41 and 0.49. So the **sign of the rear-sector result was seen before this plan was written**. The test below is therefore not a clean pre-registered test of the rear-sector hypothesis; it is written to fix definitions, the null and the decision rule before the 2,000-permutation run, and every result is labelled "orientation seen first". No other definition is tried.
- Not looked at: any p-value, interval, conditioning on a teleconnection index, or any other sector/distance/time definition.

## Question and answer format

Jason's words: for close pairs of HF lows, is the second a parent-daughter development (secondary cyclogenesis on the first's front) or do both share the same large-scale environment? Answer: parent-daughter dominates / environment only / can't tell, with the size of the most the data allow.

## Definitions (frozen)

- HF track: pipeline A track with `gust800_kt >= 71.7`, seasons 2004-05 to 2025-26 (decision 1), per basin (Atlantic 1,104, Pacific 896).
- Genesis proxy: first fix in `intensity/results/fixes_2004.csv.gz` (00/12 UTC fixes below 1010 hPa). It is the tracker's first detection, not true genesis.
- Candidate pair (P parent, D daughter), same basin: D's first fix is later than P's first fix and no later than P's last fix (P alive); |time of D's first fix - time of P's peak gust| <= 48 h; D's first fix lies within 1,500 km of P's position at that time (linear interpolation of P's fixes).
- Rear sector: bearing from P to D minus P's heading (from the two fixes bracketing that time) between 120 and 240 degrees, i.e. within 60 degrees of directly behind. This stands for the trailing-front sector. The 240-degree complement is "other".
- Observed-genesis daughters (secondary subset): D's `gen_obs` true and first fix within 12 h of the tracker start.

## Null: shared environment

Whole tracks are shifted in time by permuting peak times among the HF tracks of the same basin, season and calendar month (2,000 permutations, seed 20261008). This keeps each storm's shape and position, the seasonal cycle and the basin's storm-track geometry, and removes only the timing link between storms. It is the same logic as PR 37's Knox test. Known limitation, seen in the 40-permutation run: the null all-pairs count is larger than the observed one (the tracker, and perhaps the weather, keeps new storms from starting next to live ones), so **counts of pairs are not compared with the null**. Only the orientation of the pairs that exist is, which does not rely on the count being right.

## Tests

Primary family (2 tests, Benjamini-Hochberg, q < 0.05): rear fraction f = rear pairs / candidate pairs, Atlantic and Pacific, one-sided (rear fraction above the null mean), p from the permutation distribution of f.

Secondary (reported, not decision-bearing): the same fraction for observed-genesis daughters only (Atlantic n = 15, Pacific n = 9); observed/null ratio of the rear-sector count (biased by the count deficit above, shown only to display the deficit); two-sided p for every test.

Conditional (run only if a primary test has q < 0.05): the same fraction with the permutation stratified also by the tercile of the lagged (days -10 to -4 before D's first fix) NAO (Atlantic) or PNA (Pacific), to see how much a shared pattern accounts for.

## Power and minimum detectable effect

Planted-daughter simulation (`power.py`): observed rear count = Binomial(number of candidate pairs, null mean f) plus k extra rear-sector pairs, compared with the permutation distribution of f. Report the k with 80% power at one-sided 0.025 and what k means as a share of PR 37's spatial excess (42 pairs Atlantic, 24 Pacific: observed minus expected Knox pairs) and of HF tracks.

## Decision rules

- **Parent-daughter, at least a part**: primary q < 0.05 (rear fraction above null).
- **Not dominant (well-powered)**: q >= 0.05 and the upper 95% bound of the extra rear pairs (observed rear count minus null expectation at the observed number of pairs) is below the k with 80% power, and below half of PR 37's excess.
- **Can't tell**: anything else, in particular any mixture smaller than the stated k, and anything that rests on the 15 and 9 observed-genesis pairs.
- A hand-off to an environment-only conclusion requires a positive result of the conditional test; a null on the primary test alone is "no evidence of daughters", not "shared environment shown".

## What this study cannot show

First detection is not genesis, 16 percent of tracks start more than 12 h before their first 00/12 fix, and many "daughters" are old storms entering the domain. A true secondary on a trailing front can start while still attached to the parent's tracker track (then it is not a second track at all). The design cannot separate frontal waves from other cyclogenesis behind a storm.

## Looks at the held-out block

None. There is no fitted model and nothing scored on 2015-25; no entry is written to the looks log.

## Deviations (post hoc)

(none yet)
