# More cyclones, or a larger share reaching hurricane force? (ERA5 proxy, pipeline A)

Plan, committed before any fit: [PREREGISTRATION.md](PREREGISTRATION.md) (commit `933e68f`).
Code: [split.py](split.py), [fdr_power.py](fdr_power.py). Results: [results/](results/)
(`summary.txt` is the readable report; `results.csv` has every number; `fdr_power.txt` has FDR and detectable effects).

Reproduce (about 7 minutes, no ERA5 pull):

    python3 research/era5/freq_split/split.py all_tracks.csv.gz <cpc_indices dir> . research/era5/freq_split/results 2000
    python3 research/era5/freq_split/fdr_power.py research/era5/freq_split/results

`all_tracks.csv.gz` is `research/era5/hf_history/results/all_tracks.csv.gz` (commit `8b109b5`, hf-low PR 12).
All numbers are **pipeline A** (800 km ocean gust index, HF-equivalent at 71.7 kt) and are a **proxy**.

## Answer

Per +1 SD of the index, Oct-Apr 2004-05 to 2025-26 (22 seasons), index lagged to days -10..-4 before genesis.
Cyclones are all 9,636 Atlantic and 10,000 Pacific pipeline A tracks that start in the window.

| basin, index | RR(HF lows) | = RR(all cyclones) | x RR(share reaching HF) | share of the effect, f |
|---|---|---|---|---|
| Atlantic, NAO | 1.096 (1.03-1.16) | 0.974 (0.96-0.99) | 1.125 (1.05-1.20) | 1.29 (1.07-1.78) |
| Pacific, PNA | 1.072 (1.01-1.14) | 1.028 (1.02-1.04) | 1.043 (0.99-1.10) | 0.61 (-0.65 to 0.79) |

Intervals are 95% season-block bootstrap. Odds ratio for the share, track-level logistic: Atlantic 1.145, Pacific 1.047.

1. **Atlantic: the share, not the count.** NAO+ does not bring more Atlantic cyclones (slightly fewer, -2.6% per SD). A larger fraction of them reach 71.7 kt (+12.5%). By the rule fixed in advance (whole bootstrap interval of f above 0.5) this is "mostly share", and the count term works against it. The prediction in the pre-registration holds.
2. **Pacific: not resolved.** PNA+ brings more cyclones (+2.8%, clear) and a somewhat larger share (+4.3%, not distinguishable from 1: permutation p = 0.24). The interval for f spans 0.5, so by the pre-registered rule it is "both contribute / not resolved". The prediction (share carries more than half) is **not supported** by the gust-based test; the count term is real, the share term is not resolved.
3. **Pacific, depth version (pre-registered secondary S4).** With the HF cut replaced by a count-matched pressure cut (965.0 hPa), the Pacific share term is +11.1% per SD (f = 0.79, 0.60-0.89), and holds in 1979-80 to 1999-2000 on its own (RR(share) 1.118, f = 0.85, 0.62-1.01). So the Pacific share effect is present when "strong" means deep, not when it means high gust. We do not know which is right; this is flagged below.
4. **The lag matters.** The archive's published effects (Atlantic 1.128, Pacific 1.217) and pipeline A's own (1.18-1.21) are same-time. With the same-time index (S1) this analysis gives 1.210 and 1.186, matching them; the lagged index gives 1.096 and 1.072, less than half of the log effect (41-48%). Part of the same-time association is the storm in its own index (the point Q1 found). The same-time split is Atlantic share-dominated (f = 1.15) and Pacific mixed (f = 0.70, 0.45-0.83).
5. **Everything else is null or secondary.** ONI: no effect on either factor in either basin (RR 1.00-1.02; effective n about 22, so only RR of about 1.1 or more is detectable). Pacific NAO: null. Atlantic PNA: share RR 0.88 (p = 0.001); exploratory, one of 60 tests.

## Power

The effect per SD that each test could detect at 80% power (2.8 x clustered SE, `fdr_power.txt`): Atlantic NAO count 1.03, share 1.10, HF 1.08; Pacific PNA count 1.02, share 1.09, HF 1.10. So the Pacific share null is not strong: effects up to about 1.09 per SD (a bit more than twice the estimate) cannot be excluded. The Pacific HF-frequency effect itself (1.072, permutation p = 0.047, BH q = 0.068, Bonferroni x2 p = 0.095) is marginal under the lag. Effective n: 22 seasons x about 5.3 (NAO) and 5.7 (PNA) independent index values per season, not the cyclone count.

## Multiplicity

60 estimates (3 per analysis, 20 analyses) are in `fdr_power.txt`, 52 distinct after removing repeated estimates; Benjamini-Hochberg q is across the 52. For the two primary analyses the q values are Atlantic count 0.010, share 0.003, HF 0.010; Pacific count 0.007, share 0.30, HF 0.068. Every secondary analysis is labelled as such in `summary.txt`.

## What this does not show

- **A "cyclone" here is any pipeline A track:** a low below 1010 hPa lasting at least 24 h with at least two fixes in a basin domain. A change in the number of weak lows counts the same as a change in strong ones. S3 splits the chain at a real-cyclone stage (min pressure at most 1000 hPa and at least 48 h): the Atlantic share effect is mostly in the last step (HF given real cyclone 1.105, z 3.1; real given all 1.019); the Pacific steps are 0.993 and 1.051, neither significant.
- **Gust and depth disagree in the Pacific** (item 3). Gust-based HF is exposed to the pre-2001 drift (not an issue for the 2004+ primary window) and to land in the gust radius; depth is exposed to NAO/PNA being defined by pressure at the same centres. The lag is the guard against the second, and removes about half the effect, so some circularity may remain in the depth version. This is a disagreement between two proxies for "strong", not a resolved finding.
- **Domain and track definition.** The NAO moves Atlantic storms north (Q1: +3 degrees), and the domain is fixed. Cyclones moving in or out of the 30-67N box would change the count term; this was not checked.
- **This is the ERA5 proxy.** It is not a count of warned hurricane-force lows. The archive's own cyclone population does not exist, so the share cannot be checked against the archive.
- **Same storms as the catalog.** The HF cyclones here are the 1,003 and 825 events already used in the frequency test; the new information is the denominator.

## Departures from the plan

None changes a pre-registered estimate. Logged because they happened after the plan was committed.

1. Added Benjamini-Hochberg FDR across all tests and detectable-effect sizes (`fdr_power.py`) at the project's request (relayed 06:50Z, before results were available). The plan's Bonferroni x2 on the two primary decompositions is also reported.
2. The CPC daily files hold `-99.000` (missing) values fused to the day field on two lines; the parser was fixed to read the fixed-width columns and treat them as missing. Found while debugging the first run, before the full run.
3. ONI was lagged to the month of day -7 (the middle of the lag window); the plan did not say how. Chosen before any ONI result was seen.
4. A 20-draw debug run of the full script was done before the 2,000-draw run, and its estimates (not its p-values) are identical. The debug run was looked at, so the primary analyses were seen once before the final run; nothing was changed in response to them.
5. Resampling loops use a Newton-Raphson Poisson solver instead of statsmodels' GLM; the estimates are the same.
