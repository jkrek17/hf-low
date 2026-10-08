# Do the share results hold when only sustained HF lows count? (RA-22; ERA5 proxy, pipeline A, and the archive)

Plan, committed before any 2-fix or 3-fix outcome was run: [PREREGISTRATION.md](PREREGISTRATION.md) (`38e5cf6`, clarified in `c365311`;
corrections made afterwards are listed under "Deviations and corrections" below). Code: `labels.py`, `run_pr14.py`, `run_channels.py`,
`contrast.py`, `report.py`, `posthoc_intensity.py`. Results: [results/](results/) (`summary.txt` is the readable report; `headline.csv`,
`contrasts.csv`, `paired.csv`, `secondary.csv`, `posthoc_intensity_*.csv`, per-study folders).

Reproduce (no ERA5 pull; about 25 minutes on 4 cores; the derived track tables go to the ignored `work/`):

    python3 research/era5/sustained_hf/labels.py research/era5/sustained_hf/work
    # then, for k in 1 2 3 (and run{2,3}, gcut{2,3} for the secondary checks), with W=research/era5/sustained_hf/work:
    python3 research/era5/sustained_hf/run_pr14.py $W/all_tracks_k$k.csv.gz results/pr14_k$k.csv 2000
    python3 research/era5/sustained_hf/run_channels.py $W/all_tracks_k$k.csv.gz results/channels_k$k 2000 2000
    python3 research/era5/nao_share_barrier/barrier_share.py $W/all_tracks_k$k.csv.gz <cpc_indices> . \
        research/era5/nao_share_barrier/gh_daily.csv research/era5/highlat/gustloc_fixes.csv results/pr63_k$k 2000
    python3 research/era5/sustained_hf/contrast.py {daily|weekly|archive} results/contrast_<design>_withm.csv 2000 2000
    python3 research/era5/sustained_hf/posthoc_intensity.py {weekly|daily} results/posthoc_intensity_<design>.csv 2000 2000
    python3 research/era5/sustained_hf/report.py research/era5/sustained_hf/results

"HF_k" means pipeline A peak gust at least 71.7 kt **and** at least k in-domain 6-hourly fixes at or above 71.7 kt (archive: `hfN >= k`).
A track that fails HF_k stays a cyclone in every denominator. The originals' code runs unchanged; only which tracks count as HF changes.
At k = 1 every original number is reproduced exactly (PR 14, 41, 63, 64; see "Controls").

## Answer

**Yes, every headline result holds for sustained HF lows, and none gets weaker.** All 26 pre-registered headline tests keep their sign and
pass FDR at k = 2 and k = 3 (largest q 0.029). The point estimates are larger at 3 fixes than at 1 in all 10 PR 14 / 41 / 64 ratio comparisons
(k = 1 to k = 3), by up to 9% of the ratio (four of the twelve paired differences have p < 0.05, none passes FDR). The pattern does not favour brief lows. Whether it favours sustained ones over brief ones is
**not established**: all six sustained-versus-brief contrasts point the same way (sustained larger, by 6% to 18% in the ratio), but none passes FDR
(smallest q 0.108) and the study could detect a contrast of only about 19% to 31%. So: **no clear preference, leaning sustained.**

Per +1 SD of the index. NAO for the Atlantic and PNA for the Pacific in PR 14 and PR 63; the leave-one-season-out hemispheric pattern index of PR 41 in PR 41 and
PR 64. Oct-Apr 2004-05 to 2025-26, 22 seasons. 95% season-block bootstrap intervals in `headline.csv`.

| result | basin | 1 fix (original) | 2 fixes | 3 fixes | holds? (q) |
|---|---|---|---|---|---|
| PR 14 RR(HF lows), NAO / PNA | Atlantic | 1.096 | 1.109 | 1.137 | yes (0.026, 0.025) |
| | Pacific | 1.072 | 1.128 | 1.152 | yes (0.007, 0.007) |
| PR 14 RR(share reaching HF) | Atlantic | 1.125 | 1.139 | 1.168 | yes (0.009, 0.009) |
| | Pacific | 1.043 (p 0.24, unresolved) | 1.098 | 1.121 | yes (0.029, 0.025) |
| PR 14 RR(all cyclones), unchanged | Atl / Pac | 0.974 / 1.028 | same | same | |
| PR 64 RR(HF lows), pattern | Atlantic | 1.291 | 1.316 | 1.385 | yes (0.0009) |
| | Pacific | 1.165 | 1.219 | 1.264 | yes (0.0009) |
| PR 64 RR(share reaching HF) | Atlantic | 1.271 | 1.295 | 1.364 | yes (0.0009) |
| | Pacific | 1.129 | 1.180 | 1.224 | yes (0.0009) |
| PR 64 share channel, part of the change | Atlantic | 92% (82-100%) | 92% (83-100%) | 93% (85-100%) | |
| | Pacific | 81% (72-90%) | 85% (78-93%) | 88% (82-95%) | |
| PR 41 archive RR(HF lows), pattern | Atlantic | 1.236 | 1.268 | 1.351 | yes (0.0009) |
| | Pacific | 1.193 | 1.215 | 1.275 | yes (0.0009) |
| PR 63 P1 share without terrain-type fixes | Atlantic | 1.141 | 1.156 | 1.174 | yes (0.005, 0.006) |
| PR 63 P3 mediation, NAO share given the Greenland high / NAO share | Atlantic | 0.923 | 0.881 | 0.890 | yes (0.002, 0.003) |
| PR 63 NAO share given the Greenland high | Atlantic | 1.039 (0.93-1.15) | 1.004 (0.88-1.14) | 1.039 (0.90-1.19) | |
| PR 63 F4 Greenland high on share, per SD | Atlantic | 0.870 | 0.832 | 0.824 | yes (0.0009) |

- **PR 14 (more cyclones or a larger share).** Atlantic: the share still carries the NAO effect (f 1.29, 1.26, 1.21 at 1, 2, 3 fixes). Pacific: at 1 fix the share part
  was unresolved (f 0.61, interval -0.65 to 0.80); at 2 and 3 fixes it is clear (RR(share) 1.098 and 1.121, f 0.77 and 0.81, intervals 0.26-0.89 and 0.23-0.91). The Pacific
  single-fix lows do not respond to the PNA (RR 0.976, 0.88-1.10), which is why the share term was weak. **Part of the Pacific gust-versus-depth disagreement of PR 14 and
  PR 52 is brief gust spikes.** (Descriptive, from the pre-registered brief class; the Pacific depth version has no duration and was not repeated.)
- **PR 64 (channels).** The share channel carries the effect in both basins at every k, and the Pacific share channel rises from 81% to 88%. The cyclone-count terms do not depend on k.
- **PR 41 (archive).** The hemispheric pattern's effect on archive HF counts is not a single-fix artefact: RR per SD 1.236 to 1.351 (Atlantic) and 1.193 to 1.275 (Pacific).
  The index was fitted to the k = 1 archive counts and not refitted (the fields are not here), so this is a transfer; the rise with k is therefore, if anything, conservative.
- **PR 63.** Barrier-type fixes still carry none of the NAO share effect (P2, with-terrain-removed over original: 1.014, 1.015, 1.005; unresolved below about 4%, descriptive), and the
  Greenland-high mediation gets, if anything, stronger (P3 0.923, 0.881, 0.890; NAO share given the high 1.039, 1.004, 1.039). The Greenland high lowers the storm HF share at every k
  (0.870, 0.832, 0.824 per SD).

## Sustained versus brief (pre-registered contrast, family B, BH over six)

Ratio RR(sustained, at least 3 fixes) / RR(brief, exactly 1 fix) per SD, with the three classes' own RR.

| design | basin | brief | middle (2 fixes) | sustained | sustained / brief (95% bootstrap) | p / q | detectable at 80% |
|---|---|---|---|---|---|---|---|
| PR 14, daily NAO / PNA | Atlantic | 1.071 | 1.057 | 1.137 | 1.062 (0.88-1.27) | 0.47 / 0.47 | 1.31 |
| | Pacific | 0.976 | 1.076 | 1.152 | 1.180 (1.00-1.36) | 0.036 / 0.108 | 1.24 |
| PR 64, weekly pattern | Atlantic | 1.211 | 1.200 | 1.385 | 1.144 (0.99-1.34) | 0.067 / 0.113 | 1.23 |
| | Pacific | 1.075 | 1.118 | 1.264 | 1.176 (1.05-1.36) | 0.029 / 0.108 | 1.19 |
| PR 41, archive weekly pattern | Atlantic | 1.177 | 1.118 | 1.351 | 1.148 (0.98-1.33) | 0.076 / 0.113 | 1.25 |
| | Pacific | 1.151 | 1.083 | 1.275 | 1.108 (0.97-1.33) | 0.20 / 0.24 | 1.26 |

Overall rule (written beforehand): favours sustained needs at least 4 of 6 cells with q < 0.05 and none favouring brief. Result: 0 of 6, so **no clear preference**; all
six point estimates are above 1, so it **leans sustained**. The contrast is a comparison of two noisy classes (343 brief against 439 sustained Atlantic tracks), which is why it
cannot see anything under about 20%. The middle class (exactly 2 fixes, S3) is not monotone: it sits below the brief class in 4 of 6 cells.

Paired difference RR(HF_3)/RR(HF_1), family C (descriptive, BH over 12): 1.04 and 1.07 (PR 14 daily), 1.07 and 1.09 (PR 64), 1.09 and 1.07 (PR 41 archive), Atlantic and Pacific; four
have p < 0.05 and none q < 0.05 (smallest q 0.066). Full table: `paired.csv`.

## Duration or intensity? (secondary S2 and a post hoc check)

Sustained lows have higher peak gust by construction. S2 takes the same number of tracks per basin as HF_k but chooses them by peak gust alone (`gcut`). It gives about the same
ratios as duration: PR 64 RR(HF), 3 fixes, 1.385 / 1.264 against 1.321 / 1.311 for the gust cut (Atlantic / Pacific); PR 14 1.137 / 1.152 against 1.127 / 1.151.
**Post hoc (not in the plan):** the paired difference, duration over gust cut, is 1.049 (0.96-1.14) and 0.964 (0.90-1.03) in PR 64 at 3 fixes and 1.009 and 1.001 in PR 14; at 2 fixes it is
1.014, 0.995 (PR 64) and 0.992, 1.011 (PR 14); every interval contains 1. So the rise of the pattern effect with k is **not distinguishable from an intensity gradient**: a
higher peak-gust cut does as well as requiring a longer stay at HF. Rule 5 of the plan applies: do not describe the rise as a duration effect.

S1 (HF_k by the longest run of consecutive HF fixes) gives the same answer: PR 64 RR(HF) at run >= 3: 1.390 and 1.279; PR 14 1.153 and 1.125 (p of the Pacific PR 14 share 0.094).

## Power

Detectable RR at 80% power (exp of 2.8 x clustered SE) is in `headline.csv` for every cell: 1.10-1.19 at k = 2 to 3 for the PR 64 and PR 41 tests, 1.12-1.20 for PR 14 and PR 63 (k = 1 was 1.08-1.10). The SE of log RR(HF) from 1 to 3 fixes goes 0.029 to 0.060 (PR 14 Atlantic), 0.033 to 0.055 (PR 14 Pacific), 0.025 to 0.040 and 0.037 to 0.054 (PR 64),
0.035 to 0.036 and 0.043 to 0.061 (PR 41 archive): 1.0 to 2.1 times, against the planned 1.4 to 1.5 (sqrt of the count ratio); the daily design widens most. The cell verdicts are "holds" because the effects are larger than the original, not because intervals stayed tight. The
sustained-versus-brief contrast is the under-powered result.

## Controls (k = 1)

The unchanged originals were run on the k = 1 track table first and reproduce: PR 14 Atlantic NAO RR(HF) 1.096, RR(all) 0.974, RR(share) 1.125, f 1.29; Pacific PNA 1.072, 1.028, 1.043, 0.61.
PR 64 Atlantic T1 1.016, T2 1.271, T3 1.291, f 0.94, share channel 92.1%; Pacific 1.033, 1.129, 1.165, 0.79, 80.9%. PR 63 share 1.125, P1 1.141, NAO share given the high 1.039, F4 0.870.
PR 41 archive RR per SD 1.236 and 1.193 (PR 64 S5). `contrast.py` reproduces these point estimates before use.

## Looks, spent and logged

One look at each of the 22 gust-era seasons for pipeline A (the PR 14, 63 and 64 setups) and one at the archive's 2004-05 to 2025-26 weekly counts (the PR 41 setup), logged as look 7 at the 2015-16 to 2025-26
seasons in `hemispheric/results/heldout_looks.log`. No model was refitted and no held-out split scored (leave-one-season-out index throughout). Pre-2004 seasons were not used. Each of k = 2 and 3 was run once.
This is a robustness check of four existing findings, not a new search. The seasons are the same 22 as before: effects at k = 2 and 3 are **not independent confirmations**, since the sustained tracks
are a subset of the original ones. Effective n is 22 seasons (about 5 index values per season), not tracks or weeks.

## What this does not show

- **Not a cause**, and not new data. The proxy's HF duration is not the archive's duration (proxy: 66% of HF tracks have at least 2 fixes; archive 74%).
- The PR 41 pattern was fitted to archive counts at k = 1 and cannot be refitted without the 17 GB of fields; a pattern fitted for sustained lows might do better.
- The Pacific depth variants and the 1979-2000 within-era checks have no duration and were not repeated.
- The 60N and terrain-type rules of PR 29 and 63 are applied as in the originals, per track, not re-defined per k.
- Duration is counted in total HF fixes; the consecutive-run version (S1) agrees. Fixes are 6-hourly in both the archive and the ERA5 tracks.
- The OPC area rule is not used (its edges were recalled, not sourced).

## Deviations and corrections (none changes a test; all written after outcomes were seen, so labelled)

1. **S3 middle class** was pre-registered "in the contrast tables" but was not in the first run; `contrast.py` was extended with the class and rerun. Every number of the first run is identical to the last digit (checked, the first-run files were then deleted).
2. **Clarification text for PR 63.** The `c365311` clarification calls P1 "the difference between the share effect without terrain-type fixes and the original". That is wrong: P1 is the share effect without terrain-type fixes itself (k = 1: 1.141); the difference is P2 (k = 1: 1.014). The test that was run, and the p value used, are the original's P1; P2 is reported as descriptive.
3. **Family C** was announced as 24 paired differences; PR 14's HF and share differences are identical (RR(all) does not depend on k), so there are 12 distinct ones, and BH is over 12.
4. **Decision rule 2 (size relative to k = 1)** needs a paired difference. The PR 63 quantities have none (the original script does not produce one); for them the k = 2 and 3 estimates are compared with k = 1 by their intervals, and the table says "paired difference not computed".
5. **Post hoc:** duration against count-matched gust cut (`posthoc_intensity.py`), described above. A lead, not a finding.
6. No k = 2 or 3 outcome was run before the plan was committed, and no code was changed after those outcomes beyond items 1 and 5 above.

## Verification

A fresh Sonnet agent, given the definitions and inputs but not the code, results or READMEs, recomputed with its own implementation (statsmodels Poisson) and **matched to three decimals**:
n_hf from the catalog (4,157 events; 1,412 / 866 / 589 with 1 / 2 / 3 fixes) and the survivor counts (Atlantic 1,007 / 664 / 439, brief 343; Pacific 825 / 542 / 374, brief 283); PR 14 RR(all), RR(HF_k),
RR(share_k) at k = 1, 2, 3 in both basins; PR 64 weekly RR(all), RR(HF_k), RR(share_k); PR 41 archive RR(HF_k) and the archive event totals (902 / 675 / 453, brief 226; 811 / 609 / 430, brief 199);
brief, middle, sustained RR and the sustained/brief ratios in all three designs; PR 63 NAO share given the Greenland high, F4 and P1 at k = 1, 2, 3.

**Not independently checked:** p and q values, all bootstrap intervals, the detectable-effect figures and SEs, f intervals, the PR 64 share-channel percentages, the paired differences (family C),
secondary S1 and S2 rows, the post hoc duration-versus-intensity rows, PR 63 P2 and the PR 64 cyclone-count rows (k-independent, reproduced from the originals).
