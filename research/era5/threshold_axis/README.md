# Do the share results move with the gust threshold? (RA-23; ERA5 proxy, pipeline A)

Plan, committed before any 68 or 75 kt outcome: [PREREGISTRATION.md](PREREGISTRATION.md) (`15f66d5`). Code: `labels.py`, `thr_contrast.py`, `report.py`; the runs reuse
`sustained_hf/run_pr14.py`, `sustained_hf/run_channels.py` and `nao_share_barrier/barrier_share.py` unchanged on derived track tables. Results: [results/](results/)
(`summary.txt` readable; `headline.csv`, `paired.csv`). Reproduce (no ERA5 pull, about 15 minutes on 4 cores):

    python3 research/era5/threshold_axis/labels.py research/era5/threshold_axis/work
    # for t in t68 t717 t75: run_pr14.py, run_channels.py and barrier_share.py on work/all_tracks_$t.csv.gz (commands in sustained_hf/README.md)
    python3 research/era5/threshold_axis/thr_contrast.py {daily|weekly} results/thr_<design>.csv 2000 2000
    python3 research/era5/threshold_axis/report.py research/era5/threshold_axis/results

HF_T = pipeline A peak gust index at or above T kt (68, 71.7, 75); cyclones below T stay in every denominator. 71.7 kt is the control and reproduces PR 14, 63 and 64 exactly.
Per +1 SD of the index, Oct-Apr 2004-05 to 2025-26, 22 seasons.

## Answer

**The signs and the Atlantic-above-Pacific ranking hold at 68 and 75 kt, and the effects grow with the cut. Nothing flips. The one place the effect weakens is the low cut in the daily-count design (PR 14).**
Of 21 pre-registered headline tests, 18 hold (sign kept, BH q < 0.05). Three do not: PR 14 Pacific RR(HF) and RR(share) at 68 kt (1.034, 1.006: can't tell) and PR 14 Atlantic RR(HF) at 68 kt (1.050,
q 0.068; the preregistered rule calls this "collapse, well powered", because the detectable RR 1.055 is below the 71.7 estimate 1.096; it is an attenuation, with the same sign).
PR 64 holds at every cut in both basins.

| result | basin | 68 kt | 71.7 kt (original) | 75 kt |
|---|---|---|---|---|
| PR 14 RR(HF lows) | Atlantic NAO | 1.050 (1.01-1.09) | 1.096 | 1.124 (1.03-1.20) |
| | Pacific PNA | 1.034 (0.98-1.09) | 1.072 | 1.121 (1.03-1.22) |
| PR 14 RR(share) | Atlantic | 1.079 (1.02-1.13) | 1.125 | 1.155 (1.05-1.24) |
| | Pacific | 1.006 (0.96-1.06) | 1.043 | 1.091 (1.00-1.18) |
| PR 64 RR(HF lows) | Atlantic | 1.217 (1.18-1.26) | 1.291 | 1.288 (1.22-1.38) |
| | Pacific | 1.119 (1.07-1.18) | 1.165 | 1.217 (1.14-1.33) |
| PR 64 RR(share) | Atlantic | 1.198 (1.16-1.25) | 1.271 | 1.268 (1.19-1.37) |
| | Pacific | 1.084 (1.05-1.13) | 1.129 | 1.179 (1.11-1.28) |
| PR 64 share channel, part of change | Atl / Pac | 90% / 75% | 92% / 81% | 92% / 85% |
| PR 63 P1 share without terrain-type fixes | Atlantic | not run at 68 (see plan) | 1.141 | 1.157 (1.06-1.24) |
| PR 63 P3 mediation (share given Greenland high / share) | Atlantic | 0.916 | 0.923 | 0.891 |
| PR 63 F4 Greenland high on share | Atlantic | 0.889 | 0.870 | 0.832 |

Intervals are 95% season-block bootstrap (all in `headline.csv`). Tracks flagged HF in the daily-design window: Atlantic 1,425 / 1,003 / 703, Pacific 1,283 / 825 / 570.

- **Prediction 1 (sign and rank): holds.** Every headline ratio keeps its sign at both ends; the Atlantic share RR exceeds the Pacific one at all three cuts in both setups (6 of 6).
- **Prediction 2 (grows with the cut): mostly.** The paired ratio RR(HF at 75 kt)/RR(HF at 68 kt) (identical to the share ratio) is 1.070 and 1.084 (PR 14, Atlantic and Pacific, intervals above 1) and 1.058 and 1.087 (PR 64; the Atlantic interval
  0.999-1.123 just contains 1). After BH over the 12 paired differences only two pass (PR 64 Atlantic 68 to 71.7 and PR 64 Pacific 68 to 75, q 0.009); the rest have q 0.054 to 0.89. The Atlantic PR 64 effect plateaus from 71.7 to 75 kt (1.291, 1.288).
  So it is a lean, not an established gradient.
- **Where the effect is weakest.** At 68 kt the Pacific PNA share effect in the daily design is gone (1.006). This is the Pacific cell that was already unresolved at 71.7 kt (1.043, p 0.24) and was clear for sustained lows (RA-22).
  Lower cuts add many lows that are not strong; the Pacific PNA effect appears only when the cut selects the stronger storms.
- **PR 63 (Greenland high and terrain).** Mediation and the high's lowering of the storm HF share hold at both cuts and are, if anything, stronger at 75 kt (0.832 per SD).

## Is RA-22's lean toward sustained lows just a higher gust cut?

Largely yes, as far as these data can say. The cut that matches the number of tracks with at least 2 HF fixes is about 75.4 kt in both basins, and with at least 3 fixes 78.1 to 79.0 kt (peak gust of the N-th highest track, all months 2004-2025).
The 75 kt column here reproduces the 2-fix column of RA-22: PR 64 RR(HF) 1.288 / 1.217 against 1.316 / 1.219 for 2 fixes; PR 14 1.124 / 1.121 against 1.109 / 1.128. RA-22's post hoc comparison already found a duration-matched set indistinguishable from a count-matched gust cut.
Taken together, the rise of the pattern effect with "sustained" is the same rise as with a stronger peak gust. This does not exclude a separate duration effect (the fix-level gust table does not allow duration at 68 or 75 kt, so the two axes were not crossed).

## Power

Detectable RR at 80% power (exp of 2.8 x clustered SE) is in `headline.csv`: PR 64 1.05 to 1.12, PR 14 1.06 to 1.13, PR 63 1.06 to 1.12. At 68 kt there are 40% more events and intervals are narrower, at 75 kt 30% fewer and wider.
The 68 kt PR 14 Pacific null (RR 1.034, share 1.006, detectable 1.08) rules out effects above about 8% but not smaller ones.

## Looks

One look at the 22 gust-era seasons for pipeline A (PR 14, 63, 64 setups; leave-one-season-out index, nothing refitted), the 14th scoring of 2015-16..2025-26 by the agenda's count, logged in `hemispheric/results/heldout_looks.log`.
A robustness check, not an independent confirmation; the seasons are the same ones as before. Effective n is 22 seasons.

## What this does not show

- Not a cause. ERA5 reads low in the extremes, so a "true" HF cut may be above 71.7 kt in ERA5 units; stability from 68 to 75 kt says the mechanism acts on storm strength broadly, not that 71.7 is right. The range spans only about a tenth of the index range.
- PR 63 P1 at 68 kt was not run: `gustloc_fixes.csv` holds only fixes at 71.7 kt or above, so 68 kt tracks would be unclassified. (The script still writes a P1 row at 68; it is not valid and is not in the family.)
- P(HF) as a soft label was not done (per-cyclone P(HF) exists only for catalog events).
- Family A is BH over 21 tests; the decision rule "collapse" compares the detectable RR with the 71.7 estimate and fires for PR 14 Atlantic RR(HF) at 68 kt even though the sign and size are close to the original (1.050 against 1.096).

## Deviations (post hoc)

None to the plan. The descriptive comparison with the RA-22 gust cut equivalents (75.4 and 78-79 kt) was added when writing up; it uses counts only.

## Verification

A fresh Sonnet agent, given the definitions and inputs but not the code or results, recomputed with its own implementation and **matched to three decimals**: PR 14 RR(all), RR(HF_T), RR(share_T) at 68, 71.7, 75 kt in both basins; PR 64 weekly RR(all), RR(HF_T), RR(share_T);
PR 63 mediation ratio, F4 (0.890 / 0.870 / 0.831 against 0.889 / 0.870 / 0.832, inside tolerance) and P1 at 71.7 and 75 kt; the gust cuts that match the 2-fix and 3-fix counts (75.4 / 79.0 kt Atlantic, 75.4 / 78.1 kt Pacific); and the track counts on retained days
(Atlantic 1,425 / 1,003 / 703, Pacific 1,283 / 825 / 570). The verifier's window-wide counts before dropping days without an index are Atlantic 1,430 / 1,007 / 706 and Pacific 1,286 / 825 / 570; the README quotes the retained-day figures, which are the ones in the fits.

**Not independently checked:** p and q values, all intervals, detectable effects, the PR 64 share-channel percentages, the paired differences (family C), f, the PR 14 Atlantic and Pacific f intervals, and the verdict labels.
