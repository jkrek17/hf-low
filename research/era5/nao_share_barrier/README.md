# Is the Atlantic NAO share effect storm wind or barrier flow? (ERA5 proxy, pipeline A)

Plan, committed before any outcome was related to a predictor: [PREREGISTRATION.md](PREREGISTRATION.md) (`d642ff7`).
Code: [barrier_share.py](barrier_share.py) (analysis), [gh_pull.py](gh_pull.py) (daily Greenland high),
[gh_check.py](gh_check.py). Results: [results/](results/) (`summary.txt`, `tests.csv`, `headline.json`).
Everything is **pipeline A** (`research/era5/hf_history`, 800 km ocean gust index, HF-equivalent at 71.7 kt) and a
**proxy**. Atlantic, 1 Oct + 210 days, seasons 2004-05 to 2025-26 (22), NAO lagged to days -10..-4 before genesis,
the same setup as hf-low PR 14 (`freq_split`), which this reproduces exactly (9,636 tracks, 1,003 HF, RR(share) 1.125).

Reproduce (about 4 minutes, no ERA5 pull; `gh_daily.csv` is committed, its pull was 11 GB):

    python3 research/era5/nao_share_barrier/barrier_share.py research/era5/hf_history/results/all_tracks.csv.gz \
        <cpc_indices dir> . research/era5/nao_share_barrier/gh_daily.csv research/era5/highlat/gustloc_fixes.csv \
        research/era5/nao_share_barrier/results 2000

## Answer in plain words

**No. The barrier explanation does not account for the NAO share effect.** Taking out the Greenland barrier-type
fixes leaves all of it: the share ratio per SD of NAO goes from 1.125 (95% 1.05-1.19) to **1.141** (1.06-1.21)
when the 109 of 1,003 HF cyclones whose only HF-strength fixes are terrain-type (PR 29's classification) no
longer count as HF. The paired difference is +1.4% (95% -1.1% to +3.8%), so at most about 9% of the effect
could have been carried by those cyclones, and the test can detect a loss of 3.5% (about 30% of the effect).
The three Greenland masks agree (100 km 1.129, 300 km 1.132, storm-centre points only 1.150), and so does
dropping those cyclones from both sides (1.141). That is a well-powered null for "barrier fixes carry the effect".

Two further findings, one of which weakens the picture a little:

1. **The Greenland high does stand in for NAO, but not through barrier flow.** With the lagged Greenland-high index
   in the model the NAO share ratio falls from 1.125 to **1.039** (0.93-1.15), the pre-registered mediation
   prediction (at least half gone) is supported by its rule (r = 0.33, upper bound 0.90, q = 0.044 over all 20
   tests). The sign is the reverse of the barrier story, though: a stronger high **lowers** the share
   (0.870 per SD, 0.82-0.92). It does the same for HF cyclones with barrier fixes removed (0.862), so it works
   on storm HF, not on barrier HF. NAO and the high are correlated -0.69 in this sample (-0.59 at PR 58's case
   times), so NAO|GH is poorly determined and r itself has an interval from -1.25 to 0.90.
2. **Cutting everything north of 60N takes some of it off, not resolved.** Dropping the cyclones whose gust peak is
   at or north of 60N from both sides (36% of tracks, 36% of HF) leaves a ratio of 1.091 (0.99-1.19), q = 0.74
   of the effect (interval -0.22 to 1.16), paired difference -3.1% (-7.3% to +2.2%, p = 0.22). This is not the
   barrier route: removing barrier-type fixes changes nothing. NAO moves storms north (PR 11), so the 60N cut
   removes the storms the NAO effect lives in; the cut asks about location, not wind type.

## Numbers

Per +1 SD of lagged NAO, ratio and 95% season-block bootstrap interval (2,000 draws, 22 seasons). `q` is the
fraction of the log share effect retained, `log RR_share(variant) / log RR_share(original)`.

| HF definition | HF cyclones | RR(share) | q (95% interval) |
|---|---|---|---|
| original (PR 14) | 1,003 | 1.125 (1.05-1.19) | 1 |
| **T: not terrain-type (primary)** | 894 | **1.141 (1.06-1.21)** | **1.12 (0.91-1.56)** |
| G100: 100 km Greenland mask | 947 | 1.129 (1.05-1.19) | 1.03 (0.89-1.27) |
| G300: 300 km Greenland mask | 791 | 1.132 (1.04-1.21) | 1.05 (0.75-1.35) |
| R400: points within 400 km of the centre | 862 | 1.150 (1.07-1.22) | 1.18 (0.98-1.66) |
| D: T, reclassified tracks dropped from both sides | 894 of 9,527 | 1.141 (1.06-1.21) | 1.12 (0.93-1.52) |
| N60: tracks with gust peak at or north of 60N dropped from both sides | 642 of 6,207 | 1.091 (0.99-1.19) | 0.74 (-0.22-1.16) |

Pre-registered verdict table (q 0.75-1.25 "survives"): T and D fall inside; the upper end of the T interval
(1.56) lies in the "grows" row, so the verdict is not firm by the rule fixed in advance, but the lower end (0.91)
excludes both "partly barrier" and "mostly barrier". G100, G300, R400 and D have point estimates of 1.03-1.18
with intervals reaching from 0.75 to 1.66, so each is "survives, not firm". N60 (0.74) sits at the "partly" edge
with an interval that spans every row.

Where the NAO count effect sits (log-linear counts per SD): storm-HF cyclones (T) 1.111 (1.04-1.17), barrier-only
HF cyclones 0.979 (0.79-1.22, only 109 cyclones, minimum detectable 1.37, so **unresolved**), difference 1.135 (0.90-1.41).

The Greenland high (GH, mean MSLP over the ice sheet 64-78N, 55W-30W, 12 UTC, lagged days -10..-4; 1 SD = 9.7 hPa)
on the share: 0.870 (0.82-0.92) for original HF, 0.862 (0.81-0.92) for T, paired difference 0.990 (0.97-1.01).
NAO given GH, T definition: 1.054 (0.95-1.17).

Multiplicity: 20 pre-registered tests; **9 pass Benjamini-Hochberg q < 0.05 across all 20** (the NAO share effect
for P1, G100, G300, R400 and D; storm-HF NAO count; GH share for original and T; the mediation test P3). None of
the paired differences against the original (P2, G100, G300, R400, N60, D) passes (smallest q 0.14, R400 +2.2%).
Full table with p, family q and leave-one-season-out ranges: [results/tests.csv](results/tests.csv).

## How this relates to the Greenland-jet result (PR 58)

PR 58 found that, among times with a low near Greenland, a strong Greenland high raises the odds of a barrier-jet
gust of 71.7 kt or more about 4 times per SD, with NAO adding nothing once the high is in. This study asks a
different question (what share of all Atlantic cyclones reaches HF, by season-level lagged indices) and gets a
consistent but not identical answer: the high does absorb NAO (NAO 1.125 to 1.039), but the barrier-type HF
cyclones are only 11% of the HF cyclones (109 of 1,003), their NAO coefficient is unresolved, and the high's effect
on the share has the *opposite* sign (a weak high goes with more HF cyclones) and is the same with barrier-type
fixes removed. So the two results do not conflict: the barrier jet follows the high within a window where a
low is already present, while the cyclone-wide HF share follows the high through storm wind. The mechanism by
which a weak Greenland high raises the share (a deeper Icelandic low, a more open storm track, or both) is not
tested here, and the high is itself partly storm-made, which the 4-10 day lag reduces and does not remove.

## What this does not show

- It cannot separate proxy error from archive omission near Greenland (PR 29). No station sees HF wind there.
- The barrier flags cover only HF-strength fixes of tracks that are events (3,303 fixes, 1,103 of 1,104 Atlantic
  HF tracks 2004-05 on). One HF cyclone in the window has no fix rows and stays HF in every variant.
- Barrier-type means PR 29's positional rule (gust maximum more than 400 km from the centre and within 300 km of
  Greenland), not wind direction; some storm wind near the coast is removed by G300 and R400.
- NAO and the high are correlated -0.69, so a mediation estimate is imprecise and "NAO given GH" does not mean the
  NAO has no effect: its interval, 0.93-1.15, includes 1.125.
- No terrain mask exists in pipeline A, and Atlantic fixes north of 60N are lower-confidence throughout. 22 seasons
  is the sample; no held-out look was used. Transitioning tropical cyclones are in.

## Deviations and clarifications

None changes a pre-registered rule.

1. **Debug run seen.** A 20-draw debug run of the full script printed the point estimates before the 2,000-draw
   run; the final run changed only the intervals and p-values. Nothing was changed in response.
2. **Clarification of F2 for N60 and D** (made before any outcome was read, but after the plan was committed):
   the plan said "paired difference against the original HF fitted on the same sample". For N60 and D the
   variant is a different sample, so the comparison is against the original on the full sample; this is what the
   code does and what the table shows.
3. **Bug fixed before the full run:** a column named `T` collided with the DataFrame transpose attribute, and the
   basin subset needed `reset_index`. Both were crashes, not estimates.
4. The correlation between lagged NAO and lagged GH in this sample is -0.69 against PR 58's -0.59 (different times,
   different lag); the plan's variance-inflation figure (1.24) is therefore 1.39 here.

## Verification

GH series: matches PR 58's stage-1 GH at all 2,277 shared 12 UTC times (correlation 0.99999996, largest difference
0.005 hPa), `gh_check.py`. PR 14 numbers reproduced exactly. Independent recomputation of the numbers above: see
the section appended below by the verifier hand-off.
