# Archive counts against published OPC counts

**Question.** The archive has 46.0 Atlantic and 38.9 Pacific hurricane-force (HF) events per season (2004-05 to 2025-26). Published OPC counts are 15-23
per basin per season (Von Ahn, Sienkiewicz and Chang 2006, 2001-04), about 49 per season for both basins (Chelton et al. 2006) and about 25 per basin
(Jelenak et al. 2013, conference slides, 2000-2010). Do the two agree once definitions are matched?

**Answer.** Mostly yes. Two definitional differences, event duration and the area counted, take the archive from 46.0 and 38.9 to 23.7 and 21.8; with the
season window and tropical-cyclone removal it ends at 21.9 Atlantic and 20.8 Pacific, against 20.0 and 19.7 in Von Ahn and about 25 in Jelenak. What remains
is 10% (Atlantic) and 6% (Pacific) above Von Ahn and 12% and 17% below Jelenak. The same rules applied to ERA5 pipeline A (a proxy) in the three Von Ahn
seasons give 52 Atlantic and 47 Pacific events against 60 and 59 published (ratios 0.87 and 0.80, p = 0.51 and 0.29), which is not a detectable difference.

That agreement does not show which definition the published authors used. Several rules reach roughly 20 to 30 (see "What this cannot tell us"), the
rule chain was fixed before the results were seen but is not the only one that works, and the duration rule is an assumption about what a "warned storm" is.
So: the gap is a definition gap, not evidence that the archive is wrong, and the remaining size and the exact cause are not pinned down.

Pre-registration of the checks and the "explained" criteria: `PREREG.md` (committed first, `0efd654`). Deviations are listed at its end. Run:
`python3 -I research/era5/count_reconcile/reconcile.py` (about 6 minutes, committed files only, no ERA5 access). Output: `results/reconcile.txt`.
Archive = OPC HF archive on main; proxy = ERA5 pipeline A (`research/era5/hf_history`, 71.7 kt gust index), a proxy.

## Each check (mean events per season, 2004-05 to 2025-26, 95% season-block bootstrap)

| # | Check | Atlantic | Pacific | Share of gap (Atl / Pac) |
|---|---|---|---|---|
| B0 | at least one HF fix | 46.0 [42.1, 50.0] | 38.9 [36.0, 42.0] | gap to Von Ahn: 26.0 / 19.2 |
| 1 | unit: at least 2 HF fixes (drops one-fix, 6 h events) | 33.8 | 29.1 | 46.8% / 50.9% |
| 1 | at least 3 HF fixes (alone, not chained) | 22.4 | 20.5 | alone: 91% / 96% |
| 1 | HF or DHF fix; class `low` only | 46.0 / 44.0 | 39.0 / 38.9 | about 0 (DHF) |
| 2 | link split and duplicate IDs (11 sequential, 10 concurrent links) | 45.4 | 38.6 | 1.8% / 1.2% (chained) |
| 3 | OPC area (assumed: Atlantic west of 35W, 31-67N; Pacific east of 160E, 30-60N), first HF fix | 31.9 alone | 29.5 alone | 37.1% / 37.2% (chained) |
| 3 | same, lon edge 5 degrees outward / inward | 35.5 / 24.8 | 32.7 / 26.0 | |
| 4 | October-April only | 42.9 | 36.9 | 5.1% / 4.0% (chained) |
| 5 | drop tropical-cyclone-linked | 44.1 | 37.5 | 1.9% / 0.9% (chained) |
| 6 | analysis era | no step | no step | none |

Chain (1 at least 2 HF fixes, 2 linking, 3 area, 4 Oct-Apr, 5 TC): 46.0, 33.8, 33.4, 23.7, 22.4, 21.9 (Atlantic); 38.9, 29.1, 28.9, 21.8, 21.0, 20.8 (Pacific).
Share of the gap is measured on the step, in this order, so it depends on the order.

By the pre-set rule (a factor explains nothing if its chained share is under 10%): **explain nothing** are linking of split or duplicate IDs, the
October-April window, tropical-cyclone removal, the DHF category, and the analysis era. The archive contains no tropical-category fixes (TY, TS) at all. **Explain
the gap**: the duration rule (about half) and the area (about 37%).

- **Where the area matters.** 30.6% of Atlantic events first reach HF east of 35W and 23.9% of Pacific events west of 160E; the busiest Atlantic box (60-65N,
  40-30W) straddles 35W. The OPC edges are recalled, not read from a source here, and the 5 degree sensitivity moves the Atlantic area-only count from
  24.8 to 35.5. **This needs a source before it is quoted.**
- **Duration.** 26.3% of Atlantic and 25.1% of Pacific archive events have a single HF fix (6 h). Whether a published "warned storm" excludes them is not known.
- **Era (check 6).** Archive counts show no step at the dates recording practice changed: Pacific 39.8 before 2013-14 and 38.3 after (Welch p = 0.659),
  Atlantic 45.1 before 2017-18 and 47.2 after (p = 0.606; Benjamini-Hochberg q = 0.659 and 1.0). The share of single-fix events and of events outside the
  assumed area are also flat across seasons (exploratory). Mean events per season were 44.7 (Atlantic) and 40.3 (Pacific) in 2004-05 to 2009-10, the
  QuikSCAT era. So the gap is not a post-scatterometer or ASCAT-era change in the archive.
- **Start-up seasons.** The archive's Pacific tab begins February 2002 and the Atlantic tab September 2003, so 2002-03 (Pacific: 22) and 2003-04 (Atlantic: 11;
  Pacific: 27) are the earliest complete seasons. With no filter the archive already matches Von Ahn's Pacific total in those years (49 against 44, ratio 1.11,
  p = 0.68) and sits below the Atlantic (11 against 15, p = 0.56). That is a weak test: three basin-seasons, in the archive's start-up years.

## Reproducing the published counts

| Comparison | Rule | Archive or proxy | Published | Ratio |
|---|---|---|---|---|
| Von Ahn, Pacific 2002-03 and 2003-04, archive | none | 49 | 44 | 1.11 (p 0.68) |
| Von Ahn, Pacific 2002-03 and 2003-04, archive | area, at least 2 HF fixes | 32 | 44 | 0.73 (p 0.21) |
| Von Ahn, Atlantic 2003-04, archive | none / area + at least 2 fixes | 11 / 7 | 15 | 0.73 / 0.47 |
| Von Ahn, 2001-04 Atlantic, proxy (8b) | at least 2 fixes, area, Oct-Apr, no TC | 52 (19, 24, 9) | 60 (22, 23, 15) | 0.87 (p 0.51) |
| Von Ahn, 2001-04 Pacific, proxy (8b) | same | 47 (14, 21, 12) | 59 (15, 22, 22) | 0.80 (p 0.29) |
| Jelenak (slides), 2004-09 per season, archive | none | 44.7 / 40.3 | 27.5 / 25.6 | 1.62 / 1.57 |
| Jelenak, same | area | 28.0 / 30.7 | 27.5 / 25.6 | 1.02 / 1.20 |
| Jelenak, same | area, at least 2 fixes | 23.0 / 24.0 | 27.5 / 25.6 | 0.84 / 0.94 |
| Chelton, Oct 2003 to May 2005 (1.67 seasons, both basins), archive | none / area | 70.1 / 52.7 per season | 49 per season | 1.43 / 1.08 |

Proxy events per season without any filter are 46.7 Atlantic and 42.7 Pacific in 2001-02 to 2003-04, about twice the published counts, and 43.5 / 41.7 in
2004-05 to 2009-10. So the physical count of ERA5 HF-equivalent cyclones did not change between the Von Ahn years and the archive years; what the published
counts leave out is the same kind of thing the matched rules remove. Under the matched rules the proxy gives 14.3 Atlantic and 15.7 Pacific events per
season in 2004-05 to 2009-10, below the archive's 23.0 and 24.0 under the same rules: the proxy matches 64-68% of archive events, so it is a rough guide to
level, not a replacement for the archive.

## What this cannot tell us

1. **Which definition the published authors used.** One rule alone reaches the published range (at least 3 HF fixes: 22.4 and 20.5), others get most of the way
   (area: 31.9 and 29.5; at least 2 fixes and area: 23.7 and 21.8). The reduction to about 20 is therefore not a unique identification, and a reader
   should not take the two-factor rule as established. The choice of at least 2 fixes for the chain was made when the script was written, before any output was
   seen, but it was not in the plan (deviation 1).
2. **A direct test in overlapping years is thin.** Only Pacific 2002-03, Pacific 2003-04 and Atlantic 2003-04 are both published and archived; the archive was
   still starting. Applied to those years the matched rules undershoot (0.73 in the Pacific, 0.47 in the Atlantic) while no rule matches almost exactly. The proxy
   check (8b) is the only test that does not depend on the archive's start-up, and it is within noise but 13-20% low.
3. **Jelenak 2013 is conference slides** and its definition and exact season span are unknown (10.5 seasons assumed).
4. **The OPC area edges are recalled, not sourced**, and the duration rule is a guess about what "warned storm" means. A direct read of Von Ahn et al. 2006 for its counting
   rule would settle both; the literature review read it at abstract level, so this was not done.
5. **Not tested here:** whether the archive's HF category corresponds to an OPC warning issued at that time (the archive records the category of a fix, the published counts are
   storms for which a warning was issued); the identity of the 26% single-fix events (short-lived real storms, or HF entries that were issued and cancelled).

## Rows for Jason

No archive row was found wrong. `results/linked_ids_for_review.csv` lists 21 pairs of IDs that look like one low split in two (11) or recorded twice at the same time
(10); the build already flags 22 lows as `split`, and `data/hf_lows/collision_pairs.csv` has 39 pairs. Merging them changes the means by 0.6 (Atlantic) and 0.3 (Pacific) events per
season, so they do not matter for the reconciliation. Nothing was edited in the archive or the sheet.

## Verification

A fresh Sonnet agent, given the claims and the committed inputs but not the code, recomputed with its own code: events per season (46.0, 38.9); events with at least 2 and at
least 3 HF fixes; the share of single-fix events; the shares outside the assumed area and the area counts; October-April only; the early per-season counts; the era means and both Welch
p-values; and the proxy per-season counts under the matched rules (Atlantic 19, 24, 9; Pacific 14, 21, 12) and without filters (46.7, 42.7). All matched to the rounding shown except:
- The Atlantic single-fix share is 26.3% (README corrected from 26.4%), and the at-least-2-fixes Atlantic mean is 33.86 against 33.8 printed by the script (a fraction of one event across 22 seasons; the script
  drops duplicate timestamps within a low before counting fixes, which may account for it; not run down).
- The chain claim was given to the verifier without the ID-linking step. Without linking it gets 33.9, 24.0, 22.7 (Atlantic) and 29.1, 22.0, 21.2 (Pacific), each about 0.2-0.3 above the chain
  with linking (33.4, 23.7, 22.4; 28.9, 21.8, 21.0). That is the size of the linking step, so the difference is the brief, not a disagreement, but the **linked chain itself was not independently recomputed**.

**Not independently checked:** the linking step and the linked chain, all bootstrap intervals, the tropical-cyclone share (5.3% / 4.9%) and the TC step, the Chelton and Jelenak comparison rows,
the Von Ahn binomial p-values, the 5 degree area sensitivity, the single-fix and outside-area shares by season (exploratory), and the 8b rows for 2004-2025.
