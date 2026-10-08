# Why are our hurricane-force (HF) counts 1.5-2 times the published OPC ones?

Data-quality reconciliation, not a hypothesis test. This plan was committed before any reconciliation check was run.
Sources: the archive (OPC HF archive, `data/hf_lows/HF_Data_-_{Atl,Pac}.csv` and `docs/data/hf-lows.json`, both on main) and,
as a physical benchmark only, ERA5 **pipeline A** (`research/era5/hf_history`), a **proxy**. Pipeline A was count-matched to the
archive over 2021-26, so it cannot be an independent check of the archive's definition after 2004.

## What had been seen before this plan was written (so it is not claimed as blind)

- Atlas headline numbers: archive 46.0 Atlantic and 38.9 Pacific events per season, 2004-05 to 2025-26.
- The per-season archive and proxy counts in `hf_history/results/era5_hf_counts_by_season.csv`, including the archive's 2001-02 to 2003-04
  rows (Atlantic 0, 0, 11; Pacific 1, 22, 27) and the proxy's roughly 45 per basin in those seasons. So the archive's 2002-04
  rows are known to be incomplete (the archive was starting), and it is known that the proxy does not show 15-23 per basin in 2001-04.
- The busiest Atlantic box is 60-65N, 40-30W, which straddles 35W.

## Published benchmarks (from the literature review; Jelenak 2013 is conference slides, not peer reviewed)

| Source | Period | Count |
|---|---|---|
| Von Ahn, Sienkiewicz and Chang 2006 | 2001-02 to 2003-04 | Atlantic 22, 23, 15; Pacific 15, 22, 22 warned storms per season |
| Chelton et al. 2006 | Oct 2001 to May 2005 | 175 warnings, both basins (about 49 per season combined) |
| Jelenak et al. 2013 (slides) | 2000-2010 | 289 Atlantic, 269 Pacific (about 25 per basin per season) |

Reference values used for "gap": Atlantic 20.0 and Pacific 19.7 (Von Ahn mean), and 25 (Jelenak). The unit in each is a *warned* storm.
The OPC area of responsibility used below (Atlantic west of 35W and 31-67N; Pacific 160E eastward and 30-60N) is recalled, not read from a
source in this repository. It is a stated assumption, tested with a 5 degree sensitivity.

## Checks, in order, each with the count rule it applies

Baseline B0: an archive event is a low (ID) with at least one HF-category fix; season 1 June to 31 May by the ID's season; archive
seasons 2004-05 to 2025-26 for levels, and the raw rows for 2001-04. Count per season per basin, mean and season-block bootstrap interval.

1. **Unit.** (a) any HF or DHF fix (warned, including developing); (b) at least 2 HF fixes; (c) at least 12 h at HF (3 fixes); (d) class `low` only.
2. **Splits, merges and duplicates.** Link two IDs in one basin when one ends and the other starts within 12 h and 500 km (a low split in two), and
   merge IDs that overlap in time with fixes within 500 km of each other (duplicates). Also report the build's own `split` flag and
   `collision_pairs.csv`. Count events after linking.
3. **Domain and basin boundaries.** Fraction of events, and of first HF fixes, outside the assumed OPC area; count restricted to events whose
   first HF fix is inside it; the same with the boundary moved 5 degrees either way; positions that do not belong to their basin file.
4. **Season window.** Count events whose first HF fix is in October-April only; and November-March only.
5. **Tropical and post-tropical.** Drop archive events matched (`archive_events` in the pipeline A catalog) to a proxy event flagged as tropical-cyclone-linked
   (400 km IBTrACS proximity); drop archive rows with tropical categories (TY, TS, TYPHOON, TYPH). Archive-wide TC share is only known for the matched subset.
6. **Analysis era.** Mean archive events per season in 2004-05 to 2009-10 (QuikSCAT to 2009), 2009-10 to 2012-13, 2013-14 on in the Pacific,
   2017-18 on in the Atlantic (the dates recording practice changed). Two pre-set step tests (Pacific 2013-14, Atlantic 2017-18; Poisson
   regression on season counts) plus Benjamini-Hochberg across the checks that carry a p-value. No trend is fitted.
7. **Reproduction.** Compare archive counts under the rules above with Von Ahn, Chelton and Jelenak for the overlapping seasons. The archive is incomplete
   before 2004-05, so only seasons the archive can cover count as a test (see "Explained" below).
8. **Proxy as physical benchmark.** Proxy events per season in 2001-02 to 2003-04 against the published counts for the same seasons, and archive-matched
   versus unmatched proxy events after 2004. If the proxy shows about the archive's level in 2001-04, the published counts were *not* a lower physical count.

## What would count as "explained"

- For each factor: reduction in mean events per season, alone (marginal) and applied in the order above (sequential), as a share of the gap
  (baseline minus published reference). A factor "explains nothing" if its sequential share is under 10% of the gap (about 2 events per season).
- The gap is **reconciled** for a basin if, after the factors that are definitional (1-5), the archive's per-season mean in a window the archive covers
  lies within 20% of the published reference and a conditional-binomial (Poisson ratio) test against the pooled published count has p >= 0.05.
- If the archive after matching is still more than 20% above the reference, the **remaining gap** is reported with the candidate cause, and is stated
  as a question for Jason (a sheet fix) or as an unresolved definition (what the published authors counted), not as an error in the archive.
- All checks are reported, including those that explain nothing. Any change to the rules after seeing results is logged below as post hoc.

## Post hoc deviations

(none yet)
