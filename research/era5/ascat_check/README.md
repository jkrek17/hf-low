# ASCAT check of the ERA5 gust proxy (status: WAITING ON NETWORK ACCESS; no scatterometer result exists)

Question (Jason, 2026-10-13, item 2): does pipeline A's gust index agree with independent scatterometer winds?
**ERA5 proxy, pipeline A.** Plan: `PREREGISTRATION.md` (committed first, `5440a7b`, before any scatterometer value existed).

## State

- Every scatterometer host tried is blocked by the container's network policy (NOAA CoastWatch, NCEI, PO.DAAC/Earthdata,
  RSS, EUMETSAT, KNMI, Copernicus Marine); the AWS open-data buckets checked hold no ASCAT. Nothing has been pulled.
  Unblock: allow `coastwatch.pfeg.noaa.gov` and `data.remss.com`, then a source adapter writes the cell table below.
- Built and tested on the ERA5 side only: `era5side.py` (cache, sample, disc, interpolation), `recover_index.py` (R0),
  `collocate.py`, `analyse.py` (P1-P4, BH, season-block bootstrap), `test_synthetic.py`.
- **Not done:** the adapter that downloads and standardises a real scatterometer product (format depends on the source that
  becomes reachable), the real run, verification, PR to Q&A/agenda/ledger.

## Results so far (all ERA5-side or synthetic; none is evidence about ASCAT)

- **R0** (`results/r0_recovery_300.csv`, random 300 of 17,552 sample tracks, seed 0): of the 243 with stored index > 0,
  95.1 % recover pipeline A's index to within 1 kt from the cache (disc max over ocean cells within 800 km); the other 4.9 %
  are higher (neighbouring low's cells; median +3.3 kt). 57 of 300 have stored index 0 and are excluded (PREREGISTRATION
  deviation 1). Not independently recomputed.
- **Synthetic end-to-end run** (`results/synthetic_test.txt`): 1,638 footprints from 1,638 real storm times with a fake
  swath (wind = gust/1.3 + noise, saturating at 50 kt). Collocation ran, mean footprint cover 0.62, 1,546 pass the 40 % rule.
  It recovers rho 0.98, AUC 0.997, a +2.8 kt HF-vs-matched difference. This shows the code runs and that the
  saturating-wind construction does not break the tests; it says nothing about real winds. Its P3 standard error
  (0.08 kt) is a floor, because the fake swath is a deterministic function of ERA5. Real power must be recomputed from
  the real footprint count and the real scatter.

## Input contract for the real run

`cells.csv.gz` with columns `sensor,time,lat,lon,wind_kt,ok` (time UTC ISO; ok=1 passes the product's rain/ambiguity/
quality flags). Then `python collocate.py cells.csv.gz out` and `python analyse.py out_footprints.csv`.
Requires `numcodecs` and the shared cache at `/mnt/project-files/era5-cache`.
