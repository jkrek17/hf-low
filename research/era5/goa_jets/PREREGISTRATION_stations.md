# Gulf of Alaska jets: do coastal weather stations back the HRRR shoreline gusts?

Written and committed BEFORE any station record or HRRR/ERA5 value at a station was read for an outcome (the station list was chosen from NDBC metadata and from observation COUNTS only, no wind comparison). Asked by Jason 2026-10-17 08:05Z ("do what you think is best") after the HRRR study left open whether HRRR's extra shoreline gusts are real (`README_hrrr.md`, `README_sar.md`).
HRRR is a 3 hour MODEL forecast, ERA5 a PROXY; the station records are observations from exposed capes and islands, at different anemometer heights (not adjusted). Pipeline: this study's own zone maxima (`PREREGISTRATION_hrrr.md`), neither pipeline A nor B.

## Data and size (stated before the pull)

- Stations (NDBC C-MAN / coastal, historical standard meteorological files, GST column, 1 Oct 2018 to 30 Apr 2026): Z1/Z3 Pilot Rock `pila2` (59.74N, 149.47W; inside both boxes), Middleton Island `mdxa2` (59.44N, 146.33W; Z1; only 2024-25 on file, flagged short); Z2 Cape Decision `cdea2` (56.00N, 134.14W), Lincoln Island `lcna2` (56.06N, 132.69W), Guard Island `gixa2` (55.45N, 131.88W); Z3 East Amatuli Island `amaa2` (58.92N, 151.95W), Flat Island `fila2` (59.33N, 152.00W), Augustine Island `auga2` (59.38N, 153.36W). Cape Spencer `cspa2` (58.20N, 136.64W) lies 0.2 degree outside every zone box: used in H1 only. Chosen for being exposed outer-coast sites inside the zone boxes with at least 1,000 gust reports at 35 kt or more (counts only; Icy Bay, Pelican and inner-channel sites left out as sheltered or too few strong reports). Station heights vary and are not adjusted.
- Station files are small (under 5 MB each, about 40 MB in all). HRRR: only the GUST record of the same f03 files as the first HRRR study (0.78 MB per time x 3,394 = **about 2.7 GB, far under the 50 GB gate; the first study's fields were not kept**), read at the station's nearest 3 km cell and the maximum of its 3x3 neighbourhood. ERA5 gust: the cached 00/12 UTC file in `/mnt/project-files/era5-cache/gust_nh_12h/` (no pull), nearest 0.25 degree cell.
- Matching: each 00 and 12 UTC time takes the station record within 30 minutes of the hour (the largest GST if several). No imputation; a station-time without a record is dropped.
- Units: knots; distances in nautical miles where quoted.

## Hypotheses

**H1, calibration at the station.** For times when the station gust is 35 kt or more: per station the median HRRR/station and ERA5/station gust ratio (HRRR value = the larger of nearest-cell and 3x3 maximum is NOT used; nearest cell is the registered one, 3x3 maximum reported beside it), Spearman correlation of each with the station gust, and the number of times. Reading: "HRRR is the better coastal reference" if the HRRR median ratio is closer to 1 than ERA5's (in log) at at least 5 of the 8 stations that have 30 or more such times; otherwise "no better than ERA5". Station anemometers are not at 10 m, so a ratio near 1 is not a proof of accuracy.

**H2, are HRRR-only events seen at the coast?** Times with a zone gust of 54 kt or more in HRRR but not in ERA5 ("HRRR-only"), in both ("both"), and in neither ("neither"), using the committed zone tables (`results/hrrr_features.csv`, `stage0_zone_gust.csv`). Outcome: at least one of the zone's stations reports a gust of 45 kt or more within 30 minutes (only times when at least one zone station reported). Season-block bootstrap (8 blocks, 2,000 draws, seed 16). Registered reading by the HRRR-only rate against the neither rate, each zone, if there are at least 30 HRRR-only times with a reporting station: **supported** at 3 times or more; **partly supported** at 1.5 to 3 times; **not supported** below 1.5 times; fewer than 30 times: **cannot tell**. A few stations cannot see a jet elsewhere in the zone, so "not supported" would mean the stations show nothing, not that no jet exists. Both are one look at the same 8 seasons; the neither group is the large reference.

Tests: H1 (stations as units, no multiplicity correction needed for a count rule), H2 (3 zones; the rule is a ratio, no p-values). Minimum detectable effect is reported for each H2 zone as the rate ratio the counts could separate at 80% power.

## Reporting rules

Pre-registered verdict words are quoted as registered ("supported", "partly supported", "not supported", "cannot tell"). A null carries its interval. A fresh agent recomputes the point estimates before the PR; bootstrap intervals are flagged unrechecked if not. Deviations are logged below.

## Deviations

(none yet)
