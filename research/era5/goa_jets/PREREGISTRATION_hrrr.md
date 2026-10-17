# Gulf of Alaska jets with HRRR-Alaska (3 km model): does the 0.25 degree ERA5 proxy miss the jets?

Written and committed BEFORE any HRRR field was pulled or read for an outcome. Asked by Jason 2026-10-17 04:59Z ("ERA5 is kind of a bad way to study these smaller scale effects. what else can we do? use model data?"); he picked option 1 (HRRR-Alaska) at 05:21Z.
HRRR values are a MODEL forecast (3 hour forecast of the operational HRRR-Alaska), not observation; ERA5 is a PROXY. Pipeline: this study's own zone maxima, as in `PREREGISTRATION.md`; neither pipeline A nor B.

## Design decision (differs from the first study, on purpose)

The first study pulled ERA5 event times plus 600 random controls. That cannot be repeated with HRRR: events are defined by the HRRR gust, and an HRRR event that ERA5 did not call an event is exactly the case of interest,
so selecting times by ERA5 would hide it. **Every 00 and 12 UTC time, October to April, seasons 2018-19 to 2025-26 (HRRR-Alaska archive; 3,396 valid times)** is pulled, so rates are exact and need no control scaling.
Valid time T uses the HRRR-Alaska cycle 3 h earlier (21 UTC for 00 UTC, 09 UTC for 12 UTC) at forecast hour 3: spun up, not yet far from its start. (Cycles run every 3 h; the file at T itself is the f00 start state.)

## Data and size (stated before the pull)

`noaa-hrrr-bdp-pds`, `hrrr.YYYYMMDD/alaska/hrrr.tCCz.wrfsfcf03.ak.grib2`, byte ranges from the `.idx`: GUST surface 0.78 MB, UGRD and VGRD 10 m 1.49 MB each, MSLMA mean sea level 0.41 MB (probed from real files on 2023-01-15):
4.17 MB per time x 3,396 = **about 14.2 GB, under the 50 GB gate; no existing cache holds these fields.** Plus the LAND field once. Grid 1299 x 919 at 3 km (polar stereographic), cropped to the rows and columns covering 48-67N, 160-128W.
Missing files are counted; more than 2% in a season is reported. Decoder `grib2_hrrr.py` (the repository's packing 5.3 reader plus simple packing for the winds), geolocation and wind rotation `hrrr_grid.py`, checked against the file's own land mask and GFS winds.

## Definitions (same as the first study unless stated)

- Zones Z1, Z2, Z3 with the same boxes; **HRRR ocean cells** (LAND < 0.5 on the 3 km grid) inside the box. Zone gust = largest HRRR gust over them; zone wind = largest 10 m speed (winds rotated to east/north).
- Regimes from the zone-mean wind direction (FROM) at every time, same direction ranges (Z1, Z2 BARRIER 100-170, ONSHORE 170-260, OFFSHORE rest; Z3 GAP 300 through north to 060, INFLOW 090-200, OTHER rest).
- **Event** = zone gust at or above 54 kt, **strong** = 64 kt, hurricane force proxy 71.7 kt counted only. These are the first study's ERA5 numbers kept for comparability and are NOT calibrated to HRRR gust (a different gust diagnostic); read them as "HRRR 54 kt", and the buoy check below says what that is on a buoy.
- COASTPILE (Z1, Z2) and INTERIOR_DP (Z3) as before, from HRRR MSLMA and the HRRR land mask (distance to land from the 3 km mask). DEPTH and LOW_DIST are NOT used (the HRRR domain covers only 66% of the reference-low box); STAB is not used (no SST).
- Independent sample: the 8 seasons.

## Questions and tests (descriptive unless stated; no p-value is claimed)

**H1, how does ERA5 compare with HRRR?** At every time and zone: Spearman correlation of HRRR and ERA5 zone-maximum gust (ERA5 from `results/stage0_zone_gust.csv`); median ratio HRRR/ERA5; counts in the 2x2 table of events (HRRR 54 kt or more vs ERA5 54 kt or more), and the same for 64 kt, by zone and regime. Season-block bootstrap (8 blocks, 2,000 draws, seed 16).

**H2, where is the jet?** For HRRR events in the barrier or gap regime: distance from the peak-gust cell to the nearest land (km), median and the shares within 10, 25 and 50 km of land, beside the ERA5 figure (the first study: 29% of North Gulf barrier peaks in cells with at least 10% land).

**H3, buoys.** NDBC buoys 46001, 46061, 46080, 46076, 46082, 46083, 46077 (positions asserted against the NDBC pages, `buoy.py check`; print each buoy's HRRR nearest-cell land fraction). Times where the buoy gust at 00/12 UTC (record within 40 min) is at least 35 kt (selected on the observation, so the model cannot choose them): median and quartiles of HRRR gust / buoy gust and HRRR wind / buoy wind, and the same for ERA5 on the same times where available (`results/buoy_era5_gust.csv`, ERA5 wind from the earlier fields where pulled; else gust only). Also: at buoy gust >= 45 kt, how often HRRR zone gust is at least 54 kt. Buoy wind is at anemometer height (about 5 m) and is not adjusted.

**H4, the forecaster table redone with HRRR.** For each zone and regime: chance (%) that a 00/12 UTC time is an event or strong (events / times; exact), the ERA5 rate at the same times for comparison, and by band of COASTPILE or INTERIOR_DP (bottom half, 50th-80th, 80th-95th, top 5% of the regime's distribution) the same chance. Season-block intervals and the range across the 8 seasons. Plain rates in the write-up.

**Not done:** strong-versus-storm odds ratios (too few episodes in 8 seasons: stated, not attempted), DEPTH and LOW_DIST, stability, forecast skill (this is f03, not a forecast test), HF-class tests (counts only), anything before 2018.

## Caveats stated now

HRRR is a forecast model with its own gust diagnostic and its own biases; it assimilates little over the open Gulf; HRRR-Alaska versions changed inside 2018-2026 and are not reconstructed; 3 km still smooths the narrowest channels and fjords; 8 seasons is short, so "ERA5 reads low" is a statement about these 8 seasons and about HRRR, not a calibration. HRRR 54 kt is not ERA5 54 kt.
Reporting: name the pipeline and say HRRR is a model; nautical miles in text; plain rates; a fresh Sonnet agent recomputes every quoted number; hf-reviewer before the PR; looks logged in `research/era5/looks/goa_jets_hrrr.log`.

## Deviations

- 2026-10-17, before the full pull (commit 9d6fddb): crop south edge widened to 48N so the 300-400 km coast pile-up cells are inside the stored window; a first start of the pull with the narrower edge was stopped and its output deleted (no outcome read).
- 2026-10-17, plumbing run of `hrrr_analyse.py` on 1,113 times (which times was not recorded) to catch code errors; its tables were not read for results and were overwritten by the full run.
- 2026-10-17 (POST HOC): `hrrr_sensitivity.py` restricts the HRRR zone maximum to ocean cells at least 10 km and at least 25 km from land, because H2 showed the registered zone maximum sits within 10 km of land for most HRRR events in the North Gulf and Cook Inlet. Not registered; reported as a sensitivity, one more look at the same 8 seasons.
- 2026-10-17: 2 of 3,396 valid times have no file (2019-03-10 12 UTC and 2019-03-11 00 UTC, 404 after retries); both fall in one season (0.5% of it).
- 2026-10-17: the bootstrap used 1,000 draws (seed 16) in `hrrr_analyse.py`, not the 2,000 registered; distances are reported in kilometres, not nautical miles; the ERA5 H2 distances were computed fresh with the same land mask rather than taken from the first study.
