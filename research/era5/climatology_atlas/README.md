# Climatology atlas of hurricane-force lows

Descriptive statistics only: no hypothesis is tested, so there is no pre-registration, and no trend is estimated (decision 1).
Two sources are always shown side by side and always named:

- **archive**: the OPC hurricane-force archive on `main` (`docs/data/hf-lows.json`), seasons 2004-05 to 2025-26 (22 seasons;
  season = 1 June to 31 May). HF-window metrics only, because recording practice changed in 2013-14 (Pacific) and 2017-18
  (Atlantic) and the whole-track metrics are not comparable across seasons.
- **proxy**: ERA5 **pipeline A** (`research/era5/hf_history`, 800 km ocean gust index, HF-equivalent at 71.7 kt). A proxy:
  "a cyclone whose ERA5 fields look like the ones OPC warned for as hurricane force".

Tropical-cyclone-linked events are **included** in every count and distribution (the archive cannot separate them). An HF event is
an archive low with at least one HF-category fix, or a pipeline A event. An event's month is the month of its first HF fix.
Brackets are 95% intervals from a season-block bootstrap (seasons resampled within basin, 2,000 draws, seed 1). All numbers
below are in `results/atlas.txt` and `results/historic.txt`, and were recomputed by an independent agent (see Verification).

Reproduce (reads committed files only, no ERA5 access; about 4 minutes):

    python3 -I research/era5/climatology_atlas/atlas.py
    python3 -I research/era5/climatology_atlas/historic.py
    python3 -I research/era5/climatology_atlas/figures.py

Files: `atlas.py` (tables), `figures.py` (ten figures), `historic.py` (1979-2003 candidate list); `results/*.csv` (every table
behind a figure), `results/events_archive.csv` and `events_proxy.csv` (one row per event, the inputs to every table), `results/fig*.png`.

## What the atlas says

**Where (fig1, fig2, fig3, fig10).**
- The busiest 5° x 10° box in both sources is the Atlantic 60-65N, 40-30W (southeast of Greenland): 96.8 HF-centre hours per
  season in the archive, 88.9 in the proxy. The Pacific's is 40-45N, 160-170E (off Japan and the Kuril Islands): 43.1 and 47.7.
- Median first-HF latitude: Atlantic 55.0N (archive), 53.5N (proxy); Pacific 45.0N and 42.2N. 22.7% (archive) and 28.4% (proxy) of Atlantic events
  first reach HF at or north of 60N; in the Pacific 0.4% and 0.3%. The proxy sits about 1 to 3.5 degrees south of the archive at first HF in the Pacific, in every month.
- Archive Atlantic first-HF latitude is lowest in October (52.1N) and February (52.5N) and above 55N in November to January;
  the Pacific is at 42.5-44.0N in January-March and 46-49N in October-December (fig10).

**How many, by month (fig4, `monthly_exposure.csv`).** Per season, events by month of first HF fix (archive; Oct, Nov, Dec,
Jan, Feb, Mar, Apr): Atlantic 3.9, 4.6, 8.5, 9.4, 8.8, 6.1, 1.8; Pacific 4.0, 5.1, 8.0, 8.6, 5.7, 4.0, 1.5. Proxy is within the
intervals of the archive in nearly every month. 93-95% of events first reach HF between October and April. Atlantic HF hours per
season peak in February (178.9 [146.2, 216.5]); Pacific in January (162.8 [136.1, 191.2]). The timing is analysed in the seasonal cycle thread (PR 45).

**Season to season (fig5, `season_counts.csv`).** Archive events per season: Atlantic mean 46.0, sd 9.6, range 30-67; Pacific mean
38.9, sd 7.5, range 27-51. Most active archive seasons: Atlantic 2006-07 (67), 2025-26 (66); Pacific 2007-08 (51), 2015-16 and
2023-24 (50). The archive and the proxy agree on the ranking of seasons moderately (Spearman 0.62 Atlantic, 0.72 Pacific, 22
seasons; intervals 0.29-0.83 and 0.48-0.85). The Atlantic and Pacific archive counts are uncorrelated (-0.05, interval -0.50 to 0.39). No trend is estimated.

**How strong (fig6).** Archive median minimum pressure 963 hPa (Atlantic) and 964 hPa (Pacific); 37.8% and 34.2% below 960 hPa,
15.0% and 11.7% below 950 hPa, 3.5% and 2.5% below 940 hPa. Proxy (ERA5 MSLP): 963.1 and 963.6; 41.2% and 37.7% below 960. Proxy
maximum 24 h deepening median 1.20 Bergerons (Atlantic) and 1.41 (Pacific); 64.8% and 77.2% reach 1 Bergeron. The archive's HF-window
deepening exists for only 314 of 1,011 Atlantic and 286 of 856 Pacific events (24% and 26% of those reach 1 Bergeron): that
is a selected subset, so use the proxy for deepening. The ten deepest analysed archive pressures are 920-929 hPa
(`results/atlas.txt`); they are unverified against the sheet.

**How long and how fast (fig7, fig8).** Hours at HF, median: 12 (Atlantic archive) and 18 (Pacific archive), 12 in the proxy
in both. Distance covered while HF, archive median 337 km (Atlantic) and 412 km (Pacific; at least a quarter of events have a single HF fix, so zero distance). Mean translation speed while HF, median 23.5 and 25.5 kt (archive), 24.6 and 29.4 kt (proxy).
While HF the mean motion is toward the north-east (north in a few Labrador Sea and Irminger Sea boxes) in every box with enough steps; 61% (Atlantic) and 77% (Pacific) of archive steps head
into the north-east quadrant. Genesis-to-HF timing is PR 16 and is not repeated.

**How many at once (fig9, `simultaneous.csv`).** At least one HF low is active in the Atlantic at 14.5% of 6-hourly times from
October to April (archive; proxy 14.5%) and in the Pacific at 12.5% (proxy 12.3%); two or more at 0.9% and 0.6% (Atlantic and Pacific archive); three at the same time
happens, four never; the maximum is 3 in both basins and both sources.

**Transitioning tropical cyclones (fig9, `tc_share_by_month.csv`, proxy only).** By the 400 km proximity flag (first look, not a tightened
definition), 5.9% of Atlantic and 5.5% of Pacific events. They are confined to the early season: first HF fix in September 65% (Atlantic) and 55% (Pacific), October 24% in both, November 3% and 1%, December to April essentially none.

## Historic storms, 1979-2003 (`results/historic.txt`, `historic_storms.csv`)

Within-era candidate list from the proxy: the 15 deepest and 15 highest-gust pipeline A events per basin before 2004-05. Nothing
here is confirmed. Ranking is within the pre-2004 sample only, by two separate measures (depth, and gust within-era), because pre-2001
gust drifts upward at fixed depth. Caveats: ERA5 pressure over the open ocean in the 1980s was constrained by fewer
observations, so depths may be too shallow; the archive is incomplete in 2001-02 to 2003-04, so absence there is not evidence;
tropical-cyclone-linked events dominate the Pacific gust list (flagged). As a plausibility check, the proxy's second-deepest
Atlantic storm is 10 January 1993 (915 hPa, 61.0N 14.2W), the date of the Braer storm; that date was recalled from memory, not
checked against a source. The median ERA5 minimum pressure of events differs little between eras (Atlantic 962.3 in 1979-2000, 963.4 from 2004; Pacific 964.0 and 964.1), a descriptive check only.

## What is not here, and what the map thread could show

- Not done: genesis maps before 2004 (about 370 GB of ERA5, needs Jason's go-ahead); a sea-state climatology (wave pull not sized); the
  tightened tropical-cyclone definition (agenda question 9); any trend (gated by decision 1).
- Map layers: `hf_hours_per_season_{archive,proxy}_{atl,pac}.csv` (box density), `events_archive.csv` / `events_proxy.csv` (first-HF,
  minimum-pressure and last-HF positions), `motion_steps.csv.gz` (box-mean arrows). The map thread already draws archive tracks; these are extra layers.

## Caveats

- The proxy is not the archive. Its Atlantic domain stops near 67N and it has no terrain mask (STATUS, high-latitude assumption), so
  the Greenland box in the proxy includes barrier-type gusts; about a third of fixes north of 60N are of that type (PR 29).
- Archive first fix is where OPC began warning, not genesis; archive "first HF fix" is censored by that in about two-thirds of Atlantic events.
- Pressure statistics use archive class `low` only (tip-jet and no-centre events, 42 of 1,867, have no analysed pressure).
- 22 seasons: intervals on monthly means are about +-15-25%; season-to-season differences between sources are within those intervals.

## Verification

A fresh Sonnet agent, given the claims and the committed inputs but not the code or README, recomputed with its own code: event counts;
monthly events and HF hours; season mean, sd, range and ranking; the three Spearman correlations; archive and proxy minimum pressure
medians and threshold shares; proxy deepening; first-HF latitude and the share north of 60N; the busiest boxes; hours at HF and distance
covered; the simultaneity shares; the tropical-cyclone shares for the whole sample, September and October; and the historic Atlantic
ranks 1-2, the Pacific minimum and the era medians. All matched to the rounding shown.

**Not independently checked:** bootstrap intervals (all of them); translation speeds and the motion-field figures (including the 61% and 77% north-east
shares); the month-by-month first-HF latitudes; the share of events inside October-April; the ten deepest archive pressures; the November-to-April
tropical-cyclone shares; the historic lists beyond the entries named above, and the Braer-storm date. The figures were inspected by eye, not recomputed.
