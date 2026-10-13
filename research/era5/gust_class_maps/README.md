# ERA5 gust cache and monthly gale / storm / hurricane-force gust maps

**ERA5 PROXY (reanalysis instantaneous 10 m gust), NOT DIRECT OBSERVATION.** ERA5 gusts read low in extreme storms; near coasts, mountains and sea ice they are not reliable.

## What was kept (reuse this instead of pulling gust again)

In the project's shared folder, `/mnt/project-files/era5-cache/` (outlives containers; not in git):

- `gust_nh_12h/gust_0000.npz` ... `gust_0160.npz`: ARCO-ERA5 instantaneous 10 m gust, northern hemisphere (rows 90N to 0N, 361 x 1440, 0.25 degree, longitude 0 to 359.75), **00 and 12 UTC, 1 June 2004 to 31 May 2026** (16,070 times, 100 per file), stored as uint8 in 0.5 kt steps (`value = round(gust_kt * 2)`, clipped at 127.5 kt), zstd. 4.1 GB. Each file holds `hours_since_1900`, `shape`, `data`. Decode: `numcodecs.Blosc().decode(z["data"].tobytes())`, `np.frombuffer(..., np.uint8).reshape(z["shape"])`.
- `gust_monthly_counts.npz` (11 MB): per calendar month, cutoffs 34, 38, 50, 54, 64, 71.7, 80, 90 kt, per cell, the count of sampled times at or above the cutoff, plus times per month and the highest cache value per month. Hours per season in a month = count x 12 h / 22 seasons.
- `lsm_nh_0p25.npy`: ERA5 land-sea mask (ocean where < 0.5), same grid.

Pull cost: 51.3 GB streamed once (about 10 minutes, 12 processes), Jason's go 2026-10-13. Streaming rate was 3.2 MB per time, so any other period or 6-hourly sampling scales from that.

## Code

- `extract_cache.py` pulls and writes the cache (resumable).
- `monthly_counts.py` reads the cache and writes the monthly counts.
- `tools/build_gust_classes.py` (site, on `main`) turns the counts into `docs/data/gust-classes.js`.

## Limits

12-hourly instantaneous samples, so every hours-per-season value is a lower bound on true gust hours (each sample stands for 12 h). Sampling quantisation is 0.5 kt. 2004-05 on only; 1979-2003 is not in the cache (it would be about 25 to 30 GB more and, by decision 1, would only be a within-era view).

## Classes

Gale at or above 38 kt, storm at or above 54 kt, hurricane force at or above 71.7 kt are **nominal proxy classes** on the gust index. Only hurricane force is tied to the archive (count-matched); gale and storm are not calibrated on observed events. The earlier `gust_climo` maps use shares of times at or above 50 and 64 kt for October to April only.
