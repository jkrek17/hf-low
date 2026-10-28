# Pre-registration: pressure gradient to wind speed, western Atlantic (options A and C)

Written and committed before any ERA5 gust was matched to any pressure gradient. Requested by Jason 2026-10-28 ("link a pressure gradient to a wind speed", then "interested in seeing A, but my gut says C"; scope card: western Atlantic now, in gust, no download).

## What has been seen

- Sizing only: file layouts of the cached MSLP crop, the gust cache, the air/sea temperature crop and the ASCAT cells. No gust, wind or gradient value for this question.
- Related earlier results that shape the expectation, not this test: Gulf Stream rounds (PR 263, 269, 271) found the ASCAT wind runs about 10% above the ERA5 pressure-gradient wind when the water is 6 C or more warmer than the air, and the ratio of wind to gradient wind fell as the gradient rose (buoys, 0.84 to 0.46, post hoc). The same gradient-wind code (`gulf_stream_region/build_boxes.py::vg_field`) is reused.

## Question

How much surface wind goes with a given pressure gradient over the western Atlantic in cold season, and how does it change with latitude, with air-sea stability, and inside a storm core? Answer as a table a forecaster can read: gradient in hPa per 100 nm against the median ERA5 gust and the chance of gale, storm and hurricane-force gust.

## Data (all cached, no download)

- MSLP: `/mnt/project-files/gulf-stream/msl/msl_<year>.npz` (uint16, +85000 Pa; 0.25 degree, 55N down to 20N, 270E to 310E; 00 and 12 UTC, October to April, 2004-05 to the last winter present).
- Gust: `/mnt/project-files/era5-cache/gust_nh_12h/` (ERA5 instantaneous 10 m gust, kt, 0.5 kt steps), same times.
- Air and sea temperature: `/mnt/project-files/gulf-stream/era5temp/` (2 m temperature at 00/12 UTC, daily SST; crop 45N to 25N, 275E to 310E).
- Land mask: `era5-cache/lsm_nh_0p25.npy`.
- Secondary check only: ASCAT cells `/mnt/project-files/ascat-check/cells/` (equivalent-neutral, saturates near 45 kt, only cells within 800 km of a pipeline A fix).
- **ERA5 is a proxy.** Gust is the cached wind quantity, reads 15 to 25% above ASCAT and is inflated in unstable air. Gale, storm and hurricane-force are nominal gust classes (38, 54, 71.7 kt); only 71.7 is count-matched to the archive. ERA5 assimilates ASCAT and buoys, so ERA5 wind against ERA5 pressure is partly self-consistent.

## Definitions

- **Sample**: ocean grid cells, every 8th cell in each direction (2 degrees), 30N to 45N and 78W to 55W, at least 1 degree from land, every 00/12 UTC in October to April, winters 2004-05 onward. The region is set by the temperature crop (45N) and the 1-degree land buffer.
- **Gradient**: magnitude of the MSLP gradient by central difference over plus and minus 0.5 degree. Reported in hPa per 100 nm (1 nm = 1.852 km). **Geostrophic wind** G = gradient / (1.25 kg/m3 x f), in kt, as in `vg_field`.
- **Gradient bins (hPa per 100 nm)**: <2, 2-4, 4-6, 6-8, 8-10, 10-12, 12 and over.
- **Latitude bands**: 30-37.5N and 37.5-45N.
- **Stability**: SST minus 2 m air temperature: stable (<0 C), weakly unstable (0 to <6 C), strongly unstable (6 C or more).
- **Outcome**: ERA5 gust at the cell (kt). Per bin: n, seasons contributing, median, 10th and 90th percentile, chance of gust at or above 38, 54 and 71.7 kt, and the median of gust divided by G.
- **Storm core (option C)**: a cell within 216 nm (400 km) of the lowest MSLP in its surrounding plus and minus 3 degree box, where that minimum is below 1000 hPa and at least 3 degrees inside the crop edge. All other cells are non-core. Core status uses the pressure field alone, not the pipeline A tracks (so the table is independent of the tracker).
- **Gradient wind (option C)**: solves V^2 kappa + f V = f G, where kappa is the signed curvature of the isobars computed from the MSLP smoothed with a 1-degree box (cyclonic positive), taking the root that is real and positive; cells with no real root are marked and counted.
- **Uncertainty**: 95% intervals by a bootstrap over winters (22 winters at most), because cells and times within a winter are not independent. The independent sample is winters, and rare classes (71.7 kt) will have few. Cells with fewer than 200 points or fewer than 8 winters are not tabulated.

## Planned tests (three primary, Benjamini-Hochberg across them)

- **P1 (stability)**: at gradient 6-10 hPa per 100 nm, all latitudes pooled, the median of gust / G is higher in strongly unstable than in stable air. Statistic: difference in the winter-level medians, bootstrapped over winters. Expected direction: unstable higher.
- **P2 (core)**: at gradient of 6 hPa per 100 nm or more, the median of gust / G is lower in core cells than in non-core cells. Expected direction: lower (curvature, tight cores smoothed by 0.25 degree).
- **P3 (curvature)**: in core cells with gradient of 4 hPa per 100 nm or more, the gradient wind explains more of the variation in gust than G does. Statistic: Spearman rank correlation of gust with gradient wind minus that with G, bootstrapped over winters.
- Verdict words: "yes" if the BH-adjusted interval excludes 0 in the expected direction; "no (well powered)" only if the interval lies inside plus or minus 5% for P1 and P2, or plus or minus 0.03 for P3; otherwise "cannot tell" with the interval stated.

## Descriptive outputs, not tests

The A table (all cells) and C table (core cells) by bin, latitude band and stability; the rule-of-thumb comparison (what the textbook 4 hPa per 60 nm at 45N gives against the median gust in that bin); counts of cells with no real gradient-wind root; and a note of every empty or thin cell.

## Secondary check (labelled as such)

ASCAT S1: for ASCAT cells in the same region within 1.5 h of a 00/12 UTC time, the ratio of ASCAT wind to G by gradient bin, restricted to bins whose geostrophic wind is under 45 kt (the satellite ceiling). ASCAT is equivalent-neutral and these cells lie within 800 km of a tracked low, so this is a check of the lower end, not of storm or hurricane force. Not part of the BH family.

## Known limits, stated now

- 12-hourly samples and 0.25 degree ERA5 smooth the strongest gradients, so the upper bins describe ERA5, not the real storm.
- The gust cache is gust, not sustained wind. A sustained-wind table would need a 10 m wind pull and is not part of this round.
- Western Atlantic only (30-45N). Nothing here is tested for the Pacific or the open ocean east of 55W.
- Hurricane-force exceedances in this box will be rare; those cells report counts and "cannot tell".

## Changes after seeing results

Any change is logged below as a post-hoc deviation with the date and what had been seen.
