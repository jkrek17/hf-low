# Pre-registration: does the isallobaric wind explain part of the departure from geostrophic?

Written and committed before any isallobaric field was computed. Follows PR 282 (Jason, 2026-10-28: "a tighter gradient means a more disturbed atmosphere from geostrophic", then "your increasing isallobaric flow", then "sure"). Reuses cached pressure and the cached ASCAT/GFS tables; no download. Western Atlantic, ERA5 proxy.

## What has been seen

The gust-to-geostrophic ratio falls from about 0.8 at 4-6 hPa per 100 nm to 0.3 to 0.45 at 12 and over (PR 282). The curvature (gradient-wind) correction made the match worse (P3, -0.058). No isallobaric quantity has been computed.

## Definitions

- Pressure tendency: (MSLP(t+12 h) - MSLP(t-12 h)) / 24 h at each 00/12 UTC, ERA5 0.25 degree, smoothed with a 1-degree box. Times without both neighbours (the first and last of each winter) are dropped.
- Geostrophic wind vector: Vg = (1/(rho f)) k x grad p, rho = 1.25 kg/m3, central difference over +-0.5 degree (as in PR 282).
- Isallobaric wind vector: Va = -(1/(rho f^2)) grad (dp/dt), the gradient of the smoothed tendency over +-0.5 degree. It points toward the fastest-falling pressure.
- Corrected speed: S = |Vg + Va| (kt). A quasi-geostrophic estimate; it ignores friction and curvature.

## Planned tests (three, Benjamini-Hochberg across them)

- **I1** (match, 21 winters minus the edge times): among ocean cells with gradient of 4 hPa per 100 nm or more, Spearman(gust, S) minus Spearman(gust, |Vg|). Winter bootstrap, 2000 draws. A meaningful gain is 0.01 or more.
- **I2** (skill, 4 winters, as in `PREREGISTRATION_SKILL.md`): leave-one-winter-out gain in squared error of ASCAT wind from adding S to the GFS-wind baseline (1 + v_gfs + S) against adding only |Vg| (1 + v_gfs + |Vg|). All box-passes. A meaningful gain is 2% of MSE. Day-block bootstrap and the all-four-winters rule.
- **I3** (deepening versus steady): at gradient 6 to 10 hPa per 100 nm, the difference in median gust / G between rapidly falling pressure (tendency of -1 hPa per 12 h or faster) and rising pressure (+1 hPa per 12 h or more), winter-level medians, winter bootstrap. Two-sided (the direction is not predicted: the isallobaric wind is mostly across the isobars, so it adds to the speed only where the tendency pattern is displaced from the pressure pattern).

Verdicts as in `PREREGISTRATION.md`: "yes" if the interval excludes 0 in the stated direction and clears the meaningful size; "no (well powered)" if the interval is inside the meaningful size; otherwise "cannot tell".

## Descriptive

The size of |Va| by gradient bin, how often it is more than 5 kt, and the share of Va along the geostrophic direction versus across it.

## Limits stated now

12 hourly differences are a coarse tendency (a centred 24 h change), the 1-degree smoothing removes the small-scale pattern, and the quasi-geostrophic form is an approximation in strong flow. ERA5 gust is a proxy and the box-pass test is below about 45 kt.

## Changes after seeing results

Logged below as post hoc.
