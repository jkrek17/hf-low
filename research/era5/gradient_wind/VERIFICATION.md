# Verification

Two fresh Sonnet agents recomputed the quoted numbers from raw inputs with their own code; a third (Opus) reviewed the claims. Scratch scripts were not committed (they ran in the session scratchpad).

## Matched (first verifier)
- Gradient and geostrophic wind from raw MSLP on 40 random rows (within 0.3%, a length-convention difference of 111,000 against 111,320 m per degree); gust match exact.
- All A-table rows for both latitude bands (n, median, shares at 38, 54 and 71.7 kt, gust / geostrophic); the stability medians and shares; the 4-6 and 12+ winter-bootstrap intervals; P1 (+0.0172), P2 (+0.1066), P2b (+0.0445); the 30-37.5N and 37.5-45N 8-10 core comparisons; the ASCAT table (within three cells per bin, medians within 0.4 kt).

## Matched (second verifier)
- Skill test: 407,221 box-passes and 848 days, M1 +15.04%, by winter, >=25 kt -2.29%, M1q +5.38% and +19.92%, intervals within 0.1 point.
- Isallobaric: I1 +0.0271, |Va| medians and shares, I3 -0.0379 (n 8,672 and 5,335), I2 +0.92% and -0.34%.
- P3 -0.0581 (interval -0.0757 to -0.0385 against my -0.0782 to -0.0387); no-root shares 3.896% core and 6.258% non-core (93% of the no-root rows have a negative discriminant, the rest are small positive discriminants in a different differencing scheme); P1c +0.2448; the 8-10 core and non-core shares.

## Differences
- P1b: first run +0.022 (-0.003 to +0.044); verifier +0.012 to +0.018 depending on the minimum points per winter and bin. Restated as about +0.01 to +0.02, interval includes 0.
- P3 curvature: the verifier's scheme did not replicate the original curvature exactly (137 finite-V rows have a negative discriminant in its scheme).

## Not independently checked
- Core membership was not rebuilt from raw pressure (the stored flag was used).
- The 8-10 core / non-core intervals were matched only to within 1 point (29.6-45.6 against 29.4-44.9).
- The sustained-wind ASCAT rows "within 800 km of a tracked low" (the selection comes from the cached ASCAT product).

## Claim review changes
Verdict words aligned with the plans (P2 yes by rule, P3 no fitting word, I2 no well powered, I3 no verdict); gust-not-sustained warning added; hurricane-force counts and days added; the single-bin core lead labelled; the circularity of ERA5 with ASCAT stated; the tendency's use of future pressure stated; post hoc items logged in the plans.
