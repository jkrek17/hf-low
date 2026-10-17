# How many Sentinel-1 radar scenes cover the Gulf of Alaska storm times? (catalogue count only)

Asked by Jason 2026-10-17 ("sure", on the SAR-catalogue question). Counts of what exists in the Copernicus Data Space catalogue; no imagery was downloaded and no wind was retrieved. Script `sar_catalogue.py`, output `results/sar_catalogue.csv` and `sar_catalogue.txt`. Not pre-registered (a count, no hypothesis). The episodes are the first study's ERA5-PROXY 54 kt episodes (`results/episodes.csv`), October 2014 on (Sentinel-1 begins in 2014), control rows left out.

A scene is counted if its footprint touches the zone box (the zone boxes of the first study) and it is a Sentinel-1 wide-swath ground-range product (IW_GRDH_1S or EW_GRDM_1S).

| Zone | Episodes | Scene within 6 h of the episode time | Scene on the same UTC day | 64 kt episodes with a scene within 6 h |
|---|---|---|---|---|
| North Gulf | 181 | 36 (20%) | 61 (34%) | 7 of 28 |
| Panhandle | 172 | 34 (20%) | 61 (35%) | 4 of 26 |
| Cook Inlet | 172 | 41 (24%) | 59 (34%) | 4 of 19 |

Reading: about 1 episode in 5 has a radar scene in the zone within 6 hours, and about 1 in 3 on the same day; for the strongest (64 kt proxy) episodes 15 of 73 do. That is enough for a handful of case studies, not for a rate. A scene touching the zone box does not mean it covers the jet core near the coast, and Sentinel-1 gives a wind estimate only after a retrieval step this study has not set up (and the estimate is least reliable in the strongest winds). Counts are of scenes, not of independent storms. A fresh agent re-queried 12 of the within-6-hour counts live and all matched, and matched the per-zone totals; it did not check the same-day counts or the 64 kt row, and every Sentinel-1 product found was IW (the EW count was 0 in all checks, so the EW handling is untested).
