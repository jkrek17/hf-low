# Gust search radius around the model low (HF), GFS PROXY vs OPC HF calls

Pre-registered in `PREREG.md` (commit 842b97d, before scoring). Code `radius_sweep.py`, results `radius_sweep.json`. Data: `../results/raw_gust_fixes_96h.csv.gz`.

Result: the radius barely matters. Gust maximum within 270 nm (500 km), 432 nm (800 km) or 648 nm (1200 km) of the model low gives matched-threshold HSS within about 0.03 of each other at every lead, both basins. By the pre-registered rule the smallest radius (270 nm) is enough in 17 of 18 basin-lead cells; the exception is Pacific 60 h (432 nm, +0.021). Larger radii are slightly worse in several cells (Atlantic 24 h: HSS 0.527, 0.505, 0.475; 72 h: 0.434, 0.413, 0.396).

Not tested: radii below 270 nm and the 500 nm (926 km) Jason asked for (not stored; fields not kept). The rule is signed, so a larger radius that scores worse counts as "stopped". The truth match (OPC HF low within 600 km) is unchanged. OPC records low positions, not wind areas.

Verified by a fresh agent recomputing all 54 cells without reading the code: all POD / agree / HSS within 0.002.
