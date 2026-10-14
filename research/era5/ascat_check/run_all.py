"""Real run: RSS ASCAT cells (rss_adapter.extract output) -> collocate -> footprints -> analyse (P1-P4).
usage: run_all.py CELLS_DIR OUT_PREFIX [scat|scat+rad]"""
import sys, numpy as np, pandas as pd
import rss_adapter as R, collocate as C

cells_dir, out = sys.argv[1], sys.argv[2]
rule = sys.argv[3] if len(sys.argv) > 3 else "scat"
cells = R.load_cells(cells_dir, rain_rule=rule)
print("cells", len(cells), "ok share", round(float(cells.ok.mean()), 4), flush=True)
F, P = C.collocate(cells, with_pairs=(rule == "scat"))
F["accepted"] = F.cover >= C.MIN_COVER
F.to_csv(out + "_footprints.csv", index=False)
if len(P):
    P[["track", "basin", "season", "sensor", "lat", "lon"]].assign(
        wind_kt=P.wind_kt.astype("float32"), g=P.g.astype("float32")).to_csv(out + "_pairs.csv.gz", index=False)
print(len(F), "footprints,", int(F.accepted.sum()), "accepted")
