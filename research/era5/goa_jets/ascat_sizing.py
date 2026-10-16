"""ASCAT sizing (PREREGISTRATION_followup.md): how many cells of the existing ASCAT cache fall in the Gulf of Alaska zones, and what a fresh RSS pull for the
episode days (episode day and the day before) would cost. Writes results/ascat_sizing.txt. Reads /mnt/project-files/ascat-check (not committed)."""
import os, glob
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); A = "/mnt/project-files/ascat-check"
n = tot = 0
for f in sorted(glob.glob(A + "/cells/cells_*.npz")):
    r = np.load(f)["rows"]; lat = -89.875 + 0.25 * r[:, 0]; lon = 0.125 + 0.25 * r[:, 1]
    n += int(((lat >= 53) & (lat <= 61) & (lon >= 202) & (lon <= 230)).sum()); tot += len(r)
inv = pd.read_csv(A + "/rss_inventory.csv"); inv["d"] = inv.date.astype(int)
E = pd.read_csv(os.path.join(HERE, "results", "episodes.csv"), parse_dates=["time"])
days = set(E.time.dt.strftime("%Y%m%d").astype(int)) | set((E.time - pd.Timedelta(days=1)).dt.strftime("%Y%m%d").astype(int))
sub = inv[inv.d.isin(days)]
s = (f"ASCAT cache cells in 53-61N, 158-130W: {n} of {tot:,}\nepisodes (3 zones): {len(E)}; distinct episode days and day before: {len(days)}; RSS files on those days: {len(sub)}; "
     f"size {sub['size'].sum() / 1e9:.2f} GB\n")
open(os.path.join(HERE, "results", "ascat_sizing.txt"), "w").write(s); print(s)
