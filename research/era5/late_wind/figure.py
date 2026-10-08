"""Forest plot of the primary tests (odds ratio for LATE, season-clustered 95% interval). usage: python3 figure.py"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
t = pd.read_csv(os.path.join(R, "tests.csv"))
t = t[t.id.str.startswith("T")].iloc[::-1]
fig, ax = plt.subplots(figsize=(7.5, 3.8))
ax.errorbar(t.OR, range(len(t)), xerr=[t.OR - t.lo, t.hi - t.OR], fmt="o", color="#1f4e79", capsize=3)
ax.axvline(1, color="gray", lw=1)
ax.set_xscale("log")
ax.set_yticks(range(len(t)))
ax.set_yticklabels([f"{i}  {l}" for i, l in zip(t.id, t.label)], fontsize=8)
ax.set_xlabel("odds ratio for onset at or after the minimum-pressure fix (log scale)")
ax.set_title("RA-7, ERA5 proxy, pipeline A, 1,886 events, 2004-05 to 2025-26", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(R, "late_wind.png"), dpi=130)
