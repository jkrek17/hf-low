"""The atlas extended back to 1979, as ERA5 PROXY only (there is no archive before 2001, and it is incomplete to 2004).

usage: python3 -I research/era5/climatology_atlas/era1979.py        (writes results/era1979/ tables and figures)

WARNING that goes on every figure and table here: this is ERA5 pipeline A, a PROXY, NOT DIRECT OBSERVATION. Decision 1
(Jason, 2026-10-08): the gust index drifts upward before 2001 at fixed storm depth, so gust-based levels and trends are not
comparable across 2001 and no 1979-2025 hurricane-force trend is claimed. Nothing here estimates or tests a trend.

How the pre-2004 data are used:
  - gust-based events (pipeline A, 71.7 kt) are shown in three SEPARATE eras: 1979-80 to 2000-01 (E1, 22 seasons),
    2001-02 to 2003-04 (E2, 3 seasons; archive incomplete) and 2004-05 to 2025-26 (E3, 22 seasons). Levels are never pooled
    across the E1/E3 break. Within-era shape (month shares, latitude, depth, speed) is what is compared by eye.
  - the one continuous 1979-2025 series is DEPTH-based: tracks whose ERA5 minimum MSLP is at or below a cut. The cut is
    fixed by the rule "count-matched on 2004-05 to 2025-26": the integer hPa at which the mean annual count of all pipeline A
    tracks (all_tracks.csv.gz) equals the mean annual count of gust-based events in that period, per basin. It is chosen on
    E3 only and not revisited. Depth is not subject to the gust drift, but ERA5 pressure over the open ocean in the 1980s was
    constrained by fewer observations, and depth-only and gust-only storms are different populations (hf-low PR 52: in the
    Pacific the depth-only storms sit about 8.8 degrees farther north). Season-to-season spread is shown without a trend line.
Seasons are 1 June to 31 May, labelled by start year. Tropical-cyclone-linked events are included. Intervals: 95%
season-block bootstrap (2,000 draws, seed 1), within era.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import atlas as A_  # noqa: E402
import figures as F_  # noqa: E402

OUT = os.path.join(HERE, "results", "era1979")
ERAS = {"E1 1979-2000": (1979, 2000), "E2 2001-2003": (2001, 2003), "E3 2004-2025": (2004, 2025)}
ECOL = {"E1 1979-2000": "#6a3d9a", "E2 2001-2003": "#999999", "E3 2004-2025": "#d17a22"}
WARN = "ERA5 PROXY (pipeline A), NOT DIRECT OBSERVATION. Gust-based levels are not comparable across 2001 (decision 1); no trend is claimed."
MONTHS = A_.MONTHS
BN = A_.BN


def banner(fig, extra=""):
    fig.get_layout_engine().set(rect=(0, 0.04, 1, 0.96))
    fig.text(0.5, 0.005, WARN + (" " + extra if extra else ""), ha="center", va="bottom", fontsize=8, color="#a40000", weight="bold")


def era_of(season):
    for k, (a, b) in ERAS.items():
        if a <= season <= b:
            return k


def depth_cut(T, P3):
    """Integer hPa per basin at which the E3 mean annual count of all tracks equals the E3 mean annual proxy event count."""
    cuts = {}
    for b in A_.BASINS:
        target = (P3.basin == b).sum() / 22.0
        d = T[(T.basin == b) & T.season.between(2004, 2025)]
        best = min(range(940, 1000), key=lambda c: abs((d.minp <= c).sum() / 22.0 - target))
        cuts[b] = (best, target, (d.minp <= best).sum() / 22.0)
    return cuts


def monthly_by_era(P):
    rows = []
    rng = np.random.default_rng(1)
    for era, (a, b_) in ERAS.items():
        seas = list(range(a, b_ + 1))
        n = len(seas)
        for b in A_.BASINS:
            d = P[(P.basin == b) & P.season.between(a, b_)]
            cnt = np.zeros((n, len(MONTHS)))
            for j, m in enumerate(MONTHS):
                cnt[:, j] = d[d.mon_on == m].groupby("season").size().reindex(seas, fill_value=0).values
            idx = rng.integers(0, n, (A_.NBOOT, n))
            bc = cnt[idx].mean(1)
            share = cnt.sum(0) / cnt.sum()
            bs = cnt[idx].sum(1) / cnt[idx].sum(1).sum(1, keepdims=True)
            for j, m in enumerate(MONTHS):
                rows.append(dict(era=era, basin=b, month=m, seasons=n, events_per_season=cnt[:, j].mean(),
                                 lo=np.percentile(bc[:, j], 2.5), hi=np.percentile(bc[:, j], 97.5),
                                 share_oct_apr=share[j], share_lo=np.percentile(bs[:, j], 2.5), share_hi=np.percentile(bs[:, j], 97.5)))
    return pd.DataFrame(rows)


def era_summary(P):
    rows = []
    for era, (a, b_) in ERAS.items():
        for b in A_.BASINS:
            d = P[(P.basin == b) & P.season.between(a, b_)]
            n = b_ - a + 1
            rows.append(dict(era=era, basin=b, seasons=n, events=len(d), events_per_season=len(d) / n,
                             median_minp=d.minp.median(), share_minp_lt_960=(d.minp < 960).mean(),
                             median_gust_kt=d.gust.median(), median_lat_first_hf=d.lat_on.median(),
                             median_lat_min_pressure=d.lat_mp.median(), median_hf_hours=d.hf_h.median(),
                             median_speed_kt=d.hf_spd.median(), median_maxdeep_berg=d.maxdeep.median(),
                             tc_share=d.tc.astype(float).mean()))
    return pd.DataFrame(rows)


def depth_series(T, cuts, P):
    rows = []
    for b in A_.BASINS:
        for s in range(1979, 2026):
            d = T[(T.basin == b) & (T.season == s)]
            rows.append(dict(basin=b, season=s, era=era_of(s), depth_events=int((d.minp <= cuts[b][0]).sum()),
                             gust_events=int(((P.basin == b) & (P.season == s)).sum()), depth_cut_hPa=cuts[b][0]))
    return pd.DataFrame(rows)


def main():
    os.makedirs(OUT, exist_ok=True)
    P, PF, _, _ = A_.load_proxy(1979, 2025)
    P["era"] = P.season.map(era_of)
    PF["era"] = PF.season.map(era_of)
    T = pd.read_csv(os.path.join(A_.HFH, "all_tracks.csv.gz"))
    P3 = P[P.era == "E3 2004-2025"]
    cuts = depth_cut(T, P3)
    out = [WARN, "", "Eras (never pooled for gust-based levels): " + "; ".join(f"{k} ({v[1]-v[0]+1} seasons)" for k, v in ERAS.items()),
           f"Proxy events: " + ", ".join(f"{k}: {int((P.era==k).sum())}" for k in ERAS), ""]
    for b in A_.BASINS:
        out.append(f"Depth cut {BN[b]}: <= {cuts[b][0]} hPa minimum MSLP (E3 mean annual count {cuts[b][2]:.1f} vs gust events {cuts[b][1]:.1f}; chosen on E3 only)")
    mt = monthly_by_era(P)
    mt["caveat"] = "ERA5 proxy, not observation; gust levels not comparable across 2001"
    mt.to_csv(os.path.join(OUT, "monthly_by_era.csv"), index=False, float_format="%.3f")
    es = era_summary(P)
    es["caveat"] = "ERA5 proxy, not observation; gust levels not comparable across 2001"
    es.to_csv(os.path.join(OUT, "era_summary.csv"), index=False, float_format="%.3f")
    ds = depth_series(T, cuts, P)
    ds["caveat"] = "ERA5 proxy, not observation; depth_events is the only series meant to be read across eras, gust_events only within an era"
    ds.to_csv(os.path.join(OUT, "season_series_1979_2025.csv"), index=False)
    out.append("\n== Era summary, gust-based proxy events (compare SHAPE across eras, not levels) ==")
    for _, r in es.iterrows():
        out.append(f"{r.era} {BN[r.basin]}: {int(r.events)} events, {r.events_per_season:.1f}/season; median min MSLP {r.median_minp:.1f}; share <960 {r.share_minp_lt_960:.3f}; "
                   f"median gust index {r.median_gust_kt:.1f} kt; median lat of first HF {r.median_lat_first_hf:.1f}N; min-MSLP {r.median_lat_min_pressure:.1f}N; "
                   f"hours at HF {r.median_hf_hours:.0f}; speed {r.median_speed_kt:.1f} kt; deepening {r.median_maxdeep_berg:.2f} Bergerons; TC share {r.tc_share:.3f}")
    out.append("\n== Share of Oct-Apr events by month of first HF fix, per era (shape) ==")
    for b in A_.BASINS:
        for era in ERAS:
            s = mt[(mt.basin == b) & (mt.era == era)]
            out.append(f"{BN[b]} {era}: " + ", ".join(f"{int(m)}: {x:.3f}" for m, x in zip(s.month, s.share_oct_apr)))
    out.append("\n== Depth-based series: events per season with minimum MSLP at or below the cut ==")
    for b in A_.BASINS:
        for era in ERAS:
            v = ds[(ds.basin == b) & (ds.era == era)].depth_events
            out.append(f"{BN[b]} {era}: mean {v.mean():.1f}, sd {v.std():.1f}, range {v.min()}-{v.max()} (no trend fitted)")
        v = ds[(ds.basin == b)]
        top = v.sort_values("depth_events", ascending=False).head(5)
        out.append(f"{BN[b]} five highest seasons 1979-2025 by depth: " + ", ".join(f"{int(r.season)}-{str(int(r.season)+1)[2:]} ({r.depth_events})" for _, r in top.iterrows()))
    open(os.path.join(OUT, "era1979.txt"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))

    # ---------------- figures
    fig, axs = plt.subplots(2, 2, figsize=(12, 7.5), constrained_layout=True)
    for j, b in enumerate(A_.BASINS):
        for i, (c, lab) in enumerate((("share_oct_apr", "share of the basin's Oct-Apr events"), ("events_per_season", "events per season (within era only)"))):
            ax = axs[i, j]
            for k, era in enumerate(ERAS):
                s = mt[(mt.basin == b) & (mt.era == era)].set_index("month").loc[MONTHS]
                x = np.arange(len(MONTHS)) + (k - 1) * 0.27
                lo, hi = ("share_lo", "share_hi") if c == "share_oct_apr" else ("lo", "hi")
                ax.bar(x, s[c], 0.27, color=ECOL[era], label=era, yerr=[s[c] - s[lo], s[hi] - s[c]], capsize=1.5, error_kw=dict(lw=0.6))
            ax.set_xticks(range(len(MONTHS)))
            ax.set_xticklabels([F_.MNAME[m] for m in MONTHS])
            ax.set_ylabel(lab, fontsize=9)
            ax.set_title(BN[b], fontsize=10)
            ax.legend(fontsize=7)
    fig.suptitle("Gust-based proxy events by month and era. Compare the SHAPE (top); the bottom levels are not comparable across 2001", fontsize=11)
    banner(fig)
    fig.savefig(os.path.join(OUT, "fig11_monthly_by_era.png"), dpi=115)
    plt.close(fig)

    fig, axs = plt.subplots(2, 2, figsize=(11, 8.6), constrained_layout=True)
    for j, b in enumerate(A_.BASINS):
        vmax = 0
        for era in ("E1 1979-2000", "E3 2004-2025"):
            a, b_ = ERAS[era]
            f = PF[(PF.era == era)]
            f = f.assign(season=f.season)
            la, lo, H = A_.box_hours(f, b)
            H = H * 22 / (b_ - a + 1)
            vmax = max(vmax, H.max())
        for i, era in enumerate(("E1 1979-2000", "E3 2004-2025")):
            a, b_ = ERAS[era]
            f = PF[(PF.era == era)]
            la, lo, H = A_.box_hours(f, b)
            m = F_.pcm(axs[i, j], b, H * 22 / (b_ - a + 1), la, lo, vmax)
            axs[i, j].set_title(f"{BN[b]}: {era} (proxy; {H.sum()*22/(b_-a+1):.0f} h per season)", fontsize=10)
        fig.colorbar(m, ax=axs[:, j], shrink=0.7, label="HF-equivalent centre hours per season per box")
    fig.suptitle("Where HF-equivalent lows sit, 1979-2000 against 2004-2025: compare the PATTERN; levels differ for reasons we cannot separate from the gust drift", fontsize=10)
    banner(fig)
    fig.savefig(os.path.join(OUT, "fig12_centre_hours_E1_vs_E3.png"), dpi=115)
    plt.close(fig)

    fig, axs = plt.subplots(2, 2, figsize=(12.5, 7.6), constrained_layout=True)
    for j, b in enumerate(A_.BASINS):
        d = ds[ds.basin == b]
        ax = axs[0, j]
        ax.plot(d.season, d.depth_events, "-o", ms=3.5, color="#1b7837")
        ax.axvline(2000.5, color="k", lw=0.8, ls=":")
        ax.set_title(f"{BN[b]}: DEPTH-based, tracks with min MSLP <= {cuts[b][0]} hPa (cut matched on 2004-25 only)", fontsize=9)
        ax.set_ylabel("events per season")
        ax.set_ylim(bottom=0)
        ax = axs[1, j]
        for era, (a, b_) in ERAS.items():
            e = d[d.era == era]
            ax.plot(e.season, e.gust_events, "-o", ms=3.5, color=ECOL[era], label=era)
        ax.axvline(2000.5, color="k", lw=0.8, ls=":")
        ax.set_title(f"{BN[b]}: GUST-based (pipeline A 71.7 kt): read within each colour only", fontsize=9)
        ax.set_xlabel("season (start year)")
        ax.set_ylabel("events per season")
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=7)
    fig.suptitle("Events per season since 1979. No trend line, none claimed. Dotted line = 2001, where gust levels stop being comparable", fontsize=11)
    banner(fig, "ERA5 pressure before ~1990 is less constrained by observations.")
    fig.savefig(os.path.join(OUT, "fig13_season_series_1979_2025.png"), dpi=115)
    plt.close(fig)

    fig, axs = plt.subplots(2, 3, figsize=(15, 7.6), constrained_layout=True)
    for j, b in enumerate(A_.BASINS):
        for k, era in enumerate(ERAS):
            e = P[(P.basin == b) & (P.era == era)]
            xs, med, q1, q3 = [], [], [], []
            for i, m in enumerate(MONTHS):
                v = e[e.mon_on == m].lat_on
                xs.append(i + (k - 1) * 0.2)
                med.append(v.median()); q1.append(v.quantile(.25)); q3.append(v.quantile(.75))
            axs[j, 0].errorbar(xs, med, yerr=[np.array(med) - q1, np.array(q3) - med], fmt="o-", ms=3, color=ECOL[era], capsize=2, label=era)
            axs[j, 1].hist(e.minp, bins=np.arange(920, 1011, 4), histtype="step", density=True, color=ECOL[era], lw=1.5, label=era)
            axs[j, 2].hist(e.hf_spd.dropna(), bins=np.arange(0, 80, 3), histtype="step", density=True, color=ECOL[era], lw=1.5, label=era)
        axs[j, 0].set_xticks(range(len(MONTHS)))
        axs[j, 0].set_xticklabels([F_.MNAME[m] for m in MONTHS])
        axs[j, 0].set_ylabel(f"{BN[b]}: latitude of first HF fix (median, quartiles)", fontsize=9)
        axs[j, 1].set_xlabel("ERA5 minimum MSLP (hPa)")
        axs[j, 2].set_xlabel("speed while HF-equivalent (kt)")
        axs[j, 0].legend(fontsize=7)
    fig.suptitle("Within-era shape of proxy events: latitude by month, depth and speed, three eras kept separate", fontsize=11)
    banner(fig)
    fig.savefig(os.path.join(OUT, "fig14_shape_by_era.png"), dpi=115)
    plt.close(fig)


if __name__ == "__main__":
    main()
