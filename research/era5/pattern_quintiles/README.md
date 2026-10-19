# Gale, storm and HF lows by pattern-index quintile

**ERA5 proxy, pipeline A (`hf_history`)**, 22 seasons 2004-05 to 2025-26, Oct-Apr weeks, tropical transitions in. Plan: `PREREGISTRATION.md` (`1cf78e5`), committed before any class count met the index. Index: the PR 41 leave-one-season-out pattern index (no refit). Classes by track peak gust index: HF >= 71.7 kt, storm 54 to 71.7, gale 40 to 54 (**the 40 kt cut is new and uncalibrated**; no archive gale series exists).

## Answer in plain words
The pattern score sorts **hurricane-force** weeks, not gale or storm weeks.
- Weeks in the top fifth of the score had about **2 times** as many HF lows as the bottom fifth in both basins, after allowing for the month (Atlantic 2.1, interval 1.75-2.60; Pacific 2.0, 1.54-2.60). Rule verdict: separates, q < 0.001.
- Atlantic HF lows per week by quintile, low to high: 1.10, 1.20, 1.40, 1.56, 2.41. Pacific: 0.86, 1.13, 1.05, 1.44, 1.82. The rise is steepest in the top fifth.
- Weeks with at least one HF low: Atlantic 70% in the bottom fifth, 90% in the top; Pacific 56% and 85%. So even the quietest fifth of weeks usually has one: the score shifts how many, and a week is rarely HF-free in the top fifth.
- Storm-force lows barely change (month-adjusted ratio 1.07 Atlantic, 1.05 Pacific; lower bounds just under 1: can't tell). Gale lows do not rise at all (0.89, 0.97).
- Total tracks per week rise only from about 11.6 to 12.6 (Atlantic), 11.8 to 13.0 (Pacific): the score adds HF lows, not lows.

## Limits
- Not a fresh hold-out: all 22 seasons were scored with the index again (`looks/pattern_quintiles.log`).
- Week-level, 22 seasons: intervals are season-block bootstrap (2,000 draws). Quintile cuts and month means were fixed in the resampling.
- Index and latitude/season move together (PR 41 follow-ups); month adjustment is the only control. No position control.
- HF is the proxy class; the 40 kt gale cut is arbitrary. Weekly peak-time assignment, not storm-relative.
- This says nothing about 24 h bombs; the 24 hPa/24 h standard was not scored.

Files: `quintiles.py`, `results/ratios.csv`, `results/by_quintile.csv`.
