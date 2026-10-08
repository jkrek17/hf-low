# Status

The ledger for this repository: what is in flight, what exists, what does not, and what may not be relied on yet. Read it at the start of a session and update it at the end, in the same pull request as the work. When this file and git disagree, git is right and this file gets fixed.

Last updated: 2026-10-08, after this repository (`jkrek17/hf-low`) was split from `jkrek17/awips-tools` and three commits that landed there during the split were carried over. Commit hashes below are this repository's. Pull request numbers refer to `jkrek17/awips-tools` unless they say otherwise.

## Threads

| Branch | Head | Was | Ahead of main | What it holds |
|---|---|---|---|---|
| `main` | | | | The archive CSVs, the site, the build and publish tools, the site tests, and the session rules. No research code. |
| `claude/exciting-fermat-8vcvgq` | `5f05669` | integration branch | 30 | Pull requests 80 and 81 merged together, plus the finished ocean-masked run of ERA5 pipeline B and the finding that a gust-based criterion cannot be carried back before 2001 (`c9dc994`). Does **not** have pull request 82. |
| `claude/hf-lows-qc-mode` | `747d87d` | new during the split | 29 | A QC mode for the page that writes corrections to the spreadsheet. Branched from the integration branch as it stood at `73ef2d2`, so it lacks `c9dc994` and `5f05669`. |
| `claude/era5-hf-history` | `b964cd1` | PR 80, with PR 82 merged in | 17 | Precursor recovery, ERA5 pipeline A (`research/era5/hf_history/`), and per-cyclone P(HF). Does not have pipeline B or the commit guard. |
| `claude/era5-hf-probability` | `71013e5` | PR 82 | 16 | Same files as `claude/era5-hf-history`; that branch is this one plus a merge commit. |
| `claude/hf-lows-cleanup` | `b6e8f52` | PR 81 | 25 | Precursor recovery, ERA5 pipeline B as it stood before the masked run, collision review, the commit guard. Fully contained in the integration branch. |

**The research branches predate `CLAUDE.md` and this file.** A session started on one of them will not load the session rules. Merge `main` into a research branch before working on it; the merge adds files and changes none of the branch's own.

**No branch has everything.** The integration branch is missing the three P(HF) commits (`a3f6c8a`, `71013e5`, and the merge `b964cd1`). A trial merge of `claude/era5-hf-history` into it on 2026-10-08 (not committed) went through without conflicts. `claude/hf-lows-qc-mode` is missing the two newest integration commits.

## Decisions waiting on Jason

1. **Pipeline A before 2001.** Commit `c9dc994` finds that ERA5 gusts drift upward over 1979 to 2000 at fixed storm depth and concludes a gust-based criterion cannot be carried back before the archive. Pipeline A is a gust threshold applied from 1979, and its Atlantic trend (+1.7 per decade, p = 0.05) spans that period. The merge commit `5f05669` says the finding "is raised separately rather than acted on here". Nobody has tested pipeline A against it. Flag A's record before 2001, rebuild it on pressure depth, or run the same drift check on A's index?
2. **Two ERA5 pipelines.** A and B were built in parallel and differ in almost every definition (table below). Keep both with distinct names and purposes, or retire one?
3. **Merge order into `main`.** Nothing from the research branches is on `main` yet.
4. **Validation floor year.** Files say the proxy "cannot be validated before 2001" in some places and "before 2004" in others. Which is the statement of record?
5. **QC mode and the live sheet.** `claude/hf-lows-qc-mode` adds a path that writes to the spreadsheet. Should a work session ever use it against the real sheet, or only against the mocks in `tests/qc/`?
6. **Cutover** from `awips-tools` (section at the end).

## What exists

### Archive and site (`main`)

- Two basin CSVs under `data/hf_lows/`, exported from the spreadsheet.
- `python3 tools/build_hf_lows.py --check` on `main` reports 1,932 lows, 8,073 fixes, 25 seasons, 3 dropped rows.
- Site tests pass on `main` and on the integration branch (checked 2026-10-08 after the split): composite 39, playback 40, regress 50, teleconnect 37, plus the Python teleconnections test.
- `flat/` regenerates from `docs/` with `python3 tools/publish.py --no-fetch --deploy flat --flat --yes` and matches the committed copy apart from build stamps.
- On `claude/hf-lows-qc-mode` only: opened with `?qc`, the page lists 194 data-quality flags and lets a fix be edited. Edits post to a new `doPost` in `web/HFArchiveExport/Code.gs`, which needs its own `QC_TOKEN`, checks the cell still holds what the page showed, writes it, and logs before and after to a "QC log" tab. Tests with mocks are in `tests/qc/`.

### Precursor recovery (all four research branches)

- Fetch, parse, and track stages with tests. `data/hf_lows/precursors.csv` is committed: 6,658 rows (2,416 high confidence, 4,242 medium), covering 1,423 distinct events.
- Sidecars committed: `recovered_chain_hours.csv`, `archive_position_suspects.csv`. On `claude/hf-lows-cleanup` and the integration branch also `collision_pairs.csv` and `tools/review_collisions.py`.
- The build attaches precursors as a separate series: 696 of 1,928 hurricane-force events have a usable window.
- Validation numbers that may be quoted, and their limits, are in the `tools/track_hsf.py` docstring. Current statement: coverage 32% (16% at high confidence only); wrong-storm rate 3.5% on usable fixes; recall 50% on the fastest deepeners against about 72% for the rest.

### ERA5 pipeline A: `research/era5/hf_history/` (on `claude/era5-hf-history`, `claude/era5-hf-probability`, integration)

- A threshold on a gust index, calibrated on seasons 2021-22 to 2025-26, applied 1979 to 2025. Threshold 71.7 kt.
- Skill in `results/skill.txt`: AUC 0.990, POD 0.77, FAR 0.23, CSI 0.63, HSS 0.76, bias 0.99; leave-one-season-out threshold 71.4 to 72.3 kt.
- Catalog committed: 4,157 events and 4,154 matched null cases over 47 seasons.
- Per-cyclone P(HF) (`probability.py`, `results/probability.txt`) on `claude/era5-hf-history` and `claude/era5-hf-probability` only.
- Reproduction streams about 370 GB and needs go-ahead.

### ERA5 pipeline B: `event_fields.py`, `criterion.py`, `series.py` (integration branch; an earlier state on `claude/hf-lows-cleanup`)

- A per-moment probability from four features, fitted on early (2004-05), late (2021-25), and full (2004 to 2025) windows.
- The ocean-masked extraction is complete for all 47 seasons (`c9dc994`: 51,372 steps, none failed, about 190 GB). `hf_probability.csv.gz` on the integration branch holds 174,700 rows for 1979 to 2025.
- Skill at a count-matched cut, leave-one-season-out over 22 seasons (`c9dc994`): POD 0.727, FAR 0.274, HSS 0.705; gust alone gives HSS 0.659.
- Masking transitioning tropical cyclones (IBTrACS within 400 km) does not change skill.
- A basin term removes the pooled Atlantic and Pacific bias at no cost to skill. It is not in the committed series.
- `criterion-result.txt` holds the transfer test between windows; `series-result.txt` carries the stationarity flags.
- The copy of `hf_probability.csv.gz` on `claude/hf-lows-cleanup` is the older 13-season, pre-mask snapshot. Use the integration branch's.

### How A and B differ

| | A (`hf_history`) | B (`event_fields`) |
|---|---|---|
| Output | Event or not, by threshold; catalog with null cases | Probability per 6-hourly cyclone moment |
| Gust feature | Maximum within 800 km, points nearest this low only | Maxima within 300 and 500 km, plus area and gradient terms |
| Low detection | Below 1010 hPa | Below 1005 hPa and 4 hPa deeper than surroundings; scored at 990 hPa or below |
| Atlantic domain | 30-67N, 98W-10E | 30-70N, 80W-10E |
| Pacific domain | 27-67N, 135E-120W | 30-67N, 160E-120W |
| Months | June to May | September to May |
| Positive label | Track matched to an archive event, any fix category | Moment within 400 km of an archive HF fix of class "low" |
| Calibration | 2021-22 to 2025-26 | 2004-05, 2021-25, and 2004 to 2025 |
| Seasons committed | 47 | 47 on the integration branch |
| Ocean mask | From the start | Since `c9dc994` on the integration branch |

## What does not exist

- A single branch holding all of the research.
- Any research code or data on `main`.
- A gust-based record that may be used before 2001 (see the first gate below).
- Any test of pipeline A against the gust drift found in `c9dc994`.
- The basin term in the committed pipeline B series.
- Validation of either ERA5 pipeline before the archive begins. There is nothing to validate against.
- Repair of the archive position errors the tracker found. They are listed, not fixed.
- Tests for `tools/review_collisions.py` or anything under `research/`.

## Gates in force

Do not build on these without closing the gate or stating the dependence.

- **Backfilled statistics across seasons.** The build reports "COVERAGE TREND: Bf not safe to compare across seasons". Headline coverage must come from recovered rows alone.
- **Recovered subset.** Not a random sample of storms; a rate from it describes it, not the archive.
- **Low-confidence precursor tier.** Refused.
- **ERA5 record is a proxy.** Label it so everywhere. No validation before the archive.
- **Gust-based indices before 2001.** `c9dc994`: over 1979 to 2000 the mean sea gust at fixed storm depth rises 0.65 kt per decade (t = 2.40) and gust exceedances rise (t = 3.00), while the same quantities are flat over the archive period. Applied backwards, a gust criterion would flag too few early events and manufacture an increase. This reverses the earlier gate (`d34dd3d`) that put the start at 1979. It binds pipeline B's series and, until someone shows otherwise, pipeline A's catalog and trends before 2001. Counts based on pressure depth did not drift.
- **Pipeline B series against the archive.** `series-result.txt` flags a trend in the series minus the archive over the 22 overlap seasons (+30.4 per decade, t = 2.63 against a critical 2.09). The commit notes that Mann-Kendall does not confirm it and that part of the divergence is the archive's own labelling. Do not read the series sum as a trend estimate.
- **Pipeline B `p_full`.** In-sample for seasons from 2004.
- **Pipeline B at high latitude.** Over-predicts by about 55% at 60 to 71N and in the most land-affected quartile (`c9dc994`).
- **Pipeline A before 2006.** Bias 2.05 against the archive in 2001-05 because the archive was still starting, so those seasons do not test the threshold.

## Statements in the repository that are stale or conflict

Fix these as the files are next touched, in the same commit.

- `research/era5/README.md` (all research branches) says there is no calibrated gust index, no extended record, and that pressure gradients run 11% stronger before 1979. Pipelines A and B exist, and `stationarity2` on `claude/hf-lows-cleanup` reversed the gradient result. Its script list omits both pipelines.
- `research/era5/hf_history/README.md` presents the 1979 to 2025 record and its trend with no mention of the gust drift in `c9dc994`.
- Commit `d34dd3d` and text that follows it give 1979 as a safe start year. `c9dc994` withdraws that for gust-based indices.
- `research/era5/hf_history/track.py` docstring names two output files; the code writes one.
- Validation floor: 2001 in `criterion-result.txt`, `research/era5/README.md`, and `hf_history/README.md`; 2004 in `series.py` and `series-result.txt`.
- Event totals: 1,928 in the tracker docstring and build message, 1,932 lows in the build count and later commits. These may count different things (hurricane-force events against all lows); nobody has written down which.
- Commit `c48a8a8` says precursors are attributed to 1,860 events; the committed file has 1,423 distinct event ids.
- Track counts: 8,138 in `hf_history/results/skill.txt`, 8,144 in `results/probability.txt`, although the later commit says the calibration reproduces exactly.
- High Seas cache size: about 70 MB in `.gitignore`, about 450 MB in commit `1c70427`.
- On the integration branch, `hf_history/README.md` and `results/era5_hf_catalog.csv` are the versions from before P(HF).
- Left over from `awips-tools`: `tools/publish.py` and `.gitignore` mention `docs/cps/`, `docs/img/`, and other paths that are not in this repository; `.claude/skills/clasp/SKILL.md` describes a different Apps Script app; `docs/README.md` omits several files from its list; `flat/README.md` is a copy of `docs/README.md`; the publish manifest is still named `.awips-publish-manifest.json`; User-Agent strings still say `awips-tools`.

## Open work

- Pipeline A: test its gust index for the drift in `c9dc994`; until then treat its record before 2001 as gated.
- Pipeline B: put the basin term in the committed series; look at the over-prediction at high latitude.
- Bring P(HF) into the integration branch, and the two newest integration commits into `claude/hf-lows-qc-mode`.
- Precursors: resolve the collision worklist (39 pairs: 31 sequential, 8 concurrent); recall on fast deepeners is not at parity and tuning was stopped deliberately; tropical and post-tropical systems are not parsed.
- The unexplained +2.6 residual in the explosive share (commit `a9be84c`).
- Hard-coded `/home/user/awips-tools/...` paths in `research/era5/` pilot scripts (`compare.py`, `envdemo.py`, `matchmonth.py`, `mismatch.py`, `sample.py`, `tip.py`). They will not run from a clone of this repository as written.

## Fixes wanted in the sheet

Corrections that belong in the spreadsheet, because a CSV-only fix is overwritten by the next fetch.

- 38 suspected position errors in `data/hf_lows/archive_position_suspects.csv` (research branches), including 8 first hurricane-force fixes the tracker refuses to anchor on. The QC mode on `claude/hf-lows-qc-mode` is a route for making such fixes in the sheet once it is merged and deployed.

## Cutover from awips-tools

`awips-tools` still holds its copy of these files and still publishes `/hf-lows/` on the public site. Until the steps below are done, **`awips-tools` is the copy that publishes** and changes made only here do not reach the public site.

**Do all hurricane-force work here from now on.** Sessions were still committing to `awips-tools` while the split was under way. Three commits that landed there were carried over (`c9dc994`, `5f05669`, `747d87d`). Anything committed to `awips-tools` after 2026-10-08 05:01 UTC is not in this repository; check before relying on either copy.

1. Done 2026-10-08: the setup pull request (number 1 in this repository) is merged.
2. Here: turn on GitHub Pages by hand, Settings > Pages > Build and deployment > Source = "GitHub Actions", then re-run "Deploy archive site to Pages". Its first run failed at the step that tries to enable Pages automatically; the data rebuild, the page validation, and the flat-sync check all passed.
3. In `awips-tools`: remove `hf-lows` from `OWN_FOLDERS` and the archive build from `site_publish.yml`, so that workflow copies `/hf-lows/` forward as it does other projects' folders.
4. Here: set the secret `WEB_TOKEN` and the variable `PUBLIC_SITE_REPO`. `publish-site.yml` then publishes `/hf-lows/`. Do not do this before step 3, or two repositories write the same folder.
5. In `awips-tools`: remove the moved files and leave a pointer to this repository. Close pull requests 80, 81, and 82 there with a link to the branches here.
6. Production: the NOAA copy is made by hand with `tools/publish.py` from a checkout. Point that checkout at this repository.

One deliberate difference: `publish-site.yml` leaves `data/qc-report.txt` out of the public copy, as `tools/publish.py` does for production. The `awips-tools` workflow copied it.
