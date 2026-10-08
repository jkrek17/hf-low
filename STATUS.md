# Status

The ledger for this repository: what is in flight, what exists, what does not, and what may not be relied on yet. Read it at the start of a session and update it at the end, in the same pull request as the work. When this file and git disagree, git is right and this file gets fixed.

Last updated: 2026-10-08, when this repository (`jkrek17/hf-low`) was split from `jkrek17/awips-tools`. Commit hashes below are this repository's. Pull request numbers refer to `jkrek17/awips-tools`.

## Threads

| Branch | Head | Was | Ahead of main | What it holds |
|---|---|---|---|---|
| `main` | `1425324` | | | The archive CSVs, the site, the build and publish tools, the site tests. No research code. |
| `claude/exciting-fermat-8vcvgq` | `73ef2d2` | integration branch | 28 | Pull requests 80 and 81 merged together: precursor recovery, both ERA5 pipelines, the commit guard. Does **not** have pull request 82. |
| `claude/era5-hf-history` | `b964cd1` | PR 80, with PR 82 merged in | 17 | Precursor recovery, ERA5 pipeline A (`research/era5/hf_history/`), and per-cyclone P(HF). Does not have pipeline B or the commit guard. |
| `claude/era5-hf-probability` | `71013e5` | PR 82 | 16 | Same files as `claude/era5-hf-history`; that branch is this one plus a merge commit. |
| `claude/hf-lows-cleanup` | `b6e8f52` | PR 81 | 25 | Precursor recovery, ERA5 pipeline B, collision review, the commit guard. Fully contained in the integration branch. |
| `harness-setup` | | this setup | | `CLAUDE.md`, the `hf-harness` skill, this file, the README, the workflows. |

Every research branch is two merge commits behind `main`; those merges changed no files.

**No branch has everything.** The integration branch is missing three commits (P(HF): `a3f6c8a`, `71013e5`, and the merge `b964cd1`). The simplest route to one complete branch is to merge `claude/era5-hf-history` into `claude/exciting-fermat-8vcvgq`. A trial merge on 2026-10-08 (not committed) went through without conflicts: it adds `probability.py` and its two result files and brings `hf_history/README.md` and `results/era5_hf_catalog.csv` up to their P(HF) versions.

## Decisions waiting on Jason

1. **Two ERA5 pipelines.** A and B were built in parallel and differ in almost every definition (table below). Keep both with distinct names and purposes, or retire one?
2. **Merge order into `main`.** Nothing from the research branches is on `main` yet.
3. **Validation floor year.** Files say the proxy "cannot be validated before 2001" in some places and "before 2004" in others. Which is the statement of record?
4. **Cutover** from `awips-tools` (section at the end).

## What exists

### Archive and site (`main`)

- Two basin CSVs under `data/hf_lows/`, exported from the spreadsheet.
- `python3 tools/build_hf_lows.py --check` on `main` reports 1,932 lows, 8,073 fixes, 25 seasons, 3 dropped rows.
- Site tests pass on `main` and on the integration branch (checked 2026-10-08 after the split): composite 39, playback 40, regress 50, teleconnect 37, plus the Python teleconnections test.
- `flat/` regenerates from `docs/` with `python3 tools/publish.py --no-fetch --deploy flat --flat --yes` and matches the committed copy apart from build stamps.

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

### ERA5 pipeline B: `event_fields.py`, `criterion.py`, `series.py` (on `claude/hf-lows-cleanup`, integration)

- A per-moment probability from four features, fitted on early (2004-05), late (2021-25), and full windows.
- `criterion-result.txt` holds the transfer test between windows. `stationarity2.py` and its result are the basis for starting in 1979.
- `hf_probability.csv.gz` is committed but is a **partial, superseded snapshot**: 13 of 47 seasons, and extracted before gust features were masked to ocean points.

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
| Calibration | 2021-22 to 2025-26 | 2004-05, 2021-25, and all cached seasons from 2004 |
| Seasons committed | 47 | 13 |
| Ocean mask | From the start | Added later; committed outputs predate it |

## What does not exist

- A single branch holding all of the research.
- Any research code or data on `main`.
- A pipeline B probability series for all 47 seasons, or any pipeline B output made with the ocean mask. Commit `40c3ea5` started the masked re-extraction; no commit records it finishing. `event_fields_steps.csv` holds only seasons 2004 and 2005.
- A fit of pipeline B that includes 2006 to 2020.
- Validation of either ERA5 pipeline before the archive begins. There is nothing to validate against.
- Repair of the archive position errors the tracker found. They are listed, not fixed.
- Tests for `tools/review_collisions.py` or anything under `research/`.

## Gates in force

Do not build on these without closing the gate or stating the dependence.

- **Backfilled statistics across seasons.** The build reports "COVERAGE TREND: Bf not safe to compare across seasons". Headline coverage must come from recovered rows alone.
- **Recovered subset.** Not a random sample of storms; a rate from it describes it, not the archive.
- **Low-confidence precursor tier.** Refused.
- **ERA5 record is a proxy.** Label it so everywhere. No validation before the archive.
- **Years before 1979.** Sensitivity check only; 1958 to 1978 has not been run.
- **Pipeline B series trend.** Two borderline trends are flagged in `series-result.txt`; "the all-rows sum should not be read as a trend estimate".
- **Pipeline B `p_full`.** In-sample for seasons from 2004 and will change when 2006 to 2020 enter the fit.
- **Pipeline A before 2006.** Bias 2.05 against the archive in 2001-05 because the archive was still starting, so those seasons do not test the threshold.

## Statements in the repository that are stale or conflict

Fix these as the files are next touched, in the same commit.

- `research/era5/README.md` (all research branches) says there is no calibrated gust index, no extended record, and that pressure gradients run 11% stronger before 1979. Pipelines A and B exist, and `stationarity2` on `claude/hf-lows-cleanup` reversed the gradient result. Its script list omits both pipelines.
- `research/era5/hf_history/track.py` docstring names two output files; the code writes one.
- Validation floor: 2001 in `criterion-result.txt`, `research/era5/README.md`, and `hf_history/README.md`; 2004 in `series.py` and `series-result.txt`.
- Event totals: 1,928 in the tracker docstring and build message, 1,932 lows in the build count and later commits. These may count different things (hurricane-force events against all lows); nobody has written down which.
- Commit `c48a8a8` says precursors are attributed to 1,860 events; the committed file has 1,423 distinct event ids.
- Track counts: 8,138 in `hf_history/results/skill.txt`, 8,144 in `results/probability.txt`, although the later commit says the calibration reproduces exactly.
- High Seas cache size: about 70 MB in `.gitignore`, about 450 MB in commit `1c70427`.
- On the integration branch, `hf_history/README.md` and `results/era5_hf_catalog.csv` are the versions from before P(HF).
- Left over from `awips-tools`: `tools/publish.py` and `.gitignore` mention `docs/cps/`, `docs/img/`, and other paths that are not in this repository; `.claude/skills/clasp/SKILL.md` describes a different Apps Script app; `docs/README.md` omits several files from its list; `flat/README.md` is a copy of `docs/README.md`; the publish manifest is still named `.awips-publish-manifest.json`; User-Agent strings still say `awips-tools`.

## Open work

- Pipeline B: finish the masked re-extraction (needs go-ahead, about 190 GB), rerun `series.py`, refit with 2006 to 2020, add a basin term reported both ways, report POD, FAR, and HSS, test masking transitioning tropical cyclones.
- Precursors: resolve the collision worklist (39 pairs: 31 sequential, 8 concurrent); recall on fast deepeners is not at parity and tuning was stopped deliberately; tropical and post-tropical systems are not parsed.
- The unexplained +2.6 residual in the explosive share (commit `a9be84c`).
- Hard-coded `/home/user/awips-tools/...` paths in `research/era5/` pilot scripts (`compare.py`, `envdemo.py`, `matchmonth.py`, `mismatch.py`, `sample.py`, `tip.py`). They will not run from a clone of this repository as written.

## Fixes wanted in the sheet

Corrections that belong in the spreadsheet, because a CSV-only fix is overwritten by the next fetch.

- 38 suspected position errors in `data/hf_lows/archive_position_suspects.csv` (research branches), including 8 first hurricane-force fixes the tracker refuses to anchor on.

## Cutover from awips-tools

`awips-tools` still holds its copy of these files and still publishes `/hf-lows/` on the public site. Until the steps below are done, **`awips-tools` is the copy that publishes** and changes made only here do not reach the public site.

1. Merge `harness-setup` here. That turns on the Pages preview for this repository and the checks in `pages.yml`.
2. In `awips-tools`: remove `hf-lows` from `OWN_FOLDERS` and the archive build from `site_publish.yml`, so that workflow copies `/hf-lows/` forward as it does other projects' folders.
3. Here: set the secret `WEB_TOKEN` and the variable `PUBLIC_SITE_REPO`. `publish-site.yml` then publishes `/hf-lows/`. Do not do this before step 2, or two repositories write the same folder.
4. In `awips-tools`: remove the moved files and leave a pointer to this repository. Close pull requests 80, 81, and 82 there with a link to the branches here.
5. Production: the NOAA copy is made by hand with `tools/publish.py` from a checkout. Point that checkout at this repository.

One deliberate difference: `publish-site.yml` leaves `data/qc-report.txt` out of the public copy, as `tools/publish.py` does for production. The `awips-tools` workflow copied it.
