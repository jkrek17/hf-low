---
name: hf-harness
description: Operating procedure for work in the hurricane-force low archive repository. Load before any task beyond a small fix - archive data, site build, precursor recovery, ERA5 proxy record, teleconnection statistics, or resuming a branch in flight.
---

# HF low project harness

How to run a piece of work in this repository so that it does not collide with another session, does not publish by accident, and does not put a number into the record that the evidence will not carry.

`CLAUDE.md` has the standing rules and the two actions that need Jason's go-ahead. This file is the procedure and the checks. `STATUS.md` is the ledger of what is in flight. Current numbers live in `STATUS.md`, in result files, and in commit bodies, never in this file, so that this file does not go stale.

## Why this project needs a procedure

Three things have already gone wrong here, each caught late:

- Two sessions built separate ERA5 pipelines in parallel, with different domains, radii, masks, and calibration windows. Neither knew about the other until both were mostly done.
- A README went on saying "there is no calibrated gust index" after one had been committed beside it.
- A tracker was calibrated on the part of the archive where deepening had already finished, so it refused to follow exactly the storms the project exists to measure. A guard against one bias produced a worse one in the other direction.

Each step below exists to stop one of these from happening again.

## The loop

### 1. Orient

Do what `CLAUDE.md` says under "Every session, first". Then answer three questions before touching anything:

- **Is someone already on this?** Look in `STATUS.md` for a thread that owns the topic. If one exists, extend that branch. If you believe a second approach is warranted, write one line in `STATUS.md` saying why it is deliberate before you start. Parallel implementations are allowed; accidental ones are the problem.
- **What does the ledger say exists and does not exist?** Believe result files and commit bodies over README prose when they disagree, and note the disagreement for step 5.
- **Is anything you need gated?** A gate in `STATUS.md` ("not safe to compare across seasons", "flag, not a finding") binds your work too. Building on a gated result needs the gate closed first, with evidence, or the dependence stated in your output.

### 2. Frame

Write down, in the pull request description or the first commit body:

- The goal in one or two sentences.
- Which strand it belongs to (section "Checks by strand" below).
- The checks that will show it is done. Take them from the strand's list; add your own where the list does not reach.
- Any step that needs go-ahead under `CLAUDE.md`.

If the request is ambiguous and a wrong guess would cost a long extraction or change what a published number means, ask. Otherwise take the most reasonable reading and record it as an assumption.

### 3. Work

- Work on a branch. Open a draft pull request early so the thread is visible, and add it to `STATUS.md`.
- Commit often. A `WIP:` commit that says what holds the work back is a good commit. When a stage's own validation says its output is not fit to use, commit the code and the validation numbers and withhold the data file, saying so in the subject line.
- Before any expensive pull or fetch: estimate the volume, check whether a cache already covers it, make the job resumable, and send intermediates to an ignored `work/` directory.
- When you delegate to a fresh agent, tell it everything it needs; it cannot see this session. Give it the branch, the files, what it may change, and the exact shape of what to return.

### 4. Verify

Run the strand's checks. Then decide whether an independent check is owed:

**A number that will be quoted** in a README, a commit subject, a report, the site's Method tab, or `STATUS.md` **gets recomputed by an agent that did not produce it**, from the committed files, with no access to your reasoning. Give it the claim, the files, and the question "is this number what these files produce?" If you cannot start a fresh agent, recompute it yourself by a second method and say in the commit that the check was not independent.

Three questions to put to every result before it goes out. Each has a failure in this project behind it:

1. **Was anything chosen using the outcome?** A calibration sample, a cost weight, a filter, or a validation set picked or tuned after seeing the result is selection on the outcome. Say how many times a held-out test was looked at. Break recall down by the quantity you care about, since an average hides a failure concentrated in the cases that matter.
2. **Does an artefact correlate with position, era, or recording practice?** The teleconnection work measures where and when storms occur, so contamination that varies with location (land in a gust radius) or with era (archive lead fixes appearing in 2017, the warning category appearing in 2001) can manufacture or hide the signal. Check coverage and skill by era and by basin before pooling.
3. **What is the real sample size?** For any question about a climate index it is the number of seasons, and for a strongly autocorrelated index like ONI it is about one value per season. State n, give an interval or a signal-to-noise figure, and do not report a trend from a handful of points as a finding.

Report each check as passed, failed, or not checked with the reason. "Not checked" is acceptable. Reporting an unchecked number as established is not.

### 5. Close

- Update `STATUS.md`: heads, what exists, gates opened or closed, what is unfinished, what superseded what.
- Search for sentences your change made false: `git grep -n -i -E "there is no|does not exist|not yet|still starts|no .* exists"` over READMEs and docstrings in the area you touched. Fix them in the same commit.
- When a number replaces an earlier one, the commit body says which and why ("supersedes 16% in <hash>: the calibration sample was wrong").
- Leave the pull request as a draft unless Jason asked for it to be ready.
- Report to Jason briefly: what was produced, which checks passed and which were not run, assumptions made, and what is waiting on him.

## Checks by strand

### Archive data and site

| Check | Command |
|---|---|
| Build runs and reports counts | `python3 tools/build_hf_lows.py --check` |
| Rebuild site data | `python3 tools/build_hf_lows.py` |
| Site logic | `node tests/composite/test_composite.js`, `node tests/playback/test_playback.js`, `node tests/regress/test_regress.js`, `node tests/teleconnect/test_teleconnect.js` |
| Index data | `python3 tests/teleconnections/test_teleconnections.py` |
| Build and backfill (research branches) | `python3 tests/hf_lows/test_hf_lows.py` |
| Regenerate `flat/` | `python3 tools/publish.py --no-fetch --deploy flat --flat --yes`, then do not commit `flat/.awips-publish-manifest.json` |
| Page validity, as the downstream CI runs it | `npx html-validate@11.15.0 docs` |

- Every build stamps a new `generated` time and `build` object, so a plain diff always shows a change. To compare two builds, drop those two fields first. `tools/publish.py` prints a delta report that already does this.
- After changing `docs/`, regenerate `flat/` in the same commit. The Pages workflow fails when they drift, because a forecaster copies `flat/` to production by hand.
- `docs/data/qc-report.txt` is for the repository only. It is never deployed.
- Rows before the 2004-05 season are short-counted because the archive was still starting. They stay in the data and out of the default view; say so when a statistic includes them.
- A correction to an archive row belongs in the spreadsheet. If you must patch the CSV, log it in `STATUS.md` under "Fixes wanted in the sheet".

### Statistics on the archive and climate indices

- Before any composite or regression, confirm every season is present for every basin in the range you claim. A silent gap of four Atlantic seasons got into an earlier analysis this way.
- The MJO series is a longitude, not eight phases. Do not label it RMM or bin it into phases; `tools/build_teleconnections.py` explains why.
- Statistics from backfilled (pre-hurricane-force) data are not comparable across seasons while the build reports its coverage-trend flag. A rate computed on the recovered subset describes that subset, not the archive.
- Apply the three questions in step 4. Question 3 decides most of what can be said here.

### Precursor recovery (High Seas Forecast text)

| Stage | Command | Output |
|---|---|---|
| Fetch | `python3 tools/fetch_hsf.py` | `data/hsf_cache/`, ignored, resumable |
| Parse | `python3 tools/parse_hsf.py` | `data/hf_lows/hsf_lows.csv`, ignored, never committed |
| Track | `python3 tools/track_hsf.py` | `data/hf_lows/precursors.csv`, committed |
| Validate | `python3 tools/track_hsf.py --validate` | the numbers that may be quoted |
| Review collisions | `python3 tools/review_collisions.py` | a worklist; nothing is merged or dropped |

Tests: `python3 tests/fetch_hsf/test_fetch_hsf.py`, `python3 tests/hsf_parse/test_parse_hsf.py`, `python3 tests/track_hsf/test_track_hsf.py`.

- Read the block headed "LIMITS OF THE VALIDATION AND OF THE RESULT" in `tools/track_hsf.py` before quoting anything from this stage. In short: the main validation has been looked at several times and is no longer a clean held-out test, the recovered storms are not a random sample, and the fastest deepeners are under-represented.
- If you run the validation again after a design change, add to the count of looks in that docstring.
- `warn_cat` is the warning a low was filed under, which for a developing storm is what it is forecast to become. It is not an observation of hurricane-force wind at that time.
- The low-confidence tier is refused on purpose. Do not admit it to raise coverage.

### ERA5 proxy record

- Call it a proxy everywhere it appears: file headers, column descriptions, figures, commit subjects. An event in it is "a cyclone whose ERA5 fields look like the ones OPC warned for as hurricane force", which is well defined and is not the archive.
- It cannot be validated before the archive begins. How far back an index may be used at all is a gate in `STATUS.md`, and it has moved more than once as evidence came in. Check it before extending or quoting any pre-archive figure. Gust-based and depth-based quantities have not behaved the same way before the archive, so a result for one says nothing about the other.
- Two pipelines exist (`research/era5/hf_history/` and the `event_fields.py` / `criterion.py` / `series.py` set). They differ in domain, gust radius, low detection, months, labels, calibration window, and what they output. Name the pipeline with every number, and do not mix their outputs.
- Gust features are taken over ocean points only. The archive is about wind over water, and land or terrain inside the radius adds error that varies with position.
- Report skill in the terms a forecaster reads: POD, FAR, CSI, HSS, and bias, alongside AUC. Give leave-one-season-out ranges and the transfer to seasons outside the fit.
- A sum of probabilities is not an event count. A count fitted to match the archive in the calibration window is on the archive's scale only there.
- Transitioning tropical cyclones are a known source of false alarms. Say whether they are in or out.
- A re-extraction must reproduce the earlier numbers or explain every difference. If the track count changes, say why.
- Never commit raw fields or per-season caches. Commit the scripts and the small derived outputs. A cache made before a feature definition changed is invalid; say so where the derived file is described.

## Briefing a fresh agent

```
ROLE: <scout | worker | verifier>
GOAL: <the one thing to accomplish, and what it feeds>
CONTEXT: <branch, files, decisions already made, definitions it must use>
SCOPE: <what it may read and change>
OUT OF SCOPE: <what it must not touch, read, or decide>
DONE WHEN: <checkable conditions>
RETURN: <exact shape: table, verdict per claim, file path>
IF BLOCKED: stop and report what is missing; do not guess or widen scope
```

A verifier gets the claim and the files, is told not to read your notes or the pull request discussion, and is asked to recompute, not to review.
