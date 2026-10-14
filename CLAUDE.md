# Working in this repository

The archive of hurricane-force extratropical lows (North Atlantic and North Pacific), the site that presents it, and the research built on it. Forecasters and a public page depend on the numbers here, and several work sessions run against it in parallel, so the rules below are about not publishing by accident, not duplicating another session's work, and not letting a number outlive the evidence for it.

## Every session, first

1. Read the ledger as it stands on `main`: `git fetch --all`, then `git show origin/main:STATUS.md`. It says which branches are in flight, what exists, what does not, and which results are gated. The copy of `STATUS.md` on a research branch can be older; `main` holds the one that counts.
2. Compare the branch heads with the table in the ledger. The Head column is each branch's head when the ledger was last updated, so a difference means something has landed since. Look at what it was, and bring the ledger up to date before doing anything else. Git is right when the two disagree.
3. Run `git config core.hooksPath .githooks` if the session did not do it for you.
4. For anything beyond a small fix, load the `hf-harness` skill (`.claude/skills/hf-harness/SKILL.md`) and follow it.

## Wait for Jason's explicit go-ahead

These two never run on a general instruction, a plan approval, or your own initiative. Do the preparation, say exactly what you are about to run, and stop.

- **Publishing to production.** `tools/publish.py --deploy <path>` for any path outside this repository, and anything else that writes to the NOAA web root. (Regenerating the committed `flat/` copy with `--deploy flat --flat` is not publishing and needs no go-ahead.)
- **Large ERA5 pulls.** Any extraction expected to stream more than 50 GB. State the estimate, how you got it, and whether an existing cache already covers part of it.

## Never

- Commit raw ERA5 fields, `data/hsf_cache/`, `data/hf_lows/hsf_lows.csv`, per-season extraction caches, pip wheels, or anything over 10 MB.
- Commit `tools/publish.local.json`, the export token, or `.clasp.json`.
- Use `git add -A` or `git add .`. Stage files by name; regenerable intermediates sit next to committed data.
- Edit `docs/data/` or `flat/` by hand. Both are generated.
- Reorder fields in the site payload's wire format. Append only.

## Every session, last

- Update `STATUS.md` on `main`: branch heads, what now exists, any gate you opened or closed, anything left unfinished. When your work is on a research branch, the ledger update is a small pull request of its own against `main`.
- If your change makes a sentence in a README or docstring false ("there is no ...", "not yet ..."), fix that sentence in the same commit.
- Put the numbers in the commit body. Work that is not ready is committed as `WIP:` with what holds it back, not left in scratch space.
- Turn every correction into something durable. A correction is anything that made you change work you had treated as done: Jason said so, a verifier's number did not match, a reviewer called a claim overstated, a check failed for a real reason. Add one file under `failures/` and convert it, preferring a test or hook over a written rule. `failures/README.md` has the format.
- If that class of failure is already in `failures/`, fixing the instance is not enough. Open a pull request that changes how the work is done and title it "Upgrade for review:" so Jason sees it.
- Run `python3 tests/invariants/run.py` before you open or merge a pull request. A NEW violation gets fixed, not added to the known list.

## Things that are easy to get wrong here

- The spreadsheet is the source of the archive. A correction made only in `data/hf_lows/*.csv` is overwritten by the next fetch, so log it under "Fixes wanted in the sheet" in `STATUS.md`.
- The ERA5 record is a proxy. Call it one wherever it appears, and never present it as validated before the archive begins.
- Two ERA5 pipelines exist with different definitions. Say which one a number came from.
- A season runs 1 June to 31 May. The independent sample for any climate-index question is seasons, not storms.
- A script that opens `/home/...` or `/tmp/...` runs once. Find the repository from the script's own location, and read committed files.
