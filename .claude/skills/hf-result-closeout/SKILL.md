---
name: hf-result-closeout
description: Use when an hf-low research result is finished, to get its numbers independently verified, its PR merged under the standing authority, the STATUS.md ledger updated on main, and the result forwarded to the Science Q&A document thread.
---

# Closing out a research result (hf-low)

`hf-harness` steps 4 and 5 give the checks and the ledger rules; this skill is the order of operations and the standing permissions. Do not repeat them, follow them.

## 1. Verify with a fresh agent

- Start a Sonnet verifier (coding and verifier subagents are Sonnet) with the harness brief shape (ROLE/GOAL/CONTEXT/SCOPE/OUT OF SCOPE/DONE WHEN/RETURN/IF BLOCKED).
- Give it the committed inputs and the pre-registration, not your code, README or notes. Ask it to reimplement or recompute each quoted number.
- The verification script finds the repository from its own location (`os.path.dirname(os.path.abspath(__file__))` and up), never an absolute path such as `/home/claude/hf-low/` or `/tmp/`, and reads committed files only. If it needed an uncommitted intermediate, `VERIFICATION.md` says which and how to rebuild it. A script nobody else can run is a record that a check happened, not a check.
- Save its script and output beside the work (`verify/verify_numbers.py`, `VERIFICATION.md`). `VERIFICATION.md` lists what matched, any difference with its cause, and **what was not independently checked** (for example bootstrap intervals).
- Every quoted number is either recomputed or listed as unchecked. Tell Jason which.
- Finish the README: the answer in plain words first, pipeline A or B, proxy wording, the post hoc items labelled, and the count of held-out looks.

## 2. Review the claims with `hf-reviewer`

After the numbers are verified and the README is finished, and before any merge or ledger entry:

- Call the `hf-reviewer` subagent (Opus, high effort). Give it the topic folder, the commit subject or pull request title, and the exact `STATUS.md` text you intend to add. Do not give it your notes.
- It returns a verdict per claim: supported, overstated, unsupported, or not checkable. Change every overstated or unsupported statement, or state the disagreement in the README and tell Jason.
- A ledger entry is checked before it reaches `main`, not corrected after. This step exists because that happened twice (`failures/`, class `claim-exceeds-evidence`).
- If the reviewer or the verifier made you change something, that is a correction: log it in `failures/` (step 5).

## 3. Commit hygiene

Stage files by name, never `git add -A` or `git add .`. Nothing over 10 MB; no raw ERA5, `data/hsf_cache/`, `hsf_lows.csv` or per-season caches. Put the numbers in the commit body. A result that is not ready is `WIP:` with what holds it back. Fix any README or docstring sentence your result made false. Draft PR first; PR body starts with the project attribution block. Run `python3 tests/invariants/run.py`; the Checks workflow runs it and the test suites on the pull request.

## 4. Merge authority (standing, from Jason)

- Finished, verified and reviewed research PR into `claude/exciting-fermat-8vcvgq`, with the Checks workflow green: merge it yourself and say so in the thread. This replaces the harness line about leaving PRs as drafts. If Checks is red, fix the cause; do not merge around it.
- A pull request titled "Upgrade for review:" (a change to a skill, a check, a hook, or `CLAUDE.md` made because a failure repeated) is never self-merged. It waits for Jason.
- Research or site code into `main`, and any live deploy or write to the NOAA web root: wait for Jason. `tools/publish.py --deploy <path outside repo>` needs his explicit go-ahead.
- Large ERA5 pull (over 50 GB): ask first.

## 5. Ledger PR to main

A separate small PR containing `STATUS.md` and, when there was a correction, one new file under `failures/` (format in `failures/README.md`), and nothing else: branch head, what now exists (with the verified headline numbers), gates opened or closed, open items and any "Fixes wanted in the sheet". Self-merge once the numbers are verified. If a permission check blocks the merge, retry once, then leave the PR open or write to `/mnt/project-files/ledger-pending/`; never work around it.

## 6. Forward and report

- Send the finished result to the "Science Q&A document" thread (via the coordinator): question, pre-registered vs post hoc status, method, how to reproduce, and the result with nulls given equal prominence.
- Reply once to Jason when the goal is finished or stuck: the answer in plain words, the numbers' verification status, what the reviewer made you change, and anything waiting on him. Decisions go as one batched card with a recommendation.
- Routine choices: record the assumption in `STATUS.md` and keep going; do not ask.
