---
name: hf-result-closeout
description: Use when an hf-low research result is finished, to get its numbers independently verified, its PR merged under the standing authority, the STATUS.md ledger updated on main, and the result forwarded to the Science Q&A document thread.
---

# Closing out a research result (hf-low)

`hf-harness` steps 4 and 5 give the checks and the ledger rules; this skill is the order of operations and the standing permissions. Do not repeat them, follow them.

## 1. Verify with a fresh agent

- Start a Sonnet verifier (coding and verifier subagents are Sonnet) with the harness brief shape (ROLE/GOAL/CONTEXT/SCOPE/OUT OF SCOPE/DONE WHEN/RETURN/IF BLOCKED).
- Give it the committed inputs and the pre-registration, not your code, README or notes. Ask it to reimplement or recompute each quoted number.
- Save its script and output beside the work (`verify/verify_numbers.py`, `VERIFICATION.md`). `VERIFICATION.md` lists what matched, any difference with its cause, and **what was not independently checked** (for example bootstrap intervals).
- Every quoted number is either recomputed or listed as unchecked. Tell Jason which.
- Finish the README: the answer in plain words first, pipeline A or B, proxy wording, the post hoc items labelled, and the count of held-out looks.

## 2. Commit hygiene

Stage files by name, never `git add -A` or `git add .`. Nothing over 10 MB; no raw ERA5, `data/hsf_cache/`, `hsf_lows.csv` or per-season caches. Put the numbers in the commit body. A result that is not ready is `WIP:` with what holds it back. Fix any README or docstring sentence your result made false. Draft PR first; PR body starts with the project attribution block.

## 3. Merge authority (standing, from Jason)

- Finished and verified research PR into `claude/exciting-fermat-8vcvgq`: merge it yourself and say so in the thread. This replaces the harness line about leaving PRs as drafts.
- Research or site code into `main`, and any live deploy or write to the NOAA web root: wait for Jason. `tools/publish.py --deploy <path outside repo>` needs his explicit go-ahead.
- Large ERA5 pull (over 50 GB): ask first.

## 4. Ledger PR to main

A separate small PR containing `STATUS.md` only: branch head, what now exists (with the verified headline numbers), gates opened or closed, open items and any "Fixes wanted in the sheet". Self-merge once the numbers are verified. If a permission check blocks the merge, retry once, then leave the PR open or write to `/mnt/project-files/ledger-pending/`; never work around it.

## 5. Forward and report

- Send the finished result to the "Science Q&A document" thread (via the coordinator): question, pre-registered vs post hoc status, method, how to reproduce, and the result with nulls given equal prominence.
- Reply once to Jason when the goal is finished or stuck: the answer in plain words, the numbers' verification status, and anything waiting on him. Decisions go as one batched card with a recommendation.
- Routine choices: record the assumption in `STATUS.md` and keep going; do not ask.
