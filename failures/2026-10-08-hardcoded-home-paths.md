---
id: 2026-10-08-hardcoded-home-paths
date: 2026-10-08
class: environment-assumption
caught_by: thread
became: test
became_ref: tests/invariants/run.py (no-session-paths)
repeat_of: 
workflow_change: 
---

**What happened.** Six pilot scripts under `research/era5/` open `/home/user/awips-tools/docs/data/hf-lows.json`. They ran in the session that wrote them and in no other checkout.

**How it was caught.** Found when the repository was read end to end for the split.

**What it was turned into.** A check that fails on an absolute path into a home directory or `/tmp/` in tracked code. The six scripts are listed as known violations until they are fixed or retired.
