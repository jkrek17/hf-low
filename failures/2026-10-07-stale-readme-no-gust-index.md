---
id: 2026-10-07-stale-readme-no-gust-index
date: 2026-10-07
class: stale-statement
caught_by: thread
became: test
became_ref: tests/invariants/run.py (stale-claims) with tests/invariants/claims.json
repeat_of: 
workflow_change: 
---

**What happened.** `research/era5/README.md` went on saying "There is no calibrated gust index" and "no extended event record" after both had been committed beside it. A reader starting from the README would have concluded the work did not exist.

**How it was caught.** Found when the repository was read end to end for the split.

**What it was turned into.** First a step in `hf-harness` (search for sentences your change made false). That step depends on a session remembering to run it. The README gained a "What exists now" section on 2026-10-08 (`3817ae8`), which removed the sentence, and the sentence is now also guarded by a check: `claims.json` pairs a sentence with the file whose existence makes it false, so it cannot come back unnoticed. Add a pair there whenever a "there is no X" sentence is written about something a later commit could create.
