---
id: 2026-10-07-duplicate-era5-pipelines
date: 2026-10-07
class: parallel-duplication
caught_by: thread
became: skill
became_ref: .claude/skills/hf-harness/SKILL.md (step 1, "Is someone already on this?"); STATUS.md ("How A and B differ")
repeat_of: 
workflow_change: 
---

**What happened.** Two sessions built ERA5 proxy pipelines at the same time: pipeline A in `0ed86eb`, and pipeline B in the line of commits running from `c904636` to `c9dc994`. They differ in domain, gust radius, low detection, months, labels and calibration window. The harness skill records that neither session knew about the other until both were mostly done; no commit says so directly.

**How it was caught.** Partway through, the pipeline B session read a skill written for the same task and found its own extraction lacked an ocean mask (`20b8656`). The two pipelines were first set side by side in the ledger written at the repository split (`57e9d0d`).

**What it was turned into.** The ledger (`STATUS.md` on `main`, one row per thread, and the table "How A and B differ") and the first question in `hf-harness`. There is no automated check.
