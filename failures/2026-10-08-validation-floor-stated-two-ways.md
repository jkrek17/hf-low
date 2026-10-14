---
id: 2026-10-08-validation-floor-stated-two-ways
date: 2026-10-08
class: stale-statement
caught_by: thread
became: gate
became_ref: STATUS.md ("Windows before 2004"); .claude/skills/hf-preregistered-test/SKILL.md (section 2, "Seasons"); tests/invariants/run.py (record-start)
repeat_of: 2026-10-07-stale-readme-no-gust-index
workflow_change: hf-preregistered-test section 2 fixes "fit and test on 2004-05 onward" for every test, and the record-start check pins RECORD_START at 2004 in the build and the payload
---

**What happened.** The year before which the proxy cannot be checked against the archive was written three ways. `criterion-result.txt`, `research/era5/README.md` and `hf_history/README.md` said 2001; `series.py` and `series-result.txt` said 2004; the first ledger said pipeline A's seasons "before 2006" do not test the threshold. The build's `RECORD_START` had been 2004 throughout. Nothing was computed wrongly because of it, but a session choosing a fit window had three answers to pick from.

**How it was caught.** Listed as a conflict in the first ledger (`57e9d0d`), then settled when the pipeline B session reported at close-out that `RECORD_START = 2004` must gate every window (`db57b41`).

**What it was turned into.** One gate in the ledger, a fixed rule in `hf-preregistered-test`, and a check that fails if the build constant or the payload stops saying 2004. Logged after the fact: the gate and the rule were written on the day, before this library existed.
