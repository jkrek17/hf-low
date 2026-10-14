---
id: 2026-06-24-missing-atlantic-seasons
date: 2026-06-24
class: incomplete-sample
caught_by: jason
became: test
became_ref: tests/invariants/run.py (seasons-complete)
repeat_of: 
workflow_change: 
---

**What happened.** An ENSO, NAO and MJO compositing analysis ran with four Atlantic seasons (2022-23 to 2025-26) missing from its input. Nothing failed; the composites were simply computed on 17 seasons while the text said 21.

**How it was caught.** Jason noticed the season count during a final audit of the session. This happened in a chat session in June 2026, before the work had a repository, so nothing here can confirm the details.

**What it was turned into.** A rule in `hf-harness` (confirm every season is present for every basin before a composite or regression), and now a check that fails when any season from the record start to the latest season has no events for a basin in the site payload. The check is narrower than the failure: an analysis that assembles its own input is still covered only by the rule.
