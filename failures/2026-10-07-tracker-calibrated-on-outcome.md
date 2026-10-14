---
id: 2026-10-07-tracker-calibrated-on-outcome
date: 2026-10-07
class: selection-on-outcome
caught_by: thread
became: skill
became_ref: .claude/skills/hf-harness/SKILL.md (step 4, question 1); fix in 7f89e3d
repeat_of: 
workflow_change: 
---

**What happened.** The precursor tracker's pressure cost was calibrated on the mature hurricane-force segment of the archive, where deepening has largely finished. Recall on storms that fell 16 hPa or more into their first hurricane-force fix was 16%, against about 72% for the rest (`36e6e35`). A guard against one bias (preferring candidates that bomb) had produced a worse one, and still selection on the outcome.

**How it was caught.** By the stage's own validation, once recall was broken down by how fast the storm was deepening.

**What it was turned into.** Recalibration on the regime being tracked (`7f89e3d`, recall 50%), the withheld-output convention, and question 1 in `hf-harness`: was anything chosen using the outcome, and break recall down by the quantity that matters.
