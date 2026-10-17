---
id: 2026-10-17-goa-kept-fields-not-kept
date: 2026-10-17
class: claim-exceeds-evidence
caught_by: thread
became: skill
became_ref: .claude/skills/hf-harness/SKILL.md (claim review: "kept at" needs an ls in the same session)
repeat_of: 2026-10-17-goa-hrrr-readme-overclaims
workflow_change: any statement that data are kept at a shared-folder path is checked with ls in the session that writes it; found while planning a follow-up that needed the HRRR grids
---

**What happened.** The first Gulf of Alaska README, the follow-up ledger entry and the HRRR README said the ERA5 window, the GFS forecast fields and the HRRR fields were kept in `/mnt/project-files/goa-jets/`. Only the result tables were there. Logged after the fact: found when a follow-up study needed the HRRR grids and had to pull again (2.7 GB for the gust alone).

**How it was caught.** The thread, by listing the folder before planning a pull.

**What it was turned into.** A step in the hf-harness claim review: write "kept at" only after an `ls` in the same session. The sentences in both READMEs and the ledger are corrected.
