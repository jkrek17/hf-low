---
id: 2026-10-15-readme-figures-not-in-result-files
date: 2026-10-15
class: claim-exceeds-evidence
caught_by: reviewer
became: skill
became_ref: .claude/skills/hf-harness/SKILL.md (step 4, figures in plain-words text must be in a result file)
repeat_of: 2026-10-08-hf-vs-storm-figures-corrected
workflow_change: the harness verify step now requires every figure in a README's plain-words text to appear in a committed result file, and unit and wording checks (one SD in physical units, "cannot tell" for nulls) before the claim review
---

**What happened.** The tip jets README (hf-low PR 243) had two figures in no result file (episodes per month "36 to 48", which belonged to the unwidened classes, and "about 11 hPa" per SD of depth, which was another study's SD). A widened class was labelled as the original, nulls were written as "does nothing" with detectable effects near 2, and the proposed ledger text said "three setups" for 328 of 514 episodes.

**How it was caught.** The hf-reviewer claim review before merge. The verifier's recomputation had matched every table number; the errors were in the prose around them.

**What it was turned into.** A step in the hf-harness verify section: figures in plain-words text come from a committed result file, nulls are written as "cannot tell, detectable about N", and a widened class is named wherever it appears. Fixed in PR 243 before merge. Second claim-exceeds-evidence entry after 2026-10-08-hf-vs-storm-figures-corrected (and the ASCAT independence wording on main).
