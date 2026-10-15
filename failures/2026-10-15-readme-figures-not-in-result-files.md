---
id: 2026-10-15-readme-figures-not-in-result-files
date: 2026-10-15
class: claim-exceeds-evidence
caught_by: reviewer
became: skill
became_ref: .claude/skills/hf-harness/SKILL.md (step 4, figures in plain-words text must be in a result file)
repeat_of: 2026-10-14-ascat-independence-overstated
workflow_change: the harness verify step now requires every figure in a README's plain-words text to appear in a committed result file, and unit and wording checks (one SD in physical units, "cannot tell" for nulls) before the claim review
---

The tip jets README (hf-low PR 243) was reviewed before merge. Two figures were in no result file (episodes per month "36 to 48", which belonged to the unwidened classes, and "about 11 hPa" per SD of depth, which was another study's SD), a widened class was labelled as the original, nulls were written as "does nothing" with detectable effects near 2, and the proposed ledger text said "three setups" for 328 of 514 episodes. The verifier's recomputation matched every table number; the errors were in the prose around them. Fixed in the same PR before merge.
