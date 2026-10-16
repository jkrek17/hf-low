---
id: 2026-10-16-goa-readme-causal-headings
date: 2026-10-16
class: claim-exceeds-evidence
caught_by: reviewer
became: skill
became_ref: .claude/skills/hf-harness/SKILL.md (step 4, quoted-numbers list before the claim review; wording of association headings; placebo counts)
repeat_of: 2026-10-15-readme-figures-not-in-result-files
workflow_change: a quoted-numbers list (number, result file, column) is written before the claim review and the reviewer runs before the pull request is opened; association headings use same-time wording; shuffle results are quoted as counts
---

**What happened.** The first Gulf of Alaska jets README (`b96b45a`) had causal or superlative headings ("decides how strong it gets", "matters most"), called the 54 kt proxy class "storm force" although open-water buoys read a median 44 to 48 kt when ERA5 showed it, said "no shuffle of 1,000 matched it" for a placebo that matched once, said bootstrap ranges exceeded 100% in "three bands" (it was 11), quoted per-SD units and an ASCAT figure that were in no result file, and left out that two bolded Q2 leads fail the within-season shuffle.

**How it was caught.** The hf-reviewer claim review, run after the verifier and before the pull request was opened. The verifier's recomputation had matched every table number.

**What it was turned into.** A quoted-numbers list and reviewer-before-PR step in hf-harness, with wording rules for association headings and shuffle counts. SDs and placebo counts moved into `q2_tests.csv` and `audit.csv`; README reworded in `c8a809d`. Third claim-exceeds-evidence entry after 2026-10-15-readme-figures-not-in-result-files, whose step (figures from result files) was followed for the tables but not for units or borrowed figures.
