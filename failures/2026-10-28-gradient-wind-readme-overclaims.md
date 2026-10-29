---
id: 2026-10-28-gradient-wind-readme-overclaims
date: 2026-10-28
class: claim-exceeds-evidence
caught_by: reviewer
became: skill
became_ref: .claude/skills/hf-harness/SKILL.md (step 4, cross-table and verdict-word rule)
repeat_of: 2026-10-17-goa-station-readme-overclaims
workflow_change: the claim reviewer's checklist in hf-harness step 4 now lists four things it was catching by hand: a sentence that compares two tables names both samples; a verdict word is compared with the plan's own rule (including when the interval sits near the bound) before it is written; a table built from gust says plainly that gust is not sustained wind and is not a warning category; every post hoc item is logged under "Changes after seeing results" in the plan before the README is written
---

**What happened.** The first draft of the gradient-to-wind README (hf-low PR 282) was matched number for number by two fresh agents, then a claim reviewer found 20 sentences to fix. Examples: a ratio of 0.56 quoted beside a 34 kt and 54 kt pair that gives 0.63 (two samples mixed); "direction as expected" for a test that was a "yes" under the plan's own rule, and "not a yes" for a result outside every verdict band; a core versus non-core lead taken from one band and one bin and stated as a finding; the post hoc items labelled in the README but not logged in the plans; a 15% skill gain presented without saying it used the analysis pressure and that ERA5 assimilates the ASCAT outcome; "by construction" for something untested; a gust table with no warning that gust is not the sustained wind warnings are issued on. The verifier had also caught five wording errors in the first pass (for example "most tight gradients are in cores" when core cells were 22% at 6-8 hPa per 100 nm).

**How it was caught.** The first verifier (wording and one size claim), the second verifier (P1b restated), and the claim reviewer, all before merge.

**What it was turned into.** A line in the hf-harness skill's step 4 (the quoted-numbers paragraph), so the reviewer and the author check the same four things. No lint: these need judgement. This is the ninth entry of the class and a repeat of the 17 October README overclaims; the existing lint (`overclaim-words`) did not catch any of these, so the pull request adds the checklist and is marked "Upgrade for review" for Jason to decide whether a stronger enforcement is wanted (for example a plan-versus-README verdict-word diff).
