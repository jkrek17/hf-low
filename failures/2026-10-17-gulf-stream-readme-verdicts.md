---
id: 2026-10-17-gulf-stream-readme-verdicts
date: 2026-10-17
class: claim-exceeds-evidence
caught_by: reviewer
became: skill
became_ref: .claude/skills/hf-preregistered-test/SKILL.md (Verdict words come from the plan)
repeat_of: 2026-10-17-goa-hrrr-readme-overclaims
workflow_change: the pre-registered-test skill now fixes the verdict wording to the plan's decision rules, voids power statements when a placebo fails, and asks for a one-station-at-a-time look before a group verdict
---

**What happened.** The first README draft of the Gulf Stream buoy study said "a few knots" where the pre-registered figure is +1.0 kt (the larger number was a post hoc lead), called the pooled wind result "no effect" although it averages opposite-signed groups, gave two significant opposite-sign tests no verdict, said a 2 ft wave effect "would probably have shown" although the placebo had failed and the power run did not fit its own standard errors, and put the Stream-versus-shelf split forward as a Stream effect although one Stream buoy is -7% and two shelf buoys are +7% and +10%.

**How it was caught.** hf-reviewer, before the pull request, after a verifier matched the numbers.

**What it was turned into.** A paragraph in the `hf-preregistered-test` skill (section 3). Fourth entry of this class in two days; the lint route used for the HRRR entry does not reach this wording, so the change is to the skill and to the order of work: verifier, then reviewer, then pull request.
