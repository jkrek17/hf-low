---
id: 2026-10-16-goa-followup-verdict-words
date: 2026-10-16
class: claim-exceeds-evidence
caught_by: reviewer
became: skill
became_ref: .claude/skills/hf-harness/SKILL.md (step 4, registered verdict wording)
repeat_of: 2026-10-16-goa-readme-causal-headings
workflow_change: the README states the pre-registered reading in its registered words ("partly", "cannot tell"), with a minimum detectable effect beside any null, and every figure borrowed from a sizing step has a committed script and output
---

**What happened.** The first draft of the Gulf of Alaska follow-up README said "mostly the gust formula" where the registered rule gave "partly", answered "No." where the rule gave "cannot tell", called a forecast-gust-versus-proxy agreement "usable", gave ASCAT counts that no committed file held, wrote 22 season blocks for 21, and left several deviations out of the pre-registration log.

**How it was caught.** hf-reviewer, before the pull request, after a verifier had matched every point estimate.

**What it was turned into.** A rule in hf-harness: quote the registered verdict in its own words, put a minimum detectable effect or interval bound beside any null, and commit the sizing script and output for any borrowed figure. Same class, same day as the entry it repeats, whose quoted-numbers list covered tables but not verdict wording or sizing figures.
