---
id: 2026-10-14-ascat-independence-overstated
date: 2026-10-14
class: claim-exceeds-evidence
caught_by: reviewer
became: none
became_ref: open: no conversion yet; question for Jason (see workflow_change)
repeat_of: 2026-10-08-gust-headline-walked-back
workflow_change: not yet made. hf-result-closeout step 2 (claim review before any merge) already exists; PR 214 was merged before the reviewer ran. A mechanical guard (for example a check that a research PR touching a README with a verdict also carries a reviewer line) needs an "Upgrade for review" PR, which waits for Jason.
---

**What happened.** The ASCAT check (research PR 214, merged as `76217ce`, commit subject of `a5b8d9d` "ERA5 gust proxy tracks ASCAT winds") described ASCAT winds as independent satellite observations, said the proxy "tracks observed winds", counted P4 (a significant sensor difference) as a pass in "5 of 5", wrote "Decision 1 holds", said ASCAT "saturates" and gave "lower bounds", and called the work independently verified while the verification files were uncommitted. ERA5 assimilates ASCAT inside the same 12 h windows, so the data are a consistency check, not independent validation. The integration branch carried these statements for about 15 minutes; the ledger text on `main` had not been merged.

**How it was caught.** The advisor pass before closeout raised the independence point; the Opus claim reviewer then listed every overstated or unsupported sentence. Both ran after the merge, not before, so step 2 of `hf-result-closeout` was done out of order.

**What it was turned into.** The wording was corrected in a follow-up research PR and in the ledger entry (the commit subject cannot be edited, so the ledger says not to quote it). No self-enforcing conversion yet; this is the second `claim-exceeds-evidence` entry and the rule calls for a workflow change, which is a question for Jason (above).
