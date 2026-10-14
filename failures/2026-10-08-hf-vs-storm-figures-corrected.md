---
id: 2026-10-08-hf-vs-storm-figures-corrected
date: 2026-10-08
class: claim-exceeds-evidence
caught_by: thread
became: skill
became_ref: .claude/agents/hf-reviewer.md; .claude/skills/hf-result-closeout/SKILL.md (step 2)
repeat_of: 2026-10-08-gust-headline-walked-back
workflow_change: hf-reviewer subagent, called from hf-result-closeout step 2 before any merge or ledger entry
---

**What happened.** After the HF versus storm-force composites were merged and entered in the ledger, the depth-matched 48-kt radius figures had to be corrected in `STATUS.md`, in the README, and in the script that builds the result's page, and limits on latitude and ocean fraction added (`88edbb9`, `a93e795`). The uncorrected ledger entry was on `main` for about six minutes.

**How it was caught.** By the thread that produced the result, after its ledger entry was already on `main`.

**What it was turned into.** This is the second time a statement reached the record ahead of its evidence, so the workflow changes: a result's claims and its ledger text are read against the result files by a separate high-effort reviewer before the merge, not after.
