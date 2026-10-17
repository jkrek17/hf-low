---
id: 2026-10-17-gulf-stream-round2-verdicts-and-headlines
date: 2026-10-17
class: claim-exceeds-evidence
caught_by: reviewer
became: gate
became_ref: STATUS.md (Gulf Stream round 2 entry: a README "no" verdict needs a power line for that test in the topic's power file); PR 264 stays the process upgrade
repeat_of: 2026-10-08-gust-headline-walked-back
workflow_change: PR 264 (verdict-wording rule, open for Jason) is the upgrade. A candidate invariant, not yet built: fail when a README "Plan verdict" cell says "no" with no matching line in the topic's power file. Offered to Jason, not built without his word.
---

**What happened.** The first README for Gulf Stream round 2 had the numbers right (a fresh agent matched 15 of 15 within 0.04 points) but several words past the evidence: B2e called "no" on power borrowed from round 1 after its estimate was seen; B4e called "can't tell" where the plan's rule says "no, well powered, marginal"; the two forecaster headlines said "cold-air-outbreak effect" and "not stronger over the front than elsewhere" (the plan verdict is can't tell, opposite sign); "the difference is zero", "explains", "not a latitude artefact", "not a method artefact" and "if anything the model is too strong" were stated as findings; the registered A1m row was missing.

**How it was caught.** hf-reviewer, after the verifier. All fixed in the same pull request (hf-low PR 271) before merge.

**What it was turned into.** A gate in the ledger entry and the pending PR 264; the lint `overclaim-words` passed this README, so it does not yet catch "explains", "is zero" or "is a ... effect".
