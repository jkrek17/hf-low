---
id: 2026-10-08-gust-headline-walked-back
date: 2026-10-08
class: claim-exceeds-evidence
caught_by: thread
became: gate
became_ref: STATUS.md ("Gust-based event counts before 2001"); recomputed in c55e74d
repeat_of: 
workflow_change: 
---

**What happened.** Commit `c9dc994` is titled "A gust-based criterion cannot be carried back before 2001". Its measurements stand: gust rises at fixed storm depth before 2001. Its headline went further than they do. Within the hour the same session said the conclusion was wrong at event level, on a reconciliation that existed only in its conversation. For about 25 minutes the ledger on `main` carried the commit's conclusion as a gate and as its first decision (`4527a7b` to `db57b41`). The recomputation that followed found that neither "cannot be carried back" nor the walk-back's "counts are fine" is supported, and that the reconciliation's own test had 12 to 15% power (`c55e74d`).

**How it was caught.** The session reported it when it was closed out.

**What it was turned into.** A ledger gate that separates what is measured from what is unsettled, and a recomputation from committed files. Commit subjects cannot be edited, so the ledger entry is the correction.
