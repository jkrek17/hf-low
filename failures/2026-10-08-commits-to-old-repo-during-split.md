---
id: 2026-10-08-commits-to-old-repo-during-split
date: 2026-10-08
class: wrong-target
caught_by: thread
became: gate
became_ref: STATUS.md ("Cutover from awips-tools": "Do all hurricane-force work here from now on")
repeat_of: 
workflow_change: 
---

**What happened.** While the hurricane-force files were being moved out of `awips-tools`, two sessions kept committing there. Three commits landed after the copy was taken (`c9dc994`, `5f05669`, `747d87d` as they are numbered here).

**How it was caught.** The old repository's branch heads were compared with the copy after the move (`4527a7b`).

**What it was turned into.** The three commits were carried over with the same filter, and the ledger says to do all hurricane-force work here. That is a written instruction only. `awips-tools` still holds its copy of the files and still publishes `/hf-lows/` until the cutover steps in the ledger are done, so the same mistake remains possible until then.
