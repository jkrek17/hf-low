---
id: 2026-10-07-large-files-nearly-committed
date: 2026-10-07
class: unsafe-commit
caught_by: thread
became: hook
became_ref: .githooks/pre-commit (4f0f6b8); tests/invariants/run.py (forbidden-files)
repeat_of: 
workflow_change: 
---

**What happened.** In one session, three large files came within one `git add -A` of the history: a 67 MB parsed intermediate (`36e6e35`), a fabricated development CSV sitting at a real data path, and about 51 MB of pip wheels (`54ca46a`).

**How it was caught.** By the session itself, before committing each time.

**What it was turned into.** Ignore rules, then a pre-commit hook that refuses staged files over 10 MB and build artefacts, enabled at session start (`4f0f6b8`). The forbidden-files check applies the same limits to files that are already tracked, so a clone without the hook enabled is still caught.
