---
id: 2026-10-14-session-paths-in-gfs-test
date: 2026-10-14
class: environment-assumption
caught_by: ci
became: hook
became_ref: .githooks/pre-commit (runs tests/invariants/run.py)
repeat_of: 2026-10-08-hardcoded-home-paths
workflow_change: the pre-commit hook now runs tests/invariants/run.py, so a NEW no-session-paths violation blocks the commit locally instead of surfacing as a red CI check after merge
---

**What happened.** Three scripts of the GFS test of the storm-pattern terms (research PR 206: `extract_gfs.py`, `gfs_fixes.py`, `score_gfs.py`) defaulted their work directory to `/tmp/claude-0/...`. They were merged into the research branch, and the `invariants` check (`no-session-paths`) went red on the next pull request. The same mechanism as 2026-10-08-hardcoded-home-paths: a script that runs only in the session that wrote it.

**How it was caught.** The `invariants` CI check on research PR 210, noticed by the summary-page thread. The session that wrote the scripts did not run `tests/invariants/run.py` before merging (the CLAUDE.md rule asking for it was added after the merge).

**What it was turned into.** The scripts were fixed in research PR 212 (work directory defaults to `gfs/work/` or `GFS_WORK`; shared folder via `STORM_PATCH_DIR`). The pre-commit hook now runs the invariants (0.2 s), so a new violation fails at commit time. The check itself already existed; what changed is when it runs.
