---
id: 2026-10-08-verifier-scripts-hardcode-clone-path
date: 2026-10-08
class: environment-assumption
caught_by: thread
became: test
became_ref: tests/invariants/run.py (no-session-paths); .claude/skills/hf-result-closeout/SKILL.md (step 1)
repeat_of: 2026-10-08-hardcoded-home-paths
workflow_change: hf-result-closeout step 1 says how a verification script must find its inputs; the no-session-paths check is run by the Checks workflow added in the same pull request as this entry
---

**What happened.** Fourteen verification scripts written on 2026-10-08 under `research/era5/*/verify/` open files by an absolute path into the sandbox that ran them (`/home/claude/hf-low/...`). One also reads intermediates from an ignored `work/` folder, and one reads from `/tmp/`. The verifications were real when they ran, but the scripts cannot be run again as committed, which is most of the reason for keeping them. Three more were committed the same afternoon (`02c6454`), after this entry was written and before it or its check had been merged, so no session had been told otherwise. They are counted here, not as a further repeat.

**How it was caught.** By running the new path check by hand against the research branch while it was being written. The same mechanism had been logged once already for six pilot scripts.

**What it was turned into.** Second failure of this class, so the workflow changes. The close-out skill now says how a verification script must locate its inputs, and the check fails any new script that opens a path under a home directory or `/tmp/`. The check does not test that a script reads committed files only; that part is a rule for the author and a question for the reviewer. The existing scripts are listed in `tests/invariants/known_violations.json` for a thread to fix.
