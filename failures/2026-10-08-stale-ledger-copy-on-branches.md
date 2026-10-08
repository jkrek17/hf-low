---
id: 2026-10-08-stale-ledger-copy-on-branches
date: 2026-10-08
class: stale-statement
caught_by: thread
became: rule
became_ref: CLAUDE.md ("Every session, first", step 1)
repeat_of: 2026-10-07-stale-readme-no-gust-index
workflow_change: CLAUDE.md step 1 reads the ledger with `git show origin/main:STATUS.md` (6e1f4e6)
---

**What happened.** `STATUS.md` records the head of every research branch. Once `main` was merged into a research branch, that branch carried a copy that went stale as soon as any branch moved, including itself. A session reading the ledger from its own branch would have compared git against an old table.

**How it was caught.** While planning the first merges after the split, before any session had been misled.

**What it was turned into.** Sessions read the ledger from `origin/main` whatever branch they are on, and update it there with a ledger-only pull request (`6e1f4e6`). Logged after the fact as the second stale-statement entry: the change was made on the day, before this library existed.
