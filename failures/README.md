# Failure library

One file per thing that went wrong in this project and what it was turned into. The point is that a correction made once should not have to be made again: it should end up as a test, a hook, a step in a skill, or a gate in the ledger, and this folder is the record of which.

Read this before logging an entry. `python3 tools/failures.py tally` prints the running counts.

## When to add an entry

Add one when anything made you change work you had treated as done:

- Jason corrected it.
- A verifier's recomputation did not match, or a reviewer judged a claim overstated or unsupported.
- A check or test failed for a reason that is not a typo.
- You found, after the fact, that a committed statement or number was wrong.

A check that fails while you are still working, and that you fix before anyone relies on the work, is the system doing its job. It does not need an entry.

## The two rules

1. **Every entry is converted.** Prefer the form that enforces itself, in this order: a **test** in `tests/invariants/` or a test suite, a **hook**, a step in a **skill** (anything under `.claude/skills/` or `.claude/agents/`), a **rule** in `CLAUDE.md`, a **gate** in `STATUS.md`. A written rule depends on a session reading it; a failing check does not. `became: none` is allowed only while the conversion is in a pull request, and the tally lists it as open.
2. **The second failure of a class changes the workflow.** If the class already has an entry, fixing the instance is not enough. Open a pull request that changes how the work is done (a new check, a new step, a new role) and mark it "Upgrade for review" so Jason sees it. Set `repeat_of` to the earlier entry and `workflow_change` to what changed. `python3 tools/failures.py check` fails when a repeat has neither.

An entry may be logged after the fact, for something already fixed. Say so in it, so the record does not suggest the fix was a response to the library.

## Format

A file named `YYYY-MM-DD-short-slug.md`, so that two threads logging at once never touch the same file. The header is plain `key: value` lines between `---` markers:

```
---
id: 2026-10-08-short-slug          (same as the file name, without .md)
date: 2026-10-08                   (the day it happened, US Eastern)
class: claim-exceeds-evidence      (one of the classes below, or a new one)
caught_by: jason                   (jason | verifier | reviewer | ci | thread)
became: test                       (test | hook | skill | rule | gate | none)
became_ref: tests/invariants/run.py (seasons-complete)
repeat_of:                         (id of the first entry in this class, if any)
workflow_change:                   (what changed in the workflow, if a repeat)
---
```

Then three short paragraphs: **What happened.** **How it was caught.** **What it was turned into.** Cite commits. Write for someone who was not there.

## Classes

Reuse a class when the mechanism is the same, even if the subject differs. Add a class when none fits, and add it to this list in the same pull request.

| Class | The mechanism |
|---|---|
| `parallel-duplication` | Two threads did the same work without knowing about each other. |
| `stale-statement` | A document, ledger copy, or docstring kept saying something that had stopped being true, and was read as current. |
| `selection-on-outcome` | A sample, weight, filter, or test set was chosen or tuned in a way that depends on the result. |
| `claim-exceeds-evidence` | A headline, commit subject, README, or ledger entry said more than the result files support, or carried a figure they do not contain. |
| `position-or-era-artefact` | Something that varies with location, era, or recording practice was read as signal, or as a failure of the method. |
| `incomplete-sample` | An analysis ran on fewer seasons, basins, or events than it stated. |
| `unsafe-commit` | Something that should never enter history was committed or nearly committed. |
| `wrong-target` | Work went to the wrong repository, branch, or environment. |
| `environment-assumption` | Code depended on a path, cache, or machine state that exists only in the session that wrote it. |

## Where entries go

Entries live on `main`. A thread working on a research branch adds its entry in the same ledger pull request that updates `STATUS.md`.
