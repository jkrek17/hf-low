---
id: 2026-10-07-land-in-gust-radius
date: 2026-10-07
class: position-or-era-artefact
caught_by: thread
became: skill
became_ref: .claude/skills/hf-harness/SKILL.md (step 4, question 2; ERA5 section, ocean points only); fix in 40c3ea5
repeat_of: 
workflow_change: 
---

**What happened.** Pipeline B's gust features were drawn from a 500 km disc with no land mask. On a sample of candidates, 51% had more than 10% land in the disc and 20% had terrain above 2000 m within it (`20b8656`). The contamination varies with position, and position is what the teleconnection work measures.

**How it was caught.** The session read a skill written for the same task that specified ocean gusts, then measured its own exposure before deciding.

**What it was turned into.** A masked re-extraction of all 47 seasons (`40c3ea5`, reported finished in `c9dc994`), and question 2 in `hf-harness`.
