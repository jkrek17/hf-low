---
id: 2026-10-17-goa-hrrr-readme-overclaims
date: 2026-10-17
class: claim-exceeds-evidence
caught_by: reviewer
became: test
became_ref: tests/invariants/run.py (storm-force-label)
repeat_of: 2026-10-16-goa-followup-verdict-words
workflow_change: a lint on research READMEs rejects "storm-force" next to a zone-maximum, HRRR or ERA5 54 kt figure unless the sentence also says uncalibrated, proxy, or reads as; the reviewer's required-fix list is applied before the pull request opens
---

**What happened.** The first draft of the HRRR README called HRRR's 54 kt zone maxima "storm-force" (they are an uncalibrated model gust), said both models "read about 10 percent high" at the buoys from ratios that span 0.96 to 1.15, put "about 40 times as many cells" where the count is 44 to 48 and untested as a cause, stated the shoreline finding for all three zones when the Panhandle does not show it, left ERA5 unrestricted in the 25 km comparison without saying so, and gave the buoy gust a wrong averaging time. The pre-registration log misdescribed the crop change and the plumbing run and omitted the bootstrap draw count and the kilometre units.

**How it was caught.** hf-reviewer, before the pull request, after a verifier matched the point estimates.

**What it was turned into.** A lint (`storm-force-label`) in `tests/invariants/run.py`, the mechanical form of the label rule that the two earlier entries in this class left as prose. Third entry in the class in two days: the upgrade is the lint plus applying the reviewer's list before opening a pull request.
