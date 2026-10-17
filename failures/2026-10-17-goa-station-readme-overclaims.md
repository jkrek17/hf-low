---
id: 2026-10-17-goa-station-readme-overclaims
date: 2026-10-17
class: claim-exceeds-evidence
caught_by: reviewer
became: test
became_ref: tests/invariants/run.py (overclaim-words)
repeat_of: 2026-10-17-goa-hrrr-readme-overclaims
workflow_change: a lint rejects "really do see", "are real", "not just a gust/bias" and "under-counts" in research READMEs unless the sentence says it was not shown; and the registered verdict word is used in the headline, with any unregistered comparison labelled post hoc
---

**What happened.** The first draft of the coastal-station README headed with "Yes, partly" where the registered H2 verdict was "supported" (the "partly" came from an unregistered comparison with the "both" class), said the stations "really do see" strong gusts when 8 to 9 in 10 HRRR-only times had none at 45 kt, said HRRR "does not read higher" than the stations when the ratios were taken on station-selected times (selecting on the model instead gives 0.86 to 1.49), claimed "ERA5 under-counts" and "not just a gust-diagnostic bias" with no test, applied the shoreline finding to the Panhandle (peak cell about 26 nautical miles out), and gave the sample as running to April 2026 when the station files stop on 1 January 2026.

**How it was caught.** hf-reviewer, before the pull request, after a fresh agent had matched every point estimate.

**What it was turned into.** A lint (`overclaim-words`) on research READMEs, beside the `storm-force-label` lint of the previous entry, plus a rule in this entry's workflow change: headline in the registered word, post hoc comparisons labelled. A sample shortfall against the pre-registered window (`incomplete-sample` class, second occurrence) was also found; the mechanical check for it (comparing the pair table's last date with the window stated in the pre-registration) is not built and is raised in the pull request.
