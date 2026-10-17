---
id: 2026-10-17-gulf-stream-region-prose-and-s5
date: 2026-10-17
class: claim-exceeds-evidence
caught_by: reviewer
became: test
became_ref: research/era5/gulf_stream_region/test_convert_en.py (conversion size check); hf-reviewer list applied before the pull request opened
repeat_of: 2026-10-08-gust-headline-walked-back
workflow_change: the verdict-wording rule in PR 264 (open for Jason) is the upgrade; this entry adds that every derived figure in README prose (placebo-corrected values, knots, shares) must be printed by the analysis script into a result file, and that a secondary correction (here the satellite-to-true-wind conversion) gets a size test before it is quoted
---

**What happened.** The first README for the region-wide Gulf Stream round (verdict words were right this time) said GFS "carries most of" the Stream effect and that the bump is "biggest away from the Stream front" in bold although the plan reads that as can't tell; gave a mechanism (cold air outbreaks) that no result file covers; said the shuffle check was run for A1 only (B1 had one); said "no test here is a null" (two were); said the front is undefined east of 65W (30% and 43% of boxes are east of 65W); quoted a 7% figure that is 3.5% in absolute terms; quoted a placebo-corrected value that only existed in prose. Separately, the equivalent-neutral conversion (S5) used the bulk Richardson number as the stability parameter and so moved A1 by 0.8 points instead of about 4; with the correct size, the GFS gap in true 10 m wind goes from +2.2% to -1.9%, which reverses the forecaster-facing sentence.

**How it was caught.** hf-reviewer, after a verifier matched the point estimates. A first verifier's placebo disagreed because its own shuffle wrote values to the wrong rows; the second, told how to write them back, matched.

**What it was turned into.** A size test for the conversion (`test_convert_en.py`: a 3-point correction at 12 m/s for 7.5 versus 1.5 C), the placebo-corrected value and band shares now printed by `analyse_region.py` into the result file, and the reviewer's list applied before the pull request. Another `claim-exceeds-evidence` entry (the Gulf of Alaska entries of 2026-10-16 and 2026-10-17 are on the research branch, not yet on main); the process change is PR 264 (open for Jason).
