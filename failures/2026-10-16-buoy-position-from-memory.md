---
id: 2026-10-16-buoy-position-from-memory
date: 2026-10-16
class: position-or-era-artefact
caught_by: reviewer
became: skill
became_ref: .claude/skills/hf-harness/SKILL.md (ERA5 section, point comparisons read and assert positions and print the cell's land fraction); research/era5/goa_jets/buoy.py check_positions
repeat_of: 2026-10-07-land-in-gust-radius
workflow_change: any comparison at a named point reads its position from the source, asserts it to 0.05 degree and prints the model cell's land fraction before a ratio is computed; a "the model cannot see X" finding is checked at the true position first
---

**What happened.** The Gulf of Alaska jets study (research/era5/goa_jets, branch commit `b96b45a`) compared ERA5 with seven NDBC buoys at positions typed from memory. Buoy 46077 (Shelikof Strait) was put at 156.94W; NDBC has it at 154.211W, about 160 km east. The wrong position fell in an ERA5 land cell (land fraction 0.92), so the README said "ERA5 cannot see the strait or the inlet" and that real gap winds are "underestimated by an unknown, large amount". At the true position the cell is 14% land, ERA5 wind is 0.81 to 0.88 of the buoy's and the gust about equal.

**How it was caught.** The hf-reviewer claim review before the pull request was opened. It fetched the NDBC station page and checked the land fraction at both positions. The verifier had only been asked about the numbers, and they did match the (wrong) pairs.

**What it was turned into.** `buoy.py check_positions` (reads the NDBC station pages, asserts to 0.05 degree, prints land fraction), a bullet in the hf-harness ERA5 section, corrected results and README in `c8a809d`. Second position-or-era-artefact entry after 2026-10-07-land-in-gust-radius (a location-dependent error read as a property of the method).
