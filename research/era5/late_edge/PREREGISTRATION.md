# Is the high-latitude Atlantic "late HF" effect an artefact of the domain edge? Pre-registration (RA-26 follow-up)

Written and committed **before** the lateness contrast (RA-7's T2) was tabulated on any subset defined by edge contact in this thread.

What had been looked at: (1) PR 86's headline (T2: late share 52.7% vs 26.9%, OR 3.10, 2.85 without terrain-flagged fixes; joint OR 1.82) and the RA-26 result (PR 94); (2) **counts only**, no outcome:
Atlantic events 1,039; group A60 (onset gust maximum at or north of 60N) 283; of those 129 touch the north or east edge as defined below and 154 do not; south of 60N 756, 214 touching and 542 not.
(3) How `lifecycle.py` defines the minimum: `h_on_minp` uses the lowest MSLP over **all** fixes of the track, in or out of domain (the tracker runs 20-75N). So the lowest-pressure fix is not truncated at the domain edge; the agenda's mechanism has to act through onset instead: an HF fix counts only in domain, so a storm that is HF before it enters the domain, or whose in-domain part starts after its minimum, gets a late "onset". This is the artefact tested here.

All results are ERA5 **proxy**, **pipeline A**, seasons 2004-05 to 2025-26, non-TC events (as RA-7). Pipeline A's Atlantic domain is 30-67N and from 98W east to 10E (`hf_history/track.py`). No pull.

## Question and the answer I will give

Does PR 86's odds ratio of lateness north of 60N survive when tracks that touch the domain edge are removed? Answer in plain words: it stays (and by how much), it shrinks (and by how much), or it vanishes (edge artefact).

## Definitions (fixed now)

- Data and LATE: exactly RA-7's (`late_wind.load()`); LATE = `h_on_minp >= 0`; group indicator T2 = onset gust maximum at or north of 60N (RA-7). Atlantic only.
- **TOUCH (primary)**: any fix on the whole 6-hourly track (`era5_hf_catalog_tracks.csv`, includes out-of-domain fixes) that is outside the Atlantic domain at latitude 50N or higher, i.e. the track left through the north or east edge (or the 75N detector limit) at some time. Southern and western out-of-domain fixes (genesis at low latitude) are not edge contact for this question.
- TOUCH (buffered, secondary): TOUCH, or any in-domain fix north of 64N or between 7E and 10E. Expected to leave about 70 events in A60; reported, labelled underpowered.

## Tests (one family, BH-FDR over all seven)

Model: logistic regression, season as cluster (CR1, 22 clusters), Wald test, as RA-7. Intervals also by season-block bootstrap (2,000 draws, seed 7) as a check.

- **E1 (primary)**: OR of LATE for T2, Atlantic events that do not TOUCH. Full-sample OR (3.10) reproduced first as a check.
- E2: the same OR in TOUCH events.
- E3: difference in log OR between TOUCH and not-TOUCH (interaction T2 x TOUCH).
- E4: E1 adjusted for marginal storm (peak gust < 76.7 kt), translation speed and month group.
- E5: E1 with LATE6 (>= 6 h).
- E6: E1 with terrain-flagged fixes removed (Greenland, within 100 km).
- E7: E1 with the buffered TOUCH.

## Decision rules (thresholds fixed now)

- **Persists** if E1's OR is at least 1.76 (half the full-sample log odds ratio, ln 3.10) and its 95% interval excludes 1.
- **Vanishes (edge artefact)** if E1's interval includes 1, its upper bound is below 3.10 (the full-sample OR is excluded), and power to detect OR 2.0 is at least 80%.
- **Shrinks** if the interval excludes 1 but the OR is below 1.76.
- Otherwise **can't tell**.

## Power

Minimum detectable OR at 80% from the cluster-robust standard error (exp(2.8 SE)), and the power at OR 2.0. No planted-effect simulation.

## Looks

Uses all 22 seasons including 2015-25 once, descriptively; no fit scored on a separate block. By the agenda's count, look #17 at 2015-25 (after RA-26's #16); the agenda thread renumbers at merge. Zero pre-2001 (gust-based onset). Entry goes in `research/era5/looks/late_edge.log`.

## What this cannot show

Even a surviving effect is not explained, only not an edge artefact. TOUCH removes about half the high-latitude events, so the remaining group may differ in other ways (more slow-moving or Greenland-sector storms). No terrain mask; the proxy near Greenland has the highest false-alarm ratio (0.44).

## Deviations (post hoc)

(none yet)
