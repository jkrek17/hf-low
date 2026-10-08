# High-latitude Atlantic HF lows: what was found

Plan and hypotheses: `PLAN.md` (committed before the analysis ran, `b938ecd`, `9c916ef`).
Numbers: `results.txt` (archive and ERA5 pipeline A), `pcs-result.txt` (station).
Everything ERA5 here is **pipeline A, a proxy**, seasons 2004-05 to 2025-26, Atlantic.
Intervals are season block bootstraps (22 seasons).

## Answer in one paragraph

Two things stack up. Most of the Atlantic's high-latitude share is storms
whose own wind field reaches HF strength between Iceland and Greenland. About a
third of the high-latitude HF-strength fixes in pipeline A instead take their
gust from a point 400 km or more from the centre and within 300 km of
Greenland, mostly in northerly or north-easterly flow along the east coast
(barrier-wind type), not the westerly forward tip jet. Those terrain-type
fixes are poorly matched by the archive, so pipeline A is less trustworthy
north of 60N, not untrustworthy: the extra false alarms are 18 points of FAR
(26% to 44%), and recall is the same.

## Numbers

| Quantity | Value | 95% CI |
|---|---|---|
| Archive Atlantic HF lows with an HF fix >= 60N | 377 of 1,011 = 37.3% | 32.7-42.3 |
| Archive Pacific | 12 of 856 = 1.4% | 0.7-2.2 |
| Pipeline A Atlantic events peaking >= 60N | 34.0% (2004-05 on); 37.5% all seasons | |
| Gust maximum > 400 km from the centre: fixes >= 60N vs < 60N | 44.1% vs 19.4% | |
| Maximum within 300 km of Greenland: >= 60N vs < 60N | 62.8% vs 7.7% | |
| Terrain type (both of the above): difference >= 60N minus < 60N | +0.268 (31.4% vs 4.7%) | 0.227-0.301 |
| Index needs the Greenland 300 km band (events losing HF) | 20.4% of events; 43.7% of those peaking >= 60N | |
| Peak >= 60N share with the Greenland 300 km band removed | 24.0% (from 34.0%) | |
| Index needs the Greenland 100 km band | 5.5% of events | 4.2-6.9 |
| Index needs sea-ice points (> 0.15) | 2.3% of events (5.6% of those peaking >= 60N) | 1.6-3.0 |
| FAR >= 60N vs < 60N (peak position) | 0.443 vs 0.258, difference +0.185 | 0.122-0.244 |
| POD >= 60N vs < 60N (archive position) | 0.759 vs 0.744, difference +0.014 | -0.040 to +0.071 |
| Archive match rate, events needing the 300 km band vs not | 47.1% vs 73.3% | |
| Archive pressure at HF, Denmark Strait/Iceland minus south of 60N | -6.5 hPa | -9.5 to -3.0 |
| Archive pressure at HF, Cape Farewell/Irminger minus south | -3.5 hPa | -7.0 to 0.0 |
| Prins Christian Sund mean wind, median, when A's maximum is at Cape Farewell vs other times | 36.0 vs 13-17 kt | not bootstrapped |

## Reading, by pre-registered hypothesis

- **H1, terrain jets: supported as a description, not as the whole answer.**
  Terrain-type maxima are six times as common at fixes north of 60N (31% vs
  5%). But 69% of high-latitude HF-strength fixes are not terrain type, and
  with the Greenland 300 km band removed the Atlantic still has 24% of events
  peaking north of 60N against 0.4% in the Pacific. That removal also takes
  out genuine storm winds near the coast, so it overstates how much is
  terrain and 24% is a floor on the share that the storm's own field carries.
  Where terrain-type maxima occur, the wind is from the north or north-east in
  79% of fixes and from the west in 20%, and only 19% lie in the Cape Farewell
  box, so this is mostly barrier flow along the east coast and the Denmark
  Strait, not the classic forward tip jet. The fixes are 315 and they cluster
  in storms, so read these percentages loosely.
- **H2, weaker central pressure near Greenland: not supported.** The sign is
  the opposite of the pre-registered one: HF fixes at 60N and beyond come with
  lower pressure than those further south (-3.5 to -6.5 hPa; the Cape Farewell
  interval touches zero). The Iceland Low lowers pressure there for any
  storm, so this does not rule out orographic enhancement; it only fails to
  show it.
- **H3, proxy artefact: partly supported.** The pre-registered rule called for a
  recommendation if FAR north of 60N exceeded FAR south by more than 0.10 with
  an interval excluding zero. It did (+0.185, 0.122-0.244). The other
  trigger, more than 10% of events depending on points within 100 km of
  Greenland or on sea ice, was not met (5.5% and 2.3%). Recall does not
  differ (+0.014, interval spans zero), which is a null with moderate power:
  the interval allows a difference of about 7 points either way.
- **Station check.** At Prins Christian Sund the mean wind is 36 kt (median)
  at the 202 times when A's maximum is in the Cape Farewell box, against 13-17
  kt at the other times; 54% reach 34 kt against 11-12%. The strong wind there
  is real. It is not hurricane force (no report at 64 kt in any group), and the
  ERA5 index does not track the station within those times (Spearman -0.18,
  n = 163). A 19 m coastal station in a fjord mouth cannot confirm 64 kt.

## What "false alarm" means here

An ERA5 event with no archive match is a false alarm against the archive. At
high latitude that may also be OPC not warning (sparse observations, the
archive's area ends at 66.6N) so the 0.443 mixes proxy error with archive
omission. This study cannot separate them. A scatterometer record would.

## Does pipeline A need a terrain mask?

My recommendation, not a change: do not mask. A Greenland 300 km mask would
drop 20% of events, including real storm wind. Report pipeline A's Atlantic
events north of 60N as a separate stratum with its FAR of about 0.44, and treat
the 225 events whose index depends on the Greenland 300 km band as lower
confidence (47% archive match against 73%). A hard mask at 100 km is cheap
(5.5% of events) and could be tried as a sensitivity check.

## Deviations and limits

- Error caught by the verifier: the first run took the sea-ice maximum over
  non-blank fixes only, which counted 28 events (2.5%, CI 1.9-3.2) as losing HF
  without ice points. Blank means ice was not fetched (fix south of 55N) and the
  index there is unchanged; with that fixed it is 25 events (2.3%, CI 1.6-3.0).
  The corrected 25 was recomputed by the verifier; its interval and the 5.6%
  were not. Nothing else changed. Both versions are in git history.
- Consistent with a separate thread: the wind-structure thread (PR #27) found
  14% of Atlantic HF fixes have their gust maximum within 100 km of Greenland
  or Iceland. Here 13.4% of the 3,303 fixes do (12.8% Greenland alone). Different
  fix sets and code, so this is agreement, not a check.
- The agenda's 37.5% used 1,005 events; the current payload gives 377 of
  1,011 (37.3%). Logged in `PLAN.md`.
- Pipeline A's domain stops at 67N and the archive's latest Atlantic fix is
  66.6N, so nothing north of that is tested.
- The Pacific side was not studied directly. Why the Pacific has almost no
  high-latitude HF (its storm track ends near the Gulf of Alaska) is outside
  this data.
- The agenda's question about a trend in the high-latitude share (sea-ice
  retreat) was not run; it is not in `PLAN.md` and would have low power.
- One ERA5 extraction, 3,094 times, about 27 GB.
