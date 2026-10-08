# Is the lagged hemispheric pattern's skill upstream or local? Pre-registration

Agenda item RA-2 (`/mnt/project-files/research-agenda/AGENDA.md`). Written and committed **before any region model was
fitted or scored**. At that point nothing had been run in this thread: no field extraction, no outcome joined to a region.
Background that shaped the design (all from the PR 41 README, `research/era5/hemispheric/`): the pattern is a ridge-penalised
Poisson model on 10 EOF scores each of Z500, U250 and SST (5.6 degree, 12 UTC daily, 7-day mean over days -7..-1 before the
week); archive weekly counts, Oct-Apr, 22 seasons 2004-05..2025-26; baseline B1 = month + log(1 + previous week's count).
ERA5 fields are a **proxy** for the atmosphere; the outcome is the **archive**; pipeline A (`research/era5/hf_history`) enters
only in the secondary check R1 (depth counts 1979-2000) and in S-imprint, and says so.

Question (Jason, via the agenda): *does the pattern's week-ahead skill come from upstream fields (a remote wave train) or from
local ones (the jet and trough over the basin itself)?* Answer wanted in plain words per basin: how much of the skill is upstream,
how much local, and how sure.

## The held-out reservoir, stated first

The PR 41 held-out seasons (2015-16..2025-26) were scored four times as pre-registered (primary, swap, lag 2, concurrent) plus
the post hoc scripts, and the 22-season attribution and leave-one-season-out descriptives used all 22 seasons
(`hemispheric/results/heldout_looks.log`). **No fresh archive-count split exists**: every archive season 2004-05..2025-26 has
been in a PR 41 look, and the 2026-27 season has just begun (about 1 week of 30). The design therefore does not rely on a fresh
split for its primary answer:

- The primary answer compares **region models against the full pattern and against each other on the same weeks, under the
  same nested leave-one-season-out scheme.** Regions, variables and the decision rule are fixed here, before any score. The
  comparison is a difference between models that are all out-of-sample on each season, and none is chosen from outcomes, so
  re-using the seasons does not select anything. It does mean the absolute skill numbers are not new evidence for the pattern.
- The one source of data that no region model has met is **pipeline A depth counts in 1979-80..2000-01** (proxy outcome, a
  different definition, within-era variation only under decision 1). PR 41's S5 used it once, with the frozen full pattern. R1
  below uses it for the first time for every region model and for the full pattern refitted on 22 seasons (a second look for the
  full pattern, a first look for the others). This is the replication, and it is counted in `looks.log`.

Every score computed in this thread is appended to `results/looks.log` with the model set and the seasons.

## Regions (fixed here)

Fields on the PR 41 grid: 12 latitude rows (25.3-87.2N) x 64 longitudes (5.625 degree cell centres, 0..354.4E). A cell belongs
to a sector by its centre longitude in (-180, 180]; sector edges are half-open (a, b]. Latitude is not cut: all 12 rows. Three
sectors of 120 degrees tile the hemisphere. **Sectors per target basin** are chosen from the basin's pipeline A domain (Atlantic
98W-10E, Pacific 135E-120W), widened to the 120-degree width, not from the pattern maps:

| basin | LOC (local: the basin's own sector) | UP (upstream: the 120 degrees to its west) | DN (downstream: the remaining 120 degrees to its east) |
|---|---|---|---|
| Atlantic | (-100, 20]: 100W to 20E | (140, 180] and (-180, -100]: 140E to 100W | (20, 140]: 20E to 140E |
| Pacific | (120, 180] and (-180, -120]: 120E to 120W | (0, 120]: 0 to 120E | (-120, 0]: 120W to 0 |

Reading the Atlantic row against the PR 41 maps: the Alaska / north-east Pacific ridge (about 135W) is in UP; the North
American trough core (about 90W) sits **inside LOC**, and the jet from 90W to 20W lies wholly in LOC. That placement is the
literal reading of the basin domain and is **not** what the agenda text calls "upstream". So the trough-upstream version is
pre-registered as sensitivity variant V1 below, with the same status as the primary. Pacific: the PR 41 trough at the dateline
and the extended jet are in LOC; the central North American ridge (about 90W) is in DN, which is downstream in the usual sense.
A one-week-lag pattern has no clean upstream/downstream separation (a wave packet at 30 m/s travels about 1,600 km/day, so
DN fields at day -7 can have been shaped by the basin a few days earlier): "DN-only keeps skill" is reported as a diagnostic of
how redundant the hemisphere is, not as a mechanism.

## Models (all fitted per basin on the same outcome; Z500 + U250 + SST, K = 10 EOFs per field, as PR 41 unless stated)

Each model restricts the grid cells the EOFs see; EOFs, standardisation and the ridge penalty are otherwise identical to PR 41
(`hemlib.py`, lambda grid {1,...,3000, none}, chosen by inner leave-one-season-out). "No PCs" may win.

- **Primary (7 per basin):** FULL (the PR 41 pattern, a reproduction check); keep-only UP, LOC, DN; drop-one noUP (= LOC+DN), noLOC
  (= UP+DN), noDN (= UP+LOC).
- **V1, trough upstream (2, Atlantic only):** keep-only UP1 = (140, 180] and (-180, -70] and keep-only LOC1 = (-70, 20]. The
  North American trough core goes to UP. The Pacific has no analogous boundary question, so V1 is not run for it.
- **V2, variables (4 per basin):** keep-only UP with Z500 only; UP with U250 only; LOC with Z500 only; LOC with U250 only
  (10 PCs, no SST). Reading: "ridge/trough heights upstream" against "the local jet".
- **V3, sliding window (16 per basin):** 120-degree windows centred every 22.5 degrees, keep-only, all three fields. Descriptive
  curve of skill against window centre; not a test of anything beyond the 16 q values it shares a family with.
- **S-imprint (secondary):** FULL, UP and LOC refitted with the baseline B2 = B1 + log(1 + all pipeline A cyclones in the basin in
  the previous 7 days) (PH1 of PR 41). Reason: storms of the previous week imprint on fields over the basin, so LOC can carry
  skill that is storm clustering; UP cannot be imprinted by the basin's storms in the same way. Pipeline A counts are a proxy
  and enter only as a covariate.

## Scoring (nested leave-one-season-out, all 22 seasons)

For each held-out season s: EOFs and lambda from the other 21 seasons (inner leave-one-season-out), models refitted on those 21,
scored on season s (the scheme of `hemispheric/attribute.py`, which gave FULL +5.74% Atlantic and +3.40% Pacific; reproducing
these two numbers exactly is a **pre-condition**: if FULL does not reproduce, nothing else is reported until it does).

Skill: pooled Poisson deviance, SS(M|B1) = 1 - sum D_M / sum D_B1 over the 660 weeks per basin.
Shares: retained share f_M = SS(M|B1) / SS(FULL|B1). **Shapley share** of each sector: with the 8 subsets of {UP, LOC, DN}
(empty = B1, SS 0; singles = keep-only; pairs = the three drop-one models; all three = FULL), the average marginal gain of a
sector over the 6 orderings divided by SS(FULL). The three Shapley shares sum to 1 by construction; the primary headline is
`Shapley(UP)` against `Shapley(LOC)`, and the keep-only f values give the agenda's test.

Uncertainty: season-block bootstrap, 10,000 draws over the 22 seasons, shared across models (paired). Intervals are 90%
for shares (reported) and 95% for differences. Permutation p for SS > 0: 10,000 season-block permutations of the penalised-PC
contribution (the PR 41 `perm_p` rule), one-sided. Seed 20261008 (+ offset per stage).

## Prediction and decision rules (fixed now)

Prediction (from the agenda): **Atlantic: UP-only retains at least half of FULL's skill and LOC-only less than half.** If the
reverse holds, remote control is not supported. The same two numbers are reported for the Pacific with no directional prediction
(the agenda gave none); the Pacific pattern looks like a strong +PNA, which is a local trough with a downstream ridge, and a
local-dominant result is not surprising there.

Per basin, only if FULL has SS > 0 at q < 0.05 (Family 1):

- **REMOTE**: f_UP >= 0.5 and f_LOC < 0.5. **LOCAL**: f_LOC >= 0.5 and f_UP < 0.5. **SHARED**: both >= 0.5 (either region carries
  it; the wave train is redundant across regions and the two cannot be separated). **DISTRIBUTED**: both < 0.5 (skill needs the
  combination).
- A category is **resolved** only if the 90% bootstrap interval of both f_UP and f_LOC lies wholly on its side of 0.5.
  Otherwise the answer is "unresolved, leaning X" and says what interval overlaps. If f_FULL's own interval includes 0 for a
  basin the ratios are not interpretable and the answer is "can't tell".
- The Shapley shares must agree in direction (UP vs LOC larger) with f for the answer to be called robust to the order of
  ablation. If drop-one and keep-only disagree (for example UP-only keeps most while noUP loses little), the report says the
  skill is redundant across regions, not that one region is the cause.
- If DN-only also retains at least half of FULL's skill, regions are not separable in this record, and the answer is "the pattern
  is hemispherically redundant" whatever UP and LOC show.
- A pairwise difference of keep-only skills (UP - LOC, UP - DN, LOC - DN) is resolved only if its 95% interval excludes 0.
  The report states the smallest difference resolvable at 80% power from the bootstrap standard error (2.8 x SE).

R1 (replication on data the region models have not met): pipeline A depth counts, 1979-80..2000-01, proxy outcome, within-era
variation only (decision 1; no level or trend). Each model is refitted on all 22 archive seasons, its penalised-PC linear predictor is
standardised on that era, and the S5 test is applied (month + previous week + gamma x index; one-sided season-block permutation,
5,000). The same classification rule is applied to the ratios of deviance gains. R1 **agreeing** with the primary (same
category, or same larger of UP and LOC) is called a replication; disagreement is reported as a failed replication and the primary
answer is then "unresolved".

## Multiplicity

Per basin and tests (all one-sided unless stated):

- **Family 1** (7): SS > 0 for FULL, UP, LOC, DN, noUP, noLOC, noDN.
- **Family 2** (6): pairwise keep-only differences UP-LOC, UP-DN, LOC-DN (two-sided, bootstrap); drop-one losses
  SS(FULL) - SS(noUP), - SS(noLOC), - SS(noDN) > 0 (bootstrap one-sided).
- **Family 3, secondary**: V1 (2 SS > 0), V2 (4 SS > 0), V3 (16 SS > 0), S-imprint (3 SS > 0), R1 (7 gamma > 0).
BH q is reported within each family and across every test in the thread, per basin and over both basins. The report states how
many tests pass at q < 0.05 out of the total.

## Power

Resolution is stated, not assumed. (a) The bootstrap SE of SS(FULL) and of each difference; the smallest keep-only difference
resolvable at 80% power is 2.8 x its SE. (b) Agenda note: differences under about 1.5 skill points are treated as unresolved
unless the SE says otherwise. (c) A well-powered null on "UP-only retains half" needs the interval for f_UP to be narrower than
+/- 0.25; if the interval is wider the answer is "can't tell", not "no".
Pacific FULL skill is small (3.4%), so ratios there have wide intervals by construction.

## What this does not do

- It is not a causal test. Cutting fields by sector tells where the predictive information sits, not what drives storms.
- Sector edges are a choice; V1 and V3 show how much the answer depends on them.
- Wave-train collinearity means a null for one block does not make the other the cause.
- 22 seasons, no fresh archive split; absolute skill is not new evidence. ERA5 is a proxy; R1 uses pipeline A depth counts, a proxy
  outcome with a different definition and within-era use only.
- 5.6 degree daily fields do not resolve jets in detail.

## Deviations (post hoc)

None yet.
