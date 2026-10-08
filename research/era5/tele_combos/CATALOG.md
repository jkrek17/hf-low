# Which teleconnection combinations could change where HF lows form, track and peak? (Stage 1 catalog)

Written 2026-10-08, before any model in this directory was fitted. ERA5 **proxy**, **pipeline A**
(`research/era5/hf_history`: 800 km ocean gust index, HF-equivalent at 71.7 kt). Not observations.

## The question, in Jason's words (2026-10-08)

"Do the states of these provide constructive or destructive feedback to the atmosphere? Does a +PNA and
El Niño cause more HF lows? Does a farther-west warm anomaly of El Niño lead to more genesis over the
Kuroshio?"

Read as statistics: for two teleconnection states A and B, is the HF-low response to *both together*
different from what the two give separately, in **where** the lows peak or how many there are in a named
region? *Constructive* means the joint effect is larger than the sum of the separate effects on the log scale
(counts) or on the position scale (degrees); *destructive* means smaller or opposite. That is an
**interaction coefficient** (A x B), with the main effects held in the model. Interaction is the whole claim.

Why the earlier test could not see it: `/mnt/project-files/teleconnection-test/REPORT.md` fitted weekly
**basin-wide** counts. A change in *where* lows peak leaves a basin count unchanged. Its answer, superposition
(no interaction in basin totals), holds for totals only. It says nothing about location.

## What is and is not available (checked 2026-10-08)

| Item | Status |
|---|---|
| Daily NAO, PNA (CPC), 1979-2026 | `/mnt/project-files/teleconnection-test/cpc_indices/`, complete but for 14 missing days |
| ONI (monthly), 1978-2026 | `cpc_indices/oni_1978_2002.txt` plus `docs/data/teleconnections.json` |
| MJO as ten longitude indices, pentad, 1978-2026 | `cpc_indices/proj_norm_order.ascii`. A longitude series, not RMM phases. Sign: negative = enhanced convection. The monthly mean of the 100E-140E series correlates 0.4-0.5 with ONI (`teleconnections.json`, `conventionEvidence`), so the raw series carry ENSO. Stage 2 high-passes them. |
| AO (daily) | `teleconnections.json`, from 2001 only |
| QBO, polar vortex | Not on hand. Both are derivable from ERA5 zonal-mean wind (WeatherBench2 and ARCO are reachable; a 10 hPa, 60N daily series is under 10 GB) |
| HF-equivalent event with **peak position** (`peak_lat`, `peak_lon`) | `hf_history/results/all_tracks.csv.gz`: 75,087 tracks, 1979-2025. Position is where the 800 km gust index peaks. |
| **Genesis position** for all tracks | **Not committed** for 1979-2003. Fix-level files exist for 2004+ only (`intensity/results/fixes_2004.csv.gz`, 00/12 UTC). Rebuilding earlier fixes means re-running the extraction (~370 GB, gated). So stage 2 uses **peak position**, not genesis position. Where "genesis" appears below, read "where the low is strongest". |

Marginals (Oct-Apr, pipeline A, from `all_tracks.csv.gz`; counts only, no association with any index looked at):
2004-05..2025-26: 10,149 Pacific lows, 826 HF-equivalent (about 38 a season); 9,779 Atlantic lows, 1,007
HF-equivalent (about 46 a season). 1979-80..2025-26: 21,698 / 1,738 Pacific, 20,747 / 2,074 Atlantic.
HF-equivalent peaks sit at Pacific mean 44.5N (sd 6.6), 176E (sd 23.5); Atlantic mean 54.3N (sd 9.1), 318.5E (sd 37.1).

## How to read the power column

Effective n = seasons x independent values of the **predictor** per season, not the number of cyclones.
Values per season (Oct-Apr, 210 days), from the additive test and from lag-1 autocorrelation of each series:
ONI about 1; NAO about 5.3; PNA about 5.7; high-passed MJO about 4 (to be re-estimated in `power.py`);
polar vortex about 2-3; QBO about 0.4 (a 28-month cycle: roughly 10 independent values in 22 seasons).
An interaction needs both factors to vary independently inside the same seasons, so its effective n is
set by the **slower** factor times a joint-variation penalty: with ONI it is close to the number of
seasons, 22 (gust-based, 2004-05 on) or 47 (pressure-based, 1979 on). Cyclone-level n is large (800 to 1,000
HF-equivalent peaks per basin) but those peaks are not independent draws of the index.

`power.py` (committed before any fit) turns this into a detectable effect: season-block permutation of the
index series against the real outcome series, giving the null spread of the interaction coefficient. The
detectable effect at 80% power is 2.8 x that spread. It reads outcome *marginals* only; any association of
outcome with index is destroyed by the permutation.

## The combinations

Regions: Pacific box set P-W (27-45N, 135-165E, Kuroshio-Oyashio), P-C (35-55N, 165-200E), P-E (40-60N,
200-240E, Gulf of Alaska). Atlantic box set A-W (27-50N, 262-305E, US coast and Gulf Stream), A-C (40-62N,
305-335E), A-E (50-67N, 335-10E, Norwegian Sea and Nordic seas).

| # | Pair | Basin | Mechanism | Expected signature | Index and lag | Power (effective n) | Covered by |
|---|---|---|---|---|---|---|---|
| 1 | **ONI x PNA** | Pacific | El Niño favours +PNA (deeper Aleutian low, jet pushed south and east). If both act on one jet, +PNA with El Niño is stronger than either; if ENSO acts only through PNA the interaction is zero and the effect is a pathway. | More HF peaks in P-E and P-C, fewer in P-W; peaks displaced south and east. Count total up. | PNA days -10..-4, ONI of day -7 month | 22 seasons x ~1; low | **Thread "El Niño and PNA together"** (branch `claude/enso-pna-interaction-cy2cog`). Family member, not rerun. |
| 2 | **ENSO flavour x Kuroshio genesis** | Pacific | A warm anomaly farther west (central-Pacific El Niño) moves the tropical heating and the jet-exit region west, and changes Kuroshio-front baroclinicity. Jason: more genesis over the Kuroshio. | More HF peaks in P-W. | Flavour index (Modoki-type) against ONI | ~22; low | **Thread "El Niño flavor and Kuroshio genesis"**. Family member, not rerun. |
| 3 | **ONI x MJO** | Pacific | The MJO's convective heating over the Maritime Continent and west Pacific extends (suppressed convection west, enhanced near the dateline) or retracts the East Asian jet about a week later (Moore et al. 2010; Henderson et al. 2016; from memory, not rechecked). El Niño already extends the jet. When both push the same way the jet response adds; against each other it cancels. | Peak longitude of Pacific HF lows shifts east when both extend the jet and west when both retract it; P-C and P-E counts follow. Mean position shift is the signature, not total count. | MJO 120E and 140E series, high-passed, pentad centred days -10..-4; ONI of day -7 month | ~22 x ~4; **moderate** | Not covered. Basin-count version tested: archive q = 0.045, ERA5 1979-2025 RR 0.991 (no effect). **Location version is new.** |
| 4 | **PNA x MJO** | Pacific | The MJO in phases that force the PNA teleconnection on days 5-10 (Henderson et al. 2016) adds to or fights the slowly varying PNA already in place. | Same displacement as #1 (south and east) in aligned cases; Gulf of Alaska (P-E) count; peak latitude. | PNA days -10..-4; MJO 160E and 120W series, high-passed, same window | ~22 x 5.7 x 4; **moderate** | Not covered. Basin-count secondary: RR 1.10 in the archive (p 0.03-0.29), 1.02 in ERA5 1979-2000. Location version is new. |
| 5 | **NAO x PNA** | Atlantic | Pacific-to-Atlantic link: a +PNA ridge over western North America reinforces a trough over eastern North America, changing the Atlantic jet entrance and where cyclones deepen (synoptic-eddy feedbacks; Hoskins-type wave-train arguments; from memory). | A-W (US coast) share changes, peak longitude shifts west with +PNA under -NAO and east under +NAO. | NAO and PNA days -10..-4 | ~22 x 5.3 x 5.7; **good** (both fast) | Not covered. Basin-count version tested: no interaction. **Location version is new.** |
| 6 | **ONI x NAO** | Atlantic | ENSO modifies the NAO's effect on the Atlantic storm track (El Niño tends to a southward, eastward track; Pozo-Vazquez et al. 2001-type results; from memory). The NAO shift in storm track may be larger or smaller by ENSO state. | Peak latitude and longitude of Atlantic HF lows; A-C vs A-E. | NAO days -10..-4; ONI of day -7 month | ~22 x ~1; low | Not covered. Basin-count version tested: null. **Location version is new.** |
| 7 | **MJO x NAO (chain)** | Atlantic | The MJO precedes an NAO phase change by 7-12 days (Cassou 2008). The question is whether an MJO-primed NAO has a different effect from a free-running one. | NAO-type shift in track; mostly a pathway, an interaction only if the primed NAO differs. | MJO 70E-100E high-passed, days -20..-11; NAO days -10..-4 | ~22 x 4; moderate | Not covered. Basin count: null. **Not picked**: mostly a pathway question with a weak prior for an interaction, and it overlaps #5's NAO variation. |
| 8 | **NAO x stratospheric polar vortex** | Atlantic | A strong vortex couples to a +NAO/poleward jet over the following weeks; a weak one to -NAO (Baldwin and Dunkerton 2001; Kidston et al. 2015). Interaction: NAO already +/- *and* the vortex agrees, so the jet is reinforced; if it disagrees, the troposphere-driven NAO is damped. | Peak latitude of Atlantic HF lows (poleward under both strong); A-E count (Nordic seas). | NAO days -10..-4; 10 hPa 60N zonal-mean wind anomaly, mean of days -30..-11 | ~22 x 5.3 x ~2.5; **low to moderate** | Not covered. New index to derive; this is the only combination that adds a distinct region the other tests cannot reach. |
| 9 | **ONI x QBO** | both | Holton-Tan: westerly QBO with La Niña/neutral strengthens the vortex; easterly QBO with El Niño weakens it, acting on the NAO downstream. | Only through #8. | ONI, QBO 30 hPa or 50 hPa | ~10 independent QBO values in 22 seasons; **too low** | Not covered; **not picked**. A well-powered result is impossible. |
| 10 | **ONI x PNA x MJO (three-way)** | Pacific | Jet extension needing all three. | | | Far too low | **Not picked**: a three-way interaction on 22 seasons cannot be powered. |

Two further lists were considered and rejected as primary tests: AO (only 2001 on, so 25 seasons, redundant with NAO and
the polar vortex) and the Atlantic NAO x MJO (#7).

## Stage 2 picks

Picked on mechanism, distinctness from the two sibling threads, and the power column. Not on any result.

| Test | Pair | Why picked |
|---|---|---|
| T1 | ONI x MJO, Pacific | Jason's jet-extension idea; ENSO set the base state and the MJO moves it on a faster clock; moderate power |
| T2 | PNA x MJO, Pacific | the direct Pacific wave-train version; fast indices; moderate power |
| T3 | NAO x PNA, Atlantic | Pacific-to-Atlantic linkage, both fast; the best-powered pair |
| T4 | ONI x NAO, Atlantic | the second ENSO combination; low power, so a null here will be called inconclusive unless the interval is tight |
| T5 | NAO x polar vortex, Atlantic | the one distinct-region modulator; needs a derived 10 hPa index; contingent on that derivation passing its own check (rule in the plan) |

Tests #1 and #2 of this catalog (ONI x PNA; ENSO flavour x Kuroshio) are the two sibling threads. They are
members of the family. Their primary p-values enter a joint false-discovery correction at the end if their
results land in time. Mine are corrected on their own and with theirs.

Not picked: #7 (pathway, weak interaction prior), #9 (cannot be powered), #10 (cannot be powered).

## References (recalled from memory, **not rechecked**)

Moore, Martius and Spengler (2010, Mon. Wea. Rev.); Henderson, Maloney and Barnes (2016, J. Climate);
Cassou (2008, Nature); Baldwin and Dunkerton (2001, Science); Kidston et al. (2015, Nature Geoscience);
Holton and Tan (1980); Hoskins and Hodges (2002); Eichler and Higgins (2006). Use them to motivate the
mechanism; check each before quoting in an article.
