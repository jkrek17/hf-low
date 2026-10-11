# Vertical motion (omega) and ascent coupling as predictors of P(HF) (ERA5 proxy, pipeline A, Oct-Apr)

Question (Jason, 2026-10-09): is the search area big enough, could an HF storm sit under one jet's entrance and another's exit, and is there a better measure of vertical lift? Plan: `PREREGISTRATION.md` (committed before any omega field was read; one deviation, below). ERA5 proxy, pipeline A; fit and test 2004-05 to 2021-22, Oct-Apr, 74,009 fixes, 16,874 tracks, 18 seasons, 4,012 hf24 events; leave-one-season-out. WeatherBench2 1.5 deg `vertical_velocity` at 500 and 700 hPa, box +-5000 km. Pull: 1,946 chunks, 21.96 GB (sized 21.9 GB by HEAD request). Transitioning tropical cyclones in.

## Answer in plain words
- **Direct lift helps, but what it measures is the storm itself.** Adding the ERA5 ascent around the low to the PR 12 model raises the Brier skill for "HF within 24 h" from 0.4133 to 0.4221 (+0.0088, 90% interval +0.0055 to +0.0121, 15 of 18 seasons better, q 0.002). For rapid deepening it is +0.0127 (17 of 18 seasons). That passes the pre-registered rule (gain >= 0.005, interval above 0, q < 0.05) and is larger than anything in the jet-streak and trough test (PR 149: +0.0016 to +0.0029 for HF within 24 h).
- **The gain comes from the ascent close to the low, not from lift in the surrounding environment.** The 750-2,500 km ring alone adds nothing (-0.0001, well-powered null), so the information is the ascent within 1,000 km, which includes the storm's own latent-heating ascent. It tracks a storm that is already strong, so it is partly an outcome-side diagnostic, not a precursor. Lags barely change it (lag 0 +0.0086, with 24 h +0.0088, with 12/24/48 h +0.0090), so it is a fix-time reading.
- **Coupled ascent centres add nothing.** Whether a second distinct ascent centre sits within 2,500 km adds -0.0000 (primary rule; interval -0.0006 to +0.0005) and the single flag +0.0003. Both well-powered nulls (gains of 0.005 or more are excluded). The same holds on the other two centre rules and on rapid deepening and hf48. This does not test two jet streaks (that needs 250 hPa winds again, +40.6 GB, not pulled); it tests whether the lift field shows a second centre.
- **On top of the jet and trough features** the ascent still adds +0.0075 for HF within 24 h and +0.0103 for rapid deepening, so it is not just another view of J2/T2 (PR 149).
- **Atlantic vs Pacific:** Pacific +0.0133 (q < 0.001), Atlantic +0.0052 [-0.0001, +0.0099], can't tell (descriptive split).

## Tests
28 non-descriptive tests, 16 pass the rule: all 16 are the W and W+C groups (8 families x 2). The other 12 are 11 coupling tests (C in 8 families, `couple_w` alone, and the two other centre rules), none of which passes and all 11 of which come out as well-powered nulls on the 0.005 threshold, and the ring-only test Wx (also a null). Primary family (BH over 3, hf24, lags 0 and 24 h): W +0.0088 (adds skill), C -0.0000 (null), W+C +0.0086 (adds skill). Full table: `skill.txt`, `gains.csv`; probabilities `hf24_probs.npz`.

## Detection QC (before outcomes were joined; `qc.txt`)
- Q1 passed: an ascent maximum >= 0.2 Pa/s within 2,500 km at 100% of fixes; two or more centres at 31.0% (Atlantic 28.5%, Pacific 33.2%) of fixes under the primary rule (bounds 10 to 70%).
- Q2 passed: ascent within 1,000 km is higher in the jet-streak quadrants RE and LX than in LE and RX by +0.152 Pa/s (Atlantic) and +0.138 (Pacific), q 0.0005, so the lift field and the PR 149 quadrants agree.
- Q3 (12 cases drawn, `qc_cases.png`, `qc_cases.csv`): the detector marks real maxima, but several "centres" are single-grid-point spots at the ring edge or small convective-looking features, and in a few cases the stronger centre is a thin spot far from the storm. This is why the coupling flag may carry little information. Logged, not tuned.

## Descriptive companion (`descriptive.txt`; HF onset against the storm-force-only peak, matched on basin and month, 380 and 296 HF onsets on extracted fixes)
Ascent within 1,000 km is far higher at HF onset (Atlantic 1.21 vs 0.84 Pa/s, Pacific 1.35 vs 0.93, q 0.003); ascent 750-2,500 km out differs little (Atlantic +0.02, q 0.21; Pacific +0.05, q 0.003); the share with a second centre does not differ (Atlantic +0.04, Pacific -0.03, q 0.14 and 0.33). HF is defined by the gust index, so HF storms differ from storm-force ones in strength by construction.

## Deviations (post hoc; also in `PREREGISTRATION.md`)
1. Centre rule changed before the pull (0.4 Pa/s, 1,000 km apart, 4,000 km -> 0.8 Pa/s, 1,500 km, 3,000 km) because the first rule gave two or more centres at 100% of fixes on a 24-chunk dry run. Original and 0.6 Pa/s rules are sensitivities (`_o`, `_m`) and both are nulls.
2. The Wx sensitivity group is one variable (`a500_x750` plus lags), not a full ring replica of W.
3. The card said "about 22 GB, omega plus a wider box and a second-streak flag". The wind-based second-streak flag was not within that size (+40.6 GB) and was not done; coupling is defined from the ascent field. Reported to Jason with the results.

## Limits
Not a forecast test (perfect prognosis from ERA5). 1.5 deg omega smooths narrow frontal ascent. Omega includes the storm's own latent heating, so the W gain is partly the storm describing itself; the ring and quadrant checks are the guards. No Q-vector, no 0.25 deg, no pull after 2023-01-09. Same fixes as PR 149, so the 18 seasons are reused for a new set of tests (counted above).

## Reproduce
`ERA5_WORK=DIR python3 extract_w.py 4` (about 22 GB, 5 min), `ERA5_WORK=DIR python3 qc_w.py`, `python3 qc_cases_w.py`, `python3 skill_w.py results`, `python3 desc_w.py`. `python3 test_synthetic.py` runs the synthetic check.

## Verification
See `VERIFICATION.md` (fresh Sonnet agent recomputation). Not independently checked: extraction code (decoding, great-circle placement), the QC drawings, the 28-test LOSO tables other than the recomputed ones.
