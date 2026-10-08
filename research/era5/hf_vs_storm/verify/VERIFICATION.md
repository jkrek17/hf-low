# Verification of HF-vs-SF ERA5-proxy numbers

Independent recomputation (python3 -I verify_numbers.py). Claim 3 SF means: per-month means weighted by HF_onset non-NaN counts per month (peak-month stratum); HF = plain mean. Tolerance 5% of quoted value (claim 3), 0.5 m z500 / 0.05 m/s u250 (claim 4).

| claim | recomputed | quoted | match? |
|---|---|---|---|
| 1. atl HF n | 1104 | 1104 | YES |
| 1. atl SF n | 4032 | 4032 | YES |
| 1. pac HF n | 896 | 896 | YES |
| 1. pac SF n | 4104 | 4104 | YES |
| 2. anchor rows | 2400 | 2400 | YES |
| 2. rows per anchor | {'SF_peak': 800, 'HF_onset': 800, 'HF_peak': 800} | 800 each | YES |
| 2. distinct times (times_storm.csv) | 1894 | 1894 | YES |
| 2. distinct times (stats) | 1894 | 1894 | YES |
| 2. max |g800-g800_cat| (n=2400) | 0.1 | <=0.5 | YES |
| 2. n g800_cat>0 with |d|>0.5 | 0 | 0 | YES |
| 3. atl pc hf (nHF=400,nSF=400) | 971.05 | 971.05 | YES |
| 3. atl pc sf (nHF=400,nSF=400) | 983.5262 | 983.5262 | YES |
| 3. atl pc diff (nHF=400,nSF=400) | -12.4763 | -12.4763 | YES |
| 3. atl msl_grad hf (nHF=400,nSF=400) | 4.2804 | 4.2804 | YES |
| 3. atl msl_grad sf (nHF=400,nSF=400) | 2.6819 | 2.6819 | YES |
| 3. atl msl_grad diff (nHF=400,nSF=400) | 1.5984 | 1.5984 | YES |
| 3. atl gmax_r hf (nHF=400,nSF=400) | 302.185 | 302.185 | YES |
| 3. atl gmax_r sf (nHF=400,nSF=400) | 407.3625 | 407.3625 | YES |
| 3. atl gmax_r diff (nHF=400,nSF=400) | -105.1775 | -105.1775 | YES |
| 3. atl g48_rmax hf (nHF=400,nSF=400) | 945.0325 | 945.0325 | YES |
| 3. atl g48_rmax sf (nHF=400,nSF=400) | 786.8175 | 786.8175 | YES |
| 3. atl g48_rmax diff (nHF=400,nSF=400) | 158.215 | 158.215 | YES |
| 3. atl d2m_500 hf (nHF=400,nSF=400) | 1.7185 | 1.7185 | YES |
| 3. atl d2m_500 sf (nHF=400,nSF=400) | 0.0593 | 0.0593 | YES |
| 3. atl d2m_500 diff (nHF=400,nSF=400) | 1.6592 | 1.6592 | YES |
| 3. atl gust_factor hf (nHF=400,nSF=400) | 1.5088 | 1.5088 | YES |
| 3. atl gust_factor sf (nHF=400,nSF=400) | 1.4647 | 1.4647 | YES |
| 3. atl gust_factor diff (nHF=400,nSF=400) | 0.0441 | 0.0441 | YES |
| 3. pac pc hf (nHF=400,nSF=400) | 973.8288 | 973.8288 | YES |
| 3. pac pc sf (nHF=400,nSF=400) | 984.5025 | 984.5025 | YES |
| 3. pac pc diff (nHF=400,nSF=400) | -10.6738 | -10.6738 | YES |
| 3. pac msl_grad hf (nHF=400,nSF=400) | 4.4894 | 4.4894 | YES |
| 3. pac msl_grad sf (nHF=400,nSF=400) | 2.9208 | 2.9208 | YES |
| 3. pac msl_grad diff (nHF=400,nSF=400) | 1.5686 | 1.5686 | YES |
| 3. pac gmax_r hf (nHF=400,nSF=400) | 232.6075 | 232.6075 | YES |
| 3. pac gmax_r sf (nHF=400,nSF=400) | 333.0875 | 333.0875 | YES |
| 3. pac gmax_r diff (nHF=400,nSF=400) | -100.48 | -100.48 | YES |
| 3. pac g48_rmax hf (nHF=400,nSF=400) | 895.1775 | 895.1775 | YES |
| 3. pac g48_rmax sf (nHF=400,nSF=400) | 741.415 | 741.415 | YES |
| 3. pac g48_rmax diff (nHF=400,nSF=400) | 153.7625 | 153.7625 | YES |
| 3. pac d2m_500 hf (nHF=400,nSF=400) | 6.7707 | 6.7707 | YES |
| 3. pac d2m_500 sf (nHF=400,nSF=400) | 2.5751 | 2.5751 | YES |
| 3. pac d2m_500 diff (nHF=400,nSF=400) | 4.1956 | 4.1956 | YES |
| 3. pac gust_factor hf (nHF=400,nSF=400) | 1.5207 | 1.5207 | YES |
| 3. pac gust_factor sf (nHF=400,nSF=400) | 1.4827 | 1.4827 | YES |
| 3. pac gust_factor diff (nHF=400,nSF=400) | 0.0381 | 0.0381 | YES |
| 4. atl z500 L n=1079/3759 (quoted n 1080/3773) | diff -7.198 hf -7.898 sf -0.700 | diff -7.171 hf -7.872 sf -0.701 | YES |
| 4. atl z500 k0 n=1079/3759 (quoted n 1080/3773) | diff -21.985 hf -27.077 sf -5.092 | diff -21.953 hf -27.042 sf -5.090 | YES |
| 4. atl u250 L n=1079/3759 (quoted n 1080/3773) | diff 1.121 hf 1.001 sf -0.120 | diff 1.119 hf 1.000 sf -0.119 | YES |
| 4. atl u250 k0 n=1079/3759 (quoted n 1080/3773) | diff 1.943 hf 2.115 sf 0.172 | diff 1.942 hf 2.115 sf 0.173 | YES |
| 4. pac z500 L n=874/3874 (quoted n 875/3883) | diff -5.837 hf -6.267 sf -0.430 | diff -5.817 hf -6.240 sf -0.424 | YES |
| 4. pac z500 k0 n=874/3874 (quoted n 875/3883) | diff -9.019 hf -9.671 sf -0.652 | diff -9.090 hf -9.727 sf -0.637 | YES |
| 4. pac u250 L n=874/3874 (quoted n 875/3883) | diff 0.820 hf 0.837 sf 0.017 | diff 0.821 hf 0.837 sf 0.016 | YES |
| 4. pac u250 k0 n=874/3874 (quoted n 875/3883) | diff 1.311 hf 1.473 sf 0.162 | diff 1.314 hf 1.474 sf 0.160 | YES |
## Notes
- All claims 1-4 match. Claims 1-3 reproduce to the 4th decimal (exact).
- Claim 4: recomputed sample is slightly smaller than the quoted n (HF 1078 vs 1080 atl, 872 vs 875 pac; SF 3759 vs 3773 atl, 3874 vs 3883 pac). Cause not fully identified: 26 atl HF onsets fall in Jun-Aug (1104-26=1078, my count); the original apparently excludes by a different rule (excluding on peak-month 6-8 plus onset-month 6-7 gives 1079/874). The remaining SF gap (about 14 atl, 9 pac) is probably a different handling of missing lags (e.g. per-lag NaN-skipping in the lead window). Effect on results: diffs within 0.04 m (atl z500), 0.005 m/s (u250); worst case pac z500 k0 0.09 m. Well inside the stated tolerance.
- Pac z500 k0 recompute differs by 0.09 m (-9.09 quoted vs -9.00 here); u250 diffs <= 0.01 m/s.
- Max |g800 - g800_cat| over the 2400 anchors with g800_cat>0 is 0.1 kt (all 2400 have g800_cat>0).

## NOT checked
- Bootstrap p values, CIs (lo/hi), standard errors, BH q values, MDE (mde80, mde_std), sd/d_std
- Secondary families (P1b lags k=1,2,4,7 and mslp/sst/tcwv box tests), C2 comparison, g48_right, S1/S2/S3
- Pixel maps / composites (maps_*.npz, comp_*.npz), series.csv, figures
- Claim that recomputed g800 was re-derived from raw ERA5 (only the stored g800 vs g800_cat consistency was checked), the raw 0.25 deg extraction itself, heading/rotation
- Group size claims for 'matched' 400/400 random draw (seed) selection; sizes_storm.txt
