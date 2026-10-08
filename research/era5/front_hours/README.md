# Does SST-front strength off Japan track Pacific HF-centre hours? (RA-16, ERA5 proxy)

Plan, committed before any index or outcome: [PREREGISTRATION.md](PREREGISTRATION.md) (commit `d55ee7c`). Code: `sst_front.py`, `analysis.py`, `power.py`. Results: `results/` (`summary.txt`, `results.csv`, `season_table.csv`, `power.txt`). Verification: [VERIFICATION.md](VERIFICATION.md).
**Pipeline A (`research/era5/hf_history`), an ERA5 proxy**; ERA5 0.25 degree SST. Box 40-45N 160-170E, tropical-cyclone-linked events included.

## Answer
**Cannot tell. The point estimate is weakly positive (Spearman rho +0.22, 95% interval -0.26 to +0.61, permutation p 0.36, n = 19 seasons 2007-2025), but the study could only have detected a correlation of about 0.60 (0.63 by simulation).** A front-strength effect up to about 0.6 in rank correlation is not excluded, and only a correlation above about 0.7 would have survived the 11-test correction. In ratio terms: 1.08 times more box fixes per SD of front strength (robust 95% 0.93-1.26). 0 of 11 tests have p < 0.05 or BH q < 0.10 (smallest q 0.49).

## Numbers (autumn SON front index, lead of 0-3 months before the Dec-Mar peak)
| id | test | rho | 95% | p | q |
|---|---|---|---|---|---|
| P1 | FRONT(SON) vs hours, 2007-25 | +0.22 | -0.26, +0.61 | 0.36 | 0.49 |
| S1 | 2004-25, product dummy removed (n=22) | +0.29 | -0.16, +0.64 | 0.18 | 0.49 |
| S2 | Kuroshio Extension band only | -0.24 | -0.62, +0.24 | 0.33 | 0.49 |
| S3 | Oyashio band only | +0.20 | -0.28, +0.60 | 0.40 | 0.49 |
| S4 | concurrent DJF (contrast) | +0.10 | -0.37, +0.53 | 0.69 | 0.75 |
| S5 | deepening-phase hours | +0.21 | -0.27, +0.61 | 0.39 | 0.49 |
| S6 | mature/decay-phase hours | -0.08 | -0.51, +0.39 | 0.75 | 0.75 |
| S7 | storm count | +0.22 | -0.26, +0.61 | 0.36 | 0.49 |
| S8 | OPC archive box fixes | +0.27 | -0.22, +0.65 | 0.26 | 0.49 |
| S9 | year trend removed | +0.22 | -0.28, +0.62 | 0.37 | 0.49 |
| S10 | rho(S5) - rho(S6) | +0.29 | -0.39, +0.84 | 0.39 | 0.49 |

Box hours per season 2007-25: mean 47.1, SD 26.0, range 6-114; 5.1 storms per season. Power (planted effect, marginal of H): 80% at latent rho about 0.63; 38% at 0.4; 21% at 0.3.
The phase prediction (S5 positive, S6 about zero, S10 positive) has the right signs but S10 is not significant, so it is not supported.

## Caveats
- Autocorrelation: FRONT_SON 0.50, hours -0.50 lag-1; the front rises in recent seasons (2024 SON 3.6 against 2.5-3.2 earlier), handled only by S9.
- ERA5 SST product switched in Sept 2007; primary sample starts after it. The 2004-06 values show no visible break.
- One season can be dominated by one slow storm (max single-storm share 1.00, median 0.33); S7 and the rank statistic limit this.
- Not independently checked: see VERIFICATION.md. Post hoc: `power.py` simulation added after outcomes (marginal only). No other deviations.
- Pull: 1.62 GB of ERA5 SST, 132 months x 9 samples. Looks: 21st at 2015-25, none pre-2001 (`research/era5/looks/ra16_front.log`).

## Reproduce
    python3 research/era5/front_hours/sst_front.py work/front_monthly.csv   # 1.6 GB, committed copy in results/
    python3 -I research/era5/front_hours/analysis.py research/era5/front_hours/results/front_monthly.csv research/era5/front_hours/results
    python3 -I research/era5/front_hours/power.py research/era5/front_hours/results/season_table.csv
