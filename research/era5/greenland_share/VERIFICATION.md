# Verification

A fresh Sonnet agent, given the committed inputs and the plan but not the code, README or results, recomputed claims 1-9
(`verify/verify_numbers.py`).

**Matched:** 9,636 tracks / 1,003 HF / 9,622 with fix rows; correlations -0.693 (days), -0.692 (storms), -0.892 (season
means); share ratios GH 0.870, NAO 1.125, joint 0.893 / 1.039; count ratios (GH alone HF 0.895, all 1.029; joint 0.910,
1.020); tercile shares and cell sizes (terciles cut over the 9,622 storms with fix rows; one cell 7.1 against the
README's 7.2, rounding); first-position longitude +1.958, minimum pressure +1.001, Eady -0.016; logistic GH
coefficients -0.135 / -0.157 (116%) / -0.040 (30%); 500-draw season bootstrap of GH given NAO [-0.194, -0.019].

**Differed, then resolved:** the box logistic (README B10, -0.132). The verifier got -0.107 with 2,011 storms in the box
(the producer also has 2,011). Its script (line 71) entered the box indicator as a covariate rather than as the
outcome, so its value is not a test of the claim. A second method by the producer (statsmodels binomial GLM, same
covariates, box as outcome) gives -0.132. **This one claim is therefore checked by a second method, not independently.**

**Not independently checked:** the 2,000-draw intervals, p and q values, minimum detectable effects, leave-one-season-out
ranges, retained-fraction intervals, and the post hoc table.
