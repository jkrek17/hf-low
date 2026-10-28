"""No-network test: plant an overconfident ensemble and a reliable one; the slope test must flag the first and not the second;
isotonic calibration must remove the planted bias; BH must behave. Run: python3 test_stats.py"""
import numpy as np, stats as S

def make(N, true_of_p, seasons=6, per=40000, seed=1):
    rng = np.random.default_rng(seed); nb, yb = {}, {}
    for s in range(seasons):
        # latent frequency of the forecast fraction: most cases k=0, some higher
        k = np.minimum(rng.geometric(0.5, per) - 1, N) * (rng.random(per) < 0.12)
        k = np.where(rng.random(per) < 0.9, 0, np.minimum(rng.geometric(0.25, per), N))
        p = true_of_p(k / N); y = rng.random(per) < p
        nb[s] = np.bincount(k, minlength=N + 1); yb[s] = np.bincount(k, weights=y, minlength=N + 1).astype(int)
    return nb, yb

N = 31
rel = make(N, lambda f: np.clip(f, 0.002, 1))                 # reliable: observed rate equals fraction
over = make(N, lambda f: np.clip(0.02 + 0.5 * f, 0.002, 1))   # overconfident: high fractions verify far less often
r1 = S.summarize(*rel, N, nboot=400); r2 = S.summarize(*over, N, nboot=400)
print("reliable slope %.2f %s p=%.3f gain %.4f" % (r1["slope"], r1["slope_ci"], r1["p_overconf"], r1["cal_gain"]))
print("overconf slope %.2f %s p=%.3f gain %.4f %s" % (r2["slope"], r2["slope_ci"], r2["p_overconf"], r2["cal_gain"], r2["cal_gain_ci"]))
assert 0.85 < r1["slope"] < 1.15 and r1["p_overconf"] > 0.05
assert r2["slope"] < 0.7 and r2["p_overconf"] < 0.05 and r2["cal_gain"] > 0 and r2["cal_gain_ci"][0] > 0
assert S.bh([0.001, 0.2, 0.9]) == [True, False, False]
n = np.array([5, 5, 5]); y = np.array([3, 1, 2]); assert np.all(np.diff(S.iso_fit(n, y)) >= 0)
print("ok")
