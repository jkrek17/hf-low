"""Synthetic checks for Stage A scoring: a perfect predictor gives BSS near 1, a constant climatological predictor BSS 0,
a null G gives A1 skill about 0 and a planted G positive. python tests/test_a_score.py (no network)."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import a_score as S


def test():
    rng = np.random.default_rng(1)
    n = 5000
    seas = rng.integers(2021, 2026, n)
    base = 0.1
    y = (rng.random(n) < base).astype(int)
    pc = np.full(n, base)
    B = S.draws(5)
    perfect, _ = S.bss_stat(np.where(y == 1, 0.999, 0.001), pc, y, seas, S.SEASONS, B)
    assert perfect["bss"] > 0.99 and perfect["ci"][0] > 0.98, perfect
    const, bs = S.bss_stat(pc, pc, y, seas, S.SEASONS, B)
    assert const["bss"] == 0.0 and np.allclose(bs, 0.0), const
    assert abs(S.auc(np.where(y == 1, 1.0, 0.0), y) - 1.0) < 1e-12
    assert abs(S.auc(rng.random(n), y) - 0.5) < 0.05
    thr = S.best_hss_threshold(np.where(y == 1, 0.9, 0.1), y)
    assert abs(thr - 0.9) < 1e-12
    cs = S.cat_scores(S.contingency(np.where(y == 1, True, False), y.astype(bool)))
    assert cs["hss"] == 1.0 and cs["pod"] == 1.0 and cs["far"] == 0.0
    # A1: month harmonics + G; null G -> skill about 0, planted G -> clearly positive
    doy = rng.integers(0, 365, n)
    c, s = np.cos(2 * np.pi * doy / 365.25), np.sin(2 * np.pi * doy / 365.25)
    g = rng.normal(60, 10, n)
    sk = []
    for eff in (0.0, 0.08):
        p = 1 / (1 + np.exp(-(-2.2 + 0.5 * c + eff * (g - 60))))
        yy = (rng.random(n) < p).astype(int)
        pg, pm = S.loso_probs(g, c, s, yy, seas)
        st, _ = S.bss_stat(pg, pm, yy, seas, S.SEASONS, B)
        sk.append(st["bss"])
    assert abs(sk[0]) < 0.01 and sk[1] > 0.03, sk
    # BH
    q = S.bh({"a": 0.001, "b": 0.04, "c": 0.5})
    assert abs(q["a"] - 0.003) < 1e-12 and abs(q["b"] - 0.06) < 1e-12 and q["c"] == 0.5
    print("ok perfect BSS %.3f, constant BSS %.3f, null G %.4f, planted G %.3f" % (perfect["bss"], const["bss"], sk[0], sk[1]))


if __name__ == "__main__":
    test()
