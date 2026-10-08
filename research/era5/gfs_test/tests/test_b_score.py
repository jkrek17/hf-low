"""Synthetic check: the Stage B skill function is ~0 for a null index and positive for a planted one. python tests/test_b_score.py (no network)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "research", "era5", "gfs_test"))
import b_score as B


def make(seed, rr):
    rng = np.random.default_rng(seed)
    seas = np.repeat(np.arange(2004, 2019), 30)
    month = np.tile(np.repeat([10, 11, 12, 1, 2, 3, 4], [5, 4, 5, 4, 4, 4, 4])[:30], 15)
    z = rng.normal(size=450)
    mu = 1.3 * rr ** z
    y = rng.poisson(mu).astype(float)
    prev = np.r_[0, y[:-1]]
    return pd.DataFrame(dict(month=month, lprev=np.log1p(prev))), y, z, seas


def test():
    sk = []
    for s in range(6):
        df, y, z, seas = make(s, 1.0)
        D0, D1 = B.loso_dev(df, y, z, seas); sk.append(B.skill(D0, D1))
    assert abs(np.mean(sk)) < 0.005 and max(sk) < 0.01, sk          # null index: skill about 0 (slightly negative: one extra fitted coefficient)
    df, y, z, seas = make(0, 1.4)
    D0, D1 = B.loso_dev(df, y, z, seas)
    assert B.skill(D0, D1) > 0.03
    print("ok null mean skill %.4f, planted skill %.4f" % (np.mean(sk), B.skill(D0, D1)))


if __name__ == "__main__":
    test()
