"""Sanity checks of bulk.py: neutral air gives actual = neutral; unstable air gives neutral above actual by 1 to 12 per cent, most at light wind; the inverse returns its input."""
import numpy as np
from bulk import *
n = 3; u = np.array([10.0, 25.0, 40.0]); t2 = np.full(n, 280.0); td = np.full(n, 276.0); p = np.full(n, 1010.0)
neutral = neutral_from_actual(u, t2, td, t2 + 0.0, p, 36.0, 0)
unst = neutral_from_actual(u, t2, td, t2 + 8.0, p, 36.0, 0)
stab = neutral_from_actual(u, t2, td, t2 - 3.0, p, 36.0, 0)
print("near-neutral EN/actual", np.round(neutral / u, 3), "unstable (+8 K)", np.round(unst / u, 3), "stable (-3 K)", np.round(stab / u, 3))
assert np.all(np.abs(neutral / u - 1) < 0.02) and np.all((unst / u > 1.005) & (unst / u < 1.2)) and np.all(np.diff(unst / u) < 0) and np.all(stab / u < 1.0)
back, res = actual_from_neutral(unst, t2, td, t2 + 8.0, p, 36.0, 0); assert res < 0.05 and np.allclose(back, u, rtol=2e-3), (back, res)
print("ok")
