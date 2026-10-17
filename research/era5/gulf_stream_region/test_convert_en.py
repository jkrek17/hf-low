"""Check the equivalent-neutral conversion has the size it should: at 12 m/s the correction for dT = 7.5 C versus 1.5 C is about 3 points, not 0.5 (the first version's error)."""
import numpy as np
from region_lib import convert_en
v = 12 * 1.944; t2 = 275.0
c = lambda dT: convert_en(v, t2, t2 + dT) / v
d = 100 * (c(1.5) - c(7.5))
assert 2.0 < d < 5.0, d
print("ok", round(d, 2))
