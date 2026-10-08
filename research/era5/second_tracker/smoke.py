"""Smoke test on one chunk each (no comparison with A): detections per field and a plotted sanity count."""
import sys, time
import common as C
for var, det in (("mean_sea_level_pressure", C.detect_m), ("vorticity", C.detect_v)):
    c, _ = C.chunk_of(__import__("datetime").datetime(2015, 1, 10))
    t = time.time(); a = C.fetch(var, c); print(var, a.shape, float(a.min()), float(a.max()), round(time.time() - t, 1), "s", C.BYTES[0] / 1e6, "MB")
    t = time.time(); d = det(a[0]); print(len(d), "detections", round(time.time() - t, 2), "s", d[:3])
