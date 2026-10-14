"""Patch geometry and helpers for the storm-centred pattern test (re-uses hf_env_composites/common.py)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "hf_env_composites"))
import common as _c
from common import *          # WB2, CLIM, LEVELS, G, sample, destination, decode, zarray, fs, WORK, WB2_T0

NP, DKM = 13, 154.3333        # 13 x 13 points, edge = 6 * 154.333 = 926 km = 500 nautical miles
_a = (np.arange(NP) - NP // 2) * DKM
PX, PY = np.meshgrid(_a, _a, indexing="xy")          # PX east, PY north
PR = np.hypot(PX, PY)
PBANG = np.degrees(np.arctan2(-PY, PX))              # with rot = 90 the bearing is east for +x, north for +y
CH = ["a_z500", "a_msl", "a_u250", "a_v250", "z500", "msl"]    # anomalies first, raw (display only) last


# gcsfs does not work through this container's proxy; plain HTTPS to the public bucket does.
import time, requests
_sess = requests.Session()


def _cat(key):
    url = "https://storage.googleapis.com/" + key
    for k in range(6):
        try:
            r = _sess.get(url, timeout=120)
            if r.status_code == 200:
                return r.content
            if r.status_code == 404:
                raise FileNotFoundError(key)
        except FileNotFoundError:
            raise
        except Exception:
            pass
        time.sleep(2 ** k)
    raise RuntimeError(f"fetch failed {key}")


_c.cat = _cat
cat = _cat
