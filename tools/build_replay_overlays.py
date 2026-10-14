#!/usr/bin/env python3
"""Pull the GFS context overlays (lows, trough axes, jet cores, Coupled flags; see tools/outlook_overlays.py) for every cycle on the Outlook replay page.

    python3 tools/build_replay_overlays.py [--payload data/outlook_replay/replay_payload.json] [--out data/outlook_replay/replay_overlays.json]

Writes {cycle 'YYYYMMDD': [overlay for day 1..4 or null]}; resumable (cycles already in the output are skipped). About 44 GRIB records of
0.6-0.9 MB per cycle, roughly 4.5 GB streamed for the 126 cycles, nothing kept but the small result. tools/build_outlook_replay.py adds the result
to each case as "ov" (append-only). Standard library plus numpy."""
import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
import outlook_overlays as V


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--payload", default=os.path.join(ROOT, "data", "outlook_replay", "replay_payload.json"))
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "outlook_replay", "replay_overlays.json"))
    a = ap.parse_args()
    cycles = [c["cycle"] for c in json.load(open(a.payload))["cases"]]
    done = json.load(open(a.out)) if os.path.exists(a.out) else {}
    t0 = time.time()
    for k, cy in enumerate(cycles):
        if cy in done and all(x is not None for x in done[cy]):
            continue
        done[cy] = V.for_cycle(cy)
        json.dump(done, open(a.out, "w"), separators=(",", ":"))
        print(f"{k + 1}/{len(cycles)} {cy} {[x is not None for x in done[cy]]} {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
