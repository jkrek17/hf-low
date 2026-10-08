"""PR 64 channels.py (24 primary tests, decomposition) on a derived track table. chanlib.TRACKS is replaced; nothing else changes.
usage: run_channels.py ALL_TRACKS_VARIANT.csv.gz OUTDIR [NPERM] [NBOOT]
"""
import os, sys
ERA = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(ERA, "hem_channels"))
import chanlib as C  # noqa: E402
C.TRACKS = os.path.abspath(sys.argv[1])
sys.argv = ["channels.py"] + sys.argv[2:]
import channels  # noqa: E402
channels.main()
