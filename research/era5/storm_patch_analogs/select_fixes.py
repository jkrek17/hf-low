"""Write results/fixes_patch.csv: the 32,288 fixes of the pre-registration (seasons 2014-15..2021-22, Oct-Apr)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "intensity"))
import model as M
D = M.load(os.path.join(HERE, "..", "intensity", "results", "fixes_2004.csv.gz"), os.path.join(HERE, "..", "intensity", "results", "env_2004.csv.gz"))
S = D[(D.season >= 2014) & (D.season <= 2021) & D.month.isin([10, 11, 12, 1, 2, 3, 4])]
assert len(S) == 32288 and S.track.nunique() == 7379
S[["track", "time", "basin", "season", "lat", "lon"]].to_csv(os.path.join(HERE, "results", "fixes_patch.csv"), index=False)
print(len(S))
