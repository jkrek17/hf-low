"""Does adding the four regions (see diagrams.py) add skill to state, state+Hart or the full model? Post hoc, descriptive.
ERA5 proxy, pipeline A, leave-one-season-out, season-bootstrap 90 % intervals on the gain. Run from the repo root:
python3 research/era5/bomb_env/pair_regions.py > research/era5/bomb_env/results/pair_regions.txt
Region cuts are the upper thirds over all fixes (a mild look-ahead on two quantiles).
"""
import sys,numpy as np,pandas as pd
sys.path.insert(0,"research/era5/bomb_env"); sys.path.insert(0,"research/era5/intensity")
import score as S, model as M
rng=np.random.default_rng(1)
D=S.load()
hv,he=np.nanquantile(D.vadv500,2/3),np.nanquantile(D.eady,2/3)
D["reg"]=(D.vadv500>=hv).astype(int)*2+(D.eady>=he).astype(int)
for k in (1,2,3): D[f"r{k}"]=(D.reg==k).astype(float)
REG=["r1","r2","r3"]
S.SETS={"hart":M.SETS["hart"],"hart_reg":M.SETS["hart"]+REG,"full":M.SETS["full"],"full_reg":M.SETS["full"]+REG,"state_reg":M.STATE+REG,"state":M.STATE}
P4=S.loso(D,"gcls",4); y=D.gcls.values
L4={k:S.rps(P4[k],y,4) for k in P4}
print("outlook RPSS vs clim")
for k in L4: print(" ",k,round(S.skill(L4[k],L4["clim"]),4))
for a,b in (("hart_reg","hart"),("full_reg","full"),("state_reg","state"),("hart_reg","full")):
    print(" gain",a,"over",b,[round(float(x),4) for x in S.boot_gain(D.season.values,L4[a],L4[b],L4["clim"],rng)])
R=D[D.cls>=0].reset_index(drop=True); R["rd"]=(R.cls==4).astype(int)
P2=S.loso(R,"rd",2); L2={k:S.brier(P2[k][:,1],R.rd.values) for k in P2}
print("rd BSS")
for k in L2: print(" ",k,round(S.skill(L2[k],L2["clim"]),4))
for a,b in (("hart_reg","hart"),("full_reg","full"),("state_reg","state"),("hart_reg","full")):
    print(" gain",a,"over",b,[round(float(x),4) for x in S.boot_gain(R.season.values,L2[a],L2[b],L2["clim"],rng)])
print(M.HART)
hb=[c for c in M.HART][0]
print(D.groupby("reg")[M.HART].median().round(2))
