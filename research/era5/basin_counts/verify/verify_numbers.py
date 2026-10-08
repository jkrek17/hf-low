import numpy as np, pandas as pd, datetime as dt
R='/home/claude/hf-low/research/era5/hemispheric/results/'
d=pd.read_csv(R+'weekly_table.csv.gz')
rng=np.random.default_rng(1)
def get(a,p,s0,s1):
    x=d[(d.season>=s0)&(d.season<=s1)].sort_values(['season','week'])
    assert len(x)==660
    return x[a].values.reshape(22,30).astype(float),x[p].values.reshape(22,30).astype(float),x
def ds(M): return M-M.mean(0)
def corr(a,b): return np.corrcoef(a.ravel(),b.ravel())[0,1]
def perm_p(A,P,n=200000):
    a=A.ravel()-A.mean();a=A-A.mean();p=P-P.mean()
    G=a@p.T; den=np.sqrt((a**2).sum()*(p**2).sum()); r0=np.trace(G)/den
    cnt=0
    for _ in range(n//10000):
        pm=np.array([rng.permutation(22) for _ in range(10000)])
        rs=G[np.arange(22)[None,:],pm].sum(1)/den
        cnt+=(rs<=r0+1e-12).sum()
    return r0,(1+cnt)/(1+n)
S={'ARCH':('y_atl','y_pac',2004,2025),'PRX':('pAhf_atl','pAhf_pac',2004,2025),'DEP':('pAdepth_atl','pAdepth_pac',1979,2000)}
mon=np.array([(dt.date(2001,10,1)+dt.timedelta(days=7*k+3)).month for k in range(30)])
D=np.array([[m==mm for mm in (10,11,12,1,2,3,4)] for m in mon],float) # 30x7
Dall=np.tile(D,(22,1))
def resid(y,X): 
    b=np.linalg.lstsq(X,y,rcond=None)[0]; return y-X@b
for k,(a,p,s0,s1) in S.items():
    A,P,x=get(a,p,s0,s1); Ad,Pd=ds(A),ds(P)
    r,pv=perm_p(Ad,Pd); print(k,'primary r %.4f p %.4f'%(r,pv))
    print(' S1 %.4f'%corr(resid(A.ravel(),Dall).reshape(22,30),resid(P.ravel(),Dall).reshape(22,30)))
    print(' S2a %.4f'%np.corrcoef(A.sum(1),P.sum(1))[0,1])
    X=np.c_[np.ones(660),x.nao.values,x.pna.values]
    print(' S3 %.4f'%corr(resid(Ad.ravel(),X),resid(Pd.ravel(),X)))
    print(' halves %.4f %.4f'%(corr(Ad[:11],Pd[:11]),corr(Ad[11:],Pd[11:])))
    print(' both zero %.4f'%((A==0)&(P==0)).mean())
    for b in (5,10):
        Ab=A.reshape(22,30//b,b).sum(2);Pb=P.reshape(22,30//b,b).sum(2)
        print(' block%d %.4f'%(b,corr(ds(Ab),ds(Pb))))
    if k=='ARCH':
        print(' var/mean atl %.3f pac %.3f'%(Ad.var()*22/21/A.mean(), Pd.var()*22/21/P.mean()))
        print(' var/mean (per-cell var)', (Ad**2).sum()/(22*30-30)/A.mean(), (Pd**2).sum()/(22*30-30)/P.mean())
        ia=pd.read_csv(R+'oos_index_atl.csv').idx.values; ip=pd.read_csv(R+'oos_index_pac.csv').idx.values
        ri=np.corrcoef(ia,ip)[0,1]; print(' index r %.4f ratio %.4f'%(ri,r/ri))
