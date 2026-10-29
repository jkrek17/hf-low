"""Joint (GEFS members reaching T, ENS members reaching T) case tables for winter 2025-26, from pulled cell maxima and truth.
Per (run date, lead, basin, threshold, truth): n[kg, ke] cases and y[kg, ke] events. Same truth rules as ../build_counts.py (imported).
usage: build_joint.py OUT.npz   (ENS_WORK must hold gefs/, ens/, ascat/)"""
import os, sys, datetime as dt
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, ".."))
import ens_fetch as F, build_counts as B

NG, NE = 31, 51


def joint(kg, ke, y):
    n = np.zeros((NG + 1, NE + 1), np.int64); e = np.zeros_like(n)
    np.add.at(n, (kg, ke), 1); np.add.at(e, (kg, ke), y.astype(np.int64)); return n, e


def build(fixes, season=2025, leads=(24, 48, 72, 96)):
    keys, N, Y = [], [], []
    for ymd in F.run_dates(season):
        run = dt.datetime.strptime(ymd, "%Y%m%d")
        for lead in leads:
            pg = os.path.join(F.WORK, "gefs", "%s_%03d.npz" % (ymd, lead)); pe = os.path.join(F.WORK, "ens", "%s_%03d.npz" % (ymd, lead))
            if not (os.path.exists(pg) and os.path.exists(pe)): continue
            zg, ze = np.load(pg), np.load(pe); V = run + dt.timedelta(hours=lead) + dt.timedelta(days=int(os.environ.get('PLACEBO_DAYS', '0')))   # truth date; shifted only in the placebo run
            tz = os.path.join(F.WORK, "ascat", V.strftime("%Y%m%d") + ".npz"); a = np.load(tz) if os.path.exists(tz) else None
            fx = fixes[(fixes.t >= V - dt.timedelta(hours=6)) & (fixes.t <= V + dt.timedelta(hours=6))]
            for basin in F.BOX:
                mg, me = zg[basin], ze[basin]
                if a is not None:
                    wnd, cov, oc = a[basin + "_wind"], a[basin + "_cov"], a[basin + "_oc"]; ok = (oc >= 8) & (cov >= 0.5 * oc)
                    for T in (48, 64):      # 64 kt vs ASCAT is not scored (saturation); only 48 here
                        if T != 48: continue
                        n, e = joint((mg >= 2 * T).sum(0)[ok], (me >= 2 * T).sum(0)[ok], (wnd >= 2 * T)[ok])
                        keys.append((ymd, lead, basin, T, "ascat")); N.append(n); Y.append(e)
                cr, cc, idx, LA, LO = B.neighbours(basin); fg, fe = mg.reshape(NG, -1), me.reshape(NE, -1)
                clat, clon = LA[cr, cc], LO[cr, cc]
                oc = (a[basin + "_oc"][cr, cc] >= 8) if a is not None else np.ones(len(cr), bool)
                fb = fx[fx.basin == basin]; T = 64
                kg = np.stack([fg[:, ix].max(1) for ix in idx], 1); ke_ = np.stack([fe[:, ix].max(1) for ix in idx], 1)
                f = fb[fb.Category.isin(B.CATS[T])]
                if len(f):
                    d = B.gc(clat[:, None], clon[:, None], f.lat.values[None, :], f.lon.values[None, :]); y = (d <= B.RAD).any(1)
                else: y = np.zeros(len(cr), bool)
                n, e = joint((kg >= 2 * T).sum(0)[oc], (ke_ >= 2 * T).sum(0)[oc], y[oc])
                keys.append((ymd, lead, basin, T, "opc")); N.append(n); Y.append(e)
    return np.array([list(map(str, k)) for k in keys]), np.stack(N), np.stack(Y)


if __name__ == "__main__":
    K, N, Y = build(B.opc_fixes())
    np.savez_compressed(sys.argv[1], keys=K, n=N, y=Y); print(len(K), "case groups; cases", int(N.sum()), "events", int(Y.sum()))
