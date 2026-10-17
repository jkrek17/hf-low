"""COARE 3.6 (pycoare) wrappers: actual 10 m wind <-> stress-equivalent (neutral) 10 m wind. Winds in kt in and out. See PREREGISTRATION.md (Method > Bulk conversion)."""
import numpy as np
from pycoare import coare_36
KT = 1.94384449
def rh_from(t_c, td_c):
    """relative humidity (%) from 2 m temperature and dewpoint, Magnus formula"""
    e = lambda x: np.exp(17.625 * x / (243.04 + x))
    return np.clip(100.0 * e(td_c) / e(t_c), 1.0, 100.0)
def neutral_from_actual(u_kt, t2_k, td_k, ts_k, p_hpa, lat, jcool, chunk=250000):
    """stress-equivalent neutral 10 m wind (kt) from the actual 10 m wind, with 2 m air temperature and dewpoint, sea (or skin) temperature, pressure"""
    u = np.asarray(u_kt, float) / KT; t2 = np.asarray(t2_k, float) - 273.15; td = np.asarray(td_k, float) - 273.15; ts = np.asarray(ts_k, float) - 273.15
    p = np.asarray(p_hpa, float); la = np.broadcast_to(np.asarray(lat, float), u.shape); out = np.full(u.shape, np.nan)
    for i in range(0, len(u), chunk):
        s = slice(i, i + chunk)
        c = coare_36(u[s], t=t2[s], rh=rh_from(t2[s], td[s]), zu=10, zt=2, zq=2, ts=ts[s], p=p[s], lat=la[s], jcool=jcool, rs=100, rl=310, zi=600)
        out[s] = np.asarray(c.velocities.u_n_rf, float).reshape(-1)
    return out * KT
def actual_from_neutral(un_kt, t2_k, td_k, ts_k, p_hpa, lat, jcool, nit=4):
    """invert by fixed point: U <- Un * U / En(U); converges in a few steps because En/U varies slowly with U"""
    un = np.asarray(un_kt, float); u = un.copy()
    for _ in range(nit):
        en = neutral_from_actual(u, t2_k, td_k, ts_k, p_hpa, lat, jcool); u = un * u / en
    return u, np.nanmax(np.abs(neutral_from_actual(u, t2_k, td_k, ts_k, p_hpa, lat, jcool) - un))
def neutral_from_stress(tau, t2_k, td_k, ts_k, p_hpa, lat, jcool, nit=5, chunk=250000):
    """stress-equivalent neutral 10 m wind (kt) from a surface stress tau (N/m2): u* = sqrt(tau/rho); U_N = u* / sqrt(Cdn(U_N)), iterated, Cdn from COARE 3.6 at the current guess"""
    tau = np.asarray(tau, float); t2 = np.asarray(t2_k, float) - 273.15; td = np.asarray(td_k, float) - 273.15; ts = np.asarray(ts_k, float) - 273.15
    p = np.asarray(p_hpa, float); la = np.broadcast_to(np.asarray(lat, float), tau.shape); rho = p * 100.0 / (287.05 * (t2 + 273.15) * (1 + 0.608 * 0.622 * 6.11 * np.exp(17.625 * td / (243.04 + td)) / p))
    us = np.sqrt(np.maximum(tau, 1e-6) / rho); u = us / 0.04; out = np.full(tau.shape, np.nan)
    for i in range(0, len(u), chunk):
        s = slice(i, i + chunk); uu = u[s].copy()
        for _ in range(nit):
            c = coare_36(uu, t=t2[s], rh=rh_from(t2[s], td[s]), zu=10, zt=2, zq=2, ts=ts[s], p=p[s], lat=la[s], jcool=jcool, rs=100, rl=310, zi=600)
            uu = us[s] / np.sqrt(np.asarray(c.transfer_coefficients.cdn_rf, float).reshape(-1))
        out[s] = uu
    return out * KT
