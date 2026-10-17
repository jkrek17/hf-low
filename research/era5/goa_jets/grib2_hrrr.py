"""HRRR-Alaska variant of the repository GRIB2 reader (outlook_replay/grib2_t53.py): same packing 5.3 decoder, grid template 3.20 (polar stereographic, 1299 x 919). Returns the array as stored (rows = j, columns = i); see hrrr_grid.py for the geolocation."""
"""Minimal GRIB2 reader for the one packing the GFS 0.25 degree gust field uses: grid template 3.0 (lat/lon), data
representation template 5.3 (complex packing + spatial differencing, order 1 or 2), no missing values.
Written because eccodes cannot be installed here (the package index is not reachable). Checked against ERA5 gust, see README."""
import struct
import numpy as np


def _smag(b):
    v = int.from_bytes(b, "big"); s = 1 << (8 * len(b) - 1)
    return -(v & (s - 1)) if v & s else v


def _unpack(bits, off, widths, n_each=None):
    """bits: uint8 bit array; read len(widths) values of per-value widths starting at bit offset off."""
    widths = np.asarray(widths, dtype=np.int64); out = np.zeros(len(widths), dtype=np.int64)
    pos = off + np.concatenate([[0], np.cumsum(widths)[:-1]])
    for w in np.unique(widths):
        if w == 0: continue
        sel = np.nonzero(widths == w)[0]
        idx = pos[sel][:, None] + np.arange(w)[None, :]
        out[sel] = bits[idx].astype(np.int64) @ (1 << np.arange(w - 1, -1, -1, dtype=np.int64))
    return out, off + int(widths.sum())


def decode_message(b):
    assert b[:4] == b"GRIB" and b[7] == 2
    p = 16; sec = {}
    while b[p:p + 4] != b"7777":
        L = struct.unpack(">I", b[p:p + 4])[0]; sec[b[p + 4]] = b[p:p + L]; p += L
    s3, s5, s7 = sec[3], sec[5], sec[7]
    assert struct.unpack(">H", s3[12:14])[0] == 20, "grid template (HRRR-Alaska polar stereographic 3.20)"
    ni, nj = struct.unpack(">II", s3[30:38]); scan = s3[64]
    assert (ni, nj) == (1299, 919), (ni, nj, scan)
    npts = struct.unpack(">I", s5[5:9])[0]; tmpl = struct.unpack(">H", s5[9:11])[0]
    if tmpl == 0:        # simple packing (HRRR 10 m UGRD, VGRD)
        R = struct.unpack(">f", s5[11:15])[0]; E = _smag(s5[15:17]); D = _smag(s5[17:19]); nb = s5[19]
        bits = np.unpackbits(np.frombuffer(s7[5:], dtype=np.uint8))[:npts * nb].reshape(npts, nb)
        X = bits.astype(np.int64) @ (1 << np.arange(nb - 1, -1, -1, dtype=np.int64))
        return ((R + X.astype(np.float64) * 2.0 ** E) / 10.0 ** D).astype(np.float32).reshape(nj, ni)
    assert tmpl == 3, "packing template"
    R = struct.unpack(">f", s5[11:15])[0]; E = _smag(s5[15:17]); D = _smag(s5[17:19]); nb = s5[19]
    assert s5[22] == 0, "missing values present"
    NG = struct.unpack(">I", s5[31:35])[0]; wref, wbits = s5[35], s5[36]
    lref = struct.unpack(">I", s5[37:41])[0]; linc = s5[41]; llast = struct.unpack(">I", s5[42:46])[0]; lbits = s5[46]
    order, extra = s5[47], s5[48]
    d = s7[5:]; q = 0; init = []
    for _ in range(order): init.append(_smag(d[q:q + extra])); q += extra
    dmin = _smag(d[q:q + extra]); q += extra
    bits = np.unpackbits(np.frombuffer(d[q:], dtype=np.uint8)); off = 0
    def block(n, w):
        nonlocal off
        v, off = _unpack(bits, off, np.full(n, w, dtype=np.int64)) if w else (np.zeros(n, dtype=np.int64), off)
        off = (off + 7) // 8 * 8
        return v
    gref = block(NG, nb); gw = block(NG, wbits) + wref; gl = block(NG, lbits) * linc + lref; gl[-1] = llast
    assert int(gl.sum()) == npts, (int(gl.sum()), npts)
    vals = np.zeros(npts, dtype=np.int64)
    # per-value width = its group's width; values then offset by the group reference
    wv = np.repeat(gw, gl); raw, off = _unpack(bits, off, wv)
    z = raw + np.repeat(gref, gl)
    if order == 1:
        z[0] = init[0]; z[1:] += dmin; z = np.cumsum(z)
    elif order == 2:
        z[2:] += dmin; z[0] = 0; z[1] = init[1] - init[0]
        z = init[0] + np.cumsum(np.cumsum(z))
    return ((R + z.astype(np.float64) * 2.0 ** E) / 10.0 ** D).astype(np.float32).reshape(nj, ni)
