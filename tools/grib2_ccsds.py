"""Minimal GRIB2 reader for ECMWF open data (IFS, AIFS, ENS, AIFS ENS) 0.25 degree fields: grid template 3.0, data
representation template 5.42 (CCSDS / AEC), through the system libaec via ctypes. Written because eccodes cannot be
installed here. Grid is 1440 x 721, scanned +i, -j, starting at longitude -180 (stored as 180); decode() returns it
rolled to 0..359.75 and north-up, the same order as tools/grib2_t53.py gives for GFS and GEFS.
Ported idea: jkrek17/awips-tools marine_wind/pipeline/fetch.py (decode with eccodes). Read-only reference."""
import ctypes, ctypes.util, struct
import numpy as np

_AEC = None
def _lib():
    global _AEC
    if _AEC is None:
        name = ctypes.util.find_library("aec") or "libaec.so.0"
        _AEC = ctypes.CDLL(name)
    return _AEC

class _Stream(ctypes.Structure):
    _fields_ = [("next_in", ctypes.c_void_p), ("avail_in", ctypes.c_size_t), ("total_in", ctypes.c_size_t),
                ("next_out", ctypes.c_void_p), ("avail_out", ctypes.c_size_t), ("total_out", ctypes.c_size_t),
                ("bits_per_sample", ctypes.c_uint), ("block_size", ctypes.c_uint), ("rsi", ctypes.c_uint),
                ("flags", ctypes.c_uint), ("state", ctypes.c_void_p), ("reserved_state", ctypes.c_int)]

def _smag(b):
    v = int.from_bytes(b, "big"); s = 1 << (8 * len(b) - 1)
    return -(v & (s - 1)) if v & s else v

def _aec_decode(data, nbits, block, rsi, flags, npts):
    nbytes = (nbits + 7) // 8
    if nbytes == 3 and not flags & 2: nbytes = 4
    out = ctypes.create_string_buffer(npts * nbytes)
    src = ctypes.create_string_buffer(bytes(data), len(data))
    s = _Stream(); s.next_in = ctypes.cast(src, ctypes.c_void_p); s.avail_in = len(data)
    s.next_out = ctypes.cast(out, ctypes.c_void_p); s.avail_out = npts * nbytes
    s.bits_per_sample = nbits; s.block_size = block; s.rsi = rsi; s.flags = flags
    lib = _lib()
    rc = lib.aec_buffer_decode(ctypes.byref(s))
    if rc != 0 or s.total_out != npts * nbytes: raise ValueError("libaec rc=%d out=%d want %d" % (rc, s.total_out, npts * nbytes))
    dt = ">" if flags & 4 else "<"
    if nbytes == 3:
        a = np.frombuffer(out.raw, dtype=np.uint8).reshape(npts, 3).astype(np.uint32)
        return ((a[:, 0] << 16) | (a[:, 1] << 8) | a[:, 2]) if flags & 4 else ((a[:, 2] << 16) | (a[:, 1] << 8) | a[:, 0])
    return np.frombuffer(out.raw, dtype=np.dtype("%su%d" % (dt, nbytes))).astype(np.uint64)

def decode_message(b):
    assert b[:4] == b"GRIB" and b[7] == 2
    p = 16; sec = {}
    while b[p:p + 4] != b"7777":
        L = struct.unpack(">I", b[p:p + 4])[0]; sec[b[p + 4]] = b[p:p + L]; p += L
    s3, s5, s7 = sec[3], sec[5], sec[7]
    assert struct.unpack(">H", s3[12:14])[0] == 0, "grid template"
    ni, nj = struct.unpack(">II", s3[30:38]); scan = s3[71]
    assert (ni, nj) == (1440, 721) and scan == 0, (ni, nj, scan)
    lat0 = _smag(s3[46:50]) / 1e6; lon0 = struct.unpack(">I", s3[50:54])[0] / 1e6
    npts = struct.unpack(">I", s5[5:9])[0]; assert struct.unpack(">H", s5[9:11])[0] == 42, "packing template"
    R = struct.unpack(">f", s5[11:15])[0]; E = _smag(s5[15:17]); D = _smag(s5[17:19]); nb = s5[19]
    flags, block, rsi = s5[21], s5[22], struct.unpack(">H", s5[23:25])[0]
    if nb == 0: x = np.zeros(npts)
    else: x = _aec_decode(s7[5:], nb, block, rsi, flags, npts).astype(np.float64)
    v = ((R + x * 2.0 ** E) / 10.0 ** D).astype(np.float32).reshape(nj, ni)
    if lat0 < 0: v = v[::-1]
    if lon0 > 100: v = np.roll(v, ni // 2, axis=1)
    return v
