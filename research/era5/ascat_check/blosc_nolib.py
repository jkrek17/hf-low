"""Decode the era5-cache blobs (Blosc1 container, zstd codec, bitshuffle, typesize 1) without numcodecs.

The 2026-10-14 environment has no PyPI access, so numcodecs cannot be installed. This reads the Blosc1 header, decompresses
each block with the system libzstd through ctypes, and undoes bitshuffle in numpy. Only the layout the cache uses is
supported (typesize 1; zstd or raw blocks); anything else raises. Validated against the stored pipeline A index (R0)."""
import ctypes, ctypes.util, struct
import numpy as np

_z = None


def _zstd():
    global _z
    if _z is None:
        name = ctypes.util.find_library("zstd") or "libzstd.so.1"
        _z = ctypes.CDLL(name)
        _z.ZSTD_decompress.restype = ctypes.c_size_t
        _z.ZSTD_decompress.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_char_p, ctypes.c_size_t]
        _z.ZSTD_isError.argtypes = [ctypes.c_size_t]
    return _z


def _unbitshuffle(buf, n):
    """bitshuffle for 1-byte elements: bit-plane b holds bit b of all elements, packed little-endian."""
    m = n - n % 8
    out = np.empty(n, np.uint8)
    if m:
        planes = np.frombuffer(buf, np.uint8, m).reshape(8, m // 8)
        bits = np.unpackbits(planes, axis=1, bitorder="little")          # (8, m)
        out[:m] = (bits.astype(np.uint8) << np.arange(8, dtype=np.uint8)[:, None]).sum(0, dtype=np.uint8)
    out[m:] = np.frombuffer(buf, np.uint8, n - m, offset=m)
    return out


def decode(raw):
    raw = bytes(raw)
    ver, verlz, flags, typesize = raw[0], raw[1], raw[2], raw[3]
    nbytes, blocksize, cbytes = struct.unpack("<III", raw[4:16])
    if typesize != 1:
        raise ValueError("only typesize 1 supported")
    codec = flags >> 5
    if flags & 0x2:                                            # stored uncompressed
        return np.frombuffer(raw, np.uint8, nbytes, offset=16).copy()
    nblocks = -(-nbytes // blocksize)
    bstarts = struct.unpack(f"<{nblocks}I", raw[16:16 + 4 * nblocks])
    out = np.empty(nbytes, np.uint8)
    z = _zstd()
    for i, s in enumerate(bstarts):
        bs = min(blocksize, nbytes - i * blocksize)
        clen = struct.unpack("<I", raw[s:s + 4])[0]
        body = raw[s + 4:s + 4 + clen]
        if clen == bs:
            blk = np.frombuffer(body, np.uint8)
        else:
            if codec != 4:
                raise ValueError(f"codec {codec} not supported")
            dst = ctypes.create_string_buffer(bs)
            r = z.ZSTD_decompress(dst, bs, body, clen)
            if z.ZSTD_isError(r) or r != bs:
                raise ValueError("zstd error")
            blk = np.frombuffer(dst.raw, np.uint8)
        if flags & 0x4:
            blk = _unbitshuffle(blk.tobytes(), bs)
        elif flags & 0x1:
            raise ValueError("byteshuffle not needed for typesize 1")
        out[i * blocksize:i * blocksize + bs] = blk
    return out
