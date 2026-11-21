"""lut_synth fast search kernels (numba): vectors of at most 128 quotient bits as (lo, hi) uint64 pairs,
branch-free row elimination. fsearch = one-level Gray walk (one-sided products), fsearch2 = two-level walk
(two-sided product, both fixed parts enumerated) with re-compaction after the outer elimination."""
import numpy as np
from numba import njit

ONE = np.uint64(1)
ZERO = np.uint64(0)


@njit(cache=True, inline="always")
def _member(glo, ghi, n, tlo, thi, wlo, whi):
    for m in range(n):
        wlo[m] = glo[m]
        whi[m] = ghi[m]
    for i in range(n):
        a = wlo[i]
        b = whi[i]
        if a != ZERO:
            piv = a & (~a + ONE)
            if tlo & piv:
                tlo ^= a
                thi ^= b
            for j in range(i + 1, n):
                msk = ZERO - np.uint64((wlo[j] & piv) != ZERO)
                wlo[j] ^= a & msk
                whi[j] ^= b & msk
        elif b != ZERO:
            piv = b & (~b + ONE)
            if thi & piv:
                thi ^= b
            for j in range(i + 1, n):
                msk = ZERO - np.uint64((whi[j] & piv) != ZERO)
                whi[j] ^= b & msk
    return tlo == ZERO and thi == ZERO


DB = np.uint64(0x03F79D71B4CB0A89)
S58 = np.uint64(58)
IDX64 = np.array([0, 1, 48, 2, 57, 49, 28, 3, 61, 58, 50, 42, 38, 29, 17, 4, 62, 55, 59, 36, 53, 51, 43, 22, 45, 39,
                  33, 30, 24, 18, 12, 5, 63, 47, 56, 27, 60, 41, 37, 16, 54, 35, 52, 21, 44, 32, 23, 11, 46, 26, 40,
                  15, 34, 20, 31, 10, 25, 14, 19, 9, 13, 8, 7, 6], dtype=np.uint64)


@njit(cache=True, inline="always")
def _member1(g, n, t, w):
    """Is t in span(g[0..n))? One 64-bit word per vector; shift-and-mask elimination (no data-dependent branch
    in the inner loop)."""
    for m in range(n):
        w[m] = g[m]
    for i in range(n):
        a = w[i]
        if a != ZERO:
            s = IDX64[((a & (~a + ONE)) * DB) >> S58]
            t ^= a & (ZERO - ((t >> s) & ONE))
            for j in range(i + 1, n):
                w[j] ^= a & (ZERO - ((w[j] >> s) & ONE))
    return t == ZERO


@njit(cache=True)
def _fsearch_nb(Glo, Ghi, blo, bhi, start, count):
    """g_m(x) = sum_k x_k G[k, m]; fixed rhs b. First Gray index in [start, start + count) with b in span(g(x))."""
    K, M = Glo.shape
    glo = np.zeros(M, dtype=np.uint64)
    ghi = np.zeros(M, dtype=np.uint64)
    gray = start ^ (start >> 1)
    for k in range(K):
        if (gray >> k) & 1:
            for m in range(M):
                glo[m] ^= Glo[k, m]
                ghi[m] ^= Ghi[k, m]
    wlo = np.empty(M, dtype=np.uint64)
    whi = np.empty(M, dtype=np.uint64)
    one_word = True
    for k in range(K):
        for m in range(M):
            if Ghi[k, m] != ZERO:
                one_word = False
    if bhi != ZERO:
        one_word = False
    for it in range(count):
        if one_word:
            hit = _member1(glo, M, blo, wlo)
        else:
            hit = _member(glo, ghi, M, blo, bhi, wlo, whi)
        if hit:
            return start + it
        nxt = start + it + 1
        k = 0
        while ((nxt >> k) & 1) == 0:
            k += 1
        if k >= K:
            return -1
        for m in range(M):
            glo[m] ^= Glo[k, m]
            ghi[m] ^= Ghi[k, m]
    return -1


@njit(cache=True)
def _fsearch2_nb(Dlo, Dhi, Qlo, Qhi, blo, bhi, nbits, use0, D0lo, D0hi, Q0lo, Q0hi, b0lo, b0hi,
                 y_start, n_outer, x_start, x_count):
    """Two-sided product (u0 + u')(v0 + v'), v0 = sum_j y_j p_j (outer Gray walk), u0 = sum_k x_k p_k (inner).
    D[k, c] = proj(delta_c p_k), Q[k, j] = proj(p_k p_j), b = proj(J); nbits = quotient bits (<= 128).
    Hit: b + proj(u0 v0) in span(delta_c v0) + span(delta_c u0), and b + proj(u0 v0) not already in span(delta_c v0).
    Returns (outer index, inner index) or (-1, -1)."""
    K, H = Dlo.shape
    flo = np.zeros(H, dtype=np.uint64)
    fhi = np.zeros(H, dtype=np.uint64)
    Blo = np.zeros(K, dtype=np.uint64)
    Bhi = np.zeros(K, dtype=np.uint64)
    gy = y_start ^ (y_start >> 1)
    for j in range(K):
        if (gy >> j) & 1:
            for c in range(H):
                flo[c] ^= Dlo[j, c]
                fhi[c] ^= Dhi[j, c]
            for k in range(K):
                Blo[k] ^= Qlo[k, j]
                Bhi[k] ^= Qhi[k, j]
    Elo = np.empty(H, dtype=np.uint64)
    Ehi = np.empty(H, dtype=np.uint64)
    Ew = np.empty(H, dtype=np.int64)        # 0: pivot in lo, 1: pivot in hi
    Ep = np.empty(H, dtype=np.uint64)       # pivot bit
    keep = np.empty(128, dtype=np.int64)
    Drl = np.empty((K, H), dtype=np.uint64)
    Drh = np.empty((K, H), dtype=np.uint64)
    Brl = np.empty(K, dtype=np.uint64)
    Brh = np.empty(K, dtype=np.uint64)
    glo = np.empty(H, dtype=np.uint64)
    ghi = np.empty(H, dtype=np.uint64)
    wlo = np.empty(H, dtype=np.uint64)
    whi = np.empty(H, dtype=np.uint64)
    v0lo = np.empty(2 * H, dtype=np.uint64)
    v0hi = np.empty(2 * H, dtype=np.uint64)
    w0lo = np.empty(2 * H, dtype=np.uint64)
    w0hi = np.empty(2 * H, dtype=np.uint64)
    for oi in range(n_outer):
        # echelon of the fixed generators delta_c v0
        ne = 0
        plo = ZERO
        phi = ZERO
        for c in range(H):
            a = flo[c]
            b = fhi[c]
            for i in range(ne):
                if Ew[i] == 0:
                    if a & Ep[i]:
                        a ^= Elo[i]
                        b ^= Ehi[i]
                else:
                    if b & Ep[i]:
                        a ^= Elo[i]
                        b ^= Ehi[i]
            if a != ZERO:
                Elo[ne] = a
                Ehi[ne] = b
                Ew[ne] = 0
                Ep[ne] = a & (~a + ONE)
                plo |= Ep[ne]
                ne += 1
            elif b != ZERO:
                Elo[ne] = a
                Ehi[ne] = b
                Ew[ne] = 1
                Ep[ne] = b & (~b + ONE)
                phi |= Ep[ne]
                ne += 1
        nk = 0
        for q in range(nbits):
            if q < 64:
                if not ((plo >> np.uint64(q)) & ONE):
                    keep[nk] = q
                    nk += 1
            else:
                if not ((phi >> np.uint64(q - 64)) & ONE):
                    keep[nk] = q
                    nk += 1
        # reduce and compact the tables
        brl = ZERO
        brh = ZERO
        for k in range(K + 1):
            ncol = H if k < K else 1
            for c in range(ncol + (1 if k < K else 0)):
                if k < K and c < H:
                    a = Dlo[k, c]
                    b = Dhi[k, c]
                elif k < K:
                    a = Blo[k]
                    b = Bhi[k]
                else:
                    a = blo
                    b = bhi
                for i in range(ne):
                    if Ew[i] == 0:
                        if a & Ep[i]:
                            a ^= Elo[i]
                            b ^= Ehi[i]
                    else:
                        if b & Ep[i]:
                            a ^= Elo[i]
                            b ^= Ehi[i]
                lo = ZERO
                hi = ZERO
                for t in range(nk):
                    q = keep[t]
                    if q < 64:
                        bit = (a >> np.uint64(q)) & ONE
                    else:
                        bit = (b >> np.uint64(q - 64)) & ONE
                    if t < 64:
                        lo |= bit << np.uint64(t)
                    else:
                        hi |= bit << np.uint64(t - 64)
                if k < K and c < H:
                    Drl[k, c] = lo
                    Drh[k, c] = hi
                elif k < K:
                    Brl[k] = lo
                    Brh[k] = hi
                else:
                    brl = lo
                    brh = hi
        # inner Gray walk
        xs = x_start if oi == 0 else 1
        for c in range(H):
            glo[c] = ZERO
            ghi[c] = ZERO
        rlo = brl
        rhi = brh
        gx = xs ^ (xs >> 1)
        for k in range(K):
            if (gx >> k) & 1:
                for c in range(H):
                    glo[c] ^= Drl[k, c]
                    ghi[c] ^= Drh[k, c]
                rlo ^= Brl[k]
                rhi ^= Brh[k]
        xi = xs
        one_word = nk <= 64
        while xi < x_count:
            if one_word:
                hit = rlo != ZERO and _member1(glo, H, rlo, wlo)
            else:
                hit = (rlo != ZERO or rhi != ZERO) and _member(glo, ghi, H, rlo, rhi, wlo, whi)
            if hit:
                trivial = False
                if use0:
                    # second test in the quotient WITHOUT the targets: a hit whose right-hand side is already
                    # in span(delta v0, delta u0) involves no target (tau = 0) and is skipped
                    yy = y_start + oi
                    gy2 = yy ^ (yy >> 1)
                    gx2 = xi ^ (xi >> 1)
                    for c in range(2 * H):
                        v0lo[c] = ZERO
                        v0hi[c] = ZERO
                    r0l = b0lo
                    r0h = b0hi
                    for j in range(K):
                        if (gy2 >> j) & 1:
                            for c in range(H):
                                v0lo[c] ^= D0lo[j, c]
                                v0hi[c] ^= D0hi[j, c]
                    for k in range(K):
                        if (gx2 >> k) & 1:
                            for c in range(H):
                                v0lo[H + c] ^= D0lo[k, c]
                                v0hi[H + c] ^= D0hi[k, c]
                            for j in range(K):
                                if (gy2 >> j) & 1:
                                    r0l ^= Q0lo[k, j]
                                    r0h ^= Q0hi[k, j]
                    trivial = _member(v0lo, v0hi, 2 * H, r0l, r0h, w0lo, w0hi)
                if not trivial:
                    return y_start + oi, xi
            xi += 1
            k = 0
            while ((xi >> k) & 1) == 0:
                k += 1
            if k >= K:
                break
            for c in range(H):
                glo[c] ^= Drl[k, c]
                ghi[c] ^= Drh[k, c]
            rlo ^= Brl[k]
            rhi ^= Brh[k]
        # next outer
        nxt = y_start + oi + 1
        j = 0
        while ((nxt >> j) & 1) == 0:
            j += 1
        if j >= K:
            return -1, -1
        for c in range(H):
            flo[c] ^= Dlo[j, c]
            fhi[c] ^= Dhi[j, c]
        for k in range(K):
            Blo[k] ^= Qlo[k, j]
            Bhi[k] ^= Qhi[k, j]
    return -1, -1


def split2(arr):
    """(..., W) uint64 with W in (1, 2) -> (lo, hi) contiguous arrays."""
    if arr.shape[-1] == 1:
        lo = np.ascontiguousarray(arr[..., 0])
        return lo, np.zeros_like(lo)
    assert arr.shape[-1] == 2, "quotient wider than 128 bits"
    return np.ascontiguousarray(arr[..., 0]), np.ascontiguousarray(arr[..., 1])


# ---- C kernels (built on the box from lut_synth_kern.c); numba versions are the fallback
import ctypes
import os

_LIB = None
_SO = os.path.expanduser("~/qcc/ws/lut_synth/kern.so")
if os.path.exists(_SO) and not os.environ.get("LUT_SYNTH_NUMBA"):
    _LIB = ctypes.CDLL(_SO)
    _U = np.ctypeslib.ndpointer(dtype=np.uint64, flags="C_CONTIGUOUS")
    _LIB.fsearch_c.restype = ctypes.c_longlong
    _LIB.fsearch_c.argtypes = [_U, _U, ctypes.c_int, ctypes.c_int, ctypes.c_uint64, ctypes.c_uint64,
                               ctypes.c_int, _U, _U, ctypes.c_uint64, ctypes.c_uint64,
                               ctypes.c_longlong, ctypes.c_longlong]
    _LIB.fsearch2_c.restype = None
    _LIB.fsearch2_c.argtypes = [_U, _U, _U, _U, ctypes.c_int, ctypes.c_int, ctypes.c_uint64, ctypes.c_uint64,
                                ctypes.c_int, ctypes.c_int, _U, _U, _U, _U, ctypes.c_uint64, ctypes.c_uint64,
                                ctypes.c_longlong, ctypes.c_longlong, ctypes.c_longlong, ctypes.c_longlong,
                                ctypes.POINTER(ctypes.c_longlong)]


def _c(a):
    return np.ascontiguousarray(a, dtype=np.uint64)


def fsearch(Glo, Ghi, blo, bhi, start, count, G0lo=None, G0hi=None, b0lo=0, b0hi=0):
    """G0 / b0: the same tables in the quotient without the targets (optional); with them the C kernel skips hits
    whose right-hand side needs no target. The numba fallback ignores them (the caller verifies every hit)."""
    if _LIB is None:
        return int(_fsearch_nb(Glo, Ghi, blo, bhi, start, count))
    K, M = Glo.shape
    use0 = G0lo is not None
    if not use0:
        G0lo, G0hi = Glo, Ghi
    return int(_LIB.fsearch_c(_c(Glo), _c(Ghi), K, M, int(blo), int(bhi), int(use0), _c(G0lo), _c(G0hi),
                              int(b0lo), int(b0hi), int(start), int(count)))


def fsearch2(Dlo, Dhi, Qlo, Qhi, blo, bhi, nbits, use0, D0lo, D0hi, Q0lo, Q0hi, b0lo, b0hi,
             y_start, n_outer, x_start, x_count):
    if _LIB is None:
        r = _fsearch2_nb(Dlo, Dhi, Qlo, Qhi, blo, bhi, nbits, use0, D0lo, D0hi, Q0lo, Q0hi, b0lo, b0hi,
                         y_start, n_outer, x_start, x_count)
        return int(r[0]), int(r[1])
    K, H = Dlo.shape
    out = (ctypes.c_longlong * 2)()
    _LIB.fsearch2_c(_c(Dlo), _c(Dhi), _c(Qlo), _c(Qhi), K, H, int(blo), int(bhi), int(nbits), int(bool(use0)),
                    _c(D0lo), _c(D0hi), _c(Q0lo), _c(Q0hi), int(b0lo), int(b0hi), int(y_start), int(n_outer),
                    int(x_start), int(x_count), out)
    return int(out[0]), int(out[1])


KERNEL = "C" if _LIB is not None else "numba"
