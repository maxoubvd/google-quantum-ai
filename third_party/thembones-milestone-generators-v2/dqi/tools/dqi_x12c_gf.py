"""dqi_x12c_gf.py -- GF(2^12) arithmetic for the classical models of workstream dqi_x12c.

Field: GF(2^12) with modulus 0x1053 (the competition's a^12 + a^6 + a^4 + a + 1), primitive element 2.
Vector functions work elementwise on numpy int64 arrays (one entry per shot); scalar functions take ints.
Polynomials are lists (scalar) or 2-D arrays (shots, coefficients), lowest degree first.
"""
import numpy as np

MOD = 0x1053
ORD = 4095
I64 = np.int64

EXP = np.zeros(2 * ORD + 2, dtype=I64)
LOG = np.zeros(4096, dtype=I64)
_x = 1
for _i in range(ORD):
    EXP[_i] = _x
    LOG[_x] = _i
    _x <<= 1
    if _x & 4096:
        _x ^= MOD
assert _x == 1
EXP[ORD:2 * ORD] = EXP[:ORD]


# ---- vector arithmetic --------------------------------------------------------------------------------------
def mul(a, b):
    a = np.asarray(a, dtype=I64)
    b = np.asarray(b, dtype=I64)
    return np.where((a == 0) | (b == 0), 0, EXP[LOG[a] + LOG[b]])


def inv0(a):
    """1/a, and 0 for a = 0 (the kit's convention)."""
    a = np.asarray(a, dtype=I64)
    return np.where(a == 0, 0, EXP[(ORD - LOG[a]) % ORD])


def sq(a):
    a = np.asarray(a, dtype=I64)
    return np.where(a == 0, 0, EXP[(2 * LOG[a]) % ORD])


# ---- scalar arithmetic --------------------------------------------------------------------------------------
def smul(a, b):
    return int(EXP[LOG[a] + LOG[b]]) if a and b else 0


def sinv(a):
    assert a
    return int(EXP[(ORD - LOG[a]) % ORD])


def spow(a, n):
    if n == 0:
        return 1
    return int(EXP[(LOG[a] * n) % ORD]) if a else 0


# ---- scalar polynomials -------------------------------------------------------------------------------------
def pmul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                c[i + j] ^= smul(x, y)
    return c


def pdiv_linear(p, x):
    """Exact quotient p / (z + x) (asserts that x is a root)."""
    out = [0] * (len(p) - 1)
    carry = 0
    for i in reversed(range(1, len(p))):
        carry = p[i] ^ smul(carry, x)
        out[i - 1] = carry
    assert (p[0] ^ smul(carry, x)) == 0, 'not a root'
    return out


# ---- vector polynomials (shape (shots, n), lowest degree first) -------------------------------------------
def veval(P, x):
    """P(x) per shot; x is an int or an array of length shots."""
    r = np.zeros(P.shape[0], dtype=I64)
    for j in reversed(range(P.shape[1])):
        r = mul(r, x) ^ P[:, j]
    return r


def vderiv(P):
    """Formal derivative (characteristic 2: only odd coefficients survive)."""
    D = np.zeros_like(P)
    D[:, 0:-1:2] = P[:, 1::2]
    return D


def vmul_scalar_poly(s, P):
    """s * P with s an array (shots,) and P (shots, n)."""
    return mul(s[:, None], P)


def vpad(P, n):
    if P.shape[1] >= n:
        return P
    return np.concatenate([P, np.zeros((P.shape[0], n - P.shape[1]), dtype=I64)], axis=1)


def vdegree(P):
    nz = P != 0
    return np.where(nz.any(axis=1), P.shape[1] - 1 - nz[:, ::-1].argmax(axis=1), -1)
