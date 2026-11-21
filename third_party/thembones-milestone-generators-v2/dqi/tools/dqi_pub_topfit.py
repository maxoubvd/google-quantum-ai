#!/usr/bin/env python3
"""dqi_pub_topfit.py -- builder D1: fit the Laurent top X~_34 of our converter from a dqi_pub_convtest dump.

From the dump (A, B, C = ref output, K = a3 output; normal domain) in the square-root domain:
  X_0..X_33 by back-substitution from C'_34..C'_1 = H(B')X; true X~_34 = C'_0 + sum_{j<=33} B'_{j+1} X_j;
  the circuit's X^_34 = X~_34 + sqrt(e), e = K_34 - (M C mod B)_34 (the error is e M);
  top_raw = X^_34 + sum_{j in J} B'_{j+1} X_j (J of tools/dqi_pub_a3.py).
Fits X~_34 - top_raw = sum c_ab B'_a X_b over GF(2^12) (a = 0..35, b = 0..33) and prints the nonzero c_ab.
Usage: python dqi_pub_topfit.py <dump prefix> <count>
"""
import sys

import numpy as np

_MOD = 4179
_EXP = np.zeros(8190, dtype=np.int64)
_LOG = np.zeros(4096, dtype=np.int64)
_x = 1
for _i in range(4095):
    _EXP[_i] = _x
    _LOG[_x] = _i
    _x <<= 1
    if _x & 4096:
        _x ^= _MOD
_EXP[4095:] = _EXP[:4095]


def mul(a, b):
    a = np.asarray(a, dtype=np.int64)
    b = np.asarray(b, dtype=np.int64)
    return np.where((a == 0) | (b == 0), 0, _EXP[(_LOG[a] + _LOG[b]) % 4095])


class FT:
    @staticmethod
    def inv(a):
        return int(_EXP[(4095 - _LOG[a]) % 4095])
J = (0, 2, 4, 6, 8, 10, 12, 14, 17, 18, 20, 22, 24, 26, 28, 30, 32)


def gsqrt(v):
    out = v.copy()
    for _ in range(11):
        out = mul(out, out)
    return out


def main():
    prefix, count = sys.argv[1], int(sys.argv[2])
    rows_F, rows_t = [], []
    for i in range(count):
        z = np.load('%s_%d.npz' % (prefix, i))
        B, C, K = z['B'], z['C'], z['K']
        n = B.shape[0]
        Bs, Cs = gsqrt(B), gsqrt(C)
        Bf = np.concatenate([Bs, np.ones((n, 1), dtype=np.int64)], axis=1)
        X = np.zeros((n, 34), dtype=np.int64)
        for j in range(34):                      # C'_{34-j} = sum_{i<=j} B'_{35-j+i} X_i
            k = 34 - j
            acc = Cs[:, k].copy()
            for i in range(j):
                acc ^= mul(Bf[:, k + i + 1], X[:, i])
            X[:, j] = acc                        # coefficient B'_35 = 1
        row0 = np.zeros(n, dtype=np.int64)
        for j in range(34):
            row0 ^= mul(Bf[:, j + 1], X[:, j])
        xt = Cs[:, 0] ^ row0
        # error scalar from the a3 output
        Mn = np.array(z['M']) if 'M' in z.files else None
        e = z['E'] if 'E' in z.files else None
        if e is None:
            raise SystemExit('dump lacks E (error scalar); rerun convtest with the E field')
        xhat = xt ^ gsqrt(e)
        top_raw = xhat.copy()
        for j in J:
            top_raw ^= mul(Bf[:, j + 1], X[:, j])
        rows_t.append(xt ^ top_raw)
        rows_F.append(np.stack([mul(Bf[:, a], X[:, b]) for a in range(36) for b in range(34)], axis=1))
    F = np.concatenate(rows_F)
    t = np.concatenate(rows_t)
    np.savez('/tmp/d1_topfit.npz', F=F, t=t)
    cols = [(a, b) for a in range(36) for b in range(34)]
    hold = 64
    nf = len(t) - hold
    M = np.concatenate([F[:nf], t[:nf, None]], axis=1).copy()
    piv, r = [], 0
    for c in range(len(cols)):
        nz = np.flatnonzero(M[r:, c])
        if len(nz) == 0:
            continue
        p = r + nz[0]
        M[[r, p]] = M[[p, r]]
        M[r] = mul(M[r], FT.inv(int(M[r, c])))
        f = M[:, c].copy()
        f[r] = 0
        nzr = np.flatnonzero(f)
        if len(nzr):
            M[nzr] ^= mul(f[nzr, None], M[r][None, :])
        piv.append(c)
        r += 1
        if r == M.shape[0]:
            break
    sol = np.zeros(len(cols), dtype=np.int64)
    for i, c in enumerate(piv):
        sol[c] = M[i, len(cols)]
    incons = int(np.count_nonzero(M[len(piv):, len(cols)]))
    pred = np.zeros(hold, dtype=np.int64)
    for k in np.flatnonzero(sol):
        pred ^= mul(F[nf:, k], sol[k])
    print('shots', len(t), 'rank', len(piv), 'inconsistent', incons, 'held-out mismatches',
          int(np.count_nonzero(pred != t[nf:])), 'of', hold)
    for k in np.flatnonzero(sol):
        print('  B%d * X%d  c=%d' % (cols[k][0], cols[k][1], sol[k]))


if __name__ == '__main__':
    main()
