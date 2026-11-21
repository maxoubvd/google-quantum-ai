#!/usr/bin/env python3
"""r48_dqi_ret2_forms.py -- seat R48-DQI-RET2: D1's CRT body of un-synthesis with ONE point at infinity.

D1's a1 (tools/dqi_pub_a1.py / dqi_pub_a1model.py, read only) writes the body's quotient F = (l S_c - U) / q through
its values at t_j = 0..n-1 (n = len(T)): a free premap of the dividend U before the division plus n rank-one
products after it. Here the last point is infinity: F = sum_{j < n-1} F(t_j) L_j(x) + F_top (x^(n-1) - sum_j
t_j^(n-1) L_j(x)), with L_j the Lagrange basis of the n - 1 finite points and F_top = l_1 S_top (the leading
coefficient of l S_c; U has degree < 2 and does not reach it). Same number of products (n). What changes: the premap
is a combination of polynomials of degree <= n - 2, so it never writes the LAST target word (the leading
coefficient in the q-adic digit basis = the leader word a2 / beta of the frame). The late-flag division needs that:
on a terminated fold the raw dividend is leader x constant, and a premap that wrote the leader would rescale it.

Functions take D1's model module (pure Python GF(2^12)) as `M`. Self-test:
    python -I r48_dqi_ret2_forms.py <dir holding dqi_pub_a1model.py> [trials=300] [seed=1]
checks (1) premap + products = the kit's floor map on random inputs (D1's own check, ported), (2) pre[last] = 0,
(3) the terminated case: with U_raw = lead (c0, c1), l = 0: premap + completion constants corr = kappa + co_T give
T + lead co_T, corr[last] = 0, for the 15 flagged sites of the line.
"""
import importlib.util
import os
import random
import sys


def load_model(d):
    spec = importlib.util.spec_from_file_location('r48_a1model', os.path.join(d, 'dqi_pub_a1model.py'))
    M = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(M)
    return M


def crt_forms_inf(M, m, n, k, pts=None):
    """Forms for a body with source length m and target length n = m - 1, last point at infinity.
    Returns (pre, updates): pre[i] = (coefficient of U[0], of U[1]) in T[i]; updates = [(ac, bc, tc)].
    pts: the n - 1 finite points (distinct field elements; default 0..n-2)."""
    assert n == m - 1 and n >= 2
    mul, inv, power = M.mul, M.inv, M.power
    Qs, Qt = M.qadic_matrix(m, k), M.qadic_matrix(n, k)
    Qs_inv = M.matinv(Qs)
    pts = list(range(n - 1)) if pts is None else list(pts)
    assert len(pts) == n - 1 and len(set(pts)) == n - 1
    V = [[power(t, j) for j in range(n - 1)] for t in pts]
    W = M.matinv(V)                                  # coefficients of L_j (degree <= n - 2): column j of W
    q = [k, 1, 1]
    tcs = []
    for j, t in enumerate(pts):
        qt = M.peval(q, t)
        Lj = [W[i][j] for i in range(n - 1)] + [0]
        tcs.append([mul(c, inv(qt)) for c in M.matvec(Qt, Lj)])
    top = [0] * n                                    # x^(n-1) - sum_j t_j^(n-1) L_j
    top[n - 1] = 1
    for j, t in enumerate(pts):
        f = power(t, n - 1)
        for i in range(n - 1):
            top[i] ^= mul(f, W[i][j])
    tc_inf = M.matvec(Qt, top)
    pre = [[0, 0] for _ in range(n)]
    for j, t in enumerate(pts):
        for i in range(n):
            pre[i][0] ^= tcs[j][i]
            pre[i][1] ^= mul(tcs[j][i], t)
    updates = []
    for j, t in enumerate(pts):
        row = [power(t, i) for i in range(m)]
        bc = [0] * m
        for i in range(m):
            for c in range(m):
                bc[c] ^= mul(row[i], Qs_inv[i][c])
        updates.append(([1, t], bc, tcs[j]))
    bc = [Qs_inv[m - 1][c] for c in range(m)]        # the top coefficient of S_c as a form of the digit words
    updates.append(([0, 1], bc, tc_inf))
    return pre, updates


def constant_qadic(M, k, degree, n):
    """The kit's padded_return.constant_qadic: q-adic digits of (k + x) x^(degree - 2)."""
    mul = M.mul
    p = [0] * (degree - 2) + [k, 1]
    out = []
    while len(out) < n:
        q = [0] * max(0, len(p) - 2)
        for i in range(len(p) - 1, 1, -1):
            a = p[i]
            q[i - 2] = a
            p[i] ^= a
            p[i - 1] ^= a
            p[i - 2] ^= mul(k, a)
        out.extend((p + [0, 0])[:2])
        p = q
    assert not any(p)
    return out[:n]


def kpoints(nt, k):
    """Finite points {0, k} for a 3-word target: both are roots of the padding polynomial (k + x) x^(d-2), so the
    premap of a terminated fold's raw dividend IS the completion and no completion gate is left (corr = 0)."""
    return [0, k] if nt == 3 else None


def completion(M, k, degree, ntarget, nt, pts=None):
    """(c0, c1, corr): the dividend's padding digit and the T-constants of the late-flag completion
    (kit constant + what the infinity premap added on the raw dividend); corr[-1] must be 0."""
    mul = M.mul
    co = constant_qadic(M, k, degree, ntarget)
    assert ntarget == nt + 2 and co[-1] == 0
    pre, _ = crt_forms_inf(M, nt + 1, nt, k, pts)
    corr = [co[2 + i] ^ mul(pre[i][0], co[0]) ^ mul(pre[i][1], co[1]) for i in range(nt)]
    return co[0], co[1], corr


def check(M, trials, seed):
    mul = M.mul
    rnd = random.Random(seed)
    bad = lastbad = 0
    for trial in range(trials):
        n = rnd.choice([3, 5, 7, 9, 15, 33, 35]) - (trial % 2)
        m = n + 1
        c = 2 * rnd.randrange(17)
        k = 2048 ^ mul(c, c) ^ c
        S = [rnd.randrange(4096) for _ in range(m)]
        T = [rnd.randrange(4096) for _ in range(n)]
        l = [rnd.randrange(4096), rnd.randrange(1, 4096)]
        Qs, Qt = M.qadic_matrix(m, k), M.qadic_matrix(n, k)
        Sc = M.matvec(M.matinv(Qs), S)
        fl, rem = M.divmodp(M.pmul(l, Sc), [k, 1, 1])
        fl = (fl + [0] * n)[:n]
        kit = [a ^ b for a, b in zip(T, M.matvec(Qt, fl))]
        U = (rem + [0, 0])[:2]
        pre, ups = crt_forms_inf(M, m, n, k)
        lastbad += pre[-1] != [0, 0]
        mine = list(T)
        for i in range(n):
            mine[i] ^= mul(pre[i][0], U[0]) ^ mul(pre[i][1], U[1])
        for ac, bc, tc in ups:
            lv = mul(ac[0], l[0]) ^ mul(ac[1], l[1])
            sv = 0
            for a, s in zip(bc, S):
                sv ^= mul(a, s)
            pr = mul(lv, sv)
            for i in range(n):
                mine[i] ^= mul(tc[i], pr)
        bad += mine != kit
    return bad, lastbad


def sites(M):
    """The 15 flagged divisions of the line: (head, frame, which, n, k, degree, ntarget, nt)."""
    mul = M.mul
    out = []
    for j in range(17):
        c = 2 * j
        k = 2048 ^ mul(c, c) ^ c
        n = 35 - 2 * j
        for which, degree, ntarget, nt in (('first', n + 1, n + 2, n), ('second', n, n + 1, n - 1)):
            head = 2 * j + (which == 'second')
            if head >= 19:
                out.append((head, j, which, n, k, degree, ntarget, nt))
    return out


def check_terminated(M, seed):
    """Terminated fold on the wires: U_raw = lead (c0, c1), head 0. Late flag: premap(U_raw), no product
    (l = 0), completion T += lead corr. The line: T += lead co_T. Both must agree and leave the leader."""
    mul = M.mul
    rnd = random.Random(seed)
    bad = 0
    for head, j, which, n, k, degree, ntarget, nt in sites(M):
        c0, c1, corr = completion(M, k, degree, ntarget, nt)
        co = constant_qadic(M, k, degree, ntarget)
        pre, _ = crt_forms_inf(M, nt + 1, nt, k)
        ok = corr[-1] == 0 and pre[-1] == [0, 0]
        for _ in range(20):
            lead = rnd.randrange(1, 4096)
            T = [rnd.randrange(4096) for _ in range(nt - 1)] + [lead]      # the last target word IS the leader
            U = [mul(lead, c0), mul(lead, c1)]
            mine = [T[i] ^ mul(pre[i][0], U[0]) ^ mul(pre[i][1], U[1]) for i in range(nt)]
            mine = [mine[i] ^ mul(lead, corr[i]) for i in range(nt)]
            line = [T[i] ^ mul(lead, co[2 + i]) for i in range(nt)]
            ok = ok and mine == line and mine[-1] == lead
        print('head %2d frame %2d %-6s nt %2d k %4d  (c0, c1) (%4d, %4d)  corr nonzero %2d of %2d  corr_last %d  %s'
              % (head, j, which, nt, k, c0, c1, sum(1 for v in corr if v), nt, corr[-1], 'ok' if ok else 'BAD'))
        bad += not ok
    return bad


def check_kpoints(M, seed):
    """Target length 3 with finite points {0, k}: floor map on random inputs, and corr = 0 at head 32."""
    mul = M.mul
    rnd = random.Random(seed)
    bad = 0
    for trial in range(200):
        c = 2 * rnd.randrange(17)
        k = 2048 ^ mul(c, c) ^ c
        n, m = 3, 4
        S = [rnd.randrange(4096) for _ in range(m)]
        T = [rnd.randrange(4096) for _ in range(n)]
        l = [rnd.randrange(4096), rnd.randrange(4096)]
        Qs, Qt = M.qadic_matrix(m, k), M.qadic_matrix(n, k)
        Sc = M.matvec(M.matinv(Qs), S)
        fl, rem = M.divmodp(M.pmul(l, Sc), [k, 1, 1])
        fl = (fl + [0] * n)[:n]
        kit = [a ^ b for a, b in zip(T, M.matvec(Qt, fl))]
        U = (rem + [0, 0])[:2]
        pre, ups = crt_forms_inf(M, m, n, k, kpoints(n, k))
        mine = [T[i] ^ mul(pre[i][0], U[0]) ^ mul(pre[i][1], U[1]) for i in range(n)]
        for ac, bc, tc in ups:
            lv = mul(ac[0], l[0]) ^ mul(ac[1], l[1])
            sv = 0
            for a, s in zip(bc, S):
                sv ^= mul(a, s)
            pr = mul(lv, sv)
            mine = [mine[i] ^ mul(tc[i], pr) for i in range(n)]
        bad += mine != kit or pre[-1] != [0, 0]
    out = []
    for head, j, which, n, k, degree, ntarget, nt in sites(M):
        if nt == 3:
            out.append((head, completion(M, k, degree, ntarget, nt, kpoints(nt, k))[2]))
    return bad, out


if __name__ == '__main__':
    Mod = load_model(sys.argv[1])
    trials = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    b1, b2 = check(Mod, trials, seed)
    print('floor map: mismatches %d of %d; premap writes the last word in %d' % (b1, trials, b2))
    print('terminated folds: bad sites', check_terminated(Mod, seed))
    b3, corr = check_kpoints(Mod, seed)
    print('points {0, k} at target length 3: floor map mismatches %d of 200; (head, corr): %s' % (b3, corr))
