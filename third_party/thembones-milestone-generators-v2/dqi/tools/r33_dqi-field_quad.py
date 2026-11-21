#!/usr/bin/env python3
"""r33_dqi-field_quad.py -- exact AND count of small field products when an AND may multiply two MIXED
linear forms (seat DQI-FIELD, design round R33).

Map: one scalar x in GF(2^k) times n words y_1..y_n in GF(2^k) -> the n products (k*n output bits).
Program: L AND gates g_i = l_i(v) * l'_i(v), v = (x bits, y bits), and an output matrix W (k*n by L):
output_c = sum_i W[c][i] g_i + (linear terms, free).  The quadratic part of g_i is the alternating form
l_i ^ l'_i; the program is correct iff for every output c the alternating forms add up to the target:
    cross block (x_a, y_b): coefficient of x_a y_b in output c
    x-x block and y-y block: zero ("garbage" must cancel), unless the rule frees a block.
Rules:
    bil     l on x only, l' on y only (the bilinear rank; sanity: GF(4) 3, GF(8) 6, GF(16) 9)
    pure    mixed forms, both garbage blocks must cancel (the true quadratic AND count of the map)
    xfree   mixed forms, x-x garbage free (the scalar's own quadratic forms come from a cache)
    yfree   mixed forms, y-y garbage free, x-x must cancel: ONE word x times n CACHED scalars y_1..y_n
            (the word's own garbage is paid once inside the program and shared by the n products)
    yfree1  as yfree, but only the y-y pairs inside one scalar are free (no cross-scalar cache)
    yfree2  as yfree1 with one pair of each scalar forced to cancel (a smaller cache)
    xyfree  mixed forms, both garbage blocks free (both operands carry a cache)
    A rule with the prefix F (Fbil, Fpure, Fxfree, ...) takes the fused chain step of GF(2^K) over GF(2^(K/2)),
    (H, D) -> (H D, sqrt(Tr H) D), instead of the product (x = H, y = D; N is ignored; use mw = 0 there).
    A rule with the prefix P (Pbil, Ppure, Pxfree) takes the fold (a + b z) A(z): x = the head's two words,
    y = the N words of A, N + 1 output words (evaluation costs (N + 1) x rank).
Usage: r33_dqi-field_quad.py K N L RULE [time_limit_s=600] [seed=0] [lex=1] [mw=0] [fixw=0]
  mw > 0 adds the implied minimum-weight clauses for combinations of up to mw output rows (same answer, faster).
Prints one JSON line: {"k","n","L","rule","result": "SAT"|"UNSAT"|"UNKNOWN","seconds", "program" if SAT (verified)}.
"""
import sys
import json
import time
import itertools

MODS = {2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011, 7: 0b10000011, 8: 0b100011011,
        12: 4179}   # 12: the competition modulus x^12 + x^6 + x^4 + x + 1


def gf_mul(a, b, k):
    m = MODS[k]
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> k:
            a ^= m
    return r


def targets(k, n):
    """T[c] = set of (a, b) variable pairs (a < b) whose product appears in output c. Variables: x 0..k-1,
    word j bit b -> k + j*k + b. Outputs: c = j*k + bit."""
    N = k + k * n
    T = [set() for _ in range(k * n)]
    for a in range(k):
        for b in range(k):
            p = gf_mul(1 << a, 1 << b, k)
            for j in range(n):
                for bit in range(k):
                    if (p >> bit) & 1:
                        T[j * k + bit].add((a, k + j * k + b))
    return N, T


def gf_pow(a, e, k):
    r = 1
    while e:
        if e & 1:
            r = gf_mul(r, a, k)
        a = gf_mul(a, a, k)
        e >>= 1
    return r


def targets_fused(k):
    """The c = 0 fused chain step in the field L = GF(2^k) over K = GF(2^(k/2)) (the lens log's analogue):
    (H, D) -> (H D, sqrt(Tr_{L/K} H) D).  Variables: H bits 0..k-1, D bits k..2k-1; outputs 2k bits."""
    assert k % 2 == 0
    N = 2 * k
    T = [set() for _ in range(2 * k)]
    for a in range(k):
        h = 1 << a
        phi = gf_pow(h ^ gf_pow(h, 1 << (k // 2), k), 1 << (k - 1), k)
        for b in range(k):
            p1 = gf_mul(h, 1 << b, k)
            p2 = gf_mul(phi, 1 << b, k)
            for bit in range(k):
                if (p1 >> bit) & 1:
                    T[bit].add((a, k + b))
                if (p2 >> bit) & 1:
                    T[k + bit].add((a, k + b))
    return N, T


def targets_fold(k, n):
    """The fold of the synthesis pair in the field GF(2^k): head (a, b) times A_0..A_{n-1}, the product
    (a + b z) A(z) with its n + 1 coefficient words c_j = a A_j + b A_{j-1}.  Variables: a bits 0..k-1,
    b bits k..2k-1, word j bit t at 2k + j k + t.  Evaluation at n + 1 points costs (n + 1) x rank(k)."""
    N = 2 * k + k * n
    T = [set() for _ in range(k * (n + 1))]
    for s in range(k):
        for t in range(k):
            p = gf_mul(1 << s, 1 << t, k)
            for bit in range(k):
                if (p >> bit) & 1:
                    for j in range(n):
                        T[j * k + bit].add((s, 2 * k + j * k + t))            # a A_j -> c_j
                        T[(j + 1) * k + bit].add((k + s, 2 * k + j * k + t))  # b A_j -> c_(j+1)
    return N, T


FUSED = False
FOLD = False


def build(k, n, L, rule, lex=True, seed=0, mw=0, fixw=0):
    from pycryptosat import Solver
    N, T = targets_fused(k) if FUSED else (targets_fold(k, n) if FOLD else targets(k, n))
    nout = len(T)
    kx = 2 * k if FOLD else k          # number of x variables (the head of a fold has two words)
    s = Solver(threads=1)
    cnt = [0]

    def new():
        cnt[0] += 1
        return cnt[0]

    xs = list(range(k))
    ys = list(range(k, N))
    lv = [[None] * N for _ in range(L)]
    lpv = [[None] * N for _ in range(L)]
    for i in range(L):
        for a in range(N):
            if rule == 'bil':
                lv[i][a] = new() if a < kx else 0
                lpv[i][a] = new() if a >= kx else 0
            else:
                lv[i][a] = new()
                lpv[i][a] = new()
    W = [[new() for _ in range(L)] for _ in range(nout)]

    def AND(p, q):
        if p == 0 or q == 0:
            return 0
        r = new()
        s.add_clause([-r, p])
        s.add_clause([-r, q])
        s.add_clause([r, -p, -q])
        return r

    def XOR2(p, q):
        if p == 0:
            return q
        if q == 0:
            return p
        r = new()
        s.add_xor_clause([r, p, q], False)
        return r

    def constrained(a, b):
        ax, bx = a < kx, b < kx
        if ax != bx:
            return True
        if rule in ('bil', 'pure'):
            return True
        if rule == 'xfree':
            return not ax  # y-y must cancel, x-x free
        if rule == 'yfree':
            return ax      # x-x must cancel, every y-y pair (also across words) free
        if rule == 'yfree1':
            # x-x must cancel; y-y free only inside one word (each scalar brings its own cache, no cross terms)
            return ax or ((a - k) // k != (b - k) // k)
        if rule == 'yfree2':
            # as yfree1, and the pair of the two highest bits of each scalar must cancel too: the cache of a
            # scalar then has at most C(k,2) - 1 forms (k = 3: two forms per scalar)
            if ax or ((a - k) // k != (b - k) // k):
                return True
            return (a - k) % k == k - 2 and (b - k) % k == k - 1
        if rule == 'xyfree':
            return False
        raise ValueError(rule)

    pairs = [(a, b) for a in range(N) for b in range(a + 1, N) if constrained(a, b)]
    d = {}
    for i in range(L):
        for (a, b) in pairs:
            d[(i, a, b)] = XOR2(AND(lv[i][a], lpv[i][b]), AND(lv[i][b], lpv[i][a]))
    for c in range(nout):
        for (a, b) in pairs:
            terms = []
            for i in range(L):
                t = AND(W[c][i], d[(i, a, b)])
                if t:
                    terms.append(t)
            rhs = (a, b) in T[c]
            if terms:
                s.add_xor_clause(terms, rhs)
            elif rhs:
                s.add_clause([])  # impossible
    if fixw:
        # Single product, L = 2^(k-1) gates, minimum weight k: the output code is the [2^(k-1), k] code whose
        # columns are an affine hyperplane {Tr(kappa g) = 1}; x -> kappa x is a symmetry of the product, so the
        # columns are fixed to the elements of absolute trace 1, in increasing order (valid for k = 4, L = 8).
        assert n == 1 and L == (1 << (k - 1))
        def tr(g):
            t, v = 0, g
            for _ in range(k):
                t ^= v
                v = gf_mul(v, v, k)
            return t
        cols = [g for g in range(1, 1 << k) if tr(g) == 1]
        assert len(cols) == L, (len(cols), L)
        for i, g in enumerate(cols):
            for c in range(k):
                s.add_clause([W[c][i] if (g >> c) & 1 else -W[c][i]])
        lex = False
    if mw:
        # implied constraint (coding bound): every nonzero combination of outputs is carried by at least
        # dmin gates, dmin = k (one garbage block must cancel) or ceil(k/2) (both free).  Only combinations
        # of at most mw output rows are written out.
        dmin = (k + 1) // 2 if rule == 'xyfree' else k
        if FUSED:
            dmin = (k - 1 + 1) // 2 if rule == 'xyfree' else k - 1
        size = L - dmin + 1
        if size >= 1:
            for r in range(1, mw + 1):
                for rows in itertools.combinations(range(nout), r):
                    if r == 1:
                        z = W[rows[0]]
                    else:
                        z = []
                        for i in range(L):
                            v = new()
                            s.add_xor_clause([v] + [W[c][i] for c in rows], False)
                            z.append(v)
                    for S in itertools.combinations(range(L), size):
                        s.add_clause([z[i] for i in S])
    if lex and L > 1:
        # W columns in non-increasing lexicographic order (gate order is a symmetry)
        for i in range(L - 1):
            eq_prev = None  # literal: columns i and i+1 agree on rows < c
            for c in range(nout):
                p, q = W[c][i], W[c][i + 1]
                # if equal so far then p >= q
                if eq_prev is None:
                    s.add_clause([p, -q])
                else:
                    s.add_clause([-eq_prev, p, -q])
                if c < nout - 1:
                    e = new()  # e <-> eq_prev and (p == q)
                    x = new()
                    s.add_xor_clause([x, p, q], False)  # x = p xor q
                    if eq_prev is None:
                        s.add_clause([-e, -x])
                        s.add_clause([e, x])
                    else:
                        s.add_clause([-e, eq_prev])
                        s.add_clause([-e, -x])
                        s.add_clause([e, -eq_prev, x])
                    eq_prev = e
    return s, N, T, lv, lpv, W, cnt[0]


def verify(k, n, L, rule, N, T, lvals, lpvals, Wvals):
    """Check the alternating-form identity directly from the decoded program."""
    nout = len(T)
    for c in range(nout):
        acc = {}
        for i in range(L):
            if not Wvals[c][i]:
                continue
            for a in range(N):
                for b in range(a + 1, N):
                    v = (lvals[i][a] & lpvals[i][b]) ^ (lvals[i][b] & lpvals[i][a])
                    if v:
                        acc[(a, b)] = acc.get((a, b), 0) ^ 1
        got = {p for p, v in acc.items() if v}
        kx = 2 * k if FOLD else k
        for (a, b) in itertools.combinations(range(N), 2):
            ax, bx = a < kx, b < kx
            if ax != bx:
                need = True
            elif rule in ('bil', 'pure'):
                need = True
            elif rule == 'xfree':
                need = not ax
            elif rule == 'yfree':
                need = ax
            elif rule == 'yfree1':
                need = ax or ((a - k) // k != (b - k) // k)
            elif rule == 'yfree2':
                need = ax or ((a - k) // k != (b - k) // k) or ((a - k) % k == k - 2 and (b - k) % k == k - 1)
            else:
                need = False
            if need and (((a, b) in got) != ((a, b) in T[c])):
                return False
    return True


def main():
    global FUSED
    k, n, L = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rule = sys.argv[4]
    global FOLD
    if rule.startswith('F'):      # Fbil / Fpure / Fxfree / Fyfree / Fxyfree: the fused-step analogue (N = 0 by convention)
        FUSED = True
        rule = rule[1:]
    elif rule.startswith('P'):    # Pbil / Ppure / Pxfree: the fold (a + b z) A(z), A of N words
        FOLD = True
        rule = rule[1:]
    tl = float(sys.argv[5]) if len(sys.argv) > 5 else 600.0
    seed = int(sys.argv[6]) if len(sys.argv) > 6 else 0
    lex = (sys.argv[7] != '0') if len(sys.argv) > 7 else True
    mw = int(sys.argv[8]) if len(sys.argv) > 8 else 0
    fixw = int(sys.argv[9]) if len(sys.argv) > 9 else 0
    t0 = time.time()
    s, N, T, lv, lpv, W, nv = build(k, n, L, rule, lex, seed, mw, fixw)
    sat, sol = s.solve(time_limit=tl) if hasattr(s, 'solve') else (None, None)
    dt = time.time() - t0
    out = {'k': k, 'n': n, 'L': L, 'rule': ('F' if FUSED else 'P' if FOLD else '') + rule, 'vars': nv, 'seconds': round(dt, 1), 'lex': lex, 'mw': mw, 'fixw': fixw}
    if sat is True:
        lvals = [[int(bool(sol[v])) if v else 0 for v in row] for row in lv]
        lpvals = [[int(bool(sol[v])) if v else 0 for v in row] for row in lpv]
        Wvals = [[int(bool(sol[v])) for v in row] for row in W]
        ok = verify(k, n, L, rule, N, T, lvals, lpvals, Wvals)
        out['result'] = 'SAT' if ok else 'SAT-BUT-VERIFY-FAILED'
        out['program'] = {'l': [''.join(map(str, r)) for r in lvals], 'lp': [''.join(map(str, r)) for r in lpvals],
                          'W': [''.join(map(str, r)) for r in Wvals]}
    elif sat is False:
        out['result'] = 'UNSAT'
    else:
        out['result'] = 'UNKNOWN'
    print(json.dumps(out), flush=True)


if __name__ == '__main__':
    main()
