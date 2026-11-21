#!/usr/bin/env python3
"""dqi_pub_qmap.py -- builder D1 (holder 4): fewest ANDs for a given QUADRATIC map over GF(2) (exact, SAT).

Input: a JSON file {"anf": [[monomial, ...] per output]} (monomials as bitmasks over n variables, degree <= 2), e.g.
from tools/dqi_pub_inv27model.py. Program: L gates g_i = (l_i . x) (l'_i . x), output_c = sum_i W[c][i] g_i +
(linear terms, free). Correct iff for every output c and every pair u < v the alternating forms add up to the
target's quadratic part (diagonal products only add linear terms).
Usage: dqi_pub_qmap.py <anf.json> <L> [time_limit_s=600] [seed=0] [out.json]
"""
import itertools
import json
import sys
import time


def main():
    D = json.load(open(sys.argv[1]))
    anf = D['anf']
    L = int(sys.argv[2])
    tl = float(sys.argv[3]) if len(sys.argv) > 3 else 600
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    out_path = sys.argv[5] if len(sys.argv) > 5 else None
    n = max(m.bit_length() for mons in anf for m in mons)
    nout = len(anf)
    T = [set() for _ in range(nout)]
    for c, mons in enumerate(anf):
        for m in mons:
            bits = [i for i in range(n) if (m >> i) & 1]
            if len(bits) == 2:
                T[c].add((bits[0], bits[1]))
    from pycryptosat import Solver
    s = Solver(threads=1)
    cnt = [0]

    def new():
        cnt[0] += 1
        return cnt[0]
    l = [[new() for _ in range(n)] for _ in range(L)]
    lp = [[new() for _ in range(n)] for _ in range(L)]
    W = [[new() for _ in range(L)] for _ in range(nout)]

    def AND(p, q):
        r = new()
        s.add_clause([-r, p])
        s.add_clause([-r, q])
        s.add_clause([r, -p, -q])
        return r

    def XOR2(p, q):
        r = new()
        s.add_xor_clause([r, p, q], False)
        return r
    d = {}
    for i in range(L):
        for u in range(n):
            for v in range(u + 1, n):
                d[(i, u, v)] = XOR2(AND(l[i][u], lp[i][v]), AND(l[i][v], lp[i][u]))
    for c in range(nout):
        for u in range(n):
            for v in range(u + 1, n):
                terms = [AND(W[c][i], d[(i, u, v)]) for i in range(L)]
                s.add_xor_clause(terms, (u, v) in T[c])
    # symmetry: W columns in non-increasing lexicographic order
    for i in range(L - 1):
        eq = None
        for c in range(nout):
            p, q = W[c][i], W[c][i + 1]
            s.add_clause([p, -q] if eq is None else [-eq, p, -q])
            if c < nout - 1:
                e, x = new(), new()
                s.add_xor_clause([x, p, q], False)
                if eq is None:
                    s.add_clause([-e, -x])
                    s.add_clause([e, x])
                else:
                    s.add_clause([-e, eq])
                    s.add_clause([-e, -x])
                    s.add_clause([e, -eq, x])
                eq = e
    t0 = time.time()
    sat, sol = s.solve(time_limit=tl)
    res = dict(n=n, outputs=nout, L=L, seconds=round(time.time() - t0, 1),
               result='SAT' if sat else ('UNSAT' if sat is False else 'UNKNOWN'))
    if sat:
        lv = [[int(sol[x]) for x in row] for row in l]
        lpv = [[int(sol[x]) for x in row] for row in lp]
        Wv = [[int(sol[x]) for x in row] for row in W]
        # verify the alternating forms
        ok = True
        for c in range(nout):
            for u in range(n):
                for v in range(u + 1, n):
                    acc = 0
                    for i in range(L):
                        if Wv[c][i]:
                            acc ^= (lv[i][u] & lpv[i][v]) ^ (lv[i][v] & lpv[i][u])
                    ok &= acc == int((u, v) in T[c])
        res['verified'] = bool(ok)
        res['program'] = dict(l=[''.join(map(str, r)) for r in lv], lp=[''.join(map(str, r)) for r in lpv],
                              W=[''.join(map(str, r)) for r in Wv])
    line = json.dumps(res)
    print(line)
    if out_path:
        open(out_path, 'w').write(line + '\n')


if __name__ == '__main__':
    main()
