#!/usr/bin/env python3
"""dqi_pub_xag.py -- builder D1 (holder 4): exact XOR-AND synthesis (SAT) of the GF(16) inverse-trace leaf.

Targets: h_j(C) = Tr16(beta^j C^-1), j = 0..3, C in GF(16) = F2[beta]/(beta^4 + beta + 1) as 4 bits (C = 0 -> 0).
Circuit: k AND gates; gate g sees linear forms (no constant) of the inputs and of earlier gates; the outputs are
linear forms of inputs and gates. Then (-1)^Tr16(B / C) = prod_j (-1)^(B_j h_j(C)): CZ of B_j with an output form.
Usage: dqi_pub_xag.py <k> [time_limit_s=600] [out.json]
"""
import json
import sys
import time


def gmul(a, b):
    r = 0
    for i in range(4):
        if (b >> i) & 1:
            r ^= a << i
    for i in (6, 5, 4):
        if (r >> i) & 1:
            r ^= 0b10011 << (i - 4)
    return r


def ginv(a):
    if not a:
        return 0
    return next(x for x in range(1, 16) if gmul(a, x) == 1)


def gtr(a):
    s, v = 0, a
    for _ in range(4):
        s ^= v
        v = gmul(v, v)
    return s & 1


def targets():
    return [[gtr(gmul(1 << j, ginv(c))) for c in range(16)] for j in range(4)]


def main():
    k = int(sys.argv[1])
    tl = float(sys.argv[2]) if len(sys.argv) > 2 else 600
    out_path = sys.argv[3] if len(sys.argv) > 3 else None
    n, m = 4, 4
    T = targets()
    from pycryptosat import Solver
    s = Solver(threads=1)
    cnt = [0]

    def new():
        cnt[0] += 1
        return cnt[0]

    def AND(p, q):
        if p is True:
            return q
        if q is True:
            return p
        if p is False or q is False:
            return False
        r = new()
        s.add_clause([-r, p])
        s.add_clause([-r, q])
        s.add_clause([r, -p, -q])
        return r

    def XOR(lits, rhs_target=None):
        """returns a literal equal to the XOR of lits (True / False constants folded)."""
        const = False
        vs = []
        for x in lits:
            if x is True:
                const = not const
            elif x is False:
                continue
            else:
                vs.append(x)
        if not vs:
            return const
        r = new()
        s.add_xor_clause([r] + vs, const)       # r = xor(vs) xor const
        return r
    A = [[new() for _ in range(n + j)] for j in range(k)]
    Bs = [[new() for _ in range(n + j)] for j in range(k)]
    Cs = [[new() for _ in range(n + k)] for _ in range(m)]
    for t in range(16):
        sig = [bool((t >> i) & 1) for i in range(n)]
        for j in range(k):
            a = XOR([AND(A[j][i], sig[i]) for i in range(n + j)])
            b = XOR([AND(Bs[j][i], sig[i]) for i in range(n + j)])
            sig.append(AND(a, b))
        for o in range(m):
            terms = [AND(Cs[o][i], sig[i]) for i in range(n + k)]
            vs = [x for x in terms if x not in (True, False)]
            const = sum(1 for x in terms if x is True) & 1
            want = bool(T[o][t]) ^ bool(const)
            if vs:
                s.add_xor_clause(vs, want)
            elif want:
                s.add_clause([])
    t0 = time.time()
    sat, sol = s.solve(time_limit=tl)
    res = dict(k=k, seconds=round(time.time() - t0, 1), result='SAT' if sat else ('UNSAT' if sat is False else 'UNKNOWN'))
    if sat:
        prog = dict(A=[[int(sol[x]) for x in row] for row in A], B=[[int(sol[x]) for x in row] for row in Bs],
                    C=[[int(sol[x]) for x in row] for row in Cs])
        # verify by simulation
        ok = True
        for t in range(16):
            sig = [(t >> i) & 1 for i in range(n)]
            for j in range(k):
                a = sum(prog['A'][j][i] & sig[i] for i in range(n + j)) & 1
                b = sum(prog['B'][j][i] & sig[i] for i in range(n + j)) & 1
                sig.append(a & b)
            for o in range(m):
                ok &= (sum(prog['C'][o][i] & sig[i] for i in range(n + k)) & 1) == T[o][t]
        res['verified'] = bool(ok)
        res['program'] = prog
    line = json.dumps(res)
    print(line)
    if out_path:
        open(out_path, 'w').write(line + '\n')


if __name__ == '__main__':
    main()
