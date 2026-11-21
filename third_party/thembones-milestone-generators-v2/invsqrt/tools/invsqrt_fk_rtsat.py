"""F1 (third holder): SAT search for in-place router networks with FEW swap pairs.

n stages, stage s = a matching M_s on the 26 ring positions (free variables), class k fires the stages of its codeword
c_k (free variables, class 0 fires none). Requirement: the ordered product of the fired matchings of class k is the
cyclic shift x -> x + k (k = 0..10), i.e. the router of the endpoint in the labelling where class k is the shift by
k (w.l.o.g. when class 0 fires nothing). Caps: pairs per stage <= tmax, total pairs <= nmax.
Encoding: pos[k][s][x][y] = token x sits on y after stage s in class k (one-hot rows); stage transition
(c ? pos . M_s : pos) as 4-literal implications; final rows fixed to the shift.
    python invsqrt_fk_rtsat.py n tmax nmax [solver=cadical153] [seed=0] [out=FILE.json]
Prints SAT with the matchings and codewords, or UNSAT.
"""
import json
import sys
import time

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool, CNF
from pysat.solvers import Solver

N = 26
K = 11
W4 = [0]      # weight of every class codeword (0 = free); option w4=<w>
D4 = [0]      # minimum pairwise code distance (with w4); option d4=<d>
INV13 = [False]  # option inv13=1: stages invariant under translation by N/2
Q0FREE = [False]  # option q0free=1: class 0 may fire stages (final layouts compared up to the shift)
CODES = [None]  # fixed codewords; option codes=rm (f_k = sum_i flag_i(k) l_i on F_2^3, l = x1, x2, x3, x1+x2+x3+1)


def rm_codes(order=None):
    pts = order or list(range(8))
    def f(k, x):
        b = [(x >> i) & 1 for i in range(3)]
        ls = [b[0], b[1], b[2], b[0] ^ b[1] ^ b[2] ^ 1]
        o = 0
        for i in range(4):
            if (k >> i) & 1:
                o ^= ls[i]
        return o
    return ["".join(str(f(k, pts[s])) for s in range(8)) for k in range(K)]


def build(n, tmax, nmax):
    pool = IDPool()
    cnf = CNF()
    m = {}
    for s in range(n):
        for a in range(N):
            for b in range(a + 1, N):
                m[(s, a, b)] = pool.id(("m", s, a, b))

    def M(s, a, b):
        return m[(s, min(a, b), max(a, b))]
    u = {}
    for s in range(n):
        for y in range(N):
            u[(s, y)] = pool.id(("u", s, y))
            lits = [M(s, y, w) for w in range(N) if w != y]
            # u <-> OR lits ; and at most one pair per position
            cnf.append([-u[(s, y)]] + lits)
            for l in lits:
                cnf.append([-l, u[(s, y)]])
            for cl in CardEnc.atmost(lits, 1, vpool=pool, encoding=EncType.pairwise).clauses:
                cnf.append(cl)
        allm = [m[(s, a, b)] for a in range(N) for b in range(a + 1, N)]
        for cl in CardEnc.atmost(allm, tmax, vpool=pool, encoding=EncType.seqcounter).clauses:
            cnf.append(cl)
    allm = list(m.values())
    for cl in CardEnc.atmost(allm, nmax, vpool=pool, encoding=EncType.totalizer).clauses:
        cnf.append(cl)
    c = {}
    for k in range(K):
        for s in range(n):
            c[(k, s)] = pool.id(("c", k, s))
            if k == 0 and not Q0FREE[0]:
                cnf.append([-c[(k, s)]])
    if INV13[0]:                                # every stage invariant under x -> x + N/2 (commutes with shifts)
        h = N // 2
        for s in range(n):
            for a in range(N):
                for b in range(a + 1, N):
                    a2, b2 = (a + h) % N, (b + h) % N
                    l1, l2 = M(s, a, b), M(s, a2, b2)
                    if l1 != l2:
                        cnf.extend([[-l1, l2], [l1, -l2]])
    # symmetry (translations commute with every shift): position 0 is matched in SOME stage. (Before 06:00Z this
    # clause forced stage 0 itself to match position 0, which loses networks whose stage 0 is empty.)
    cnf.append([u[(s, 0)] for s in range(n)])
    if CODES[0]:                                # fixed codewords (list of strings, class 0 = zeros)
        for k in range(K):
            for s in range(n):
                cnf.append([c[(k, s)] if CODES[0][k][s] == "1" else -c[(k, s)]])
    if W4[0]:                                   # RM(1,3)-like codes: weight W4 per class >= 1, pairwise distance >= D4
        for k in range(1, K):
            for cl in CardEnc.equals([c[(k, s)] for s in range(n)], W4[0], vpool=pool,
                                     encoding=EncType.seqcounter).clauses:
                cnf.append(cl)
        for k in range(1, K):
            for k2 in range(k + 1, K):
                d = []
                for s in range(n):
                    e = pool.id(("d", k, k2, s))
                    a, b = c[(k, s)], c[(k2, s)]
                    cnf.extend([[-e, a, b], [-e, -a, -b], [e, -a, b], [e, a, -b]])
                    d.append(e)
                for cl in CardEnc.atleast(d, D4[0], vpool=pool, encoding=EncType.seqcounter).clauses:
                    cnf.append(cl)
    pos = {}
    for k in range(K):
        for x in range(N):
            for y in range(N):
                pos[(k, 0, x, y)] = None          # constant [x == y]
        for s in range(n):
            for x in range(N):
                for z in range(N):
                    pos[(k, s + 1, x, z)] = pool.id(("p", k, s + 1, x, z))
            for x in range(N):
                row = [pos[(k, s + 1, x, z)] for z in range(N)]
                for cl in CardEnc.equals(row, 1, vpool=pool, encoding=EncType.seqcounter).clauses:
                    cnf.append(cl)
            cs = c[(k, s)]
            for x in range(N):
                for y in range(N):
                    py = pos[(k, s, x, y)]
                    # py is either a literal or a constant
                    if py is None:
                        if x != y:
                            continue
                        pre = []                       # true
                    else:
                        pre = [-py]
                    # not fired: pos' = pos
                    cnf.append(pre + [cs, pos[(k, s + 1, x, y)]])
                    # fired: y -> z with M(y,z), or stays if unmatched
                    cnf.append(pre + [-cs, u[(s, y)], pos[(k, s + 1, x, y)]])
                    for z in range(N):
                        if z != y:
                            cnf.append(pre + [-cs, -M(s, y, z), pos[(k, s + 1, x, z)]])
        if not Q0FREE[0]:
            for x in range(N):
                cnf.append([pos[(k, n, x, (x + k) % N)]])
    if Q0FREE[0]:                               # general: token x of class k ends where token x + k of class 0 ends
        for k in range(1, K):                   # (P_k(x) = L(x + k) = P_0(x + k); the word layout is the labelling)
            for x in range(N):
                for z in range(N):
                    a, b = pos[(k, n, x, z)], pos[(0, n, (x + k) % N, z)]
                    cnf.extend([[-a, b], [a, -b]])
    return cnf, pool, m, c


def main(n, tmax, nmax, *opts):
    n, tmax, nmax = int(n), int(tmax), int(nmax)
    kw = dict(o.split("=") for o in opts)
    global N, K
    N = int(kw.get("ring", N))
    K = int(kw.get("classes", K))
    W4[0] = int(kw.get("w4", 0))
    D4[0] = int(kw.get("d4", 0))
    Q0FREE[0] = kw.get("q0free", "0") == "1"
    INV13[0] = kw.get("inv13", "0") == "1"
    if kw.get("codes") == "rm":
        order = [int(x) for x in kw["order"].split(",")] if "order" in kw else None
        CODES[0] = rm_codes(order)
        print("codes", CODES[0])
    t0 = time.time()
    cnf, pool, m, c = build(n, tmax, nmax)
    print(f"n {n} tmax {tmax} nmax {nmax}: {pool.top} vars, {len(cnf.clauses)} clauses ({time.time() - t0:.1f}s)",
          flush=True)
    with Solver(name=kw.get("solver", "cadical153"), bootstrap_with=cnf.clauses) as S:
        ok = S.solve()
        print("SAT" if ok else "UNSAT", f"({time.time() - t0:.1f}s)", flush=True)
        if ok:
            model = set(l for l in S.get_model() if l > 0)
            stages = [[(a, b) for a in range(N) for b in range(a + 1, N) if m[(s, a, b)] in model] for s in range(n)]
            codes = ["".join("1" if c[(k, s)] in model else "0" for s in range(n)) for k in range(K)]
            print("pairs per stage", [len(x) for x in stages], "total", sum(len(x) for x in stages))
            print("codes", codes)
            res = dict(n=n, tmax=tmax, nmax=nmax, stages=stages, codes=codes)
            if "out" in kw:
                json.dump(res, open(kw["out"], "w"))
            print(json.dumps(res))


if __name__ == "__main__":
    main(*sys.argv[1:])
