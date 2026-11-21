#!/usr/bin/env python3
"""invsqrt_fc_s18cubrank.py -- FC eighteenth holder: is a cubic partial function ONE gated CCZ away from Clifford?
Question: do affine forms L1, L2, L3 of the inputs exist with  target = Q + L1 L2 L3  (Q of degree <= 2) on every reached
pattern?  If yes, the phase (-1)^(m target) of an X-measured wire costs ONE CCZ under m (0.5 a shot) plus Cliffords.
SAT (pysat CaDiCaL): unknowns = the three forms' coefficients (inputs + constant) and Q's coefficients (1 + n + n(n-1)/2);
one XOR equation per reached pattern.  Symmetry is not broken (the answer is all that matters).
  invsqrt_fc_s18cubrank.py SPEC.json [K3=1] [K4=0]   (F3 / F1 'rows' format: {"inputs": [...], "rows": [[pattern, value], ...]})
  K3 products of three forms and K4 products of four forms (a quartic product = one gated AND into a temporary + one
  gated CCZ: 1.0 a shot; a cubic product = one gated CCZ: 0.5 a shot).
Safe to run: reads only, prints only."""
import itertools
import json
import sys
import time

from pysat.formula import IDPool
from pysat.solvers import Solver


def main():
    d = json.load(open(sys.argv[1]))
    names = d["inputs"]
    rows = d["rows"]
    n = len(names)
    pool = IDPool()
    cls = []
    T = pool.id("T")
    cls.append([T])

    def xor_chain(lits):
        if not lits:
            return -T
        acc = lits[0]
        for x in lits[1:]:
            c = pool.id()
            cls.extend([[-acc, -x, -c], [acc, x, -c], [acc, -x, c], [-acc, x, c]])
            acc = c
        return acc

    def and2(a, b):
        c = pool.id()
        cls.extend([[-c, a], [-c, b], [c, -a, -b]])
        return c
    K3 = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    K4 = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    shapes = [3] * K3 + [4] * K4
    L = [[[pool.id(("L", t, k, j)) for j in range(n + 1)] for k in range(sz)] for t, sz in enumerate(shapes)]
    quads = [()] + [(j,) for j in range(n)] + list(itertools.combinations(range(n), 2))
    Qc = [pool.id(("Q", i)) for i in range(len(quads))]
    for p, y in rows:
        p = int(p)
        bits = [(p >> j) & 1 for j in range(n)]
        prods = []
        for t, sz in enumerate(shapes):
            lv = []
            for k in range(sz):
                lits = [L[t][k][j] for j in range(n) if bits[j]] + [L[t][k][n]]
                lv.append(xor_chain(lits))
            pr = lv[0]
            for x in lv[1:]:
                pr = and2(pr, x)
            prods.append(pr)
        qlits = [Qc[i] for i, q in enumerate(quads) if all(bits[j] for j in q)]
        tot = xor_chain(qlits + prods)
        cls.append([tot] if int(y) else [-tot])
    t0 = time.time()
    with Solver(name="cadical153", bootstrap_with=cls) as S:
        ok = S.solve()
        print("%s: %d inputs, %d rows: Q + %d cubic + %d quartic products (%.1f a shot) %s (%.1fs)" % (
            sys.argv[1].rsplit("/", 1)[-1], n, len(rows), K3, K4, 0.5 * K3 + 1.0 * K4,
            "EXISTS" if ok else "impossible (UNSAT)", time.time() - t0), flush=True)
        if ok:
            m = set(x for x in S.get_model() if x > 0)
            for t, sz in enumerate(shapes):
                for k in range(sz):
                    print("  P%d L%d = %s" % (t + 1, k + 1, " ^ ".join([names[j] for j in range(n) if L[t][k][j] in m] +
                                                                  (["1"] if L[t][k][n] in m else [])) or "0"))


if __name__ == "__main__":
    main()
