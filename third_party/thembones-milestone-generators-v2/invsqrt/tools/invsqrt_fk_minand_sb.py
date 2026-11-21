#!/usr/bin/env python3
"""F1: exact minimum-AND search with symmetry breaking (copy of builder F2's encoding in tools/invsqrt_fl_minand.py,
multi-output as tools/invsqrt_fk_minandm.py). Two sound reductions:
 (1) no constants inside the factors: (L^a)(R^b) = L R ^ b L ^ a R ^ a b, and the affine rest folds into every later
     factor and output (which are affine anyway), so a program with constants maps to one without, same AND count;
 (2) each factor pair (L, R) in reduced echelon form with respect to the variable order (inputs, then earlier gates):
     pivot(L) < pivot(R), R is 0 at pivot(L), L is 0 at pivot(R). Any basis {L', R'} of span{L, R} gives
     L' R' = L R ^ (a linear form of L, R), again absorbed by later forms; L, R non-zero and different.
Outputs keep their constant (it is free). Prints the program like F2's tool; out=prog.json as minandm.
    python invsqrt_fk_minand_sb.py <kmax> <p1.json> [<p2.json> ...] [kmin=K] [solver=NAME] [out=prog.json]
"""
import json
import sys
import time

from pysat.formula import IDPool
from pysat.solvers import Solver

sys.path.insert(0, "<HOME>/Quantum/tools")


def load(paths):
    ds = [json.load(open(p)) for p in paths]
    rows = int(ds[0]["rows"])
    names = [x["name"] for x in ds[0]["inputs"]]
    ins = [int(x["tt"], 16) for x in ds[0]["inputs"]]
    tg = [int(d["target"], 16) for d in ds]
    pats = {}
    for r in range(rows):
        p = tuple((t >> r) & 1 for t in ins)
        y = tuple((t >> r) & 1 for t in tg)
        assert pats.setdefault(p, y) == y
    return names, pats, len(tg)


def encode(pats, n, m, k, ret_rows=False):
    pool = IDPool()
    cls = []
    T = pool.id("T")
    cls.append([T])
    F = -T

    def xor2(a, b):
        c = pool.id()
        cls.extend([[-a, -b, -c], [a, b, -c], [a, -b, c], [-a, b, c]])
        return c

    def and2(a, b):
        c = pool.id()
        cls.extend([[-c, a], [-c, b], [c, -a, -b]])
        return c

    def xor_list(lits):
        if not lits:
            return F
        acc = lits[0]
        for x in lits[1:]:
            acc = xor2(acc, x)
        return acc
    coef = {}
    for t in range(k):
        for side in "LR":
            for j in range(n + t):                   # no constant (reduction 1)
                coef[(t, side, j)] = pool.id(("c", t, side, j))
    for o in range(m):
        for j in range(1 + n + k):
            coef[("o", o, j)] = pool.id(("o", o, j))
    # reduction 2: reduced echelon pair per gate
    for t in range(k):
        M = n + t
        piv = {}
        for side in "LR":
            seen = None                              # literal: some 1 strictly before j
            for j in range(M):
                x = coef[(t, side, j)]
                p = pool.id(("p", t, side, j))
                # p <-> x & ~seen
                if seen is None:
                    cls.extend([[-p, x], [p, -x]])
                else:
                    cls.extend([[-p, x], [-p, -seen], [p, -x, seen]])
                ns = pool.id(("s", t, side, j))      # ns <-> seen | x
                if seen is None:
                    cls.extend([[-ns, x], [ns, -x]])
                else:
                    cls.extend([[-ns, seen, x], [ns, -seen], [ns, -x]])
                seen = ns
                piv[(side, j)] = p
            cls.append([seen])                       # non-zero
        sL = None
        for j in range(M):
            pl, pr = piv[("L", j)], piv[("R", j)]
            # pivot(R) > pivot(L): pr -> some pl before j
            if sL is None:
                cls.append([-pr])
            else:
                cls.append([-pr, sL])
            cls.append([-pl, -coef[(t, "R", j)]])    # R zero at L's pivot
            cls.append([-pr, -coef[(t, "L", j)]])    # L zero at R's pivot
            nsl = pool.id(("sl", t, j))
            if sL is None:
                cls.extend([[-nsl, pl], [nsl, -pl]])
            else:
                cls.extend([[-nsl, sL, pl], [nsl, -sL], [nsl, -pl]])
            sL = nsl
    def row_clauses(p, ys):
        start = len(cls)
        g = []
        for t in range(k):
            vals = []
            for side in "LR":
                terms = [coef[(t, side, i)] for i in range(n) if p[i]]
                terms += [and2(coef[(t, side, n + s)], g[s]) for s in range(t)]
                vals.append(xor_list(terms))
            g.append(and2(vals[0], vals[1]))
        for o in range(m):
            terms = [coef[("o", o, 0)]] + [coef[("o", o, 1 + i)] for i in range(n) if p[i]]
            terms += [and2(coef[("o", o, 1 + n + s)], g[s]) for s in range(k)]
            v = xor_list(terms)
            cls.append([v] if ys[o] else [-v])
        out = cls[start:]
        del cls[start:]
        return out
    if ret_rows:
        return cls, coef, row_clauses
    for p, ys in pats.items():
        cls.extend(row_clauses(p, ys))
    return cls, coef


def evaluate(prog, n, p):
    vals = list(p)
    for g in prog["gates"]:
        a = 0
        for j in g["left"]:
            a ^= vals[j]
        b = 0
        for j in g["right"]:
            b ^= vals[j]
        vals.append(a & b)
    outs = []
    for lst in prog["outputs"]:
        o = 0
        for j in lst:
            o ^= 1 if j == -1 else vals[j]
        outs.append(o)
    return tuple(outs)


def extract(mdl, coef, n, m, k):
    prog = dict(gates=[], outputs=[])
    for t in range(k):
        sides = [[j for j in range(n + t) if coef[(t, side, j)] in mdl] for side in "LR"]
        prog["gates"].append(dict(left=sides[0], right=sides[1]))
    for o in range(m):
        lst = ([-1] if coef[("o", o, 0)] in mdl else [])
        lst += [i for i in range(n) if coef[("o", o, 1 + i)] in mdl]
        lst += [n + s_ for s_ in range(k) if coef[("o", o, 1 + n + s_)] in mdl]
        prog["outputs"].append(lst)
    return prog


def solve_cegar(pats, n, m, k, solver, start=24, add=8, seed=1):
    """Counterexample-guided: rows are added incrementally; UNSAT on a subset is UNSAT on all rows."""
    import random
    rnd = random.Random(seed)
    keys = list(pats)
    rnd.shuffle(keys)
    base, coef, row_clauses = encode(pats, n, m, k, ret_rows=True)
    used = set()
    with Solver(name=solver, bootstrap_with=base) as s:
        for p in keys[:start]:
            s.append_formula(row_clauses(p, pats[p])); used.add(p)
        while True:
            if not s.solve():
                return None, len(used)
            mdl = set(l for l in s.get_model() if l > 0)
            prog = extract(mdl, coef, n, m, k)
            bad = [p for p in keys if p not in used and evaluate(prog, n, p) != pats[p]]
            if not bad:
                return prog, len(used)
            for p in bad[:add]:
                s.append_formula(row_clauses(p, pats[p])); used.add(p)


def main():
    args = [a for a in sys.argv[1:] if "=" not in a]
    opt = dict(a.split("=", 1) for a in sys.argv[1:] if "=" in a)
    kmax = int(args[0])
    names, pats, m = load(args[1:])
    n = len(names)
    nm = lambda j: names[j] if j < n else f"g{j - n}"
    solver = opt.get("solver", "cadical153")
    for k in range(int(opt.get("kmin", "1")), kmax + 1):
        t0 = time.time()
        if opt.get("cegar", "0") == "1":     # measured slower on these problems (s14 k = 3: 201 s vs 3.7 s)
            prog, used = solve_cegar(pats, n, m, k, solver, seed=int(opt.get("seed", "1")))
            ok = prog is not None
            info = f"{used} of {len(pats)} patterns used"
        else:
            cls, coef = encode(pats, n, m, k)
            with Solver(name=solver, bootstrap_with=cls) as s:
                ok = s.solve()
                prog = extract(set(l for l in s.get_model() if l > 0), coef, n, m, k) if ok else None
            info = f"{len(pats)} patterns, {len(cls)} clauses"
        print(f"k = {k}: {'SAT' if ok else 'UNSAT'} ({info}, {m} outputs, {time.time() - t0:.1f}s)", flush=True)
        if ok:
            assert all(evaluate(prog, n, p) == y for p, y in pats.items())
            prog["inputs"] = names
            for t, g in enumerate(prog["gates"]):
                print(f"  g{t} = ({' ^ '.join(map(nm, g['left']))}) & ({' ^ '.join(map(nm, g['right']))})")
            for o, lst in enumerate(prog["outputs"]):
                print(f"  out{'' if m == 1 else o} = " + " ^ ".join("1" if j == -1 else nm(j) for j in lst))
            if "out" in opt:
                json.dump(prog, open(opt["out"], "w"), indent=1)
            return


if __name__ == "__main__":
    main()
