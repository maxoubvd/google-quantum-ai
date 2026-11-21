#!/usr/bin/env python3
"""F1: exact minimum-AND search for SEVERAL output bits sharing one straight-line program (multi-output copy of
builder F2's tools/invsqrt_fl_minand.py encoding). Program with k ANDs: g_t = (affine over inputs and g_1..g_{t-1})
* (affine over the same); each output = affine over inputs and g_1..g_k. Problems: F2's JSON format, one per output,
with the same inputs and rows (tools/invsqrt_fk_signsat.py with the same input list writes them so).
    python invsqrt_fk_minandm.py <kmax> <p1.json> <p2.json> [...] [kmin=K] [solver=NAME] [out=prog.json]
Writes the first SAT program as JSON {inputs, outputs: [names], gates: [{left, right}] (index lists, constants as
-1), outputs_forms: [[idx...]]} when out= is given.
"""
import json
import sys
import time

from pysat.formula import IDPool
from pysat.solvers import Solver


def load(paths):
    ds = [json.load(open(p)) for p in paths]
    rows = int(ds[0]["rows"])
    names = [x["name"] for x in ds[0]["inputs"]]
    ins = [int(x["tt"], 16) for x in ds[0]["inputs"]]
    for d in ds[1:]:
        assert int(d["rows"]) == rows and [x["name"] for x in d["inputs"]] == names
        assert [int(x["tt"], 16) for x in d["inputs"]] == ins
    tg = [int(d["target"], 16) for d in ds]
    pats = {}
    for r in range(rows):
        p = tuple((t >> r) & 1 for t in ins)
        y = tuple((t >> r) & 1 for t in tg)
        assert pats.setdefault(p, y) == y
    return names, pats, len(tg)


def encode(pats, n, m, k):
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
            for j in range(1 + n + t):
                coef[(t, side, j)] = pool.id(("c", t, side, j))
    for o in range(m):
        for j in range(1 + n + k):
            coef[("o", o, j)] = pool.id(("o", o, j))
    for p, ys in pats.items():
        g = []
        for t in range(k):
            vals = []
            for side in "LR":
                terms = [coef[(t, side, 0)]]
                terms += [coef[(t, side, 1 + i)] for i in range(n) if p[i]]
                terms += [and2(coef[(t, side, 1 + n + s)], g[s]) for s in range(t)]
                vals.append(xor_list(terms))
            g.append(and2(vals[0], vals[1]))
        for o in range(m):
            terms = [coef[("o", o, 0)]] + [coef[("o", o, 1 + i)] for i in range(n) if p[i]]
            terms += [and2(coef[("o", o, 1 + n + s)], g[s]) for s in range(k)]
            v = xor_list(terms)
            cls.append([v] if ys[o] else [-v])
    return cls, coef


def main():
    args = [a for a in sys.argv[1:] if "=" not in a]
    opt = dict(a.split("=", 1) for a in sys.argv[1:] if "=" in a)
    kmax = int(args[0])
    names, pats, m = load(args[1:])
    n = len(names)
    for k in range(int(opt.get("kmin", "1")), kmax + 1):
        cls, coef = encode(pats, n, m, k)
        t0 = time.time()
        with Solver(name=opt.get("solver", "cadical153"), bootstrap_with=cls) as s:
            ok = s.solve()
            mdl = set(l for l in s.get_model() if l > 0) if ok else None
        print(f"k = {k}: {'SAT' if ok else 'UNSAT'} ({len(pats)} patterns, {m} outputs, {len(cls)} clauses, "
              f"{time.time() - t0:.1f}s)", flush=True)
        if ok:
            def form(key0):
                idx = [-1] if coef[key0(0)] in mdl else []
                return idx
            prog = dict(inputs=names, gates=[], outputs=[])
            for t in range(k):
                sides = []
                for side in "LR":
                    lst = ([-1] if coef[(t, side, 0)] in mdl else [])
                    lst += [i for i in range(n) if coef[(t, side, 1 + i)] in mdl]
                    lst += [n + s_ for s_ in range(t) if coef[(t, side, 1 + n + s_)] in mdl]
                    sides.append(lst)
                prog["gates"].append(dict(left=sides[0], right=sides[1]))
                nm = lambda j: "1" if j == -1 else (names[j] if j < n else f"g{j - n}")
                print(f"  g{t} = ({' ^ '.join(map(nm, sides[0])) or '0'}) & ({' ^ '.join(map(nm, sides[1])) or '0'})")
            for o in range(m):
                lst = ([-1] if coef[("o", o, 0)] in mdl else [])
                lst += [i for i in range(n) if coef[("o", o, 1 + i)] in mdl]
                lst += [n + s_ for s_ in range(k) if coef[("o", o, 1 + n + s_)] in mdl]
                prog["outputs"].append(lst)
                nm = lambda j: "1" if j == -1 else (names[j] if j < n else f"g{j - n}")
                print(f"  out{o} = {' ^ '.join(map(nm, lst))}")
            if "out" in opt:
                json.dump(prog, open(opt["out"], "w"), indent=1)
            return


if __name__ == "__main__":
    main()
