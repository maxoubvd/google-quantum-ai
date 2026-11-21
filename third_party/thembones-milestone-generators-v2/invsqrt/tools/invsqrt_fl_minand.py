#!/usr/bin/env python3
"""invsqrt_fl_minand.py -- builder F2: exact minimum-AND search (SAT, pysat) for one output bit over given input
truth tables. Program with k ANDs: g_t = (affine over inputs and g_1..g_{t-1}) * (affine over the same);
output = affine over inputs and g_1..g_k. Rows = distinct input patterns (the target must be a function of them).

  invsqrt_fl_minand.py h2 <model.json> <second spec128.json> <set> <k> [timeout=600]
  invsqrt_fl_minand.py json <problem.json> <kmax>      generic: {"rows", "inputs": [{"name", "tt"}], "target"[, "mask"]}
  sets: s6 (D9..D12, H0, H1: the public patch 8 inputs), s8 (+ R, a5), s10 (+ D7, D8), s12 (+ D5, D6)
"""
import json
import sys
import time

from pysat.formula import IDPool
from pysat.solvers import Solver

NV = 4096


def rows_of(inputs, target):
    pats = {}
    for r in range(NV):
        p = tuple((t >> r) & 1 for t in inputs)
        y = (target >> r) & 1
        if pats.setdefault(p, y) != y:
            return None
    return pats


def encode(pats, n, k):
    pool = IDPool()
    cls = []
    T = pool.id("T")
    cls.append([T])                      # constant true literal
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
    # coefficient variables: for AND t, two forms over 1 + n + (t-1) terms; output form over 1 + n + k terms
    coef = {}
    for t in range(k):
        for side in "LR":
            for j in range(1 + n + t):
                coef[(t, side, j)] = pool.id(("c", t, side, j))
    for j in range(1 + n + k):
        coef[("o", j)] = pool.id(("o", j))
    for p, y in pats.items():
        g = []
        for t in range(k):
            vals = []
            for side in "LR":
                terms = [coef[(t, side, 0)]]
                terms += [coef[(t, side, 1 + i)] for i in range(n) if p[i]]
                terms += [and2(coef[(t, side, 1 + n + s)], g[s]) for s in range(t)]
                vals.append(xor_list(terms))
            g.append(and2(vals[0], vals[1]))
        terms = [coef[("o", 0)]] + [coef[("o", 1 + i)] for i in range(n) if p[i]]
        terms += [and2(coef[("o", 1 + n + s)], g[s]) for s in range(k)]
        o = xor_list(terms)
        cls.append([o] if y else [-o])
    return cls, coef, pool


SOLVER = "cadical153"     # json mode option solver=<pysat name> (glucose4, lingeling, maplechrono, minisat22, ...)


def generic(path, kmax, tmo, kmin=1):
    """JSON {"rows": N, "inputs": [{"name", "tt"}], "target": tt} (tt = hex, bit r = value at row r; only rows
    listed in an optional "mask" tt count). Tries k = kmin..kmax and prints the first SAT program."""
    d = json.load(open(path))
    n_rows = int(d["rows"])
    ins = [int(x["tt"], 16) for x in d["inputs"]]
    names = [x["name"] for x in d["inputs"]]
    tgt = int(d["target"], 16)
    mask = int(d.get("mask", "0"), 16) or (1 << n_rows) - 1
    pats = {}
    for r in range(n_rows):
        if not (mask >> r) & 1:
            continue
        p = tuple((t >> r) & 1 for t in ins)
        y = (tgt >> r) & 1
        if pats.setdefault(p, y) != y:
            print("target is not a function of these inputs")
            return
    for k in range(kmin, kmax + 1):
        cls, coef, pool = encode(pats, len(ins), k)
        t0 = time.time()
        with Solver(name=SOLVER, bootstrap_with=cls) as s:
            ok = s.solve()
            mdl = s.get_model() if ok else None
        print(f"k = {k}: {'SAT' if ok else 'UNSAT'} ({len(pats)} patterns, {len(cls)} clauses, "
              f"{time.time() - t0:.1f}s)", flush=True)
        if ok:
            pos = set(l for l in mdl if l > 0)
            for t in range(k):
                fs = []
                for side in "LR":
                    terms = (["1"] if coef[(t, side, 0)] in pos else [])
                    terms += [names[i] for i in range(len(ins)) if coef[(t, side, 1 + i)] in pos]
                    terms += [f"g{s_}" for s_ in range(t) if coef[(t, side, 1 + len(ins) + s_)] in pos]
                    fs.append(" ^ ".join(terms) or "0")
                print(f"  g{t} = ({fs[0]}) & ({fs[1]})")
            terms = (["1"] if coef[("o", 0)] in pos else [])
            terms += [names[i] for i in range(len(ins)) if coef[("o", 1 + i)] in pos]
            terms += [f"g{s_}" for s_ in range(k) if coef[("o", 1 + len(ins) + s_)] in pos]
            print("  out = " + " ^ ".join(terms))
            return


def main():
    if sys.argv[1] == "json":
        opt = dict(a.split("=", 1) for a in sys.argv[4:])
        global SOLVER
        SOLVER = opt.get("solver", SOLVER)
        generic(sys.argv[2], int(sys.argv[3]), float(opt.get("timeout", "600")), int(opt.get("kmin", "1")))
        return
    model = json.load(open(sys.argv[2]))
    sec = json.load(open(sys.argv[3]))
    which, k = sys.argv[4], int(sys.argv[5])
    tmo = float(dict(a.split("=", 1) for a in sys.argv[6:]).get("timeout", "600"))
    R128 = int(sec["inputs"][0]["tt"], 16)
    D, SS = model["D"], model["stored_slope"]

    def tt(f):
        return sum(1 << v for v in range(NV) if f(v))
    Db = [tt(lambda v, b=b: (D[v] >> b) & 1) for b in range(13)]
    H = [tt(lambda v, b=b: (SS[v >> 5] >> b) & 1) for b in range(5)]
    R = tt(lambda v: (R128 >> (v >> 5)) & 1)
    a5 = tt(lambda v: (v >> 10) & 1)
    sets = dict(s6=[Db[9], Db[10], Db[11], Db[12], H[0], H[1]])
    sets["s8"] = sets["s6"] + [R, a5]
    sets["s10"] = sets["s8"] + [Db[7], Db[8]]
    sets["s12"] = sets["s10"] + [Db[5], Db[6]]
    names = dict(s6="D9 D10 D11 D12 H0 H1", s8="+R a5", s10="+D7 D8", s12="+D5 D6")
    ins = sets[which]
    pats = rows_of(ins, H[2])
    if pats is None:
        print(f"{which}: H2 is not a function of these inputs")
        return
    cls, coef, pool = encode(pats, len(ins), k)
    t0 = time.time()
    with Solver(name="cadical153", bootstrap_with=cls) as s:
        ok = s.solve()
        model_ = s.get_model() if ok else None
    print(f"{which} ({names[which]}): {len(pats)} input patterns, k = {k}: {'SAT' if ok else 'UNSAT'} "
          f"({len(cls)} clauses, {time.time() - t0:.1f}s)")
    if ok:
        pos = set(l for l in model_ if l > 0)
        prog = {str(key): (coef[key] in pos) for key in coef}
        print(json.dumps({k2: v for k2, v in prog.items() if v}))


if __name__ == "__main__":
    main()
