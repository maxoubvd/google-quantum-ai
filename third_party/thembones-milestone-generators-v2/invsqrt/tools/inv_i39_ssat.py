#!/usr/bin/env python3
"""inv_i39_ssat.py -- I39 (builder, second holder): minimum-AND rebuild programs in GHOST's STAGED shape, searched
as ONE SAT problem: TEMPLATE "k0:I0,k1:I1,..." = stages run one after the other; stage s has k_s gates, its first I_s
on temporaries (each reads only earlier temporaries of its own stage), the rest CCX straight into the target (they
read only their stage's temporaries and are read by nobody); the target is the XOR of every terminal gate, any
temporaries and inputs.  Width = the largest I_s (the stages do not overlap), cost = the sum of k_s.  The answer is
written in GHOST's staged format (r63_inv_ghost_state.py replay / the staged emitter) and checked on every row.
F1's exact encoding (tools/invsqrt_fk_minand_sb.py, imported read-only) plus unit clauses for the shape.
  python3 inv_i39_ssat.py TEMPLATE PROBLEM.json [solver=NAME] [out=PROG.json] [name=q18]
  python3 inv_i39_ssat.py TEMPLATE PROBLEM.json dump=X.cnf             (DIMACS only, for kissat / cadical binaries)
  python3 inv_i39_ssat.py TEMPLATE PROBLEM.json model=SOLVER.out out=PROG.json   (decode a binary's model, check it)
Safe to run: reads the problem, writes out= only."""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_fk_minand_sb as MB  # noqa: E402
from pysat.solvers import Solver  # noqa: E402


def layout(template):
    st = []
    for part in template.split(","):
        k, i = (int(x) for x in part.split(":"))
        assert 1 <= k and 0 <= i <= k
        st.append((k, i))
    gates = []                                       # global gate t -> (stage, local j, internal?)
    for s, (k, i) in enumerate(st):
        for j in range(k):
            gates.append((s, j, j < i))
    return st, gates


def main():
    template, path = sys.argv[1], sys.argv[2]
    opts = dict(a.split("=", 1) for a in sys.argv[3:] if "=" in a)
    st, gates = layout(template)
    K = len(gates)
    names, pats, m = MB.load([path])
    assert m == 1
    n = len(names)
    t0 = time.time()
    cls, coef = MB.encode(pats, n, m, K)
    first = {}
    for t, (s, j, _) in enumerate(gates):
        first.setdefault(s, t)
    for t2, (s2, j2, int2) in enumerate(gates):
        for t, (s, j, intl) in enumerate(gates[:t2]):
            ok = (s == s2) and intl and (int2 is False or j < j2)
            if not ok:
                for side in "LR":
                    cls.append([-coef[(t2, side, n + t)]])
        if not int2:                                 # a terminal gate goes into the target
            cls.append([coef[("o", 0, 1 + n + t2)]])
    if "g0in" in opts:                               # (third holder) gate 0 must equal one of the listed products
        top = max(abs(x) for c in cls for x in c)
        sel = []
        for ln in open(opts["g0in"]):
            if "|" not in ln:
                continue
            sides = []
            for s in ln.split("|"):
                items = [x.strip() for x in s.split("^") if x.strip()]
                if any(x not in names for x in items):
                    sides = None
                    break
                sides.append(set(names.index(x) for x in items))
            if not sides:
                continue
            top += 1
            v = top
            sel.append(v)
            for side, js in zip("LR", sides):
                for j in range(n):
                    cls.append([-v, coef[(0, side, j)] if j in js else -coef[(0, side, j)]])
        print("g0in: %d candidate products" % len(sel), flush=True)
        cls.append(sel)
    if "dump" in opts:                               # DIMACS for a borrowed box (solver input file only)
        nv = max(abs(x) for c in cls for x in c)
        with open(opts["dump"], "w") as f:
            f.write("p cnf %d %d\n" % (nv, len(cls)))
            for c in cls:
                f.write(" ".join(map(str, c)) + " 0\n")
        print("dumped", opts["dump"], nv, "vars", len(cls), "clauses")
        return
    if "model" in opts:                              # decode a solver's output (lines 'v ...' / 's SATISFIABLE')
        txt = open(opts["model"]).read()
        if "s SATISFIABLE" not in txt:
            print("no model in", opts["model"], "(UNSAT or unfinished)")
            return
        lits = [int(x) for ln in txt.splitlines() if ln.startswith("v ") for x in ln[2:].split()]
        model = set(x for x in lits if x > 0)
        print("template %s: model read from %s" % (template, opts["model"]))
    else:
        with Solver(name=opts.get("solver", "cadical153"), bootstrap_with=cls) as S:
            ok = S.solve()
            print("template %s (K %d, width %d): %s (%d patterns, %d clauses, %.1fs)" % (
                template, K, max(i for _, i in st), "SAT" if ok else "UNSAT", len(pats), len(cls), time.time() - t0),
                flush=True)
            if not ok:
                return
            model = set(v for v in S.get_model() if v > 0)
    prog = MB.extract(model, coef, n, m, K)
    assert all(MB.evaluate(prog, n, p) == y for p, y in pats.items())
    # GHOST's staged format: per stage, gate indices = inputs 0..n-1 then the stage's own gates n + j
    out_all = prog["outputs"][0]
    stages = []
    for s, (k, i) in enumerate(st):
        ts = [t for t, g in enumerate(gates) if g[0] == s]
        loc = {n + t: n + gates[t][1] for t in ts}
        sg = []
        for t in ts:
            g = prog["gates"][t]
            sg.append(dict(left=[x if x < n else loc[x] for x in g["left"]],
                           right=[x if x < n else loc[x] for x in g["right"]]))
        out = [loc[x] for x in out_all if x >= n and x in loc and gates[x - n][2]]
        if s == 0:
            out = [x for x in out_all if x < n] + out          # inputs and the constant (-1) once
        stages.append(dict(gates=sg, internal=i, out=out))
    sprog = dict(inputs=names, internal=0, stages=stages, outputs=[[]], target_name=opts.get("name", "q18"),
                 note="I39 inv_i39_ssat.py template %s, %d ANDs, width %d" % (template, K, max(i for _, i in st)))
    # replay the staged form on every pattern (the staged emitter's reading: stage = out + its terminal gates)
    for p, y in pats.items():
        tot = 0
        for stg in sprog["stages"]:
            vals = list(p)
            for g in stg["gates"]:
                a = 0
                for x in g["left"]:
                    a ^= 1 if x == -1 else vals[x]
                b = 0
                for x in g["right"]:
                    b ^= 1 if x == -1 else vals[x]
                vals.append(a & b)
            o = 0
            for x in list(stg["out"]) + [n + j for j in range(stg["internal"], len(stg["gates"]))]:
                o ^= 1 if x == -1 else vals[x]
            tot ^= o
        assert tot == y[0], "staged replay differs"
    print("staged replay: exact on all %d patterns" % len(pats))
    if "out" in opts:
        json.dump(sprog, open(opts["out"], "w"), indent=1)
        print("written", opts["out"])


if __name__ == "__main__":
    main()
