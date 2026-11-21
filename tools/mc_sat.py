#!/usr/bin/env python3
"""Exact multiplicative complexity (minimum AND count of an XOR-AND graph) by SAT.

Plan step S2.24 (screens of DQI-B13 and CHEM-B19). A Boolean function F: {0,1}^n -> {0,1}^m has an XAG with k ANDs
when there are gates g_1..g_k, g_j = L_j^0 * L_j^1, each L a XOR of inputs and earlier gates, and each output a XOR
of a constant, inputs and gates. In a reversible kickmix circuit such an XAG costs k Toffolis forward and 0 to
uncompute (each AND target is erased by HMR plus a CZ fix-up), so for a phase oracle (-1)^f the price is MC(f).

Encoding (one SAT instance per budget k): selector bits a (inputs) and b (earlier gates) for both operands of every
gate, output selector bits; one value variable per gate per truth-table row; XORs by Tseitin chains; products
selector * gate value by 3-clause ANDs. Without loss of generality the gate operands carry no constant (a constant
on an operand adds an affine term that later operands and the outputs can absorb), operand 0 of each gate is
lexicographically below operand 1 is NOT imposed (kept simple), and for a single output the last gate must feed the
output (otherwise it is useless and k is not minimal).

Size: O(2^n * k * (n + k)) clauses. Practical for n <= 8 inputs and k <= ~10; above that use upper bounds from
constructions. Time and memory are reported per instance; instances for several k run in parallel (fork).

    ~/venvs/sat/bin/python tools/mc_sat.py trinv N [--kmin K] [--kmax K] [--timeout S]
    ~/venvs/sat/bin/python tools/mc_sat.py tt HEXTABLE NVARS [...]          # single-output truth table
"""
from __future__ import annotations

import argparse
import itertools
import json
import multiprocessing as mp
import sys
import time

from pysat.formula import IDPool
from pysat.solvers import Solver

# ---------------------------------------------------------------------------------------------- finite fields
IRRED = {2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011, 8: 0b100011011, 12: 0b1000001010011}


def gf_mul(x: int, y: int, n: int) -> int:
    poly, r = IRRED[n], 0
    while y:
        if y & 1:
            r ^= x
        y >>= 1
        x <<= 1
        if x >> n:
            x ^= poly
    return r


def gf_inv(x: int, n: int) -> int:
    """x^(2^n - 2); 0 -> 0."""
    r, e, b = 1, (1 << n) - 2, x
    while e:
        if e & 1:
            r = gf_mul(r, b, n)
        b = gf_mul(b, b, n)
        e >>= 1
    return r if x else 0


def gf_trace(x: int, n: int) -> int:
    t, y = 0, x
    for _ in range(n):
        t ^= y
        y = gf_mul(y, y, n)
    assert t in (0, 1)
    return t


def table_trinv(n: int) -> list[int]:
    """Truth table of Tr(a / b) on GF(2^n)^2, input index = a | b << n (2n variables)."""
    return [gf_trace(gf_mul(i & ((1 << n) - 1), gf_inv(i >> n, n), n), n) for i in range(1 << (2 * n))]


# ------------------------------------------------------------------------------------------------ SAT model
def _xor_chain(clauses, pool, lits, target_val=None):
    """Return a literal equal to XOR(lits) (lits may contain the constants True/False)."""
    const, vs = False, []
    for l in lits:
        if l is True:
            const = not const
        elif l is False:
            continue
        else:
            vs.append(l)
    if not vs:
        return const
    acc = vs[0]
    for l in vs[1:]:
        z = pool.id()
        clauses += [[-z, acc, l], [-z, -acc, -l], [z, -acc, l], [z, acc, -l]]
        acc = z
    return -acc if const else acc


def _and(clauses, pool, x, y):
    if x is False or y is False:
        return False
    if x is True:
        return y
    if y is True:
        return x
    z = pool.id()
    clauses += [[-z, x], [-z, y], [z, -x, -y]]
    return z


def build(tables: list[list[int]], nvars: int, k: int):
    """CNF for 'the m outputs in `tables` have an XAG with k ANDs'. Returns (clauses, pool, selectors)."""
    pool, cl = IDPool(), []
    rows = 1 << nvars
    A = [[[pool.id(("a", j, s, i)) for i in range(nvars)] for s in range(2)] for j in range(k)]
    B = [[[pool.id(("b", j, s, l)) for l in range(j)] for s in range(2)] for j in range(k)]
    m = len(tables)
    OC = [pool.id(("oc", o)) for o in range(m)]
    OX = [[pool.id(("ox", o, i)) for i in range(nvars)] for o in range(m)]
    OG = [[pool.id(("og", o, j)) for j in range(k)] for o in range(m)]
    # each operand non-trivial (some selector set); for single output, the last gate feeds the output
    for j in range(k):
        for s in range(2):
            cl.append(A[j][s] + B[j][s])
    if m == 1 and k:
        cl.append([OG[0][k - 1]])
    else:  # each gate used somewhere (later operand or an output)
        for j in range(k):
            cl.append([B[l][s][j] for l in range(j + 1, k) for s in range(2)] + [OG[o][j] for o in range(m)])
    for t in range(rows):
        x = [(t >> i) & 1 for i in range(nvars)]
        g = []
        for j in range(k):
            ops = []
            for s in range(2):
                lits = [A[j][s][i] for i in range(nvars) if x[i]]
                lits += [_and(cl, pool, B[j][s][l], g[l]) for l in range(j)]
                ops.append(_xor_chain(cl, pool, lits))
            g.append(_and(cl, pool, ops[0], ops[1]))
        for o in range(m):
            lits = [OC[o]] + [OX[o][i] for i in range(nvars) if x[i]] + [_and(cl, pool, OG[o][j], g[j]) for j in range(k)]
            out = _xor_chain(cl, pool, lits)
            want = tables[o][t]
            if out is True or out is False:
                if out != bool(want):
                    cl.append([])
            else:
                cl.append([out] if want else [-out])
    return cl, pool, (A, B, OC, OX, OG)


def decode(model, pool, sel, nvars, k):
    s = set(l for l in model if l > 0)
    A, B, OC, OX, OG = sel
    gates = [[([i for i in range(nvars) if A[j][q][i] in s], [l for l in range(j) if B[j][q][l] in s]) for q in range(2)]
             for j in range(k)]
    outs = [dict(c=int(OC[o] in s), x=[i for i in range(nvars) if OX[o][i] in s], g=[j for j in range(k) if OG[o][j] in s])
            for o in range(len(OC))]
    return dict(gates=gates, outputs=outs)


def check(circ, tables, nvars):
    """Evaluate a decoded XAG on every row; True if it reproduces all outputs."""
    for t in range(1 << nvars):
        x = [(t >> i) & 1 for i in range(nvars)]
        g = []
        for (a0, b0), (a1, b1) in circ["gates"]:
            v0 = (sum(x[i] for i in a0) + sum(g[l] for l in b0)) & 1
            v1 = (sum(x[i] for i in a1) + sum(g[l] for l in b1)) & 1
            g.append(v0 & v1)
        for o, out in enumerate(circ["outputs"]):
            v = (out["c"] + sum(x[i] for i in out["x"]) + sum(g[j] for j in out["g"])) & 1
            if v != tables[o][t]:
                return False
    return True


def solve(args):
    tables, nvars, k, timeout = args
    t0 = time.perf_counter()
    cl, pool, sel = build(tables, nvars, k)
    with Solver(name="cadical195", bootstrap_with=cl) as s:
        res = s.solve()  # the time limit is enforced by the parent (hard kill of the worker pool)
        out = dict(k=k, result={True: "SAT", False: "UNSAT", None: "TIMEOUT"}[res], vars=pool.top, clauses=len(cl),
                   seconds=round(time.perf_counter() - t0, 2))
        if res:
            circ = decode(s.get_model(), pool, sel, nvars, k)
            out["verified"] = check(circ, tables, nvars)
            out["circuit"] = circ
    return out


def degree(table, nvars):
    """Algebraic degree via the Moebius transform (lower bound: MC >= degree - 1)."""
    a = list(table)
    for i in range(nvars):
        for t in range(1 << nvars):
            if t >> i & 1:
                a[t] ^= a[t ^ (1 << i)]
    return max((bin(t).count("1") for t in range(1 << nvars) if a[t]), default=0)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("trinv")
    p.add_argument("n", type=int)
    p = sub.add_parser("tt")
    p.add_argument("hex")
    p.add_argument("nvars", type=int)
    for p in sub.choices.values():
        p.add_argument("--kmin", type=int)
        p.add_argument("--kmax", type=int, default=12)
        p.add_argument("--timeout", type=float, default=0)
        p.add_argument("--workers", type=int, default=6)
        p.add_argument("--json")
    a = ap.parse_args()
    if a.cmd == "trinv":
        nvars, table = 2 * a.n, table_trinv(a.n)
    else:
        nvars = a.nvars
        v = int(a.hex, 16)
        table = [(v >> t) & 1 for t in range(1 << nvars)]
    deg = degree(table, nvars)
    kmin = a.kmin if a.kmin is not None else max(deg - 1, 0)
    print(json.dumps(dict(cmd=a.cmd, nvars=nvars, degree=deg, lower_bound=max(deg - 1, 0))), flush=True)
    ks = list(range(kmin, a.kmax + 1))
    results = []
    # One worker per budget k, all started together; each k gets the same wall-clock deadline. An in-solver
    # interrupt did not stop CaDiCaL (S2.24 first run hung 56 min), so the pool is terminated instead.
    pool = mp.get_context("fork").Pool(min(a.workers, len(ks)))
    t0 = time.perf_counter()
    jobs = [(k, pool.apply_async(solve, (([table], nvars, k, a.timeout),))) for k in ks]
    try:
        for k, job in jobs:
            left = None if not a.timeout else max(1.0, a.timeout - (time.perf_counter() - t0))
            try:
                r = job.get(timeout=left)
            except mp.TimeoutError:
                r = dict(k=k, result="TIMEOUT", seconds=round(time.perf_counter() - t0, 1))
            print(json.dumps({x: r[x] for x in r if x != "circuit"}), flush=True)
            results.append(r)
            if r["result"] == "SAT":
                break
    finally:
        pool.terminate()
        pool.join()
    if a.json:
        json.dump(results, open(a.json, "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
