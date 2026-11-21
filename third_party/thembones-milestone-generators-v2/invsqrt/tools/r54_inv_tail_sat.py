#!/usr/bin/env python3
"""r54_inv_tail_sat.py -- R54-INV-TAIL (designer): RELATIONAL minimum-AND search for the last j rows of a divider
tail, with the tolerance's don't-cares.  The gate encoding and its two symmetry reductions are builder F1's
(tools/invsqrt_fk_minand_sb.py: no constants in the factors; each factor pair in reduced echelon form), copied; the
output condition is new: per class of the map's inputs the output WORD has to be one of the allowed words (every
disallowed word of a reachable class is blocked by one clause); unreachable patterns carry no clause at all.

  r54_inv_tail_sat.py TABLE.npz K M0 [solver=cadical153] [out=prog.json] [keep0=1] [dump=file.cnf]
      TABLE.npz  from tools/r54_inv_tail_dc.py (minslack per class and word)
      K          the number of ANDs;   M0: a word is allowed on a class iff minslack >= M0 ("margin" = today's margin,
                 "m/1.35", "m/2", or a number; a negative number allows words that today's angles miss by that much)
      keep0=1    the first digit (the row's own sign wire) is kept as it is (as the rows do)
      exact=1    only today's word is allowed (the exact map: control)
      inj=1      (free gates) the map stays one-to-one on the window: classes with the same v bits and the same
                 lowest window bit get different words (necessary for an in-place tail with no held wire)
      fam=1      THE ROW FAMILY instead of free gates (K is ignored): today's rows (framed adds, same wires, same
                 carries, so the same AND count forward and in the mirror, in place by construction) with every
                 carry-in and operand literal a free affine form of the v bits the class reads (digits=1: and of the
                 sub-map's digit wires); fixops=1: only the carry-in literals free; taper=r[,r]: row r (0-based
                 among the j rows) one target narrower (one AND less)
Prints `k = K: SAT / UNSAT` and the program.  Safe to run: reads the table, writes only out= / dump=."""
import json
import sys
import time

import numpy as np
from pysat.formula import IDPool
from pysat.solvers import Solver


def encode_gates(n, m, k):
    pool = IDPool()
    cls = []
    T = pool.id("T")
    cls.append([T])
    coef = {}
    for t in range(k):
        for side in "LR":
            for j in range(n + t):
                coef[(t, side, j)] = pool.id(("c", t, side, j))
    for o in range(m):
        for j in range(1 + n + k):
            coef[("o", o, j)] = pool.id(("o", o, j))
    for t in range(k):
        M = n + t
        piv = {}
        for side in "LR":
            seen = None
            for j in range(M):
                x = coef[(t, side, j)]
                p = pool.id(("p", t, side, j))
                if seen is None:
                    cls.extend([[-p, x], [p, -x]])
                else:
                    cls.extend([[-p, x], [-p, -seen], [p, -x, seen]])
                ns = pool.id(("s", t, side, j))
                if seen is None:
                    cls.extend([[-ns, x], [ns, -x]])
                else:
                    cls.extend([[-ns, seen, x], [ns, -seen], [ns, -x]])
                seen = ns
                piv[(side, j)] = p
            cls.append([seen])
        sL = None
        for j in range(M):
            pl, pr = piv[("L", j)], piv[("R", j)]
            if sL is None:
                cls.append([-pr])
            else:
                cls.append([-pr, sL])
            cls.append([-pl, -coef[(t, "R", j)]])
            cls.append([-pr, -coef[(t, "L", j)]])
            nsl = pool.id(("sl", t, j))
            if sL is None:
                cls.extend([[-nsl, pl], [nsl, -pl]])
            else:
                cls.extend([[-nsl, sL, pl], [nsl, -sL], [nsl, -pl]])
            sL = nsl
    return pool, cls, coef, T


def add_pattern(pool, cls, coef, T, n, m, k, p, allowed):
    """p: tuple of input bits; allowed: set of output words (tuples of m bits)."""
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
    g = []
    for t in range(k):
        vals = []
        for side in "LR":
            terms = [coef[(t, side, i)] for i in range(n) if p[i]]
            terms += [and2(coef[(t, side, n + s)], g[s]) for s in range(t)]
            vals.append(xor_list(terms))
        g.append(and2(vals[0], vals[1]))
    outs = []
    for o in range(m):
        terms = [coef[("o", o, 0)]] + [coef[("o", o, 1 + i)] for i in range(n) if p[i]]
        terms += [and2(coef[("o", o, 1 + n + s)], g[s]) for s in range(k)]
        outs.append(xor_list(terms))
    if len(allowed) == 1:
        (w,) = tuple(allowed)
        for o in range(m):
            cls.append([outs[o]] if w[o] else [-outs[o]])
        return outs
    for wi in range(1 << m):
        w = tuple((wi >> (m - 1 - o)) & 1 for o in range(m))
        if w not in allowed:
            cls.append([(-outs[o] if w[o] else outs[o]) for o in range(m)])
    return outs


def evaluate(prog, p):
    vals = list(p)
    for g in prog["gates"]:
        a = 0
        for j in g["left"]:
            a ^= 1 if j == -1 else vals[j]
        b = 0
        for j in g["right"]:
            b ^= 1 if j == -1 else vals[j]
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


def problem(tab, M0s, keep0, exact):
    ms = tab["ms"]
    mg = float(tab["margin"])
    if M0s == "margin":
        M0 = mg - 1
    elif M0s.startswith("m/"):
        M0 = mg / float(M0s[2:])
    else:
        M0 = float(M0s)
    nH, nv, j = int(tab["nH"]), int(tab["nv"]), int(tab["j"])
    m = j + 1
    n = nv + nH
    pats = {}
    for c in range(ms.shape[0]):
        H, vt, wt = int(tab["Hval"][c]), int(tab["vtop"][c]), int(tab["wtoday"][c])
        p = tuple((vt >> i) & 1 for i in range(nv)) + tuple((H >> i) & 1 for i in range(nH))
        al = set()
        for wi in range(1 << m):
            w = tuple((wi >> (m - 1 - o)) & 1 for o in range(m))
            ok = (wi == wt) if exact else (ms[c, wi] >= M0)
            if keep0 and w[0] != ((wt >> (m - 1)) & 1):
                ok = False
            if ok:
                al.add(w)
        pats[p] = al
    names = ["v%d" % (7 - nv + i) for i in range(nv)] + ["H%d" % i for i in range(nH)]
    return names, pats, n, m, M0


# ------------------------------------------------------------------ the ROW FAMILY (in place by construction)
S1 = [0, 0, 0, 2, 3, 4, 5, 6]
S2 = [0, 0, 1, 2, 3, 4, 6, 6]
BT = 6


def fam_rows(which, j, shift, taper=()):
    """per row of the last j rows: (q, kk): the row adds / subtracts its operand to the window bits q .. n-2 (n = the
    window's bits, top = the sign that becomes the digit); kk = today's taper (names today's operand literals)."""
    if which == "MAIN":                     # k9f: the last row (b = 1) runs at taper 6 on two wires (R41_XT)
        rows = [[0, BT - b] for b in range(j, 1, -1)] + [[1, 6]]
    else:
        sched = S1 if which == "DIV1" else S2
        k0 = 9 - j
        rows = [[sched[k - 1] + 1 - (shift + k - k0), sched[k - 1]] for k in range(k0, 9)]
    for r in taper:
        rows[r][0] += 1
    return [tuple(x) for x in rows]


def today_lits(kk, cnt, kkmin):
    """today's literals of a row with taper kk and cnt targets, as ctx index lists (ctx 0 = const, 1 + i = v bit
    kkmin + i): cin = v_kk, operand bit i = Vb[kk + 1 + i] (v, then the constant 1 at index 7, then nothing)."""
    cin = [1 + kk - kkmin]
    ops = []
    for i in range(cnt):
        b = kk + 1 + i
        ops.append([1 + b - kkmin] if b <= 6 else ([0] if b == 7 else []))
    return cin, ops


def fam_program(rows, lits, nv, nH, pol=1):
    """the family member with the literals `lits` (per row: (cin, [operand literals]); a literal = list of ctx
    indices: 0 const, 1 + i = input i (a v bit), 1 + nv + r = the digit of sub-map row r) as an XOR-AND program."""
    n_in = nv + nH
    gates = []
    M = [frozenset([nv + i]) for i in range(nH)]
    digits = []

    def form(lit):
        f = frozenset()
        for c in lit:
            if c == 0:
                f ^= frozenset([-1])
            elif c <= nv:
                f ^= frozenset([c - 1])
            else:
                f ^= digits[c - 1 - nv]
        return f
    for (q, kk), (cin, ops) in zip(rows, lits):
        n = len(M)
        su = M[n - 1]
        digits.append(su)
        e = su ^ frozenset([-1]) if pol else su
        c = form(cin)
        for i, pos in enumerate(range(q, n - 1)):
            t, o = M[pos], form(ops[i])
            if pos < n - 2:
                gates.append(dict(left=sorted(t ^ e ^ c), right=sorted(o ^ c)))
                cn = c ^ frozenset([n_in + len(gates) - 1])
            M[pos] = t ^ o ^ c
            if pos < n - 2:
                c = cn
        M = M[:n - 1]
    outs = [sorted(d) for d in digits] + [sorted(M[-1])]
    return dict(gates=gates, outputs=outs, leftover=[sorted(x) for x in M[:-1]])


def fam_solve(tab, pats, names, n, m, M0, opt):
    which, j, shift, kkmin = str(tab["which"]), int(tab["j"]), int(tab["shift"]), int(tab["kkmin"])
    nH, nv = int(tab["nH"]), int(tab["nv"])
    taper = tuple(int(x) for x in opt["taper"].split(",")) if opt.get("taper") else ()
    rows = fam_rows(which, j, shift, taper)
    use_dig = int(opt.get("digits", 0))
    shapes, today = [], []
    nn = nH
    for (q, kk) in rows:
        cnt = nn - 1 - q
        shapes.append(cnt)
        today.append(today_lits(kk, cnt, kkmin))
        nn -= 1
    ands = sum(max(0, c - 1) for c in shapes)
    pol = int(opt.get("pol", 1))
    if not taper:                                   # self-test: today's literals give today's word on every class
        ok_pol = None
        for pl in (1, 0):
            prog = fam_program(rows, today, nv, nH, pl)
            good = True
            for c in range(tab["ms"].shape[0]):
                H, vt, wt = int(tab["Hval"][c]), int(tab["vtop"][c]), int(tab["wtoday"][c])
                p = tuple((vt >> i) & 1 for i in range(nv)) + tuple((H >> i) & 1 for i in range(nH))
                w = evaluate(prog, p)
                if sum(b << (m - 1 - o) for o, b in enumerate(w)) != wt:
                    good = False
                    break
            if good:
                ok_pol = pl
                break
        print("row family: rows (q, kk) %s, targets %s, %d ANDs; self-test with today's literals: %s"
              % (rows, shapes, ands, "PASS (pol %d)" % ok_pol if ok_pol is not None else "FAIL"), flush=True)
        if ok_pol is None:
            return None
        pol = ok_pol
    else:
        print("row family: rows (q, kk) %s, targets %s, %d ANDs (tapered rows %s)" % (rows, shapes, ands, taper),
              flush=True)
    pool = IDPool()
    cls = []
    T = pool.id("T")
    cls.append([T])
    F = -T

    def xor2(a, b):
        if a == F:
            return b
        if b == F:
            return a
        if a == T:
            return -b
        if b == T:
            return -a
        c = pool.id()
        cls.extend([[-a, -b, -c], [a, b, -c], [a, -b, c], [-a, b, c]])
        return c

    def and2(a, b):
        c = pool.id()
        cls.extend([[-c, a], [-c, b], [c, -a, -b]])
        return c

    def maj3(a, b, c_):
        d = pool.id()
        cls.extend([[-a, -b, d], [-a, -c_, d], [-b, -c_, d], [a, b, -d], [a, c_, -d], [b, c_, -d]])
        return d
    nctx = 1 + nv
    coef = {}
    for r, cnt in enumerate(shapes):
        for l in range(cnt + 1):                    # l = 0: cin; l = 1 + i: operand bit i
            for cidx in range(nctx + (r + 1 if use_dig else 0)):
                coef[(r, l, cidx)] = pool.id(("k", r, l, cidx))
    if opt.get("fixops"):                           # only the carry-in literals are free
        for r, cnt in enumerate(shapes):
            for i in range(cnt):
                for cidx in range(nctx + (r + 1 if use_dig else 0)):
                    v = coef[(r, 1 + i, cidx)]
                    cls.append([v] if cidx in today[r][1][i] else [-v])
    cache = {}

    def litval(r, l, vbits, digs):
        key = (r, l, vbits)
        if key not in cache:
            acc = coef[(r, l, 0)]
            for i in range(nv):
                if vbits[i]:
                    acc = xor2(acc, coef[(r, l, 1 + i)])
            cache[key] = acc
        acc = cache[key]
        if use_dig:
            for t, dl in enumerate(digs):
                acc = xor2(acc, and2(coef[(r, l, nctx + t)], dl))
        return acc
    for p, al in pats.items():
        vbits = p[:nv]
        M = [T if b else F for b in p[nv:]]
        digs = []
        for r, ((q, kk), cnt) in enumerate(zip(rows, shapes)):
            nb = len(M)
            su = M[nb - 1]
            digs.append(su)
            e = -su if pol else su
            c = litval(r, 0, vbits, digs)
            for i, pos in enumerate(range(q, nb - 1)):
                t = M[pos]
                o = litval(r, 1 + i, vbits, digs)
                if pos < nb - 2:
                    cn = maj3(xor2(t, e), o, c)
                M[pos] = xor2(xor2(t, o), c)
                if pos < nb - 2:
                    c = cn
            M = M[:nb - 1]
        outs = digs + [M[-1]]
        for wi in range(1 << m):
            w = tuple((wi >> (m - 1 - o)) & 1 for o in range(m))
            if w not in al:
                cls.append([(-outs[o] if w[o] else outs[o]) for o in range(m)])
    solver = opt.get("solver", "cadical153")
    t0 = time.time()
    with Solver(name=solver, bootstrap_with=cls) as s_:
        ok = s_.solve()
        mdl = set(l for l in s_.get_model() if l > 0) if ok else None
    print("k = %d: %s (row family, %d classes, %d clauses, %d outputs, solver %s, %.1f s)"
          % (ands, "SAT" if ok else "UNSAT", len(pats), len(cls), m, solver, time.time() - t0), flush=True)
    if not ok:
        return None
    lits = []
    cname = ["1"] + names[:nv] + ["dig%d" % t for t in range(len(rows))]
    for r, cnt in enumerate(shapes):
        row = []
        for l in range(cnt + 1):
            row.append([cidx for cidx in range(nctx + (r + 1 if use_dig else 0)) if coef[(r, l, cidx)] in mdl])
        lits.append((row[0], row[1:]))
        td = today[r]
        print("  row %d (q %d, taper %d): cin = %s (today %s); operand bits %s (today %s)"
              % (r, rows[r][0], rows[r][1], "^".join(cname[c] for c in row[0]) or "0",
                 "^".join(cname[c] for c in td[0]) or "0",
                 [("^".join(cname[c] for c in o) or "0") for o in row[1:]],
                 [("^".join(cname[c] for c in o) or "0") for o in td[1]]))
    prog = fam_program(rows, lits, nv, nH, pol)
    bad = [p for p, al in pats.items() if evaluate(prog, p) not in al]
    assert not bad, "family program disagrees with the encoding on %d classes" % len(bad)
    prog["inputs"] = names
    prog["family"] = dict(rows=rows, lits=lits, pol=pol, ands=ands, names=cname)
    return prog


def main():
    args = [a for a in sys.argv[1:] if "=" not in a]
    opt = dict(a.split("=", 1) for a in sys.argv[1:] if "=" in a)
    tab = np.load(args[0], allow_pickle=True)
    k = int(args[1])
    names, pats, n, m, M0 = problem(tab, args[2], int(opt.get("keep0", 0)), int(opt.get("exact", 0)))
    empty = sum(1 for a in pats.values() if not a)
    free = sum(1 for a in pats.values() if len(a) > 1)
    tot = sum(len(a) for a in pats.values())
    print("%s last %d rows: %d inputs %s, %d outputs, %d reachable classes (of %d patterns), M0 %.0f: allowed words "
          "per class %.2f, classes with a choice %d, classes with NO allowed word %d"
          % (str(tab["which"]), int(tab["j"]), n, names, m, len(pats), 1 << n, M0, tot / len(pats), free, empty),
          flush=True)
    if empty:
        print("k = %d: UNSAT (a class has no allowed word at this M0)" % k)
        return
    if int(opt.get("fam", 0)):
        prog = fam_solve(tab, pats, names, n, m, M0, opt)
        if prog is not None and "out" in opt:
            json.dump(prog, open(opt["out"], "w"), indent=1)
        return
    t0 = time.time()
    pool, cls, coef, T = encode_gates(n, m, k)
    outs_by_p = {}
    for p, al in pats.items():
        outs_by_p[p] = add_pattern(pool, cls, coef, T, n, m, k, p, al)
    if int(opt.get("inj", 0)):
        # THE WINDOW STAYS ONE-TO-ONE: reachable classes with the same v bits and the same lowest window bit (that
        # wire is the leftover: affine in known wires, cleared by a relation) must get different words; otherwise the
        # tail cannot be a reversible map of its own wires and a wire is held
        nv_ = int(tab["nv"])
        groups = {}
        for p, outs in outs_by_p.items():
            groups.setdefault(p[:nv_] + (p[nv_],), []).append(outs)
        npair = 0
        cap2 = str(tab["which"]) == "MAIN"          # MAIN's window ends with TWO leftover wires: a word may be
        for g in groups.values():                   # taken by two classes of a group, never by three
            eq = {}
            for a_ in range(len(g)):
                for b_ in range(a_ + 1, len(g)):
                    diff = []
                    for o in range(m):
                        x = pool.id()
                        A, Bq = g[a_][o], g[b_][o]
                        cls.extend([[-A, -Bq, -x], [A, Bq, -x], [A, -Bq, x], [-A, Bq, x]])
                        diff.append(x)
                    if not cap2:
                        cls.append(diff)
                    else:
                        e = pool.id()                # e <- the two words are equal
                        cls.append([e] + diff)
                        eq[(a_, b_)] = e
                    npair += 1
            if cap2:
                for a_ in range(len(g)):
                    for b_ in range(a_ + 1, len(g)):
                        for c_ in range(b_ + 1, len(g)):
                            cls.append([-eq[(a_, b_)], -eq[(b_, c_)]])
        print("inj: %d groups (v bits, lowest window bit), %d pairs %s" % (len(groups), npair, "(no word three "
              "times in a group)" if cap2 else "kept apart"), flush=True)
    if "dump" in opt:
        from pysat.formula import CNF
        CNF(from_clauses=cls).to_file(opt["dump"])
        print("dumped %d clauses, %d vars" % (len(cls), pool.top), flush=True)
        if opt.get("solver") == "none":
            return
    solver = opt.get("solver", "cadical153")
    with Solver(name=solver, bootstrap_with=cls) as s:
        ok = s.solve()
        prog = extract(set(l for l in s.get_model() if l > 0), coef, n, m, k) if ok else None
    print("k = %d: %s (%d classes, %d clauses, %d outputs, solver %s, %.1f s)"
          % (k, "SAT" if ok else "UNSAT", len(pats), len(cls), m, solver, time.time() - t0), flush=True)
    if ok:
        nm = lambda j_: names[j_] if j_ < n else "g%d" % (j_ - n)
        bad = [p for p, al in pats.items() if evaluate(prog, p) not in al]
        assert not bad
        prog["inputs"] = names
        for t, g in enumerate(prog["gates"]):
            print("  g%d = (%s) & (%s)" % (t, " ^ ".join(map(nm, g["left"])), " ^ ".join(map(nm, g["right"]))))
        for o, lst in enumerate(prog["outputs"]):
            print("  out%d = %s" % (o, " ^ ".join("1" if j_ == -1 else nm(j_) for j_ in lst)))
        if "out" in opt:
            json.dump(prog, open(opt["out"], "w"), indent=1)


if __name__ == "__main__":
    main()
