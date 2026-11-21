#!/usr/bin/env python3
"""invsqrt_fc2_patch_h1c.py -- FC2 patch "h1c" for FC's g8 (own code): the higher-order stage H1' for any root split
m (nb = 17 - m folded bits of T: m 8 -> 9, m 7 -> 10), with the arithmetic of the R37 designer's integer reference
tools/r37_inv_40_ref.py (read as a model):
  T = y >> 9 (nb + 1 bits), L = y & 511, s = sign, Tf_i = T_i xor s (i < nb).
  A   = lin Tf - sum w_ij (Tf_i xor Tf_j)   half units, lin = 3 * 2^(nb - 9);
        w_ij = 3 * 2^(i+j-9) (i + j >= 9; w_45 = 2), 1 for i + j = 8 (unless in drop8), else 0
  cub = sum g_k P_k(Tf)  units, P_k = AND of folded bits;  h = (A + K2) >> 1 + (1 - 2 s) cub;  U = y + h + off.
Circuit: register R = [Gd, L_0..L_8, c0 (, c1, c2 if nc 3)] in half units (Gd = guard, fresh 0); every word is added /
subtracted into R by FC's adder A.add (narrow block + explicit undo), mod 2^len(R):
  lin      Tf at positions i + nb - 9 and i + nb - 8 (two adds)
  row j    pair parts of (i, j) as in-place literals x_ij = Tf_i xor Tf_j on T[i]: part 2 at i + j - 8,
           part 1 at i + j - 9; merged rows j (option merge): part 2 of row j + part 1 of row j + 1 as ONE word
           (e = Tf_j xor Tf_j+1; {0: e, i + j - 7: e OR x_ij}: R36's merge, OR temps)
  cubic    (1 - 2 s) g P = g (P xor s) - g s: the word {pos(g): P xor s} added; the -15 s units either as a word
           (subtracted) or (cubfold) as s at position 5 of the lowest merged word (-32 s) + two carry-ins s (+2 s)
  const    K2 + 2 off (half units) as carry-ins / borrow-ins ONE
Gd is the XOR of the position-0 literals (cleared by CNOTs). nc 1: C >> 1 = c1 + 2 c2 from FC2's threshold rule
(emit_chi, 17 ANDs) on the folded T before the increment; nc 3: in place. Then unfold, T += C; c1, c2 are
X-measured and the mirror's explicit undo rebuilds them by the rule on T' (repays the debts); c0 is held.
  cfg["h1c"] = {"nb": 10, "cub": [[g, [bits]], ...], "K2": 0, "off": -1, "drop8": [], "nc": 1, "merge": [5, 6, 7],
                "cubfold": 1}
  invsqrt_fc2_patch_h1c.py check SRC N CFG.json 'OVERRIDES'   model vs circuit (word, x, ancillas) on sampled inputs
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

DEF = dict(nb=10, cub=[[4, [9]], [8, [9, 8]], [2, [9, 7]], [1, [9, 6]]], K2=0, off=-1, drop8=[], nc=3, merge=[],
           cubfold=0)


def opts(cfg):
    o = dict(DEF)
    o.update(cfg.get("h1c", {}))
    return o


def weights(o):
    """[(i, j, position)] of every pair part the stage subtracts (half units)."""
    nb = int(o["nb"])
    drop8 = {tuple(p) for p in o.get("drop8", [])}
    out = []
    for j in range(nb):
        for i in range(j):
            if i + j == 8:
                if (i, j) not in drop8:
                    out.append((i, j, 0))
            elif (i, j) == (4, 5):
                out.append((i, j, 1))
            elif i + j >= 9:
                out.append((i, j, i + j - 9))
                out.append((i, j, i + j - 8))
    return out


def active_parts(o):
    """weights(o) minus option drops: "dropparts" [[i, j, pos]] and "mlow" {j: imin} (merged row j starts at pair
    imin: part 2 of (i, j) and part 1 of (i, j + 1) for i < imin are dropped; both sit at position i + j - 8)."""
    drop = {tuple(int(t) for t in q) for q in o.get("dropparts", [])}
    for j, imin in (o.get("mlow") or {}).items():
        j, imin = int(j), int(imin)
        for i in range(imin):
            drop.add((i, j, i + j - 8))
            drop.add((i, j + 1, i + j - 8))
    return [q for q in weights(o) if q not in drop]


def table_prog(k, vals, nout=5):
    """F1's minand format program for the bits of vals[t] (t = k input bits z0..z(k-1), z0 lowest) from the ANF:
    every monomial of degree >= 2 an AND chain (shared prefixes); outputs = XOR of monomials (-1 = constant 1)."""
    n = k
    gates, gid = [], {}

    def prod(ms):                                      # ms: sorted tuple of input indices (len >= 2)
        if ms in gid:
            return gid[ms]
        left = ms[0] if len(ms) == 2 else prod(ms[:-1])
        gates.append({"left": [left], "right": [ms[-1]]})
        gid[ms] = n + len(gates) - 1
        return gid[ms]
    outs = []
    for bit in range(nout):
        anf = [(int(vals[t]) >> bit) & 1 for t in range(1 << k)]
        for b in range(k):
            for z in range(1 << k):
                if z >> b & 1:
                    anf[z] ^= anf[z ^ (1 << b)]
        terms = []
        for z in range(1 << k):
            if not anf[z]:
                continue
            ms = tuple(i for i in range(k) if z >> i & 1)
            if len(ms) == 0:
                terms.append(-1)
            elif len(ms) == 1:
                terms.append(ms[0])
            else:
                terms.append(prod(ms))
        outs.append(terms)
    return {"gates": gates, "outputs": outs}


_TABPROG = {"prog": None}


def emit_table(c, k, vals, zin, outs, uncompute=False):
    """outs ^= the table bits (outs zero, or holding the table when uncompute): F3's six_a programs when cfg h1c
    "cubprog" is "f3" (CUB8 at >= 5 temporaries, else the staged 11-AND form), a named F1-format program, or the
    trie."""
    if _TABPROG["prog"] == "f3":
        import invsqrt_fw_chiprog as CP
        return CP.emit_cub_auto(c, list(zin), list(outs), uncompute=uncompute)
    if isinstance(_TABPROG["prog"], list):             # staged [[bit, given, prog], ..] (R37 re-pricer's order file)
        import invsqrt_fw_chiprog as CP
        st = _TABPROG["prog"]
        k_ = 0
        for b_, given, prog in (reversed(st) if uncompute else st):
            k_ += CP.emit(c, prog, list(zin) + [outs[d] for d in given], [outs[b_]], max(1, c.cap - c.live))
        return k_
    if _TABPROG["prog"] is not None:
        import invsqrt_fw_chiprog as CP
        return CP.emit(c, _TABPROG["prog"], list(zin), list(outs), max(1, c.cap - c.live))
    return emit_table_trie(c, k, vals, zin, outs)


def emit_phase_cond(c, prog, ins, g):
    """(-1)^(g f(ins)), f = output 0 of an F1-format program, g a measured bit: a gate read by later gates is
    computed into a temporary by CCX if g (CANDC: 1/2 expected) and erased by measurement (CANDE: free); a gate read
    only by the output acts as CZ if g between its two factors (free); input / temporary terms as Z if g, the
    constant as NEG if g. Returns the number of conditioned ANDs."""
    n = len(ins)
    used = [False] * len(prog["gates"])
    for gt in prog["gates"]:
        for x in gt["left"] + gt["right"]:
            if x >= n:
                used[x - n] = True
    wire = {i: ins[i] for i in range(n)}

    def form(lst, other):
        cand = [x for x in lst if x not in other and x in wire]
        assert cand, ("no carrier", lst, other)
        car = cand[0]
        rest = [x for x in lst if x != car]
        for x in rest:
            c.cx(wire[x], wire[car])
        return car, rest

    def unform(car, rest):
        for x in reversed(rest):
            c.cx(wire[x], wire[car])

    def factors(j):
        gt = prog["gates"][j]
        cl, rl = form(gt["left"], gt["right"])
        cr, rr = form(gt["right"], gt["left"])
        return (cl, rl), (cr, rr)

    temps, k = [], 0
    for j in range(len(prog["gates"])):
        if not used[j]:
            continue
        (cl, rl), (cr, rr) = factors(j)
        t = c.alloc(1)[0]
        c.cand_compute(wire[cl], wire[cr], t, g)
        unform(cr, rr)
        unform(cl, rl)
        wire[n + j] = t
        temps.append(j)
        k += 1
    for x in prog["outputs"][0]:
        if x == -1:
            c.neg(g)
        elif x < n or used[x - n]:
            c.z(wire[x], g)
        else:
            (cl, rl), (cr, rr) = factors(x - n)
            c.cz(wire[cl], wire[cr], g)
            unform(cr, rr)
            unform(cl, rl)
    for j in reversed(temps):
        (cl, rl), (cr, rr) = factors(j)
        c.cand_erase(wire[cl], wire[cr], wire[n + j], g)
        unform(cr, rr)
        unform(cl, rl)
        c.free_([wire[n + j]])
        del wire[n + j]
    return k


def emit_table_trie(c, k, vals, zin, outs):
    """outs[bit] ^= bit of vals[t] (t on the wires zin, zin[0] lowest): depth-first over the trie of the ANF
    monomials: one AND per distinct prefix (computed, XORed into every output bit whose ANF holds it, its children
    done, erased free); at most k - 1 temporaries live; linear terms and the constant by CNOT / X."""
    nout = len(outs)
    anfs = []
    for bit in range(nout):
        a = [(int(vals[t]) >> bit) & 1 for t in range(1 << k)]
        for b in range(k):
            for z in range(1 << k):
                if z >> b & 1:
                    a[z] ^= a[z ^ (1 << b)]
        anfs.append(a)
    mons = {}
    for bit in range(nout):
        for z in range(1 << k):
            if anfs[bit][z]:
                mons.setdefault(tuple(i for i in range(k) if z >> i & 1), []).append(bit)
    for ms, bits in mons.items():                      # constant and linear terms
        if len(ms) == 0:
            for b_ in bits:
                c.x(outs[b_])
        elif len(ms) == 1:
            for b_ in bits:
                c.cx(zin[ms[0]], outs[b_])
    prefixes = set()
    for ms in mons:
        for L_ in range(2, len(ms) + 1):
            prefixes.add(ms[:L_])
    ands = [0]

    def leaf(ch):                                      # a leaf that would be the 4th live temporary: CCX instead
        return len(ch) >= 5 and not any(ch + (r,) in prefixes for r in range(ch[-1] + 1, k))

    def visit(pref, wire):
        for b_ in mons.get(pref, []):
            c.cx(wire, outs[b_])
        for q in range(pref[-1] + 1, k):
            ch = pref + (q,)
            if ch in prefixes:
                if leaf(ch):                           # a leaf: CCX into each output that holds it (no temporary)
                    for b_ in mons.get(ch, []):
                        c.ccx(wire, zin[q], outs[b_])
                        ands[0] += 1
                    continue
                t = c.and_compute(wire, zin[q])
                ands[0] += 1
                visit(ch, t)
                c.and_erase(wire, zin[q], t)
    for i in range(k):
        for q in range(i + 1, k):
            ch = (i, q)
            if ch in prefixes:
                t = c.and_compute(zin[i], zin[q])
                ands[0] += 1
                visit(ch, t)
                c.and_erase(zin[i], zin[q], t)
    return ands[0]


def form_literal(state, target, avoid):
    """a carrier wire holding the XOR of the folded bits in `target` (a frozenset) from the T wires' current values
    state {wire: (frozenset, neg)}: returns (pre ops, carrier) or None; the carrier is not in `avoid`."""
    ws = [w for w in state if w not in avoid]
    allw = list(state)
    for w in ws:                                       # one wire
        st, ng = state[w]
        if st == target:
            return ([("x", w)] if ng else []), w
    for w in ws:                                       # carrier + one more
        for v in allw:
            if v == w:
                continue
            if state[w][0] ^ state[v][0] == target:
                ng = state[w][1] ^ state[v][1]
                return ([("cx", v, w)] + ([("x", w)] if ng else [])), w
    for w in ws:                                       # carrier + two more
        for i, v in enumerate(allw):
            for u in allw[i + 1:]:
                if w in (v, u):
                    continue
                if state[w][0] ^ state[v][0] ^ state[u][0] == target:
                    ng = state[w][1] ^ state[v][1] ^ state[u][1]
                    return ([("cx", v, w), ("cx", u, w)] + ([("x", w)] if ng else [])), w
    return None


def _do(c, op, t=None):
    """one self-inverse op: ("x", w) or ("cx", a, b); the name "T" stands for the temporary t."""
    if op[0] == "x":
        c.x(t if op[1] == "T" else op[1])
    else:
        c.cx(op[1], t if op[2] == "T" else op[2])


class _Ands:
    """temporary ANDs with recorded erasure (measurement-based, 0 Toffolis), erased in reverse order."""
    def __init__(self, c):
        self.c, self.er = c, []

    def g(self, a, b, pre=(), post=()):
        c = self.c
        for op in pre:
            _do(c, op)
        t = c.and_compute(a, b)
        for op in reversed(pre):
            _do(c, op)
        for op in post:
            _do(c, op, t)

        def er():
            for op in reversed(post):
                _do(c, op, t)
            for op in pre:
                _do(c, op)
            c.and_erase(a, b, t)
            for op in reversed(pre):
                _do(c, op)
        self.er.append(er)
        return t

    def erase(self):
        for er in reversed(self.er):
            er()
        self.er = []


_CHI = {"mode": "fc2"}


def emit_chi(c, z, l0, l1, c0, out1, out2):
    """out1 ^= c1, out2 ^= c2, C >> 1 = c1 + 2 c2 by FC2's threshold rule (exact on all y of m 7 for DEF7, before
    and after T += C; tools/invsqrt_fc2_chi.py): z = the top 4 folded bits of T (z[0] LSB), l1 l0 = L''_8 L''_7,
      Chi = [z >= t0] + [z >= t1] + [z >= t2],  t0 = 6 + 3 c0 + l1 l0,  t1 = 11 + 2 c0 + (c0 ? l1 l0 : l1),
      t2 = c0 ? 16 : 14 + l1;  so c2 = [z >= t1] and c1 = the parity of the three.  17 ANDs, all erased free."""
    if _CHI["mode"] == "f3six":                                      # F3's exact six_a programs (4 + 3 ANDs)
        import invsqrt_fw_chiprog as CP
        CP.emit_staged(c, [z[0], z[1], z[2], z[3], l0, l1, c0], out1, out2, c.cap - c.live, spec="six_a")
        return
    if _CHI["mode"] == "f3six6":                                     # F3's joint six_a program (6 ANDs, <= 5 temps)
        import invsqrt_fw_chiprog as CP
        CP.emit(c, CP.CHI6_K6, [z[0], z[1], z[2], z[3], l0, l1, c0], [out1, out2], c.cap - c.live)
        return
    if _CHI["mode"] == "f3k5":                                       # F3's chooser: CHI6_K5 (5 ANDs at 5 temps),
        import invsqrt_fw_chiprog as CP                              # CHI6_K6, or the staged form
        CP.emit_chi6(c, z, l0, l1, c0, out1, out2, outs_zero=bool(_CHI.get("zero", True)))
        return
    z0, z1, z2, z3 = z
    a = _Ands(c)                                                     # i2 = NOT c0 z3 z2 z1 (z0 OR NOT l1) -> out1
    U1 = a.g(z3, z2)
    U2 = a.g(U1, z1)
    V = a.g(z0, l1, pre=[("x", z0)])                                 # NOT z0 AND l1 = NOT (z0 OR NOT l1)
    W = a.g(U2, V, pre=[("x", V)])
    c.x(c0)
    c.ccx(W, c0, out1)
    c.x(c0)
    a.erase()
    shP = _Ands(c)
    P = shP.g(l1, l0)                                                # p = l1 l0
    shG = _Ands(c)
    G = shG.g(c0, P, pre=[("cx", l1, P)], post=[("cx", l1, "T")])    # q = l1 xor c0 (p xor l1)
    a = _Ands(c)                                                     # i1 = carry of z + (q, c0^q, ~c0~q, 0) + 1
    K1 = a.g(z0, G, pre=[("x", z0), ("x", G)], post=[("x", "T")])    # z0 OR q
    K2 = a.g(z1, G, pre=[("cx", c0, G), ("cx", K1, z1), ("cx", K1, G)], post=[("cx", K1, "T")])  # MAJ(z1, c0^q, K1)
    Rr = a.g(c0, G, pre=[("x", c0), ("x", G)])                      # NOT c0 AND NOT q
    K3 = a.g(z2, Rr, pre=[("cx", K2, z2), ("cx", K2, Rr)], post=[("cx", K2, "T")])
    c.ccx(z3, K3, out2)                                              # out2 = i1 (fresh zero before)
    c.cx(out2, out1)
    a.erase()
    shG.erase()
    a = _Ands(c)                                                     # i0 = carry of z + (~(p^c0), c0~p, c0, ~c0) + 1
    J1 = a.g(z0, P, pre=[("cx", c0, P), ("x", z0)], post=[("x", "T")])   # z0 OR NOT (p xor c0)
    M = a.g(c0, P, pre=[("x", P)])                                   # c0 AND NOT p
    J2 = a.g(z1, M, pre=[("cx", J1, z1), ("cx", J1, M)], post=[("cx", J1, "T")])
    J3 = a.g(z2, c0, pre=[("cx", J2, z2), ("cx", J2, c0)], post=[("cx", J2, "T")])
    c.x(c0)                                                          # out1 ^= MAJ(z3, NOT c0, J3)
    c.cx(J3, z3)
    c.cx(J3, c0)
    c.ccx(z3, c0, out1)
    c.cx(J3, c0)
    c.cx(J3, z3)
    c.x(c0)
    c.cx(J3, out1)
    a.erase()
    shP.erase()


def _mux_add(F8, c, T, R, s, nb, cubfold, s1=False):
    """FC2 mux (nb 10, after h1b's mux8): lin 1 + lin 2 - part 2 of row nb-1 as ONE add into R (11 wires, half
    units): with f = Tf_9 and Tl = Tf_0..8 the three words sum to f ? 8 Tl + 2 : 4 Tl (mod 2^11), i.e. the word
        {1: f, 2: NOT f AND Tl_0, p = 3..10: f ? Tl_(p-3) : Tl_(p-2)}   (+ cubfold: s at position 0, carry-in s).
    Forward: each mux bit one AND into a temporary made right before FC's narrow ripple reads it, erased after.
    Undo (explicit): Z2 = NOT f Tl_0, Z3 = f ? Tl_0 : Tl_1 (2 ANDs), Tl_1..8 rotated one place up under f
    (7 CSWAP: wire k = f ? Tl_(k-1) : Tl_k), FC's narrow undo generator, rotation back."""
    import invsqrt_fc_narrow as NW
    from invsqrt_fc2_patch_h1b import _fwd_stream
    assert nb == 10 and len(R) == 11
    f = T[9]
    hold = {}

    def plain(w):
        return (lambda: w, lambda: None)

    def mk2():
        c.x(f)
        t = c.and_compute(f, T[0])
        c.x(f)
        hold[2] = t
        return t

    def rl2():
        c.x(f)
        c.and_erase(f, T[0], hold.pop(2))
        c.x(f)

    def mkp(p):                                        # f ? Tl_(p-3) : Tl_(p-2) = b xor f (a xor b)
        def mk():
            a, b = T[p - 3], T[p - 2]
            c.cx(a, b)
            t = c.and_compute(f, b)
            c.cx(a, b)
            c.cx(b, t)
            hold[p] = t
            return t

        def rl():
            a, b = T[p - 3], T[p - 2]
            t = hold.pop(p)
            c.cx(b, t)
            c.cx(a, b)
            c.and_erase(f, b, t)
            c.cx(a, b)
        return mk, rl

    w0 = s if (cubfold and not s1) else None           # cubfold: +s here, the other +s is the cubic add's carry-in
    lits = [plain(w0), plain(f), (mk2, rl2)] + [mkp(p) for p in range(3, 11)]
    cin = None
    if s1:                                             # s1: both +s go to the table add as carry-in s at position 1;
        lits, R = lits[1:], list(R[1:])                # this add starts at position 1 (the guard untouched)
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    debts = _fwd_stream(c, NW, lits, R, NW._lit(cin))
    c.ops.append(("BEND", bid))
    assert not hold
    save = c.ops
    c.ops = []
    Z2, Z3 = c.alloc(2)
    c.x(f)
    c.and_compute(f, T[0], Z2)
    c.x(f)
    c.cx(T[0], T[1])
    c.and_compute(f, T[1], Z3)                         # f (Tl_0 xor Tl_1)
    c.cx(T[0], T[1])
    c.cx(T[1], Z3)                                     # Z3 = Tl_1 xor f (Tl_0 xor Tl_1)
    chain = [(T[k], T[k + 1]) for k in range(7, 0, -1)]
    for a_, b_ in chain:
        c.cswap(f, a_, b_)
    srcu = [w0, f, Z2, Z3] + [T[k] for k in range(2, 9)]
    if s1:
        srcu = srcu[1:]
    n = len(R)
    C = c.cap - c.live
    for w in R:
        c.x(w)
    if NW.plan(max(0, n - 2 - NW.free_prefix(srcu, n, cin)), C)[0] is not None:
        cost = NW.game_add(c, srcu, R, cin, debts, C)
    elif C > NW.cuccaro_need(srcu, R, cin):
        cost = float(NW.hybrid_add(c, srcu, R, cin, debts, C))
    else:
        cost = float(NW.cuccaro_add(c, srcu, R, cin, debts))
    for w in R:
        c.x(w)
    for a_, b_ in reversed(chain):
        c.cswap(f, a_, b_)
    c.cx(T[1], Z3)
    c.cx(T[0], T[1])
    c.and_erase(f, T[1], Z3, free=False)
    c.cx(T[0], T[1])
    c.x(f)
    c.and_erase(f, T[0], Z2, free=False)
    c.x(f)
    c.free_([Z3, Z2])
    c.blocks[bid] = c.ops
    c.ops = save
    F8.__dict__.setdefault("_fc2_mux_log", []).append((C, cost))


def stage(F8, c, U, cfg):
    """U: 27 - m wires, y in two's complement (LSB first). Leaves U = y + h + off; returns c0 (held)."""
    A = F8.A
    o = opts(cfg)
    nb, nc = int(o["nb"]), int(o["nc"])
    merge = sorted(int(j) for j in o.get("merge", []))
    cubfold = bool(o.get("cubfold"))
    s1 = bool(o.get("s1")) and cubfold and bool(o.get("mux")) and bool(o.get("cubtab")) and bool(o.get("cubdebt"))
    assert len(U) == 9 + nb + 1, ("h1c: U has %d wires, nb %d" % (len(U), nb))
    assert nc in (1, 3)
    L, T = list(U[:9]), list(U[9:])
    s = T[nb]
    gdlate = bool(o.get("gdlate")) and s1 and not o.get("gdinit")   # FC2 (fifth holder): the s1 mux add never
    if gdlate:                                         # touches position 0: the guard is allocated after it
        cw = c.alloc(nc)
        Gd = None
    else:
        Gd = c.alloc(1)[0]
        cw = c.alloc(nc)
    R = [Gd] + L + list(cw)
    nR = len(R)
    for i in range(nb):
        c.cx(s, T[i])                                  # T[i] = Tf_i
    par = set()                                        # XOR terms of R[0] (the guard): "Tf<i>", "s" or "one"

    def tog(*vs):
        for v_ in vs:
            par.symmetric_difference_update({v_})

    gi = o.get("gdinit")                               # the guard's start value: one free half unit times a literal
    if gi:                                             # "one", "t<k>" (Tf_k), "n<k>" (NOT Tf_k), "x<i>_<j>", "y<i>_<j>"
        if gi == "one":                                # (NOT x_ij); cleared with the parity at the stage's end
            c.x(Gd)
            tog("one")
        else:
            ks = [int(q) for q in gi[1:].split("_")]
            for k_ in ks:
                c.cx(T[k_], Gd)
                tog("Tf%d" % k_)
            if gi[0] in "ny":
                c.x(Gd)
                tog("one")

    kst = {"k": int(o["K2"]) + 2 * int(o["off"])}      # constant still to place (half units): carry-ins (+1 each on
                                                       # an addition) or borrow-ins (-1 each on a subtraction)

    _CHI["mode"] = o.get("chi", "fc2")
    _TABPROG["prog"] = (o["cubprog"] if o.get("cubprog") == "f3" else
                        json.load(open(os.path.expanduser(o["cubprog"]))) if o.get("cubprog") else None)
    lits = [frozenset([int(b)]) for b in o.get("cin_bits", [])]
    lits += [frozenset([int(a_), int(b_)]) for a_, b_ in o.get("cin_pairs", [])]

    def tstate(row=None, kind=None):
        """current folded-bit content of T[0..nb-1] (as a set) and its negation, during a row's subtraction."""
        st = {T[k]: (frozenset([k]), 0) for k in range(nb)}
        if row is not None:
            j, idx = row
            for i in idx:
                st[T[i]] = (frozenset([i, j]), 1 if kind == "merged" else 0)
            if kind == "merged":
                st[T[j + 1]] = (frozenset([j, j + 1]), 0)
        return st

    def addw(word, sign, cin=None, cvar=None, row=None, kind=None):
        """R += word (sign +1) or R -= word (sign -1); word = {pos: wire or ONE} (positions >= nR drop: mod 2^nR);
        cin: a given carry-in wire (cvar its parity name), else (subtractions) a borrow-in literal of the stage
        (cin_bits / cin_pairs) formed by CNOTs on a T wire that is not a source, else one unit of the constant."""
        src = [word.get(p) for p in range(nR)]
        while src and src[-1] is None:
            src.pop()
        pre = []
        if cin is not None:
            tog(cvar)
        elif sign < 0 and lits and kst["k"] == 0:
            st = tstate(row, kind)
            for L_ in list(lits):
                f_ = form_literal(st, L_, set(w for w in src if isinstance(w, int)))
                if f_ is not None:
                    pre, cin = f_
                    lits.remove(L_)
                    tog(*["Tf%d" % q for q in L_])
                    break
        elif (sign > 0 and kst["k"] > 0) or (sign < 0 and kst["k"] < 0):
            cin = F8.ONE
            kst["k"] -= sign
            tog("one")
        for op in pre:
            _do(c, op)
        if os.environ.get("FC2DBG"):
            print("DBG addw sign %+d row %s kind %s live %d cin %s word %s" % (
                sign, row[0] if row else None, kind, c.live, "lit" if (cin is not None and cin != F8.ONE) else cin,
                "".join("." if w_ is None else ("1" if w_ == F8.ONE else "x") for w_ in src)))
        if sign < 0:
            for t in R:
                c.x(t)
        A.add(c, src, R, cin=cin)
        if sign < 0:
            for t in R:
                c.x(t)
        for op in reversed(pre):
            _do(c, op)

    mux = bool(o.get("mux"))
    if mux:                                            # lin 1 + lin 2 - part 2 of the top row as one mux add
        assert nc == 1
        _mux_add(F8, c, T, R, s, nb, cubfold, s1)      # (cubfold: s at position 0; s1: none here)
        if cubfold and not s1:
            tog("s")
    else:
        lin1 = {i + nb - 9: T[i] for i in range(nb) if 0 <= i + nb - 9 < nR}
        lin2 = {i + nb - 8: T[i] for i in range(nb) if 0 <= i + nb - 8 < nR}
        if 0 in lin1:
            tog("Tf%d" % (9 - nb))
        if cubfold:                                    # cubfold: +2 s half units as two carry-ins s
            addw(lin1, +1, s, "s")
        else:
            addw(lin1, +1)
        if 0 in lin2:
            tog("Tf%d" % (8 - nb))
        if cubfold:
            addw(lin2, +1, s, "s")
        else:
            addw(lin2, +1)
    if gdlate:                                         # FC2: the guard comes to life for the first row word
        Gd = c.alloc(1)[0]
        R[0] = Gd

    # pair rows, plain words: part 2 of row j unless j is merged; part 1 of row j unless j - 1 is merged
    parts = active_parts(o)
    mlow = {int(k_): int(v_) for k_, v_ in (o.get("mlow") or {}).items()}
    for j in range(nb - 1, 0, -1):
        mine = [(i, p) for (i, jj, p) in parts if jj == j and not ((p - i == j - 8 and j in merge) or
                                                                  (p - i == j - 9 and (j - 1) in merge) or
                                                                  (mux and j == nb - 1 and p - i == j - 8))]
        if not mine:
            continue
        idx = sorted({i for i, _ in mine})
        for i in idx:
            c.cx(T[j], T[i])                           # T[i] = x_ij
        groups = {}
        for i, p in mine:
            groups.setdefault(p - i, {})[p] = T[i]     # one word per offset (p - i): part 2 / part 1 / special
        for off_, word in sorted(groups.items(), reverse=True):
            if 0 in word:
                i0 = [i for i, p in mine if p == 0 and p - i == off_][0]
                tog("Tf%d" % i0, "Tf%d" % j)
            addw(word, -1, row=(j, idx), kind="plain")
        for i in idx:
            c.cx(T[j], T[i])

    # merged rows j: part 2 of row j + part 1 of row j + 1 = (e ? ones : 2 X_j), e = Tf_j xor Tf_j+1 (R36's merge)
    drop8 = {tuple(p) for p in o.get("drop8", [])}
    for j in sorted(merge, reverse=True):
        assert 5 <= j <= nb - 2 and not any(q[1] == j for q in drop8)
        idx = [i for i in range(j) if i + j >= 8 and i >= mlow.get(j, 0)]
        p0 = idx[0] + j - 8                            # position of e (0 unless mlow drops the low pairs)
        for i in idx:
            c.cx(T[j], T[i])                           # T[i] = x_ij
        c.cx(T[j], T[j + 1])                           # T[j + 1] = e_j
        e = T[j + 1]
        c.x(e)
        ors = []
        for i in idx:
            c.x(T[i])
            t = c.and_compute(e, T[i])
            c.x(t)                                     # t = e OR x_ij
            ors.append((i, t))
        c.x(e)
        word = {p0: e}
        for i, t in ors:
            word[i + j - 7] = t
        if p0 == 0:
            tog("Tf%d" % j, "Tf%d" % (j + 1))
        if cubfold and j == min(merge):                # cubfold: -2^(top+1) s (half units) as s at that position
            ps = 6 if o.get("cubtab") else 5
            assert ps not in word
            word[ps] = s
        addw(word, -1, row=(j, idx), kind="merged")
        c.x(e)
        for i, t in reversed(ors):
            c.x(t)
            c.and_erase(e, T[i], t)
            c.x(T[i])
        c.x(e)
        c.cx(T[j], T[j + 1])
        for i in idx:
            c.cx(T[j], T[i])

    # cubic: g (P xor s) added at pos(g); g s subtracted at pos(g) (or, cubfold, -32 s + 2 s placed above)
    cub = o.get("cub", [])
    if cub:
        temps, wp, wm = [], {}, {}
        order = sorted(cub, key=lambda gb: -len(gb[1]))   # AND pieces first (they read the folded bits), then
        for g, bits in order:                             # the single-bit pieces (unfolded in place)
            g = int(g)
            p = (2 * g).bit_length() - 1
            assert (1 << p) == 2 * g and p not in wp, ("cubic weights must be distinct powers of two", g)
            if len(bits) == 1:
                b = T[int(bits[0])]
                c.cx(s, b)                             # Tf_b xor s = the unfolded bit
                wp[p] = b
                temps.append(("lin", b))
            else:
                assert len(bits) == 2
                t = c.and_compute(T[int(bits[0])], T[int(bits[1])])
                c.cx(s, t)
                wp[p] = t
                temps.append(("and", t, int(bits[0]), int(bits[1])))
            wm[p] = s
        assert 0 not in wp
        if cubfold:
            assert sorted(wp) == [1, 2, 3, 4] and merge, "cubfold needs the pieces at 1, 2, 3, 4 (sum 15) and a merge"
        if cubfold and mux:
            addw(wp, +1, s, "s")                       # the second +s (half units) of cubfold
        else:
            addw(wp, +1)
        if not cubfold:
            addw(wm, -1)
        for tk in reversed(temps):
            if tk[0] == "lin":
                c.cx(s, tk[1])
            else:
                c.cx(s, tk[1])
                c.and_erase(T[tk[2]], T[tk[3]], tk[1])

    if o.get("cubtab"):                                # cubic as a table of the top k folded bits (signed by s)
        kt, vals = o["cubtab"]
        kt = int(kt)
        nbit = max(int(v).bit_length() for v in vals)
        assert cubfold and mux and nbit == 5, "the table form expects 5 value bits, cubfold and mux"
        zin = [T[nb - kt + i] for i in range(kt)]
        gdfree = bool(o.get("gdearly")) and s1 and bool(o.get("cubdebt"))
        if gdfree:                                     # FC2 (fifth holder): nothing after this adds at position 0 (s1),
            for v_ in sorted(par):                     # so the guard (= the XOR of par) is cleared and freed now: the
                if v_ == "one":                        # table block runs with one more clean wire
                    c.x(Gd)
                elif v_ == "s":
                    c.cx(s, Gd)
                else:
                    c.cx(T[int(v_[2:])], Gd)
            c.free_([Gd])
        bw = c.alloc(nbit)
        if o.get("cubdebt"):
            # forward: table, (b xor s), FC's narrow ripple (carries measured: debts), (b xor s) back, the table bits
            # X-measured (debts). Explicit undo: table again, debts repaid, FC's narrow undo generator, uncompute.
            import invsqrt_fc_narrow as NW
            from invsqrt_fc2_patch_h1b import _fwd_stream
            bid = len(c.blocks)
            c.blocks[bid] = None
            c.ops.append(("BBEG", bid))
            _dbg = os.environ.get("FC2DBG")
            if _dbg:
                print("DBG table fwd: live with bw %d cap %d par %s" % (c.live, c.cap, sorted(par)))
            emit_table(c, kt, vals, zin, bw)
            for w in bw:
                c.cx(s, w)
            if s1:                                     # s1: no guard position; carry-in s at position 1 (= both +s)
                srcw = list(bw) + [None] * (nR - 1 - nbit)
                Rt = list(R[1:])
            else:
                srcw = [None] + list(bw) + [None] * (nR - 1 - nbit)
                Rt = list(R)
                tog("s")
            lits_ = [((lambda w_=w_: w_), (lambda: None)) for w_ in srcw]
            debts = _fwd_stream(c, NW, lits_, Rt, NW._lit(s))
            for w in bw:
                c.cx(s, w)
            mbits = [c.hmr(w) for w in bw]
            c.free_(bw)
            c.ops.append(("BEND", bid))
            save = c.ops
            c.ops = []
            bw2 = c.alloc(nbit)
            assert list(bw2) == list(bw), ("table wires", bw2, bw)
            emit_table(c, kt, vals, zin, bw)
            for w, mb in zip(bw, mbits):
                c.z(w, mb)
            for w in bw:
                c.cx(s, w)
            Cc = c.cap - c.live
            if _dbg:
                print("DBG table undo: live %d Cc %d nR %d" % (c.live, Cc, nR))
            nt = len(Rt)
            for w in Rt:
                c.x(w)
            if NW.plan(max(0, nt - 2 - NW.free_prefix(srcw, nt, s)), Cc)[0] is not None:
                NW.game_add(c, srcw, Rt, s, debts, Cc)
            elif s1 and Cc >= 2:                       # F3's mixed form: Cuccaro on the table positions (carries in
                import invsqrt_fw_c1add as CA          # the table wires), the R32 game on the constant top part
                if CA.NW is None:
                    CA.NW = NW
                CA.mixed_add(c, srcw, Rt, s, debts, Cc, nbit - 1, int(os.environ.get("C1CCZ", "1")))
            elif Cc > NW.cuccaro_need(srcw, Rt, s):
                NW.hybrid_add(c, srcw, Rt, s, debts, Cc)
            else:
                NW.cuccaro_add(c, srcw, Rt, s, debts)
            for w in Rt:
                c.x(w)
            for w in bw:
                c.cx(s, w)
            if o.get("cubphase") == "staged" and isinstance(_TABPROG["prog"], list):
                for b_, given, prog in reversed(_TABPROG["prog"]):   # staged measured clear: the last-computed
                    g_ = c.hmr(bw[b_])                 # bit first, its phase repaid by its own stage program
                    emit_phase_cond(c, prog, list(zin) + [bw[d] for d in given], g_)   # (it reads the bits
                c.free_(bw)                            # still live); 1/2 per AND read by later gates
            elif o.get("cubphase"):                    # clear the table by measurement; repay each bit's phase by
                mb2 = [c.hmr(w) for w in bw]           # its own program, conditioned on the outcome (1/2 per AND
                c.free_(bw)                            # read by later gates; output-only ANDs act as CZ)
                for k_ in range(nbit):
                    pr = json.load(open(os.path.expanduser(o["cubphase"] % k_)))
                    emit_phase_cond(c, pr, zin, mb2[k_])
            else:
                emit_table(c, kt, vals, zin, bw, uncompute=True)
                c.free_(bw)
            c.blocks[bid] = c.ops
            c.ops = save
        else:
            emit_table(c, kt, vals, zin, bw)
            for w in bw:
                c.cx(s, w)                             # (b xor s): the signed table, -62 s placed above
            addw({1 + k_: bw[k_] for k_ in range(nbit)}, +1, s, "s")
            for w in bw:
                c.cx(s, w)
            emit_table(c, kt, vals, zin, bw, uncompute=True)   # uncompute
            c.free_(bw)

    assert not lits, ("borrow-in literals not placed", lits)
    assert kst["k"] == 0, ("constant not placed", kst)
    gd_done = bool(o.get("gdearly")) and s1 and bool(o.get("cubdebt")) and bool(o.get("cubtab"))
    if not gd_done:
        for v_ in sorted(par):                         # the guard bit is the XOR of par: clear it
            if v_ == "one":
                c.x(Gd)
            elif v_ == "s":
                c.cx(s, Gd)
            else:
                c.cx(T[int(v_[2:])], Gd)
    F8.__dict__["_fc2_h1c_dbg"] = dict(Gd=Gd, R=list(R), T=list(T), L=list(L), op=len(c.ops))
    if not gd_done:
        c.free_([Gd])
    c0 = cw[0]
    if nc == 1:                                        # C >> 1 from the threshold rule on the folded T (before T += C)
        c1, c2 = c.alloc(2)
        if o.get("chiphase"):                          # its mirror is a measured clear: c1 (F3's CHI6_C1G reads c2),
            import invsqrt_fw_chiprog as CP            # then c2 (CHI6_C2), phases repaid if the outcome is 1
            ins7 = [T[nb - 4], T[nb - 3], T[nb - 2], T[nb - 1], L[7], L[8], c0]
            bidc = len(c.blocks)
            c.blocks[bidc] = None
            c.ops.append(("BBEG", bidc))
            emit_chi(c, T[nb - 4:nb], L[7], L[8], c0, c1, c2)
            c.ops.append(("BEND", bidc))
            save_ = c.ops
            c.ops = []
            g1 = c.hmr(c1)
            emit_phase_cond(c, CP.CHI6_C1G, ins7 + [c2], g1)
            g2 = c.hmr(c2)
            emit_phase_cond(c, CP.CHI6_C2, ins7, g2)
            c.blocks[bidc] = c.ops
            c.ops = save_
        else:
            emit_chi(c, T[nb - 4:nb], L[7], L[8], c0, c1, c2)
    else:
        c1, c2 = cw[1], cw[2]
    for i in range(nb):
        c.cx(s, T[i])                                  # unfold
    A.add(c, [c0, c1, c2], T)                          # T += C (mod 2^(nb + 1))
    # c1, c2 erased by X-measurement; the mirror's explicit undo rebuilds them from T' (threshold rule after the
    # increment, T' folded by its own sign) and repays the phase debts
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    m1 = c.hmr(c1)
    m2 = c.hmr(c2)
    c.free_([c1, c2])
    c.ops.append(("BEND", bid))
    save = c.ops
    c.ops = []
    r1, r2 = c.alloc(2)
    assert (r1, r2) == (c1, c2), ("explicit undo needs the same wires", (r1, r2), (c1, c2))
    s2 = T[nb]
    for i in range(nb - 4, nb):
        c.cx(s2, T[i])
    emit_chi(c, T[nb - 4:nb], L[7], L[8], c0, c1, c2)
    for i in range(nb - 4, nb):
        c.cx(s2, T[i])
    c.z(c1, m1)
    c.z(c2, m2)
    c.blocks[bid] = c.ops
    c.ops = save
    c.free_([c1, c2])
    return c0


def apply(F8, cfg):
    done = F8.__dict__.setdefault("_fc2_done", set())
    if "h1c" in done:
        return
    orig = F8.st_main

    def st_main_h1c(c, U, v, cfg_):
        c.tag("H1")
        c0 = stage(F8, c, U, cfg_)
        F8.__dict__["_fc2_h1c_held"] = [c0]
        dp = (cfg_.get("h1c") or {}).get("dump")
        if dp is not None:                             # FC2 (fifth holder): R36's DUMP on route 2: the held bit c0
            cfg_["_pre_lz"] = {int(dp): c0}            # stands in for MAIN's low zero at position dp (= m: half a
        c.tag("MAIN")                                  # unit of U), so it needs no wire of its own from row dp on
        return orig(c, U, v, cfg_)

    F8.st_main = st_main_h1c
    done.add("h1c")


def model_U(y, o):
    """the stage's integer function (= r37_inv_40_ref.stage_ref + off), here for any nb."""
    sys.path.insert(0, HERE)
    import r37_inv_40_ref as R37
    if not o.get("dropparts") and not o.get("mlow"):
        st = dict(nb=int(o["nb"]), cub=o.get("cub", []), K2=int(o["K2"]), drop8=o.get("drop8", []))
        for k_ in ("cubtab", "cin_pairs", "cin_bits"):
            if o.get(k_):
                st[k_] = o[k_]
        U, C = R37.stage_ref(y, st)
        return U + int(o["off"]), C
    return h_circ(y, o)


def h_circ(y, o):
    """the circuit's integer function for any option set: U = y + h + off, C = (L + h) >> 9."""
    nb = int(o["nb"])
    y = np.asarray(y, dtype=np.int64)
    s = (y < 0).astype(np.int64)
    T = y >> 9
    Tf = np.where(s == 1, ~T, T) & ((1 << nb) - 1)
    A = (3 << (nb - 9)) * Tf
    for i, j, p in active_parts(o):
        A = A - ((((Tf >> i) ^ (Tf >> j)) & 1) << p)
    cub = np.zeros_like(Tf)
    for g, bits in o.get("cub", []):
        pr = np.ones_like(Tf)
        for b in bits:
            pr = pr & ((Tf >> int(b)) & 1)
        cub = cub + int(g) * pr
    h = ((A + int(o["K2"])) >> 1) + np.where(s == 1, -cub, cub)
    return y + h + int(o["off"]), ((y & 511) + h) >> 9


def check(src, nsamp, cfgfile, over, seed=11):
    """model vs circuit on sampled inputs + edges: word = FC's ref_main on the model's U; x restored; ancillas 0."""
    src = os.path.abspath(os.path.expanduser(src))
    sys.path.insert(0, src)
    import invsqrt_fw_fcpolrun as W
    F8, K, _ = W.load(src)
    import invsqrt_fc_gen as G
    cfg_in = json.load(open(os.path.expanduser(cfgfile)))
    for k_, v_ in over.items():
        if v_ is None:
            cfg_in.pop(k_, None)
        else:
            cfg_in[k_] = v_
    cfg = W.cfg_for(F8, src, cfg_in)
    o = opts(cfg)
    rng = np.random.default_rng(seed)
    if nsamp <= 0:                                     # every input
        x = np.arange(1, F8.MAX_X + 1, dtype=np.int64)
    else:
        x = np.unique(np.concatenate([rng.integers(1, F8.MAX_X + 1, nsamp), np.arange(1, 300),
                                      np.arange(F8.MAX_X - 300, F8.MAX_X + 1)]))
    c, regs = F8.build_full(dict(cfg))
    sim = K.Sim(len(x), c.next + 2)
    sim.set_word(regs["X"], x)
    nf = regs["nfwd"]
    sim.run(c.ops, 0, nf)
    m = cfg["m"]
    n, k, Q, d, delta = F8.ref_chain(x, m)
    _, q1 = G.ref_div(delta.astype(np.int64) << cfg["F1"], d, F8.a_top(cfg), F8.sched1_of(cfg))
    _, P = G.ref_div(q1 + cfg["c2"] * d, d, F8.a_top(cfg), cfg["sched2"])
    y = -((P + 1) >> 1)
    U, C = model_U(y, o)
    Y = F8.ref_main(U, d, m, cfg["bt"])
    Yc = sim.get_word(regs["word"])
    bad = np.flatnonzero(Yc != Y)
    if len(bad):
        T_ = y >> 9
        print("BAD", len(bad), json.dumps([dict(x=int(x[i]), y=int(y[i]), U=int(U[i]), C=int(C[i]), T=int(T_[i]),
                                               L=int(y[i] & 511), dY=int(Yc[i]) - int(Y[i])) for i in bad[:20]]),
              flush=True)
    sim.run(c.ops, nf, None)
    others = [w for w in range(c.next) if w not in regs["X"]]
    return dict(n=len(x), Y_bad=int((Yc != Y).sum()), x_ok=bool(np.array_equal(sim.get_word(regs["X"]), x)),
                anc_ok=bool(not sim.q[others].any()), exp=c.toffolis_exp(), static=c.toffolis(), peak=c.peak,
                stages=c.stage_counts(), stage_peaks=c.stage_peaks, U_range=[int(U.min()), int(U.max())],
                C=[int(C.min()), int(C.max())])


def prof(src, cfgfile, over):
    """H1' blocks: forward / explicit-undo Toffolis of every narrow or explicit block, and the plain ops between."""
    src = os.path.abspath(os.path.expanduser(src))
    sys.path.insert(0, src)
    import invsqrt_fw_fcpolrun as W
    F8, K, _ = W.load(src)
    cfg_in = json.load(open(os.path.expanduser(cfgfile)))
    cfg_in.update(over)
    cfg = W.cfg_for(F8, src, cfg_in)
    c, regs = F8.build_full(dict(cfg))
    marks = c.tags + [(len(c.ops), "_end", c.live)]
    for (i0, nm, lv), (i1, _, _) in zip(marks[:-1], marks[1:]):
        if nm != "H1":
            continue
        j, plain, rows = i0, [], []
        while j < i1:
            if c.ops[j][0] == "BBEG":
                bid = c.ops[j][1]
                k = j
                while not (c.ops[k][0] == "BEND" and c.ops[k][1] == bid):
                    k += 1
                rows.append((c.toffolis_exp(plain), c.toffolis_exp(c.ops[j:k]), c.toffolis_exp(c.blocks[bid] or [])))
                plain = []
                j = k + 1
                continue
            plain.append(c.ops[j])
            j += 1
        tf = tu = 0.0
        for pl, f, u in rows:
            print("  plain-before %5.1f | block fwd %5.1f undo %5.1f" % (pl, f, u))
            tf += pl + f
            tu += pl + u
        print("  plain-after %5.1f" % c.toffolis_exp(plain))
        print("H1 total fwd %.1f undo %.1f (plain ops mirrored at equal cost)" % (tf + c.toffolis_exp(plain), tu + c.toffolis_exp(plain)))


if __name__ == "__main__":
    if sys.argv[1] == "check":
        print(json.dumps(check(sys.argv[2], int(sys.argv[3]), sys.argv[4], json.loads(sys.argv[5]) if len(sys.argv) > 5 else {})))
    elif sys.argv[1] == "prof":
        prof(sys.argv[2], sys.argv[3], json.loads(sys.argv[4]))
