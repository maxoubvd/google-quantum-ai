#!/usr/bin/env python3
"""invsqrt_fc2_patch_h1b.py -- FC2 patch "h1b" for FC's g8: the higher-order stage of family C in the form built and
verified by the R36 designer (tools/r36_inv_corr_patch.py, variant v00: sign-folded top Tf, 1.5 Tf as two plain
adds with the free cubic T8u - s as carry/borrow-ins, row pairs merged through e = Tf_j xor Tf_j+1 and OR gates,
c0 held to the mirror) -- copied here (the designer's file is read-only to builders) with FC2's options:
  cfg["h1b"] = {"jtop": 7, "drops": [[i, j, part], ...], "delta": d, "form": "hold" | "dump"}
  drops (FC2's t'-notation of tools/invsqrt_fc2_h1m2.py: part 1 = the half of the factor 3 at i + j - 10):
    [0, 9, 1]          bit 0 of the second linear add (Tf at half units)
    [8 - j, j + 1, 1]  for a merged pair j (5 <= j <= jtop): the half-unit part at position 0; the merged word
                       becomes {0: x_(8-j, j), 1: e, 2 + k: e OR x_k} (one OR gate less)
  delta: constant in half units on top of the form's own (K = 1533, or 1532 with [0, 9, 1] dropped, in FC2's
    model): +1 as the second linear add's carry-in, -1 per borrow-in (cin ONE) on subtractions without one.
Exact arithmetic = tools/invsqrt_fc2_h1m2.stage_word(y, g 1, koff, ctau [256], drop_parts = drops + [[4,5,1]])
with K = 1533 (- 1 if [0,9,1] dropped) + delta; checked by `check` below on sampled inputs.
  invsqrt_fc2_patch_h1b.py check N 'CFG_OVERRIDES'
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def opts(cfg):
    o = dict(jtop=7, drops=[], delta=0, form="hold")
    o.update(cfg.get("h1b", {}))
    o["drops"] = [tuple(d) for d in o["drops"]]
    return o


def model_K(o):
    if o.get("cinmap"):                                # mux8 + cinmap: A alone is K 1533 with (0,9,1) dropped;
        vals = list(o["cinmap"].values())              # each "one" borrow-in -1, a guard start value +1
        return 1533 - vals.count("one") + (1 if o.get("gdinit") else 0) + int(o["delta"])
    back = {tuple(q) for q in o.get("cinparts", [])}
    return 1533 - (1 if ((0, 9, 1) in o["drops"] and (0, 9, 1) not in back) else 0) + int(o["delta"])


def restored_parts(o):
    """dropped parts that mux8's borrow-ins / guard start value put back (for the model's drop list)."""
    names = list((o.get("cinmap") or {}).values()) + ([o["gdinit"]] if o.get("gdinit") else [])
    out = set()
    for nm in names:
        if nm in ("n09", "p09"):
            out.add((0, 9, 1))
        elif nm != "one":
            a, b = int(nm[1]), int(nm[2])
            out.add((a, b, 1))
    out |= {tuple(q) for q in o.get("cinparts", [])}
    return out


def kept_parts(o):
    """(i, j, part, acc position) of every part the stage subtracts (acc unit = half a unit of U)."""
    drops = set(o["drops"])
    if not o.get("keep451"):
        drops.add((4, 5, 1))
    out = []
    for i in range(10):
        for j in range(i + 1, 10):
            for part, ap in ((1, i + j - 9), (2, i + j - 8)):
                if ap >= 0 and (i, j, part) not in drops:
                    out.append((i, j, part, ap))
    return out


def g0_form(o):
    """D42: G0 (acc bit 1 after the rows, U's unit position) as an exact function of the folded wires z_0..z_8 =
    Tf_0..Tf_8, z_9 = s: (const, linear wire list, products [(la, lb)]) with G0 = const + sum lin + sum
    parity(la) parity(lb) (Dickson form of the quadratic part). Checked on all 1,024 points."""
    K = model_K(o)
    parts = kept_parts(o)
    tt = np.zeros(1024, dtype=np.int64)
    for zi in range(1024):
        sgn = (zi >> 9) & 1
        tf = zi & 511
        T = (-1 - tf) if sgn else tf
        tp = T + 512
        acc = K + 2 * (((T >> 8) & 1) - ((T >> 9) & 1))          # cubic T_8 - T_9 in U units
        for i, j, part, ap in parts:
            acc -= (((tp >> i) ^ (tp >> j)) & 1) << ap
        tt[zi] = (acc >> 1) & 1
    anf = tt.copy()                                    # Moebius transform -> algebraic normal form
    for b in range(10):
        for zi in range(1024):
            if zi >> b & 1:
                anf[zi] ^= anf[zi ^ (1 << b)]
    mons = [zi for zi in range(1024) if anf[zi]]
    assert all(bin(m_).count("1") <= 2 for m_ in mons), "G0 is not quadratic"
    const = int(anf[0])
    lin = set(b for b in range(10) if anf[1 << b])
    Q = set()
    for m_ in mons:
        if bin(m_).count("1") == 2:
            a, b = [q for q in range(10) if m_ >> q & 1]
            Q.add((a, b))
    prods = []
    while Q:
        a, b = min(Q)
        nb = lambda x, y: {c_ for p_ in Q for c_ in p_ if x in p_ and c_ not in (x, y)}
        A_, B_ = nb(a, b), nb(b, a)
        R = {p_ for p_ in Q if a not in p_ and b not in p_}
        for c_ in A_:
            for d_ in B_:
                if c_ == d_:
                    lin ^= {c_}
                else:
                    R ^= {(min(c_, d_), max(c_, d_))}
        prods.append(([a] + sorted(B_), [b] + sorted(A_)))
        Q = R
    for zi in range(1024):                             # exact check of the form
        z = [(zi >> q) & 1 for q in range(10)]
        v = const ^ (sum(z[q] for q in lin) & 1)
        for la, lb in prods:
            v ^= (sum(z[q] for q in la) & 1) & (sum(z[q] for q in lb) & 1)
        assert v == tt[zi], "Dickson form wrong"
    return const, sorted(lin), prods


def _fwd_stream(c, NW, lits, tgt, cin_lit):
    """FC's narrow forward (tools/invsqrt_fc_narrow.fwd, copied: FC's file is read-only) with streamed sources:
    lits[i] = (make, release): make() returns the wire holding source bit i, called right before step i reads it;
    release() erases it right after step i (the top step reads bits n - 2 and n - 1). cin_lit: a kit literal
    ("w", wire, flip). Returns the debts {carry index: HMR bit} for FC's undo generators."""
    n = len(tgt)
    debts = {}
    cur = cin_lit
    cur_fresh = False
    deferred = None
    for i in range(n - 1):
        a = lits[i][0]()
        top_a = lits[n - 1][0]() if i == n - 2 else None
        bw = tgt[i]
        f = NW._form_b(a, bw, cur)
        nxt, nfresh = None, False
        if i < n - 2:
            assert f[0] == "and"
            _, xw, yw, ops = f
            NW._emit(c, ops)
            t = c.alloc(1)[0]
            c.and_compute(xw, yw, t)
            NW._emit(c, ops[::-1])
            NW._xor_lit(c, cur, t)
            nxt, nfresh = ("w", t, 0), True
        else:
            top = tgt[n - 1]
            assert f[0] == "and"
            _, xw, yw, ops = f
            NW._emit(c, ops)
            c.ccx(xw, yw, top)
            NW._emit(c, ops[::-1])
            NW._xor_lit(c, cur, top)
            NW._xor_src(c, top_a, top)
        NW._xor_src(c, a, bw)
        NW._xor_lit(c, cur, bw)
        if cur_fresh:
            debts[i] = c.hmr(cur[1])
            c.free_([cur[1]])
        if i == n - 2:
            lits[n - 1][1]()
        lits[i][1]()
        cur, cur_fresh = nxt, nfresh
    return debts


def _mux8_add(F8, c, T, R):
    """FC2 mux8 (17:45Z): lin 1 + lin 2 + row 8 + both sign borrow-ins as ONE add into R (11 wires, half units):
    with f = Tf_8, Tl = Tf_0..7, s = T[9] the three words sum to
        Tl_0 + NOT f + 2 f Tl_0 + 4 (f ? NOT s : Tl_1) + sum_p=3..8 2^p (f ? Tl_(p-2) : Tl_(p-1))
        + 2^9 f NOT Tl_7 + 2^10 f Tl_7  - 1            (exact on all 1,024 (f, s, Tl); the -1 = one sub's ONE)
    (the f = 1 branch is Tl shifted one place further: 2 Tl against 4 Tl + const). Forward: each mux bit one AND
    into a temporary made right before the narrow ripple reads it and erased after (HMR + CZ): 9 + 10.
    Undo (explicit): the mux word in place by a rotation of (NOT s, Tl_1..7, E) one place up under f (8 CSWAP),
    Z_9 = f + E, Z_1 = f Tl_0 (1 AND), FC's narrow undo generator, rotation back (8 CSWAP)."""
    import invsqrt_fc_narrow as NW
    f, s = T[8], T[9]
    hold = {}

    def plain(w):
        return (lambda: w, lambda: None)

    def mk1():
        t = c.and_compute(f, T[0])
        hold[1] = t
        return t

    def rl1():
        c.and_erase(f, T[0], hold.pop(1))

    def mk2():                                         # f ? NOT s : Tl_1 = Tl_1 + f (Tl_1 + NOT s)
        c.cx(s, T[1]); c.x(T[1])
        t = c.and_compute(f, T[1])
        c.x(T[1]); c.cx(s, T[1])
        c.cx(T[1], t)
        hold[2] = t
        return t

    def rl2():
        t = hold.pop(2)
        c.cx(T[1], t)
        c.cx(s, T[1]); c.x(T[1])
        c.and_erase(f, T[1], t)
        c.x(T[1]); c.cx(s, T[1])

    def mkp(p):                                        # f ? Tl_(p-2) : Tl_(p-1)
        def mk():
            a, b = T[p - 2], T[p - 1]
            c.cx(a, b)
            t = c.and_compute(f, b)
            c.cx(a, b)
            c.cx(b, t)
            hold[p] = t
            return t

        def rl():
            a, b = T[p - 2], T[p - 1]
            t = hold.pop(p)
            c.cx(b, t)
            c.cx(a, b)
            c.and_erase(f, b, t)
            c.cx(a, b)
        return mk, rl

    def mk9():                                         # Z_9 = f NOT Tl_7 = f + g, Z_10 = g = f Tl_7
        g = c.and_compute(f, T[7])
        z9 = c.alloc(1)[0]
        c.cx(f, z9); c.cx(g, z9)
        hold[9], hold[10] = z9, g
        return z9

    def rl9():
        z9, g = hold.pop(9), hold.pop(10)
        c.cx(g, z9); c.cx(f, z9)
        c.free_([z9])
        c.and_erase(f, T[7], g)

    lits = [plain(T[0]), (mk1, rl1), (mk2, rl2)] + [mkp(p) for p in range(3, 9)] + \
        [(mk9, rl9), (lambda: hold[10], lambda: None)]
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    debts = _fwd_stream(c, NW, lits, R, ("w", f, 1))
    c.ops.append(("BEND", bid))
    assert not hold
    save = c.ops
    c.ops = []
    cyc = int(F8.__dict__.get("_fc2_mux8_form", 1)) == 2
    if cyc:
        # form 2: Tl_1..7 rotated cyclically under f (6 CSWAP): T1w = f ? Tl_7 : Tl_1, T_kw = Z_(k+1) (k = 2..7);
        # Z_2 = T1w + f (T1w + NOT s), Z_10 = g = f T1w, Z_9 = f + g, Z_1 = f Tl_0: 3 ANDs, 4 temporaries
        Z1, Z2, G, E9 = c.alloc(4)
        chain = [(T[k], T[k + 1]) for k in range(6, 0, -1)]
        for a, b in chain:
            c.cswap(f, a, b)
        c.cx(s, T[1]); c.x(T[1])
        c.and_compute(f, T[1], Z2)
        c.x(T[1]); c.cx(s, T[1])
        c.cx(T[1], Z2)
        c.and_compute(f, T[1], G)
        c.cx(f, E9); c.cx(G, E9)
        c.and_compute(f, T[0], Z1)
        c.x(f)
        srcu = [T[0], Z1, Z2] + [T[k] for k in range(2, 8)] + [E9, G]
    else:
        E, E9, Z1 = c.alloc(3)
        chain = [(T[7], E)] + [(T[k], T[k + 1]) for k in range(6, 0, -1)] + [(s, T[1])]
        c.x(s)
        for a, b in chain:
            c.cswap(f, a, b)
        c.cx(f, E9); c.cx(E, E9)
        c.and_compute(f, T[0], Z1)
        c.x(f)                                             # cin = NOT f
        srcu = [T[0], Z1] + [T[k] for k in range(1, 8)] + [E9, E]
    n = len(R)
    C = c.cap - c.live
    for w in R:
        c.x(w)
    if NW.plan(max(0, n - 2 - NW.free_prefix(srcu, n, f)), C)[0] is not None:
        cost = NW.game_add(c, srcu, R, f, debts, C)
    elif C > NW.cuccaro_need(srcu, R, f):
        cost = float(NW.hybrid_add(c, srcu, R, f, debts, C))
    else:
        cost = float(NW.cuccaro_add(c, srcu, R, f, debts))
    for w in R:
        c.x(w)
    c.x(f)
    c.and_erase(f, T[0], Z1, free=False)
    if cyc:
        c.cx(G, E9); c.cx(f, E9)
        c.and_erase(f, T[1], G, free=False)
        c.cx(T[1], Z2)
        c.cx(s, T[1]); c.x(T[1])
        c.and_erase(f, T[1], Z2, free=False)
        c.x(T[1]); c.cx(s, T[1])
        for a, b in reversed(chain):
            c.cswap(f, a, b)
        c.free_([E9, G, Z2, Z1])
    else:
        c.cx(E, E9); c.cx(f, E9)
        for a, b in reversed(chain):
            c.cswap(f, a, b)
        c.x(s)
        c.free_([Z1, E9, E])
    c.blocks[bid] = c.ops
    c.ops = save
    F8.__dict__.setdefault("_fc2_mux8_log", []).append((C, cost))


def stage(F8, c, U, cfg):
    """U: 19 wires, y in two's complement (LSB first). Leaves U = y + h; returns the wire c0 (held or dumped)."""
    A = F8.A
    o = opts(cfg)
    jtop = int(o["jtop"])
    drops = set(o["drops"])
    for d in drops:
        assert d == (0, 9, 1) or (d[2] == 1 and d[0] + d[1] == 9 and 6 <= d[1] <= jtop + 1), ("unsupported drop", d)
    d42 = bool(o.get("d42", 0))
    if d42:
        # D42 (FC2, 16:50Z): U holds y2 (18 wires: DIV2 one tail row shorter); U's unit position is a guard G0
        # (fresh zero, erased after the rows by its exact quadratic form); the carry bit c0 becomes U's own LSB
        assert len(U) == 18 and o.get("acc11", 0)
        G0 = c.alloc(1)[0]
        L, T = [G0] + list(U[:8]), list(U[8:])
    else:
        assert len(U) == 19
        L, T = list(U[:9]), list(U[9:])
    s = T[9]
    acc11 = bool(o.get("acc11", 0))                    # FC2: c1 out of the rows' register, computed after them
    if acc11:
        Gd, c0 = c.alloc(2)
        R = [Gd] + L + [c0]
    else:
        Gd, c0, c1 = c.alloc(3)
        R = [Gd] + L + [c0, c1]
    for i in range(9):
        c.cx(s, T[i])                                  # T[i] = Tf_i
    par = set()                                        # base variables of the guard bit (position 0), XOR

    def tog(*vs):
        for v_ in vs:
            par.symmetric_difference_update({v_})
    delta = int(o["delta"])
    mux8 = bool(o.get("mux8"))
    if mux8:                                           # FC2: lin 1 + lin 2 + row 8 + the two s borrow-ins as one add
        assert (0, 9, 1) in drops and delta == 0 and acc11 and not d42 and not o.get("keep451") and jtop == 7
        F8.__dict__["_fc2_mux8_form"] = int(o["mux8"])
        gi = o.get("gdinit")                           # FC2: the guard's start value = one more free bit at position 0
        if gi == "one":
            c.x(Gd)
            tog("one")
        elif gi == "p09":                              # + Tf_0 = - (part (0,9,1) = NOT Tf_0) + 1
            c.cx(T[0], Gd)
            tog("Tf0")
        elif gi:                                       # "nab": NOT x_ab = 1 - x_ab
            a_, b_ = int(gi[1]), int(gi[2])
            c.cx(T[a_], Gd)
            c.cx(T[b_], Gd)
            c.x(Gd)
            tog("Tf%d" % a_, "Tf%d" % b_, "one")
        _mux8_add(F8, c, T, R)
        tog("Tf0", "Tf8", "one")                       # position 0: Tl_0 + NOT f
    # lin 1: Tf << 1 with +T8u (the unfolded bit 8) at the unit position: the cubic's +kappa / -kappa part
    if not mux8:
        _lin12(F8, c, T, R, s, drops, delta, tog)
        if delta > 0:
            delta = 0
    rows = []
    for j in range(8, jtop, -1):
        if mux8:
            break
        idx = [i for i in range(j) if i + j >= 8]
        words = [{i + j - 8: (T[i], (i, j)) for i in idx}]
        if j >= jtop + 2:
            words.append({i + j - 9: (T[i], (i, j)) for i in idx if i + j - 9 >= 0})
        rows.append(("plain", j, idx, words))
    for j in range(5, jtop + 1):
        rows.append(("merged", j, list(range(8 - j, j)), None))
    st = {"k": 2, "neg": 0 if o.get("cinmap") else 1} if mux8 else {"k": 0, "neg": -delta}   # mux8: no s
    # borrow-ins left, one ONE (the -1) unless a cinmap places every borrow-in
    return _stage_rows(F8, c, o, T, L, R, s, Gd, c0, rows, st, drops, d42, acc11, par, tog)


def _lin12(F8, c, T, R, s, drops, delta, tog):
    A = F8.A
    t8 = c.alloc(1)[0]
    c.cx(T[8], t8)
    c.cx(s, t8)
    A.add(c, [T[i] for i in range(9)], R[1:11], cin=t8)
    c.cx(s, t8)
    c.cx(T[8], t8)
    c.free_([t8])
    # lin 2: Tf at half units (bit 0 dropped on request), carry-in +1 if delta > 0
    src2 = [T[i] for i in range(9)]
    if (0, 9, 1) in drops:
        src2[0] = None
    else:
        tog("Tf0")
    cin2 = None
    if delta > 0:
        assert delta == 1, "delta > 1 not supported"
        cin2 = F8.ONE
        tog("one")
        delta = 0
    A.add(c, src2, R, cin=cin2)


def _stage_rows(F8, c, o, T, L, R, s, Gd, c0, rows, st, drops, d42, acc11, par, tog):
    """the rows (plain and merged subtractions), the guard bit, D42's G0, the carry bits: the rest of stage()."""
    A = F8.A
    G0 = L[0] if d42 else None

    cinp = {tuple(q) for q in o.get("cinparts", [])}  # FC2 (mux8): dropped half-unit parts back as borrow-ins
    cmap = {int(k_): v_ for k_, v_ in (o.get("cinmap") or {}).items()}   # FC2 (mux8): row -> borrow-in literal

    def lit_ops(j, name):
        """(pre ops, cin wire, guard variables) of a borrow-in literal during merged row j's subtraction (the row's
        in-place state: T[i] = x_ij for i in idx, NOT x_ij for idx[1:], T[j+1] = e_j). n09 = NOT Tf_0."""
        if name == "one":
            return [], F8.ONE, ["one"]
        if name == "n09":
            return [("x", T[0])], T[0], ["Tf0", "one"]
        a, b = int(name[1]), int(name[2])
        idx = list(range(8 - j, j))
        neg = set(idx[1:])                             # wires holding NOT x_(i, j) during the sub
        held = set(idx[:1]) | {j + 1}                  # word sources: not to be touched
        if a in neg and b == j:                        # T[a] holds NOT x_(a, j) = NOT x_ab
            return [("x", T[a])], T[a], ["Tf%d" % a, "Tf%d" % b]
        if a in neg and b in neg:                      # NOT x_aj + NOT x_bj = x_ab
            return [("cx", T[b], T[a])], T[a], ["Tf%d" % a, "Tf%d" % b]
        if a not in idx and b not in idx and a not in held and b not in held:
            return [("cx", T[b], T[a])], T[a], ["Tf%d" % a, "Tf%d" % b]   # both untouched: Tf_a + Tf_b
        raise AssertionError("literal %s not formable in row %d" % (name, j))

    def sub(word, wvars, j=None):
        """R -= word (+ borrow-in): the first two take the sign wire (-s: the cubic's other half), then -1 each
        while delta < 0 asks for it. mux8 + cinparts: the merged rows' borrow-ins carry dropped parts (0 ANDs):
        row 5 NOT Tf_0 (= the -1 and part (0,9,1)), row 6 x_18 (part (1,8,1)), row 7 x_27 or x_36.
        mux8 + cinmap {row: literal}: any literal per row (lit_ops), the -1 only where the map says "one"."""
        k = st["k"]
        cin = None
        pre, post = [], []
        keep451 = bool(o.get("keep451", 0))
        if cmap and j in cmap:
            pre, cin, tv = lit_ops(j, cmap[j])
            post = list(reversed(pre))
            tog(*tv)
        elif cmap:
            pass
        elif j == 5 and (0, 9, 1) in cinp:
            assert st["neg"] == 1
            st["neg"] -= 1
            pre, post = [("x", T[0])], [("x", T[0])]
            cin = T[0]
            tog("Tf0", "one")
        elif j == 6 and (1, 8, 1) in cinp:
            pre, post = [("cx", T[8], T[1])], [("cx", T[8], T[1])]
            cin = T[1]
            tog("Tf1", "Tf8")
        elif j == 7 and (2, 7, 1) in cinp:            # T[2] holds NOT x_27 during row 7 (the OR step)
            pre, post = [("x", T[2])], [("x", T[2])]
            cin = T[2]
            tog("Tf2", "Tf7")
        elif j == 7 and (3, 6, 1) in cinp:            # NOT x_37 + NOT x_67 = x_36
            pre, post = [("cx", T[6], T[3])], [("cx", T[6], T[3])]
            cin = T[3]
            tog("Tf3", "Tf6")
        elif keep451 and k == 1:                       # FC2: the part (4,5,1) as a borrow-in: T[4] holds x_45
            cin = T[4]                                 # during the merged pair j = 5 (this sub)
            tog("Tf4", "Tf5")
        elif k < 2 or (keep451 and k == 2):
            cin = s
        elif st["neg"] > 0:
            cin = F8.ONE
            st["neg"] -= 1
            tog("one")
        st["k"] = k + 1
        if 0 in word:
            tog(*wvars)
        src = [word.get(pos) for pos in range(len(R))]
        while src and src[-1] is None:
            src.pop()
        if os.environ.get("FC2DBG"):                   # FC2 (fifth holder): the row words and their wires
            print("DBG h1b sub row %s live %d cin %s word %s" % (
                j, c.live, "y" if cin is not None else "n",
                "".join("." if w_ is None else ("1" if w_ == F8.ONE else "x") for w_ in src)))
        for op in pre:
            getattr(c, op[0])(*op[1:])
        for t in R:
            c.x(t)
        A.add(c, src, R, cin=cin)
        for t in R:
            c.x(t)
        for op in post:
            getattr(c, op[0])(*op[1:])

    for kind, j, idx, words in rows:
        for i in idx:
            c.cx(T[j], T[i])                           # T[i] = x_ij
        if kind == "plain":
            for w in words:
                w0 = w.get(0)
                sub({p: t for p, (t, _) in w.items()}, ["Tf%d" % w0[1][0], "Tf%d" % w0[1][1]] if w0 else [])
        else:                                          # M_j = e ? all ones : 2 X_j  (designer's merge)
            dropped0 = (8 - j, j + 1, 1) in drops
            c.cx(T[j], T[j + 1])                       # T[j + 1] = e_j
            e = T[j + 1]
            c.x(e)
            ors = []
            for i in idx[(1 if dropped0 else 0):]:
                c.x(T[i])
                t = c.and_compute(e, T[i])
                c.x(t)                                 # t = e OR x_ij
                ors.append((i, t))
            c.x(e)
            if dropped0:                               # positions 0, 1 carry x_0 + 2 e
                word = {0: T[idx[0]], 1: e}
                wv = ["Tf%d" % idx[0], "Tf%d" % j]
            else:
                word = {0: e}
                wv = ["Tf%d" % j, "Tf%d" % (j + 1)]
            for i, t in ors:
                word[1 + (i - idx[0])] = t
            k451 = bool(o.get("keep451")) and j == 5 and st["k"] == 1
            if k451:
                c.x(T[4])                              # the OR step left NOT x_45 here: the borrow-in needs x_45
            sub(word, wv, j)
            if k451:
                c.x(T[4])
            c.x(e)
            for i, t in reversed(ors):
                c.x(t)
                c.and_erase(e, T[i], t)
                c.x(T[i])
            c.x(e)
            c.cx(T[j], T[j + 1])
        for i in idx:
            c.cx(T[j], T[i])
    assert st["k"] >= 2 and st["neg"] == 0, ("borrow-ins not all placed", st)
    for v_ in sorted(par):                             # the guard bit is the XOR of par: clear it
        if v_ == "one":
            c.x(Gd)
        else:
            c.cx(T[int(v_[2:])], Gd)
    F8.__dict__["_fc2_dbg"] = dict(Gd=Gd, op=len(c.ops), T=list(T), L=list(L), c0=c0)   # (FC2 debug probe)
    c.free_([Gd])
    if d42:                                            # G0 = U's unit bit of the sum = a quadratic form of the
        const, lin, prods = g0_form(o)                 # folded T wires (z_i = T[i]): emitted as CCX on parities
        c.tag("G0E")
        F8.__dict__["_fc2_g0dbg"] = (G0, list(T), Gd)
        for la, lb in prods:
            a, b = la[0], lb[0]
            for z in la[1:]:
                c.cx(T[z], T[a])
            for z in lb[1:]:
                c.cx(T[z], T[b])
            c.ccx(T[a], T[b], G0)
            for z in lb[1:]:
                c.cx(T[z], T[b])
            for z in la[1:]:
                c.cx(T[z], T[a])
        for z in lin:
            c.cx(T[z], G0)
        if const:
            c.x(G0)
        c.free_([G0])
    if acc11:                                          # c1 = [C = 2] = NOT c0 AND NOT L'8 AND Tf8 AND Tf7 (T as it was)
        c1 = c.alloc(1)[0]
        cc = bool(o.get("c1cond"))
        if cc:                                         # FC2 c1cond: the mirror erases c1 by X-measurement + a CZ
            bid = len(c.blocks)                        # of conditioned ANDs (1 expected) instead of recomputing 3
            c.blocks[bid] = None
            c.ops.append(("BBEG", bid))
        c.x(c0)
        c.x(L[8])
        u1 = c.and_compute(T[8], T[7])
        u2 = c.and_compute(c0, L[8])
        c.ccx(u1, u2, c1)
        c.and_erase(c0, L[8], u2)
        c.and_erase(T[8], T[7], u1)
        c.x(L[8])
        c.x(c0)
        if cc:
            c.ops.append(("BEND", bid))
            mb = c.bit()
            pre = [("X", c0, None), ("X", L[8], None)]
            c.blocks[bid] = (pre + [("HMR", c1, mb), ("CANDC", T[8], T[7], u1, mb), ("CANDC", c0, L[8], u2, mb),
                                    ("CZ", u1, u2, None), ("CANDE", c0, L[8], u2, mb),
                                    ("CANDE", T[8], T[7], u1, mb)] + list(reversed(pre)))
    for i in range(9):
        c.cx(s, T[i])                                  # unfold
    A.add(c, [c0, c1], T)                              # T += C (mod 2^10)
    if o.get("c1debt"):
        # FC2 c1debt: c1 X-measured with a phase debt (-1)^(m c1); the mirror's explicit undo rebuilds c1 (3) and
        # repays the debt there by Z(c1) if m: 0 Toffolis forward (the narrow adders' debt mechanism, FC's kit)
        bid = len(c.blocks)
        c.blocks[bid] = None
        c.ops.append(("BBEG", bid))
        mb = c.hmr(c1)
        lits = [("CX", s, T[8], None), ("CX", s, T[7], None), ("X", c0, None), ("X", L[8], None)]
        u1, u2 = c.alloc(2)                            # wire ids for the undo's temporaries (free here)
        c.free_([u1, u2])
        c.free_([c1])
        c.ops.append(("BEND", bid))
        c.blocks[bid] = (lits + [("ANDC", T[8], T[7], u1), ("ANDC", c0, L[8], u2), ("CCX", u1, u2, c1),
                                 ("ANDE", c0, L[8], u2), ("ANDE", T[8], T[7], u1)] + list(reversed(lits))
                         + [("Z", c1, mb)])
        return c0
    if o.get("c1cond"):
        # FC2 c1cond: c1 erased by X-measurement, the phase (-1)^(m c1) repaid by CZ(u1, u2) with u1 = Tf'8 Tf'7 and
        # u2 = NOT c0 NOT L'8 computed only when m = 1 (conditioned ANDs, 1/2 expected each); the mirror's explicit
        # undo rebuilds c1 with 3.
        bid = len(c.blocks)
        c.blocks[bid] = None
        c.ops.append(("BBEG", bid))
        lits = [("CX", s, T[8], None), ("CX", s, T[7], None), ("X", c0, None), ("X", L[8], None)]
        for op in lits:
            c.ops.append(op)
        u1, u2 = c.alloc(2)
        mb = c.hmr(c1)
        c.cand_compute(T[8], T[7], u1, mb)
        c.cand_compute(c0, L[8], u2, mb)
        c.cz(u1, u2)
        c.cand_erase(c0, L[8], u2, mb)
        c.cand_erase(T[8], T[7], u1, mb)
        c.free_([u1, u2])
        c.free_([c1])
        for op in reversed(lits):
            c.ops.append(op)
        c.ops.append(("BEND", bid))
        c.blocks[bid] = (lits + [("ANDC", T[8], T[7], u1), ("ANDC", c0, L[8], u2), ("CCX", u1, u2, c1),
                                 ("ANDE", c0, L[8], u2), ("ANDE", T[8], T[7], u1)] + list(reversed(lits)))
        return c0
    if o.get("c1meas"):
        # FC2: c1 erased by X-measurement; the phase (-1)^(m c1), c1 = (Tf'8 Tf'7)(NOT c0 NOT L'8), repaid by a CZ
        # of two temporary ANDs (2 instead of 3). The mirror runs the explicit undo: c1 rebuilt by 3.
        bid = len(c.blocks)
        c.blocks[bid] = None
        c.ops.append(("BBEG", bid))
        lits = [("CX", s, T[8], None), ("CX", s, T[7], None), ("X", c0, None), ("X", L[8], None)]
        for op in lits:
            c.ops.append(op)
        mb = c.hmr(c1)
        u1 = c.and_compute(T[8], T[7])                 # (c1 stays allocated until the temps are gone: the undo
        u2 = c.and_compute(c0, L[8])                   # needs distinct wires for u1, u2, c1)
        c.cz(u1, u2, cond=mb)
        c.and_erase(c0, L[8], u2)
        c.and_erase(T[8], T[7], u1)
        c.free_([c1])
        for op in reversed(lits):
            c.ops.append(op)
        c.ops.append(("BEND", bid))
        c.blocks[bid] = (lits + [("ANDC", T[8], T[7], u1), ("ANDC", c0, L[8], u2), ("CCX", u1, u2, c1),
                                 ("ANDE", c0, L[8], u2), ("ANDE", T[8], T[7], u1)] + list(reversed(lits)))
        return c0
    c.cx(s, T[8])                                      # c1 = NOT c0 AND Tf'8 AND Tf'7 AND NOT L'8
    c.cx(s, T[7])
    c.x(c0)
    c.x(L[8])
    u1 = c.and_compute(T[8], T[7])
    u2 = c.and_compute(c0, L[8])
    u3 = c.and_compute(u1, u2)
    c.cx(u3, c1)
    c.and_erase(u1, u2, u3)
    c.and_erase(c0, L[8], u2)
    c.and_erase(T[8], T[7], u1)
    c.x(L[8])
    c.x(c0)
    c.cx(s, T[7])
    c.cx(s, T[8])
    c.free_([c1])
    return c0


def _make_rebuild_w25_c(F8):
    """FC's _rebuild_w25 (MAIN's word bit 25 rebuilt from F1's relation w25 = 1 + q21 + q20 w24) with a cheaper
    explicit undo (FC2, 17:00Z): the mirror erases w25 by an X-measurement and repays (-1)^(m q20 w24) by a
    conditioned CZ (0 Toffolis) instead of the CCX. Forward identical to FC's."""
    def rebuild(c, w24, v, ret):
        _, mbit = ret
        q20, q21 = v[6], v[7]
        bid = len(c.blocks)
        c.blocks[bid] = None
        c.ops.append(("BBEG", bid))
        t = c.alloc(1)[0]
        c.x(t)
        c.cx(q21, t)
        c.ccx(q20, w24, t)
        c.z(t, mbit)
        c.neg(mbit)
        c.ops.append(("BEND", bid))
        mu = c.bit()
        c.blocks[bid] = [("CX", q21, t, None), ("X", t, None), ("HMR", t, mu), ("CZ", q20, w24, mu)]
        return t
    return rebuild


def apply(F8, cfg):
    done = F8.__dict__.setdefault("_fc2_done", set())
    if "h1b" in done:
        return
    if cfg.get("h1b", {}).get("w25c"):
        F8._rebuild_w25 = _make_rebuild_w25_c(F8)
    orig = F8.st_main

    def st_main_h1b(c, U, v, cfg_):
        c.tag("H1")
        c0 = stage(F8, c, U, cfg_)
        if opts(cfg_).get("d42"):                      # D42: c0 is U's own LSB (no held wire)
            c.tag("MAIN")
            F8.__dict__["_fc2_u19"] = [c0] + list(U)
            return orig(c, [c0] + list(U), v, cfg_)
        if opts(cfg_)["form"] == "dump":
            cfg_["_pre_lz"] = {8: c0}
        c.tag("MAIN")
        return orig(c, U, v, cfg_)

    F8.st_main = st_main_h1b
    done.add("h1b")


def check(nsamp, over, seed=11):
    """model vs circuit on sampled inputs + edges: word = FC's ref_main on the model's U (held c0: no dump bit;
    dump: ref_main_g with c0), x restored, ancillas clean."""
    import invsqrt_fc_g8 as F8
    import invsqrt_fc_gen as G
    import invsqrt_fc_kit as K
    from invsqrt_fc2_h1lp import ref_main_g
    from invsqrt_fc2_h1m2 import stage_word
    cfg = dict(F8.CFG)
    cfg.update(json.load(open(os.path.expanduser("~/qcc/ws/invsqrt_fc/best/nw44wa/cfg.json"))))
    cfg.update({"nocorr": 1, "hx": 0, "fe": 0, "rowopt": 0, "u8late": 0, "bt": 7, "br": 8, "cap": 43,
                "fc2": ["norm", "h1b"], "fc2_src": None})
    cfg.update(over)
    o = opts(cfg)
    rng = np.random.default_rng(seed)
    x = np.unique(np.concatenate([rng.integers(1, F8.MAX_X + 1, nsamp), np.arange(1, 400),
                                  np.arange(F8.MAX_X - 400, F8.MAX_X + 1)]))
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
    dp = [list(t) for t in o["drops"]] + ([] if o.get("keep451") else [[4, 5, 1]])
    if o.get("d42"):                                   # FC2's D42 model (tools/invsqrt_fc2_d42m.py), K in acc units
        from invsqrt_fc2_d42m import stage_d42
        _, P2 = G.ref_div(q1 + cfg["c2"] * d, d, F8.a_top(cfg), cfg["sched2"])
        y2 = -((P2 + 1) >> 1)
        _, _, _, K0 = stage_d42(y2, 2, 0, 1, 1, dp)
        U, C, kept, Km = stage_d42(y2, 2, model_K(o) - K0, 1, 1, dp)
        Y = F8.ref_main(U, d, m, cfg["bt"])
        y = y2 * 2                                     # (for the c1 rule grid below: units of U)
    else:
        # FC2's model with K fixed: koff = K_target - round(knat)
        U0, C0, kept, K0 = stage_word(y, 1, 0, 0, 0, [256], (), dp)
        U, C, kept, Km = stage_word(y, 1, model_K(o) - K0, 0, 0, [256], (), dp)
        gb = (C & 1) if o["form"] == "dump" else np.zeros_like(C)
        Y = ref_main_g(U, gb, d, m, cfg["bt"])
    # the c1 erase rule on EVERY (T, L): c1 = NOT c0 AND Tf'8 AND Tf'7 AND NOT L'8
    if o.get("d42"):                                   # D42: (T, L2) grid, L2 8 bits, L'8 = L2's top bit
        Tg = np.repeat(np.arange(-512, 512, dtype=np.int64), 256)
        yg = (Tg << 8) + np.tile(np.arange(256, dtype=np.int64), 1024)
        Ug, Cg, _, _ = stage_d42(yg, 2, model_K(o) - K0, 1, 1, dp)
        Ug = Ug >> 1
        Tq, Lq, ltop = Ug >> 8, Ug & 255, 7
    else:
        Tg = np.repeat(np.arange(-512, 512, dtype=np.int64), 512)
        yg = (Tg << 9) + np.tile(np.arange(512, dtype=np.int64), 1024)
        Ug, Cg, _, _ = stage_word(yg, 1, model_K(o) - K0, 0, 0, [256], (), dp)
        Tq, Lq, ltop = Ug >> 9, Ug & 511, 8
    okr = (Cg >= 0) & (Cg <= 2) & (Tg + Cg <= 511)
    wq = np.where(Tq < 0, -1 - Tq, Tq)
    rule = (1 - (Cg & 1)) * (1 - ((Lq >> ltop) & 1)) * ((wq >> 8) & 1) * ((wq >> 7) & 1)
    c1bad = int((okr & (rule != (Cg >> 1))).sum())
    w0 = np.where(Tg < 0, -1 - Tg, Tg)                 # the compute rule with T before T += C (acc11)
    rule0 = (1 - (Cg & 1)) * (1 - ((Lq >> ltop) & 1)) * ((w0 >> 8) & 1) * ((w0 >> 7) & 1)
    c1bad += int((okr & (rule0 != (Cg >> 1))).sum())
    Yc = sim.get_word(regs["word"])
    sim.run(c.ops, nf, None)
    others = [w for w in range(c.next) if w not in regs["X"]]
    return dict(n=len(x), K=Km, Y_bad=int((Yc != Y).sum()), x_ok=bool(np.array_equal(sim.get_word(regs["X"]), x)),
                anc_ok=bool(not sim.q[others].any()), exp=c.toffolis_exp(), static=c.toffolis(), peak=c.peak,
                stages=c.stage_counts(), stage_peaks=c.stage_peaks, U_range=[int(U.min()), int(U.max())], C=[int(C.min()), int(C.max())],
                c1_rule_bad=c1bad, c1_rule_neg=int(((Cg < 0) & (Tg >= -512)).sum()))


if __name__ == "__main__":
    if sys.argv[1] == "check":
        print(json.dumps(check(int(sys.argv[2]), json.loads(sys.argv[3]) if len(sys.argv) > 3 else {})))
