#!/usr/bin/env python3
"""invsqrt_fc17_skipe.py -- FC seventeenth holder: R65-INV-FILE38's "skip E" on the hold frame (CH1's row 24 = P6 + P9's
digit-23 repay, FC's row 23 = invsqrt_fc16_r24n.block23).  Copied into a source snapshot; active only with cfg
"_fc17skipe": 1 (and "_ch1r24": 1, "_fc16r23n": CAP).

TODAY: row 23's undo re-makes w' (= row 24's output bit 8, o8') by ROW9's 3-AND form E and repays its forward's debt
(Z under mt); row 24's undo then X-measures o8' at once (bit m8) and runs its narrow game with o8' as a phase target.
SKIP E: o8' is never re-made.  Row 23's forward writes its X-measurement of o8' into a bit m8 RESERVED by row 24's undo
(row 24's undo is generated first but executed last), so the pending phase at row 24's undo is (-1)^(m8 o8') as today.
Row 24's undo erases digit 23 (o9') by X-measurement (bit m9, P9): o9' = A0 ^ o8' B0 with
  A0 = P Q ^ v6 o7',  P = o7'^v4^v6^d,  Q = v4^v5^v6^d^1,  B0 = v4^v5^d^1^o7'                     [P: algebra of CH1's form]
so (-1)^(m9 o9') = Cliffords (A0) times (-1)^(m9 B0 o8').  o8' is absent; with the game's toggles of the phase target
o8'(start) = w7' ^ [x y ^ c7 ^ a8] (the top step) where w7' is the restored input bit 8 (CH1's form A of v, d, w6').
At the top step every term is a form of live wires (s-form: x = a7 ^ c7, y = a7 ^ s7; o7' = s7 ^ 1 ^ a7 ^ c7;
w6' = s7 ^ 1), so the m9 part is paid there:  m9 (B0 x y ^ B0 (c7 ^ a8 ^ lin(w7')) ^ B0 F1 F2),  B0 x y = G x y with
G = v4^v5^d^a7 [P]: two CCZ under m9 (0.5 a shot each, kit op XCCZ of F3's fcgen) and CZs.  Saved: E's 3 ANDs.
"""
import sys
from itertools import permutations

ONE = "1"


# ------------------------------------------------------------------ affine forms over wires: (frozenset, const)
def F(*items):
    """XOR of items: a wire id (int), the constant "1", a form (frozenset, const) or None (nothing)."""
    ws, k = set(), 0
    for it in items:
        if it is None:
            continue
        if isinstance(it, tuple):
            W, kk = it
            ws ^= set(W)
            k ^= kk
        elif it == ONE:
            k ^= 1
        elif isinstance(it, int):
            ws ^= {it}
        else:
            raise ValueError(it)
    return (frozenset(ws), k)


def _lit_src(a):
    if isinstance(a, int):
        return F(a)
    return (frozenset(), 1 if a == ONE else 0)


def _lit_carry(cr):
    if cr[0] == "k":
        return (frozenset(), cr[1])
    return (frozenset({cr[1]}), cr[2])


def _xor(a, b):
    return (a[0] ^ b[0], a[1] ^ b[1])


def _hosts(fs):
    """a private wire per form (in no other form), or None."""
    out = []
    for i, (W, _) in enumerate(fs):
        others = set()
        for j, (W2, _) in enumerate(fs):
            if j != i:
                others |= W2
        p = sorted(W - others)
        if not p:
            return None
        out.append(p[0])
    return out


def _clean(fs):
    """drop constant-1 factors, merge equal factors; None if the product is identically 0."""
    out = []
    for W, k in fs:
        if not W:
            if k == 0:
                return None
            continue
        if (W, k) in out:
            continue
        if (W, k ^ 1) in out:
            return None                                 # L (L ^ 1) = 0
        out.append((W, k))
    return out


def plan_prod(fs, depth=4):
    """terms (lists of forms with private hosts) whose phases multiply to (-1)^(prod fs); cost = (#3, #2) minimal."""
    fs = _clean(fs)
    if fs is None:
        return [], (0, 0)
    if len(fs) <= 1 or _hosts(fs) is not None:
        return [fs], ((1 if len(fs) == 3 else 0), (1 if len(fs) == 2 else 0))
    if depth == 0:
        raise ValueError("plan_prod: no host assignment")
    best = None
    for i, j in permutations(range(len(fs)), 2):
        if not (fs[i][0] & fs[j][0]):
            continue
        # Li Lj R = (Li ^ Lj) Lj R ^ Lj R
        a = list(fs)
        a[i] = _xor(fs[i], fs[j])
        b = [f for t, f in enumerate(fs) if t != i]
        try:
            ta, ca = plan_prod(a, depth - 1)
            tb, cb = plan_prod(b, depth - 1)
        except ValueError:
            continue
        cost = (ca[0] + cb[0], ca[1] + cb[1])
        if best is None or cost < best[1]:
            best = (ta + tb, cost)
    if best is None:
        raise ValueError("plan_prod: no plan")
    return best


def emit_prod(c, fs, g):
    """apply (-1)^(g * prod of the affine forms fs) (g a classical bit), 1 to 3 forms; CNOT folds onto private hosts."""
    terms, cost = plan_prod(fs)
    n3 = 0
    for t in terms:
        if not t:
            c.neg(g)
            continue
        if len(t) == 1:
            W, k = t[0]
            for w in sorted(W):
                c.z(w, g)
            if k:
                c.neg(g)
            continue
        hs = _hosts(t)
        fold = []
        for (W, k), h in zip(t, hs):
            for w in sorted(W - {h}):
                fold.append(("CX", w, h))
            if k:
                fold.append(("X", h))
        for o in fold:
            c.cx(o[1], o[2]) if o[0] == "CX" else c.x(o[1])
        if len(t) == 2:
            c.cz(hs[0], hs[1], g)
        else:
            assert len(t) == 3
            c.ops.append(("XCCZ", hs[0], hs[1], hs[2], g))
            n3 += 1
        for o in reversed(fold):
            c.cx(o[1], o[2]) if o[0] == "CX" else c.x(o[1])
    return n3


# ------------------------------------------------------------------ the narrow game with a hook at its top step
def game_add_hooked(c, src, tgt, cin, debts, C, top_hook):
    """invsqrt_fc_narrow.game_add (FC's narrow undo, copied) + top_hook(c, src, tgt, carry, n) after the top step."""
    NW = sys.modules["invsqrt_fc_narrow"]
    _lit, _emit, _xor_lit, _xor_src, _form_s = NW._lit, NW._emit, NW._xor_lit, NW._xor_src, NW._form_s
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    assert n >= 3
    carry = {0: _lit(cin)}

    def write_sum(r):
        _xor_src(c, src[r], tgt[r])
        _xor_lit(c, carry[r], tgt[r])

    write_sum(0)
    Fp = 0
    while Fp + 1 <= n - 2:
        f = _form_s(src[Fp], tgt[Fp], carry[Fp])
        if f[0] == "and":
            break
        carry[Fp + 1] = ("k", f[1]) if f[0] == "const" else ("w", f[1], f[2])
        Fp += 1
        write_sum(Fp)
    L = n - 2 - Fp
    fresh = set()
    hooked = [0]

    def place(r, first, cond=None):
        f = _form_s(src[r - 1], tgt[r - 1], carry[r - 1])
        assert f[0] == "and"
        _, xw, yw, ops = f
        t = c.alloc(1)[0]
        _emit(c, ops)
        if cond is None:
            c.and_compute(xw, yw, t)
        else:
            c.cand_compute(xw, yw, t, cond)
        _emit(c, ops[::-1])
        _xor_lit(c, carry[r - 1], t, cond)
        carry[r] = ("w", t, 0)
        fresh.add(r)
        if first:
            write_sum(r)
            if r in debts:
                c.z(t, debts[r])

    def remove(r, cond=None):
        t = carry[r][1]
        _xor_lit(c, carry[r - 1], t, cond)
        f = _form_s(src[r - 1], tgt[r - 1], carry[r - 1])
        _, xw, yw, ops = f
        _emit(c, ops)
        if cond is None:
            c.and_erase(xw, yw, t)
        else:
            c.cand_erase(xw, yw, t, cond)
            c.free_([t])
        _emit(c, ops[::-1])
        del carry[r]
        fresh.discard(r)

    def ghost(r):
        t = carry[r][1]
        g = c.hmr(t)
        c.free_([t])
        del carry[r]
        fresh.discard(r)
        return g

    def repay(r, g):
        lo = carry[r - 1]
        if lo[0] == "k":
            if lo[1]:
                c.neg(g)
        else:
            c.z(lo[1], g)
            if lo[2]:
                c.neg(g)
        f = _form_s(src[r - 1], tgt[r - 1], carry[r - 1])
        _, xw, yw, ops = f
        _emit(c, ops)
        c.cz(xw, yw, g)
        _emit(c, ops[::-1])

    def top():
        tw = tgt[n - 1]
        f = _form_s(src[n - 2], tgt[n - 2], carry[n - 2])
        if f[0] == "const":
            if f[1]:
                c.x(tw)
        elif f[0] == "alias":
            c.cx(f[1], tw)
            if f[2]:
                c.x(tw)
        else:
            _, xw, yw, ops = f
            _emit(c, ops)
            c.ccx(xw, yw, tw)
            _emit(c, ops[::-1])
            _xor_lit(c, carry[n - 2], tw)
        _xor_src(c, src[n - 1], tw)
        top_hook(c, src, tgt, carry, n, f)
        hooked[0] += 1

    moves, cost = NW.plan(L, C)
    if moves is None:
        raise ValueError("skipE: no narrow schedule for L=%d at C=%d" % (L, C))
    cb = {"place": lambda k, first, cond: place(Fp + k, first, cond),
          "remove": lambda k, cond: remove(Fp + k, cond),
          "ghost": lambda k: ghost(Fp + k),
          "repay": lambda k, b: repay(Fp + k, b),
          "top": top,
          "mark": lambda m: c.ops.append(m)}
    got = NW.run_moves(L, C, list(moves), cb)
    assert got is not None and abs(got - cost) < 1e-9
    assert not fresh
    assert hooked[0] == 1, ("skipE: the top step ran %d times" % hooked[0])
    return cost


def hybrid_add_hooked(c, src, tgt, cin, debts, C, top_hook):
    """invsqrt_fc_narrow.hybrid_add (copied) + top_hook right after the downward step at position n-2: there the
    source bit a_(n-2), the carry c_(n-2) and the sum s_(n-2) sit on their wires (the game's s-form at its top)."""
    NW = sys.modules["invsqrt_fc_narrow"]
    _w, _lit, _xor_lit, _xor_src = NW._w, NW._lit, NW._xor_lit, NW._xor_src
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    assert n >= 3
    alias = NW._cuccaro_alias(src, n, cin)
    i0 = 1 if alias else 0
    budget = C
    mats = []
    if alias:
        c0 = tgt[0]
    elif _w(cin):
        c0 = cin
    else:
        c0 = c.alloc(1)[0]
        if cin == ONE:
            c.x(c0)
        mats.append((c0, cin))
        budget -= 1
    pos = list(range(i0, n - 1))
    clean = {i for i in pos if not _w(src[i])}
    budget -= len(clean)
    assert budget >= 0, "hybrid_add: constants need more clean wires than C"
    for i in sorted((i for i in pos if _w(src[i])), reverse=True):
        if budget <= 0:
            break
        clean.add(i)
        budget -= 1
    chain = {i0: c0}
    cnt = 0
    for i in pos:
        x, y = chain[i], tgt[i]
        if i in clean:
            a = src[i]
            c.cx(x, y)
            if _w(a):
                c.cx(x, a)
                t = c.and_compute(a, y)
                c.cx(x, a)
            elif a == ONE:
                c.x(x)
                t = c.and_compute(x, y)
                c.x(x)
            else:
                t = c.and_compute(x, y)
            c.cx(x, t)
            cnt += 1
            chain[i + 1] = t
        else:
            z = src[i]
            c.cx(z, y)
            c.cx(z, x)
            c.ccx(x, y, z)
            cnt += 1
            chain[i + 1] = z
        if (i + 1) in debts:
            c.z(chain[i + 1], debts[i + 1])
    _xor_src(c, src[n - 1], tgt[n - 1])
    c.cx(chain[n - 1], tgt[n - 1])
    hooked = 0
    for i in reversed(pos):
        x, y, nxt = chain[i], tgt[i], chain[i + 1]
        if i in clean:
            a = src[i]
            c.cx(x, nxt)
            if _w(a):
                c.cx(x, a)
                c.and_erase(a, y, nxt)
                c.cx(x, a)
                c.cx(a, y)
            elif a == ONE:
                c.x(x)
                c.and_erase(x, y, nxt)
                c.x(x)
                c.x(y)
            else:
                c.and_erase(x, y, nxt)
        else:
            z = nxt
            c.ccx(x, y, z)
            c.cx(z, x)
            c.cx(x, y)
            cnt += 1
        if i == n - 2:
            top_hook(c, src, tgt, {n - 2: ("w", chain[n - 2], 0)}, n, ("and",))
            hooked += 1
    assert hooked == 1
    if alias:
        c.x(tgt[0])
    for w, xv in reversed(mats):
        if xv == ONE:
            c.x(w)
        c.free_([w])
    return cnt


def cuccaro_add_hooked(c, src, tgt, cin, debts, top_hook):
    """invsqrt_fc_narrow.cuccaro_add (copied) + top_hook right after UMA at position n-2 (s-form state there)."""
    NW = sys.modules["invsqrt_fc_narrow"]
    _w, _lit, _xor_lit, _xor_src = NW._w, NW._lit, NW._xor_lit, NW._xor_src
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    assert n >= 3
    mats = []

    def mat(x):
        if _w(x):
            return x
        w = c.alloc(1)[0]
        if x == ONE:
            c.x(w)
        mats.append((w, x))
        return w

    i0 = 0
    if NW._cuccaro_alias(src, n, cin):
        i0 = 1
        c0 = tgt[0]
        A_ = [None] + [mat(src[i]) for i in range(1, n - 1)]
    else:
        c0 = mat(cin)
        A_ = [mat(src[i]) for i in range(n - 1)]
    chain = [c0] + A_[:]
    if i0 == 1:
        chain[1] = tgt[0]
    cnt = 0
    for i in range(i0, n - 1):
        x, y, z = chain[i], tgt[i], A_[i]
        c.cx(z, y)
        c.cx(z, x)
        c.ccx(x, y, z)
        cnt += 1
        if (i + 1) in debts:
            c.z(z, debts[i + 1])
    _xor_src(c, src[n - 1], tgt[n - 1])
    c.cx(A_[n - 2] if n - 2 >= i0 else chain[n - 1], tgt[n - 1])
    hooked = 0
    for i in range(n - 2, i0 - 1, -1):
        x, y, z = chain[i], tgt[i], A_[i]
        c.ccx(x, y, z)
        c.cx(z, x)
        c.cx(x, y)
        cnt += 1
        if i == n - 2:
            top_hook(c, src, tgt, {n - 2: ("w", chain[n - 2], 0)}, n, ("and",))
            hooked += 1
    assert hooked == 1
    if i0 == 1:
        c.x(tgt[0])
    for w, x in reversed(mats):
        if x == ONE:
            c.x(w)
        c.free_([w])
    return cnt


_HOOKED = {}


def _hooked_gen_add(hook):
    """F3's fcgen.gen_add (the R32 game's adaptive policy) rebuilt from its own source with ONE line added at the end
    of op_top(g, carry): hook(c, src, tgt, carry, n, f, g) (g = the path's predicate bit, None at the root)."""
    import inspect
    G = sys.modules["invsqrt_fw_fcgen"]
    src_txt = inspect.getsource(G.gen_add)
    old = "        xsrc(src[n - 1], tw, g)\n\n    placed = set()"
    assert src_txt.count(old) == 1, "fcgen.gen_add: op_top anchor not found"
    src_txt = src_txt.replace(old, "        xsrc(src[n - 1], tw, g)\n        _FC17_HOOK(c, src, tgt, carry, n, f, g)\n\n"
                                   "    placed = set()")
    ns = dict(G.__dict__)
    ns["_FC17_HOOK"] = hook
    exec(compile(src_txt, "<fc17 gen_add>", "exec"), ns)
    return ns["gen_add"]


def mixed_add_hooked(c, src, tgt, cin, debts, C, k, ccz, hook):
    """F3's c1add.mixed_add rebuilt from its own source with its game part (fcgen.gen_add) replaced by the hooked
    copy; the hook sees the top step of every policy path."""
    import inspect
    C1A = sys.modules["invsqrt_fw_c1add"]
    gen = _hooked_gen_add(hook)
    src_txt = inspect.getsource(C1A.mixed_add)
    assert src_txt.count("G.gen_add(") == 1, "c1add.mixed_add: gen_add call not found"
    src_txt = src_txt.replace("G.gen_add(", "_FC17_GEN(")
    ns = dict(C1A.__dict__)
    ns["_FC17_GEN"] = gen
    exec(compile(src_txt, "<fc17 mixed_add>", "exec"), ns)
    return ns["mixed_add"](c, src, tgt, cin, debts, C, k, ccz)


# ------------------------------------------------------------------ row 24 (CH1's P6 forward, skip-E undo)
LOG = []


def r24_block(c, src, tgt, info):
    K4 = sys.modules["invsqrt_fc2_k4fuse"]
    W = sys.modules["invsqrt_fk_w24ret"]
    NW = sys.modules["invsqrt_fc_narrow"]
    install_uret_cap(c)                            # debug only (FC17_CAP + FC17_URET)
    n = len(tgt)
    assert n == 10, n
    src = list(src) + [None] * (n - len(src))
    top, t8, d, v = tgt[9], tgt[8], info["ctl"], info["v"]
    assert len(v) == 7
    v3, v4, v5, v6 = v[3], v[4], v[5], v[6]
    if not hasattr(c, "blocks"):
        c.blocks = {}
    bid = len(c.blocks)
    c.blocks[bid] = None
    # ---- forward: CH1's P6, unchanged
    c.ops.append(("BBEG", bid))
    W._r24_form(c, tgt, d, v, top, erase=True)
    c.free.insert(0, top)
    c.live -= 1
    K4._ch1_w7(c, tgt, d, v, t8, erase=True)
    c.free.insert(0, t8)
    c.live -= 1

    def loader(c_, tg, tp, lz):
        K4._ch1_w7(c_, tgt, d, v, tp, erase=False)
    debts = W._fwd_lazy(c, NW, src[:9], tgt[:9], t8, info, loader=loader)
    assert t8 not in c.free, "row 24 at 38: w7 was not taken"
    c.free.remove(top)
    c.live += 1
    c.peak = max(c.peak, c.live)
    cap_mark(c, 3, list(v) + [d] + list(tgt[:9]))   # debug (FC17_CAP): the inputs of digit 23's birth
    K4._ch1_o9(c, tgt, d, v, top, erase=False)
    cap_mark(c, 4, [top])                           # debug (FC17_CAP): digit 23 (framed) on the top wire
    c.ops.append(("BEND", bid))
    # ---- undo (skip E): executed after row 23's undo, which leaves o8' (wire t8) at |0> with (-1)^(m8 o8') pending
    save = c.ops
    c.ops = []
    o7 = tgt[7]
    m8 = c.bit()                                   # RESERVED: row 23's forward X-measures o8' into this bit
    # Between row 23's undo and this undo the mirror runs row 23's unframe (X, CX(d23)) and row 24's frame (X,
    # CX(ctl24)) on every window wire, and the top becomes 1 ^ d23 ^ ctl24 (three X's on it); so the o8' wire, |0>
    # after row 23's skip-E undo, holds V = d23 ^ ctl24 = top ^ 1 here [M: the op list, FC17_DBG].  Clear it, and
    # split the pending phase (-1)^(m8 W) = (-1)^(m8 o8') (-1)^(m8 V)  (o8' = W ^ V is this undo's game target).
    c.cx(top, t8)
    c.x(t8)
    c.z(top, m8)
    c.neg(m8)
    m9 = c.hmr(top)                                # digit 23 X-measured (P9)
    emit_prod(c, [F(o7, v4, v6, d), F(v4, v5, v6, d, ONE)], m9)        # A0 = P Q ^ v6 o7'
    emit_prod(c, [F(v6), F(o7)], m9)
    c.free.insert(0, top)
    c.live -= 1
    c.free.insert(0, t8)                           # o8' absent (row 23's undo did not re-make it)
    c.live -= 1
    C = c.cap - c.live
    print("FC17 r24 undo: cap %d live %d C %d" % (c.cap, c.live, C), file=sys.stderr)
    assert C >= 2, ("row 24 skipE: undo needs two clean wires", C)
    SENT = -1000003
    tg2 = list(tgt[:8]) + [SENT]
    s9 = src[:9]
    nccz = [0]

    def hook(c_, src_, tgt_, carry, n_, f, g=None):
        cond = m9 if g is None else c_.band(((g, 1), (m9, 1)))   # fcgen's policy paths: under the path's bit too
        a7, s7, a8 = src_[n_ - 2], tgt_[n_ - 2], src_[n_ - 1]
        A7, A8, C7, S7 = _lit_src(a7), _lit_src(a8), _lit_carry(carry[n_ - 2]), F(s7)
        B0 = F(v4, v5, d, S7, A7, C7)              # = v4^v5^d^1^o7', o7' = s7^1^a7^c7
        w6 = F(S7, ONE)                            # restored input bit 7 (the frame's X still on)
        lin = F(v3, v4, v5, v6, w6)                # w7' = lin ^ F1 F2 (CH1's form A)
        F1, F2 = F(v3, d, w6), F(v3, v5)
        if f[0] == "and":                          # top toggles: x y ^ c7 ^ a8
            G = F(v4, v5, d, A7)                   # B0 x y = G x y  [P]
            nccz[0] += emit_prod(c_, [G, F(A7, C7), F(A7, S7)], cond)
            rest = F(C7, A8, lin)
        elif f[0] == "const":
            rest = F(A8, lin, (frozenset(), f[1]))
        else:                                      # alias: wire f[1] with flip f[2]
            rest = F(A8, lin, F(f[1]), (frozenset(), f[2]))
        emit_prod(c_, [B0, rest], cond)
        nccz[0] += emit_prod(c_, [B0, F1, F2], cond)

    import os
    i0 = len(c.ops)
    for w in tg2:
        c.x(w)
    C1A = sys.modules.get("invsqrt_fw_c1add")
    ccz = int(os.environ.get("C1CCZ", "1"))
    mp = None
    if (C1A is not None and C1A.NW is not None and os.environ.get("C1UNDO", "0") == "1"
            and len(tg2) >= 3 and (len(tg2) - 1) not in (debts or {}) and C1A.c1_need(s9, tg2, None) <= C
            and os.environ.get("C1MIX", "0") == "1"):
        mp = C1A.mixed_plan(s9, tg2, None, C, ccz)     # F3's undo route, as NW.hybrid_add (C1UNDO) takes it
    if NW.plan(max(0, 9 - 2 - NW.free_prefix(s9, 9, None)), C)[0] is not None:
        cost, kind = game_add_hooked(c, s9, tg2, None, debts, C, hook), "game"
    elif mp is not None:
        cost, kind = float(mixed_add_hooked(c, s9, tg2, None, debts, C, mp[0], ccz, hook)), "mixed k=%d" % mp[0]
    elif C > NW.cuccaro_need(s9, tg2, None):
        cost, kind = float(hybrid_add_hooked(c, s9, tg2, None, debts, C, hook)), "hybrid"
    else:
        cost, kind = float(cuccaro_add_hooked(c, s9, tg2, None, debts, hook)), "cuccaro"
    for w in tg2:
        c.x(w)
    body = c.ops[i0:]
    del c.ops[i0:]
    cost -= K4._phase_rewrite(c, body, SENT, m8)
    c.free.remove(t8)
    c.live += 1
    K4._ch1_w7(c, tgt, d, v, t8, erase=False)      # w7' re-made, then the m8 part of the debt
    c.z(t8, m8)
    c.free.remove(top)
    c.live += 1
    W._r24_form(c, tgt, d, v, top, erase=False)    # w8' re-made
    cost += 2.0 + 0.5 * nccz[0]
    undo_ops = c.ops
    c.blocks[bid] = undo_ops
    c.ops = save
    if not hasattr(c, "block_log"):
        c.block_log = []
    c.block_log.append((n, C, cost))
    LOG.append(dict(fc17r24=int(top), C_undo=C, kind=kind, undo_cost=cost, xccz=nccz[0], m8=m8, m9=m9))
    print("FC17 skipE row 24: undo kind %s C %d xccz %d" % (kind, C, nccz[0]), file=sys.stderr)
    c._fc17 = dict(m8=m8, wire=t8)
    install_r23()
    return None


# ------------------------------------------------------------------ row 23 (FC's block23 with skip E)
def install_r23():
    N = sys.modules["invsqrt_fc16_r24n"]
    if getattr(N.block23, "_fc17", False):
        return
    orig = N.block23

    def block23(c, src, tgt, info, cap):
        st = getattr(c, "_fc17", None)
        if not st:
            return orig(c, src, tgt, info, cap)
        c._fc17 = None
        return block23_skipe(c, src, tgt, info, cap, st)
    block23._fc17 = True
    N.block23 = block23


def block23_skipe(c, src, tgt, info, cap, st):
    """invsqrt_fc16_r24n.block23 with (1) the forward's X-measurement of w' written into row 24's reserved bit and
    (2) the undo ending after its game: w' (= o8') is not re-made and its debt is not repaid here (row 24's undo
    carries it as its phase target)."""
    N = sys.modules["invsqrt_fc16_r24n"]
    NW = sys.modules["invsqrt_fc_narrow"]
    _cx_all = N._cx_all
    n = len(tgt)
    assert n == 10 and list(tgt[1:]) == info["expect"], "r23n: row 23 not where expected"
    src = list(src) + [None] * (n - len(src))
    d, v = info["ctl"], info["v"]
    v2, v3, v4, v5, v6 = v[2], v[3], v[4], v[5], v[6]
    top, t8 = tgt[9], tgt[8]
    assert top == st["wire"], ("skipE: row 23's top is not row 24's o8' wire", top, st)
    peak = [c.live]

    def chk():
        peak[0] = max(peak[0], c.live)
        assert c.live <= cap, ("r23n: live %d above the cap %d" % (c.live, cap))
    a_orig = c.alloc

    def alloc_chk(k=1):
        out = a_orig(k)
        chk()
        return out
    chk()
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    i0 = len(c.ops)
    c.ops.append(("HMR", top, st["m8"]))            # the debt (-1)^(m8 w'), carried by row 24's undo
    c.free.insert(0, top)
    c.live -= 1
    c.alloc = alloc_chk
    try:
        debts = NW.fwd(c, src[:9], list(tgt[:9]), None)
    finally:
        c.alloc = a_orig
    c.free.remove(top)
    c.live += 1
    chk()
    cap_mark(c, 1, list(v) + [d] + list(tgt[:9]))   # debug (FC17_CAP): the inputs of digit 22's birth
    h = c.alloc(1)[0]
    chk()
    _cx_all(c, (v5, v6, t8), v3)
    c.cx(v6, v4); c.cx(d, v4); c.x(v4)
    c.and_compute(v3, v4, h)
    c.x(v4); c.cx(d, v4); c.cx(v6, v4)
    _cx_all(c, (t8, v6, v5), v3)
    _cx_all(c, (v3, v5, v6, t8, d), h)
    c.cx(v5, v3); c.cx(v6, v3); c.x(v3)
    c.and_compute(v3, h, top)
    c.x(v3); c.cx(v6, v3); c.cx(v5, v3)
    _cx_all(c, (d, t8, v6, v5, v3), h)
    c.cx(t8, v5); c.x(v5)
    c.cx(v6, h); c.cx(t8, h)
    c.ccx(v5, h, top)
    c.cx(t8, h); c.cx(v6, h)
    c.x(v5); c.cx(t8, v5)
    c.cx(h, top)
    _cx_all(c, (v3, v4, v5, t8, d), top); c.x(top)
    _cx_all(c, (v5, v6, t8), v3)
    c.cx(v6, v4); c.cx(d, v4); c.x(v4)
    c.and_erase(v3, v4, h)
    c.x(v4); c.cx(d, v4); c.cx(v6, v4)
    _cx_all(c, (t8, v6, v5), v3)
    cap_mark(c, 2, [top])                           # debug (FC17_CAP): digit 22 (framed) on the top wire
    fwd_cost = float(c.toffolis_exp(c.ops[i0:]))
    c.ops.append(("BEND", bid))
    # ---------------- explicit undo: o9' erased, the game; NO re-make of w' (skip E)
    save = c.ops
    c.ops = []
    _cx_all(c, (v3, v4, v5, t8, d), top); c.x(top)
    m9 = c.hmr(top)
    c.free_([top])
    h = c.alloc(1)[0]
    chk()
    _cx_all(c, (v5, v6, t8), v3)
    c.cx(v6, v4); c.cx(d, v4); c.x(v4)
    c.cz(v3, v4, m9)
    c.cand_compute(v3, v4, h, m9)
    c.x(v4); c.cx(d, v4); c.cx(v6, v4)
    _cx_all(c, (t8, v6, v5), v3)
    _cx_all(c, (v3, v5, v6, t8, d), h)
    c.cx(v5, v3); c.cx(v6, v3); c.x(v3)
    c.cz(v3, h, m9)
    c.x(v3); c.cx(v6, v3); c.cx(v5, v3)
    _cx_all(c, (d, t8, v6, v5, v3), h)
    c.cx(t8, v5); c.x(v5)
    c.cx(v6, h); c.cx(t8, h)
    c.cz(v5, h, m9)
    c.cx(t8, h); c.cx(v6, h)
    c.x(v5); c.cx(t8, v5)
    _cx_all(c, (v5, v6, t8), v3)
    c.cx(v6, v4); c.cx(d, v4); c.x(v4)
    c.cand_erase(v3, v4, h, m9)
    c.x(v4); c.cx(d, v4); c.cx(v6, v4)
    _cx_all(c, (t8, v6, v5), v3)
    c.free_([h])
    C = cap - c.live
    gcost, kind = N._game(c, NW, src[:9], list(tgt[:9]), debts, C, alloc_chk)
    undo_cost = float(c.toffolis_exp(c.ops))
    c.blocks[bid] = c.ops
    c.ops = save
    c.free.remove(top)                              # the forward goes on with o9' on the top wire
    c.live += 1
    c.block_log.append((n, C, undo_cost))
    LOG.append(dict(fc17r23=int(top), cap=cap, C_undo=C, game=kind, game_cost=gcost, fwd=fwd_cost, undo=undo_cost,
                    peak=peak[0], debts=len(debts)))


# ------------------------------------------------------------------ debug: every block's undo count at exit
_DUMP = {}


def register_dump(c):
    """FC17_BLK=path: at exit, write one line per block (bid, undo Toffoli-exp, undo ops) of the LAST build."""
    import os
    import atexit
    p = os.environ.get("FC17_BLK")
    if not p:
        return
    _DUMP["c"] = c
    if _DUMP.get("reg"):
        return
    _DUMP["reg"] = True

    def dump():
        cc = _DUMP["c"]
        with open(p, "w") as f:
            for bid in sorted(k for k in cc.blocks if isinstance(cc.blocks[k], list)):
                ops = cc.blocks[bid]
                f.write("%d %.3f %d\n" % (bid, float(cc.toffolis_exp(ops)), len(ops)))
    atexit.register(dump)


def register_opdump(c):
    """FC17_OPS=path:bid,bid: at exit, write those blocks' undo op lists (debug)."""
    import os
    import atexit
    spec = os.environ.get("FC17_OPS")
    if not spec or _DUMP.get("opreg"):
        _DUMP["c"] = c
        return
    _DUMP["opreg"] = True
    _DUMP["c"] = c
    path, bids = spec.split(":")

    def dump():
        cc = _DUMP["c"]
        with open(path, "w") as f:
            for b in bids.split(","):
                ops = cc.blocks.get(int(b))
                f.write("### block %s: %d ops, Toffoli-exp %.3f\n" % (b, len(ops or []), float(cc.toffolis_exp(ops or []))))
                for o in ops or []:
                    f.write(repr(o) + "\n")
    atexit.register(dump)


# ------------------------------------------------------------------ debug: capture wire values at marker ops
_CAP = {"marks": {}, "snap": {}, "done": False}


def cap_mark(c, tag, wires):
    """FC17_CAP=path: put a marker op (BBEG, -9000 - tag) in the op list and remember the wires to snapshot there."""
    import os
    if not os.environ.get("FC17_CAP"):
        return
    c.ops.append(("BBEG", -9000 - tag))
    _CAP["marks"][-9000 - tag] = list(wires)
    _install_cap()


def _install_cap():
    import os
    import atexit
    if _CAP.get("inst"):
        return
    _CAP["inst"] = True
    K = sys.modules["invsqrt_fc_kit"]
    orig = K.Sim.run

    def run(self, ops, start=0, stop=None):
        stop_ = len(ops) if stop is None else stop
        import os as _o
        if _CAP["done"] or not _CAP["marks"] or getattr(self, "n", 0) < int(_o.environ.get("FC17_CAPN", "1000000")):
            return orig(self, ops, start, stop)
        cut = [i for i in range(start, stop_) if ops[i][0] == "BBEG" and ops[i][1] in _CAP["marks"]]
        a = start
        for i in cut:
            orig(self, ops, a, i)
            key = ops[i][1]
            if key not in _CAP["snap"]:
                _CAP["snap"][key] = {w: np.asarray(self.q[w]).copy() for w in _CAP["marks"][key]}
            a = i + 1
        res = orig(self, ops, a, stop)
        if len(_CAP["snap"]) == len(_CAP["marks"]):
            _CAP["done"] = True
        return res
    K.Sim.run = run

    def save():
        out = {}
        for key, d in _CAP["snap"].items():
            for w, arr in d.items():
                out["m%d_q%d" % (-key - 9000, w)] = arr
        np.savez_compressed(os.environ["FC17_CAP"], **out)
    atexit.register(save)


import numpy as np  # noqa: E402


def register_logdump(c):
    """FC17_LOG=path: at exit, write c.block_log (n, C, cost) in build order (debug)."""
    import os
    import atexit
    p = os.environ.get("FC17_LOG")
    _DUMP["c"] = c
    if not p or _DUMP.get("logreg"):
        return
    _DUMP["logreg"] = True

    def dump():
        cc = _DUMP["c"]
        with open(p, "w") as f:
            for i, e in enumerate(getattr(cc, "block_log", [])):
                f.write("%d %s\n" % (i, repr(e)))
    atexit.register(dump)


def register_nwlog(c):
    """FC17_NWLOG=1: print every narrow block call (n, C, src, tgt, cin) to stderr (debug)."""
    import os
    if not os.environ.get("FC17_NWLOG") or _DUMP.get("nwreg"):
        return
    _DUMP["nwreg"] = True
    NW = sys.modules["invsqrt_fc_narrow"]
    orig = NW.block

    def block(c_, src, tgt, cin):
        print("FC17NW n %d C %d src %s tgt %s cin %s" % (len(tgt), c_.cap - c_.live, list(src), list(tgt), cin),
              file=sys.stderr)
        return orig(c_, src, tgt, cin)
    NW.block = block


def install_uret_cap(c):
    """FC17_CAP + FC17_URET=1: markers around the undo-only re-creation of digit 21 (F1's uret block, tag ('w22', 21)),
    snapshotting EVERY wire before (tag 5) and after (tag 6): exact specs for a shorter re-creation (debug)."""
    import os
    if not (os.environ.get("FC17_CAP") and os.environ.get("FC17_URET")) or _CAP.get("uret"):
        return
    _CAP["uret"] = True
    W = sys.modules["invsqrt_fk_w24ret"]
    orig = W._uret_retire

    def wrapped(F8, K, c_, wire, prog, tag, frugal=False):
        nb = len(c_.blocks)
        res = orig(F8, K, c_, wire, prog, tag, frugal)
        if tuple(tag) == ("w22", 21):
            for bid in range(nb, len(c_.blocks)):
                ops = c_.blocks.get(bid)
                if isinstance(ops, list) and ops and not any(o[0] == "BBEG" and o[1] in (-9005, -9006) for o in ops):
                    allw = list(range(c_.next))
                    _CAP["marks"][-9005] = allw
                    _CAP["marks"][-9006] = allw
                    c_.blocks[bid] = [("BBEG", -9005)] + list(ops) + [("BBEG", -9006)]
                    print("FC17 uret cap: block %d (%d ops), wire %d, prog inputs %s" % (bid, len(ops), wire,
                                                                                     prog.get("inputs")), file=sys.stderr)
            _install_cap()
        return res
    W._uret_retire = wrapped
