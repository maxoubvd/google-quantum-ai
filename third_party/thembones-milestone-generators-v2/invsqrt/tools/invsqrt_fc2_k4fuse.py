#!/usr/bin/env python3
"""invsqrt_fc2_k4fuse.py -- FC2 fourth holder (2026-11-15): MAIN's first add FUSED with row 25 (cfg "_k4fuse": 1).

Today (w26r) MAIN's first add (win -= D on the window (a, s, s, ..., s), D = 1 + 2v + 2^(m+1)) costs 8 after P1's
pass, and row 25 (an unconditional subtraction: its sign is the constant 0) runs at ONE clean wire in the in-place
Cuccaro / c1 form: 16 forward + 16 undo [M tools/invsqrt_fc2_k4prof.py].

Here the first add is F1's AND ladder (tools/invsqrt_fk_w24ret.py `_ladder_first_add`, F1's design) in an s-free
form, and row 25 runs as FC's narrow block forward (one AND per carry, carries X-measured; FC's
invsqrt_fc_narrow.fwd logic, copied) INTERLEAVED with it, one position behind the ladder. Row 25's two carry wires
are the clean wire K and the window's top E wire, which the ladder never fills: the ladder's top bit t_m is made
in K after row 25's last fresh carry, and K becomes position 28+m (the caller's E[m-1]); the old top wire is freed.
  s-free ladder [P]: t_1 = p AND y_0 (p = a ^ s, y_0 = v_0 ^ NOT s); for i >= 1, t_i = 1 forces p = 1 and y_0 = 1,
  i.e. v_0 = s, so t_(i+1) = t_i AND y_i = t_i AND NOT(v_i ^ v_0): after t_1 nothing reads a or s, and positions
  26 / 27 take their final values (NOT a, NOT v_0 ^ p) before row 25 reaches them.
Forward: m ANDs (ladder) + 9 (row 25: 8 carries + the top CCX). Undo (explicit block): row 25's undo = add D with the
carries' debts repaid (R_26 = R_25 + D; FC's block chooser at C = cap - live), then the ladder's mirror (its ANDs
erased by measurement, 0 Toffolis).
"""
import sys


def fused_top(c, cfg, A, E, win, w25, dsrc, v, m):
    NW = sys.modules["invsqrt_fc_narrow"]
    K_ = sys.modules["invsqrt_fc_kit"]
    assert len(win) == m + 2 and list(win[2:]) == list(E) and len(E) == m and len(v) == m
    a, s = win[0], win[1]
    E0 = list(E)                                   # the window's E wires as they are now (position 28 + j)
    src = list(dsrc) + [None]                      # row 25: D = [ONE] + v + [ONE], top source bit 0
    n = m + 3
    assert len(src) == n
    K = c.alloc(1)[0]                              # the clean wire (row 25's carries, then the ladder's top bit)
    if not hasattr(c, "blocks"):
        c.blocks = {}
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    lad = []                                       # the ladder's ops, in order (their mirror is the ladder's undo)

    def rec(fn):
        i0 = len(c.ops)
        fn()
        lad.extend(c.ops[i0:])

    # ---- ladder prologue: copies of s cleared, p = a ^ s, t_1, positions 26 / 27 final ----
    def lad0():
        for e in E0:
            c.cx(s, e)                             # E -> 0
        c.cx(a, s)                                 # s-wire = p
        c.cx(s, v[0])                              # v_0 -> y_0 = v_0 ^ p ^ a ^ 1 = v_0 ^ NOT s
        c.cx(a, v[0])
        c.x(v[0])
        c.and_compute(s, v[0], E0[0])              # t_1 into E[0]
        c.x(v[0])
        c.cx(a, v[0])
        c.cx(s, v[0])
        c.x(a)                                     # position 26 = R_26 bit 0 = NOT a
        c.cx(v[0], s)                              # position 27 = R_26 bit 1 = NOT v_0 ^ p
        c.x(s)
    rec(lad0)
    tdone = {1: E0[0]}                             # t_j -> wire

    def T(j):                                      # t_(j+1) = t_j AND NOT(v_j ^ v_0) into E[j] (j < m-1) or K (j = m-1)
        dst = E0[j] if j < m - 1 else K

        def f():
            c.cx(v[0], v[j])
            c.x(v[j])
            c.and_compute(tdone[j], v[j], dst)
            c.x(v[j])
            c.cx(v[0], v[j])
        rec(f)
        tdone[j + 1] = dst

    def F(j):                                      # E[j] = t_(j+1) -> R_26 bit j+2 = NOT v_(j+1) ^ t_(j+1)
        def f():
            c.cx(v[j + 1], E0[j])
            c.x(E0[j])
        rec(f)

    tgt = [w25, a, s] + E0[:m - 1] + [K]           # row 25's targets, positions 25..35 (35 = K at the end)
    framed = set()

    def frame(k):                                  # row 25 subtracts: its targets are complemented around the add
        if k not in framed:
            c.x(tgt[k])
            framed.add(k)

    def ensure(i):                                 # positions 25..25+i final (R_26) and complemented
        if i >= 3:
            j = i - 3                              # position 28 + j = E[j]
            if (j + 2) not in tdone:
                T(j + 1)
            F(j)
        frame(i)
        if i == n - 2:                             # the top step also writes position 35 (the ladder's t_m in K)
            if m not in tdone:
                T(m - 1)
            frame(n - 1)

    # ---- row 25's block forward (FC's NW.fwd, copied; carries alternate between K and the old top E wire) ----
    # m fresh carries (positions 1..m); the last must sit on the old top wire (K takes t_m at the top step), so the
    # alternation starts on K for even m (m 8, 6: as before) and on the old top wire for odd m (m 7) [FC2 fifth]
    pool = [K, E0[m - 1]] if m % 2 == 0 else [E0[m - 1], K]
    nfresh_ = [0]

    def fresh():
        w = pool[nfresh_[0] % 2]
        nfresh_[0] += 1
        return w
    debts = {}
    cur = NW._lit(None)
    cur_fresh = False
    deferred = None
    for i in range(n - 1):
        if i == n - 2:
            assert cur_fresh and cur[1] == E0[m - 1], "row 25's last carry must sit on the old top wire"
        ensure(i)
        a_, bw = src[i], tgt[i]
        f = NW._form_b(a_, bw, cur)
        nxt, nfresh = None, False
        if i < n - 2:
            if f[0] == "const":
                nxt = ("k", f[1])
            elif f[0] == "alias":
                nxt = ("w", f[1], f[2])
            else:
                _, xw, yw, ops = f
                NW._emit(c, ops)
                t = fresh()
                c.and_compute(xw, yw, t)
                NW._emit(c, ops[::-1])
                NW._xor_lit(c, cur, t)
                nxt, nfresh = ("w", t, 0), True
        else:
            top = tgt[n - 1]
            if f[0] == "const":
                if f[1]:
                    c.x(top)
            elif f[0] == "alias":
                c.cx(f[1], top)
                if f[2]:
                    c.x(top)
            else:
                _, xw, yw, ops = f
                NW._emit(c, ops)
                c.ccx(xw, yw, top)
                NW._emit(c, ops[::-1])
                NW._xor_lit(c, cur, top)
            NW._xor_src(c, src[n - 1], top)
        if f[0] == "alias" and i < n - 2:
            new_def = i
        else:
            new_def = None
            NW._xor_src(c, a_, bw)
            NW._xor_lit(c, cur, bw)
        if cur_fresh:
            debts[i] = c.hmr(cur[1])               # (the wire is reset to 0; it stays in this op's pool)
        elif deferred is not None:
            c.x(tgt[deferred])
        deferred = new_def
        cur, cur_fresh = nxt, nfresh
    assert deferred is None and not cur_fresh
    for k in range(n):
        c.x(tgt[k])                                # row 25's closing complement: R_25 = R26' - D
    c.ops.append(("BEND", bid))
    old_top = E0[m - 1]
    E[m - 1] = K                                   # position 28+m is K from here on
    c.free_([old_top])
    if cfg.get("_r24lazy"):                        # F1's lazy row 24 (tools/invsqrt_fk_w24ret.py): its hook takes row
        c._fk_r24 = dict(ctl=K, v=list(v), m=m)    # 24's add when this hand-off is set (F1's lazy row 25 sets it)
        if cfg.get("_k4r24u"):                     # FC2: row 24's undo on digit 23's wire (r24_lazy_block_k4 below)
            c._k4r24u = True
            _install_r24u()
    # ---- explicit undo: R_26 = R_25 + D (the carries of R_25 + D are the forward's: debts repaid), then the
    # ladder's mirror ----
    save = c.ops
    c.ops = []
    if int(cfg.get("_k4fuse", 1)) >= 2:
        # "_k4fuse": 2 -- row 25's undo at C 2 on digit 24's wire. In the undo, digit 24 (= R_25's top bit, on K)
        # has its last reader in row 24's undo; it is X-measured (bit ms) and K is a carry wire of row 25's undo.
        # The undo's top step only WROTE the top bit (top ^= c'_10 = MAJ(1, R_25[9], c'_9)): its writes become
        # phases on ms (X -> NEG, CX(a) -> Z(a), CCX(a, b) -> CZ(a, b); Clifford, 0 Toffolis). With
        # digit 24 = t_8 ^ c'_10 (R_26's top bit is t_8), the rest of the debt is (-1)^(ms t_8) = CZ(t_7, y'_7) if
        # ms, placed where the ladder's mirror would erase t_8 (K never holds t_8 now).
        ms = c.hmr(K)
        c.free_([K])
        SENT = -1000001                                 # the top target: never a wire
        tfs = tgt[:n - 1] + [SENT]
        C = c.cap - c.live
        assert C >= 2, ("row 25's undo on digit 24's wire needs C 2", C, c.live, c.cap)
        i0 = len(c.ops)
        if NW.plan(max(0, n - 2 - NW.free_prefix(src, n, None)), C)[0] is not None:   # FC's block chooser
            cost = NW.game_add(c, src, tfs, None, debts, C)
        elif C > NW.cuccaro_need(src, tfs, None):
            cost = float(NW.hybrid_add(c, src, tfs, None, debts, C))   # (C1UNDO: F3's mixed / c1 forms)
        else:
            cost = float(NW.cuccaro_add(c, src, tfs, None, debts))
        body = c.ops[i0:]
        del c.ops[i0:]
        def cnd(b):                                     # the phase's condition: ms, or ms AND the write's own bit
            return ms if b is None else c.band(((ms, 1), (b, 1)))
        for o in body:
            g = o[0]
            if g == "X" and o[1] == SENT:
                c.neg(cnd(o[2]))
            elif g == "CX" and o[2] == SENT:
                c.z(o[1], cnd(o[3]))
            elif g == "CCX" and o[3] == SENT:
                c.cz(o[1], o[2], ms)
                cost -= 1.0
            elif g == "XCCX" and o[3] == SENT:          # F3's gated CCX (fcgen): t ^= a b if bit g
                c.cz(o[1], o[2], cnd(o[4]))
            else:
                assert SENT not in o[1:], o
                c.ops.append(o)
        assert K in c.free
        c.free.remove(K)                                # K is position 35 again for the rest of the forward
        c.live += 1
        mir = K_.Circ.mirror(lad)
        hit = [k_ for k_, o in enumerate(mir) if o[0] == "ANDE" and o[3] == K]
        assert len(hit) == 1, hit
        o = mir[hit[0]]
        mir[hit[0]] = ("CZ", o[1], o[2], ms)
        c.ops.extend(mir)
    else:
        C = c.cap - c.live
        if NW.plan(max(0, n - 2 - NW.free_prefix(src, n, None)), C)[0] is not None:
            cost = NW.game_add(c, src, tgt, None, debts, C)
        elif C > NW.cuccaro_need(src, tgt, None):
            cost = float(NW.hybrid_add(c, src, tgt, None, debts, C))
        else:
            cost = float(NW.cuccaro_add(c, src, tgt, None, debts))
        c.ops.extend(K_.Circ.mirror(lad))
    c.blocks[bid] = c.ops
    c.ops = save
    if not hasattr(c, "block_log"):
        c.block_log = []
    c.block_log.append((n, C, cost))
    LOG.append(dict(C=C, undo_cost=cost, debts=len(debts)))


LOG = []
_ORIG = {}


def _phase_rewrite(c, body, SENT, ms):
    """append `body` to c.ops with every write to SENT turned into a phase on bit ms: X -> NEG, CX(a) -> Z(a),
    CCX(a, b) -> CZ(a, b), gated CCX (XCCX, F3 / fcgen) -> CZ under (ms AND gate). Returns the Toffolis removed."""
    gone = 0.0

    def cnd(b):
        return ms if b is None else c.band(((ms, 1), (b, 1)))
    for o in body:
        g = o[0]
        if g == "X" and o[1] == SENT:
            c.neg(cnd(o[2]))
        elif g == "CX" and o[2] == SENT:
            c.z(o[1], cnd(o[3]))
        elif g == "CCX" and o[3] == SENT:
            c.cz(o[1], o[2], ms)
            gone += 1.0
        elif g == "XCCX" and o[3] == SENT:
            c.cz(o[1], o[2], cnd(o[4]))
        else:
            assert SENT not in o[1:], o
            c.ops.append(o)
    return gone


def _install_r24u():
    """route F1's lazy row-24 block (invsqrt_fk_w24ret._r24_lazy_block, F1's code, read-only) through
    r24_lazy_block_k4 for builds that set c._k4r24u (others keep F1's function)."""
    W = sys.modules["invsqrt_fk_w24ret"]
    if "r24" in _ORIG:
        return
    _ORIG["r24"] = W._r24_lazy_block

    def route(c, src, tgt, info):
        if getattr(c, "_k4r24u", False):
            return r24_lazy_block_k4(c, src, tgt, info)
        return _ORIG["r24"](c, src, tgt, info)
    W._r24_lazy_block = route


def r24_lazy_block_k4(c, src, tgt, info):
    """F1's lazy row 24 (forward copied from invsqrt_fk_w24ret._r24_lazy_block, F1's design) with FC2's undo: in the
    undo, row 24's top target holds digit 23 (R_24's sign, framed), whose last reader is row 23's undo; it is
    X-measured (bit md) and its wire is the undo add's second carry wire (C 2). The undo add's writes to the top
    become phases on md; then the top is re-created as F1's framed R_25 bit (1 AND, F1's _r24_form) and Z(top) if md
    repays the rest: (-1)^(md (T''_top ^ W)) = (-1)^(md T'_top)."""
    W = sys.modules["invsqrt_fk_w24ret"]
    NW = sys.modules["invsqrt_fc_narrow"]
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    top, d, v = tgt[-1], info["ctl"], info["v"]
    if not hasattr(c, "blocks"):
        c.blocks = {}
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    W._r24_form(c, tgt, d, v, top, erase=True)
    c.free.insert(0, top)
    c.live -= 1

    def loader(c_, tg, tp, lz):
        W._r24_form(c_, tg, d, v, tp, erase=False)
    debts = W._fwd_lazy(c, NW, src, tgt, top, info, loader=loader)
    c.ops.append(("BEND", bid))
    assert top not in c.free, "lazy top (row 24): the top wire was not taken"
    save = c.ops
    c.ops = []
    md = c.hmr(top)
    c.free_([top])
    C = c.cap - c.live
    assert C >= 2, ("row 24's undo on digit 23's wire needs C 2", C)
    SENT = -1000002
    tg2 = list(tgt[:-1]) + [SENT]
    i0 = len(c.ops)
    for w in tg2:
        c.x(w)
    if NW.plan(max(0, n - 2 - NW.free_prefix(src, n, None)), C)[0] is not None:
        cost = NW.game_add(c, src, tg2, None, debts, C)
    elif C > NW.cuccaro_need(src, tg2, None):
        cost = float(NW.hybrid_add(c, src, tg2, None, debts, C))
    else:
        cost = float(NW.cuccaro_add(c, src, tg2, None, debts))
    for w in tg2:
        c.x(w)
    body = c.ops[i0:]
    del c.ops[i0:]
    cost -= _phase_rewrite(c, body, SENT, md)
    assert top in c.free
    c.free.remove(top)
    c.live += 1
    W._r24_form(c, tgt, d, v, top, erase=False)    # top = F1's framed R_25 bit (1 AND)
    c.z(top, md)
    cost += 1.0
    c.blocks[bid] = c.ops
    c.ops = save
    if not hasattr(c, "block_log"):
        c.block_log = []
    c.block_log.append((n, C, cost))
    LOG.append(dict(r24u=int(top), C_undo=C, undo_cost=cost))
