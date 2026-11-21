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
    if cfg.get("_r63chain"):                       # R63-P1 (r63_inv_row9_p1.py): the first block as one chain
        return _r63_chain_top(c, cfg, E, win, w25, v, m)
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



# R63-P1 (r63_inv_row9_p1.py) ------------------------------------------------------------------------------------
# (target bit, linear terms, left factor, right factor); "w" = the wire below, an int j = v_j, "1" = the constant
_R63_CHAIN = [
    (4, [1, 3, "w"], [0, 2], [2, "w", "1"]),
    (5, [2, 4, "w"], [3, "w", "1"], [1, 3]),
    (6, [3, 5, "w"], [2, "w"], [4, "w", "1"]),
    (7, [3, 4, 5, 6, "w"], [3, "w", "1"], [3, 5]),
    (8, [5, "w", "1"], [4, 6], [4, "w"]),
    (9, [6], ["w"], [5, "w", "1"]),
]


def _r63_chain_top(c, cfg, E, win, w25, v, m):
    """MAIN's first add + row 25 as ONE CHAIN: R25 = 253 - 6 v + u' on [w25, a, s] + E (10 wires), 8 ANDs forward,
    explicit undo = the mirror (1 AND).  Same hand-off as fused_top (row 24's control = the window's top wire)."""
    K_ = sys.modules["invsqrt_fc_kit"]
    assert m == 7 and len(win) == m + 2 and list(win[2:]) == list(E) and len(E) == m and len(v) == m
    a, s = win[0], win[1]
    E0 = list(E)
    top = E0[m - 1]
    if not hasattr(c, "blocks"):
        c.blocks = {}
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    i0 = len(c.ops)
    for e in E0:
        c.cx(s, e)                                 # the sign copies: E -> 0
    tmp = top                                      # c2's temporary (zero again before w9 takes the wire)
    # c2 = MAJ(v0, a, w25) = w25 ^ (v0 ^ w25)(a ^ w25)
    c.cx(w25, v[0])
    c.cx(w25, a)
    c.and_compute(v[0], a, tmp)
    c.cx(w25, a)
    c.cx(w25, v[0])
    c.cx(w25, tmp)
    # c3 = MAJ(1 ^ v1, s, c2) = c2 ^ (1 ^ v1 ^ c2)(s ^ c2), into the wire of bit 3
    c.cx(tmp, s)
    c.cx(tmp, v[1])
    c.x(v[1])
    c.and_compute(v[1], s, E0[0])
    c.x(v[1])
    c.cx(tmp, v[1])
    c.cx(tmp, s)
    c.cx(tmp, E0[0])
    # n3 = 1 ^ v0 ^ v1 ^ v2 ^ s ^ c3 ; n2 = 1 ^ v1 ^ s ^ c2 (in place on s)
    for w in (s, v[0], v[1], v[2]):
        c.cx(w, E0[0])
    c.x(E0[0])
    c.cx(tmp, s)
    c.cx(v[1], s)
    c.x(s)
    # c2 erased by measurement (a, w25 still hold U's bits)
    c.cx(w25, tmp)
    c.cx(w25, v[0])
    c.cx(w25, a)
    c.ops.append(("ANDE", v[0], a, tmp))
    c.cx(w25, a)
    c.cx(w25, v[0])
    # n1 = v0 ^ a ^ w25 (in place on a) ; n0 = w25 ^ 1
    c.cx(v[0], a)
    c.cx(w25, a)
    c.x(w25)
    # the chain: bits 4..9, each one AND straight into its own wire
    wire = {3: E0[0], 4: E0[1], 5: E0[2], 6: E0[3], 7: E0[4], 8: E0[5], 9: top}
    _m18 = {}                                    # FC18-R25C: where each bit's forward ops start
    for bit, lt, lf, rf in _R63_CHAIN:
        _m18[bit] = len(c.ops)
        below, tgt = wire[bit - 1], wire[bit]

        def wires_of(form):
            return [below if t == "w" else v[t] for t in form if t != "1"]

        def fold(form, host):
            """the ops that XOR the form onto its wire `host`."""
            ops = [("CX", w_, host, None) for w_ in wires_of(form) if w_ != host]
            if "1" in form:
                ops.append(("X", host, None))
            return ops
        wl, wr = wires_of(lf), wires_of(rf)
        pick = None
        for rank in (0, 1):                        # 0: no host is read by the other factor; 1: one is (ordered folds)
            for hl in wl:
                for hr in wr:
                    if hl == hr or pick is not None:
                        continue
                    a_in, b_in = hl in wr, hr in wl
                    if rank == 0 and not a_in and not b_in:
                        pick = (hl, hr, "LR")
                    elif rank == 1 and a_in and not b_in:
                        pick = (hl, hr, "RL")       # the right factor reads hl: fold it first
                    elif rank == 1 and b_in and not a_in:
                        pick = (hl, hr, "LR")       # the left factor reads hr: fold it first
        assert pick is not None, (bit, "no pair of fold wires")
        hl, hr, order = pick
        ol, o_r = fold(lf, hl), fold(rf, hr)
        first, second = (ol, o_r) if order == "LR" else (o_r, ol)
        c.ops.extend(first)
        c.ops.extend(second)
        c.and_compute(hl, hr, tgt)
        c.ops.extend(second[::-1])
        c.ops.extend(first[::-1])
        for t in lt:
            if t == "1":
                c.x(tgt)
            else:
                c.cx(below if t == "w" else v[t], tgt)
    fwd = c.ops[i0:]
    c.ops.append(("BEND", bid))
    c.blocks[bid] = K_.Circ.mirror(fwd)
    if cfg.get("_fc18r25c"):                     # FC18-R25C (invsqrt_fc18_r25c.py): bits 7, 8 under g only
        import invsqrt_fc18_r25c as _R25C
        c.blocks[bid] = _R25C.chain_mirror(c, K_, c.ops[i0:_m18[7]], c.ops[_m18[7]:_m18[8]],
                                           c.ops[_m18[8]:_m18[9]], c.ops[_m18[9]:i0 + len(fwd)], wire)
    cost = float(c.toffolis_exp(c.blocks[bid]))
    if cfg.get("_r24lazy"):                        # F1's lazy row 24: its control is the window's top wire
        c._fk_r24 = dict(ctl=top, v=list(v), m=m)
        if cfg.get("_k4r24u"):
            c._k4r24u = True
            _install_r24u()
    if not hasattr(c, "block_log"):
        c.block_log = []
    c.block_log.append((m + 3, c.cap - c.live, cost))
    LOG.append(dict(r63chain=1, C=c.cap - c.live, undo_cost=cost, debts=0))


# CH1 P6 (tools/inv_ch1_p6.py) -------------------------------------------------------------------------------------
def _ch1_and_into(c, left, right, tgt, erase=False):
    """tgt ^= AND(XOR left, XOR right) (forms: lists of wires and the string "1"), folding each form onto one of its
    own wires that the other form does not read; erase=True: tgt holds exactly that AND -> erased by measurement."""
    lw = [w for w in left if w != "1"]
    rw = [w for w in right if w != "1"]
    hl = [w for w in lw if w not in rw]
    hr = [w for w in rw if w not in lw]
    assert hl and hr, ("no host wires", left, right)
    hl, hr = hl[0], hr[0]
    ops = []
    for w in lw:
        if w != hl:
            ops.append(("CX", w, hl, None))
    if "1" in left:
        ops.append(("X", hl, None))
    for w in rw:
        if w != hr:
            ops.append(("CX", w, hr, None))
    if "1" in right:
        ops.append(("X", hr, None))
    c.ops.extend(ops)
    if erase:
        c.ops.append(("ANDE", hl, hr, tgt))
    else:
        c.ccx(hl, hr, tgt)
    c.ops.extend(ops[::-1])


def _ch1_w7(c, tgt, d, v, t8, erase):
    """w7' = v3^v4^v5^v6^w6' ^ (v3^d^w6')(v3^v5) on the wire t8 (w6' = tgt[7]); erase: t8 holds w7' -> 0."""
    w6 = tgt[7]
    lin = [v[3], v[4], v[5], v[6], w6]
    if erase:
        for w in lin:
            c.cx(w, t8)
        _ch1_and_into(c, [v[3], d, w6], [v[3], v[5]], t8, erase=True)
    else:
        _ch1_and_into(c, [v[3], d, w6], [v[3], v[5]], t8)
        for w in lin:
            c.cx(w, t8)


def _ch1_o9(c, tgt, d, v, top, erase):
    """o9' = v6 ^ (v4^v6^d^o8'^o7')(v4^v5^v6^d^1) ^ (v6^o8')(v6^o7') on the wire top (o8' = tgt[8], o7' = tgt[7]);
    born: 2 ANDs; erase: 1 AND + a measured erase."""
    o7, o8 = tgt[7], tgt[8]
    A = ([o7, v[4], v[6], d, o8], [v[5], v[4], v[6], d, "1"])
    B = ([o8, v[6]], [o7, v[6]])
    if erase:
        c.cx(v[6], top)
        _ch1_and_into(c, B[0], B[1], top)
        _ch1_and_into(c, A[0], A[1], top, erase=True)
    else:
        _ch1_and_into(c, A[0], A[1], top)
        _ch1_and_into(c, B[0], B[1], top)
        c.cx(v[6], top)


def _ch1_r24_block38(c, src, tgt, info):
    W = sys.modules["invsqrt_fk_w24ret"]
    NW = sys.modules["invsqrt_fc_narrow"]
    n = len(tgt)
    assert n == 10, n
    src = list(src) + [None] * (n - len(src))
    top, t8, d, v = tgt[9], tgt[8], info["ctl"], info["v"]
    assert len(v) == 7
    if not hasattr(c, "blocks"):
        c.blocks = {}
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    W._r24_form(c, tgt, d, v, top, erase=True)     # w8' -> 0 (reads w7')
    c.free.insert(0, top)
    c.live -= 1
    _ch1_w7(c, tgt, d, v, t8, erase=True)          # w7' -> 0 (reads w6')
    c.free.insert(0, t8)
    c.live -= 1

    def loader(c_, tg, tp, lz):
        _ch1_w7(c_, tgt, d, v, tp, erase=False)
    debts = W._fwd_lazy(c, NW, src[:9], tgt[:9], t8, info, loader=loader)
    assert t8 not in c.free, "row 24 at 38: w7 was not taken"
    c.free.remove(top)
    c.live += 1
    c.peak = max(c.peak, c.live)
    _ch1_o9(c, tgt, d, v, top, erase=False)        # digit 23 (framed) born by its flat form
    c.ops.append(("BEND", bid))
    save = c.ops
    c.ops = []
    _ch1_o9(c, tgt, d, v, top, erase=True)         # undo: the top erased by its form
    c.free.insert(0, top)
    c.live -= 1
    m8 = c.hmr(t8)
    c.free.insert(0, t8)
    c.live -= 1
    C = c.cap - c.live
    assert C >= 2, ("row 24 at 38: undo needs two clean wires", C)
    SENT = -1000003
    tg2 = list(tgt[:8]) + [SENT]
    s9 = src[:9]
    i0 = len(c.ops)
    for w in tg2:
        c.x(w)
    if NW.plan(max(0, 9 - 2 - NW.free_prefix(s9, 9, None)), C)[0] is not None:
        cost = NW.game_add(c, s9, tg2, None, debts, C)
    elif C > NW.cuccaro_need(s9, tg2, None):
        cost = float(NW.hybrid_add(c, s9, tg2, None, debts, C))
    else:
        cost = float(NW.cuccaro_add(c, s9, tg2, None, debts))
    for w in tg2:
        c.x(w)
    body = c.ops[i0:]
    del c.ops[i0:]
    cost -= _phase_rewrite(c, body, SENT, m8)
    c.free.remove(t8)
    c.live += 1
    _ch1_w7(c, tgt, d, v, t8, erase=False)         # w7' re-made, then the rest of S_8's debt
    c.z(t8, m8)
    c.free.remove(top)
    c.live += 1
    W._r24_form(c, tgt, d, v, top, erase=False)    # w8' re-made
    cost += 3.0
    c.blocks[bid] = c.ops
    c.ops = save
    if not hasattr(c, "block_log"):
        c.block_log = []
    c.block_log.append((n, C, cost))
    LOG.append(dict(ch1r24=int(top), C_undo=C, undo_cost=cost))


# CH1 P7 (tools/inv_ch1_p7.py) -------------------------------------------------------------------------------------
def _ch1_prod(c, left, right, tgt, mode):
    """one product of two XOR forms (wires and "1"): the left form is folded onto a wire of its own that the right
    form does not read, the right form onto its first wire; mode 'and' (tgt |0>), 'ccx' (tgt ^= product), 'ande'
    (tgt == product: measured erase), ('cz', bit) (the phase (-1)^(bit * product))."""
    lw = [w for w in left if w != "1"]
    rw = [w for w in right if w != "1"]
    hl = [w for w in lw if w not in rw]
    assert hl, ("no host for the left form", left, right)
    hl = hl[0]
    hr = [w for w in rw if w not in lw] or rw
    hr = hr[0]
    ul = [("CX", w, hl, None) for w in lw if w != hl] + ([("X", hl, None)] if "1" in left else [])
    ur = [("CX", w, hr, None) for w in rw if w != hr] + ([("X", hr, None)] if "1" in right else [])
    c.ops.extend(ul)
    c.ops.extend(ur)
    if mode == "and":
        c.and_compute(hl, hr, tgt)
    elif mode == "ccx":
        c.ccx(hl, hr, tgt)
    elif mode == "ande":
        c.ops.append(("ANDE", hl, hr, tgt))
    else:
        c.cz(hl, hr, mode[1])
    c.ops.extend(ur[::-1])
    c.ops.extend(ul[::-1])


def _ch1_r23_top(c, d, v, o8, top, mode):
    """digit 22 (framed) o9' = 1^d^g0^g1^g2^o8'^v3^v4^v5 with g0 = (o8'^v3^v5^v6)(1^d^v4^v6),
    g1 = (1^v3^v5^v6)(d^g0^o8'), g2 = (1^o8'^v5)(g0^o8'^v6).  mode 'make': into the |0> wire top (3 ANDs);
    ('phase', md): repay (-1)^(md o9') of a measured top (1 AND: g0; the rest Z / CZ under md)."""
    g0 = c.alloc(1)[0]
    F0 = ([o8, v[3], v[5], v[6]], ["1", d, v[4], v[6]])
    F1 = (["1", v[3], v[5], v[6]], [d, g0, o8])
    F2 = (["1", v[5], o8], [g0, o8, v[6]])
    _ch1_prod(c, F0[0], F0[1], g0, "and")
    lin = (d, g0, o8, v[3], v[4], v[5])
    if mode == "make":
        _ch1_prod(c, F1[0], F1[1], top, "ccx")
        _ch1_prod(c, F2[0], F2[1], top, "ccx")
        for w in lin:
            c.cx(w, top)
        c.x(top)
    else:
        md = mode[1]
        _ch1_prod(c, F1[0], F1[1], None, ("cz", md))
        _ch1_prod(c, F2[0], F2[1], None, ("cz", md))
        for w in lin:
            c.z(w, md)
        c.neg(md)
    _ch1_prod(c, F0[0], F0[1], g0, "ande")
    c.free_([g0])


def _ch1_r23_w8(c, d, v, w7, top):
    """R24 bit 8 (framed) into the |0> wire top: w8' = 1^d^h1^h2^v4^v5 with h0 = (v2^v3^v6^w7')(d^v2^v5^w7'),
    h1 = h0 (v2^v3^v5^w7'), h2 = (d^h0^v4)(1^d^v4^w7'): 3 ANDs, h0 on a helper (erased by measurement)."""
    h0 = c.alloc(1)[0]
    H0 = ([v[3], v[2], v[6], w7], [d, v[2], v[5], w7])
    _ch1_prod(c, H0[0], H0[1], h0, "and")
    _ch1_prod(c, [h0], [v[2], v[3], v[5], w7], top, "ccx")
    _ch1_prod(c, [h0, d, v[4]], ["1", d, v[4], w7], top, "ccx")
    for w in (d, v[4], v[5]):
        c.cx(w, top)
    c.x(top)
    _ch1_prod(c, H0[0], H0[1], h0, "ande")
    c.free_([h0])


def _ch1_r23_block38(c, d, src, tgt, v):
    NW = sys.modules["invsqrt_fc_narrow"]
    n = len(tgt)
    assert n == 10, n
    src = list(src) + [None] * (n - len(src))
    assert src[9] is None
    top, t8 = tgt[9], tgt[8]
    if not hasattr(c, "blocks"):
        c.blocks = {}
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    mt = c.hmr(top)                                # R24 bit 8 (framed) measured: debt (-1)^(mt w8')
    c.free.insert(0, top)
    c.live -= 1
    debts = NW.fwd(c, src[:9], tgt[:9], None)
    c.free.remove(top)
    c.live += 1
    c.peak = max(c.peak, c.live)
    _ch1_r23_top(c, d, v, t8, top, "make")         # digit 22 (framed) born by its form
    c.ops.append(("BEND", bid))
    save = c.ops
    c.ops = []
    md = c.hmr(top)                                # undo: digit 22 measured, its phase repaid through the form
    c.free.insert(0, top)
    c.live -= 1
    _ch1_r23_top(c, d, v, t8, top, ("phase", md))
    C = c.cap - c.live
    assert C >= 2, ("row 23 at 38: undo needs two clean wires", C)
    s9, t9 = src[:9], list(tgt[:9])
    for w in t9:
        c.x(w)
    if NW.plan(max(0, 9 - 2 - NW.free_prefix(s9, 9, None)), C)[0] is not None:
        cost = NW.game_add(c, s9, t9, None, debts, C)
    elif C > NW.cuccaro_need(s9, t9, None):
        cost = float(NW.hybrid_add(c, s9, t9, None, debts, C))
    else:
        cost = float(NW.cuccaro_add(c, s9, t9, None, debts))
    for w in t9:
        c.x(w)
    c.free.remove(top)
    c.live += 1
    _ch1_r23_w8(c, d, v, t8, top)                  # R24 bit 8 (framed) re-made, then the forward's debt
    c.z(top, mt)
    cost += 4.0
    c.blocks[bid] = c.ops
    c.ops = save
    if not hasattr(c, "block_log"):
        c.block_log = []
    c.block_log.append((n, C, cost))
    LOG.append(dict(ch1r23=int(top), C_undo=C, undo_cost=cost))

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
        _cfg16 = W._STATE.get("cfg") or {}                       # FC16-R24N2: rows 24 / 23 at one wire less
        _cap16 = _cfg16.get("_fc16r24n")
        if _cap16:
            import invsqrt_fc16_r24n as _N16
            return _N16.block(c, src, tgt, info, int(_cap16))
        if getattr(c, "_k4r24u", False):
            _r16 = r24_lazy_block_k4(c, src, tgt, info)
        else:
            _r16 = _ORIG["r24"](c, src, tgt, info)
        if _cfg16.get("_fc16r23n"):                              # FC's row 23 after any row-24 block
            import invsqrt_fc16_r24n as _N16
            _N16.arm23(c, tgt, info, int(_cfg16["_fc16r23n"]))
        return _r16
    W._r24_lazy_block = route


def r24_lazy_block_k4(c, src, tgt, info):
    """F1's lazy row 24 (forward copied from invsqrt_fk_w24ret._r24_lazy_block, F1's design) with FC2's undo: in the
    undo, row 24's top target holds digit 23 (R_24's sign, framed), whose last reader is row 23's undo; it is
    X-measured (bit md) and its wire is the undo add's second carry wire (C 2). The undo add's writes to the top
    become phases on md; then the top is re-created as F1's framed R_25 bit (1 AND, F1's _r24_form) and Z(top) if md
    repays the rest: (-1)^(md (T''_top ^ W)) = (-1)^(md T'_top)."""
    if (sys.modules["invsqrt_fk_w24ret"]._STATE.get("cfg") or {}).get("_ch1r24"):   # CH1 P6
        if (sys.modules["invsqrt_fk_w24ret"]._STATE.get("cfg") or {}).get("_fc17skipe"):   # FC17: skip E
            import invsqrt_fc17_skipe as _S17
            return _S17.r24_block(c, src, tgt, info)
        return _ch1_r24_block38(c, src, tgt, info)
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
    assert C >= 2 or __import__("os").environ.get("R64_K4C1"), ("row 24's undo on digit 23's wire needs C 2", C)   # R64-INV-STACK probe knob (my snapshot only)
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
