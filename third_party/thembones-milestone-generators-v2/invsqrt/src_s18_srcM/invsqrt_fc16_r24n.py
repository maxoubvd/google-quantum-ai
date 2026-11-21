#!/usr/bin/env python3
"""invsqrt_fc16_r24n.py -- FC sixteenth holder: MAIN's ROW 24 AT ONE WIRE LESS (R63-INV-ROW9's forms, R64-INV-STACK's
open step 1e), for FC's source copies.

Row 24 adds +-d 2^24 into the window R25[24..33] under the control ctl = R25's sign (fadd: the targets are framed,
t' = t xor NOT ctl, inside the add). Today (F1's lazy top + FC2's k4 undo) its forward needs 2 clean wires after its
top input is erased. This block erases BOTH top inputs and makes both top results from short forms, so the whole
block lives at (entry live) at most, i.e. one wire below today's:
  forward  (1) top' (R25 bit 8) erased by F1's form, (2) T' (R25 bit 7) erased by ROW9's chain form
           T' = v3^v4^v5^v6^w6' ^ (v3^w6'^ctl)(v3^v5)   [framed; checked on all inputs, unframed form 0 differ]
           (3) the narrow ripple on positions 0..7 only (FC's NW.fwd: carries 2..6 measured, c7 into position 7)
           (4) o8' (framed R24 bit 8) by ROW9's 3-AND form with one helper (g1 on the helper, g0 + g2 + g1 into
           the wire, helper erased free):  g0 = (v3^v6)(v5^o7'), g1 = (v2^v5^ctl^1)(v2^v3^v6^o7'^ctl),
           g2 = (v3^v5^1)(v5^o7'^ctl^1^g1), o8' = v4^v6^ctl^g0^g1^g2
           (5) o9' (framed R24 bit 9 = digit 23 framed) by ROW9's flat 2-AND form:
           o9' = v6 ^ (v4^v6^ctl^o8'^o7')(v4^v5^v6^ctl^1) ^ (v6^o8')(v6^o7')
  undo     o9' erased by its form (free), o8' erased by its form (g2's phase needs g1: one conditioned AND),
           the narrow undo game on positions 0..7 (FC's NW.game_add, the forward's debts repaid), T' re-made
           (1 AND), top' re-made by F1's form (1 AND).
cfg key "_fc16r24n": the cap the block is planned at (e.g. 39 on a 40-wire build); absent = today's block.
Every live count inside the block is asserted <= that cap. Count-only claims need the base-vs-variant check."""
import sys

LOG = []


def _cx_all(c, srcs, t):
    for s in srcs:
        c.cx(s, t)


def block(c, src, tgt, info, cap):
    W = sys.modules["invsqrt_fk_w24ret"]
    NW = sys.modules["invsqrt_fc_narrow"]
    n = len(tgt)
    assert n == 10 and info["m"] == 7, "r24n: written for m 7 (window of 10)"
    src = list(src) + [None] * (n - len(src))
    d, v = info["ctl"], info["v"]
    v2, v3, v4, v5, v6 = v[2], v[3], v[4], v[5], v[6]
    top, t8, t7 = tgt[9], tgt[8], tgt[7]
    peak = [c.live]

    def chk():
        peak[0] = max(peak[0], c.live)
        assert c.live <= cap, ("r24n: live %d above the cap %d" % (c.live, cap))
    chk()
    if not hasattr(c, "blocks"):
        c.blocks = {}
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    i0 = len(c.ops)
    # (1) top' erased by F1's form (reads T' = tgt[8], still present)
    W._r24_form(c, tgt, d, v, top, erase=True)
    # (2) T' erased by its form: T' ^ (v3^v4^v5^v6^w6') = A B, A = v3^w6'^ctl (on w6'), B = v3^v5 (on v5)
    _cx_all(c, (v3, v4, v5, v6, t7), t8)
    c.cx(v3, t7); c.cx(d, t7)
    c.cx(v3, v5)
    c.ops.append(("ANDE", t7, v5, t8))
    c.cx(v3, v5)
    c.cx(d, t7); c.cx(v3, t7)
    for w in (top, t8):                              # both wires free (taken last)
        c.free.insert(0, w)
        c.live -= 1
    # (3) the narrow ripple on positions 0..7 (carries 2..6 measured: debts; c7 into position 7 by one CCX)
    alloc0 = c.alloc

    def alloc_chk(k=1):
        out = alloc0(k)
        chk()
        return out
    c.alloc = alloc_chk
    try:
        debts = NW.fwd(c, src[:8], list(tgt[:8]), None)
    finally:
        c.alloc = alloc0
    chk()
    # (4) o8' into tgt[8]'s wire; helper = g1
    c.free.remove(t8); c.live += 1
    h = c.alloc(1)[0]
    chk()
    # helper <- g1 = (v2^v5^ctl^1)(v2^v3^v6^o7'^ctl): X1 on v5, Y1 on t7
    c.cx(v2, v5); c.cx(d, v5); c.x(v5)
    _cx_all(c, (v2, v3, v6, d), t7)
    c.and_compute(v5, t7, h)
    _cx_all(c, (d, v6, v3, v2), t7)
    c.x(v5); c.cx(d, v5); c.cx(v2, v5)
    # t8 <- g0 = (v3^v6)(v5^o7'): X0 on v3, Y0 on t7
    c.cx(v6, v3); c.cx(v5, t7)
    c.and_compute(v3, t7, t8)
    c.cx(v5, t7); c.cx(v6, v3)
    # t8 ^= g2 = (v3^v5^1)(v5^o7'^ctl^1^g1): X2 on v3, Y2 on the helper
    c.cx(v5, v3); c.x(v3)
    c.cx(v5, h); c.cx(t7, h); c.cx(d, h); c.x(h)
    c.ccx(v3, h, t8)
    c.x(h); c.cx(d, h); c.cx(t7, h); c.cx(v5, h)
    c.x(v3); c.cx(v5, v3)
    # t8 ^= g1 ^ v4 ^ v6 ^ ctl
    c.cx(h, t8)
    _cx_all(c, (v4, v6, d), t8)
    # helper (g1) erased free: X1 on v5, Y1 on t7
    c.cx(v2, v5); c.cx(d, v5); c.x(v5)
    _cx_all(c, (v2, v3, v6, d), t7)
    c.and_erase(v5, t7, h)
    _cx_all(c, (d, v6, v3, v2), t7)
    c.x(v5); c.cx(d, v5); c.cx(v2, v5)
    # (5) o9' into the top's wire: A1 B1 (A1 = v4^v6^ctl^o8'^o7' on t8, B1 = v4^v5^v6^ctl^1 on v5), then A2 B2
    c.free.remove(top); c.live += 1
    chk()
    _cx_all(c, (v4, v6, d, t7), t8)
    _cx_all(c, (v4, v6, d), v5); c.x(v5)
    c.and_compute(t8, v5, top)
    c.x(v5); _cx_all(c, (d, v6, v4), v5)
    _cx_all(c, (t7, d, v6, v4), t8)
    c.cx(v6, t8); c.cx(v6, t7)                       # A2 = v6^o8' on t8, B2 = v6^o7' on t7
    c.ccx(t8, t7, top)
    c.cx(v6, t7); c.cx(v6, t8)
    c.cx(v6, top)
    fwd_cost = float(c.toffolis_exp(c.ops[i0:]))
    c.ops.append(("BEND", bid))
    # ---------------- explicit undo (state: o0'..o7' on tgt[0..7], o8' on t8, o9' on top)
    save = c.ops
    c.ops = []
    # o9' erased: top ^= v6, HMR, CZ(A1, B1), CZ(A2, B2) under the outcome
    c.cx(v6, top)
    m9 = c.hmr(top)
    c.free_([top])
    _cx_all(c, (v4, v6, d, t7), t8)
    _cx_all(c, (v4, v6, d), v5); c.x(v5)
    c.cz(t8, v5, m9)
    c.x(v5); _cx_all(c, (d, v6, v4), v5)
    _cx_all(c, (t7, d, v6, v4), t8)
    c.cx(v6, t8); c.cx(v6, t7)
    c.cz(t8, t7, m9)
    c.cx(v6, t7); c.cx(v6, t8)
    # o8' erased: t8 ^= v4^v6^ctl -> g0^g1^g2, HMR, then CZ for g0, g1 and (with g1 re-made under the outcome) g2
    _cx_all(c, (v4, v6, d), t8)
    m8 = c.hmr(t8)
    c.free_([t8])
    c.cx(v6, v3); c.cx(v5, t7)
    c.cz(v3, t7, m8)                                 # g0
    c.cx(v5, t7); c.cx(v6, v3)
    h = c.alloc(1)[0]
    chk()
    c.cx(v2, v5); c.cx(d, v5); c.x(v5)
    _cx_all(c, (v2, v3, v6, d), t7)
    c.cz(v5, t7, m8)                                 # g1
    c.cand_compute(v5, t7, h, m8)                    # h = g1 if m8
    _cx_all(c, (d, v6, v3, v2), t7)
    c.x(v5); c.cx(d, v5); c.cx(v2, v5)
    c.cx(v5, v3); c.x(v3)
    c.cx(v5, h); c.cx(t7, h); c.cx(d, h); c.x(h)
    c.cz(v3, h, m8)                                  # g2 (h holds Y2 when m8 = 1)
    c.x(h); c.cx(d, h); c.cx(t7, h); c.cx(v5, h)
    c.x(v3); c.cx(v5, v3)
    c.cx(v2, v5); c.cx(d, v5); c.x(v5)
    _cx_all(c, (v2, v3, v6, d), t7)
    c.cand_erase(v5, t7, h, m8)
    _cx_all(c, (d, v6, v3, v2), t7)
    c.x(v5); c.cx(d, v5); c.cx(v2, v5)
    c.free_([h])
    # the narrow undo game on positions 0..7 (X^8 . add . X^8), the forward's debts repaid
    C = cap - c.live
    low = list(tgt[:8])
    for w in low:
        c.x(w)
    a0 = c.alloc
    c.alloc = alloc_chk
    try:
        if NW.plan(max(0, 8 - 2 - NW.free_prefix(src[:8], 8, None)), C)[0] is not None:
            gcost = float(NW.game_add(c, src[:8], low, None, debts, C))
            kind = "game"
        elif C > NW.cuccaro_need(src[:8], low, None):
            gcost = float(NW.hybrid_add(c, src[:8], low, None, debts, C))
            kind = "hybrid"
        else:
            gcost = float(NW.cuccaro_add(c, src[:8], low, None, debts))
            kind = "cuccaro"
    finally:
        c.alloc = a0
    for w in low:
        c.x(w)
    # T' re-made (w6' = tgt[7] restored): t8 = v3^v4^v5^v6^w6' ^ (v3^w6'^ctl)(v3^v5)
    c.free.remove(t8); c.live += 1
    chk()
    c.cx(v3, t7); c.cx(d, t7); c.cx(v3, v5)
    c.and_compute(t7, v5, t8)
    c.cx(v3, v5); c.cx(d, t7); c.cx(v3, t7)
    _cx_all(c, (v3, v4, v5, v6, t7), t8)
    # top' re-made by F1's form (reads T')
    c.free.remove(top); c.live += 1
    chk()
    W._r24_form(c, tgt, d, v, top, erase=False)
    undo_cost = float(c.toffolis_exp(c.ops))
    c.blocks[bid] = c.ops
    c.ops = save
    if not hasattr(c, "block_log"):
        c.block_log = []
    c.block_log.append((n, C, undo_cost))
    LOG.append(dict(r24n=int(top), cap=cap, C_undo=C, game=kind, game_cost=gcost, fwd=fwd_cost, undo=undo_cost,
                    peak=peak[0], debts=len(debts)))
    cap23 = (W._STATE.get("cfg") or {}).get("_fc16r23n")
    if cap23:                                        # row 23 is the next inner add of MAIN: route it to block23
        c._fc16_r23 = dict(ctl=top, v=list(v), m=7, cap=int(cap23), expect=list(tgt[:9]))
        install_add_hook()
    return debts


def _game(c, NW, src, low, debts, C, chk_alloc):
    for w in low:
        c.x(w)
    a0 = c.alloc
    c.alloc = chk_alloc
    try:
        n_ = len(low)
        if NW.plan(max(0, n_ - 2 - NW.free_prefix(src, n_, None)), C)[0] is not None:
            g, kind = float(NW.game_add(c, src, low, None, debts, C)), "game"
        elif C > NW.cuccaro_need(src, low, None):
            g, kind = float(NW.hybrid_add(c, src, low, None, debts, C)), "hybrid"
        else:
            g, kind = float(NW.cuccaro_add(c, src, low, None, debts)), "cuccaro"
    finally:
        c.alloc = a0
    for w in low:
        c.x(w)
    return g, kind


def block23(c, src, tgt, info, cap):
    """MAIN row 23 at one wire less (R63-INV-ROW9's route): window R24[23..32] framed, control d = R24[33].
    forward  the top input w' (R24[32] framed) X-measured (a debt, repaid when the undo re-makes it); the narrow
             ripple on positions 0..8; the top result o9' (R23[32] framed) by ROW9's 3-AND form with one helper:
             h0 = (v3^v5^v6^o8')(v4^v6^d^1) [helper], h1 = (v3^v5^v6^1)(v3^v5^v6^o8'^d^h0),
             h2 = (v5^o8'^1)(v6^o8'^h0), o9' = v3^v4^v5^o8'^d^1^h0^h1^h2          [framed; 0 differ on all inputs]
    undo     o9' erased by its form (h0 re-made under the outcome), the game on positions 0..8, w' re-made by
             ROW9's 3-AND form g0 = (v2^v3^v6^w7')(v2^v5^w7'^d) [helper], g1 = g0 (v2^v3^v5^w7'),
             g2 = (v4^d^g0)(v4^w7'^d^1), w' = v4^v5^g1^g2^d^1 [0 differ], then Z(w') repays the forward's debt."""
    NW = sys.modules["invsqrt_fc_narrow"]
    n = len(tgt)
    assert n == 10 and list(tgt[1:]) == info["expect"], "r23n: row 23 not where expected"
    src = list(src) + [None] * (n - len(src))
    d, v = info["ctl"], info["v"]
    v2, v3, v4, v5, v6 = v[2], v[3], v[4], v[5], v[6]
    top, t8 = tgt[9], tgt[8]
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
    mt = c.hmr(top)                                  # the debt (-1)^(mt w')
    c.free.insert(0, top)
    c.live -= 1
    c.alloc = alloc_chk
    try:
        debts = NW.fwd(c, src[:9], list(tgt[:9]), None)
    finally:
        c.alloc = a_orig
    c.free.remove(top); c.live += 1
    chk()
    h = c.alloc(1)[0]
    chk()
    # h <- h0 = X0 Y0, X0 = v3^v5^v6^o8' (on v3), Y0 = v4^v6^d^1 (on v4)
    _cx_all(c, (v5, v6, t8), v3)
    c.cx(v6, v4); c.cx(d, v4); c.x(v4)
    c.and_compute(v3, v4, h)
    c.x(v4); c.cx(d, v4); c.cx(v6, v4)
    _cx_all(c, (t8, v6, v5), v3)
    # top <- h1 = (v3^v5^v6^1)(v3^v5^v6^o8'^d^h0): Y1 on h (plain v's), then X1 on v3
    _cx_all(c, (v3, v5, v6, t8, d), h)               # h = Y1
    c.cx(v5, v3); c.cx(v6, v3); c.x(v3)              # v3 = X1
    c.and_compute(v3, h, top)
    c.x(v3); c.cx(v6, v3); c.cx(v5, v3)              # v3 plain
    _cx_all(c, (d, t8, v6, v5, v3), h)               # h = h0
    # top ^= h2 = (v5^o8'^1)(v6^o8'^h0): X2 on v5, Y2 on h
    c.cx(t8, v5); c.x(v5)
    c.cx(v6, h); c.cx(t8, h)
    c.ccx(v5, h, top)
    c.cx(t8, h); c.cx(v6, h)
    c.x(v5); c.cx(t8, v5)
    # top ^= h0 ^ v3 ^ v4 ^ v5 ^ o8' ^ d ^ 1
    c.cx(h, top)
    _cx_all(c, (v3, v4, v5, t8, d), top); c.x(top)
    # helper h0 erased free
    _cx_all(c, (v5, v6, t8), v3)
    c.cx(v6, v4); c.cx(d, v4); c.x(v4)
    c.and_erase(v3, v4, h)
    c.x(v4); c.cx(d, v4); c.cx(v6, v4)
    _cx_all(c, (t8, v6, v5), v3)
    fwd_cost = float(c.toffolis_exp(c.ops[i0:]))
    c.ops.append(("BEND", bid))
    # ---------------- explicit undo
    save = c.ops
    c.ops = []
    # o9' erased: strip the linear part, HMR, CZ for h0, h1, h2 (h0 re-made under the outcome for h1, h2)
    _cx_all(c, (v3, v4, v5, t8, d), top); c.x(top)
    m9 = c.hmr(top)
    c.free_([top])
    h = c.alloc(1)[0]
    chk()
    _cx_all(c, (v5, v6, t8), v3)
    c.cx(v6, v4); c.cx(d, v4); c.x(v4)
    c.cz(v3, v4, m9)                                 # h0's phase
    c.cand_compute(v3, v4, h, m9)                    # h = h0 if m9
    c.x(v4); c.cx(d, v4); c.cx(v6, v4)
    _cx_all(c, (t8, v6, v5), v3)
    _cx_all(c, (v3, v5, v6, t8, d), h)               # h = Y1 (when m9)
    c.cx(v5, v3); c.cx(v6, v3); c.x(v3)              # v3 = X1
    c.cz(v3, h, m9)                                  # h1's phase
    c.x(v3); c.cx(v6, v3); c.cx(v5, v3)
    _cx_all(c, (d, t8, v6, v5, v3), h)               # h = h0
    c.cx(t8, v5); c.x(v5)
    c.cx(v6, h); c.cx(t8, h)
    c.cz(v5, h, m9)                                  # h2's phase
    c.cx(t8, h); c.cx(v6, h)
    c.x(v5); c.cx(t8, v5)
    _cx_all(c, (v5, v6, t8), v3)
    c.cx(v6, v4); c.cx(d, v4); c.x(v4)
    c.cand_erase(v3, v4, h, m9)
    c.x(v4); c.cx(d, v4); c.cx(v6, v4)
    _cx_all(c, (t8, v6, v5), v3)
    c.free_([h])
    C = cap - c.live
    gcost, kind = _game(c, NW, src[:9], list(tgt[:9]), debts, C, alloc_chk)
    # w' re-made on the top wire from (v, d, w7' = tgt[8] restored), then the forward's debt repaid
    c.free.remove(top); c.live += 1
    chk()
    h = c.alloc(1)[0]
    chk()
    _cx_all(c, (v2, v3, v6), t8)                     # X = v2^v3^v6^w7'
    c.cx(v2, v5); c.cx(d, v5)                        # Y = v2^v5^d^w7'... w7' needed in Y too: use t8's copy
    c.cx(t8, v5)                                     # v5 = v5^v2^d^(v2^v3^v6^w7') = v5^d^v3^v6^w7'
    c.cx(v3, v5); c.cx(v6, v5); c.cx(v2, v5)         # v5 = v2^v5^d^w7'  (Y)
    c.and_compute(t8, v5, h)                         # h = g0
    c.cx(v2, v5); c.cx(v6, v5); c.cx(v3, v5); c.cx(t8, v5); c.cx(d, v5); c.cx(v2, v5)
    _cx_all(c, (v6, v3, v2), t8)                     # t8 = w7' again
    # top <- g1 = g0 (v2^v3^v5^w7'): the form on t8
    _cx_all(c, (v2, v3, v5), t8)
    c.and_compute(h, t8, top)
    _cx_all(c, (v5, v3, v2), t8)
    # top ^= g2 = (v4^d^g0)(v4^w7'^d^1): X on h, Y on t8
    c.cx(v4, h); c.cx(d, h)
    c.cx(v4, t8); c.cx(d, t8); c.x(t8)
    c.ccx(h, t8, top)
    c.x(t8); c.cx(d, t8); c.cx(v4, t8)
    c.cx(d, h); c.cx(v4, h)
    _cx_all(c, (v4, v5, d), top); c.x(top)
    # g0's helper erased free
    _cx_all(c, (v2, v3, v6), t8)
    c.cx(v2, v5); c.cx(d, v5); c.cx(t8, v5); c.cx(v3, v5); c.cx(v6, v5); c.cx(v2, v5)
    c.and_erase(t8, v5, h)
    c.cx(v2, v5); c.cx(v6, v5); c.cx(v3, v5); c.cx(t8, v5); c.cx(d, v5); c.cx(v2, v5)
    _cx_all(c, (v6, v3, v2), t8)
    c.z(top, mt)                                     # the forward's debt
    undo_cost = float(c.toffolis_exp(c.ops))
    c.blocks[bid] = c.ops
    c.ops = save
    c.block_log.append((n, C, undo_cost))
    LOG.append(dict(r23n=int(top), cap=cap, C_undo=C, game=kind, game_cost=gcost, fwd=fwd_cost, undo=undo_cost,
                    peak=peak[0], debts=len(debts)))


def arm23(c, tgt, info, cap):
    """after ANY row-24 block (FC2's k4, CH1's P6): route row 23's inner add to block23 (cfg "_fc16r23n" = its cap)."""
    c._fc16_r23 = dict(ctl=tgt[9], v=list(info["v"]), m=7, cap=int(cap), expect=list(tgt[:9]))
    install_add_hook()


def install_add_hook():
    """route the inner add of MAIN's row 23 to block23 (flag c._fc16_r23 set by block() right after row 24)."""
    A = sys.modules["invsqrt_fc_arith"]
    if getattr(A.add, "_fc16r23", False):
        return
    orig = A.add

    def add(c, src, tgt, cin=None, force=None):
        info = getattr(c, "_fc16_r23", None)
        if (info is not None and cin is None and force is None and len(tgt) == 10
                and list(tgt)[1:] == info["expect"]):
            c._fc16_r23 = None
            return block23(c, list(src), list(tgt), info, info["cap"])
        return orig(c, src, tgt, cin, force)
    add._fc16r23 = True
    add._fk_ladder = getattr(orig, "_fk_ladder", False)
    A.add = add


def install():
    """route F1's lazy row-24 block through block() when the cfg carries "_fc16r24n" (FC2's / F1's otherwise)."""
    W = sys.modules["invsqrt_fk_w24ret"]
    if getattr(W._r24_lazy_block, "_fc16r24n", False):
        return
    orig = W._r24_lazy_block

    def route(c, src, tgt, info):
        cap = (W._STATE.get("cfg") or {}).get("_fc16r24n")
        if cap:
            return block(c, src, tgt, info, int(cap))
        return orig(c, src, tgt, info)
    route._fc16r24n = True
    W._r24_lazy_block = route
