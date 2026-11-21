#!/usr/bin/env python3
"""invsqrt_fw_fcpol.py -- F3 (ninth holder, 2026-11-14): R32-game-optimal ADAPTIVE undo for FC's narrow blocks.

Why: FC's narrow undo (tools/invsqrt_fc_narrow.py, read-only here) runs R32's main line with STATIC repairs (each
ghost re-placed from its anchor and cleared again). For the shape (L 8 game carries, C 3 clean wires) that prices
12.5 ANDs; the optimal adaptive policy of the same game (tools/r32_inv_clean_rp_game.py, solve(10, 3, feed 0)) is
12.0: a re-placed low carry stays on its wire for the next ghost's repair. MAIN's undo in family C has 16 such
blocks per shot (FC's block profile of c_W42h1b4): -8 Toffolis per shot.

The policy (game indices 1..8 = carries F+1..F+8 after FC's free prefix F), enumerated on all 16 outcome paths:
  prefix: place 1, 2; remove 1; place 3, 4; remove 3; place 5; GHOST 2 (b2); place 6; remove 5; place 7;
          GHOST 6 (b6); place 8; top; remove 8                                         (9 ANDs)
  if b6: place 5 -> S; repay 6                                                          (0.5)
  GHOST 7 (b7)
  if d67 = b7 AND NOT b6: place 5 -> S                                                  (0.25)
  if b7: place 6 -> T (from S); repay 7; remove 6                                       (0.5)
  if e67 = b6 OR b7: remove 5 (S)
  if b2: place 1 -> U; repay 2                                                          (0.5)
  GHOST 4 (b4)
  if d24 = b4 AND NOT b2: place 1 -> U                                                  (0.25)
  if b4: place 2 -> V (from U); place 3 -> X; repay 4; remove 3; remove 2               (1.0)
  if e24 = b2 OR b4: remove 1 (U)
  expected 12.0 at 3 wires (peak 3). The conditions d / e are classical bits derived from the outcomes (kickmix:
  nested PUSH_CONDITION + BIT_INVERT; the official verifier charges a gate gated by an AND of two coins on the shots
  where it fires [M tools/audit_micro_invsqrt.py, research/ws_scorer_audit.md]).

Kit extension (monkeypatch of invsqrt_fc_kit, in memory): ops ("BAND", d, ((bit, pol), ...)) d = AND of literals
(at most one negative literal) and ("BXOR", d, (bits...)); Circ.band / Circ.bxor; Sim runs them; kmx_text writes
them; toffolis_exp weights a CANDC gated by a derived bit by that bit's probability.

  python3 invsqrt_fw_fcpol.py test [shots]     random-value check of the new undo on FC's block shapes (main kind)
  apply()                                      installs the kit extension and the undo for (L 8, C 3)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_fc_kit as K      # noqa: E402  (FC's own code, read-only)
import invsqrt_fc_narrow as NW  # noqa: E402  (FC's own code, read-only)

_APPLIED = {}


# ---------------------------------------------------------------- kit extension: derived classical bits
def _band(self, lits):
    lits = tuple((int(b), int(p)) for b, p in lits)
    assert sum(1 for _, p in lits if not p) <= 1 and any(p for _, p in lits)
    d = self.bit()
    self.ops.append(("BAND", d, lits))
    if not hasattr(self, "bprob"):
        self.bprob = {}
    pr = 1.0
    for b, _ in lits:
        pr *= self.bprob.get(b, 0.5) if b in getattr(self, "bprob", {}) else 0.5
    self.bprob[d] = pr
    return d


def _bxor(self, bits):
    bits = tuple(int(b) for b in bits)
    d = self.bit()
    self.ops.append(("BXOR", d, bits))
    return d


def _xand_erase(self, a, b, t, g):
    """erase t = a b ONLY where bit g is 1 (HMR and its CZ fix both inside PUSH_CONDITION g): t may hold another
    value where g is 0 (FC's CANDE measures t unconditionally, so it needs t = 0 there)."""
    assert len({a, b, t}) == 3
    m = self.bit()
    self.ops.append(("XANDE", a, b, t, g, m))


def _toffolis_exp(self, ops=None):
    ops = self.ops if ops is None else ops
    bp = getattr(self, "bprob", {})
    return (sum(1 for o in ops if o[0] in ("CCX", "CCZ", "CSWAP", "ANDC"))
            + sum(bp.get(o[4], 0.5) for o in ops if o[0] == "CANDC"))


def _make_run(orig_run):
    def run(self, ops, start=0, stop=None):
        seg = ops[start:stop]
        i0 = 0
        for i, o in enumerate(seg):
            if o[0] == "CANDC":            # same op; the zero-target check holds where the condition fires
                if i > i0:
                    orig_run(self, seg[i0:i])
                _, a, b, t, gb = o
                g = self.bits[gb]
                if self.check and (self.q[t] & g).any():
                    raise AssertionError("CANDC onto a nonzero target q%d where its condition holds" % t)
                self.q[t] ^= self.q[a] & self.q[b] & g
                i0 = i + 1
                continue
            if o[0] == "XANDE":
                if i > i0:
                    orig_run(self, seg[i0:i])
                _, a, b, t, gb, mb = o
                g = self.bits[gb]
                ab = self.q[a] & self.q[b]
                if self.check and ((self.q[t] != ab) & g).any():
                    raise AssertionError("XANDE: target q%d is not q%d q%d where its condition holds" % (t, a, b))
                m = self.rng.integers(0, 2, self.n).astype(bool) & g
                self.bits[mb] = m
                self.sign ^= self.q[t] & m
                self.q[t] = self.q[t] & ~g
                self.sign ^= ab & m
                i0 = i + 1
                continue
            if o[0] in ("BAND", "BXOR"):
                if i > i0:
                    orig_run(self, seg[i0:i])
                if o[0] == "BAND":
                    v = None
                    for b, p in o[2]:
                        x = self.bits[b] if p else ~self.bits[b]
                        v = x.copy() if v is None else (v & x)
                    self.bits[o[1]] = v
                else:
                    v = None
                    for b in o[2]:
                        v = self.bits[b].copy() if v is None else (v ^ self.bits[b])
                    self.bits[o[1]] = v
                i0 = i + 1
        if i0 < len(seg):
            orig_run(self, seg[i0:])
    return run


def _make_kmx(orig_kmx):
    """wrap FC's writer: each BAND / BXOR op becomes a marker the original skips (CBEG/CEND are skipped when
    KMX_CBLOCKS is off), then the marker lines are replaced by the BIT_INVERT text."""
    def kmx_text(circ, x_wires, qmax=None):
        assert not K.KMX_CBLOCKS, "fcpol: derived bits need KMX_CBLOCKS off"
        saved = circ.ops
        texts = []
        ops2 = []
        for o in saved:
            if o[0] == "XANDE":
                _, a, b, t, g, m = o
                texts.append("\n".join([f"PUSH_CONDITION if b{g}", f"HMR q{t} b{m}", f"CZ q{a} q{b} if b{m}",
                                        "POP_CONDITION"]))
                ops2.append(("X", -1000 - (len(texts) - 1), None))
                continue
            if o[0] in ("BAND", "BXOR"):
                d = o[1]
                if o[0] == "BAND":
                    pos = [b for b, p in o[2] if p]
                    neg = [b for b, p in o[2] if not p]
                    t = [f"PUSH_CONDITION if b{b}" for b in pos] + [f"BIT_INVERT b{d}"]
                    t += [f"BIT_INVERT b{d} if b{b}" for b in neg] + ["POP_CONDITION"] * len(pos)
                else:
                    t = [f"BIT_INVERT b{d} if b{b}" for b in o[2]]
                texts.append("\n".join(t))
                ops2.append(("X", -1000 - (len(texts) - 1), None))     # placeholder line "X q-1000-k"
            else:
                ops2.append(o)
        circ.ops = ops2
        try:
            txt = orig_kmx(circ, x_wires, qmax)
        finally:
            circ.ops = saved
        for k, t in enumerate(texts):
            txt = txt.replace(f"X q{-1000 - k}\n", t + "\n", 1)
        return txt
    return kmx_text


def apply_kit():
    if _APPLIED.get("kit"):
        return
    K.Circ.band = _band
    K.Circ.bxor = _bxor
    K.Circ.xand_erase = _xand_erase
    K.Circ.toffolis_exp = _toffolis_exp
    K.Sim.run = _make_run(K.Sim.run)
    K.kmx_text = _make_kmx(K.kmx_text)
    _APPLIED["kit"] = True


# ---------------------------------------------------------------- the adaptive undo for (L 8, C 3)
def pol83_add(c, src, tgt, cin, debts, C):
    """tgt += src + cin re-creating the carries by the optimal adaptive policy at (L 8, C 3); repays the debts.
    Same entry and conventions as FC's game_add (s-form carries, free prefix F, first placement writes the sum and
    repays the forward debt)."""
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    carry = {0: NW._lit(cin)}

    def write_sum(r):
        NW._xor_src(c, src[r], tgt[r])
        NW._xor_lit(c, carry[r], tgt[r])

    write_sum(0)
    F = 0
    while F + 1 <= n - 2:
        f = NW._form_s(src[F], tgt[F], carry[F])
        if f[0] == "and":
            break
        carry[F + 1] = ("k", f[1]) if f[0] == "const" else ("w", f[1], f[2])
        F += 1
        write_sum(F)
    L = n - 2 - F
    assert L == 8 and C >= 3, (L, C)

    def place(r, first, cond=None, into=None):
        R = F + r
        f = NW._form_s(src[R - 1], tgt[R - 1], carry[R - 1])
        assert f[0] == "and"
        _, xw, yw, ops = f
        t = c.alloc(1)[0] if into is None else into
        NW._emit(c, ops)
        if cond is None:
            c.and_compute(xw, yw, t)
        else:
            c.cand_compute(xw, yw, t, cond)
        NW._emit(c, ops[::-1])
        NW._xor_lit(c, carry[R - 1], t, cond)
        carry[R] = ("w", t, 0)
        if first:
            write_sum(R)
            if R in debts:
                c.z(t, debts[R])
        return t

    def remove(r, cond=None):
        R = F + r
        t = carry[R][1]
        NW._xor_lit(c, carry[R - 1], t, cond)
        f = NW._form_s(src[R - 1], tgt[R - 1], carry[R - 1])
        _, xw, yw, ops = f
        NW._emit(c, ops)
        if cond is None:
            c.and_erase(xw, yw, t)
        else:
            c.cand_erase(xw, yw, t, cond)
            c.free_([t])
        NW._emit(c, ops[::-1])
        del carry[R]

    def ghost(r):
        R = F + r
        t = carry[R][1]
        g = c.hmr(t)
        c.free_([t])
        del carry[R]
        return g

    def repay(r, g):
        R = F + r
        lo = carry[R - 1]
        if lo[0] == "k":
            if lo[1]:
                c.neg(g)
        else:
            c.z(lo[1], g)
            if lo[2]:
                c.neg(g)
        f = NW._form_s(src[R - 1], tgt[R - 1], carry[R - 1])
        _, xw, yw, ops = f
        NW._emit(c, ops)
        c.cz(xw, yw, g)
        NW._emit(c, ops[::-1])

    def top():
        tw = tgt[n - 1]
        f = NW._form_s(src[n - 2], tgt[n - 2], carry[n - 2])
        if f[0] == "const":
            if f[1]:
                c.x(tw)
        elif f[0] == "alias":
            c.cx(f[1], tw)
            if f[2]:
                c.x(tw)
        else:
            _, xw, yw, ops = f
            NW._emit(c, ops)
            c.ccx(xw, yw, tw)
            NW._emit(c, ops[::-1])
            NW._xor_lit(c, carry[n - 2], tw)
        NW._xor_src(c, src[n - 1], tw)

    # prefix (main line, unconditional): 8 first placements + top = 9
    place(1, True); place(2, True); remove(1); place(3, True); place(4, True); remove(3); place(5, True)
    b2 = ghost(2)
    place(6, True); remove(5); place(7, True)
    b6 = ghost(6)
    place(8, True); top(); remove(8)
    # ghosts 6 and 7
    S = place(5, False, b6)
    repay(6, b6)
    b7 = ghost(7)
    d67 = c.band([(b7, 1), (b6, 0)])
    place(5, False, d67, into=S)
    place(6, False, b7)
    repay(7, b7)
    remove(6, b7)
    e67 = c.bxor([b6, d67])
    remove(5, e67)
    # ghosts 2 and 4
    U = place(1, False, b2)
    repay(2, b2)
    b4 = ghost(4)
    d24 = c.band([(b4, 1), (b2, 0)])
    place(1, False, d24, into=U)
    place(2, False, b4)
    place(3, False, b4)
    repay(4, b4)
    remove(3, b4)
    remove(2, b4)
    e24 = c.bxor([b2, d24])
    remove(1, e24)
    assert set(carry) <= set(range(0, F + 1)), sorted(carry)
    return 12.0


def _ctx(c, src, tgt, cin, debts):
    """FC's game_add set-up (s-form carries, free prefix) and its move primitives, with `place` able to target a given
    zero wire and `write_sum` only on first placements (as in FC)."""
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    carry = {0: NW._lit(cin)}

    def write_sum(r):
        NW._xor_src(c, src[r], tgt[r])
        NW._xor_lit(c, carry[r], tgt[r])

    write_sum(0)
    F = 0
    while F + 1 <= n - 2:
        f = NW._form_s(src[F], tgt[F], carry[F])
        if f[0] == "and":
            break
        carry[F + 1] = ("k", f[1]) if f[0] == "const" else ("w", f[1], f[2])
        F += 1
        write_sum(F)

    def place(r, first, cond=None, into=None):
        R = F + r
        f = NW._form_s(src[R - 1], tgt[R - 1], carry[R - 1])
        assert f[0] == "and"
        _, xw, yw, ops = f
        t = c.alloc(1)[0] if into is None else into
        NW._emit(c, ops)
        if cond is None:
            c.and_compute(xw, yw, t)
        else:
            c.cand_compute(xw, yw, t, cond)
        NW._emit(c, ops[::-1])
        NW._xor_lit(c, carry[R - 1], t, cond)
        carry[R] = ("w", t, 0)
        if first:
            write_sum(R)
            if R in debts:
                c.z(t, debts[R])
        return t

    def remove(r, cond=None, exclusive=False):
        R = F + r
        t = carry[R][1]
        NW._xor_lit(c, carry[R - 1], t, cond)
        f = NW._form_s(src[R - 1], tgt[R - 1], carry[R - 1])
        _, xw, yw, ops = f
        NW._emit(c, ops)
        if cond is None:
            c.and_erase(xw, yw, t)
        elif exclusive:
            c.xand_erase(xw, yw, t, cond)
            c.free_([t])
        else:
            c.cand_erase(xw, yw, t, cond)
            c.free_([t])
        NW._emit(c, ops[::-1])
        del carry[R]

    def ghost(r):
        R = F + r
        t = carry[R][1]
        g = c.hmr(t)
        c.free_([t])
        del carry[R]
        return g

    def repay(r, g):
        R = F + r
        lo = carry[R - 1]
        if lo[0] == "k":
            if lo[1]:
                c.neg(g)
        else:
            c.z(lo[1], g)
            if lo[2]:
                c.neg(g)
        f = NW._form_s(src[R - 1], tgt[R - 1], carry[R - 1])
        _, xw, yw, ops = f
        NW._emit(c, ops)
        c.cz(xw, yw, g)
        NW._emit(c, ops[::-1])

    def top():
        tw = tgt[n - 1]
        f = NW._form_s(src[n - 2], tgt[n - 2], carry[n - 2])
        if f[0] == "const":
            if f[1]:
                c.x(tw)
        elif f[0] == "alias":
            c.cx(f[1], tw)
            if f[2]:
                c.x(tw)
        else:
            _, xw, yw, ops = f
            NW._emit(c, ops)
            c.ccx(xw, yw, tw)
            NW._emit(c, ops[::-1])
            NW._xor_lit(c, carry[n - 2], tw)
        NW._xor_src(c, src[n - 1], tw)

    return n - 2 - F, F, carry, place, remove, ghost, repay, top


def pol93_add(c, src, tgt, cin, debts, C):
    """(L 9, C 3): the R32 game's optimal adaptive policy (solve(11, 3, feed 0) = 14.0; FC's plan 15.0).
      prefix: place 1, 2, 3; remove 2, 1; place 4, 5; remove 4; place 6; GHOST 3 (b3); place 7; remove 6; place 8;
              GHOST 7 (b7); place 9; top; remove 9                                                (10)
      b7: place 6 -> S6; repay 7 | GHOST 8 (b8) | d78 = b8 & ~b7: place 6 -> S6 | b8: place 7 (from S6); repay 8;
      remove 7 | e78 = b7 | b8: remove 6                                                          (1.25)
      b3: place 1 -> S1; place 2 -> S2; repay 3 | GHOST 5 (b5) | f35 = b3 & ~b5: remove 2, remove 1 |
      g35 = b3 & b5: remove 1 | d35 = b5 & ~b3: place 1 -> S1, place 2 -> S2, remove 1 | b5: place 3 (from S2),
      place 4, repay 5, remove 4, remove 3 | GHOST 2 (b2', S2 holds c2 where b5, 0 elsewhere) | h = b5 & b2':
      place 1, repay 2, remove 1                                                                  (2.75)"""
    L, F, carry, place, remove, ghost, repay, top = _ctx(c, src, tgt, cin, debts)
    assert L == 9 and C >= 3, (L, C)
    place(1, True); place(2, True); place(3, True); remove(2); remove(1); place(4, True); place(5, True); remove(4)
    place(6, True)
    b3 = ghost(3)
    place(7, True); remove(6); place(8, True)
    b7 = ghost(7)
    place(9, True); top(); remove(9)
    S6 = place(6, False, b7)
    repay(7, b7)
    b8 = ghost(8)
    d78 = c.band([(b8, 1), (b7, 0)])
    place(6, False, d78, into=S6)
    place(7, False, b8)
    repay(8, b8)
    remove(7, b8)
    e78 = c.bxor([b7, d78])
    remove(6, e78)
    S1 = place(1, False, b3)
    S2 = place(2, False, b3)
    repay(3, b3)
    b5 = ghost(5)
    f35 = c.band([(b3, 1), (b5, 0)])
    g35 = c.band([(b3, 1), (b5, 1)])
    d35 = c.band([(b5, 1), (b3, 0)])
    # b3 & ~b5: clear c2, c1 (FC's remove frees the wire; re-claim the slots for the other branches below)
    remove(2, f35, exclusive=True)
    remove(1, f35, exclusive=True)
    # re-establish the slot bookkeeping: S1 / S2 still hold c1 / c2 in the b3 & b5 shots
    _reclaim(c, [S1, S2])
    carry[F + 1] = ("w", S1, 0)
    carry[F + 2] = ("w", S2, 0)
    remove(1, g35)
    _reclaim(c, [S1])
    carry[F + 1] = ("w", S1, 0)
    place(1, False, d35, into=S1)
    place(2, False, d35, into=S2)
    remove(1, d35)
    place(3, False, b5)
    place(4, False, b5)
    repay(5, b5)
    remove(4, b5)
    remove(3, b5)
    b2r = ghost(2)
    h = c.band([(b5, 1), (b2r, 1)])
    place(1, False, h)
    repay(2, h)
    remove(1, h)
    assert set(carry) <= set(range(0, F + 1)), sorted(carry)
    return 14.0


def _reclaim(c, wires):
    """take back wires that a conditioned remove freed although other branches still hold values on them: remove
    them from the free list and count them live again (they stay allocated until their last conditioned remove)."""
    for w in wires:
        assert w in c.free, w
        c.free.remove(w)
        c.live += 1
        c.peak = max(c.peak, c.live)


def apply():
    """install the kit extension and route FC's (L 8, C 3) and (L 9, C 3) undo blocks to the adaptive policies."""
    apply_kit()
    if _APPLIED.get("undo"):
        return
    orig = NW.game_add

    def game_add(c, src, tgt, cin, debts, C):
        n = len(tgt)
        L = n - 2 - NW.free_prefix(src, n, cin) if n >= 3 else -1
        if L == 8 and C == 3:
            return pol83_add(c, src, tgt, cin, debts, C)
        if L == 9 and C == 3:
            return pol93_add(c, src, tgt, cin, debts, C)
        return orig(c, src, tgt, cin, debts, C)
    NW.game_add = game_add
    _APPLIED["undo"] = True


def _test(nsh=4000):
    """FC's narrow test (main / wires / tail / constcin kinds) restricted to the shapes the new undo takes."""
    import numpy as np
    apply()
    rng = np.random.default_rng(11)
    nbad = ntest = 0
    for n in range(9, 15):
        for kind in ("main", "wires", "tail", "constcin"):
            c = K.Circ(n_fixed=0)
            tgt = c.alloc(n)
            if kind == "wires":
                srcw = c.alloc(n)
                src = list(srcw)
                cin = c.alloc(1)[0]
            elif kind == "main":
                v = c.alloc(max(0, n - 3))
                src = ([NW.ONE] + list(v) + [NW.ONE] + [None])[:n]
                cin = None
            elif kind == "tail":
                v = c.alloc(max(0, n - 1))
                src = (list(v) + [NW.ONE])[:n]
                cin = c.alloc(1)[0]
            else:
                v = c.alloc(max(0, n - 2))
                src = ([NW.ONE, None] + list(v))[:n]
                cin = NW.ONE
            if n - 2 - NW.free_prefix(src, n, cin) not in (8, 9):
                continue
            live0 = c.live
            C = 3
            c.cap = live0 + C
            c.blocks = {}
            fwd_start = len(c.ops)
            NW.block(c, src, tgt, cin)
            fwd_ops = c.ops[fwd_start:]
            undo = K.Circ.mirror_blocks(fwd_ops, c.blocks)
            sim = K.Sim(nsh, c.next + 4, seed=n * 7 + 1)
            vals = {}
            for w in range(live0):
                vals[w] = rng.integers(0, 2, nsh).astype(bool)
                sim.q[w] = vals[w]
            T = sim.get_word(tgt)
            S = np.zeros(nsh, dtype=np.int64)
            for k, a in enumerate(src):
                if NW._w(a):
                    S |= sim.q[a].astype(np.int64) << k
                elif a == NW.ONE:
                    S |= np.int64(1) << k
            ci = sim.q[cin].astype(np.int64) if NW._w(cin) else (1 if cin == NW.ONE else 0)
            sim.run(fwd_ops)
            ok1 = np.array_equal(sim.get_word(tgt), (T + S + ci) % (1 << n))
            sim.run(undo)
            ok2 = np.array_equal(sim.get_word(tgt), T) and all(np.array_equal(sim.q[w], vals[w]) for w in range(live0))
            others = list(range(live0, c.next))
            ok3 = not sim.q[others].any() if others else True
            ok4 = not sim.sign.any()
            exp = c.toffolis_exp(undo)
            ntest += 1
            good = ok1 and ok2 and ok3 and ok4 and c.peak <= c.cap
            nbad += 0 if good else 1
            print(f"n {n} {kind:8s} C {C}: fwd {ok1} undo {ok2} anc0 {ok3} phase {ok4} peak {c.peak}/{c.cap} "
                  f"undo expected ANDs {exp:.3f} (log {c.block_log[-1]})")
    print(f"fcpol test: {ntest} shapes, bad {nbad}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        _test(int(sys.argv[2]) if len(sys.argv) > 2 else 4000)
