#!/usr/bin/env python3
"""invsqrt_fw_m7mux.py -- F3 (twelfth holder, 2026-11-15): the explicit undo of FC2's mux add (h1c, s1 form) in the
cheapest form found by tools/invsqrt_fw_m7muxopt.py: every mux position p = 2..10 is
  rot   rotated in place (a run p1..p2 is one chain over T[p1-3..p2-2], CSWAP under f, undone after the game),
  perm  a temporary made before the game (1 AND), erased after it (free),
  trans made by the game at its first read and erased after its last read (1 AND; a wire only in that window).
The game is a copy of FC's narrow game_add (tools/invsqrt_fc_narrow.py, read-only) with lazy sources.
Default spec (24.0 against 26.0 today and 25.0 for muxsplit [6, 8] [M prof]):
  {"rot": [3, 4, 5, 6], "perm": [2, 7, 8], "trans": [9, 10], "C": 4}
"""

DEFAULT = {"rot": [3, 4, 5, 6], "perm": [2, 7, 8], "trans": [9, 10], "C": 4}
# "knowntop": 1 -- the mux add is the first add into c0 (0 before it on every input [M]): the undo game's top step
# X-measures c0 and repays the phase by Z / CZ / NEG (no CCX), and position 10's mux bit is never made (its phase is
# Clifford on T7, T8, f): -2 per shot


def game_add_lazy(c, NW, src, tgt, cin, debts, C, lazy, knowntop=None, top_phase=None):
    """FC's game_add (copy) where src[k] for k in `lazy` is made by lazy[k][0]() at its first read and erased by
    lazy[k][1](wire) after its last read (read counts from a dry run of the same plan)."""
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    total = {k: 0 for k in lazy}
    count = {k: 0 for k in lazy}
    held = {}
    dry = [True]

    def S(k):
        if k not in lazy:
            return src[k]
        if dry[0]:
            return -1000 - k                           # a dummy wire index (never emitted in the dry run)
        if k not in held:
            held[k] = lazy[k][0]()
        return held[k]

    def done(k):
        if k not in lazy:
            return
        if dry[0]:
            total[k] += 1
            return
        count[k] += 1
        if count[k] == total[k]:
            lazy[k][1](held.pop(k))

    def run(c_):
        carry = {0: NW._lit(cin)}
        fresh = set()

        def write_sum(r):
            a = S(r)
            NW._xor_src(c_, a, tgt[r])
            done(r)
            NW._xor_lit(c_, carry[r], tgt[r])

        write_sum(0)
        F = 0
        while F + 1 <= n - 2:
            f = NW._form_s(S(F), tgt[F], carry[F])
            done(F)
            if f[0] == "and":
                break
            carry[F + 1] = ("k", f[1]) if f[0] == "const" else ("w", f[1], f[2])
            F += 1
            write_sum(F)
        L = n - 2 - F

        def place(r, first, cond=None):
            f = NW._form_s(S(r - 1), tgt[r - 1], carry[r - 1])
            assert f[0] == "and"
            _, xw, yw, ops = f
            t = c_.alloc(1)[0]
            NW._emit(c_, ops)
            if cond is None:
                c_.and_compute(xw, yw, t)
            else:
                c_.cand_compute(xw, yw, t, cond)
            NW._emit(c_, ops[::-1])
            done(r - 1)
            NW._xor_lit(c_, carry[r - 1], t, cond)
            carry[r] = ("w", t, 0)
            fresh.add(r)
            if first:
                write_sum(r)
                if r in debts:
                    c_.z(t, debts[r])

        def remove(r, cond=None):
            t = carry[r][1]
            NW._xor_lit(c_, carry[r - 1], t, cond)
            f = NW._form_s(S(r - 1), tgt[r - 1], carry[r - 1])
            _, xw, yw, ops = f
            NW._emit(c_, ops)
            if cond is None:
                c_.and_erase(xw, yw, t)
            else:
                c_.cand_erase(xw, yw, t, cond)
                c_.free_([t])
            NW._emit(c_, ops[::-1])
            done(r - 1)
            del carry[r]
            fresh.discard(r)

        def ghost(r):
            t = carry[r][1]
            g = c_.hmr(t)
            c_.free_([t])
            del carry[r]
            fresh.discard(r)
            return g

        def repay(r, g):
            lo = carry[r - 1]
            if lo[0] == "k":
                if lo[1]:
                    c_.neg(g)
            else:
                c_.z(lo[1], g)
                if lo[2]:
                    c_.neg(g)
            f = NW._form_s(S(r - 1), tgt[r - 1], carry[r - 1])
            _, xw, yw, ops = f
            NW._emit(c_, ops)
            c_.cz(xw, yw, g)
            NW._emit(c_, ops[::-1])
            done(r - 1)

        def top():
            tw = tgt[n - 1]
            if knowntop == "c0m":                     # FC8 c0m: normal top into the phantom, the
                import invsqrt_fc_s8m3 as _M3t            # top literal as its Clifford phase
                f = NW._form_s(S(n - 2), tgt[n - 2], carry[n - 2])
                if f[0] == "const":
                    if f[1]:
                        c_.x(tw)
                elif f[0] == "alias":
                    c_.cx(f[1], tw)
                    if f[2]:
                        c_.x(tw)
                else:
                    _, xw, yw, ops = f
                    NW._emit(c_, ops)
                    c_.ccx(xw, yw, tw)
                    NW._emit(c_, ops[::-1])
                    NW._xor_lit(c_, carry[n - 2], tw)
                done(n - 2)
                top_phase(c_, _M3t.PH["m0"])
                done(n - 1)
                return
            if knowntop is not None:                   # the top bit's final value is the constant knowntop: X-measure
                g = c_.hmr(tw)                         # it (reset to 0), set it, repay (-1)^(g v), v = its old value
                if knowntop:                           # = NOT knowntop xor a_top xor c_top (Clifford: Z / CZ / NEG)
                    c_.x(tw)
                if not knowntop:
                    pass
                else:
                    c_.neg(g)
                top_phase(c_, g)                       # (-1)^(g a_top) on the wires a_top is made from
                f = NW._form_s(S(n - 2), tgt[n - 2], carry[n - 2])
                if f[0] == "and":                      # (-1)^(g c_top), c_top = c_(n-2) xor x y (s-form)
                    lo = carry[n - 2]
                    if lo[0] == "k":
                        if lo[1]:
                            c_.neg(g)
                    else:
                        c_.z(lo[1], g)
                        if lo[2]:
                            c_.neg(g)
                    _, xw, yw, ops = f
                    NW._emit(c_, ops)
                    c_.cz(xw, yw, g)
                    NW._emit(c_, ops[::-1])
                elif f[0] == "alias":
                    c_.z(f[1], g)
                    if f[2]:
                        c_.neg(g)
                elif f[1]:
                    c_.neg(g)
                done(n - 2)
                return
            f = NW._form_s(S(n - 2), tgt[n - 2], carry[n - 2])
            if f[0] == "const":
                if f[1]:
                    c_.x(tw)
            elif f[0] == "alias":
                c_.cx(f[1], tw)
                if f[2]:
                    c_.x(tw)
            else:
                _, xw, yw, ops = f
                NW._emit(c_, ops)
                c_.ccx(xw, yw, tw)
                NW._emit(c_, ops[::-1])
                NW._xor_lit(c_, carry[n - 2], tw)
            done(n - 2)
            NW._xor_src(c_, S(n - 1), tw)
            done(n - 1)

        moves, cost = NW.plan(L, C)
        if moves is None:
            raise ValueError("narrow undo: no schedule for L=%d at C=%d" % (L, C))
        cb = {"place": lambda k, first, cond: place(F + k, first, cond),
              "remove": lambda k, cond: remove(F + k, cond),
              "ghost": lambda k: ghost(F + k),
              "repay": lambda k, b: repay(F + k, b),
              "top": top,
              "mark": lambda m: c_.ops.append(m)}
        got = NW.run_moves(L, C, list(moves), cb)
        assert got is not None and abs(got - cost) < 1e-9
        assert not fresh
        return cost

    # dry run on a scratch copy of the allocator state (ops discarded) to count the reads of the lazy sources
    import copy
    c_dry = copy.copy(c)
    c_dry.ops = []
    c_dry.free = list(c.free)
    c_dry.nbits = c.nbits
    run(c_dry)
    dry[0] = False
    cost = run(c)
    assert not held, held
    return cost


def mux_undo(c, NW, T, f, R, w0, s1, debts, spec=None, cin=None):   # I39-SLOT10 (inv_i39_patch_slot10.py)
    """emit the explicit undo of the mux add: R -= mux word, debts repaid. Returns the game cost."""
    sp = dict(DEFAULT)
    sp.update(spec or {})
    rot, perm, trans = sorted(sp["rot"]), sorted(sp["perm"]), sorted(sp["trans"])
    lo2 = bool(sp.get("lo2")) and bool(s1)             # R48-INV-S7C: position 2 has no mux bit
    if lo2:
        perm = [p for p in perm if p != 2]
        assert sorted(rot + perm + trans) == list(range(3, 11))
    else:
        assert sorted(rot + perm + trans) == list(range(2, 11)) and 2 not in rot
    runs = []
    for p in rot:
        if runs and runs[-1][1] == p - 1:
            runs[-1][1] = p
        else:
            runs.append([p, p])
    chain_wires = set()
    for p1, p2 in runs:
        chain_wires |= set(range(p1 - 3, p2 - 1))
    for p in trans:
        assert p >= 3 and (p - 3) not in chain_wires and (p - 2) not in chain_wires, ("trans needs unrotated wires", p)

    def mk(p):
        """compute the mux bit of position p into a fresh wire (from unrotated T wires)."""
        z = c.alloc(1)[0]
        if p == 2:
            c.x(f)
            c.and_compute(f, T[0], z)
            c.x(f)
        else:
            a, b = T[p - 3], T[p - 2]
            c.cx(a, b)
            c.and_compute(f, b, z)
            c.cx(a, b)
            c.cx(b, z)
        return z

    def er(p, z):
        if p == 2:
            c.x(f)
            c.and_erase(f, T[0], z, free=False)
            c.x(f)
        else:
            a, b = T[p - 3], T[p - 2]
            c.cx(b, z)
            c.cx(a, b)
            c.and_erase(f, b, z, free=False)
            c.cx(a, b)
        c.free_([z])

    Zp = {p: mk(p) for p in perm}
    kap = None
    if lo2:                                            # X frame first; kappa = f ? R1 : Tl_0 = the carry out of
        for w in R:                                    # position 1 (sources Tl_0 and Tl_0 xor f), T unrotated
            c.x(w)
        kap = c.alloc(1)[0]
        c.cx(T[0], R[0])
        c.and_compute(f, R[0], kap)
        c.cx(T[0], R[0])
        c.cx(T[0], kap)
        if 1 in debts:
            c.z(kap, debts[1])
        c.cx(f, R[0])                                  # the sum bit of position 1: Tl_0 xor (Tl_0 xor f)
    chain = []
    for p1, p2 in runs:
        chain += [(T[k], T[k + 1]) for k in range(p2 - 3, p1 - 4, -1)]
    for a_, b_ in chain:
        c.cswap(f, a_, b_)
    srcu = [w0, f]
    lazy = {}
    for p in range(2, 11):
        if p in perm:
            srcu.append(Zp[p])
        elif p in rot:
            srcu.append(T[p - 2])
        else:
            srcu.append(None)                          # placeholder; made lazily
    if s1:
        srcu = srcu[1:]
    off = 1 if s1 else 0
    kt = bool(sp.get("knowntop"))
    ktc = sp.get("knowntop") == "c0m"                # FC8 c0m
    for p in trans:
        if kt and p == 10:
            continue                                   # position 10 is only a phase under knowntop
        lazy[p - off] = ((lambda p=p: mk(p)), (lambda z, p=p: er(p, z)))

    def top_phase(cc, g):
        """(-1)^(g m10), m10 = f ? Tl_7 : Tl_8 = T8 xor f (T7 xor T8) on unrotated T7, T8 (T9 = f)."""
        cc.z(T[8], g)
        cc.cx(T[7], T[8])
        cc.cz(f, T[8], g)
        cc.cx(T[7], T[8])
    if kt:
        assert 10 in trans and 8 not in chain_wires and 7 not in chain_wires
        srcu[10 - off] = None
    if lo2:
        srcg = list(srcu[1:])
        srcg[0] = None                                 # position 2: no literal
        cost = game_add_lazy(c, NW, srcg, R[1:], kap, {i - 1: g for i, g in debts.items() if i >= 2},
                             int(sp["C"]), {k - 1: v for k, v in lazy.items()},
                             knowntop=("c0m" if ktc else (1 if kt else None)),
                             top_phase=(top_phase if kt else None))
        for a_, b_ in reversed(chain):
            c.cswap(f, a_, b_)
        for p in reversed(perm):
            er(p, Zp[p])
        c.cx(T[0], kap)                                # kappa xor Tl_0 = f (R1_before xor Tl_0), R1_before = R1 xor f
        c.cx(T[0], R[0])
        c.cx(f, R[0])
        c.and_erase(f, R[0], kap, free=False)
        c.cx(f, R[0])
        c.cx(T[0], R[0])
        c.free_([kap])
        for w in R:
            c.x(w)
        return cost + 1.0
    for w in R:
        c.x(w)
    assert cin is None or not lo2, "a carry-in with lo2"   # I39-SLOT10 (inv_i39_patch_slot10.py)
    cost = game_add_lazy(c, NW, srcu, R, cin, debts, int(sp["C"]), lazy,
                         knowntop=("c0m" if ktc else (1 if kt else None)),   # FC8 c0m knowntop
                         top_phase=(top_phase if kt else None))
    for w in R:
        c.x(w)
    for a_, b_ in reversed(chain):
        c.cswap(f, a_, b_)
    for p in reversed(perm):
        er(p, Zp[p])
    return cost
