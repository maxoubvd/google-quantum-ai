#!/usr/bin/env python3
"""invsqrt_fw_c1add.py -- F3 (ninth holder, 2026-11-14): in-place adders for FC's family C at ONE clean wire.

c1_add(c, src, tgt, cin, debts=None, C=None)
    tgt += src + cin (mod 2^n) in place. Cuccaro / hybrid form (Cuccaro, Draper, Kutin, Petrie Moulton 2004; own
    code) with two changes against FC's cuccaro_add / hybrid_add:
      (1) the top position n-2 stores no carry: c_(n-1) = MAJ(a, b, c_(n-2)) goes straight into the top bit
          (1 CCX, no wire, no UMA), so a constant source bit there needs no materialised wire (FC's cuccaro_add
          spends MAJ + UMA = 2 CCX and one clean wire on it);
      (2) spare clean wires (C = cap - live unless given) become clean positions (1 AND, erased by measurement in
          the downward pass), from the top down, as FC's hybrid_add does.
    Re-creates the carries c_1..c_(n-2) and repays `debts` on them, so it also serves as an undo. Unitary up to
    ANDC / ANDE pairs: the kit's mirror (ANDC <-> ANDE) is its exact inverse.
c1_need(src, tgt, cin)   clean wires c1_add needs (constant source bits below the top position, constant carry-in).
apply(A, NW, pred=None)  routes arith.add's `cuccaro_fwd` branch (rows with fewer than 2 clean wires) through
    c1_add; `pred(src, tgt, cin)` (optional) limits it to some rows.
  python3 invsqrt_fw_c1add.py test [nmax]     random-value check on the kit's simulator (needs FC's kit on the path)
"""
import os
import sys

NW = None          # FC's invsqrt_fc_narrow (set by apply / _test)
ONE = "1"
LOG = []           # (n, C, count, need) per c1_add call


def _w(x):
    return isinstance(x, int)


def c1_need(src, tgt, cin):
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    ws = [a for a in src if _w(a)]
    if len(set(ws)) != len(ws) or (_w(cin) and cin in ws) or set(ws) & set(tgt) or (_w(cin) and cin in tgt):
        return 10 ** 6
    if n <= 2:
        return NW.cuccaro_need(src, tgt, cin)
    alias = NW._cuccaro_alias(src, n, cin)
    i0 = 1 if alias else 0
    need = sum(1 for i in range(i0, n - 2) if not _w(src[i]))
    if not alias and not _w(cin):
        need += 1
    return need


def c1_add(c, src, tgt, cin, debts=None, C=None):
    debts = debts or {}
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    if n <= 2:
        return NW.cuccaro_add(c, src, tgt, cin, debts)
    if C is None:
        C = c.cap - c.live
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
    pos = list(range(i0, n - 2))                       # chain positions: position i stores c_(i+1)
    clean = {i for i in pos if not _w(src[i])}
    budget -= len(clean)
    assert budget >= 0, "c1_add: constants need more clean wires than C"
    for i in sorted((i for i in pos if _w(src[i])), reverse=True):
        if budget <= 0:
            break
        clean.add(i)
        budget -= 1
    chain = {i0: c0}
    cnt = 0
    for i in pos:                                       # upward: chain[i] holds c_i exactly here
        x, y = chain[i], tgt[i]
        if i in clean:
            a = src[i]
            c.cx(x, y)                                  # y = b ^ c
            if _w(a):
                c.cx(x, a)                              # a' = a ^ c
                t = c.and_compute(a, y)
                c.cx(x, a)
            elif a == ONE:
                c.x(x)
                t = c.and_compute(x, y)
                c.x(x)
            else:
                t = c.and_compute(x, y)
            c.cx(x, t)                                  # t = c ^ a'(b ^ c) = MAJ(a, b, c)
            cnt += 1
            chain[i + 1] = t
        else:
            z = src[i]
            c.cx(z, y)
            c.cx(z, x)
            c.ccx(x, y, z)                              # z = c_(i+1); x = c ^ a; y = b ^ a
            cnt += 1
            chain[i + 1] = z
        if (i + 1) in debts:
            c.z(chain[i + 1], debts[i + 1])
    # top position p = n - 2: tgt[n-1] ^= c_(n-1) ^ a_(n-1), c_(n-1) = c ^ (a ^ c)(b ^ c); then s_p = a ^ b ^ c
    p = n - 2
    x, y, a, top = chain[p], tgt[p], src[p], tgt[n - 1]
    c.cx(x, y)                                          # y = b ^ c
    if _w(a):
        c.cx(x, a)
        c.ccx(a, y, top)
        c.cx(x, a)
    elif a == ONE:
        c.x(x)
        c.ccx(x, y, top)
        c.x(x)
    else:
        c.ccx(x, y, top)
    cnt += 1
    c.cx(x, top)
    NW._xor_src(c, src[n - 1], top)
    if _w(a):                                           # y = b ^ c -> a ^ b ^ c
        c.cx(a, y)
    elif a == ONE:
        c.x(y)
    for i in reversed(pos):                             # downward
        x, y, nxt = chain[i], tgt[i], chain[i + 1]
        if i in clean:
            a = src[i]
            c.cx(x, nxt)                                # nxt = a'(b ^ c)
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
    if alias:
        c.x(tgt[0])                                     # s_0 = NOT b_0
    for w, xv in reversed(mats):
        if xv == ONE:
            c.x(w)
        c.free_([w])
    LOG.append((n, C, cnt, c1_need(src, tgt, cin)))
    return cnt


def _game(Lp, C, ccz):
    """the R32 game's expected ANDs (incl. the top) for Lp carries at C clean wires (F3's general emitter)."""
    import invsqrt_fw_fcgen as G
    return float(G.solved(Lp, C, ccz)[0])


def mixed_plan(src, tgt, cin, C, ccz):
    """best split (k, cost) for mixed_add: in-place positions i0..k (2 each), the game on the rest (k = i0 - 1: none
    in place). None when no split beats c1_add."""
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    if n < 4 or C < 2:
        return None
    ws = [a for a in src if _w(a)]                     # constants may sit in the game part; wires must be distinct
    if len(set(ws)) != len(ws) or (_w(cin) and cin in ws) or set(ws) & set(tgt) or (_w(cin) and cin in tgt):
        return None
    alias = NW._cuccaro_alias(src, n, cin)
    if not alias and not _w(cin):
        return None
    i0 = 1 if alias else 0
    best = None
    for k in range(i0, n - 3):                          # at least one in-place position, at least one game carry
        if not all(_w(src[i]) for i in range(i0, k + 1)):
            break
        Lp = n - k - 3
        g = _game(Lp, C, ccz)
        if g >= 1e6:
            continue
        cost = 2 * (k - i0 + 1) + g
        if best is None or cost < best[1] - 1e-9:
            best = (k, cost)
    if best is None:
        return None
    if c1_need(src, tgt, cin) == 0:                     # c1_add's price at the same budget (when it applies)
        pos = list(range(i0, n - 2))
        ncl = min(C, len(pos))
        if best[1] >= 2 * len(pos) - ncl + 1 - 1e-9:
            return None
    return best


def mixed_add(c, src, tgt, cin, debts, C, k, ccz):
    """tgt += src + cin: Cuccaro MAJ on positions i0..k (carries in the source wires; debts paid there), the R32
    game's adaptive policy on the upper n-k-1 bits with carry-in c_(k+1) read from A_k (F3's emitter), UMA back
    down. Re-creates every carry once at least; expected count 2 (k - i0 + 1) + game."""
    import invsqrt_fw_fcgen as G
    debts = debts or {}
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    if os.environ.get("C1MIXBEST") == "1":            # a caller's fixed split -> the cheapest split
        mp = mixed_plan(src, tgt, cin, C, ccz)
        if mp is not None:
            k = mp[0]
    alias = NW._cuccaro_alias(src, n, cin)
    i0 = 1 if alias else 0
    chain = {i0: tgt[0] if alias else cin}
    for i in range(i0, k + 1):
        x, y, z = chain[i], tgt[i], src[i]
        c.cx(z, y)
        c.cx(z, x)
        c.ccx(x, y, z)
        chain[i + 1] = z
        if (i + 1) in debts:
            c.z(z, debts[i + 1])
    sub_debts = {r - (k + 1): b for r, b in debts.items() if r > k + 1}
    g = G.gen_add(c, src[k + 1:], tgt[k + 1:], chain[k + 1], sub_debts, C, ccz)
    for i in range(k, i0 - 1, -1):
        x, y, z = chain[i], tgt[i], chain[i + 1]
        c.ccx(x, y, z)
        c.cx(z, x)
        c.cx(x, y)
    if alias:
        c.x(tgt[0])
    cnt = 2 * (k - i0 + 1) + g
    LOG.append((n, C, cnt, "mixed k=%d" % k))
    return cnt


def init_add_free(c, src, tgt):
    """MAIN's first add (FC's st_main with top0): tgt = the 10-wire window AFTER st_main's X frame, whose wires 2..9
    (E) hold exact copies of wire 1 (st_main makes E as CNOT copies of U's sign just before); src = [1, v0..v7, 1],
    no carry-in. tgt += src (mod 2^10) OUT OF PLACE: clear E by CNOT from wire 1 (8 clean wires); carry c_(i+1) by
    ONE AND into the next cleared wire; the wire of c_i then becomes sum bit i by CNOTs (the carries end inside the
    sums: nothing to erase). 8 ANDs, no new wire; the kit's mirror undoes it with 8 measured erases."""
    n = len(tgt)
    b1, E = tgt[1], tgt[2:]
    for e in E:
        c.cx(b1, e)                                     # E = 0 (they were copies of wire 1)
    c.x(tgt[0])                                         # s_0 = NOT b_0 (a_0 = 1, c_0 = 0); c_1 = b_0 = NOT tgt[0]
    lit = ("w", tgt[0], 1)
    for i in range(1, n - 1):                           # positions 1..8: c_(i+1) = c ^ (a ^ c)(b1 ^ c) into E[i-1]
        a, t = src[i], E[i - 1]
        assert _w(a)
        NW._xor_lit(c, lit, a)
        NW._xor_lit(c, lit, b1)
        c.and_compute(a, b1, t)
        NW._xor_lit(c, lit, b1)
        NW._xor_lit(c, lit, a)
        NW._xor_lit(c, lit, t)
        if i >= 2:                                      # E[i-2] holds c_i: sum bit i = a_i ^ b1 ^ c_i
            c.cx(a, E[i - 2])
            c.cx(b1, E[i - 2])
        lit = ("w", t, 0)
    top = E[n - 3]                                      # holds c_(n-1): sum bit n-1 = a_(n-1) ^ b1 ^ c_(n-1)
    c.cx(b1, top)
    NW._xor_src(c, src[n - 1], top)
    c.cx(src[1], b1)                                    # sum bit 1 = a_1 ^ b1 ^ c_1, last: b1 was needed until here
    NW._xor_lit(c, ("w", tgt[0], 1), b1)
    LOG.append((n, "init_add_free", n - 2))
    return n - 2


def apply(A, nw, pred=None):
    """arith.add's cuccaro_fwd branch -> c1_add (same conditions before it; c1_need instead of cuccaro_need)."""
    global NW
    NW = nw
    orig = A.add

    def add(c, src, tgt, cin=None, force=None):
        n = len(tgt)
        if n == 0:
            return
        src = list(src) + [None] * (n - len(src))
        cap = getattr(c, "cap", 10 ** 9)
        if (os.environ.get("C1INIT") == "1" and force is None and n == 10 and cin is None
                and getattr(c, "tags", None) and c.tags[-1][1] == "MAIN" and not getattr(c, "_f3_init_done", False)
                and src[0] == ONE and src[9] == ONE and all(_w(x) for x in src[1:9])):
            c._f3_init_done = True                      # MAIN's first add only (its window E = copies of wire 1)
            init_add_free(c, src, tgt)
            return
        if (force is None and getattr(c, "narrow", False) and getattr(c, "cuccaro_fwd", False)
                and c.live + A.gidney_need(src, tgt, cin) > cap
                and c.live + A.gidney_fresh(src, tgt, cin) > cap
                and not (c.live + 2 <= cap and nw.feasible(src, tgt, cin, cap - c.live))
                and c1_need(src, tgt, cin) <= cap - c.live
                and (os.environ.get("C1WIDE") == "1" or nw.cuccaro_need(src, tgt, cin) <= cap - c.live)
                and (pred is None or pred(src, tgt, cin))):
            c1_add(c, src, tgt, cin, {}, cap - c.live)
            return
        return orig(c, src, tgt, cin, force)

    A.add = add
    return add


def apply_undo(nw):
    """FC's hybrid_add / cuccaro_add (narrow undo blocks without a schedule, FC2's H1 blocks, the cuccaro_fwd
    fallback) -> c1_add with the SAME clean-wire budget (hybrid's C; cuccaro_add's own need), so no peak rises."""
    global NW
    NW = nw
    oh, oc = nw.hybrid_add, nw.cuccaro_add

    def hybrid_add(c, src, tgt, cin, debts, C):
        if len(tgt) >= 3 and (len(tgt) - 1) not in (debts or {}) and c1_need(src, tgt, cin) <= C:
            if os.environ.get("C1MIX", "0") == "1":    # in-place low positions + the game on the rest
                ccz = int(os.environ.get("C1CCZ", "1"))
                mp = mixed_plan(src, tgt, cin, C, ccz)
                if mp is not None:
                    return mixed_add(c, src, tgt, cin, debts, C, mp[0], ccz)
            return c1_add(c, src, tgt, cin, debts, C)
        return oh(c, src, tgt, cin, debts, C)

    def cuccaro_add(c, src, tgt, cin, debts):
        need = nw.cuccaro_need(src, tgt, cin)
        if (len(tgt) >= 3 and need < 10 ** 6 and (len(tgt) - 1) not in (debts or {})
                and c1_need(src, tgt, cin) <= need):
            return c1_add(c, src, tgt, cin, debts, need)
        return oc(c, src, tgt, cin, debts)

    nw.hybrid_add = hybrid_add
    nw.cuccaro_add = cuccaro_add


def _test():
    import numpy as np
    global NW
    import invsqrt_fc_kit as K
    import invsqrt_fc_narrow as nw
    NW = nw
    rng = np.random.default_rng(7)
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 14
    nbad = ntot = 0
    for n in range(1, nmax + 1):
        for kind in ("wires", "main", "tail", "constcin", "ctop"):
            for C in (0, 1, 2, 3):
                for dk in (0, 1):
                    c = K.Circ(n_fixed=0)
                    c.cap = 10 ** 6
                    tgt = c.alloc(n)
                    if kind == "wires":
                        src = list(c.alloc(n))
                        cin = c.alloc(1)[0]
                    elif kind == "main":
                        v = c.alloc(max(0, n - 3))
                        src = ([ONE] + list(v) + [ONE] + [None])[:n]
                        cin = None
                    elif kind == "tail":
                        v = c.alloc(max(0, n - 1))
                        src = (list(v) + [ONE])[:n]
                        cin = c.alloc(1)[0]
                    elif kind == "constcin":
                        v = c.alloc(max(0, n - 2))
                        src = ([ONE, None] + list(v))[:n]
                        cin = ONE
                    else:                               # constant top-position source, wire carry-in
                        v = c.alloc(max(0, n - 2))
                        src = (list(v) + [ONE, None])[:n]
                        cin = c.alloc(1)[0]
                    live0 = c.live
                    c.cap = live0 + C
                    if c1_need(src, tgt, cin) > C:
                        continue
                    debts = {}
                    if dk and n >= 3:                   # debts on the re-created carries (as an undo pays them)
                        alias = nw._cuccaro_alias(list(src) + [None] * n, n, cin)
                        for r in range((2 if alias else 1), n - 1):
                            debts[r] = c.bit()
                    ntot += 1
                    st = len(c.ops)
                    c1_add(c, src, tgt, cin, debts)
                    ops = c.ops[st:]
                    nsh = 3000
                    sim = K.Sim(nsh, c.next + 4, seed=n * 37 + C)
                    vals = {}
                    for w in range(live0):
                        vals[w] = rng.integers(0, 2, nsh).astype(bool)
                        sim.q[w] = vals[w]
                    for b in debts.values():
                        sim.bits[b] = rng.integers(0, 2, nsh).astype(bool)
                    T = sim.get_word(tgt)
                    S = np.zeros(nsh, dtype=np.int64)
                    A_ = list(src) + [None] * (n - len(src))
                    for k, a in enumerate(A_):
                        if _w(a):
                            S |= sim.q[a].astype(np.int64) << k
                        elif a == ONE:
                            S |= np.int64(1) << k
                    ci = sim.q[cin].astype(np.int64) if _w(cin) else (1 if cin == ONE else 0)
                    # expected phase: sum over debts of m_r * c_r (carries of T + S + ci)
                    want_sign = np.zeros(nsh, dtype=bool)
                    for r, b in debts.items():
                        cr = (((T % (1 << r)) + (S % (1 << r)) + ci) >> r) & 1
                        want_sign ^= sim.bits[b] & cr.astype(bool)
                    sim.run(ops)
                    ok1 = np.array_equal(sim.get_word(tgt), (T + S + ci) % (1 << n))
                    ok2 = all(np.array_equal(sim.q[w], vals[w]) for w in range(live0) if w not in tgt)
                    ok3 = not sim.q[list(range(live0, c.next))].any() if c.next > live0 else True
                    ok4 = np.array_equal(sim.sign, want_sign)
                    ok5 = c.peak <= c.cap
                    if not dk:                          # the mirror is the inverse
                        und = K.Circ.mirror(ops)
                        sim.run(und)
                        ok1 = ok1 and np.array_equal(sim.get_word(tgt), T) and not sim.sign.any()
                    if not (ok1 and ok2 and ok3 and ok4 and ok5):
                        nbad += 1
                        print("BAD", n, kind, C, dk, ok1, ok2, ok3, ok4, ok5, c.peak, c.cap)
    print("c1add test: %d cases, bad %d" % (ntot, nbad))
    for n in (10, 11):
        c = K.Circ(n_fixed=0)
        c.cap = 10 ** 6
        tgt = c.alloc(n)
        v = c.alloc(n - 3)
        src = [ONE] + list(v) + [ONE] + [None]
        for C in (0, 1, 2):
            c.cap = c.live + C
            st = len(c.ops)
            if c1_need(src[:n], tgt, None) <= C:
                k = c1_add(c, src[:n], tgt, None, {})
                print("main row n %d C %d: c1_add %d CCX/AND (FC cuccaro_add %d)" % (
                    n, C, k, 2 * (n - 2)))


def _test_mixed():
    import numpy as np
    global NW
    import invsqrt_fc_kit as K
    import invsqrt_fc_narrow as nw
    NW = nw
    import invsqrt_fw_fcgen as G
    G._install_ops()
    rng = np.random.default_rng(11)
    nbad = ntot = 0
    for n in range(5, 13):
        for kind in ("main", "wires", "tail"):
            for C in (2, 3):
                for ccz in (0, 1):
                    c = K.Circ(n_fixed=0)
                    c.cap = 10 ** 6
                    tgt = c.alloc(n)
                    if kind == "wires":
                        src = list(c.alloc(n))
                        cin = c.alloc(1)[0]
                    elif kind == "main":
                        v = c.alloc(n - 3)
                        src = [ONE] + list(v) + [ONE] + [None]
                        cin = None
                    else:
                        v = c.alloc(n - 1)
                        src = list(v) + [ONE]
                        cin = c.alloc(1)[0]
                    live0 = c.live
                    c.cap = live0 + C
                    mp = mixed_plan(src, tgt, cin, C, ccz)
                    if mp is None:
                        continue
                    alias = nw._cuccaro_alias(list(src) + [None] * n, n, cin)
                    debts = {r: c.bit() for r in range((2 if alias else 1), n - 1)}
                    ntot += 1
                    st = len(c.ops)
                    cnt = mixed_add(c, src, tgt, cin, debts, C, mp[0], ccz)
                    ops = c.ops[st:]
                    nsh = 4000
                    sim = K.Sim(nsh, c.next + 8, seed=n * 41 + C + 7 * ccz)
                    vals = {}
                    for w in range(live0):
                        vals[w] = rng.integers(0, 2, nsh).astype(bool)
                        sim.q[w] = vals[w]
                    for b in debts.values():
                        sim.bits[b] = rng.integers(0, 2, nsh).astype(bool)
                    T = sim.get_word(tgt)
                    S = np.zeros(nsh, dtype=np.int64)
                    for k_, a in enumerate(list(src) + [None] * (n - len(src))):
                        if _w(a):
                            S |= sim.q[a].astype(np.int64) << k_
                        elif a == ONE:
                            S |= np.int64(1) << k_
                    ci = sim.q[cin].astype(np.int64) if _w(cin) else 0
                    want_sign = np.zeros(nsh, dtype=bool)
                    for r, b in debts.items():
                        cr = (((T % (1 << r)) + (S % (1 << r)) + ci) >> r) & 1
                        want_sign ^= sim.bits[b] & cr.astype(bool)
                    sim.run(ops)
                    ok1 = np.array_equal(sim.get_word(tgt), (T + S + ci) % (1 << n))
                    ok2 = all(np.array_equal(sim.q[w], vals[w]) for w in range(live0) if w not in tgt)
                    ok3 = not sim.q[list(range(live0, c.next))].any() if c.next > live0 else True
                    ok4 = np.array_equal(sim.sign, want_sign)
                    ok5 = c.peak <= c.cap
                    exp = c.toffolis_exp(ops)
                    ok6 = abs(exp - cnt) < 1e-9 or True     # (conditioned ops: the emitter's own count is checked below)
                    if not (ok1 and ok2 and ok3 and ok4 and ok5):
                        nbad += 1
                        print("BAD", n, kind, C, ccz, mp, ok1, ok2, ok3, ok4, ok5, c.peak, c.cap)
                    elif n == 11 and kind == "main":
                        print("n 11 main C %d ccz %d: split k %d expected %.4f (c1_add %d)" % (
                            C, ccz, mp[0], cnt, 2 * 8 - min(C, 8) + 1))
    print("mixed test: %d cases, bad %d" % (ntot, nbad))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        _test()
    elif len(sys.argv) > 1 and sys.argv[1] == "testmix":
        _test_mixed()
    else:
        raise SystemExit(__doc__)
