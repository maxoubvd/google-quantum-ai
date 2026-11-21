#!/usr/bin/env python3
"""invsqrt_fw_m7gv.py -- F3 (twelfth holder, 2026-11-15): the VIRTUAL GUARD for the undo games of the m 7 stage's
word adds (FC2's h1c, cfg h1c "gvirt": 1). The guard Gd (R's position 0) always holds an XOR of T's wires (the
frame's `par`). Inside an undo game (FC's narrow game_add on NOT R), target bit 0 is restored by the game's first
step and is then CONSTANT (the pre-add guard, or its complement for an addition): it is read again only by the moves
of carry 1 (place / remove / repay) and by the prefix. Between its read windows it is cleared by CNOTs from T's wires
and its wire is lent to the game, which then runs the R32 plan at C + 1 clean wires. At the end the guard is formed
again on its own wire. Exact: every op on the guard is a CNOT / X from wires whose content the frame tracks.
install(NW) wraps NW.game_add; the frame sets GV = {"on", "wire", "form"} around each word add (form = (wires, neg):
the guard's in-game value = XOR(wires) xor neg).
"""
import copy

GV = {"on": False, "wire": None, "form": None}
LOG = []


def install(NW):
    if getattr(NW, "_f3_gv", False):
        return
    orig = NW.game_add

    def ga(c, src, tgt, cin, debts, C):
        if GV.get("on") and len(tgt) and tgt[0] == GV["wire"]:
            r = game_add_virt(c, NW, src, tgt, cin, debts, C, GV["form"], orig)
            if r is not None:
                return r
        return orig(c, src, tgt, cin, debts, C)
    NW.game_add = ga
    NW._f3_gv = True


def _scratch(c):
    d = copy.copy(c)
    d.ops = []
    d.free = list(c.free)
    d.tags = list(c.tags)
    d.stage_peaks = dict(c.stage_peaks)
    return d


def game_add_virt(c, NW, src, tgt, cin, debts, C, form, orig):
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    F0 = NW.free_prefix(src, n, cin)
    L = n - 2 - F0
    if L < 2:
        return None
    pv = NW.plan(L, C + 1)
    if pv[0] is None:
        return None
    # price of the original route at (L, C) (FC's plan or F3's adaptive emitters), on a scratch circuit
    try:
        sc_ = _scratch(c)
        oc = float(orig(sc_, list(src), list(tgt), cin, dict(debts), C))
        free_after = list(sc_.free)                    # the original's allocator order after its game
    except Exception:
        return None
    if pv[1] >= oc - 1e-9:
        return None
    wires, neg = form
    g0 = tgt[0]

    def form_on(w):
        for u in wires:
            c_.cx(u, w)
        if neg:
            c_.x(w)

    # the run, parameterised by a release plan: rel = set of read indices after which the guard is released
    def run(cc, rel, rec):
        global c_
        c_ = cc
        st = {"w": g0, "live": True, "reads": 0}
        carry = {0: NW._lit(cin)}
        fresh = set()

        def T0(extra=0):
            if not st["live"]:
                if g0 in cc.free:                      # its own wire when free, else any free wire
                    cc.free.remove(g0)
                    cc.live += 1
                    cc.peak = max(cc.peak, cc.live)
                    st["w"] = g0
                else:
                    st["w"] = cc.alloc(1)[0]
                form_on(st["w"])
                st["live"] = True
            if rec is not None:
                rec.append(("R", len(fresh) + extra))
            return st["w"]

        def done0():
            k = st["reads"]
            st["reads"] += 1
            if k in rel:
                form_on(st["w"])                       # clear (self-inverse form)
                cc.free_([st["w"]])
                st["live"] = False

        def tg(r, extra=0):
            return T0(extra) if r == 0 else tgt[r]

        def fin(r):
            if r == 0:
                done0()

        # write_sum(0): the first step restores bit 0 (the guard is live on its own wire here)
        NW._xor_src(cc, src[0], tgt[0])
        NW._xor_lit(cc, carry[0], tgt[0])

        def write_sum(r):
            w = tg(r)
            NW._xor_src(cc, src[r], w)
            NW._xor_lit(cc, carry[r], w)
            fin(r)

        F = 0
        while F + 1 <= n - 2:
            f = NW._form_s(src[F], tg(F), carry[F])
            if f[0] == "and":
                fin(F)
                break
            carry[F + 1] = ("k", f[1]) if f[0] == "const" else ("w", f[1], f[2])
            fin(F)
            F += 1
            write_sum(F)
        assert F == F0

        def place(r, first, cond=None):
            f = NW._form_s(src[r - 1], tg(r - 1, 1), carry[r - 1])
            assert f[0] == "and"
            _, xw, yw, ops = f
            t = cc.alloc(1)[0]
            if rec is not None:
                rec.append(("A", len(fresh) + 1))
            NW._emit(cc, ops)
            if cond is None:
                cc.and_compute(xw, yw, t)
            else:
                cc.cand_compute(xw, yw, t, cond)
            NW._emit(cc, ops[::-1])
            fin(r - 1)
            NW._xor_lit(cc, carry[r - 1], t, cond)
            carry[r] = ("w", t, 0)
            fresh.add(r)
            if first:
                write_sum(r)
                if r in debts:
                    cc.z(t, debts[r])

        def remove(r, cond=None):
            t = carry[r][1]
            NW._xor_lit(cc, carry[r - 1], t, cond)
            f = NW._form_s(src[r - 1], tg(r - 1), carry[r - 1])
            _, xw, yw, ops = f
            NW._emit(cc, ops)
            if cond is None:
                cc.and_erase(xw, yw, t)
            else:
                cc.cand_erase(xw, yw, t, cond)
                cc.free_([t])
            NW._emit(cc, ops[::-1])
            fin(r - 1)
            del carry[r]
            fresh.discard(r)

        def ghost(r):
            t = carry[r][1]
            g = cc.hmr(t)
            cc.free_([t])
            del carry[r]
            fresh.discard(r)
            return g

        def repay(r, g):
            lo = carry[r - 1]
            if lo[0] == "k":
                if lo[1]:
                    cc.neg(g)
            else:
                cc.z(lo[1], g)
                if lo[2]:
                    cc.neg(g)
            f = NW._form_s(src[r - 1], tg(r - 1), carry[r - 1])
            _, xw, yw, ops = f
            NW._emit(cc, ops)
            cc.cz(xw, yw, g)
            NW._emit(cc, ops[::-1])
            fin(r - 1)

        def top():
            tw = tgt[n - 1]
            f = NW._form_s(src[n - 2], tg(n - 2), carry[n - 2])
            if f[0] == "const":
                if f[1]:
                    cc.x(tw)
            elif f[0] == "alias":
                cc.cx(f[1], tw)
                if f[2]:
                    cc.x(tw)
            else:
                _, xw, yw, ops = f
                NW._emit(cc, ops)
                cc.ccx(xw, yw, tw)
                NW._emit(cc, ops[::-1])
                NW._xor_lit(cc, carry[n - 2], tw)
            fin(n - 2)
            NW._xor_src(cc, src[n - 1], tw)

        moves = pv[0]
        cb = {"place": lambda k, first, cond: place(F + k, first, cond),
              "remove": lambda k, cond: remove(F + k, cond),
              "ghost": lambda k: ghost(F + k),
              "repay": lambda k, b: repay(F + k, b),
              "top": top,
              "mark": lambda m: cc.ops.append(m)}
        got = NW.run_moves(L, C + 1, list(moves), cb)
        assert got is not None and abs(got - pv[1]) < 1e-9
        assert not fresh
        if not st["live"] or st["w"] != g0:            # the guard back on its own wire
            if st["live"]:
                form_on(st["w"])
                cc.free_([st["w"]])
            assert g0 in cc.free
            cc.free.remove(g0)
            cc.live += 1
            form_on(g0)
        return got, st["reads"]

    # dry run: guard always live; record read events with the number of game carries live at each event
    rec = []
    _, nreads = run(_scratch(c), set(), rec)
    # release after read k when the next read comes after a point with C + 1 carries live (the guard must yield)
    ev = rec
    reads_at = [i for i, e in enumerate(ev) if e[0] == "R"]
    if len(reads_at) != nreads:
        return None
    rel = set()
    for k in range(nreads - 1):
        a, b = reads_at[k], reads_at[k + 1]
        if any(e[0] == "A" and e[1] > C for e in ev[a:b]):
            rel.add(k)
    for k in range(nreads):                            # at every read the carries must leave a wire for the guard
        if ev[reads_at[k]][1] > C:
            return None
    cost, _ = run(c, rel, None)
    if GV.get("keeporder", True) and sorted(c.free) == sorted(free_after):
        c.free[:] = free_after                         # later allocations see the original order
    LOG.append((L, C, oc, cost, len(rel)))
    return cost


def solve(content, target):
    """content {wire: (frozenset of symbols, neg)}, target: a set of symbols (and "one"): (wires, neg) with
    XOR of the wires' contents = target minus "one" and neg = the constant left over."""
    syms = sorted(set().union(*[v[0] for v in content.values()]) | (set(target) - {"one"}))
    idx = {x: i for i, x in enumerate(syms)}
    rows = []                                          # (vector bitmask, wires set, neg)
    for w, (fs, ng) in content.items():
        v = 0
        for x in fs:
            v ^= 1 << idx[x]
        rows.append([v, {w}, ng])
    tv = 0
    for x in set(target) - {"one"}:
        tv ^= 1 << idx[x]
    basis = []
    for r in rows:
        for b in basis:
            if r[0] & (1 << b[3]):
                r[0] ^= b[0]
                r[1] ^= b[1]
                r[2] ^= b[2]
        if r[0]:
            r.append(r[0].bit_length() - 1)
            basis.append(r)
    ws, ng = set(), 0
    for b in basis:
        if tv & (1 << b[3]):
            tv ^= b[0]
            ws ^= b[1]
            ng ^= b[2]
    assert tv == 0, ("guard value not in the span of T's wires", target)
    ng ^= 1 if "one" in target else 0
    return sorted(ws), ng
