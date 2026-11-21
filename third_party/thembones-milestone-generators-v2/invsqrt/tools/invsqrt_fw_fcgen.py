#!/usr/bin/env python3
"""invsqrt_fw_fcgen.py -- F3 (ninth holder, 2026-11-14): GENERAL emitter of the R32 game's optimal adaptive policy
for FC's narrow undo blocks (any game size L, clean wires C, with or without the CCZ repayment move).

The policy is the re-pricer's game (tools/r32_inv_clean_rp_game.py: solve(n = L + 2, C, feed 0, ccz)), followed
move by move as its own `trace` picks them. Emission is a depth-first walk of the policy tree: every GHOST is an
X-measurement conditioned on the current node's predicate bit P; its two children get the derived bits P & ~b and
P & b; every op of a node is conditioned on that node's bit (unconditioned at the root). The C clean wires are
allocated once as slots; each path keeps its own carry -> slot map (paths exclude each other, so slots are shared).
Expected ANDs = the game's value (a conditioned AND fires on its path's shots only).

CCZ repayment of a pending carry c_k while c_(k-2) is on a wire (FC's s-form, c_k = c_(k-1) + x y, x = a + c_(k-1),
y = a + s at position k-1):  (-1)^(c_k) = Z(c_(k-2)) CZ(x', y') CZ(a_(k-1) + c_(k-2), y) CCZ(x', y', y), with
x' = a_(k-2) + c_(k-2), y' = a_(k-2) + s_(k-2), y = a_(k-1) + s_(k-1).

  python3 invsqrt_fw_fcgen.py test L C ccz [shots]     exact random-value test on FC's block kinds
  apply(shapes)                                         route FC's game_add for the listed (L, C, ccz) to this emitter
"""
import functools
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_fw_fcpol as P   # noqa: E402  (kit extension: BAND / BXOR / XANDE; FC's modules)

K, NW = P.K, P.NW
_RP = {}


def _rp():
    if "m" not in _RP:
        spec = importlib.util.spec_from_file_location("r32rp", os.path.join(HERE, "r32_inv_clean_rp_game.py"))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _RP["m"] = mod
    return _RP["m"]


@functools.lru_cache(maxsize=None)
def solved(L, C, ccz):
    val, V, tabs = _rp().solve(L + 2, C, 0, ccz)
    return val, V, tabs


# ---------------------------------------------------------------- more kit ops (conditioned HMR / CCX / CCZ)
def _xhmr(self, q, g):
    b = self.bit()
    self.ops.append(("XHMR", q, b, g))
    return b


def _xccx(self, a, b, t, g):
    self.ops.append(("XCCX", a, b, t, g))


def _xccz(self, a, b, t, g):
    self.ops.append(("XCCZ", a, b, t, g))


def _install_ops():
    if getattr(K.Circ, "xhmr", None) is not None:
        return
    P.apply_kit()
    K.Circ.xhmr = _xhmr
    K.Circ.xccx = _xccx
    K.Circ.xccz = _xccz
    orig_run = K.Sim.run

    def run(self, ops, start=0, stop=None):
        seg = ops[start:stop]
        i0 = 0
        for i, o in enumerate(seg):
            if o[0] in ("XHMR", "XCCX", "XCCZ"):
                if i > i0:
                    orig_run(self, seg[i0:i])
                if o[0] == "XHMR":
                    _, q, b, g = o
                    gv = self.bits[g]
                    m = self.rng.integers(0, 2, self.n).astype(bool) & gv
                    self.bits[b] = m
                    self.sign ^= self.q[q] & m
                    self.q[q] = self.q[q] & ~gv
                elif o[0] == "XCCX":
                    _, a, b, t, g = o
                    self.q[t] ^= self.q[a] & self.q[b] & self.bits[g]
                else:
                    _, a, b, t, g = o
                    self.sign ^= self.q[a] & self.q[b] & self.q[t] & self.bits[g]
                i0 = i + 1
        if i0 < len(seg):
            orig_run(self, seg[i0:])
    K.Sim.run = run
    orig_kmx = K.kmx_text

    def kmx_text(circ, x_wires, qmax=None):
        saved = circ.ops
        texts, ops2 = [], []
        for o in saved:
            if o[0] in ("XHMR", "XCCX", "XCCZ"):
                if o[0] == "XHMR":
                    t = f"PUSH_CONDITION if b{o[3]}\nHMR q{o[1]} b{o[2]}\nPOP_CONDITION"
                elif o[0] == "XCCX":
                    t = f"CCX q{o[1]} q{o[2]} q{o[3]} if b{o[4]}"
                else:
                    t = f"CCZ q{o[1]} q{o[2]} q{o[3]} if b{o[4]}"
                texts.append(t)
                ops2.append(("X", -500000 - (len(texts) - 1), None))
            else:
                ops2.append(o)
        circ.ops = ops2
        try:
            txt = orig_kmx(circ, x_wires, qmax)
        finally:
            circ.ops = saved
        for k, t in enumerate(texts):
            txt = txt.replace(f"X q{-500000 - k}\n", t + "\n", 1)
        return txt
    K.kmx_text = kmx_text
    orig_exp = K.Circ.toffolis_exp

    def toffolis_exp(self, ops=None):
        ops_ = self.ops if ops is None else ops
        bp = getattr(self, "bprob", {})
        return orig_exp(self, ops) + sum(bp.get(o[4], 0.5) for o in ops_ if o[0] in ("XCCX", "XCCZ"))
    K.Circ.toffolis_exp = toffolis_exp


def _bnode(c, Pb, lit):
    """derived bit for (P AND literal); P None = the root."""
    b, pol = lit
    if Pb is None:
        if pol:
            return b
        d = c.bit()
        c.ops.append(("BAND", d, ((b, 0),)))
        c.bprob = getattr(c, "bprob", {})
        c.bprob[d] = 0.5
        return d
    return c.band([(Pb, 1), (b, pol)])


# ---------------------------------------------------------------- the emitter
def gen_add(c, src, tgt, cin, debts, C, ccz):
    _install_ops()
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    base = {0: NW._lit(cin)}

    def xsrc(a, w, g):
        if NW._w(a):
            c.cx(a, w, g)
        elif a == NW.ONE:
            c.x(w, g)

    def write_sum(r, carry, g):
        xsrc(src[r], tgt[r], g)
        NW._xor_lit(c, carry[r], tgt[r], g)

    write_sum(0, base, None)
    F = 0
    while F + 1 <= n - 2:
        f = NW._form_s(src[F], tgt[F], base[F])
        if f[0] == "and":
            break
        base[F + 1] = ("k", f[1]) if f[0] == "const" else ("w", f[1], f[2])
        F += 1
        write_sum(F, base, None)
    L = n - 2 - F
    val, V, tabs = solved(L, C, ccz)
    assert val < 1e6, ("no policy", L, C, ccz)
    rp = _rp()
    _, pow3, off = tabs
    slots = c.alloc(C)

    def form(R, carry):
        return NW._form_s(src[R - 1], tgt[R - 1], carry[R - 1])

    def zlit(lit, g):
        if lit[0] == "k":
            if lit[1]:
                c.neg(g)
        else:
            c.z(lit[1], g)
            if lit[2]:
                c.neg(g)

    def op_place(k, first, g, carry, smap):
        R = F + k
        s = next(x for x in slots if x not in smap.values())
        f = form(R, carry)
        assert f[0] == "and"
        _, xw, yw, ops = f
        NW._emit(c, ops)
        if g is None:
            c.and_compute(xw, yw, s)
        else:
            c.cand_compute(xw, yw, s, g)
        NW._emit(c, ops[::-1])
        NW._xor_lit(c, carry[R - 1], s, g)
        carry[R] = ("w", s, 0)
        smap[k] = s
        if first:
            write_sum(R, carry, g)
            if R in debts:
                c.z(s, debts[R] if g is None else c.band([(g, 1), (debts[R], 1)]))

    def op_remove(k, g, carry, smap):
        R = F + k
        s = smap.pop(k)
        NW._xor_lit(c, carry[R - 1], s, g)
        f = form(R, carry)
        _, xw, yw, ops = f
        NW._emit(c, ops)
        if g is None:
            c.and_erase(xw, yw, s, free=False)
        else:
            c.xand_erase(xw, yw, s, g)
        NW._emit(c, ops[::-1])
        del carry[R]

    def op_ghost(k, g, carry, smap):
        R = F + k
        s = smap.pop(k)
        b = c.hmr(s) if g is None else c.xhmr(s, g)
        del carry[R]
        return b

    def op_repay(k, g, carry):
        R = F + k
        zlit(carry[R - 1], g)
        f = form(R, carry)
        _, xw, yw, ops = f
        NW._emit(c, ops)
        c.cz(xw, yw, g)
        NW._emit(c, ops[::-1])

    def op_cczrepay(k, g, carry):
        """(-1)^(c_k) = Z(c_(k-2)) CZ(x', y') CZ(a_(k-1) + c_(k-2), y) CCZ(x', y', y)."""
        R = F + k
        lo2 = carry[R - 2]
        zlit(lo2, g)
        f2 = NW._form_s(src[R - 2], tgt[R - 2], lo2)
        assert f2[0] == "and", f2
        _, xw2, yw2, ops2 = f2
        a1, s1 = src[R - 1], tgt[R - 1]
        # y = a_(k-1) + s_(k-1) on wire s1
        xsrc(a1, s1, None)
        NW._emit(c, ops2)
        c.cz(xw2, yw2, g)
        if g is None:
            c.ccz(xw2, yw2, s1)
        else:
            c.xccz(xw2, yw2, s1, g)
        NW._emit(c, ops2[::-1])
        # CZ(a_(k-1) + c_(k-2), y)
        if NW._w(a1):
            NW._xor_lit(c, lo2, a1, None)
            c.cz(a1, s1, g)
            NW._xor_lit(c, lo2, a1, None)
        else:
            ka = 1 if a1 == NW.ONE else 0
            if lo2[0] == "k":
                if ka ^ lo2[1]:
                    c.z(s1, g)
            else:
                c.cz(lo2[1], s1, g)
                if ka ^ lo2[2]:
                    c.z(s1, g)
        xsrc(a1, s1, None)

    def op_top(g, carry):
        tw = tgt[n - 1]
        f = NW._form_s(src[n - 2], tgt[n - 2], carry[n - 2])
        if f[0] == "const":
            if f[1]:
                c.x(tw, g)
        elif f[0] == "alias":
            c.cx(f[1], tw, g)
            if f[2]:
                c.x(tw, g)
        else:
            _, xw, yw, ops = f
            NW._emit(c, ops)
            if g is None:
                c.ccx(xw, yw, tw)
            else:
                c.xccx(xw, yw, tw, g)
            NW._emit(c, ops[::-1])
            NW._xor_lit(c, carry[n - 2], tw, g)
        xsrc(src[n - 1], tw, g)

    placed = set()

    def walk(state, g, carry, smap, cost_acc):
        """follow the policy from `state` (m, cfg, f, t) under predicate bit g (None = root) until the end."""
        m, cfg, fz, t = state
        for _ in range(4000):
            d = [(cfg // int(pow3[i])) % 3 for i in range(m)]
            wires = sum(1 for x in d if x == 1)
            bse = int(off[m]) + cfg
            if m == L and t == 1 and cfg == 0:
                return
            cands = []
            if m == L and t == 0 and L > 0 and d[L - 1] == 1:
                cands.append((1.0 + V[bse * 4 + fz * 2 + 1], "top", 0, (m, cfg, fz, 1)))
            if m < L and wires < C and ((fz == 1) if m == 0 else d[m - 1] == 1):
                cands.append((1.0 + V[(int(off[m + 1]) + cfg + int(pow3[m])) * 4 + fz * 2], "place", m + 1,
                              (m + 1, cfg + int(pow3[m]), fz, 0)))
            for i in range(m):
                pred = (fz == 1) if i == 0 else (d[i - 1] == 1)
                p3 = int(pow3[i])
                if d[i] == 1:
                    if pred:
                        cands.append((V[(bse - p3) * 4 + fz * 2 + t], "remove", i + 1, (m, cfg - p3, fz, t)))
                    else:
                        e = 0.5 * V[(bse - p3) * 4 + fz * 2 + t] + 0.5 * V[(bse + p3) * 4 + fz * 2 + t]
                        cands.append((e, "ghost", i + 1, ((m, cfg - p3, fz, t), (m, cfg + p3, fz, t))))
                elif d[i] == 2:
                    if pred:
                        cands.append((V[(bse - 2 * p3) * 4 + fz * 2 + t], "repay", i + 1, (m, cfg - 2 * p3, fz, t)))
                    elif ccz:
                        anc = True if i == 0 else ((fz == 1) if i == 1 else d[i - 2] == 1)
                        if anc:
                            cands.append((1.0 + V[(bse - 2 * p3) * 4 + fz * 2 + t], "cczrepay", i + 1,
                                          (m, cfg - 2 * p3, fz, t)))
                else:
                    if pred and wires < C:
                        cands.append((1.0 + V[(bse + p3) * 4 + fz * 2 + t], "replace", i + 1, (m, cfg + p3, fz, t)))
            cands.sort(key=lambda z: (round(z[0], 9), 0 if (z[1] in ("remove", "repay")) else 1))
            _, mv, k, nxt = cands[0]
            if mv == "top":
                op_top(g, carry)
            elif mv in ("place", "replace"):
                first = k not in placed
                if first:
                    placed.add(k)
                op_place(k, first, g, carry, smap)
            elif mv == "remove":
                op_remove(k, g, carry, smap)
            elif mv == "repay":
                op_repay(k, g, carry)
            elif mv == "cczrepay":
                assert k >= 2, "CCZ repayment of c1 (carry-in form) is not implemented"
                op_cczrepay(k, g, carry)
            else:   # ghost: two children
                b = op_ghost(k, g, carry, smap)
                g0 = _bnode(c, g, (b, 0))
                g1 = _bnode(c, g, (b, 1))
                snap_placed = set(placed)
                walk(nxt[0], g0, dict(carry), dict(smap), cost_acc)
                placed.clear()
                placed.update(snap_placed)
                walk(nxt[1], g1, dict(carry), dict(smap), cost_acc)
                return
            m, cfg, fz, t = nxt
        raise AssertionError("policy walk did not end")

    walk((0, 0, 1, 0), None, dict(base), {}, [0.0])
    c.free_(slots)
    return float(val)


def apply(shapes):
    """route FC's game_add for every (L, C) in `shapes` (dict (L, C) -> ccz flag) to gen_add; FC's `plan` reports
    the game's value there, so that FC's `block` takes game_add even where it had no plan (C 2: Cuccaro before)."""
    _install_ops()
    orig = NW.game_add
    orig_plan = NW.plan

    def plan(L, C):
        if (L, C) in shapes:
            return (("gen",), solved(L, C, shapes[(L, C)])[0])
        return orig_plan(L, C)
    NW.plan = plan

    def game_add(c, src, tgt, cin, debts, C):
        n = len(tgt)
        L = n - 2 - NW.free_prefix(src, n, cin) if n >= 3 else -1
        if (L, C) in shapes:
            return gen_add(c, src, tgt, cin, debts, C, shapes[(L, C)])
        return orig(c, src, tgt, cin, debts, C)
    NW.game_add = game_add


def _test(L, C, ccz, nsh=4000):
    import numpy as np
    _install_ops()
    rng = np.random.default_rng(13)
    ntest = nbad = 0
    for n in range(L + 2, L + 6):
        for kind in ("main", "wires", "tail", "constcin"):
            c = K.Circ(n_fixed=0)
            tgt = c.alloc(n)
            if kind == "wires":
                src = list(c.alloc(n)); cin = c.alloc(1)[0]
            elif kind == "main":
                v = c.alloc(max(0, n - 3)); src = ([NW.ONE] + list(v) + [NW.ONE] + [None])[:n]; cin = None
            elif kind == "tail":
                v = c.alloc(max(0, n - 1)); src = (list(v) + [NW.ONE])[:n]; cin = c.alloc(1)[0]
            else:
                v = c.alloc(max(0, n - 2)); src = ([NW.ONE, None] + list(v))[:n]; cin = NW.ONE
            if n - 2 - NW.free_prefix(src, n, cin) != L:
                continue
            live0 = c.live
            c.cap = live0 + C
            debts = NW.fwd(c, src, tgt, cin)
            fwd_ops = list(c.ops)
            mark = len(c.ops)
            for w in tgt:
                c.x(w)
            cost = gen_add(c, src, tgt, cin, debts, C, ccz)
            for w in tgt:
                c.x(w)
            undo = c.ops[mark:]
            sim = K.Sim(nsh, c.next + 4, seed=n * 5 + C)
            vals = {}
            for w in range(live0):
                vals[w] = rng.integers(0, 2, nsh).astype(bool)
                sim.q[w] = vals[w]
            T = sim.get_word(tgt)
            S = np.zeros(nsh, dtype=np.int64)
            for kk, a in enumerate(src):
                if NW._w(a):
                    S |= sim.q[a].astype(np.int64) << kk
                elif a == NW.ONE:
                    S |= np.int64(1) << kk
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
            print(f"L {L} C {C} ccz {ccz} n {n} {kind:8s}: fwd {ok1} undo {ok2} anc0 {ok3} phase {ok4} peak "
                  f"{c.peak}/{c.cap} expected {exp:.4f} (game {cost:.4f}) ops {len(undo)}", flush=True)
    print(f"fcgen test L {L} C {C} ccz {ccz}: {ntest} shapes, bad {nbad}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        _test(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]) if len(sys.argv) > 5 else 4000)
