"""F1 (third holder): the invsqrt endpoint (26-slot ring, five controlled reflection stages, path-local walk, wrap
features, slot rotations) as a module for FC's gate kit (tools/invsqrt_fc_kit.py, imported read-only), with the hold
schedule chosen by F1's exact router model (tools/invsqrt_fk_router.py) for any number of free wires.

Algorithm credit: the ring endpoint ("walk8" core, agents/composer/p/w2core_e2.py of the public invsqrt generator,
Copyright 2026 Syed Asad Ali, Apache-2.0). This file is our own code written from that design (parse only).

    emit_endpoint(c, word, flags, slot_ang, wrap_ang, cap, spares=()) -> dict
  word     : 26 wire ids; word[i] holds token i (i < 25: a quotient bit, as a bit); word[25] must be |0> (the
             constant one is made and erased here);
  flags    : 4 wire ids, k = 0..10 in binary (f0 = LSB);
  slot_ang : 26 angles (numerators over 2^62 turns) for base positions p = 0..25;
  wrap_ang : 7 angles for the wrap features j = 0..6 (j = 7: 1/8 turn, 8: 1/4, 9: 1/2 are fixed);
  cap      : wires the endpoint may hold at once besides word / flags (spares first, then fresh allocations);
  spares   : wire ids that are |0> and may be used (returned as |0>).
  Phase: sum_i token_i * (slot_ang[(i + k) % 26] + [i + k >= 26] * F[i + k - 26]), F = wrap_ang + (1/8, 1/4, 1/2).
  Every other wire is returned unchanged; all helpers return to |0>; HMR debts are repaid inside.
Self-test (exact, random tokens, all 11 classes, kit simulator):  python invsqrt_fk_epkit.py [cap ...]
"""
import contextlib
import random
import sys

sys.path.insert(0, "<HOME>/Quantum/tools")
import invsqrt_fc_kit as K  # noqa: E402
import invsqrt_fk_router as R  # noqa: E402

N = 26
CENTRES = (2, 6, 24, 0, 25)
CONTROLS = ((0, 4), (8,), (0, 2, 4, 8), (0, 1, 2), (0, 1))   # affine terms over [f0, f1, f2, f3, p]; 0 = 1
SLOTS = [(p + 21) % N for p in range(N)]                     # ring position of base slot p
POS_SLOT = {SLOTS[p]: p for p in range(N)}
ONE_TURN = 1 << K.ANG_BITS
# quadratic forms as (products [(A, B)], linear L, const); variables by name
ENC = {0: ([("f1 f2 f3 p", "s f0")], "f0", 0), 1: ([("f1 f2 f3 p", "s")], "", 0),
       2: ([("f0", "f1 p"), ("f2 f3", "s")], "", 0), 3: ([("f2 f3", "s")], "", 0),
       6: ([("f0", "p"), ("f3", "s")], "", 0), 7: ([("f3", "s")], "", 0)}
THR = {5: ([("f0 f1", "f2 p")], "f3 p", 0), 6: ([("f1", "f2")], "f3", 0), 9: ([("f0 f1", "f3")], "", 0)}
P_FORM = ([("f1", "f2 f3")], "", 0)


def pairs(c):
    return [(x, (c - x) % N) for x in range(N) if x < (c - x) % N]


def fl(k):
    f = [(k >> i) & 1 for i in range(4)]
    return dict(f0=f[0], f1=f[1], f2=f[2], f3=f[3], p=f[1] & (f[2] ^ f[3]))


def form_value(form, vals):
    prods, lin, const = form
    o = const
    for A, B in prods:
        a = 0
        for v in A.split():
            a ^= vals[v]
        b = 0
        for v in B.split():
            b ^= vals[v]
        o ^= a & b
    for v in lin.split():
        o ^= vals[v]
    return o


def _check_forms():
    for j, form in ENC.items():
        for k in range(11):
            for s in ((1,) if k == j + 1 else (0, 1)):
                vals = dict(fl(k), s=s)
                assert form_value(form, vals) == (int(k >= j + 1) & s), ("ENC", j, k, s)
    for t, form in THR.items():
        for k in range(11):
            assert form_value(form, fl(k)) == int(k >= t), ("THR", t, k)
    for k in range(11):
        v = fl(k)
        assert form_value(P_FORM, v) == v["p"]


_check_forms()


class Lanes:
    def __init__(self, c, spares, cap):
        self.c, self.spare, self.cap, self.used, self.peak = c, list(spares), cap, 0, 0
        self.fresh = []

    def get(self):
        self.used += 1
        self.peak = max(self.peak, self.used)
        assert self.used <= self.cap, ("endpoint lanes", self.used, self.cap)
        if self.spare:
            return self.spare.pop()
        return self.c.alloc(1)[0]

    def put(self, q, spare_ids):
        self.used -= 1
        if q in spare_ids:
            self.spare.append(q)
        else:
            self.c.free_([q])


class Endpoint:
    def __init__(self, c, word, flags, slot_ang, wrap_ang, cap, spares):
        self.c, self.word, self.slot_ang, self.wrap_ang = c, list(word), list(slot_ang), list(wrap_ang)
        self.var = dict(f0=flags[0], f1=flags[1], f2=flags[2], f3=flags[3], p=None)
        self.vlist = list(flags) + [None]
        self.spare_ids = set(spares)
        self.lanes = Lanes(c, spares, cap)
        self.featured = set()
        self.rotated = set()

    # ---- affine controls and forms
    @contextlib.contextmanager
    def affine(self, terms):
        lin = [t.bit_length() - 1 for t in terms if t]
        piv = self.vlist[lin[0]]
        for i in lin[1:]:
            self.c.cx(self.vlist[i], piv)
        if 0 in terms:
            self.c.x(piv)
        try:
            yield piv
        finally:
            if 0 in terms:
                self.c.x(piv)
            for i in reversed(lin[1:]):
                self.c.cx(self.vlist[i], piv)

    def _wires(self, names, s=None):
        out = []
        for v in names.split():
            out.append(s if v == "s" else self.var[v])
        return out

    def _pair(self, A, B):
        """two distinct wires holding XOR(A) and XOR(B) in place; returns (a, b, undo list of cx)."""
        A = list(dict.fromkeys(A))
        B = list(dict.fromkeys(B))
        ops = []
        a = A[0]
        for w in A[1:]:
            self.c.cx(w, a)
            ops.append((w, a))
        T = set(B) - {a}
        if a in B:
            T ^= set(A) - {a}
            T.add(a)                        # a now holds XOR(A); B's a-term = a ^ XOR(A \ {a})
            T.symmetric_difference_update(set())
        cand = sorted(T - {a})
        assert cand, "degenerate product"
        bp = cand[0]
        for w in sorted(T - {bp}):
            self.c.cx(w, bp)
            ops.append((w, bp))
        return a, bp, ops

    def _undo(self, ops):
        for cc, t in reversed(ops):
            self.c.cx(cc, t)

    def form_into(self, form, h, s=None):
        prods, lin, const = form
        for A, B in prods:
            a, b, ops = self._pair(self._wires(A, s), self._wires(B, s))
            self.c.ccx(a, b, h)
            self._undo(ops)
        for w in self._wires(lin, s):
            self.c.cx(w, h)
        if const:
            self.c.x(h)

    def form_erase(self, form, h, s=None):
        prods, lin, const = form
        m = self.c.hmr(h)
        for A, B in prods:
            a, b, ops = self._pair(self._wires(A, s), self._wires(B, s))
            self.c.cz(a, b, cond=m)
            self._undo(ops)
        for w in self._wires(lin, s):
            self.c.z(w, cond=m)
        if const:
            self.c.neg(cond=m)

    # ---- swaps
    def cswap(self, stage, i, j):
        with self.affine(CONTROLS[stage]) as ctl:
            self.c.cswap(ctl, self.word[i], self.word[j])

    def hold(self, stage, i, j):
        a, b = self.word[i], self.word[j]
        anc = self.lanes.get()
        with self.affine(CONTROLS[stage]) as ctl:
            self.c.cx(a, b)
            self.c.ccx(ctl, b, anc)
            self.c.cx(a, b)
        self.c.cx(anc, a)
        self.c.cx(anc, b)
        return anc

    def release(self, stage, i, j, anc):
        a, b = self.word[i], self.word[j]
        self.c.cx(anc, a)
        self.c.cx(anc, b)
        m = self.c.hmr(anc)
        with self.affine(CONTROLS[stage]) as ctl:
            self.c.cz(ctl, a, cond=m)
            self.c.cz(ctl, b, cond=m)
        self.lanes.put(anc, self.spare_ids)

    # ---- rotations and features
    def visit(self, pos):
        p = POS_SLOT[pos]
        assert pos not in self.rotated
        self.rotated.add(pos)
        self.c.zp(self.word[pos], self.slot_ang[p])
        if p <= 8:
            self.feature(p, self.word[pos])

    def feature(self, j, s):
        assert j not in self.featured
        self.featured.add(j)
        if j in ENC:
            h = self.lanes.get()
            self.form_into(ENC[j], h, s)
            self.c.zp(h, self.wrap_ang[j] if j <= 6 else ONE_TURN // 8)
            self.form_erase(ENC[j], h, s)
            self.lanes.put(h, self.spare_ids)
            return
        t = j + 1
        q = self.lanes.get()
        self.form_into(THR[t], q)
        if j <= 6:                                   # generic angle on thr * s: one AND
            h = self.lanes.get()
            self.c.ccx(q, s, h)
            self.c.zp(h, self.wrap_ang[j])
            m = self.c.hmr(h)
            self.c.cz(q, s, cond=m)
            self.lanes.put(h, self.spare_ids)
        else:                                        # j = 8: S on thr * s as three T gates
            T8 = ONE_TURN // 8
            self.c.zp(q, T8)
            self.c.zp(s, T8)
            self.c.cx(q, s)
            self.c.zp(s, -T8)
            self.c.cx(q, s)
        self.form_erase(THR[t], q)
        self.lanes.put(q, self.spare_ids)


def schedule(cap):
    """hold sets from the exact model: S3 / S4 pairs and S2 pairs with windows, earliest-deadline packing under cap
    with the feature helpers and p as fixed demands (S0, S1 never held here)."""
    order = R.public_walk()
    win = R.windows(CENTRES, order)
    base = R.fixed_demand(order)
    cand = [key for key in win if key[0] >= 2]
    free = [cap - d for d in base]
    held = set()
    for key in sorted(cand, key=lambda kk: (win[kk][1], win[kk][0])):
        a, b = win[key]
        if all(free[t] > 0 for t in range(a, b + 1)):
            for t in range(a, b + 1):
                free[t] -= 1
            held.add(key)
    e_hold = {k for k in range(1, 13) if (3, tuple(sorted((13 - k, 13 + k)))) in held}
    f_hold = {k for k in range(13) if (4, tuple(sorted((13 + k, (12 - k) % N)))) in held}
    r2 = {}
    for (s, pr) in held:
        if s == 2:
            a, b = win[(s, pr)]
            r2[pr] = (a // 2, b // 2)
    return e_hold, f_hold, r2


def emit_endpoint(c, word, flags, slot_ang, wrap_ang, cap, spares=(), hold_sets=None):
    E = Endpoint(c, word, flags, slot_ang, wrap_ang, cap, spares)
    e_hold, f_hold, r2 = hold_sets if hold_sets is not None else schedule(cap)
    w = E.word
    t0 = c.toffolis()
    c.x(w[25])                                       # the constant one
    late2 = [pr for pr in pairs(CENTRES[2]) if pr not in r2]
    for pr in pairs(CENTRES[0]):
        E.cswap(0, *pr)
    for pr in pairs(CENTRES[1]):
        E.cswap(1, *pr)
    for pr in late2:
        E.cswap(2, *pr)
    t_glob = c.toffolis() - t0
    e = {k: (13 - k, 13 + k) for k in range(1, 13)}
    f = {k: (13 + k, (12 - k) % N) for k in range(0, 13)}
    held_e, held_r2 = {}, {}

    def e_forward(k):
        if k in e_hold:
            held_e[k] = E.hold(3, *e[k])
        else:
            E.cswap(3, *e[k])

    def e_undo(k):
        if k in e_hold:
            E.release(3, *e[k], held_e.pop(k))
        else:
            E.cswap(3, *e[k])

    def f_round(k):
        i, j = f[k]
        if k in f_hold:
            anc = E.hold(4, i, j)
            E.visit(i)
            E.visit(j)
            E.release(4, i, j, anc)
        else:
            E.cswap(4, i, j)
            E.visit(i)
            E.visit(j)
            E.cswap(4, i, j)

    def r2_open(step):
        for pr, (first, last) in sorted(r2.items()):
            if first == step:
                held_r2[pr] = E.hold(2, *pr)

    def r2_close(step):
        for pr, (first, last) in sorted(r2.items()):
            if last == step:
                E.release(2, *pr, held_r2.pop(pr))

    pw = None
    for k in range(0, 13):
        if k == 8:                                   # p before the first feature (slot 0 is visited at step 8)
            pw = E.lanes.get()
            E.form_into(P_FORM, pw)
            E.var["p"] = pw
            E.vlist[4] = pw
        r2_open(k)
        if k == 0:
            e_forward(1)
            f_round(0)
        else:
            if k + 1 <= 12:
                e_forward(k + 1)
            f_round(k)
            e_undo(k)
        r2_close(k)
    assert not held_e and not held_r2 and len(E.rotated) == N and E.featured == set(range(9))
    c.cz(E.var["f1"], E.var["f3"])                    # j = 9: Z on the one at k = 10 (f1 f3)
    for pr in reversed(late2):
        E.cswap(2, *pr)
    for pr in reversed(pairs(CENTRES[1])):
        E.cswap(1, *pr)
    for pr in reversed(pairs(CENTRES[0])):
        if 25 in pr:                                 # the one-token trick: its last reverse swap costs no AND
            x = w[pr[0] if pr[1] == 25 else pr[1]]
            one = w[25]
            with E.affine(CONTROLS[0]) as ctl:
                c.cx(one, x)
                c.x(one)
                m = c.hmr(one)
                c.cz(ctl, x, cond=m)
                c.x(one)
                c.cx(one, x)
        else:
            E.cswap(0, *pr)
    E.form_erase(P_FORM, pw)
    E.var["p"] = None
    E.vlist[4] = None
    E.lanes.put(pw, E.spare_ids)
    c.x(w[25])
    return dict(toffolis=c.toffolis() - t0, global_fwd=t_glob, peak_lanes=E.lanes.peak, e_hold=sorted(e_hold),
                f_hold=sorted(f_hold), r2=r2)


def selftest(cap, shots=4096, seed=7):
    rnd = random.Random(seed)
    c = K.Circ(n_fixed=26 + 4)
    word = list(range(26))
    flags = [26, 27, 28, 29]
    spares = c.alloc(3)                              # three |0> spares, like the line's retired root bits
    slot_ang = [rnd.randrange(ONE_TURN) for _ in range(N)]
    wrap_ang = [rnd.randrange(ONE_TURN) for _ in range(7)]
    info = emit_endpoint(c, word, flags, slot_ang, wrap_ang, cap, spares=spares)
    S = K.Sim(shots, max(c.next, 40), seed=seed)
    import numpy as np
    rng = np.random.default_rng(seed)
    ks = rng.integers(0, 11, shots)
    toks = rng.integers(0, 2, (25, shots)).astype(bool)
    for i in range(25):
        S.q[word[i]] = toks[i]
    for b in range(4):
        S.q[flags[b]] = ((ks >> b) & 1).astype(bool)
    init = S.q.copy()
    S.run(c.ops)
    F = wrap_ang + [ONE_TURN // 8, ONE_TURN // 4, ONE_TURN // 2]
    exp = np.zeros(shots, dtype=object)
    for n in range(shots):
        k = int(ks[n])
        tot = 0
        for i in range(26):
            t = 1 if i == 25 else int(toks[i][n])
            if t:
                tot += slot_ang[(i + k) % N]
                if i + k >= N:
                    tot += F[i + k - N]
        exp[n] = tot % ONE_TURN
    got = [(int(S.ang[n]) + (ONE_TURN // 2 if S.sign[n] else 0)) % ONE_TURN for n in range(shots)]
    bad_phase = sum(1 for n in range(shots) if got[n] != exp[n])
    bad_wires = int((S.q != init).any(axis=1).sum())
    return info, bad_phase, bad_wires


if __name__ == "__main__":
    caps = [int(x) for x in sys.argv[1:]] or [5]
    for cap in caps:
        info, bp, bw = selftest(cap)
        print(f"cap {cap}: Toffolis {info['toffolis']} (global fwd {info['global_fwd']}), peak lanes "
              f"{info['peak_lanes']}, holds e {len(info['e_hold'])} f {len(info['f_hold'])} S2 {len(info['r2'])}; "
              f"phase mismatches {bp}, wires changed {bw}")
