"""F1 (third holder): GENERIC endpoint emitter for any router network (stages = matchings on the 26 ring positions,
applied in order, each with an affine control over (f0..f3, p)), any visit order and any hold set, for FC's gate kit
(tools/invsqrt_fc_kit.py, read-only). Features, rotations and the one-token trick as in tools/invsqrt_fk_epkit.py.

Scheduling: every pair is opened right before the first visit of its cone (outer stages first) and closed right after
the last visit of its cone (inner stages first); a HELD pair keeps its AND on a free wire (1 Toffoli), an unheld pair
is a CSWAP at each end (2). Slot p (base position p) sits at ring position slots[p]; word[i] (token i) starts at ring
position i; class k must send token i to slots[(i + k) % 26].
    emit_generic(c, word, flags, net, order, held, slot_ang, wrap_ang, cap, spares=()) -> info
    net = dict(stages=[[(a, b), ...], ...], controls=[terms, ...], slots=[26 ring positions])
Self-test on the public network: python invsqrt_fk_epgen.py [cap]
"""
import sys

sys.path.insert(0, "<HOME>/Quantum/tools")
import invsqrt_fc_kit as K  # noqa: E402
import invsqrt_fk_epkit as E  # noqa: E402
import invsqrt_fk_router as R  # noqa: E402
import invsqrt_fk_rteval as V  # noqa: E402

N = 26


def emit_generic(c, word, flags, net, order, held, slot_ang, wrap_ang, cap, spares=()):
    stages = [[tuple(p) for p in st] for st in net["stages"]]
    controls = net["controls"]
    slots = net["slots"]
    pos_slot = {slots[p]: p for p in range(N)}
    ep = E.Endpoint(c, word, flags, slot_ang, wrap_ang, cap, spares)
    w = ep.word
    t0 = c.toffolis()
    c.x(w[25])
    win = V.windows(stages, order)
    vis = {x: i for i, x in enumerate(order)}
    fvis = [vis[slots[j]] for j in range(9)]
    p_on, p_off = min(fvis), max(fvis)
    one_home = 25
    last_one = max((key for key in win if one_home in key[1]), key=lambda kk: (win[kk][1], -kk[0]))

    def ctl_ctx(s):
        return ep.affine(controls[s])

    def open_pair(s, pr):
        if (s, pr) in held:
            return ep.hold_ctl(controls[s], *pr)
        with ctl_ctx(s) as q:
            c.cswap(q, w[pr[0]], w[pr[1]])
        return None

    def close_pair(s, pr, anc):
        if anc is not None:
            ep.release_ctl(controls[s], *pr, anc)
            return
        if (s, pr) == last_one:                      # one-token trick (the one returns home here for good)
            x = w[pr[0] if pr[1] == one_home else pr[1]]
            one = w[one_home]
            with ctl_ctx(s) as q:
                c.cx(one, x)
                c.x(one)
                m = c.hmr(one)
                c.cz(q, x, cond=m)
                c.x(one)
                c.cx(one, x)
            return
        with ctl_ctx(s) as q:
            c.cswap(q, w[pr[0]], w[pr[1]])

    opens = {}
    for key, (a, b) in win.items():
        opens.setdefault(a, []).append(key)
    closes = {}
    for key, (a, b) in win.items():
        closes.setdefault(b, []).append(key)
    live = {}
    pw = None
    for t in range(N):
        if t == p_on:
            pw = ep.lanes.get()
            ep.form_into(E.P_FORM, pw)
            ep.var["p"] = pw
            ep.vlist[4] = pw
        for key in sorted(opens.get(t, []), key=lambda kk: kk[0]):          # outer stages first
            live[key] = open_pair(*key)
        x = order[t]
        ep.rotated.add(x)
        p = pos_slot[x]
        c.zp(w[x], slot_ang[p])
        if p <= 8:
            ep.feature(p, w[x])
        for key in sorted(closes.get(t, []), key=lambda kk: -kk[0]):        # inner stages first
            close_pair(*key, live.pop(key))
        if t == p_off:
            ep.form_erase(E.P_FORM, pw)
            ep.var["p"] = None
            ep.vlist[4] = None
            ep.lanes.put(pw, ep.spare_ids)
    assert not live and len(ep.rotated) == N and ep.featured == set(range(9))
    c.cz(ep.var["f1"], ep.var["f3"])
    c.x(w[25])
    return dict(toffolis=c.toffolis() - t0, peak_lanes=ep.lanes.peak)


def _hold_ctl(self, terms, i, j):
    a, b = self.word[i], self.word[j]
    anc = self.lanes.get()
    with self.affine(terms) as ctl:
        self.c.cx(a, b)
        self.c.ccx(ctl, b, anc)
        self.c.cx(a, b)
    self.c.cx(anc, a)
    self.c.cx(anc, b)
    return anc


def _release_ctl(self, terms, i, j, anc):
    a, b = self.word[i], self.word[j]
    self.c.cx(anc, a)
    self.c.cx(anc, b)
    m = self.c.hmr(anc)
    with self.affine(terms) as ctl:
        self.c.cz(ctl, a, cond=m)
        self.c.cz(ctl, b, cond=m)
    self.lanes.put(anc, self.spare_ids)


E.Endpoint.hold_ctl = _hold_ctl
E.Endpoint.release_ctl = _release_ctl


def best_holds(stages, order, slots, cap):
    """exact ILP packing (scipy) -> held set."""
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    win = V.windows(stages, order)
    base = V.fixed_demand(order, slots)
    keys = list(win)
    A = np.zeros((N, len(keys)))
    for jj, kk in enumerate(keys):
        a, b = win[kk]
        A[a:b + 1, jj] = 1
    ub = np.array([cap - x for x in base], dtype=float)
    res = milp(-np.ones(len(keys)), constraints=LinearConstraint(A, -np.inf, ub), integrality=np.ones(len(keys)),
               bounds=Bounds(0, 1))
    return {keys[jj] for jj in range(len(keys)) if res.x[jj] > 0.5}


def selftest_public(cap=5, shots=4096, seed=11):
    import random
    import numpy as np
    stages = [E.pairs(cc) for cc in E.CENTRES]
    net = dict(stages=stages, controls=list(E.CONTROLS), slots=list(E.SLOTS))
    order = R.public_walk()
    held = best_holds(stages, order, net["slots"], cap)
    rnd = random.Random(seed)
    c = K.Circ(n_fixed=30)
    word, flags = list(range(26)), [26, 27, 28, 29]
    spares = c.alloc(3)
    slot_ang = [rnd.randrange(E.ONE_TURN) for _ in range(N)]
    wrap_ang = [rnd.randrange(E.ONE_TURN) for _ in range(7)]
    info = emit_generic(c, word, flags, net, order, held, slot_ang, wrap_ang, cap, spares=spares)
    S = K.Sim(shots, max(c.next, 40), seed=seed)
    rng = np.random.default_rng(seed)
    ks = rng.integers(0, 11, shots)
    toks = rng.integers(0, 2, (25, shots)).astype(bool)
    for i in range(25):
        S.q[word[i]] = toks[i]
    for b in range(4):
        S.q[flags[b]] = ((ks >> b) & 1).astype(bool)
    init = S.q.copy()
    S.run(c.ops)
    F = wrap_ang + [E.ONE_TURN // 8, E.ONE_TURN // 4, E.ONE_TURN // 2]
    bad = 0
    for n in range(shots):
        k = int(ks[n])
        tot = 0
        for i in range(26):
            t = 1 if i == 25 else int(toks[i][n])
            if t:
                tot += slot_ang[(i + k) % N] + (F[i + k - N] if i + k >= N else 0)
        got = (int(S.ang[n]) + (E.ONE_TURN // 2 if S.sign[n] else 0)) % E.ONE_TURN
        bad += got != tot % E.ONE_TURN
    print(f"generic emitter, public network, cap {cap}: held {len(held)}, Toffolis {info['toffolis']}, peak lanes "
          f"{info['peak_lanes']}; phase mismatches {bad}, wires changed {int((S.q != init).any(axis=1).sum())}")


if __name__ == "__main__":
    selftest_public(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
