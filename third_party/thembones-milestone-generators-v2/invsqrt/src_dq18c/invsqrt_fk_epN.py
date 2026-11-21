"""F1 (fourth holder): the invsqrt endpoint for a ring of ANY size N (26 = the line, 27 = FC's family-C word) as a module
for FC's gate kit (tools/invsqrt_fc_kit.py, read-only): five controlled reflection stages, walk, held ANDs on free
wires, slot rotations and the ten wrap features. Network and holds come from tools/invsqrt_fk_ringN.py.
Algorithm credit: the ring endpoint ("walk8" core of the public invsqrt generator, Copyright 2026 Syed Asad Ali,
Apache-2.0); our own code from that design and from F1's router model.

    emit_ring(c, word, flags, net, slot_ang, wrap_ang, cap, spares=(), const_token=None) -> info
  word      : N wire ids, word[i] = token i (a bit). If const_token is given, word[const_token] must be |0>: the
              constant one is made and erased here (one-token trick on its last reverse swap).
  flags     : 4 wire ids, k = 0..10 in binary (f0 = LSB).
  net       : dict(N, centres (5), controls (5 affine term tuples over [f0, f1, f2, f3, p]; 0 = constant 1),
              L (ring position of slot p, p = 0..N-1), order (walk: ring positions in visit order), held (set of
              (stage, pair)))
  slot_ang  : N angles (numerators over 2^64 turns) for slots p = 0..N-1;
  wrap_ang  : 10 generic angles for features j = 0..9, or 7 (then j = 7, 8, 9 are T, S, Z);
  cap       : wires the endpoint may hold at once (spares first, then fresh allocations).
  Phase: sum_i token_i * (slot_ang[(i + k) % N] + [i + k >= N] * F[i + k - N]), F = the feature angles.
  Every other wire is returned unchanged; helpers return to |0>; HMR debts are repaid inside.
Self-test (exact, kit simulator, random tokens, all 11 classes):
    python invsqrt_fk_epN.py N cap [generic|clifford] [const]
"""
import os
import sys

sys.path.insert(0, "<HOME>/Quantum/tools")
import invsqrt_fc_kit as K  # noqa: E402
import invsqrt_fk_epkit as E  # noqa: E402
import invsqrt_fk_epgen  # noqa: E402,F401  (adds Endpoint.hold_ctl / release_ctl)
import invsqrt_fk_ringN as G  # noqa: E402

ONE = E.ONE_TURN
# [k >= t] over (f0, f1, f2, f3, p) on k = 0..10: affine (one product with s) or one product (A AND B) + L, found by
# exhaustive search (fourth holder, 07:55Z); valid for ANY token at the slot (no constant needed)
TFORM = {1: ([("f0", "f1 f2 f3 p")], "f0 f1 f2 f3 p", 0), 2: ([], "f1 f2 f3 p", 0),
         3: ([("f0", "f1 p")], "f2 f3", 0), 4: ([], "f2 f3", 0), 5: ([("f0 f1", "f2 p")], "f3 p", 0),
         6: ([("f1", "f2")], "f3", 0), 7: ([("f0", "p")], "f3", 0), 8: ([], "f3", 0),
         9: ([("f0 f1", "f3")], "", 0), 10: ([("f1", "f2")], "p", 0)}


def _check_tforms():
    for t, form in TFORM.items():
        for k in range(11):
            assert E.form_value(form, E.fl(k)) == int(k >= t), ("TFORM", t, k)


_check_tforms()


def feature_noconst(ep, j, s, wrap_ang):
    """FC profile: any token at slot j. Feature = [k >= j + 1] AND s, rotated (generic) or T / S / Z (j = 7, 8, 9
    when wrap_ang has 7 entries)."""
    c = ep.c
    generic = len(wrap_ang) == 10 or j <= 6
    t = j + 1
    prods, lin, _ = TFORM[t]
    if not generic and j == 9:                             # Z on s * f1 * f3 (k = 10)
        c.ccz(s, ep.var["f1"], ep.var["f3"])
        return
    ang = wrap_ang[j] if generic else ONE // 8
    if not prods:                                          # affine threshold: one AND with s
        h = ep.lanes.get()
        form = ([(lin, "s")], "", 0)
        ep.form_into(form, h, s)
        c.zp(h, ang)
        ep.form_erase(form, h, s)
        ep.lanes.put(h, ep.spare_ids)
        return
    q = ep.lanes.get()
    ep.form_into(TFORM[t], q)
    if generic:
        h = ep.lanes.get()
        c.ccx(q, s, h)
        c.zp(h, ang)
        m = c.hmr(h)
        c.cz(q, s, cond=m)
        ep.lanes.put(h, ep.spare_ids)
    else:                                                  # j = 8: S on q * s as three T gates
        T8 = ONE // 8
        c.zp(q, T8)
        c.zp(s, T8)
        c.cx(q, s)
        c.zp(s, -T8)
        c.cx(q, s)
    ep.form_erase(TFORM[t], q)
    ep.lanes.put(q, ep.spare_ids)


def feature(ep, j, s, wrap_ang):
    """[k >= j + 1] AND s, rotated by its angle. j <= 8 use the line's forms (ENC: one form, THR: threshold + AND)."""
    generic = len(wrap_ang) == 10 or j <= 6
    ang = wrap_ang[j] if generic else (ONE // 8, ONE // 4, ONE // 2)[j - 7]
    c = ep.c
    if j in E.ENC and (generic or j == 7):
        h = ep.lanes.get()
        ep.form_into(E.ENC[j], h, s)
        c.zp(h, ang)
        ep.form_erase(E.ENC[j], h, s)
        ep.lanes.put(h, ep.spare_ids)
        return
    if j == 9:
        f1, f3 = ep.var["f1"], ep.var["f3"]
        if not generic:                                    # Z on s * f1 * f3
            c.ccz(s, f1, f3)
            return
        q = ep.lanes.get()
        c.ccx(f1, f3, q)
        h = ep.lanes.get()
        c.ccx(q, s, h)
        c.zp(h, ang)
        m = c.hmr(h)
        c.cz(q, s, cond=m)
        ep.lanes.put(h, ep.spare_ids)
        m = c.hmr(q)
        c.cz(f1, f3, cond=m)
        ep.lanes.put(q, ep.spare_ids)
        return
    t = j + 1
    q = ep.lanes.get()
    ep.form_into(E.THR[t], q)
    if generic:
        h = ep.lanes.get()
        c.ccx(q, s, h)
        c.zp(h, ang)
        m = c.hmr(h)
        c.cz(q, s, cond=m)
        ep.lanes.put(h, ep.spare_ids)
    else:                                                  # j = 8: S on q * s as three T gates
        T8 = ONE // 8
        c.zp(q, T8)
        c.zp(s, T8)
        c.cx(q, s)
        c.zp(s, -T8)
        c.cx(q, s)
    ep.form_erase(E.THR[t], q)
    ep.lanes.put(q, ep.spare_ids)


EARLY_LOG = []                                             # FC sixth holder: (t_first, early pairs) per build


def emit_ring(c, word, flags, net, slot_ang, wrap_ang, cap, spares=(), const_token=None, const_one=False,
              retire=(), const_lane=False, retire_slots=False, late_w25=None, const_early=False):
    """const_one: the constant token's wire already holds 1 (FC's top word bit w_26 = 1 on every input) and is
    left holding 1; otherwise it must be |0> and the one is made and erased here.
    retire: list of (wire, rebuild): each wire is X-measured before the routing and lent to the walk as a spare
    lane (count it in cap); after the walk, in list order, rebuild(c, wire) must write its value back into the |0>
    wire from wires that are home again (e.g. a minimum-AND program), then the measurement is repaid (Z if the
    outcome was 1). The line does this with v9..v11; FC's q21 / q20 / q19 rebuild in 1 / 2 / 4 ANDs."""
    retire = list(retire)
    mbits = [c.hmr(wq) for (wq, _) in retire]
    spares = list(spares) + [wq for (wq, _) in retire]
    N, L, order, held = net["N"], net["L"], net["order"], set(net["held"])
    controls = net["controls"]
    stages = [G.pairs(cc, N) for cc in net["centres"]]
    win = G.windows(net["centres"], order, N)
    nfeat = 10
    pos_slot = {L[p]: p for p in range(N)}
    r = net.get("stride", 1)                               # token i sits at ring position r * i (r = +-1)
    ringw = [None] * N
    for i in range(N):
        ringw[(r * i) % N] = word[i]
    if const_token is not None:
        const_token = (r * const_token) % N               # from here on: the constant's ring position
    ep = E.Endpoint(c, ringw, flags, slot_ang, wrap_ang, cap, spares)
    w = ep.word
    t0 = c.toffolis()
    t_cl = None
    if const_token is not None and const_lane:
        # F1 fifth holder: after the constant token's last pair it sits at home, constant and untouched to the end of
        # the walk (FC's relwin measurement): from the next visit it is X'd to 0 and lent as one more lane; the held
        # pairs are re-packed with that extra capacity (pairs ending earlier keep the same choices).
        vis0 = {x: i for i, x in enumerate(order)}
        t_cl = 1 + max(win[key][1] for key in win if const_token in key[1])
        t_25 = None
        if late_w25 is not None:                           # word bit 25 lent after its last pair AND its home
            h25 = (r * 25) % N                             # slot's rotation; rebuilt after the walk (1 AND)
            t_25 = max(1 + max(win[key][1] for key in win if h25 in key[1]), vis0[h25] + 1)
            if t_25 >= N:
                t_25 = None
        if t_cl < N or t_25 is not None:
            nf_m = 9 if net.get("profile_line", True) else 10
            base = [0] * N
            for j in range(nf_m):
                base[vis0[L[j]]] += G.HELP.get(j, 0)
            fv = [vis0[L[j]] for j in G.PFEATS if j < nf_m]
            for t in range(min(fv), max(fv) + 1):
                base[t] += 1
            base = [b - (1 if t >= t_cl else 0) - (1 if t_25 is not None and t >= t_25 else 0)
                    for t, b in enumerate(base)]
            h2, held2 = G.pack(win, base, cap)
            if h2 >= len(held):
                held = set(held2)
        if t_cl >= N:
            t_cl = None
    else:
        t_25 = None
    if const_token is not None:
        if not const_one:
            c.x(w[const_token])
        cands =[key for key in win if const_token in key[1] and key not in held]
        last_one = max(cands, key=lambda kk: (win[kk][1], -kk[0])) if cands else None
    else:
        last_one = None
    vis = {x: i for i, x in enumerate(order)}
    early = set()                                          # FC sixth holder: the constant's wire as a lane before
    t_first = None                                         # its first pair (R45-INV-END's packing note)
    if const_token is not None and const_early and const_one:
        t_first = min(win[key][0] for key in win if const_token in key[1])
        busy = -1
        for key in sorted(win, key=lambda kk: (win[kk][1], win[kk][0])):
            a_, b_ = win[key]
            if key in held or b_ >= t_first or a_ <= busy:
                continue
            early.add(key)
            busy = b_
        if not early:
            t_first = None
    EARLY_LOG.append((t_first, sorted(early)))
    def reads_p(form):
        return "p" in " ".join([a + " " + b for a, b in form[0]] + [form[1]]).split()

    if const_token is not None:                            # the line's forms (ENC j, else THR j + 1)
        pform_feats = [j for j in range(9) if reads_p(E.ENC[j] if j in E.ENC else E.THR[j + 1])]
        feat_fn = feature
    else:
        pform_feats = [j for j in range(nfeat) if reads_p(TFORM[j + 1]) and (len(wrap_ang) == 10 or j != 9)]
        feat_fn = feature_noconst
    pv = [vis[L[j]] for j in pform_feats]
    p_on, p_off = min(pv), max(pv)

    def open_pair(s, pr):
        if (s, pr) in early:                               # FC: on the constant's wire (|0> here), alone on it
            sp0 = list(ep.lanes.spare)
            ep.lanes.spare = [w[const_token]]
            ep.lanes.cap += 1
            anc = ep.hold_ctl(controls[s], *pr)
            assert anc == w[const_token], ("early lane", anc)
            ep.lanes.spare = sp0
            return anc
        if (s, pr) in held:
            return ep.hold_ctl(controls[s], *pr)
        with ep.affine(controls[s]) as q:
            c.cswap(q, w[pr[0]], w[pr[1]])
        return None

    def close_pair(s, pr, anc):
        if (s, pr) in early:                               # FC: the constant's wire back to |0>, out of the pool
            ep.release_ctl(controls[s], *pr, anc)
            assert w[const_token] in ep.lanes.spare
            ep.lanes.spare.remove(w[const_token])
            ep.lanes.cap -= 1
            return
        if anc is not None:
            ep.release_ctl(controls[s], *pr, anc)
            return
        if (s, pr) == last_one:                            # one-token trick: the one returns home here for good
            home = const_token
            x = w[pr[0] if pr[1] == home else pr[1]]
            one = w[home]
            with ep.affine(controls[s]) as q:
                c.cx(one, x)
                c.x(one)
                m = c.hmr(one)
                c.cz(q, x, cond=m)
                c.x(one)
                c.cx(one, x)
            return
        with ep.affine(controls[s]) as q:
            c.cswap(q, w[pr[0]], w[pr[1]])

    opens, closes = {}, {}
    for key, (a, b) in win.items():
        opens.setdefault(a, []).append(key)
        closes.setdefault(b, []).append(key)
    live = {}
    pw = None
    zlog = []                                              # (op index, 'slot' / 'feat', index, angle) for LPs

    def zrec(kind, idx, n0):
        for oi in range(n0, len(c.ops)):
            if c.ops[oi][0] == "ZP":
                zlog.append((oi, kind, idx, c.ops[oi][2]))
    if t_first is not None:                                # FC: 1 -> 0, a lane for the early pairs only
        c.x(w[const_token])
        ep.spare_ids.add(w[const_token])
    for t in range(N):
        if t_first is not None and t == t_first:           # FC: back to the constant 1 before its first pair
            assert w[const_token] not in ep.lanes.spare
            ep.spare_ids.discard(w[const_token])
            c.x(w[const_token])
        if t_cl is not None and t == t_cl:                 # the constant's wire becomes a lane (holds 0 from here)
            c.x(w[const_token])
            ep.lanes.cap += 1
            ep.lanes.spare.append(w[const_token])
            ep.spare_ids.add(w[const_token])
        if t_25 is not None and t == t_25:                 # w25 is home for good and rotated: X-measure, lend
            m25 = c.hmr(w[h25])
            ep.lanes.cap += 1
            ep.lanes.spare.append(w[h25])
            ep.spare_ids.add(w[h25])
        if t == p_on:
            pw = ep.lanes.get()
            ep.form_into(E.P_FORM, pw)
            ep.var["p"] = pw
            ep.vlist[4] = pw
        for key in sorted(opens.get(t, []), key=lambda kk: kk[0]):
            live[key] = open_pair(*key)
        x = order[t]
        p = pos_slot[x]
        n0 = len(c.ops)
        c.zp(w[x], slot_ang[p])
        zrec("slot", p, n0)
        if p < nfeat and not (p == 9 and const_token is not None):
            n0 = len(c.ops)
            feat_fn(ep, p, w[x], wrap_ang)
            zrec("feat", p, n0)
        for key in sorted(closes.get(t, []), key=lambda kk: -kk[0]):
            close_pair(*key, live.pop(key))
        if t == p_off:
            ep.form_erase(E.P_FORM, pw)
            ep.var["p"] = None
            ep.vlist[4] = None
            ep.lanes.put(pw, ep.spare_ids)
    assert not live
    if t_cl is not None:                                   # take the constant's wire back and make it 1 again
        assert w[const_token] in ep.lanes.spare, "const lane still in use after the walk"
        ep.lanes.spare.remove(w[const_token])
        ep.lanes.cap -= 1
        ep.spare_ids.discard(w[const_token])
        c.x(w[const_token])
    if t_25 is not None:                                   # take w25's wire back (|0>); rebuilt after the retires
        assert w[h25] in ep.lanes.spare, "w25 lane still in use after the walk"
        ep.lanes.spare.remove(w[h25])
        ep.lanes.cap -= 1
        ep.spare_ids.discard(w[h25])
    if const_token is not None:
        if len(wrap_ang) == 10:                           # the one sits at slot 9 only in class 10: its feature
            q = ep.lanes.get()
            c.ccx(ep.var["f1"], ep.var["f3"], q)
            zlog.append((len(c.ops), "feat", 9, wrap_ang[9] % ONE))
            c.zp(q, wrap_ang[9])
            m = c.hmr(q)
            c.cz(ep.var["f1"], ep.var["f3"], cond=m)
            ep.lanes.put(q, ep.spare_ids)
        else:
            c.cz(ep.var["f1"], ep.var["f3"])
        if not const_one:
            c.x(w[const_token])
    for idx, ((wq, rebuild), m) in enumerate(zip(retire, mbits)):
        others = [r_[0] for r_ in retire[idx + 1:]] if retire_slots else []
        if others:                                         # F1: the |0> slots of bits rebuilt later serve as temps
            c.free_(others)
        rebuild(c, wq)
        for o in others:
            c.free.remove(o)
            c.live += 1
        c.z(wq, cond=m)
    if t_25 is not None:                                   # w25 = 1 + q21 + q20 w24 (FC's w25ret relation), then
        q20w, q21w = late_w25                              # the X-measurement's kick repaid
        c.ccx(q20w, w[(r * 24) % N], w[h25])
        c.cx(q21w, w[h25])
        c.x(w[h25])
        c.z(w[h25], cond=m25)
    return dict(toffolis=c.toffolis() - t0, peak_lanes=ep.lanes.peak, zps=zlog)


BEST_CENTRES = {26: [(2, 6, 24, 0, 25)],                              # ringN winners, p-free controls
                27: [(2, 6, 2, 0, 26), (2, 6, 25, 0, 26), (2, 0, 4, 0, 26)]}



def _r55_net(path, cap):
    """R55-INV-ROUTER (own copy of the src only): a GENERAL layered swap network with one linear control per pair,
    read from json {layers: [[(u, v, mask)]], order: [slots]} in slot coordinates (class 0 fires nothing: L = identity).
    Every pair is its own "stage" (numbered in layer order); G.windows is replaced by the closure cones."""
    import json as _json
    d = _json.load(open(os.path.expanduser(path)))
    layers = [[tuple(e) for e in ll] for ll in d["layers"] if ll]
    order = list(d["order"])
    N = 26
    keys, controls = [], []
    s = 0
    for l, lay in enumerate(layers):
        for (u, v, f) in lay:
            cone = {u, v}
            for l2 in range(l + 1, len(layers)):
                add = set()
                for (a, b, _) in layers[l2]:
                    if a in cone or b in cone:
                        add |= {a, b}
                cone |= add
            keys.append((s, (min(u, v), max(u, v)), frozenset(cone)))
            controls.append(tuple(1 << b for b in range(4) if f >> b & 1))
            s += 1

    def windows(_centres, order_, _n):
        vis = {x: i for i, x in enumerate(order_)}
        out = {}
        for (s_, pr, cone) in keys:
            ts = [vis[x] for x in cone]
            out[(s_, pr)] = (min(ts), max(ts))
        return out

    G.windows = windows
    G.pairs = lambda c, n: []
    win = windows(None, order, N)
    vis = {x: i for i, x in enumerate(order)}
    L = list(range(N))
    base = [0] * N
    for j in range(9):
        base[vis[L[j]]] += G.HELP.get(j, 0)
    fv = [vis[L[j]] for j in G.PFEATS if j < 9]
    for t in range(min(fv), max(fv) + 1):
        base[t] += 1
    h, held = G.pack(win, base, cap)
    print("R55 net: pairs %d, cap %d, model cost %d, held %d" % (len(win), cap, 2 * len(win) - h, h), flush=True)
    return dict(N=N, centres=tuple(range(len(keys))), controls=controls, L=L, order=order, held=held,
                model_cost=2 * len(win) - h, stride=1)


def best_net(N, cap, rank=0, profile="generic", centres=None):
    G.use_profile(profile)
    if os.environ.get("R55_NET") and N == 26:
        return _r55_net(os.environ["R55_NET"], cap)
    nfeat = 10 if profile != "line" else 9
    cl = [centres] if centres is not None else BEST_CENTRES.get(N)
    if cl:
        res = []
        for cen in cl:
            order = G.walk_order(cen[3], cen[4], N)
            res += [(r[0], cen, order, r) for r in G.evaluate(cen, N, cap, order, nfeat)]
        res.sort(key=lambda t: t[0])
    else:
        res, _ = G.search(N, cap, top=rank + 1, nfeat=nfeat)
    if len(res) <= rank:
        raise ValueError(f"no endpoint network fits {cap} free wires on ring {N} (profile {profile}): the wrap "
                         f"features alone need 2 helpers + p; give the endpoint at least 3 lanes")
    cost, centres, order, r = res[rank]
    _, npairs, h, kind, rr, d0, ctl, held, L = r
    controls = []
    for s in range(5):
        mask = G.AFF[ctl[s]]
        terms = tuple(([0] if mask & 1 else []) + [1 << b for b in range(5) if mask >> (b + 1) & 1])
        controls.append(terms)
    return dict(N=N, centres=centres, controls=controls, L=L, order=order, held=held, model_cost=cost, stride=rr)


def selftest(N, cap, mode="generic", const=None, shots=4096, seed=5, net=None, const_one=False):
    import random
    import numpy as np
    profile = "line" if const is not None else ("generic" if mode == "generic" else "clifford")
    net = net or best_net(N, cap, profile=profile)
    rnd = random.Random(seed)
    c = K.Circ(n_fixed=N + 4)
    word, flags = list(range(N)), list(range(N, N + 4))
    spares = c.alloc(min(cap, 3))
    slot_ang = [rnd.randrange(ONE) for _ in range(N)]
    wrap_ang = [rnd.randrange(ONE) for _ in range(10 if mode == "generic" else 7)]
    info = emit_ring(c, word, flags, net, slot_ang, wrap_ang, cap, spares=spares, const_token=const,
                     const_one=const_one)
    S = K.Sim(shots, max(c.next, N + 8), seed=seed)
    rng = np.random.default_rng(seed)
    ks = rng.integers(0, 11, shots)
    toks = rng.integers(0, 2, (N, shots)).astype(bool)
    if const is not None:
        toks[const] = const_one
    for i in range(N):
        S.q[word[i]] = toks[i]
    for b in range(4):
        S.q[flags[b]] = ((ks >> b) & 1).astype(bool)
    init = S.q.copy()
    S.run(c.ops)
    F = wrap_ang if len(wrap_ang) == 10 else wrap_ang + [ONE // 8, ONE // 4, ONE // 2]
    bad = 0
    for n in range(shots):
        k = int(ks[n])
        tot = 0
        for i in range(N):
            t = 1 if i == const else int(toks[i][n])
            if t:
                tot += slot_ang[(i + k) % N] + (F[i + k - N] if i + k >= N else 0)
        got = (int(S.ang[n]) + (ONE // 2 if S.sign[n] else 0)) % ONE
        bad += got != tot % ONE
    changed = int((S.q != init).any(axis=1).sum())
    print(f"ring {N} cap {cap} {mode} const {const}: model {net['model_cost']}, held {len(net['held'])}, "
          f"Toffolis {info['toffolis']}, peak lanes {info['peak_lanes']}; phase mismatches {bad}, wires changed {changed}")
    return info, bad, changed


def selftest_retire(N=27, cap=3, nret=1, const=26, shots=4096, seed=9):
    """nret extra wires, each the AND of two data tokens, retired across the walk and rebuilt with one Toffoli each;
    the endpoint gets cap + nret lanes. Checks the phase on every shot and that every wire comes back."""
    import random
    import numpy as np
    net = best_net(N, cap + nret, profile="line")
    rnd = random.Random(seed)
    c = K.Circ(n_fixed=N + 4 + nret)
    word, flags = list(range(N)), list(range(N, N + 4))
    rw = list(range(N + 4, N + 4 + nret))
    pairs_ = [(1 + 2 * i, 2 + 2 * i) for i in range(nret)]
    slot_ang = [rnd.randrange(ONE) for _ in range(N)]
    wrap_ang = [rnd.randrange(ONE) for _ in range(7)]
    retire = [(rw[i], (lambda cc, q, a=pairs_[i][0], b=pairs_[i][1]: cc.ccx(word[a], word[b], q)))
              for i in range(nret)]
    info = emit_ring(c, word, flags, net, slot_ang, wrap_ang, cap + nret, spares=(), const_token=const,
                     const_one=True, retire=retire)
    S = K.Sim(shots, max(c.next, N + 8 + nret), seed=seed)
    rng = np.random.default_rng(seed)
    ks = rng.integers(0, 11, shots)
    toks = rng.integers(0, 2, (N, shots)).astype(bool)
    toks[const] = True
    for i in range(N):
        S.q[word[i]] = toks[i]
    for b in range(4):
        S.q[flags[b]] = ((ks >> b) & 1).astype(bool)
    for i in range(nret):
        a, b = pairs_[i]
        S.q[rw[i]] = toks[a] & toks[b]
    init = S.q.copy()
    S.run(c.ops)
    F = wrap_ang + [ONE // 8, ONE // 4, ONE // 2]
    bad = 0
    for n in range(shots):
        k = int(ks[n])
        tot = 0
        for i in range(N):
            if toks[i][n]:
                tot += slot_ang[(i + k) % N] + (F[i + k - N] if i + k >= N else 0)
        got = (int(S.ang[n]) + (ONE // 2 if S.sign[n] else 0)) % ONE
        bad += got != tot % ONE
    changed = int((S.q != init).any(axis=1).sum())
    print(f"retire test ring {N} cap {cap} + {nret} retired: Toffolis {info['toffolis']} (incl. {nret} rebuild ANDs), "
          f"peak lanes {info['peak_lanes']}; phase mismatches {bad}, wires changed {changed}")


if __name__ == "__main__":
    if sys.argv[1] == "retire":
        selftest_retire(cap=int(sys.argv[2]), nret=int(sys.argv[3]))
        sys.exit(0)
    N, cap = int(sys.argv[1]), int(sys.argv[2])
    mode = sys.argv[3] if len(sys.argv) > 3 else "generic"
    const = int(sys.argv[4]) if len(sys.argv) > 4 else None
    if len(sys.argv) > 5 and sys.argv[5] == "one":
        selftest(N, cap, mode, const, const_one=True)
        sys.exit(0)
    selftest(N, cap, mode, const)
