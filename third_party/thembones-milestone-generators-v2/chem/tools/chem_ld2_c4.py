"""chem_ld2_c4: the 68-AND high cache (item C-4) for our own loader compiler (builder C2, seventh holder).

Mechanism (read from the public 4,038 loaders, research/ws_chem_ld2.md, seventh holder): on loader k >= 1 the
delta words D_k = W_{k-1} XOR W_k admit a block set V (77 high blocks; |V| even, XOR of indices 0, block 76 in V)
whose V-sum of every payload bit is the same on every low row l (rows 1220..1231 = block 76, l >= 4, are free and
chosen to match).  Then every feature's high part lies in H = {f : sum_V f = 0} (76 dims incl. 1 and the 7 address
bits), so the high cache needs 68 ANDs instead of 69, provided the plan's tau (XOR of sigma = sum_V X Y over the
products feeding a wire) equals the table's kappa (V-sum at l = 0) on every wire.

This module:
  vset(words, k)             -> V as a 77-bit int (None if the table has no such V for loader k)
  pad_table(words, Vs)       -> words extended to 1,232 rows: padding rows hold words whose deltas are V-consistent
                                for every k with a V (bit 7 = bit 6 XOR the slot's carrier function, as in the table)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_ld2_c4pub as cp  # noqa: E402
import chem_loaders_core as core  # noqa: E402
import chem_ld2_hb as hbm  # noqa: E402  (tenth holder: block count from CHEM_HB, default 77)

J = 5
NFULL = hbm.NFULL
_LAST = hbm.LAST


def vset(words, k):
    return cp.vset(words, k)


def carrier_funcs(words):
    """Per slot j: the pure low function g_j(l) with bit7 ^ bit6 == g_j(l mod 16) on every row of every batch, or
    None if that slot has no such carrier (then bit 7 is treated as an ordinary payload bit)."""
    A = words.shape[0]
    l = np.arange(A) & 15
    out = []
    tb = core.BROT - 1                                 # tenth holder: top bit of a word (7 at 8 bits, 6 at 7 bits)
    for j in range(J):
        x = ((words[:, j::J] >> tb) ^ (words[:, j::J] >> (tb - 1))) & 1      # rows x batches
        g = np.zeros(16, np.int64)
        ok = True
        for ll in range(16):
            v = x[l == ll]
            if v.size and not (np.all(v == v.flat[0])):
                ok = False
                break
            g[ll] = v.flat[0] if v.size else 0
        out.append(g if ok else None)
    return out


def pad_table(words, Vs):
    """words (A x 75) -> (1232 x 75). Vs: dict k -> V (77-bit int) for the loaders that use the 68-AND cache.
    Batch 0's padding words hold 0 in bits 0-6; batch k's padding = batch k-1's XOR the V-consistent delta."""
    A, ncol = words.shape
    K = ncol // J
    assert A <= NFULL and A > NFULL - 16
    g = carrier_funcs(words)
    out = np.zeros((NFULL, ncol), np.int64)
    out[:A] = words
    pl = np.arange(A, NFULL) & 15                      # low row of each padding row (block 76)
    # bits 0-6 (and 7 where a slot has no carrier) of batch 0's padding: 0, or V_0-consistent when loader 0 has C-4
    if 0 in Vs and Vs[0] is not None:
        V0 = Vs[0]
        Vp0 = [p for p in range(_LAST) if (V0 >> p) & 1]
        cur0 = out[:, 0:J]
        for j in range(J):
            nb = core.BROT - 1 if g[j] is not None else core.BROT     # payload bits (7 / 8 at 8-bit words)
            for b in range(nb):
                grid = np.zeros((hbm.HB, 16), np.int64)
                grid.reshape(-1)[:A] = (words[:, j] >> b) & 1
                kappa = int(np.bitwise_xor.reduce(grid[Vp0 + [_LAST], 0])) & 1
                s = np.bitwise_xor.reduce(grid[Vp0, :], axis=0) & 1
                want76 = (kappa ^ s) & 1
                for r, ll in zip(range(A, NFULL), pl):
                    cur0[r, j] = (cur0[r, j] & ~(1 << b)) | (int(want76[ll]) << b)
    for k in range(1, K):
        cur = out[:, J * k:J * k + J]
        prev = out[:, J * k - J:J * k]
        if k not in Vs or Vs[k] is None:
            cur[A:] = prev[A:]                         # delta 0 on the padding rows (any choice is valid)
            continue
        V = Vs[k]
        Vp = [p for p in range(_LAST) if (V >> p) & 1]     # V without the last block
        d = (prev[:A] ^ cur[:A])
        for j in range(J):
            nb = core.BROT - 1 if g[j] is not None else core.BROT
            for b in range(nb):
                db = (d[:, j] >> b) & 1
                grid = np.zeros((hbm.HB, 16), np.int64)
                grid.reshape(-1)[:A] = db
                kappa = int(np.bitwise_xor.reduce(grid[Vp + [_LAST], 0])) & 1
                s = np.bitwise_xor.reduce(grid[Vp, :], axis=0) & 1      # V-sum without block 76, per l
                want76 = (kappa ^ s) & 1                                  # block 76's delta bit per l
                for r, ll in zip(range(A, NFULL), pl):
                    if want76[ll]:
                        cur[r, j] = (cur[r, j] & ~(1 << b)) | ((prev[r, j] ^ (1 << b)) & (1 << b))
                    else:
                        cur[r, j] = (cur[r, j] & ~(1 << b)) | (prev[r, j] & (1 << b))
    # carrier bit 7 on every padding word
    for j in range(J):
        if g[j] is None:
            continue
        for k in range(K):
            c = out[A:, J * k + j]
            tb = core.BROT - 1                         # top bit = (bit below it) xor the slot's carrier
            b6 = (c >> (tb - 1)) & 1
            out[A:, J * k + j] = (c & ((1 << tb) - 1)) | (((b6 ^ g[j][pl]) & 1) << tb)
    return out


def check_padded(words_p, Vs):
    """For every k in Vs: the V-sum of every bit of D_k is constant over all 16 low rows (incl. padding)."""
    bad = {}
    for k, V in Vs.items():
        if V is None:
            continue
        Vl = [p for p in range(hbm.HB) if (V >> p) & 1]
        d = words_p[:, 0:J].copy() if k == 0 else words_p[:, J * k - J:J * k] ^ words_p[:, J * k:J * k + J]
        nbad = 0
        for j in range(J):
            for b in range(core.BROT):
                grid = ((d[:, j] >> b) & 1).reshape(hbm.HB, 16)
                s = np.bitwise_xor.reduce(grid[Vl, :], axis=0) & 1
                nbad += int(np.any(s != s[0]))
        bad[k] = nbad
    return bad


# ------------------------------------------------------------------------------------------------ high functions
HB = hbm.HB
_HT =[sum(1 << p for p in range(HB) if (p & M) == M) for M in range(HB)]     # truth of x^M over the 77 blocks


def hmask_truth(mask):
    """77-bit truth table of the high function with monomial mask (bit M-1 = monomial M), core's convention."""
    t = 0
    while mask:
        b = mask & -mask
        t ^= _HT[b.bit_length()]
        mask ^= b
    return t


def bform(x, y, V):
    return bin(x & y & V).count("1") & 1


def kappa_wires(an, V):
    """Bit q = V-sum of wire q's delta at low row 0 (rows 16p)."""
    Vl = [p for p in range(HB) if (V >> p) & 1]
    out = 0
    for q, d in enumerate(an.delta_int):
        if sum((d >> (16 * p)) & 1 for p in Vl) & 1:
            out |= 1 << q
    return out


def gmap(an, s):
    """Wire-space vector of a basis-output mask s (wire q is fed iff G[q] & s is odd)."""
    out = 0
    for q, g in enumerate(an.G):
        if (g & s).bit_count() & 1:
            out |= 1 << q
    return out


# ------------------------------------------------------------------------------------------------ sector closure
def normalize(prods):
    """Every product rewritten with the canonical basis (l1 < l2) of its low plane: (l1, B1, l2, B2, s, phase) with
    features L(l1) H(B1) + L(l2) H(B2), i.e. the product (L(l1) + H(B2)) (L(l2) + H(B1))."""
    out = []
    for (a1, b1, a2, b2, s, ph) in prods:
        if not a1 or not a2 or a1 == a2:               # a single (rank-1 features): kept as it is, in no sector
            out.append([a1, b1, a2, b2, s, ph, "single"])
            continue
        els = sorted({a1, a2, a1 ^ a2})
        l1, l2 = els[0], els[1]

        def co(x):
            return (1, 0) if x == l1 else (0, 1) if x == l2 else (1, 1)
        al, ga = co(a1), co(a2)
        B1 = (b1 if al[0] else 0) ^ (b2 if ga[0] else 0)
        B2 = (b1 if al[1] else 0) ^ (b2 if ga[1] else 0)
        out.append([l1, B1, l2, B2, s, ph, "pair"])
    return out


def tau_state(an, P, V):
    t = 0
    for p in P:
        l1, B1, l2, B2, s = p[:5]
        if bform(hmask_truth(B2), hmask_truth(B1), V):
            t ^= s
    return gmap(an, t)


def solve_gf2(vecs, target):
    """Mask c with the XOR of vecs[i] over the set bits of c equal to target, or None."""
    piv = {}
    for i, v in enumerate(vecs):
        tag = 1 << i
        while v:
            hb = v.bit_length() - 1
            if hb not in piv:
                piv[hb] = (v, tag)
                break
            pv, pt = piv[hb]
            v ^= pv
            tag ^= pt
    t, tag = target, 0
    while t:
        hb = t.bit_length() - 1
        if hb not in piv:
            return None
        pv, pt = piv[hb]
        t ^= pv
        tag ^= pt
    return tag


def _move(P, x, y):
    """Elementary move inside a sector: s_y += s_x, X_x += X_y, Y_x += Y_y (features unchanged)."""
    P[y][4] ^= P[x][4]
    P[x][3] ^= P[y][3]
    P[x][1] ^= P[y][1]


def closure(an, prods, V, seed=0, trials=400):
    """Sector closure: elementary moves inside sectors (same low plane and phase) until tau == kappa on every wire.
    Returns (new products in core's tuple form, info)."""
    import random
    rng = random.Random(seed)
    P = normalize(prods)
    kap = kappa_wires(an, V)
    sectors = {}
    for i, p in enumerate(P):
        if p[6] == "pair":
            sectors.setdefault((p[0], p[2], p[5]), []).append(i)
    multi = [v for v in sectors.values() if len(v) >= 2]
    mu0 = tau_state(an, P, V) ^ kap
    info = dict(sectors=len(sectors), multi=len(multi), mismatch0=bin(mu0).count("1"), moves=0, trials=0,
                multi_products=sum(len(v) for v in multi))
    X = [hmask_truth(p[3]) for p in P]
    Y = [hmask_truth(p[1]) for p in P]
    gv = [gmap(an, p[4]) for p in P]
    mu = mu0
    w = bin(mu).count("1")
    steps = trials * 25
    for t in range(steps):
        if mu == 0:
            break
        info["trials"] = t + 1
        # every move (x -> y) with omega_xy = 1 toggles tau by product x's wire vector; pick the best one
        best, bw = [], w
        for idx in multi:
            for x in idx:
                nw_ = bin(mu ^ gv[x]).count("1")
                if nw_ > bw:
                    continue
                ys = [y for y in idx if y != x and P[x][4] != P[y][4] and (bform(X[x], Y[y], V) ^ bform(X[y], Y[x], V))]
                if not ys:
                    continue
                if nw_ < bw:
                    best, bw = [(x, y) for y in ys], nw_
                elif nw_ == bw and bw < w:
                    best += [(x, y) for y in ys]
        if best:
            x, y = rng.choice(best)
            mu ^= gv[x]
        else:                                            # no improving move: a random move (it may be neutral)
            idx = rng.choice(multi)
            x, y = rng.sample(idx, 2)
            if P[x][4] == P[y][4]:
                continue
            if bform(X[x], Y[y], V) ^ bform(X[y], Y[x], V):
                mu ^= gv[x]
        _move(P, x, y)
        X[x] ^= X[y]
        Y[x] ^= Y[y]
        gv[y] ^= gv[x]
        w = bin(mu).count("1")
        info["moves"] += 1
    mu_inc = mu
    mu = tau_state(an, P, V) ^ kap
    info["formula_ok"] = bool(mu == mu_inc)              # the predicted change omega_xy s_x, checked exactly
    info["mismatch"] = bin(mu).count("1")
    info["factors_outside_H"] = sum(bform(hmask_truth(p[1]), (1 << HB) - 1, V) + bform(hmask_truth(p[3]), (1 << HB) - 1, V)
                                    for p in P)
    info["singles"] = sum(1 for p in P if p[6] == "single")
    return [tuple(p[:6]) for p in P], info


def _rank(vecs):
    piv = {}
    for v in vecs:
        while v:
            hb = v.bit_length() - 1
            if hb not in piv:
                piv[hb] = v
                break
            v ^= piv[hb]
    return len(piv)


def closure2(an, prods, V, seed=0, walk=600, trials=3000):
    """Sector closure by reachable sets: for every sector (same low plane and phase, >= 2 products) a random walk of
    elementary moves records each reachable change of tau (wire space) with the move sequence that reaches it from
    the start; then one change per sector is chosen by randomized GF(2) elimination so that the changes sum to the
    mismatch. Returns (products, info); info['mutable_rank'] is the rank of all reachable changes."""
    import random
    rng = random.Random(seed)
    P = normalize(prods)
    kap = kappa_wires(an, V)
    sectors = {}
    for i, p in enumerate(P):
        if p[6] == "pair":
            sectors.setdefault((p[0], p[2], p[5]), []).append(i)
    multi = [v for v in sectors.values() if len(v) >= 2]
    mu = tau_state(an, P, V) ^ kap
    info = dict(sectors=len(sectors), multi=len(multi), multi_products=sum(len(v) for v in multi),
                mismatch0=bin(mu).count("1"))
    reach = []                                          # per sector: {delta: move sequence}
    for idx in multi:
        loc = {i: [P[i][1], P[i][3], P[i][4]] for i in idx}      # B1, B2, s of the sector's products
        X = {i: hmask_truth(loc[i][1]) for i in idx}
        Y = {i: hmask_truth(loc[i][0]) for i in idx}
        gv = {i: gmap(an, loc[i][2]) for i in idx}
        seen = {0: []}
        seq, delta = [], 0
        for t in range(walk):
            x, y = rng.sample(idx, 2)
            if gv[x] == gv[y] or loc[x][2] == loc[y][2]:
                continue
            if bform(X[x], Y[y], V) ^ bform(X[y], Y[x], V):
                delta ^= gv[x]
            loc[y][2] ^= loc[x][2]
            loc[x][1] ^= loc[y][1]
            loc[x][0] ^= loc[y][0]
            X[x] ^= X[y]
            Y[x] ^= Y[y]
            gv[y] ^= gv[x]
            seq.append((x, y))
            if delta not in seen:
                seen[delta] = list(seq)
        reach.append(seen)
    allv = [d for seen in reach for d in seen if d]
    info["mutable_rank"] = _rank(allv)
    info["reach_sizes"] = sorted((len(s) for s in reach), reverse=True)[:12]
    sol = None
    if mu == 0:
        sol = [0] * len(multi)
    for t in range(trials if mu else 0):
        picks = [rng.choice([d for d in seen if d]) if len(seen) > 1 else 0 for seen in reach]
        c = solve_gf2(picks, mu)
        if c is not None:
            sol = [picks[i] if (c >> i) & 1 else 0 for i in range(len(multi))]
            info["trials"] = t + 1
            break
    if sol is None:
        info["mismatch"] = bin(mu).count("1")
        info["solved"] = False
        return [tuple(p[:6]) for p in P], info
    for seen, d in zip(reach, sol):
        for (x, y) in seen[d]:
            _move(P, x, y)
    mu = tau_state(an, P, V) ^ kap
    info["mismatch"] = bin(mu).count("1")
    info["solved"] = mu == 0
    info["factors_outside_H"] = sum(bform(hmask_truth(p[1]), (1 << HB) - 1, V) + bform(hmask_truth(p[3]), (1 << HB) - 1, V)
                                    for p in P)
    return [tuple(p[:6]) for p in P], info


def wire_weights(an, base=4.0):
    """Weight of a mismatched wire = price scale of its kappa toggle: base ** bit (bit of the wire's (slot, bit))."""
    return [base ** b for (j, b) in an.out_wires]


def closure3(an, prods, V, seed=0, walk=800, trials=2000, stages=60, mask=None):
    """Multi-stage sector closure. Each stage: reachable tau changes per sector from the CURRENT state (random walks
    of elementary moves), then the one-change-per-sector GF(2) solve (randomized) on the wires of `mask` (all if
    None); if it fails, the single sector change that lowers the weighted mismatch most is applied and the next
    stage starts from there. Wires outside `mask` are left for kappa toggles (cheap bits)."""
    import random
    rng = random.Random(seed)
    nwire = len(an.out_wires)
    if mask is None:
        mask = (1 << nwire) - 1
    wts = wire_weights(an)

    def wcost(v):
        return sum(wts[q] for q in range(nwire) if (v >> q) & 1)
    P = normalize(prods)
    kap = kappa_wires(an, V)
    sectors = {}
    for i, p in enumerate(P):
        if p[6] == "pair":
            sectors.setdefault((p[0], p[2], p[5]), []).append(i)
    multi = [v for v in sectors.values() if len(v) >= 2]
    mu = tau_state(an, P, V) ^ kap
    info = dict(sectors=len(sectors), multi=len(multi), multi_products=sum(len(v) for v in multi),
                mismatch0=bin(mu).count("1"), stages=0, solved=False)

    def reach_of(idx):
        loc = {i: [P[i][1], P[i][3], P[i][4]] for i in idx}
        X = {i: hmask_truth(loc[i][1]) for i in idx}
        Y = {i: hmask_truth(loc[i][0]) for i in idx}
        gv = {i: gmap(an, loc[i][2]) for i in idx}
        seen = {0: []}
        seq, delta = [], 0
        for t in range(walk):
            x, y = rng.sample(idx, 2)
            if loc[x][2] == loc[y][2]:
                continue
            if bform(X[x], Y[y], V) ^ bform(X[y], Y[x], V):
                delta ^= gv[x]
            loc[y][2] ^= loc[x][2]
            loc[x][1] ^= loc[y][1]
            loc[x][0] ^= loc[y][0]
            X[x] ^= X[y]
            Y[x] ^= Y[y]
            gv[y] ^= gv[x]
            seq.append((x, y))
            if delta not in seen:
                seen[delta] = list(seq)
        return seen

    for st in range(stages):
        if mu & mask == 0:
            break
        info["stages"] = st + 1
        reach = [reach_of(idx) for idx in multi]
        if st == 0:
            info["mutable_rank"] = _rank([d & mask for seen in reach for d in seen if d & mask])
        sol = None
        for t in range(trials):
            picks = [rng.choice([d for d in seen if d]) if len(seen) > 1 else 0 for seen in reach]
            c = solve_gf2([d & mask for d in picks], mu & mask)
            if c is not None:
                sol = [picks[i] if (c >> i) & 1 else 0 for i in range(len(multi))]
                break
        if sol is not None:
            for seen, d in zip(reach, sol):
                for (x, y) in seen[d]:
                    _move(P, x, y)
            mu = tau_state(an, P, V) ^ kap
            break
        best = None                                      # greedy: the sector change with the lightest weighted residual
        cur = wcost(mu)
        for si, seen in enumerate(reach):
            for d in seen:
                if d:
                    w = wcost(mu ^ d)
                    if w < cur and (best is None or w < best[0]):
                        best = (w, si, d)
        if best is None:
            break
        _, si, d = best
        for (x, y) in reach[si][d]:
            _move(P, x, y)
        mu = tau_state(an, P, V) ^ kap
    info["mismatch"] = bin(mu).count("1")
    info["mismatch_masked"] = bin(mu & mask).count("1")
    info["residual"] = mu
    info["residual_cost_units"] = wcost(mu)
    info["solved"] = mu == 0
    info["factors_outside_H"] = sum(bform(hmask_truth(p[1]), (1 << HB) - 1, V) + bform(hmask_truth(p[3]), (1 << HB) - 1, V)
                                    for p in P)
    return [tuple(p[:6]) for p in P], info


def h_basis_masks(V):
    """Monomial masks (bit M-1 = monomial M) of a basis of the non-constant part of H = {f : sum_V f = 0}."""
    sig = {M: bin(_HT[M] & V).count("1") & 1 for M in range(1, HB)}
    star = next(M for M in range(1, HB) if sig[M])
    out = []
    for M in range(1, HB):
        if M == star:
            continue
        m = 1 << (M - 1)
        if sig[M]:
            m |= 1 << (star - 1)
        out.append(m)
    return out


DEBUG = []


def deficient_products(E, seed=0, rounds=4, tries=4):
    """For an engine state with deficiency 1 (one single in the parity basis): the products with the single completed
    into a pair by a candidate line whose second feature depends on the basis, and the relation: the list of
    (product index, tuple position 1 or 3) whose high coefficient can take + g for any g in H without changing a
    feature. Returns (prods, rel) or None when no completing line exists."""
    import importlib
    eng = importlib.import_module(type(E).__mro__[-2].__module__)
    E._sync()
    r = E.r
    an, binfo, points = E.an, E.binfo, E.points
    for attempt in range(tries):
        xr = np.random.default_rng(1009 * attempt + 17 * an.k + 5 + 7919 * seed)
        n = eng._lines(E.pt_blk, E.pt_ymask, E.pt_outs, E.pt_phase, E.ncoord, E.Rm, E.ecount, E.nb, E.li, E.lj,
                       E.ls, E.deg)
        lx = xr.integers(1, 1 << eng.GF_BITS, eng.MAXLINES).astype(np.uint16)
        Y = eng._build_y(n, E.li, E.lj, E.ls, lx, E.pt_blk, E.pt_ymask, E.ncoord, E.cellid, E.freeidx, E.nb, r, r)
        nu = eng.gf_rank(Y.copy(), eng.GF_EXP, eng.GF_LOG) // 2
        d = r - 2 * nu
        if d != 1:
            DEBUG.append(f"attempt {attempt}: deficiency {d} (r {r}, nu {nu})")
            return None
        ntot = eng._singles(E.pt_blk, E.pt_ymask, E.pt_outs, E.ncoord, E.Rm, E.ecount, E.nb, d, n, E.li, E.lj, E.ls)
        Y = eng._build_y(ntot, E.li, E.lj, E.ls, lx, E.pt_blk, E.pt_ymask, E.ncoord, E.cellid, E.freeidx, E.nb, r, r + d)
        ok, N = eng.gf_inverse(Y, eng.GF_EXP, eng.GF_LOG)
        if not ok:
            continue
        order = xr.permutation(ntot).astype(np.int64)
        keep = eng._extract(N, ntot, E.li, E.lj, E.ls, lx, order, E.pt_blk, E.pt_ymask, E.ncoord, E.cellid,
                            E.freeidx, E.nb, r, eng.GF_EXP, eng.GF_LOG)
        kept = np.flatnonzero(keep)
        pairs = [(int(E.li[q]), int(E.lj[q]), int(E.ls[q])) for q in kept if E.lj[q] >= 0]
        singles = [(int(E.li[q]), int(E.ls[q])) for q in kept if E.lj[q] < 0]
        if len(singles) != 1 or 2 * len(pairs) + 1 != r:
            DEBUG.append(f"attempt {attempt}: kept {len(pairs)} pairs + {len(singles)} singles for r {r}")
            continue
        si, ss = singles[0]
        partners = [int(E.lj[q]) if int(E.li[q]) == si else int(E.li[q]) for q in range(n)
                    if E.lj[q] >= 0 and int(E.ls[q]) == ss and si in (int(E.li[q]), int(E.lj[q]))]
        if not partners:
            DEBUG.append(f"attempt {attempt}: single ({si}, {ss}) has no candidate line ({n} lines)")
            continue
        feats = []
        for m, (i, j, s) in enumerate(pairs):
            feats.append((E._flat(i, s), ("p", m, 0)))
            feats.append((E._flat(j, s), ("p", m, 1)))
        feats.append((E._flat(si, ss), ("s", 0, 0)))
        fs = core.Span()
        for idx, (v, _) in enumerate(feats):
            assert fs.add(v, 1 << idx), "basis features dependent"
        coef = [0] * len(feats)
        for bi in binfo:
            for M, v in enumerate(bi["rows"]):
                rr, c = fs.reduce(v << bi["offset"])
                if rr:
                    break
                for idx in core.bits(c):
                    coef[idx] |= 1 << M
            else:
                continue
            break
        else:
            sj = partners[0]
            rr, tag = fs.reduce(E._flat(sj, ss))
            if rr:
                continue
            A = {key: coef[idx] for idx, (_, key) in enumerate(feats)}

            def phase_of(*pts):
                common = None
                for p in pts:
                    ph = set(binfo[points[p][0]]["phases"])
                    common = ph if common is None else common & ph
                return min(common) if common else None
            prods, where = [], {}
            for m, (i, j, s) in enumerate(pairs):
                where[("p", m, 0)] = (len(prods), 1)
                where[("p", m, 1)] = (len(prods), 3)
                prods.append([points[i][4], A[("p", m, 0)], points[j][4], A[("p", m, 1)], s, phase_of(i, j)])
            ph = phase_of(si, sj)
            if ph is None:
                continue
            where[("s", 0, 0)] = (len(prods), 1)
            prods.append([points[si][4], A[("s", 0, 0)], points[sj][4], 0, ss, ph])
            rel = [where[feats[idx][1]] for idx in core.bits(tag)] + [(len(prods) - 1, 3)]
            return prods, rel
    return None


def relation_closure(an, prods, rel, V):
    """tau == kappa by g in H: every relation line's high coefficient takes + g (features unchanged); tau(g) is
    linear in g [P: B is alternating on H]. Returns (products, info)."""
    kap = kappa_wires(an, V)
    basis = h_basis_masks(V)

    def with_g(gmask):
        P = [list(p) for p in prods]
        for (pi, pos) in rel:
            P[pi][pos] ^= gmask
        return P
    t0 = tau_state(an, [p + ["pair"] for p in with_g(0)], V)
    cols = [tau_state(an, [p + ["pair"] for p in with_g(g)], V) ^ t0 for g in basis]
    mu = t0 ^ kap
    c = solve_gf2(cols, mu)
    info = dict(mismatch0=bin(mu).count("1"), rank=_rank(cols), solved=c is not None)
    if c is None:
        return [tuple(p) for p in with_g(0)], info
    g = 0
    for i in range(len(basis)):
        if (c >> i) & 1:
            g ^= basis[i]
    P = with_g(g)
    info["mismatch"] = bin(tau_state(an, [p + ["pair"] for p in P], V) ^ kap).count("1")
    return [tuple(p) for p in P], info


def extract_sectored(E, rounds=4, seed=0, tries=4):
    """Products of engine E (chem_own_engine.Engine on its current spans) from a parity basis that favours repeated
    point pairs (large sectors): the extraction deletes candidate lines in the given order while Y stays
    nonsingular, so lines late in the order are kept; lines are ordered by how many kept lines share their point
    pair in the previous round (round 0: by how many candidate lines share it). Returns (prods, sector sizes)."""
    import importlib
    from collections import Counter
    eng = importlib.import_module(type(E).__mro__[-2].__module__)   # the module that defines the base Engine
    E._sync()
    r = E.r
    best = None
    for attempt in range(tries):
        xr = np.random.default_rng(1009 * attempt + 17 * E.an.k + 5 + 7919 * seed)
        n = eng._lines(E.pt_blk, E.pt_ymask, E.pt_outs, E.pt_phase, E.ncoord, E.Rm, E.ecount, E.nb, E.li, E.lj,
                       E.ls, E.deg)
        lx = xr.integers(1, 1 << eng.GF_BITS, eng.MAXLINES).astype(np.uint16)
        Y = eng._build_y(n, E.li, E.lj, E.ls, lx, E.pt_blk, E.pt_ymask, E.ncoord, E.cellid, E.freeidx, E.nb, r, r)
        nu = eng.gf_rank(Y.copy(), eng.GF_EXP, eng.GF_LOG) // 2
        d = r - 2 * nu
        ntot = n
        if d:
            ntot = eng._singles(E.pt_blk, E.pt_ymask, E.pt_outs, E.ncoord, E.Rm, E.ecount, E.nb, d, n, E.li, E.lj, E.ls)
        Y = eng._build_y(ntot, E.li, E.lj, E.ls, lx, E.pt_blk, E.pt_ymask, E.ncoord, E.cellid, E.freeidx, E.nb, r,
                         r + d)
        key = [(min(int(E.li[q]), int(E.lj[q])), max(int(E.li[q]), int(E.lj[q]))) if E.lj[q] >= 0 else None
               for q in range(ntot)]
        score = Counter(k for k in key if k is not None)
        for rd in range(rounds):
            ok, N = eng.gf_inverse(Y.copy(), eng.GF_EXP, eng.GF_LOG)
            if not ok:
                break
            tie = xr.random(ntot)
            order = np.array(sorted(range(ntot), key=lambda q: (-1 if key[q] is None else score[key[q]], tie[q])),
                             np.int64)
            keep = eng._extract(N, ntot, E.li, E.lj, E.ls, lx, order, E.pt_blk, E.pt_ymask, E.ncoord, E.cellid,
                                E.freeidx, E.nb, r, eng.GF_EXP, eng.GF_LOG)
            kept = np.flatnonzero(keep)
            pairs = [(int(E.li[q]), int(E.lj[q]), int(E.ls[q])) for q in kept if E.lj[q] >= 0]
            singles = [(int(E.li[q]), int(E.ls[q])) for q in kept if E.lj[q] < 0]
            if 2 * len(pairs) + len(singles) != r:
                continue
            prods = E._products(pairs, singles)
            if prods is None:
                continue
            kc = Counter((min(i, j), max(i, j)) for i, j, _ in pairs)
            multi = sum(c for c in kc.values() if c > 1)
            if best is None or multi > best[1]:
                best = (prods, multi, sorted(kc.values(), reverse=True))
            score = Counter({k: 1000 * kc.get(k, 0) + score.get(k, 0) for k in score})
    if best is None:
        raise RuntimeError("extract_sectored: no plan extracted")
    return best[0], best[2]


# ------------------------------------------------------------------------------------------------ the 68-AND cache
class TSpan:
    """Span of truth tables (ints) with tags (ints, XOR-combined)."""

    def __init__(self):
        self.piv = {}

    def reduce(self, v, tag=0):
        while v:
            hb = v.bit_length() - 1
            if hb not in self.piv:
                return v, tag
            pv, pt = self.piv[hb]
            v ^= pv
            tag ^= pt
        return 0, tag

    def add(self, v, tag):
        r, t = self.reduce(v, tag)
        if r:
            self.piv[r.bit_length() - 1] = (r, t)
            return True
        return False


def build_hcache(V, hbits=7, need=None):
    if need is None:                                   # tenth holder: HB - 2 - hbits (68 at 77 blocks, 64 at 73)
        need = HB - 2 - hbits
    """Greedy cache of `need` ANDs whose outputs, with 1 and the address bits, span H = {f : sum_V f = 0}.
    Wires: index i < hbits = address bit i; index hbits + j = cache AND j. Returns (ops, truths): ops[j] =
    (terms1, terms2), the wire lists whose XORs are the two controls; truths[w] = 77-bit truth of wire w."""
    truths = [sum(1 << p for p in range(HB) if (p >> i) & 1) for i in range(hbits)]
    S = TSpan()
    S.add((1 << HB) - 1, 0)
    for w, f in enumerate(truths):
        S.add(f, 1 << w)
    ops = []
    rep = {1 << i: i for i in range(hbits)}            # monomial -> wire that carries it (up to lower terms)

    def try_pair(t1, t2):
        f1 = 0
        for w in t1:
            f1 ^= truths[w]
        f2 = 0
        for w in t2:
            f2 ^= truths[w]
        if bform(f1, f2, V):
            return None
        g = f1 & f2
        if S.reduce(g)[0] == 0:
            return None
        return g

    def accept(t1, t2, g):
        w = len(truths)
        truths.append(g)
        S.add(g, 1 << w)
        ops.append((list(t1), list(t2)))
        return w

    for M in [M for M in range(1, HB) if M & (M - 1)]:
        if len(ops) >= need:
            break
        top = M.bit_length() - 1
        Mp = M ^ (1 << top)
        if Mp not in rep:
            continue
        a, b = rep[Mp], top
        cand = [([a], [b])] + [([a], [b, c]) for c in range(len(truths)) if c not in (a, b)] + \
            [([a, c], [b]) for c in range(len(truths)) if c not in (a, b)]
        for t1, t2 in cand:
            g = try_pair(t1, t2)
            if g is not None:
                rep[M] = accept(t1, t2, g)
                break
    while len(ops) < need:                             # top up with any pair that still adds a dimension
        n = len(truths)
        found = None
        for x in range(n):
            for y in range(x + 1, n):
                g = try_pair([x], [y])
                if g is not None:
                    found = ([x], [y], g)
                    break
            if found:
                break
        if found is None:                              # controls may be any XOR of wires: random small combinations
            import random
            rr = random.Random(len(ops) * 7919 + n)
            for t in range(200000):
                t1 = rr.sample(range(n), rr.randint(1, 4))
                t2 = rr.sample(range(n), rr.randint(1, 4))
                if set(t1) == set(t2):
                    continue
                g = try_pair(t1, t2)
                if g is not None:
                    found = (t1, t2, g)
                    break
        if found is None:
            raise ValueError(f"68-AND cache: stuck at {len(ops)} ANDs")
        accept(*found)
    assert len(S.piv) == 1 + hbits + need == HB - 1, f"cache spans {len(S.piv)} dimensions"
    assert all(bform(f, (1 << HB) - 1, V) == 0 for f in truths), "a cache function outside H"
    return ops, truths


# ------------------------------------------------------------------------------------------------ the C-4 emitter
def emit_c4(an, prods, V):
    """kmx text and stats of loader k with the 68-AND high cache (copy of chem_loaders_core.emit; the high cache,
    the high parts of the product controls and the pure-high corrections are expressed in the cache basis of H).
    Needs tau == kappa on every wire (chem_ld2_c4.closure); fails loudly otherwise."""
    bits, BROT, J, QUBIT_CAP = core.bits, core.BROT, core.J, core.QUBIT_CAP
    A, lb, L, nl, n_sel = an.A, an.lb, an.L, an.nl, an.n_sel
    assert an.H == HB and an.hbits == 7
    a = np.arange(A)
    lidx, hidx = a & (L - 1), a >> lb
    ALL = (1 << A) - 1
    LBf = []
    for bs in an.basis:
        tt = np.zeros(L, np.uint8)
        for m in bits(bs["anf"]):
            tt[m] = 1
        tt = core.mobius(tt, lb)
        LBf.append(core.to_int(tt[lidx]))

    def addr_truth(g):
        bv = np.array([(g >> p) & 1 for p in range(HB)], np.uint8)
        return core.to_int(bv[hidx])
    HMf = [None] + [core.to_int((hidx & M) == M) for M in range(1, an.H)]

    ops_c, truths = build_hcache(V)
    nwc = len(truths)                                  # 7 address bits + 68 cache wires
    CT = TSpan()                                       # cache basis over the 77 blocks; tag bit w = wire w, bit nwc = 1
    CT.add((1 << HB) - 1, 1 << nwc)
    for w, f in enumerate(truths):
        CT.add(f, 1 << w)

    def hwires(mask):
        r, tag = CT.reduce(hmask_truth(mask))
        assert r == 0, "a product high part is outside H"
        return tag

    def wire_q(w):
        return n_sel + BROT * w[0] + w[1]

    ops = []
    free = sorted(wire_q(c["wire"]) for c in an.carriers if c["clean"]) + list(range(n_sel + J * BROT, QUBIT_CAP))
    used_max = n_sel + J * BROT - 1

    def take():
        nonlocal used_max
        if not free:
            raise ValueError("out of qubits")
        q = free.pop(0)
        used_max = max(used_max, q)
        return q

    def give(q):
        free.append(q)
        free.sort()

    lowq = [None] * nl
    for t, bs in enumerate(an.basis):
        if bs["kind"] == "lin":
            lowq[t] = bs["bit"]
        elif bs["kind"] == "car":
            lowq[t] = wire_q(bs["wire"])
    prep = []
    for c in an.carriers:
        q = wire_q(c["wire"])
        if c["conj"]:
            prep.append(f"CX q{q - 1} q{q}")
        if c["const"]:
            prep.append(f"X q{q}")
        if c["clean"]:
            for t in c["combo"]:
                prep.append(f"CX q{lowq[t]} q{q}")
    ops += prep

    def controls(s1, s2, const1=False, const2=False):
        """Pre-gates that put XOR(s1) (+1) on pivot p1 and XOR(s2) (+1) on pivot p2 (the core product logic)."""
        s1, s2 = set(s1), set(s2)
        assert s1 and s2 and s1 != s2, "degenerate controls"
        p1 = min(s1 - s2) if s1 - s2 else min(s1)
        pre = [f"CX q{q} q{p1}" for q in sorted(s1) if q != p1]
        if p1 in s2:
            s2 = (s1 ^ s2) | {p1}
            const2 = const2 ^ const1
        if const1:
            pre.append(f"X q{p1}")
        p2 = min(q for q in s2 if q != p1)
        pre += [f"CX q{q} q{p2}" for q in sorted(s2) if q != p2]
        if const2:
            pre.append(f"X q{p2}")
        return pre, p1, p2

    cq = [lb + i for i in range(7)]                    # qubit of each cache-basis wire
    made_high = []
    for (t1, t2) in ops_c:
        pre, p1, p2 = controls([cq[w] for w in t1], [cq[w] for w in t2])
        t = take()
        ops.extend(pre)
        ops.append(f"CCX q{p1} q{p2} q{t}")
        ops.extend(reversed(pre))
        cq.append(t)
        made_high.append((pre, p1, p2, t))

    def low_and(t):
        par, i = an.recipe[t]
        qs = [lowq[x] for x in par]
        piv = next(q for q in qs if q != i)
        pre = [f"CX q{q} q{piv}" for q in qs if q != piv]
        tq = take()
        ops.extend(pre)
        ops.append(f"CCX q{piv} q{i} q{tq}")
        ops.extend(reversed(pre))
        lowq[t] = tq
        return (pre, piv, i, tq, t)

    def low_del(rec):
        pre, piv, i, tq, t = rec
        ops.extend(pre)
        ops.append(f"HMR q{tq} b0")
        ops.append(f"CZ q{piv} q{i} if b0")
        ops.extend(reversed(pre))
        lowq[t] = None
        give(tq)

    nw = len(an.out_wires)
    wq = [wire_q(w) for w in an.out_wires]

    def func(low_mask, high_mask):
        f = 0
        for t in bits(low_mask):
            f ^= LBf[t]
        for M in bits(high_mask):
            f ^= HMf[M + 1]
        return f

    acc = [0] * nw
    targets_of = []
    for a1, b1, a2, b2, s, _ in prods:
        P = func(a1, b2) & func(a2, b1)
        tg = [q for q in range(nw) if (an.G[q] & s).bit_count() & 1]
        assert tg, "product without a target"
        targets_of.append(tg)
        for q in tg:
            acc[q] ^= P
    pure = core.Span()
    pure.add(ALL, 1)
    for t in range(nl):
        pure.add(LBf[t], 1 << (1 + t))
    for w, f in enumerate(truths):
        pure.add(addr_truth(f), 1 << (1 + nl + w))
    corr = []
    for q in range(nw):
        r, c = pure.reduce(an.delta_int[q] ^ acc[q])
        assert r == 0, f"wire {an.out_wires[q]}: pure-high part outside H (tau != kappa) on circuit {an.k}"
        corr.append(c)

    core_made = [low_and(t) for t in an.core]
    by_phase = [[] for _ in an.excl]
    for i, p in enumerate(prods):
        by_phase[p[5]].append(i)
    n_emitted = 0

    def do_product(i):
        nonlocal n_emitted
        a1, b1, a2, b2, s, _ = prods[i]
        h2, h1 = hwires(b2), hwires(b1)
        s1 = [lowq[t] for t in bits(a1)] + [cq[w] for w in bits(h2 & ((1 << nwc) - 1))]
        s2 = [lowq[t] for t in bits(a2)] + [cq[w] for w in bits(h1 & ((1 << nwc) - 1))]
        assert None not in s1 and None not in s2, "product uses a low function that is not cached in its phase"
        pre, p1, p2 = controls(s1, s2, bool(h2 >> nwc & 1), bool(h1 >> nwc & 1))
        tg = targets_of[i]
        first = wq[tg[0]]
        fan = [f"CX q{first} q{wq[q]}" for q in tg[1:]]
        ops.extend(pre)
        ops.extend(fan)
        ops.append(f"CCX q{p1} q{p2} q{first}")
        ops.extend(reversed(fan))
        ops.extend(reversed(pre))
        n_emitted += 1

    def low_corrections(ts):
        for q in range(nw):
            for t in ts:
                if corr[q] >> (1 + t) & 1:
                    ops.append(f"CX q{lowq[t]} q{wq[q]}")

    excl_all = {t for ex in an.excl for t in ex}
    for ph, ex in enumerate(an.excl):
        made = [low_and(t) for t in ex]
        for i in by_phase[ph]:
            do_product(i)
        low_corrections(ex)
        if ph == len(an.excl) - 1:
            low_corrections([t for t in range(nl) if t not in excl_all])
            for q in range(nw):
                c = corr[q]
                if c & 1:
                    ops.append(f"X q{wq[q]}")
                for w in range(nwc):
                    if c >> (1 + nl + w) & 1:
                        ops.append(f"CX q{cq[w]} q{wq[q]}")
        for rec in reversed(made):
            low_del(rec)
    for rec in reversed(core_made):
        low_del(rec)
    for pre, p1, p2, t in reversed(made_high):
        ops.extend(pre)
        ops.append(f"HMR q{t} b0")
        ops.append(f"CZ q{p1} q{p2} if b0")
        ops.extend(reversed(pre))
    ops.extend(reversed(prep))

    head = []
    for i in range(n_sel):
        head.append(f"APPEND_TO_REGISTER q{i} r0")
    head.append('REGISTER r0 "selection"')
    for j in range(J):
        for b in range(BROT):
            head.append(f"APPEND_TO_REGISTER q{n_sel + BROT * j + b} r{j + 1}")
        head.append(f'REGISTER r{j + 1} "word_{j}"')
    infra = len(ops_c) + len(an.comp)
    stats = dict(k=an.k, infra=infra, products=n_emitted, toffoli=infra + n_emitted, qubits=used_max + 1,
                 phases=len(an.excl), outputs=len(an.out_level), ops=len(ops), cache="c4-68")
    return "\n".join(head + ops) + "\n", stats
