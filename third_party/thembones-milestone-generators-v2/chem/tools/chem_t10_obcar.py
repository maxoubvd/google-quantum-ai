#!/usr/bin/env python
"""chem_t10_obcar (T10 seat 3, Nov 18; ledger row O14): the 76 ONE-BODY rows under a relation file, ROUNDED DIRECTLY IN
C2's CARRIER FORM (R53-CHEM-GAUGE's rule for the one-body rows), so that every relation holds on the stored words and
every carrier is exact after the table step.

C2's carriers ask, for the row at address a and column t: bit6 XOR bit5 of the stored word = c_t(a & 15) (the slot
monomial, r53_chem_gauge_car.MONO). A row in carrier form is word_t = x_t + 64 c_t with x_t in [-32, 31]; bits 0..5 of
the stored word are those of x_t. So: RK7's staged rounding (chem_rk7_ob.progressive_rel: RelRank's rank rule) and its
relation-keeping polish (polish_rel) run on x in the signed range with the relations on the bits of x (= the stored
bits), and the error_1 objective decodes the TRUE words x + 64 c (half turns c). The rows' addresses are FIXED per try
(RelRank's matching off: c depends on the address); tries draw random placements as table7c's do. Start angles: the
signed-range angles z of the eigenvectors, gauge-moved: x = (-1)^n z, n = N c (r53_chem_gauge_car.gauge_lists; the
decoded vector is unchanged: checked on every row). After the table step's obsigned, C2's unchanged carrier gauge
(chem_ld2_e11car7.py) returns these words (the carrier form of a row is unique).
Copies of chem_rk7_ob.progressive_rel / polish_rel (read-only source) with the offset; everything else imported.

    python chem_t10_obcar.py ROUNDED.npz OUT_DIR --rel @rel.txt [table7c options] [--ob-tries N] [--ob-workers W]
                             [--ob-seed S] [--ob-fine 3] [--rel-out FILE]
    python chem_t10_obcar.py --selftest
"""
import multiprocessing as mp
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_onebody import STAGES, angles_and_mass, decode_jac, decode_pi, error1  # noqa: E402
from chem_rk7_rel import RelRank, check, fmt, parse, ttable  # noqa: E402
from chem_rk7_ob import BITS32, RelFilter, dp_reparity, fine_stages, kernel_sets  # noqa: E402
from r53_chem_gauge_car import MONO, gauge_lists  # noqa: E402

OPT = {"--rel": "", "--ob-tries": "1", "--ob-workers": "1", "--ob-seed": "0", "--ob-fine": "3", "--rel-out": "",
       "--ob-sweeps": "0", "--ob-match": "samel", "--ob-btab": "1", "--ob-mode": "ge"}
LO, HI = -32, 31


def ttable_x(t, addrs):
    """chem_rk7_rel.ttable plus the kind ('FB', base, x): the base target XOR a free constant per block on blocks 0..3
    (bit b of x = the constant of block b; addresses 0..63 hold one-body rows only, so the loader's family F takes it),
    and R65's 'G/<base>/<hex>' targets."""
    if isinstance(t, str) and t.startswith("G/"):
        from r65_chem_all15_relg import ttable_g
        return ttable_g(t, addrs)
    if isinstance(t, tuple) and t and t[0] == "FB":
        a_ = np.asarray(addrs, np.int64)
        blk = a_ >> 4
        return ttable(t[1], a_) ^ np.where(blk < 4, (int(t[2]) >> np.minimum(blk, 3)) & 1, 0).astype(np.int64)
    return ttable(t, addrs)


def btab_rels(rels, rng):
    """every relation's target XOR a random constant on each of the blocks 0..3 (one draw per try)."""
    return [(k, o, ("FB", m, int(rng.integers(0, 16)))) for (k, o, m) in rels]


def fmt_x(rels):
    """relations with ('FB', base, x) targets in R65's G/ form (block constants hex; blocks 0..3 from x)."""
    out = []
    for (k, o, m) in rels:
        if isinstance(m, tuple):
            base, x = m[1], int(m[2])
            out.append((k, o, f"G/{base}/{x:x}" if x else base))
        else:
            out.append((k, o, m))
    return fmt(out)


def carriers_at(addrs, ncol=75):
    low = np.asarray(addrs, np.int64) % 16
    return np.stack([((low & MONO[t % 5]) == MONO[t % 5]).astype(np.int64) for t in range(ncol)], 1)


def carrier_start(z, edges, addrs):
    """signed-range angles z (rows, 75; units of pi) at addresses -> (x, c): x = (-1)^n z (units of pi), c the half
    turns (0/1), n = N c: the same vector (the gauge of r53_chem_gauge_car)."""
    c = carriers_at(addrs, z.shape[1])
    deps = gauge_lists(edges)
    n = np.zeros_like(c)
    for t in range(z.shape[1]):
        if deps[t]:
            n[:, t] = c[:, deps[t]].sum(1) & 1
    return np.where(n == 1, -z, z), c


def check_x(words, addrs, rels):
    """chem_rk7_rel.check with ttable_x targets: (relation, row) pairs broken (rows in address order)."""
    w = np.asarray(words, np.int64)
    bad = 0
    for (k, o, m) in rels:
        acc = np.zeros(len(w), np.int64)
        for (j, b) in o:
            acc ^= ((w[:, 5 * k + j] ^ w[:, 5 * k - 5 + j]) >> b) & 1
        bad += int(np.sum(acc != ttable_x(m, addrs)))
    return bad


def progressive_car(th0, C, edges, lam_at, H, scale, iters, rels, stages=STAGES, match="samel"):
    """chem_rk7_ob.progressive_rel on x in [-32, 31] with the half turns C (units of pi) added in every decode. The
    address of a row moves only among the addresses with its own low four bits (match 'samel': the same carriers, so
    the same x), or not at all (match 'none')."""
    nrow = th0.shape[0]
    base = [(k, o, (m[1] if isinstance(m, tuple) else m)) for (k, o, m) in rels]
    dyn = RelRank(base, np.arange(nrow), nrow, 0)
    dyn.tabs = [ttable_x(m, dyn.addrs).astype(bool) for (_k, _o, m) in rels]       # t10: block constants on 0..3
    dyn.zero = not any(np.any(t_) for t_ in dyn.tabs)
    dyn._L = {}
    lowa = np.arange(dyn.na) & 15
    if match == "samel":
        dyn.allowed = (lowa[None, :] == (np.arange(nrow) & 15)[:, None])
    else:
        dyn.allowed = np.eye(nrow, dyn.na, dtype=bool)
    theta = th0.copy()
    frozen = np.zeros(theta.shape, bool)
    qv = np.zeros(theta.shape)
    lam2 = lam_at[:, None] ** 2
    ncol = theta.shape[1]

    def fun(x, mask, base):
        t = base.copy()
        t[mask] = x
        V, J = decode_jac(t + C, edges)
        D = (V.T * lam_at) @ V - H
        gV = 4.0 * lam_at[:, None] * (V @ D)
        return float(np.sum(D * D)), np.einsum("jn,jkn->jk", gV, J)[mask]

    for frac in stages:
        need = int(round(frac * theta.size)) - int(frozen.sum())
        if need > 0:
            now = np.where(frozen, qv, theta)
            mass = angles_and_mass(decode_pi(now + C, edges), edges)[1]
            xs = now * scale
            q = np.clip(np.round(xs), LO, HI)
            impact = lam2 * mass * (q / scale - now) ** 2
            impact[frozen] = -np.inf
            pick = np.argsort(-impact, axis=None)[:need]
            for flat in pick:
                j, c = divmod(int(flat), ncol)
                zc = int(q[j, c])
                v = dyn.choose(j, c, float(xs[j, c]), LO, HI)
                if v is not None:
                    zc = int(v)
                dyn.commit(j, c, zc)
                frozen[j, c] = True
                qv[j, c] = zc / scale
        free = ~frozen
        if free.any():
            base = np.where(frozen, qv, theta)
            r = minimize(fun, theta[free], args=(free, base), jac=True, method="L-BFGS-B",
                         options={"maxiter": iters, "maxcor": 20, "ftol": 1e-14, "gtol": 1e-10})
            theta = base.copy()
            theta[free] = r.x
        else:
            theta = qv.copy()
    return np.round(qv * scale).astype(np.int64), dyn


def polish_car(wq, C, edges, H, brot, sweeps, rels, addrs, inner=25, maxsize=4, dp=True):
    """chem_rk7_ob.polish_rel (non-full range) on x in [-32, 31] with the half turns 64 C added in every decode."""
    import itertools
    scale = float(1 << (brot - 1))
    N = edges.shape[0] + 1
    ncol = wq.shape[1]
    wq = np.array(wq, np.int64)
    C64 = 64 * np.asarray(C, np.int64)
    filt = RelFilter(rels, ncol)
    target = np.stack([ttable_x(m, addrs) for (_, _, m) in rels], axis=1) if rels else np.zeros((len(wq), 0), np.int64)
    free, ksets = kernel_sets(rels, ncol, maxsize)
    dp = dp and any(all(b == 0 for (_, b) in o) for (_, o, _) in rels)
    by_k = {}
    for ri, (k, o, _) in enumerate(rels):
        if all(b == 0 for (_, b) in o):
            by_k.setdefault(k, []).append((sum(1 << j for (j, _) in o), ri))
    hib = np.zeros(ncol, np.int64)
    for (k, o, _) in rels:
        for (j, b) in o:
            for c_ in (5 * k + j, 5 * k - 5 + j):
                hib[c_] = max(hib[c_], b)
    keepm = ((1 << hib) - 1)
    sdel = np.array([d for d in range(-3, 4) if d], np.int64)
    dpn = dpa = 0
    moves = []
    for t in range(ncol):
        for d in (1, -1, 2, -2, 4, -4):
            moves.append(((t,), (d,)))
    for S in ksets:
        for sg in itertools.product((1, -1), repeat=len(S)):
            moves.append((S, sg))
    # t10: BIT-1 KERNEL MOVES: +-2 on two (or three) words whose bit-1 syndromes cancel keep every bit-1 relation (and
    # bit 0): without them the polish cannot touch bit 1 of a relation word at all
    syn1 = [0] * ncol
    for i_, (k_, o_, _m) in enumerate(rels):
        for (j_, b_) in o_:
            if b_ == 1:
                syn1[5 * k_ + j_] ^= 1 << i_
                syn1[5 * k_ - 5 + j_] ^= 1 << i_
    cols1 = [c_ for c_ in range(ncol) if syn1[c_]]
    by1 = {}
    for c_ in cols1:
        by1.setdefault(syn1[c_], []).append(c_)
    n1 = 0
    for grp in by1.values():
        for a_, b_ in itertools.combinations(grp, 2):
            for sg in itertools.product((2, -2), repeat=2):
                moves.append(((a_, b_), sg))
                n1 += 1
    if os.environ.get("T10_OB_K1TRI", "1") == "1":
        seen3 = set()
        for a_, b_ in itertools.combinations(cols1, 2):
            for c_ in by1.get(syn1[a_] ^ syn1[b_], []):
                t3 = tuple(sorted((a_, b_, c_)))
                if c_ in (a_, b_) or t3 in seen3:
                    continue
                seen3.add(t3)
                for sg in itertools.product((2, -2), repeat=3):
                    moves.append((t3, sg))
                    n1 += 1
    Dm = np.zeros((len(moves), ncol), np.int64)
    for i, (S, sg) in enumerate(moves):
        for c, d in zip(S, sg):
            Dm[i, c] += d
    V = decode_pi(np.mod(wq + C64, 1 << brot) / scale, edges)
    e1, perm, D, lam = error1(V, H)
    best_e, best_wq = e1, wq.copy()
    stats = dict(kernel_sets=len(ksets), free_cols=len(free), moves=len(moves), moved=0, sweeps=0, e_in=float(e1))

    def rowval(Vc, base, lj):
        quad = np.einsum("ip,pq,iq->i", Vc, base, Vc)
        n4 = np.sum(Vc * Vc, axis=1) ** 2
        return 2 * lj * quad + lj * lj * n4

    for sweep in range(sweeps):
        moved = 0
        for j in range(N):
            lj = lam[perm[j]]
            base = D - lj * np.outer(V[j], V[j])
            cj = C64[j][None, :]
            if dp:
                v0 = V[j]
                cur = 2 * lj * float(v0 @ base @ v0) + lj * lj * float(v0 @ v0) ** 2
                Cs = np.repeat(wq[j][None, :], ncol * len(sdel), axis=0)
                ci = np.repeat(np.arange(ncol), len(sdel))
                Cs[np.arange(len(Cs)), ci] += np.tile(sdel, ncol)
                okc = (Cs[np.arange(len(Cs)), ci] >= LO) & (Cs[np.arange(len(Cs)), ci] <= HI)
                vals = np.full(len(Cs), np.inf)
                vals[okc] = rowval(decode_pi(np.mod(Cs[okc] + cj, 1 << brot) / scale, edges), base, lj) - cur
                vals = vals.reshape(ncol, len(sdel))
                newv = Cs[np.arange(len(Cs)), ci].reshape(ncol, len(sdel))
                if np.any(hib):
                    sameh = (((np.mod(newv, 128) >> 1) ^ (np.mod(wq[j], 128)[:, None] >> 1)) & keepm[:, None]) == 0
                    vals = np.where(sameh, vals, np.inf)
                par = np.mod(newv, 2)
                cost = np.zeros((ncol, 2))
                pick = np.zeros((ncol, 2), np.int64)
                p0 = np.mod(wq[j], 2)
                for p in (0, 1):
                    vv = np.where(par == p, vals, np.inf)
                    a_ = np.argmin(vv, axis=1)
                    cost[:, p] = vv[np.arange(ncol), a_]
                    pick[:, p] = newv[np.arange(ncol), a_]
                    same = p0 == p
                    better = same & (cost[:, p] > 0)
                    cost[better, p] = 0.0
                    pick[better, p] = wq[j][better]
                pats, tc = dp_reparity(cost, by_k, target[j])
                newbits = BITS32[pats].reshape(-1)[:ncol]
                if np.any(newbits != p0) and tc < -1e-13:
                    dpn += 1
                    c = wq[j].copy()
                    ch = newbits != p0
                    c[ch] = pick[ch, newbits[ch]]
                    if not rels or np.all(filt.parities(c[None, :])[0] == target[j]):
                        Vn = decode_pi(np.mod(c[None, :] + cj, 1 << brot) / scale, edges)
                        if rowval(Vn, base, lj)[0] < cur - 1e-13:
                            wq[j] = c
                            V[j] = Vn[0]
                            dpa += 1
                            moved += 1
            for _ in range(inner):
                Cc = wq[j][None, :] + Dm
                ok = np.all((Cc >= LO) & (Cc <= HI), axis=1)
                if rels:
                    ok &= np.all(filt.parities(Cc) == target[j][None, :], axis=1)
                Cc = Cc[ok]
                if not len(Cc):
                    break
                Vc = decode_pi(np.mod(Cc + cj, 1 << brot) / scale, edges)
                v0 = V[j]
                val = rowval(Vc, base, lj)
                cur = 2 * lj * float(v0 @ base @ v0) + lj * lj * float(v0 @ v0) ** 2
                i = int(np.argmin(val))
                if val[i] >= cur - 1e-13:
                    break
                wq[j] = Cc[i]
                V[j] = Vc[i]
                moved += 1
            D = base + lj * np.outer(V[j], V[j])
        e_prev = e1
        e1, perm, D, lam = error1(V, H)
        stats["moved"] += moved
        stats["sweeps"] = sweep + 1
        if e1 < best_e:
            best_e, best_wq = e1, wq.copy()
        if not moved or e_prev - e1 < 1e-7 * e1:
            break
    stats.update(dp_tried=dpn, dp_taken=dpa, e_out=float(best_e))
    return best_wq, best_e, stats


def one_try_car(args):
    """(H, edges, brot, iters, sweeps, rels, t, seed, fine[, match, btab]) -> (error_1, STORED words in address order,
    info). btab: tries t > 0 draw a free constant per (relation, block 0..3); match: 'samel' or 'none'."""
    H, edges, brot, iters, sweeps, rels, t, seed, fine = args[:9]
    match = args[9] if len(args) > 9 else "samel"
    btab = bool(args[10]) if len(args) > 10 else True
    from r48_chem_q7c_obsigned import signed_angles
    if btab and t > 0:
        rels = btab_rels(rels, np.random.default_rng([int(seed), int(t), 91]))
    scale = float(1 << (brot - 1))
    lam_e, vec = np.linalg.eigh(H)
    zs, _ = signed_angles(vec.T, edges)
    n = zs.shape[0]
    perm = np.arange(n) if t == 0 else np.random.default_rng(1000 * seed + t).permutation(n)
    x0, C = carrier_start(zs[perm], edges, np.arange(n))
    v_a = decode_pi(zs[perm], edges)
    v_b = decode_pi(x0 + C, edges)
    dv = float(np.max(np.minimum(np.abs(v_a - v_b).max(1), np.abs(v_a + v_b).max(1))))
    assert dv < 1e-9, f"the carrier start changed a vector ({dv})"
    st = fine_stages(STAGES, fine) if fine > 1 else STAGES
    t0 = time.time()
    xq, dyn = progressive_car(x0, C, edges, lam_e[perm], H, scale, iters, rels, st, match)
    assert np.array_equal(np.arange(n) & 15, dyn.mrow & 15), "a one-body row moved to another low address"
    xa = np.zeros_like(xq)
    xa[dyn.mrow] = xq                                  # row -> its address
    CA = carriers_at(np.arange(n), xq.shape[1])        # = C of each row at its new address (same low four bits)
    stored0 = np.mod(xa + 64 * CA, 1 << brot)
    bad0 = check_x(stored0, np.arange(n), rels)
    assert bad0 == 0, f"{bad0} relations broken after the staged rounding"
    e_pre = error1(decode_pi(stored0 / scale, edges), H)[0]
    xa, _, st_ = polish_car(xa, CA, edges, H, brot, max(sweeps, 1), rels, np.arange(n))
    stored = np.mod(xa + 64 * CA, 1 << brot)
    bad = check_x(stored, np.arange(n), rels)
    assert bad == 0, f"{bad} relations broken on the one-body rows"
    car = ((stored >> 6) ^ (stored >> 5)) & 1
    assert np.array_equal(car, CA), "a carrier is wrong"
    e = error1(decode_pi(stored / scale, edges), H)[0]
    info = dict(free=int(dyn.free_ok), forced=int(dyn.forced_n), completions=int(dyn.ncons), e_staged=float(e_pre),
                moved_rows=int((dyn.mrow != np.arange(n)).sum()), seconds=round(time.time() - t0, 1), **st_)
    info["rels"] = fmt_x(rels)
    return float(e), stored, info


def progressive_ge(th0, edges, lam_at, H, scale, iters, rels, stages=STAGES):
    """MODE ge: the one-body rows' staged rounding (as progressive_rel: stiffest first, L-BFGS on error_1 between
    stages) in the SIGNED range z, with the relations held by RelPlacerGE (tools/chem_t10_relge.py) as ONE group of 76
    places = the 76 one-body addresses: each row's address is chosen by the placer's assignment (any address), every
    relation is asked of y = eps x z (eps per address, exact recursion), and the block constants of every relation on
    blocks 0..3 are fitted (block 4 pinned: it is shared with two-body rows)."""
    from chem_t10_relge import RelPlacerGE
    nrow, ncol = th0.shape
    fp = RelPlacerGE(1, nrow, [(k, tuple(o), m) for (k, o, m) in rels], edges, n_one=0, ncols=ncol, rounds=1,
                     pin_blocks=(4,))
    theta = th0.copy()
    frozen = np.zeros(theta.shape, bool)
    qv = np.zeros(theta.shape)
    lam2 = lam_at[:, None] ** 2

    def fun(x, mask, base):
        t = base.copy()
        t[mask] = x
        V, J = decode_jac(t, edges)
        D = (V.T * lam_at) @ V - H
        gV = 4.0 * lam_at[:, None] * (V @ D)
        return float(np.sum(D * D)), np.einsum("jn,jkn->jk", gV, J)[mask]

    forced = 0
    REGAUGE = os.environ.get("T10_OB_REGAUGE", "0") == "1"      # measured 10:55Z: 8 broken after the rounding (1 without): off
    PLUSPOLE = os.environ.get("T10_OB_PLUSPOLE", "0") == "1"    # measured 11:10Z on the real h1': no valid try (fine 1): off
    pluspole = [0]
    deps_ = gauge_lists(edges)
    governed = [[s_ for s_ in range(ncol) if t_ in deps_[s_]] for t_ in range(ncol)]    # a half turn on t negates these
    regauged = [0]
    if REGAUGE:
        from r65_chem_all15_frelg import regauge as _regauge     # read-only import: the rounding tool's exact re-gauging
    POLE_FIRST = float(os.environ.get("T10_OB_POLE", "0"))     # words within this many steps of -32 freeze first, at -32
    if POLE_FIRST > 0:
        xs0 = theta * scale
        pole = (xs0 <= LO + POLE_FIRST) | (xs0 >= HI + 1 - POLE_FIRST)
        # a word at +32 - d is the same angle as -32 + (64 - ...) only through a half turn: take the pole -32 only from the
        # negative side; on the positive side the nearest in range (31) is kept, the half turn left to the gauge
        neg = pole & (xs0 < 0)
        frozen[neg] = True
        qv[neg] = LO / scale
        forced += int(neg.sum())
    for frac in stages:
        need = int(round(frac * theta.size)) - int(frozen.sum())
        if need > 0:
            now = np.where(frozen, qv, theta)
            mass = angles_and_mass(decode_pi(now, edges), edges)[1]
            xs = now * scale
            q = np.clip(np.round(xs), LO, HI)
            impact = lam2 * mass * (q / scale - now) ** 2
            impact[frozen] = -np.inf
            pick = np.argsort(-impact, axis=None)[:need]
            cst = lam2 * mass + 1e-30
            r, m = fp.residues(xs, cst, frozen, np.round(qv * scale).astype(np.int64))
            fp.infeasible_max = max(getattr(fp, "infeasible_max", 0), fp.infeasible)     # a late pole: counted, the try is judged at its end
            i_, j_ = np.unravel_index(pick, theta.shape)
            mm = m[j_]
            zz = q[i_, j_].copy()
            sel = mm >= 2
            zraw = np.where(sel, r[i_, j_] + mm * np.round((xs[i_, j_] - r[i_, j_]) / np.maximum(mm, 1)), np.round(xs[i_, j_]))
            zr = np.where(zraw > HI, zraw - mm, np.where(zraw < LO, zraw + mm, zraw))
            zz[sel] = zr[sel]
            up = np.flatnonzero(zraw == HI + 1) if PLUSPOLE else np.zeros(0, np.int64)
            zz[up] = LO
            forced += int(np.sum(zz != q[i_, j_]))
            frozen[i_, j_] = True
            qv[i_, j_] = zz / scale
            # t10: a best value of +32 is taken as the pole by an EXACT gauge move (half turn: -32, and the words it governs
            # negated, frozen or free, this stage's picks included): C2's stored words keep their bits, no clip to 31
            for u_ in up:
                rr_, tt_ = int(i_[u_]), int(j_[u_])
                for s_ in governed[tt_]:
                    qv[rr_, s_] = -qv[rr_, s_]
                    theta[rr_, s_] = -theta[rr_, s_]
                pluspole[0] += 1
        free = ~frozen
        if free.any():
            base = np.where(frozen, qv, theta)
            res_ = minimize(fun, theta[free], args=(free, base), jac=True, method="L-BFGS-B",
                            options={"maxiter": iters, "maxcor": 20, "ftol": 1e-14, "gtol": 1e-10})
            theta = base.copy()
            theta[free] = res_.x
            if REGAUGE:     # t10: free angles that drifted past +-pi/2 back into the signed range by exact gauge moves
                allw = np.where(frozen, qv, theta)          # (half turns + the negations they imply; C2's stored words,
                ones = np.ones(allw.shape, bool)            # hence every relation on them, are unchanged)
                out_, nfl, nng, ncl = _regauge(allw, frozen, ones, ones, edges, -0.5, 0.5)
                qv = np.where(frozen, out_, qv)
                theta = np.where(frozen, theta, out_)
                regauged[0] += nfl
        else:
            theta = qv.copy()
    fp.regauged = regauged[0]
    fp.pluspole = pluspole[0]
    return np.round(qv * scale).astype(np.int64), fp, forced


def repair_rows(xA, CA, edges, H, brot, rels, addrs):
    """rows whose relations are broken (a late pole flipped an eps): per row the move set of least row error_1 that makes
    EVERY relation of the row hold: single moves +-1, +-2, +-4 of any word, then pairs of +-2 moves on the words of its
    broken relations; x stays in [-32, 31] (carriers kept). Returns (x, rows repaired, rows left)."""
    import itertools
    scale = float(1 << (brot - 1))
    x = np.array(xA, np.int64)
    C64 = 64 * np.asarray(CA, np.int64)
    ncol = x.shape[1]
    filt = RelFilter(rels, ncol)
    target = np.stack([ttable_x(m, addrs) for (_, _, m) in rels], axis=1)
    V = decode_pi(np.mod(x + C64, 1 << brot) / scale, edges)
    e1, perm, D, lam = error1(V, H)
    bad = np.flatnonzero(np.any(filt.parities(x) != target, axis=1))
    fixed, left = 0, 0
    b1cols = sorted({5 * k + jj for (k, o, _m) in rels for (jj, b) in o if b == 1} |
                    {5 * k - 5 + jj for (k, o, _m) in rels for (jj, b) in o if b == 1})
    for j in bad:
        supp = b1cols                                      # pairs over every bit-1 relation column of the row
        cands = []
        for c in range(ncol):
            for d in (1, -1, 2, -2, 4, -4):
                v = x[j].copy()
                v[c] += d
                cands.append(v)
        for c1, c2 in itertools.combinations(supp, 2):
            for d1 in (2, -2):
                for d2 in (2, -2):
                    v = x[j].copy()
                    v[c1] += d1
                    v[c2] += d2
                    cands.append(v)
        Cc = np.array(cands)
        ok = np.all((Cc >= LO) & (Cc <= HI), axis=1) & np.all(filt.parities(Cc) == target[j][None, :], axis=1)
        if not ok.any():
            left += 1
            continue
        Cc = Cc[ok]
        lj = lam[perm[j]]
        base = D - lj * np.outer(V[j], V[j])
        Vc = decode_pi(np.mod(Cc + C64[j][None, :], 1 << brot) / scale, edges)
        val = 2 * lj * np.einsum("ip,pq,iq->i", Vc, base, Vc) + lj * lj * np.sum(Vc * Vc, axis=1) ** 2
        i = int(np.argmin(val))
        x[j] = Cc[i]
        V[j] = Vc[i]
        D = base + lj * np.outer(V[j], V[j])
        fixed += 1
    return x, fixed, left


def one_try_ge(args):
    """MODE ge: (H, edges, brot, iters, sweeps, rels, t, seed, fine) -> (error_1, STORED words in address order, info)."""
    H, edges, brot, iters, sweeps, rels, t, seed, fine = args[:9]
    from r48_chem_q7c_obsigned import signed_angles
    from chem_t10_relge import stored_words
    scale = float(1 << (brot - 1))
    lam_e, vec = np.linalg.eigh(H)
    zs, _ = signed_angles(vec.T, edges)
    n = zs.shape[0]
    perm = np.arange(n) if t == 0 else np.random.default_rng(1000 * seed + t).permutation(n)
    st = fine_stages(STAGES, fine) if fine > 1 else STAGES
    t0 = time.time()
    zq, fp, forced = progressive_ge(zs[perm], edges, lam_e[perm], H, scale, iters, rels, st)
    od = fp.order(np.arange(n))                        # od[address] = row
    zA = zq[od]
    stored0 = stored_words(zA, edges, np.arange(n), brot)
    used = fp.rels_used()
    bad0 = check_x(stored0, np.arange(n), used)
    CA = carriers_at(np.arange(n), zq.shape[1])
    xA = np.mod(stored0 - 64 * CA + 64, 128) - 64
    assert np.all((xA >= LO) & (xA <= HI)), "a stored word outside its carrier class"
    e_pre = error1(decode_pi(stored0 / scale, edges), H)[0]
    rep = (0, 0)
    if bad0:      # rows whose relations a late pole broke: the cheapest repair by the row's own error_1, before the polish
        xA, f_, l_ = repair_rows(xA, CA, edges, H, brot, used, np.arange(n))
        rep = (f_, l_)
    xA, _, st_ = polish_car(xA, CA, edges, H, brot, max(sweeps, 1), used, np.arange(n))
    stored = np.mod(xA + 64 * CA, 1 << brot)
    bad = check_x(stored, np.arange(n), used)
    car = ((stored >> 6) ^ (stored >> 5)) & 1
    assert np.array_equal(car, CA), "a carrier is wrong"
    e = error1(decode_pi(stored / scale, edges), H)[0]
    info = dict(forced=int(forced), broken_after_rounding=int(bad0), repaired_rows=rep[0], unrepaired_rows=rep[1],
                broken=int(bad), e_staged=float(e_pre),
                moved_rows=int((od != np.arange(n)).sum()), gflips=int(fp.gflips), pluspole=int(getattr(fp, "pluspole", 0)),
                seconds=round(time.time() - t0, 1), **st_)
    from r65_chem_all15_relg import fmt_g
    info["rels"] = fmt_g(used)
    if bad:
        return float("inf"), stored, info
    return float(e), stored, info


def build_rows7w_car(H, edges, brot, iters, sweeps):
    import r50_chem_lbl_table7c as T7
    rels = parse(OPT["--rel"])
    n = int(OPT["--ob-tries"])
    if int(OPT["--ob-sweeps"]):
        sweeps = int(OPT["--ob-sweeps"])
    if iters == 0:      # table7c's provisional rows (replaced by round 1): one quick try, no polish sweeps beyond one
        n, iters, sweeps = 1, 40, 1
    args = [(H, edges, brot, iters, sweeps, rels, t, int(OPT["--ob-seed"]), int(OPT["--ob-fine"]), OPT["--ob-match"],
             int(OPT["--ob-btab"])) for t in range(n)]
    nw = max(1, min(int(OPT["--ob-workers"]), n))
    fn = one_try_ge if OPT["--ob-mode"] == "ge" else one_try_car
    if nw > 1:
        with mp.get_context("fork").Pool(nw) as pool:
            res = pool.map(fn, args)
    else:
        res = [fn(x) for x in args]
    es = np.array([r[0] for r in res])
    assert np.isfinite(es).any(), "no one-body try kept every relation"
    b = int(np.argmin(es))
    T7.log(f"[table7c] one-body rows IN CARRIER FORM (t10 obcar, mode {OPT['--ob-mode']}) under {len(rels)} relations, {n} tries (fine "
           f"{OPT['--ob-fine']}): error_1 best {1e3 * es[b]:.3f} mHa (try {b}), median {1e3 * float(np.median(es)):.3f}; "
           f"{ {k_: v_ for k_, v_ in res[b][2].items() if k_ != 'rels'} }")
    if OPT["--rel-out"]:
        open(OPT["--rel-out"], "w").write("# relations kept by the one-body rows (chem_t10_obcar, carrier form)\n" + res[b][2]["rels"])
    return res[b][1], res[b][0]


def selftest():
    """the carrier start keeps every vector; a one-body rounding of a random H under a bit-1 layout: relations on the
    stored words and carriers exact; C2's move law (as r53_chem_gauge_car's selftest) leaves the rows unchanged."""
    from r53_chem_gauge_car import carrier_gauge
    rng = np.random.default_rng(3)
    N = 76
    edges = np.array([(int(rng.integers(0, t + 1)), t + 1) for t in range(N - 1)], np.int64)
    A = rng.normal(size=(N, N))
    H = (A + A.T) / 2
    rels = parse("6 | 0:1 1:1 2:1 3:1 4:1 | F:6150;8 | 0:1 1:1 2:1 3:1 4:1 | F:8150;9 | 0:1 1:1 2:1 3:1 4:1 | F:9150;"
                 "8 | 0:0 | 101;8 | 1:0 | 110;9 | 2:0 | 4")
    e, st, info = one_try_car((H, edges, 7, 60, 1, rels, 0, 0, 1))
    assert check(st, np.arange(N), rels) == 0
    want = carriers_at(np.arange(N))
    st2, n2, h2 = carrier_gauge(st, edges, want)
    print(f"selftest: error_1 {1e3 * e:.3f} (random H), relations broken on the stored words 0, carriers exact; "
          f"C2's carrier form of the result equal word for word: {bool(np.array_equal(st2, st))}; {info['seconds']} s")
    assert np.array_equal(st2, st)


def main():
    if "--selftest" in sys.argv:
        selftest()
        print("selftest passed")
        return
    argv = []
    it = iter(sys.argv[1:])
    for x in it:
        if x in OPT:
            OPT[x] = next(it)
        else:
            argv.append(x)
    assert OPT["--rel"], "--rel @file is required"
    import r50_chem_lbl_table7c as T7
    sys.argv = [sys.argv[0]] + argv
    T7.WIRES.append(("relations", OPT["--rel"]))
    T7.build_rows7w = build_rows7w_car
    relax = float(os.environ.get("T10_OB_RELAX", "0"))
    if relax > 0:       # t10: the FIRST refit with the rows fixed (table7c round 1) takes error_1 <= t1 + RELAX, so that round 2
        _refit = T7.cor.refit      # re-rounds the rows on an h1' already moved toward them; every later refit is exact
        ncall = [0]

        def refit_relaxed(prob, w, W, t2, t1, *ra, **kw):
            if t1 < 1e5:
                ncall[0] += 1
                if ncall[0] == 1:
                    T7.log(f"[table7c] t10: round-1 refit relaxed to error_1 <= {t1 + relax:.4f} mHa")
                    return _refit(prob, w, W, t2, t1 + relax, *ra, **kw)
            return _refit(prob, w, W, t2, t1, *ra, **kw)
        T7.cor.refit = refit_relaxed
    free_t2 = float(os.environ.get("T10_FREE_T2", "0"))
    if free_t2 > 0:     # t10 seat 6 (Nov 19 08:50Z): the FREE refit (no error_1 limit) held to error_2 <= FREE_T2 mHa. On a
        _refit0 = T7.cor.refit      # table with error_2 room (y46a: 19.766) the free refit spent all of it on lambda and drove
                                    # error_1 to 156 mHa: the rows were then rounded on a far-moved h1' (24.2 staged, pB1q 21.8)

        def refit_free_capped(prob, w, W, t2, t1, *ra, **kw):
            if t1 >= 1e5 and free_t2 < t2:
                T7.log(f"[table7c] t10: the free refit held to error_2 <= {free_t2:.4f} mHa (room kept for the rows)")
                return _refit0(prob, w, W, free_t2, t1, *ra, **kw)
            return _refit0(prob, w, W, t2, t1, *ra, **kw)
        T7.cor.refit = refit_free_capped
    T7.main()


if __name__ == "__main__":
    main()
