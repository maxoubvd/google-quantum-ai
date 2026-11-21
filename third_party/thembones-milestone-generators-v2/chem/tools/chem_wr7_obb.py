"""chem_wr7_obb (WR7): the 76 one-body rows under delta wires on bits 0 AND 1.

progressive_dynb: the designer's progressive_dyn (r50_chem_lbl_obchain.py) with chem_wr7_chainb.DynChainB: a wired
  word takes the nearest word whose bits keep the rows matchable to the addresses 0..75 (inside [lo, hi]).
polishb: exact coordinate descent on error_1 (polish_coded's objective) with the moves that keep every relation:
  an unwired word +-1; a wired word +-2 (bit-0 wires only on its column) or +-4 (bit 1 too); pair moves (dx, dy) in
  -3..3 on both words of a wired pair, kept when every wired bit of the pair keeps its XOR.
one_try_b: one placement (random initial row-address assignment, as chem_wr7_table7c.one_try), fine schedule, polish.
"""
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_onebody import STAGES, angles_and_mass, decode_jac, decode_pi, error1  # noqa: E402
from chem_wr7_chainb import DynChainB, check  # noqa: E402
from chem_wr7_table7c import fine_stages  # noqa: E402


def progressive_dynb(th0, edges, lam_at, H, scale, iters, wires, stages=STAGES, lo=-64, hi=63):
    nrow = th0.shape[0]
    dyn = DynChainB(wires, np.arange(nrow), nrow)
    theta = th0.copy()
    frozen = np.zeros(theta.shape, bool)
    qv = np.zeros(theta.shape)
    lam2 = lam_at[:, None] ** 2
    ncol = theta.shape[1]
    gauge = hi - lo + 1 < 128      # signed range: the rounding tool's exact regauge after every refit (fix round 1)
    if gauge:
        from r50_chem_lbl_round import negatable, neighbour_matrix, regauge
        NB = neighbour_matrix(edges)
        c1 = np.zeros(ncol, bool)
        for (k_, j_, b_, _) in wires:
            if b_ >= 1:
                c1[[5 * k_ - 5 + j_, 5 * k_ + j_]] = True
        glo, ghi = -0.5 - 0.5 / scale, 0.5 - 0.5 / scale      # the tool's range in units of pi (rounds into lo..hi)

    def fun(x, mask, base):
        t = base.copy()
        t[mask] = x
        V, J = decode_jac(t, edges)
        D = (V.T * lam_at) @ V - H
        gV = 4.0 * lam_at[:, None] * (V @ D)
        return float(np.sum(D * D)), np.einsum("jn,jkn->jk", gV, J)[mask]

    for frac in stages:
        need = int(round(frac * theta.size)) - int(frozen.sum())
        if need > 0:
            now = np.where(frozen, qv, theta)
            mass = angles_and_mass(decode_pi(now, edges), edges)[1]
            xs = now * scale
            if hi - lo + 1 >= 128:      # the full circle wraps
                q = np.mod(np.round(xs) - lo, hi - lo + 1) + lo
            else:
                q = np.clip(np.round(xs), lo, hi)
            impact = lam2 * mass * (q / scale - now) ** 2
            impact[frozen] = -np.inf
            pick = np.argsort(-impact, axis=None)[:need]
            for flat in pick:
                j, c = divmod(int(flat), ncol)
                zc = int(q[j, c])
                v = dyn.choose(j, c, float(xs[j, c]), lo, hi)
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
        if gauge:       # exact half turns into the signed range; odd frozen words of bit-1 columns never negated
            wv = np.round(qv * scale).astype(np.int64)
            prot = frozen & c1[None, :] & (np.mod(wv, 2) != 0)
            negat, flip_ok = negatable(frozen, wv, prot, NB)
            theta, _, _, _ = regauge(np.where(frozen, qv, theta), frozen, flip_ok, negat, edges, glo, ghi)
            qv = np.where(frozen, theta, qv)
    return np.round(qv * scale).astype(np.int64), dyn


def polishb(wq, edges, H, brot, sweeps, wires, inner=25, lo=-64, hi=63):
    scale = float(1 << (brot - 1))
    N = edges.shape[0] + 1
    ncol = wq.shape[1]
    wq = np.array(wq, np.int64)
    bits_of = {}
    pairs = {}
    for (k_, j_, b_, _) in wires:
        e_, f_ = 5 * k_ - 5 + j_, 5 * k_ + j_
        bits_of[e_] = max(bits_of.get(e_, -1), b_)
        bits_of[f_] = max(bits_of.get(f_, -1), b_)
        pairs.setdefault((e_, f_), set()).add(b_)
    single = []
    for t in range(ncol):
        s = 1 if t not in bits_of else (2 if bits_of[t] == 0 else 4)
        for d in (s, -s):
            single.append((t, d))
    pl = [(e_, f_, sorted(bs)) for (e_, f_), bs in pairs.items()]
    dd = [(dx, dy) for dx in range(-3, 4) for dy in range(-3, 4) if dx and dy]
    V = decode_pi(np.mod(wq, 1 << brot) / scale, edges)
    e1, perm, D, lam = error1(V, H)
    for sweep in range(sweeps):
        moved = 0
        for j in range(N):
            lj = lam[perm[j]]
            base = D - lj * np.outer(V[j], V[j])
            for _ in range(inner):
                cands = []
                for (t, d) in single:
                    v = wq[j, t] + d
                    if lo <= v <= hi:
                        c = wq[j].copy()
                        c[t] = v
                        cands.append(c)
                for (e_, f_, bs) in pl:
                    we, wf = int(wq[j, e_]), int(wq[j, f_])
                    x0 = we ^ wf
                    for (dx, dy) in dd:
                        a_, b2 = we + dx, wf + dy
                        if not (lo <= a_ <= hi and lo <= b2 <= hi):
                            continue
                        x1 = a_ ^ b2
                        if all(((x1 >> b) & 1) == ((x0 >> b) & 1) for b in bs):
                            c = wq[j].copy()
                            c[e_], c[f_] = a_, b2
                            cands.append(c)
                C = np.array(cands, np.int64)
                Vc = decode_pi(np.mod(C, 1 << brot) / scale, edges)
                v0 = V[j]
                quad = np.einsum("ip,pq,iq->i", Vc, base, Vc)
                n4 = np.sum(Vc * Vc, axis=1) ** 2
                val = 2 * lj * quad + lj * lj * n4
                cur = 2 * lj * float(v0 @ base @ v0) + lj * lj * float(v0 @ v0) ** 2
                i = int(np.argmin(val))
                if val[i] >= cur - 1e-13:
                    break
                wq[j] = C[i]
                V[j] = Vc[i]
                moved += 1
            D = base + lj * np.outer(V[j], V[j])
        e_prev = e1
        e1, perm, D, lam = error1(V, H)
        if not moved or e_prev - e1 < 1e-7 * e1:
            break
    return wq, e1


def one_try_b(args):
    """(H, edges, brot, iters, sweeps, wires, t, seed, fine, lo, hi) -> (error_1, words in address order, free, forced)"""
    H, edges, brot, iters, sweeps, wires, t, seed, fine, lo, hi = args
    from chem_wr7_table7c import T7
    scale = float(1 << (brot - 1))
    lam_e, vec = np.linalg.eigh(H)
    if hi < (1 << (brot - 1)) - 1:      # signed range asked: start from the signed gauge (obsigned's), theta in [-pi/2, pi/2]
        from r48_chem_q7c_obsigned import signed_angles
        th0, _ = signed_angles(vec.T, edges)
    else:
        th0, _ = angles_and_mass(vec.T, edges)
    n = th0.shape[0]
    perm = np.arange(n) if t == 0 else np.random.default_rng(1000 * seed + t).permutation(n)
    st = fine_stages(STAGES, fine) if fine > 1 else STAGES
    wq, dyn = progressive_dynb(th0[perm], edges, lam_e[perm], H, scale, iters, wires, st, lo, hi)
    wq, _ = polishb(wq, edges, H, brot, max(sweeps, 1), wires, lo=lo, hi=hi)
    out = np.zeros_like(wq)
    out[dyn.mrow] = wq
    bad = check(out, np.arange(n), wires)
    assert bad == 0, f"{bad} relations broken on the one-body rows"
    out = np.mod(out, 1 << brot)         # the table's storage: 0 .. 2^brot - 1 (bits 0 and 1 unchanged)
    e = error1(T7.reconstruct_u_vectors(out / scale, edges), H)[0]
    return float(e), out, int(dyn.free_ok), int(dyn.forced_n)
