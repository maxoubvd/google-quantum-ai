"""chem_cd7_oblay (CD7, 17:10Z): the ONE-BODY rows under an R51-CHEM-BIT1 LAYOUT (coset chains on bits 0 and 1) for my
finishing chain (tools/chem_cd7_table7.py --spec). A copy of R51's progressive_layout (tools/r51_chem_bit1_obround.py:
26 stages, impact order, the rounding tool's exact residues per stage, L-BFGS refit of the free angles on error_1^2)
with three changes:
  (1) rows 76-79 (two-body, block 4, low addresses 12-15) enter the decoder as FROZEN pseudo-rows: their bits are
      fixed (moving them broke cd7_r1's first chain), block 4's one-body rows complete every chain around them;
  (2) every word stays in the SIGNED RANGE -32..31 (the decoder's range set to 32: coded words move by 2^(coded bits),
      uncoded words are clamped; the polish rejects moves out of it): obsigned then has nothing to move, and its
      negations, which flip bit 1 of an odd word, never touch a coded bit 1;
  (3) placements as my chain's (light4 / heavy5), chosen by the caller.
"""
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r51_chem_bit1_round as RD  # noqa: E402
from r48_chem_q7c_onebody import STAGES, angles_and_mass, decode_jac, decode_pi, error1  # noqa: E402
from r51_chem_bit1_codemodel import parse_spec  # noqa: E402

CW = np.array([[(a0 ^ (bin(a & l).count("1") & 1)) for l in range(16)] for a0 in (0, 1) for a in range(16)], np.int64)
LO, HI = -32, 31


def setup(spec_text, ncols=75):
    """Layout -> (codes, deg, step, coded columns); sets the tool's chains (RD.CHAIN) and range (RD.HALF = 32)."""
    _, codes = parse_spec(spec_text, np.arange(ncols))
    deg, chain = RD.layout_arrays([[b_, 1, ch] for b_, d_, chains in codes for ch in chains], (ncols, 8))
    RD.CHAIN = chain
    RD.HALF = 32
    step = np.ones(ncols, np.int64)
    coded = set()
    for bit, _, chains in codes:
        for ch in chains:
            for c in ch:
                step[c] = max(step[c], 1 << (bit + 1))
                coded.add(int(c))
    return codes, deg, step, coded


def signed(w):
    w = np.mod(np.asarray(w, np.int64), 128)
    return np.where(w >= 64, w - 128, w)


def progressive_layout_fixed(th0, edges, lam_at, H, scale, iters, deg, w76):
    n, nc = th0.shape
    theta = th0.copy()
    frozen = np.zeros(theta.shape, bool)
    qv = np.zeros(theta.shape)
    val = np.zeros(theta.shape, np.int64)
    lam2 = lam_at[:, None] ** 2
    blk_e = np.concatenate([np.arange(n) // 16, np.full(4, 4)])
    pos_e = np.concatenate([np.arange(n) % 16, np.arange(12, 16)])
    xdeg = np.zeros((15, 8), np.int64)
    w76s = signed(w76)
    moved = 0

    def fun(x, mask, base):
        t = base.copy()
        t[mask] = x
        V, J = decode_jac(t, edges)
        D = (V.T * lam_at) @ V - H
        gV = 4.0 * lam_at[:, None] * (V @ D)
        return float(np.sum(D * D)), np.einsum("jn,jkn->jk", gV, J)[mask]

    for frac in STAGES:
        need = int(round(frac * theta.size)) - int(frozen.sum())
        if need > 0:
            now = np.where(frozen, qv, theta)
            mass = angles_and_mass(decode_pi(now, edges), edges)[1]
            q = np.round(now * scale) / scale
            impact = lam2 * mass * (q - now) ** 2
            impact[frozen] = -np.inf
            pick = np.argsort(-impact, axis=None)[:need]
            t_e = np.vstack([now * scale, w76s.astype(np.float64)])
            c_e = np.vstack([lam2 * mass + 1e-30, np.ones((4, nc))])
            f_e = np.vstack([frozen, np.ones((4, nc), bool)])
            v_e = np.vstack([val, w76s])
            r_, m_ = RD.residues(t_e, c_e, blk_e, pos_e, deg, xdeg, f_e, v_e)
            for fl in pick:
                j, t = divmod(int(fl), nc)
                xs = now[j, t] * scale
                z = int(RD.into_range(RD.residue_round(xs, r_[j, t], m_[t]), m_[t]))
                moved += int(z != int(np.round(xs)))
                frozen[j, t] = True
                qv[j, t] = z / scale
                val[j, t] = z
        free = ~frozen
        if free.any():
            base = np.where(frozen, qv, theta)
            r = minimize(fun, theta[free], args=(free, base), jac=True, method="L-BFGS-B",
                         options={"maxiter": iters, "maxcor": 20, "ftol": 1e-14, "gtol": 1e-10})
            theta = base.copy()
            theta[free] = r.x
        else:
            theta = qv.copy()
    return np.round(qv * scale).astype(np.int64), moved


def polish_range(wq, edges, H, brot, sweeps, step, inner=25):
    """R51's exact coordinate descent on error_1 (step 2^(coded bits) per column), moves out of -32..31 rejected."""
    scale = float(1 << (brot - 1))
    N = edges.shape[0] + 1
    ncol = wq.shape[1]
    wq = np.array(wq, np.int64)
    V = decode_pi(wq / scale, edges)
    e1, perm, D, lam = error1(V, H)
    ar = np.arange(ncol)
    for sweep in range(sweeps):
        moved = 0
        for j in range(N):
            lj = lam[perm[j]]
            base = D - lj * np.outer(V[j], V[j])
            for _ in range(inner):
                v0 = V[j]
                cand = np.repeat(wq[j][None, :], 2 * ncol, axis=0)
                cand[ar, ar] += step
                cand[ncol + ar, ar] -= step
                ok = np.concatenate([wq[j] + step <= HI, wq[j] - step >= LO])
                Vc = decode_pi(cand / scale, edges)
                quad = np.einsum("ip,pq,iq->i", Vc, base, Vc)
                n4 = np.sum(Vc * Vc, axis=1) ** 2
                v_ = np.where(ok, 2 * lj * quad + lj * lj * n4, np.inf)
                cur = 2 * lj * float(v0 @ base @ v0) + lj * lj * float(v0 @ v0) ** 2
                i = int(np.argmin(v_))
                if v_[i] >= cur - 1e-13:
                    break
                wq[j] = cand[i]
                V[j] = Vc[i]
                moved += 1
            D = base + lj * np.outer(V[j], V[j])
        e_prev = e1
        e1, perm, D, lam = error1(V, H)
        if not moved or e_prev - e1 < 1e-7 * e1:
            break
    return wq, e1


def violations_full(wq, w76, codes):
    """Chain links not affine on blocks 0-4, block 4 with rows 76-79 at low 12-15."""
    w = np.vstack([np.mod(wq, 128), np.mod(np.asarray(w76, np.int64), 128)])
    bad = 0
    for bit, _, chains in codes:
        for ch in chains:
            for c0, c1 in zip(ch[:-1], ch[1:]):
                for blk in range(5):
                    rows = np.arange(16 * blk, 16 * blk + 16) if blk < 4 else np.arange(64, 80)
                    d = ((w[rows, c0] ^ w[rows, c1]) >> bit) & 1
                    bad += int(not (CW == d[None, :]).all(axis=1).any())
    return bad
