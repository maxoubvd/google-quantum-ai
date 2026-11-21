"""r49_chem_free_codemodel: the damage of low-bit CODES at 7 bits in the linear model of the staged rounding
(R49-CHEM-FREE). One group at a time, the model of tools/r48_chem_q7c_linmodel.py (copied): S = Gauss-Newton matrix of
the group's angles with its weights eliminated; the staged procedure freezes the entries in the order of S_kk (the
real tool's --rank sdyn), every free entry sitting at its conditional mean (exact Gaussian conditioning).

New here: a code on bit i of the words of a column, per 16-address block of the table: the 16 bits must be a codeword
of RM(d, 4). In the staged order a word's bit is FREE while its address is outside the span (in the code's monomial
space) of the addresses already frozen in its (block, column), and FORCED afterwards; a forced word goes to the
nearest integer with the forced low bits. So the stiff words choose and the soft words pay: sensitivity-weighted
decoding by construction. The addresses are the real table's (rows matched by overlap: out/cost_perm.npy of
r49_chem_free_costmap.py), blocks cut at multiples of 16; rows of other groups in a boundary block are ignored.

    python r49_chem_free_codemodel.py POINT.npz --edges TABLE.npz --perm out/cost_perm.npy --group 0
           --codes "none;0:1:15-74;0:1:0-74;0:1:15-74,1:2:15-74" [--rank sdiag]
Prints the staged damage (mHa^2) of every code set for the group.
"""
import argparse
import os
import sys
import time

import numpy as np
import scipy.linalg as sl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_round import LONG, N, Target, angles_and_mass  # noqa: E402
from r48_chem_q7c_absorb import tangents  # noqa: E402

STEP = np.pi / 64.0


def build_S(r, U, W, A, I, J, s, edges):
    R, B, C = W.shape
    u, V = tangents(U[r], edges)
    T_ = V.shape[1]
    p = (u[:, I] * u[:, J]) * s
    q = ((V[:, :, I] * u[:, None, J] + u[:, None, I] * V[:, :, J]) * s).reshape(B * T_, -1)
    g = (W[r] @ A[r]) * s
    ar = A[r] * s
    bk = np.repeat(np.arange(B), T_)
    Gqq = q @ q.T
    Gqg = q @ g.T
    Ggg = g @ g.T
    Gqp, Gqa = q @ p.T, q @ ar.T
    Ggp, Gga = g @ p.T, g @ ar.T
    Gpp, Gaa, Gap = p @ p.T, ar @ ar.T, ar @ p.T
    nw = B * C
    X = Gqg[:, bk]
    Haa = 2.0 * (Ggg[np.ix_(bk, bk)] * Gqq + X.T * X)
    del Gqq, X
    wb, wc = np.repeat(np.arange(B), C), np.tile(np.arange(C), B)
    Haw = 2.0 * (Gga[np.ix_(bk, wc)] * Gqp[:, wb] + Ggp[np.ix_(bk, wb)] * Gqa[:, wc])
    Hww = 2.0 * (Gaa[np.ix_(wc, wc)] * Gpp[np.ix_(wb, wb)] + Gap[np.ix_(wc, wb)] * Gap[np.ix_(wc, wb)].T)
    lam = 1e-10 * np.trace(Hww) / nw
    S = Haa - Haw @ sl.cho_solve(sl.cho_factor(Hww + lam * np.eye(nw)), Haw.T)
    return 0.5 * (S + S.T) * (STEP * STEP * 1e6)


class Coder:
    """Per (block, column, bit): the frozen addresses' monomial vectors in reduced form, with their bits."""

    def __init__(self, codes):
        self.codes = codes              # list of (bit, degree, set of columns)
        self.monos = {}
        for bit, deg, _ in codes:
            self.monos[deg] = [m for m in range(16) if bin(m).count("1") <= deg]
        self.basis = {}                 # (ci, block, col) -> list of (vec, val)

    def vec(self, deg, l):
        v = 0
        for i, m in enumerate(self.monos[deg]):
            if (l & m) == m:
                v |= 1 << i
        return v

    def forced(self, block, col, l):
        """{bit: value} forced for the word at low address l of (block, col); None entries = free."""
        out = {}
        for ci, (bit, deg, cols) in enumerate(self.codes):
            if col not in cols:
                continue
            v, val = self.vec(deg, l), 0
            for bv, bval in self.basis.get((ci, block, col), []):
                if v & (bv & -bv):
                    v ^= bv
                    val ^= bval
            out[bit] = (val if v == 0 else None, v, val)
        return out

    def commit(self, block, col, l, z, info):
        for ci, (bit, deg, cols) in enumerate(self.codes):
            if col not in cols:
                continue
            fv, v, val = info[bit]
            if fv is None:              # a new independent position: its bit becomes a basis row
                b = (int(z) >> bit) & 1
                self.basis.setdefault((ci, block, col), []).append((v, val ^ b))


class CosetCoder:
    """C2's weaker requirement for bit 0: per (block, slot) the columns of the slot lie in ONE coset of RM(1,4): the
    parity of the word at low address l in column number c of its slot is v_l + a_c0 + a_c.l (16 + 5 per column
    unknowns). First-come: a word is free while its equation is independent of those already frozen."""

    def __init__(self, cols):
        self.cols, self.c0 = cols, min(cols)
        self.basis = {}

    def forced(self, block, col, l):
        if col not in self.cols:
            return {}
        v = 1 << int(l)
        cidx = (int(col) - self.c0) // 5
        for i, m in enumerate((0, 1, 2, 4, 8)):
            if (int(l) & m) == m:
                v |= 1 << (16 + 5 * cidx + i)
        val = 0
        for bv, bval in self.basis.get((int(block), int(col) % 5), []):
            if v & (bv & -bv):
                v ^= bv
                val ^= bval
        return {0: (val if v == 0 else None, v, val)}

    def commit(self, block, col, l, z, info):
        fv, v, val = info[0]
        if fv is None:
            self.basis.setdefault((int(block), int(col) % 5), []).append((v, val ^ (int(z) & 1)))


CW1 = np.array([[(a0 ^ (bin(a & l).count("1") & 1)) for l in range(16)] for a0 in (0, 1) for a in range(16)], np.int64)


def staged_coded(S, x, fracs, sk, coder, block_of, low_of, col_of, weighted=False):
    """weighted: a word whose bit 0 is still FREE in its (block, column) takes the parity of the best codeword of
    RM(1,4) for the whole set: among the codewords that agree with the bits already frozen there, the one of least
    sum over the unfrozen words of (conditional stiffness 1 / P_jj) x (distance of the present mean to the nearest
    integer of that parity)^2. Plain: it takes its own nearest integer."""
    n = len(x)
    members = {}
    if weighted and coder is not None:
        for k in range(n):
            members.setdefault((int(block_of[k]), int(col_of[k])), []).append(k)
    frozen_bits = {}
    posarr = np.zeros(n, np.int64)
    P = sl.cho_solve(sl.cho_factor(S + 1e-9 * np.trace(S) / n * np.eye(n)), np.eye(n))
    mu = x.copy()
    free = np.arange(n)
    E, done, nforced, nflip = 0.0, 0, 0, 0
    for f in fracs:
        need = int(round(f * n)) - done
        if need <= 0:
            continue
        m = mu[free]
        score = sk[free]
        pick = np.argsort(-score, kind="stable")[:need]
        if weighted:
            posarr[:] = -1
            posarr[free] = np.arange(len(free))
            dP = np.diag(P)
        keep = np.setdiff1d(np.arange(len(free)), pick)
        zt = np.round(m[pick])
        if coder is not None:
            for ii, pk in enumerate(pick):
                k = free[pk]
                info = coder.forced(block_of[k], col_of[k], low_of[k])
                if info:
                    cons = {b: fv for b, (fv, _, _) in info.items() if fv is not None}
                    key = (int(block_of[k]), int(col_of[k]))
                    if weighted and 0 in info and info[0][0] is None:
                        fb = frozen_bits.get(key, {})
                        ok = np.ones(len(CW1), bool)
                        for l_, b_ in fb.items():
                            ok &= CW1[:, l_] == b_
                        cost_cw = np.zeros(len(CW1))
                        for j in members[key]:
                            pj = posarr[j]
                            if pj < 0 or int(low_of[j]) in fb:
                                continue
                            mj, wj = m[pj], 1.0 / max(dP[pj], 1e-300)
                            pen = [(2.0 * np.round((mj - par) / 2.0) + par - mj) ** 2 for par in (0, 1)]
                            cost_cw += wj * np.where(CW1[:, int(low_of[j])] == 0, pen[0], pen[1])
                        cost_cw[~ok] = np.inf
                        cons[0] = int(CW1[int(np.argmin(cost_cw)), int(low_of[k])])
                    if cons:
                        base = int(np.round(m[pk]))
                        best = None
                        for dz in (0, 1, -1, 2, -2, 3, -3, 4, -4):
                            zc = base + dz
                            if all(((zc >> b) & 1) == v for b, v in cons.items()):
                                if best is None or abs(zc - m[pk]) < abs(best - m[pk]):
                                    best = zc
                        nforced += 1
                        nflip += int(best != base)
                        zt[ii] = best
                    coder.commit(block_of[k], col_of[k], low_of[k], zt[ii], info)
                    if weighted and 0 in info:
                        frozen_bits.setdefault(key, {})[int(low_of[k])] = int(zt[ii]) & 1
        e = zt - m[pick]
        Pdd = P[np.ix_(pick, pick)]
        cf = sl.cho_factor(Pdd + 1e-14 * np.trace(Pdd) / len(pick) * np.eye(len(pick)))
        sol = sl.cho_solve(cf, e)
        E += float(e @ sol)
        Pud = P[np.ix_(keep, pick)]
        mu[free[keep]] += Pud @ sol
        P = P[np.ix_(keep, keep)] - Pud @ sl.cho_solve(cf, Pud.T)
        free = free[keep]
        done += len(pick)
    return E, nforced, nflip


def parse(spec):
    if spec == "none":
        return None
    out = []
    for part in spec.split(","):
        b, d, rng = part.split(":")
        lo, hi = (int(v) for v in rng.split("-"))
        out.append((int(b), int(d), set(range(lo, hi + 1))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("init")
    ap.add_argument("--edges", required=True)
    ap.add_argument("--perm", required=True)
    ap.add_argument("--group", type=int, default=0)
    ap.add_argument("--codes", default="none;0:1:15-74")
    a = ap.parse_args()
    t0 = time.time()
    z = np.load(a.init)
    U, W = np.array(z["U"], np.float64), np.array(z["W"], np.float64)
    R, B, C = W.shape
    U /= np.linalg.norm(U, axis=2, keepdims=True)
    edges = np.array(np.load(a.edges)["tree_edges"], np.int64)
    perm = np.load(a.perm)
    tgt = Target()
    I, J = tgt.i, tgt.j
    s = np.sqrt(np.where(I == J, 1.0, 2.0))
    A = np.einsum("rbc,rbp->rcp", W, U[:, :, I] * U[:, :, J])
    r = a.group
    S = build_S(r, U, W, A, I, J, s, edges)
    T_ = N - 1
    th, _ = angles_and_mass(U[r], edges)
    x = (th * 64.0).reshape(-1)
    addr = np.zeros(B, np.int64)                       # address of point row b
    addr[perm[r]] = 76 + r * B + np.arange(B)
    block_of = np.repeat(addr // 16, T_)
    low_of = np.repeat(addr % 16, T_)
    col_of = np.tile(np.arange(T_), B)
    sk = np.diag(S).copy()
    print(f"group {r}: S built ({time.time() - t0:.0f} s); {B * T_} words", flush=True)
    base = None
    for spec in a.codes.split(";"):
        wtd = spec.startswith("w")
        if spec.startswith("c"):                      # c0:1:15-74 = the coset form of the bit-0 code
            lo_, hi_ = (int(v) for v in spec.split(":")[2].split("-"))
            cd, coder = True, CosetCoder(set(range(lo_, hi_ + 1)))
        else:
            cd = parse(spec[1:] if wtd else spec)
            coder = Coder(cd) if cd else None
        E, nf, nfl = staged_coded(S, x, LONG, sk, coder, block_of, low_of, col_of, weighted=wtd)
        if base is None:
            base = E
        print(f"  group {r} codes [{spec}]: staged damage {E:.2f} mHa^2 ({E - base:+.2f}, x{E / base:.3f}); words with a forced bit {nf}, "
              f"moved off their nearest value {nfl} ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
