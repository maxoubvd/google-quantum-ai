"""r51_chem_bit1_codemodel: the damage of code LAYOUTS (any coded bit, RM(1,4) or RM(2,4) cosets, explicit chains of
columns, two row deals) at 7 bits in the linear model of the staged rounding (R51-CHEM-BIT1).

The model is R49-CHEM-FREE's (tools/r49_chem_free_codemodel.py, itself the model of tools/r48_chem_q7c_linmodel.py):
one group at a time, S = the Gauss-Newton matrix of the group's angles with its weights eliminated; the staged
procedure freezes the words in the order of S_kk, every free word at its conditional mean (exact Gaussian
conditioning). build_S is copied unchanged.

New here (GenCoder): a code = (bit, degree d, chains). A CHAIN is a list of columns whose bit `bit` must lie in ONE
coset of RM(d,4) on every 16-address block (the loaders' requirement for the columns of a slot in consecutive
batches: a chain of m columns = m - 1 coded (loader, slot) pairs). Equation of the word at low address l in the
q-th column of its chain: bit = v_l + sum_mono a_{q,mono} mono(l) (16 + nmono per column unknowns). First-come in
the staged order: a word's coded bit is FREE while its equation is independent of those already frozen in its
(block, chain), FORCED afterwards; a word with forced bits goes to the nearest integer carrying them.
Deals: plain = the table's row order (--perm); heavy5 = inside every block the group's rows re-placed, strongest
first, on low addresses taken in affine-independence-first order (0, 1, 2, 4, 8, then the rest).

Shorthand of a spec (specs separated by ';'):  [h:]comp,comp,...   h: = heavy5 deal (hs: = heavy5 by S_kk row sums)
  comp = bit:deg:B<lo>-<hi>      batches lo..hi, chains by slot (the tool's --codes bit:11:lo-hi / bit:12:lo-hi)
       | bit:deg:C<n>x<q>        the n cheapest columns of --colcost, dealt by rank into q chains
       | bit:deg:L<c>+<c>+..|<c>+..   explicit chains

    python r51_chem_bit1_codemodel.py POINT.npz --edges TABLE.npz --perm cost_perm.npy --colcost colcost.npy
           --group 2 --specs "none;0:1:B7-14;0:1:B7-14,1:1:B7-14;h:0:1:C30x5,1:1:C30x5" --json out.json
"""
import argparse
import json
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


def monos(deg):
    return [m for m in range(16) if bin(m).count("1") <= deg]


class GenCoder:
    """codes: list of (bit, deg, chains); chains = list of lists of columns."""

    def __init__(self, codes):
        self.codes = codes
        self.where = {}                       # column -> list of (code index, chain index, position in chain)
        self.mon = {}
        for ci, (bit, deg, chains) in enumerate(codes):
            self.mon[ci] = monos(deg)
            for hi, ch in enumerate(chains):
                for q, c in enumerate(ch):
                    self.where.setdefault(int(c), []).append((ci, hi, q))
        self.basis = {}
        self.bits_of = {c: sorted({codes[ci][0] for ci, _, _ in v}) for c, v in self.where.items()}

    def forced(self, block, col, l):
        """list of (code index, key, bit, forced value or None, reduced vector, value)"""
        out = []
        for ci, hi, q in self.where.get(int(col), ()):
            mon = self.mon[ci]
            v = 1 << int(l)
            base = 16 + len(mon) * q
            for i, m in enumerate(mon):
                if (int(l) & m) == m:
                    v |= 1 << (base + i)
            val = 0
            key = (ci, int(block), hi)
            for bv, bval in self.basis.get(key, ()):
                if v & (bv & -bv):
                    v ^= bv
                    val ^= bval
            out.append((key, self.codes[ci][0], val if v == 0 else None, v, val))
        return out

    def commit(self, info, z):
        for key, bit, fv, v, val in info:
            if fv is None:
                self.basis.setdefault(key, []).append((v, val ^ ((int(z) >> bit) & 1)))


def staged(S, P0, x, fracs, sk, coder, block_of, low_of, col_of, dyn=None):
    """dyn (n,) array or None: filled with 1 / P_jj of every word at the stage it is frozen (its stiffness then)."""
    n = len(x)
    P = P0.copy()
    mu = x.copy()
    free = np.arange(n)
    E, done = 0.0, 0
    nforced, nmoved, sq_forced = 0, 0, 0.0
    for f in fracs:
        need = int(round(f * n)) - done
        if need <= 0:
            continue
        m = mu[free]
        pick = np.argsort(-sk[free], kind="stable")[:need]
        keep = np.setdiff1d(np.arange(len(free)), pick)
        zt = np.round(m[pick])
        if coder is not None:
            for ii, pk in enumerate(pick):
                k = free[pk]
                if int(col_of[k]) not in coder.where:
                    continue
                info = coder.forced(block_of[k], col_of[k], low_of[k])
                cons = {}
                ok = True
                for key, bit, fv, v, val in info:
                    if fv is not None:
                        if bit in cons and cons[bit] != fv:
                            ok = False
                        cons[bit] = fv
                assert ok, "two codes force one bit two ways"
                if cons:
                    span = 1 << (max(cons) + 1)
                    base = int(np.round(m[pk]))
                    best = None
                    for dz in range(-span, span + 1):
                        zc = base + dz
                        if all(((zc >> b) & 1) == v for b, v in cons.items()):
                            if best is None or abs(zc - m[pk]) < abs(best - m[pk]):
                                best = zc
                    nforced += 1
                    nmoved += int(best != base)
                    sq_forced += float((best - m[pk]) ** 2)
                    zt[ii] = best
                coder.commit(info, zt[ii])
        e = zt - m[pick]
        Pdd = P[np.ix_(pick, pick)]
        cf = sl.cho_factor(Pdd + 1e-14 * np.trace(Pdd) / len(pick) * np.eye(len(pick)))
        if dyn is not None:      # the price of a unit error on this word alone, the stage's other words frozen too
            dyn[free[pick]] = np.diag(sl.cho_solve(cf, np.eye(len(pick))))
        sol = sl.cho_solve(cf, e)
        E += float(e @ sol)
        Pud = P[np.ix_(keep, pick)]
        mu[free[keep]] += Pud @ sol
        P = P[np.ix_(keep, keep)] - Pud @ sl.cho_solve(cf, Pud.T)
        free = free[keep]
        done += len(pick)
    return E, nforced, nmoved, sq_forced


def staged_tool(P0, x, fracs, sk, codes, addr, B, T_, wmode):
    """The staged procedure with the picked words set by the BUILDER'S decoder (tools/r51_chem_bit1_round.py residues
    with the layout's chains: per stage the cheapest residues of all words given the frozen ones), as CD7's kill test
    (a) did for the uniform code (tools/chem_cd7_codemodel.py staged_tool, copied). wmode S: weights S_kk; P: 1 / P_jj."""
    import r51_chem_bit1_round as RD
    RD.HALF = 1 << 20                                  # the model's angles are not gauged into the signed range
    deg, chain = RD.layout_arrays([[b_, 1, ch] for b_, d_, chains in codes for ch in chains], (T_, 8))
    RD.CHAIN = chain
    xdeg = np.zeros((15, 8), np.int64)
    blk_row = addr // 16
    blk_row = blk_row - blk_row.min()
    pos_row = addr % 16
    n = len(x)
    P = P0.copy()
    mu = x.copy()
    free = np.arange(n)
    frozen = np.zeros(n, bool)
    val = np.zeros(n, np.int64)
    E, done, nmoved = 0.0, 0, 0
    cP = np.zeros(n)
    for f in fracs:
        need = int(round(f * n)) - done
        if need <= 0:
            continue
        pick = np.argsort(-sk[free], kind="stable")[:need]
        keep = np.setdiff1d(np.arange(len(free)), pick)
        t = np.where(frozen, val.astype(np.float64), mu)
        if wmode == "S":
            c = sk
        else:
            cP[:] = 0.0
            cP[free] = 1.0 / np.maximum(np.diag(P), 1e-300)
            c = np.where(frozen, 1.0, cP)
        r_, m_ = RD.residues(t.reshape(B, T_), c.reshape(B, T_) + 1e-30, blk_row, pos_row, deg, xdeg,
                             frozen.reshape(B, T_), val.reshape(B, T_))
        allowed = RD.residue_round(t.reshape(B, T_), r_, m_).reshape(-1)
        k_ = free[pick]
        zt = allowed[k_]
        nmoved += int(np.sum(zt != np.round(mu[k_])))
        e = zt - mu[k_]
        Pdd = P[np.ix_(pick, pick)]
        cf = sl.cho_factor(Pdd + 1e-14 * np.trace(Pdd) / len(pick) * np.eye(len(pick)))
        sol = sl.cho_solve(cf, e)
        E += float(e @ sol)
        Pud = P[np.ix_(keep, pick)]
        mu[free[keep]] += Pud @ sol
        P = P[np.ix_(keep, keep)] - Pud @ sl.cho_solve(cf, Pud.T)
        frozen[k_] = True
        val[k_] = zt.astype(np.int64)
        mu[k_] = zt
        free = free[keep]
        done += len(pick)
    # the layout must hold on the finished words (rows of this group only)
    v = val.reshape(B, T_)
    CW = np.array([[(a0 ^ (bin(a_ & l).count("1") & 1)) for l in range(16)] for a0 in (0, 1) for a_ in range(16)])
    bad = 0
    for bit, d_, chains in codes:
        for ch in chains:
            for c0, c1 in zip(ch[:-1], ch[1:]):
                for b_ in np.unique(blk_row):
                    rr = np.flatnonzero(blk_row == b_)
                    bits = ((v[rr, c0] ^ v[rr, c1]) >> bit) & 1
                    bad += int(not (CW[:, pos_row[rr]] == bits[None, :]).all(axis=1).any())
    return E, bad, nmoved, 0.0


def indep_first(lows):
    """the low addresses `lows` in affine-independence-first order (greedy, ascending)."""
    lows = sorted(int(v) for v in lows)
    chosen, rest, basis = [], [], []

    def red(v):
        for b in basis:
            if v & (b & -b):
                v ^= b
        return v
    for l in lows:
        v = red((1 << 4) | l)                 # affine: the constant coordinate + the four bits
        if v:
            basis.append(v)
            chosen.append(l)
        else:
            rest.append(l)
    return chosen + rest


def heavy_deal(addr, strength):
    """addr (B,) addresses of the group's rows; inside every block the rows are re-placed strongest first on the
    block's addresses (of this group) in independence-first order."""
    new = addr.copy()
    for b in np.unique(addr // 16):
        rows = np.flatnonzero(addr // 16 == b)
        order = rows[np.argsort(-strength[rows], kind="stable")]
        lows = indep_first(addr[rows] % 16)
        new[order] = 16 * b + np.array(lows, np.int64)
    return new


def parse_spec(spec, colcost):
    deal = "plain"
    if spec.startswith("tp:"):                # the builder's decoder, weights 1 / P_jj (plain deal)
        deal, spec = "toolP", spec[3:]
    elif spec.startswith("t:"):               # the builder's decoder, weights S_kk (plain deal)
        deal, spec = "toolS", spec[2:]
    elif spec.startswith("hs:"):
        deal, spec = "heavy5s", spec[3:]
    elif spec.startswith("h:"):
        deal, spec = "heavy5", spec[2:]
    if spec == "none":
        return deal, []
    rank = np.argsort(colcost, kind="stable")
    codes = []
    for comp in spec.split(","):
        bit, deg, what = comp.split(":")
        if what[0] == "B":
            lo, hi = (int(v) for v in what[1:].split("-"))
            chains = [[5 * k + j for k in range(lo, hi + 1)] for j in range(5)]
        elif what[0] == "C":
            n_, q_ = (int(v) for v in what[1:].split("x"))
            chains = [[int(c) for c in rank[i:n_:q_]] for i in range(q_)]
        elif what[0] == "G":                 # G<m>@<lo>-<hi>: batches lo..hi cut into runs of m, chains by slot
            m_, rng = what[1:].split("@")
            lo, hi = (int(v) for v in rng.split("-"))
            runs = [list(range(b, min(b + int(m_), hi + 1))) for b in range(lo, hi + 1, int(m_))]
            chains = [[5 * k + j for k in run] for run in runs if len(run) > 1 for j in range(5)]
        elif what[0] == "L":
            chains = [[int(c) for c in ch.split("+")] for ch in what[1:].split("|")]
        else:
            raise ValueError(comp)
        codes.append((int(bit), int(deg), chains))
    return deal, codes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("init")
    ap.add_argument("--edges", required=True)
    ap.add_argument("--perm", default=None)
    ap.add_argument("--order", default=None, help="instead of --perm: the rounding tool's order file (order[address - 76] = flat row)")
    ap.add_argument("--colcost", required=True)
    ap.add_argument("--group", type=int, default=0)
    ap.add_argument("--specs", required=True, help="specs separated by ';' or @file (one per line)")
    ap.add_argument("--json", default=None)
    ap.add_argument("--dump", default=None, help="npz: per word the stiffness at its freeze stage in the uncoded run, S_kk, addresses")
    a = ap.parse_args()
    t0 = time.time()
    z = np.load(a.init)
    U, W = np.array(z["U"], np.float64), np.array(z["W"], np.float64)
    R, B, C = W.shape
    U /= np.linalg.norm(U, axis=2, keepdims=True)
    ze = np.load(a.edges)
    edges = np.array(ze["tree_edges" if "tree_edges" in ze.files else "edges"], np.int64)
    perm = np.load(a.perm) if a.perm else None
    colcost = np.load(a.colcost)
    tgt = Target()
    I, J = tgt.i, tgt.j
    s = np.sqrt(np.where(I == J, 1.0, 2.0))
    A = np.einsum("rbc,rbp->rcp", W, U[:, :, I] * U[:, :, J])
    r = a.group
    S = build_S(r, U, W, A, I, J, s, edges)
    T_ = N - 1
    th, _ = angles_and_mass(U[r], edges)
    x = (th * 64.0).reshape(-1)
    addr0 = np.zeros(B, np.int64)                       # address of point row b
    if a.order:                                         # the rounding tool's order file: order[address - 76] = flat row
        order = np.load(a.order)
        pos = np.empty(R * B, np.int64)
        pos[order] = np.arange(R * B)
        addr0 = 76 + pos[r * B + np.arange(B)]
        assert np.all((addr0 - 76) // B == r), "the order moves rows between groups"
    else:
        addr0[perm[r]] = 76 + r * B + np.arange(B)
    col_of = np.tile(np.arange(T_), B)
    sk = np.diag(S).copy()
    n = len(x)
    P0 = sl.cho_solve(sl.cho_factor(S + 1e-9 * np.trace(S) / n * np.eye(n)), np.eye(n))
    del S
    w2 = np.sum(W[r] ** 2, axis=1)
    skrow = sk.reshape(B, T_).sum(axis=1)
    print(f"group {r}: S and P built ({time.time() - t0:.0f} s); {n} words", flush=True)
    specs = open(a.specs[1:]).read().split() if a.specs.startswith("@") else a.specs.split(";")
    res, base = [], None
    for spec in specs:
        t1 = time.time()
        spec0 = spec
        deal, codes = parse_spec(spec, colcost)
        if deal in ("toolS", "toolP"):
            t1 = time.time()
            E, nf, nm, sq = staged_tool(P0, x, LONG, sk, codes, addr0, B, T_, deal[-1])
            links = sum(len(ch) - 1 for _, _, chains in codes for ch in chains)
            res.append(dict(group=r, spec=spec0, damage=E, base=base, forced=nf, moved=nm, sq_forced=sq, links=links))
            print(f"  g{r} [{spec0[:70]}]: {E:.2f} ({E - base:+.2f}, x{E / base:.4f}); links {links}; violations {nf}, moved {nm} "
                  f"({time.time() - t1:.0f} s)", flush=True)
            if a.json:
                json.dump(res, open(a.json, "w"), indent=1)
            continue
        addr = addr0 if deal == "plain" else heavy_deal(addr0, w2 if deal == "heavy5" else skrow)
        block_of, low_of = np.repeat(addr // 16, T_), np.repeat(addr % 16, T_)
        coder = GenCoder(codes) if codes else None
        dyn = np.zeros(n) if (a.dump and not codes) else None
        E, nf, nm, sq = staged(None, P0, x, LONG, sk, coder, block_of, low_of, col_of, dyn=dyn)
        if dyn is not None:
            np.savez(a.dump, dyn=dyn.reshape(B, T_), sk=sk.reshape(B, T_), addr=addr0, w2=w2, group=r, damage=E)
        if base is None:
            base = E
        links = sum(len(ch) - 1 for _, _, chains in codes for ch in chains)
        res.append(dict(group=r, spec=spec, damage=E, base=base, forced=nf, moved=nm, sq_forced=sq, links=links))
        print(f"  g{r} [{spec}]: {E:.2f} ({E - base:+.2f}, x{E / base:.4f}); links {links}; forced {nf}, moved {nm}, "
              f"sum sq of forced moves {sq:.0f} ({time.time() - t1:.0f} s)", flush=True)
        if a.json:
            json.dump(res, open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
