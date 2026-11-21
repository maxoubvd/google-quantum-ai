"""r50_chem_lbl_place: THE BIT-0 COSET CODE RM(1,4) PER 16-ADDRESS BLOCK WITH THE ROW'S BLOCK AND PLACE DECIDED AT
FREEZE TIME, in the linear model of the staged rounding (R50-CHEM-LBL, second hold, attack 2).

The model is R51-CHEM-BIT1's `staged_tool` (tools/r51_chem_bit1_codemodel.py: R49-CHEM-FREE's staged model with the
BUILDER'S decoder, the exact soft coset decoder `shared_coset_bits` of the rounding tool), copied; build_S imported.
New here (Placer): before every stage's decode the rows of the group are RE-PLACED on the group's addresses by an
alternating minimisation of the decoder's own cost: (1) decode every (block, chain) on the present places (exact);
(2) with the decoded bit patterns held, a linear assignment of rows to places (cost of row i at place p = the sum over
the coded columns of the cost of the bit the pattern asks there; a frozen word of the wrong parity is prohibitive);
repeat until the total stops falling. Frozen words never break: the present places always cost zero on them.

Spec = <geom>:<mode>:<lo>-<hi>[:<wires>]   (specs separated by ';')
  geom  asis    the group's addresses 76 + 127 r + (0..126): two partial blocks (the model sees the neighbours' rows
                as absent: R51's convention, the reference numbers)
        al      the group's rows on addresses 0..126: seven full blocks and one of 15 (what a block is in truth)
  mode  fixed   the deal of --order (CD7's)                      place  every row may move at every stage
        int     only rows dealt to full blocks move (asis)       once   placed at the first stage only, then fixed
        sw      place + annealed pair swaps with the exact cost of the blocks touched (--n0 at the first stage, --ns
                at every later one)                              sw1    the swaps at the first stage only, then fixed
  lo-hi the coded batches (chains by slot), or 'none'
  wires (attack 3) a file of bit-0 delta wires k:j:0:m (r50_chem_lbl_chain.py): columns joined by a wire are ONE
        effective column of the code, so "3-14" with wires on k4, k6, .., k14 = wires on the even loaders and the code
        on the odd ones; 'none' with wires = the wires alone (places by this tool's soft decoder, not DynChain)

    python r50_chem_lbl_place.py POINT.npz --edges TREE.npz --order order.npy --group 2 --specs "asis:fixed:5-14;asis:place:5-14"
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import scipy.linalg as sl
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_round import LONG, N, Target, angles_and_mass  # noqa: E402
from r51_chem_bit1_codemodel import build_S  # noqa: E402
import r51_chem_bit1_round as RD  # noqa: E402


_H = RD._shared_checks(11)
_HM = (_H * (1 << np.arange(_H.shape[0]))[:, None]).sum(axis=0)
_NS = 1 << _H.shape[0]
_PERM = [np.arange(_NS) ^ int(_HM[i]) for i in range(16)]


def block_cost(t0, t1, groups):
    """The decoder's cost of ONE block (exact: the DP of `shared_coset_bits` without the backtrack). t0, t1 (16, nc):
    the cost of bit 0 / 1 at every place; groups: the column index lists of the chains."""
    V = np.full((t0.shape[1], _NS), np.inf)
    V[:, 0] = 0.0
    for i in range(16):
        if not (t0[i].any() or t1[i].any()):            # an empty place: either bit is free
            V = np.minimum(V, V[:, _PERM[i]])
            continue
        V = np.minimum(V + t0[i][:, None], V[:, _PERM[i]] + t1[i][:, None])
    return float(sum(V[g].sum(axis=0).min() for g in groups))


class Placer:
    """Rows of ONE group against its places. pblk / ppos: block and low address of every place."""

    def __init__(self, chains, pblk, ppos, paddr=None, wires=()):
        """chains: the code's chains (lists of columns: one coset of RM(1,4) per block). wires (attack 3): (k, j, 0, m)
        = bit 0 of word 5k+j XOR bit 0 of word 5k-5+j = parity(address AND m) in every row. Columns joined by wires
        form a RUN = ONE effective column of the code (its bits are those of the run's first column XOR the parity of
        the cumulative mask at the row's address); the coset code then links the effective columns of a chain."""
        wcols = {5 * k + j for k, j, _, _ in wires} | {5 * k - 5 + j for k, j, _, _ in wires}
        self.cc = np.array(sorted({int(c) for ch in chains for c in ch} | wcols), np.int64)
        idx = {int(c): i for i, c in enumerate(self.cc)}
        root, cm = {int(c): int(c) for c in self.cc}, {int(c): 0 for c in self.cc}
        for k, j, _, m in sorted(wires):
            c0, c1 = 5 * k - 5 + j, 5 * k + j
            root[c1], cm[c1] = root[c0], cm[c0] ^ int(m)
        roots = sorted(set(root.values()))
        self.e_of = np.array([roots.index(root[int(c)]) for c in self.cc], np.int64)      # effective column of each column
        self.cmask = np.array([cm[int(c)] for c in self.cc], np.int64)
        ech = np.array([len(chains) + e for e in range(len(roots))], np.int64)             # a run outside the code: alone
        for h, ch in enumerate(chains):
            for c in ch:
                ech[self.e_of[idx[int(c)]]] = h
        self.chain = ech                                                                   # chain id of each effective column
        self.wired = bool(wires)
        self.pblk, self.ppos = pblk, ppos
        self.nb = int(pblk.max()) + 1
        self.trace, self.strace = [], []
        self.groups_e = [np.flatnonzero(ech == h) for h in np.unique(ech)]
        self.groups = [np.flatnonzero(ech[self.e_of] == h) for h in np.unique(ech)]       # columns of each chain
        if self.wired:
            x = paddr[:, None] & self.cmask[None, :]
            self.Wp = np.zeros_like(x)
            while np.any(x):
                self.Wp ^= x & 1
                x = x >> 1
            self.Emat = (self.e_of[:, None] == np.arange(len(roots))[None, :]).astype(np.float64)

    def eff(self, cost, place, rows=None):
        """the cost of bit 0 / 1 of every effective column for the rows at their places"""
        if not self.wired:
            return cost if rows is None else [cost[0][rows], cost[1][rows]]
        rows = np.arange(len(place)) if rows is None else rows
        w = self.Wp[place[rows]] == 0
        c0, c1 = cost[0][rows], cost[1][rows]
        return [np.where(w, c0, c1) @ self.Emat, np.where(w, c1, c0) @ self.Emat]

    def costs(self, t, c, frozen, val):
        tc, cw = t[:, self.cc], c[:, self.cc]
        cost = [cw * (RD.residue_round(tc, v, 2) - tc) ** 2 for v in (0, 1)]
        have = val[:, self.cc] & 1
        big = 1e4 * max(float(x.max()) for x in cost) + 1e-300
        fz = frozen[:, self.cc]
        return [np.where(fz, np.where(have == v, 0.0, big), cost[v]) for v in (0, 1)]

    def patterns(self, cost, place):
        ec = self.eff(cost, place)
        tab = [np.zeros((self.nb, 16, len(self.chain))) for _ in (0, 1)]
        for v in (0, 1):
            tab[v][self.pblk[place], self.ppos[place]] = ec[v]
        pat = RD.shared_coset_bits(tab[0], tab[1], self.cc, 11, chain=self.chain)[self.pblk, self.ppos]
        return (pat[:, self.e_of] ^ self.Wp) if self.wired else pat      # (places, coded columns)

    def step(self, t, c, frozen, val, place, movable, maxit=40):
        cost = self.costs(t, c, frozen, val)
        rows = np.flatnonzero(movable)
        ar = np.arange(len(place))
        first = None
        for it in range(maxit):
            patp = self.patterns(cost, place).astype(np.float64)
            # a row keeps its own free bit per chain (the coset's v): at a place it pays the cheaper of the chain's
            # pattern and its complement (without this the patterns are fitted to their own rows and nothing moves)
            M = np.zeros((len(place), len(place)))
            for g in self.groups:
                d = (cost[1] - cost[0])[:, g] @ patp[:, g].T
                M += np.minimum(cost[0][:, g].sum(axis=1)[:, None] + d, cost[1][:, g].sum(axis=1)[:, None] - d)
            cur = float(M[ar, place].sum())
            if first is None:
                first = cur
            pl = place[rows]
            ri, ci = linear_sum_assignment(M[np.ix_(rows, pl)])
            new = place.copy()
            new[rows[ri]] = pl[ci]
            tot = float(M[ar, new].sum())
            if not tot < cur - 1e-12 * abs(cur) or np.array_equal(new, place):
                break
            place = new
        self.trace.append((first, cur, it))
        return place

    def blk_cost(self, cost, b, place):
        t = [np.zeros((16, len(self.chain))) for _ in (0, 1)]
        rows = np.flatnonzero(self.pblk[place] == b)
        ec = self.eff(cost, place, rows)
        for v in (0, 1):
            t[v][self.ppos[place[rows]]] = ec[v]
        return block_cost(t[0], t[1], self.groups_e)

    def own(self, cost, place):
        """every row's cost at its own place under the decoded patterns (the rows' free bits kept)"""
        patp = self.patterns(cost, place).astype(np.float64)[place]
        out = np.zeros(len(place))
        for g in self.groups:
            d = ((cost[1] - cost[0])[:, g] * patp[:, g]).sum(axis=1)
            out += np.minimum(cost[0][:, g].sum(axis=1) + d, cost[1][:, g].sum(axis=1) - d)
        return out

    def swaps(self, t, c, frozen, val, place, movable, nprop, tfrac, rng, refresh=200):
        """Annealed pair swaps with the EXACT cost of the blocks touched (what the alternating assignment misses: the
        codewords re-fitted after the move). The first row of a pair is drawn by its present excess."""
        cost = self.costs(t, c, frozen, val)
        base = np.minimum(cost[0], cost[1]).sum(axis=1)
        eb = np.array([self.blk_cost(cost, b, place) for b in range(self.nb)])
        start = float(eb.sum() - base.sum())
        rows = np.flatnonzero(movable)
        acc, T0, prob = 0, 0.0, None
        best, bestp = float(eb.sum()), place.copy()
        for it in range(nprop):
            if it % refresh == 0:
                e = np.maximum(self.own(cost, place) - base, 0.0)[rows]
                tot = float(e.sum())
                if tot <= 1e-9 * (abs(start) + 1e-300):
                    break
                prob = e + 0.02 * tot / len(rows)
                prob /= prob.sum()
                if it == 0:
                    T0 = tfrac * tot / len(rows)
            T = T0 * (1.0 - it / nprop) ** 2 + 1e-300
            i, j = rows[rng.choice(len(rows), p=prob)], rows[rng.integers(len(rows))]
            if i == j:
                continue
            bi, bj = self.pblk[place[i]], self.pblk[place[j]]
            new = place.copy()
            new[i], new[j] = place[j], place[i]
            ni = self.blk_cost(cost, bi, new)
            nj = self.blk_cost(cost, bj, new) if bj != bi else ni
            d = ni - eb[bi] + ((nj - eb[bj]) if bj != bi else 0.0)
            if d < 0 or rng.random() < np.exp(-d / T):
                place, acc = new, acc + 1
                eb[bi], eb[bj] = ni, nj
                if eb.sum() < best:
                    best, bestp = float(eb.sum()), place.copy()
        self.strace.append((start, best - float(base.sum()), acc))
        return bestp


def violations(v, chains, blk_row, pos_row, wires=(), addr=None):
    CW = np.array([[(a0 ^ (bin(a_ & l).count("1") & 1)) for l in range(16)] for a0 in (0, 1) for a_ in range(16)])
    bad = 0
    for ch in chains:
        for c0, c1 in zip(ch[:-1], ch[1:]):
            for b_ in np.unique(blk_row):
                rr = np.flatnonzero(blk_row == b_)
                bits = (v[rr, c0] ^ v[rr, c1]) & 1
                bad += int(not (CW[:, pos_row[rr]] == bits[None, :]).all(axis=1).any())
    for k, j, _, m in wires:
        par = np.array([bin(int(a_) & int(m)).count("1") & 1 for a_ in addr])
        bad += int(np.sum(((v[:, 5 * k + j] ^ v[:, 5 * k - 5 + j]) & 1) != par))
    return bad


def staged_place(P0, x, fracs, sk, chains, pblk, ppos, place0, B, T_, wmode="S", mode="fixed", movable=None, dyn=None,
                 n0=8000, ns=400, seed=0, wires=(), paddr=None, cw=None, record=None):
    """R51's staged_tool with the places re-decided before every stage's decode (mode place / int / once).
    dyn: a StagedChain of r50_chem_lbl_chain (bit-0 delta wires, address decided at freeze time) for attack 3."""
    n = len(x)
    coded = bool(chains) or bool(wires)
    if coded:
        RD.HALF = 1 << 20
        if not wires:
            deg, chain = RD.layout_arrays([[0, 1, ch] for ch in chains], (T_, 8))
            RD.CHAIN = chain
            xdeg = np.zeros((15, 8), np.int64)
        placer = Placer(chains, pblk, ppos, paddr, wires)
    place = place0.copy()
    rng = np.random.default_rng(seed)
    if movable is None:
        movable = np.ones(B, bool)
    P = P0.copy()
    mu = x.copy()
    free = np.arange(n)
    frozen = np.zeros(n, bool)
    val = np.zeros(n, np.int64)
    E, done, nmoved, stage = 0.0, 0, 0, 0
    cP = np.zeros(n)
    for f in fracs:
        need = int(round(f * n)) - done
        if need <= 0:
            continue
        pick = np.argsort(-sk[free], kind="stable")[:need]
        keep = np.setdiff1d(np.arange(len(free)), pick)
        k_ = free[pick]
        if coded:
            t = np.where(frozen, val.astype(np.float64), mu)
            if wmode == "S":
                c = sk
            elif wmode == "D":          # the word's stiffness at its own freeze stage in the uncoded run (a dry run's record)
                c = cw
            else:
                cP[:] = 0.0
                cP[free] = 1.0 / np.maximum(np.diag(P), 1e-300)
                c = np.where(frozen, 1.0, cP)
            t2, c2 = t.reshape(B, T_), c.reshape(B, T_) + 1e-30
            f2, v2 = frozen.reshape(B, T_), val.reshape(B, T_)
            if mode in ("place", "int", "sw") or (mode in ("once", "sw1") and stage == 0):
                place = placer.step(t2, c2, f2, v2, place, movable)
            if mode in ("sw", "sw1") and (stage == 0 or mode == "sw") and (n0 if stage == 0 else ns) > 0:
                place = placer.swaps(t2, c2, f2, v2, place, movable, n0 if stage == 0 else ns, 0.3 if stage == 0 else 0.05, rng)
                if stage == 0:
                    place = placer.step(t2, c2, f2, v2, place, movable)
            if wires:                                   # my decode: the bits the patterns ask at the rows' places
                r_, m_ = np.zeros((B, T_), np.int64), np.ones(T_, np.int64)
                r_[:, placer.cc] = placer.patterns(placer.costs(t2, c2, f2, v2), place)[place]
                m_[placer.cc] = 2
            else:
                r_, m_ = RD.residues(t2, c2, pblk[place], ppos[place], deg, xdeg, f2, v2)
            zt = RD.residue_round(t2, r_, m_).reshape(-1)[k_]
        else:
            zt = np.round(mu[k_])
        nmoved += int(np.sum(zt != np.round(mu[k_])))
        e = zt - mu[k_]
        Pdd = P[np.ix_(pick, pick)]
        cf = sl.cho_factor(Pdd + 1e-14 * np.trace(Pdd) / len(pick) * np.eye(len(pick)))
        if record is not None:
            record[k_] = np.diag(sl.cho_solve(cf, np.eye(len(pick))))
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
        stage += 1
    bad = violations(val.reshape(B, T_), chains, pblk[place], ppos[place], wires, None if paddr is None else paddr[place]) if coded else 0
    info = dict(bad=bad, nmoved=nmoved, rows_moved=int((place != place0).sum()),
                rows_new_block=int((pblk[place] != pblk[place0]).sum()) if coded else 0)
    if coded and mode != "fixed":
        tr = placer.trace
        info["stage0_cost"] = [round(tr[0][0], 3), round(tr[0][1], 3), tr[0][2]]
        info["restages"] = int(sum(1 for a_, b_, _ in tr[1:] if b_ < a_ - 1e-12 * abs(a_)))
        if placer.strace:
            st = placer.strace
            info["stage0_swaps"] = [round(st[0][0], 3), round(st[0][1], 3), st[0][2]]
            info["later_swaps_accepted"] = int(sum(v[2] for v in st[1:]))
    return E, info, place, val.reshape(B, T_)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("init")
    ap.add_argument("--edges", required=True)
    ap.add_argument("--order", required=True, help="the rounding tool's order file (order[address - 76] = flat row)")
    ap.add_argument("--group", type=int, default=0)
    ap.add_argument("--specs", required=True)
    ap.add_argument("--wmode", default="S")
    ap.add_argument("--json", default=None)
    ap.add_argument("--base", type=float, default=None, help="the uncoded damage of the group (skips the 'none' run)")
    ap.add_argument("--n0", type=int, default=8000, help="mode sw: swap proposals at the first stage")
    ap.add_argument("--ns", type=int, default=400, help="mode sw: swap proposals at every later stage")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--drift", action="store_true")
    ap.add_argument("--save", default=None, help="prefix: the final places and words of every placed spec")
    a = ap.parse_args()
    t0 = time.time()
    z = np.load(a.init)
    U, W = np.array(z["U"], np.float64), np.array(z["W"], np.float64)
    R, B, C = W.shape
    U /= np.linalg.norm(U, axis=2, keepdims=True)
    ze = np.load(a.edges)
    edges = np.array(ze["tree_edges" if "tree_edges" in ze.files else "edges"], np.int64)
    tgt = Target()
    I, J = tgt.i, tgt.j
    s = np.sqrt(np.where(I == J, 1.0, 2.0))
    A = np.einsum("rbc,rbp->rcp", W, U[:, :, I] * U[:, :, J])
    r = a.group
    S = build_S(r, U, W, A, I, J, s, edges)
    T_ = N - 1
    th, _ = angles_and_mass(U[r], edges)
    x = (th * 64.0).reshape(-1)
    order = np.load(a.order)
    pos = np.empty(R * B, np.int64)
    pos[order] = np.arange(R * B)
    place0 = pos[r * B + np.arange(B)] - r * B              # place index of point row b inside its group
    assert np.all((place0 >= 0) & (place0 < B)), "the order moves rows between groups"
    sk = np.diag(S).copy()
    n = len(x)
    P0 = sl.cho_solve(sl.cho_factor(S + 1e-9 * np.trace(S) / n * np.eye(n)), np.eye(n))
    del S
    print(f"group {r}: S and P built ({time.time() - t0:.0f} s); {n} words", flush=True)
    res, base = [], a.base
    cw = None
    if a.wmode == "D":
        cw = np.zeros(n)
        E0, _, _, _ = staged_place(P0, x, LONG, sk, [], None, None, place0, B, T_, record=cw)
        base = E0 if base is None else base
        print(f"  dry run (uncoded) {E0:.2f}: the words' stiffness at their freeze stage recorded ({time.time() - t0:.0f} s)", flush=True)
    if a.drift:      # how well a word's parity at the START predicts the parity it is frozen with (uncoded run)
        E0, _, _, wd = staged_place(P0, x, LONG, sk, [], None, None, place0, B, T_)
        chg = ((np.round(x).astype(np.int64) ^ wd.reshape(-1)) & 1).reshape(B, T_)
        rk = np.empty(n, np.int64)
        rk[np.argsort(-sk, kind="stable")] = np.arange(n)
        rk = rk.reshape(B, T_)
        for lo_, hi_ in ((25, 75), (15, 75)):
            sel = np.zeros((B, T_), bool)
            sel[:, lo_:hi_] = True
            out = []
            for q0, q1 in ((0, .05), (.05, .1), (.1, .2), (.2, .4), (.4, .6), (.6, .8), (.8, 1.0001)):
                m_ = sel & (rk >= q0 * n) & (rk < q1 * n)
                out.append(f"{q0:.2f}-{min(q1, 1):.2f}: {chg[m_].mean():.3f} ({int(m_.sum())})")
            print(f"  DRIFT g{r} columns {lo_}-{hi_ - 1}: share of words frozen with the other parity than at the start, by freeze rank: " + "; ".join(out), flush=True)
    for spec in a.specs.split(";"):
        t1 = time.time()
        parts = spec.split(":")
        geom, mode, rng = parts[0], parts[1], parts[2]
        paddr = (76 + B * r + np.arange(B)) if geom == "asis" else np.arange(B)
        pblk, ppos = paddr // 16 - paddr.min() // 16, paddr % 16
        chains, wires = [], ()
        if rng != "none":
            for part in rng.split("+"):                     # "7-10+11-14": separate chains (no code between them)
                lo, hi = (int(v) for v in part.split("-"))
                chains += [[5 * k + j for k in range(lo, hi + 1)] for j in range(5)]
        if len(parts) > 3:                                  # wires: a file of k:j:0:m, comma separated
            wires = [tuple(int(v) for v in t.split(":")) for t in open(parts[3]).read().strip().split(",") if t]
        movable = None
        if mode == "int":
            full = np.bincount(pblk, minlength=pblk.max() + 1) == 16
            movable = full[pblk[place0]]
        E, info, place, words = staged_place(P0, x, LONG, sk, chains, pblk, ppos, place0, B, T_, a.wmode, mode, movable,
                                             n0=a.n0, ns=a.ns, seed=a.seed, wires=wires, paddr=paddr, cw=cw)
        if base is None:
            base = E
        res.append(dict(group=r, spec=spec, damage=E, base=base, **info))
        print(f"  g{r} [{spec}]: {E:.2f} ({E - base:+.2f}, x{E / base:.4f}); {json.dumps(info)} ({time.time() - t1:.0f} s)", flush=True)
        if a.json:
            json.dump(res, open(a.json, "w"), indent=1)
        if a.save and mode != "fixed" and (chains or wires):
            np.savez(f"{a.save}_{parts[0]}_{parts[1]}_{parts[2]}{'_w' if wires else ''}_g{r}.npz", place=place, place0=place0, words=words, paddr=paddr)


if __name__ == "__main__":
    main()
