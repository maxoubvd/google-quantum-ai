"""r50_chem_lbl_hook: THE ROWS' PLACES DECIDED AT FREEZE TIME, AS A HOOK FOR A STAGED ROUNDING OF THE WHOLE TABLE
(R50-CHEM-LBL, second hold): the bit-0 coset code RM(1,4) per 16-address block (CD7's code 11) and / or bit-0 delta
wires (k:j:0:m), all groups at once, a row moving only inside its group. It is the Placer of r50_chem_lbl_place.py
(tested there in the staged model, one group at a time) wrapped for a tool whose rows stand in address order at the
start (address = n_one + row).

    fp = FreezePlacer(R, B, chains, wires)            # chains: lists of columns (one coset per block); wires: (k, j, 0, m)
    at every freeze, BEFORE the stage's words are set:
        r, m = fp.residues(t, c, frozen, value)       # t (rows, 75) in steps, c the weights, value the frozen words
        allowed = r + m * np.round((t - r) / m)       # the tool's residue_round; the stage's words take `allowed`
      (code only: `blk, pos = fp.addr() // 16, fp.addr() % 16` and the tool's own residues() give the same bits)
    at the end:  order_out = fp.order(order_in)        # order[address - n_one] = row of INIT; write words AND weights so
    a resumed run: fp.place = saved place (row -> address - n_one), nothing else is kept between calls.

--selftest: a mock staged freeze on random angles (3 groups): no relation broken, the decoder's cost against the fixed
order.
"""
import argparse
import os
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r51_chem_bit1_round as RD  # noqa: E402
from r50_chem_lbl_place import Placer, violations  # noqa: E402


class FreezePlacer:
    def __init__(self, R, B, chains, wires=(), n_one=76, ncols=75, fixed=False):
        self.R, self.B, self.n_one, self.ncols, self.fixed = R, B, n_one, ncols, fixed
        paddr = n_one + np.arange(R * B)
        self.paddr = paddr
        self.pl = Placer(chains, paddr // 16 - n_one // 16, paddr % 16, paddr, wires)
        self.group = np.arange(R * B) // B                 # of a row; a place p belongs to group p // B
        self.place = np.arange(R * B)
        self.moves = 0

    def addr(self):
        return self.paddr[self.place]

    def step(self, cost, maxit=40):
        pl, place = self.pl, self.place
        for _ in range(maxit):
            patp = pl.patterns(cost, place).astype(np.float64)
            new, gain = place.copy(), 0.0
            for g in range(self.R):
                rows = np.flatnonzero(self.group == g)
                pls = place[rows]
                M = np.zeros((len(rows), len(rows)))
                for grp in pl.groups:
                    d = (cost[1] - cost[0])[np.ix_(rows, grp)] @ patp[np.ix_(pls, grp)].T
                    M += np.minimum(cost[0][np.ix_(rows, grp)].sum(axis=1)[:, None] + d,
                                    cost[1][np.ix_(rows, grp)].sum(axis=1)[:, None] - d)
                ri, ci = linear_sum_assignment(M)
                cur, tot = float(np.trace(M)), float(M[ri, ci].sum())
                if tot < cur - 1e-12 * abs(cur):
                    new[rows[ri]] = pls[ci]
                    gain += cur - tot
            if gain <= 0.0:
                break
            self.moves += int((new != place).sum())
            place = new
        self.place = place

    def residues(self, t, c, frozen, value):
        cost = self.pl.costs(t, c, frozen, value)
        if not self.fixed:
            self.step(cost)
        r = np.zeros(t.shape, np.int64)
        m = np.ones(t.shape[1], np.int64)
        r[:, self.pl.cc] = self.pl.patterns(cost, self.place)[self.place]
        m[self.pl.cc] = 2
        return r, m

    def order(self, order_in):
        out = np.empty_like(order_in)
        out[self.place] = order_in
        return out


def selftest(seed=0):
    rng = np.random.default_rng(seed)
    R, B, T_ = 3, 127, 75
    RD.HALF = 1 << 20
    chains = [[5 * k + j for k in range(9, 15)] for j in range(5)]
    wires = [(k, j, 0, int(rng.integers(1, 128))) for k in (10, 12, 14) for j in range(5)]
    for name, ch, ws in (("code 9-14", chains, ()), ("wires k10, k12, k14", [], wires), ("wires + code on k11, k13", chains, wires)):
        out = []
        for fixed in (True, False):
            r_ = np.random.default_rng(seed + 1)
            x = r_.uniform(-30, 30, (R * B, T_))
            c = np.exp(r_.normal(0, 2.5, (R * B, T_))) * np.exp(r_.normal(0, 1.5, (R * B, 1)))
            fp = FreezePlacer(R, B, ch, ws, fixed=fixed)
            frozen = np.zeros(x.shape, bool)
            val = np.zeros(x.shape, np.int64)
            rank = np.argsort(-c, axis=None)
            cost = 0.0
            for st in np.array_split(rank, 12):
                t = np.where(frozen, val, x + (0.0 if not frozen.any() else 0.05 * r_.normal(size=x.shape)))
                r, m = fp.residues(t, c, frozen, val)
                allowed = (r + m * np.round((t - r) / m)).astype(np.int64)
                i, j = np.unravel_index(st, x.shape)
                val[i, j] = allowed[i, j]
                cost += float((c[i, j] * ((allowed[i, j] - t[i, j]) ** 2 - (np.round(t[i, j]) - t[i, j]) ** 2)).sum())
                frozen[i, j] = True
            a = fp.addr()
            bad = violations(val, ch, a // 16, a % 16, ws, a)
            od = fp.order(np.arange(R * B))
            assert bad == 0 and sorted(od.tolist()) == list(range(R * B)) and np.all(od // B == np.arange(R * B) // B)
            out.append(cost)
        print(f"selftest {name}: relations broken 0; rows kept in their groups; excess cost fixed order {out[0]:.1f}, "
              f"places at freeze time {out[1]:.1f} (x{out[1] / out[0]:.3f})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
