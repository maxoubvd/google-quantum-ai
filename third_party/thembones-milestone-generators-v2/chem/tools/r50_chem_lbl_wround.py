"""r50_chem_lbl_wround: DELTA WIRES ON ADDRESS BITS in the linear model of the staged rounding (R50-CHEM-LBL,
construction 9). The relation of a wire k:j:b:m (bit b of word 5k+j XOR bit b of word 5k-5+j = parity(the row's
address AND m)) is imposed the way a code is: the first of the two words to be frozen is free, the second takes the bit
that the relation asks for (its nearest integer with that bit). The rows' addresses come from r50_chem_lbl_wcost.py
(--addr; default: the present order). build_S and staged_coded are FREE's, imported unchanged.

    python r50_chem_lbl_wround.py POINT.npz TREE.npz --group 2 --wires "14:0:0:0,..." [--addr ADDR.npz] [--jitter SEED:AMP]
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_round import LONG, N, Target, angles_and_mass  # noqa: E402
from r49_chem_free_codemodel import build_S, staged_coded  # noqa: E402


class Wires:
    def __init__(self, wires, addr_of_row):
        self.by_col = {}
        for (k, j, b, i) in wires:
            e, f = 5 * k - 5 + j, 5 * k + j
            assert e not in self.by_col and f not in self.by_col, "a word in two relations: use a chain-free layout"
            self.by_col.setdefault(e, []).append((f, b, i))
            self.by_col.setdefault(f, []).append((e, b, i))
        self.addr = addr_of_row
        self.z = {}

    def forced(self, row, col, low):
        out = {}
        for (p, b, i) in self.by_col.get(int(col), ()):
            zp = self.z.get((int(row), p))
            if zp is None:
                out.setdefault(b, (None, None, None))
            else:
                out[b] = (((zp >> b) & 1) ^ (bin(int(self.addr[int(row)]) & i).count("1") & 1), None, None)
        return out

    def commit(self, row, col, low, z, info):
        self.z[(int(row), int(col))] = int(z)


class DynWires:
    """THE ADDRESS DECIDED AT FREEZE TIME. A row's address is free until its relations are complete: when the second
    word of a pair is frozen it takes its own nearest word, and the parity it leaves becomes a known bit of the row's
    address, as long as the rows of the group can still be given distinct addresses (a bipartite matching, re-solved
    by one augmenting path); only when they cannot is the word forced."""

    def __init__(self, wires, addrs, nrows):
        self.pair = {}
        self.L = []
        for w, (k, j, b, m) in enumerate(wires):
            assert b == 0
            e, f = 5 * k - 5 + j, 5 * k + j
            assert e not in self.pair and f not in self.pair, "a word in two relations: use a chain-free layout"
            self.pair[e] = (f, w)
            self.pair[f] = (e, w)
            self.L.append(np.array([bin(int(x) & m).count("1") & 1 for x in addrs], np.int64))
        self.na = len(addrs)
        self.allowed = np.ones((nrows, self.na), bool)
        self.mrow = np.arange(nrows)                    # row -> address index
        self.maddr = np.full(self.na, -1, np.int64)
        self.maddr[:nrows] = np.arange(nrows)
        self.z = {}
        self.free_ok = self.forced_n = 0

    def _augment(self, row, allowed_row, apply):
        """row loses its address and must be re-matched inside allowed_row; True if a perfect matching exists."""
        mrow, maddr = self.mrow.copy(), self.maddr.copy()
        maddr[mrow[row]] = -1
        mrow[row] = -1
        allowed = self.allowed
        seen = np.zeros(self.na, bool)

        def dfs(r, al):
            for a_ in np.where(al & ~seen)[0]:
                if seen[a_]:
                    continue
                seen[a_] = True
                o = maddr[a_]
                if o < 0 or dfs(o, allowed[o]):
                    maddr[a_] = r
                    mrow[r] = a_
                    return True
            return False

        ok = dfs(row, allowed_row)
        if ok and apply:
            self.mrow, self.maddr = mrow, maddr
        return ok

    def _restrict(self, row, w, v, apply):
        al = self.allowed[row] & (self.L[w] == v)
        if al[self.mrow[row]]:
            ok = True
        else:
            ok = self._augment(row, al, apply)
        if ok and apply:
            self.allowed[row] = al
        return ok

    def forced(self, row, col, low):
        row, col = int(row), int(col)
        if col not in self.pair:
            return {}
        p, w = self.pair[col]
        zp = self.z.get((row, p))
        if zp is None:
            return {0: (None, None, None)}
        vcur = int(self.L[w][self.mrow[row]])
        if self._restrict(row, w, 1 - vcur, False):
            self.free_ok += 1
            return {0: (None, None, None)}
        self.forced_n += 1
        return {0: ((zp & 1) ^ vcur, None, None)}

    def commit(self, row, col, low, z, info):
        row, col = int(row), int(col)
        self.z[(row, col)] = int(z)
        p, w = self.pair[col]
        zp = self.z.get((row, p))
        if zp is not None:
            assert self._restrict(row, w, (int(z) & 1) ^ (zp & 1), True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("point")
    ap.add_argument("tree")
    ap.add_argument("--group", type=int, default=2)
    ap.add_argument("--wires", required=True)
    ap.add_argument("--addr", default="")
    ap.add_argument("--jitter", default="")
    ap.add_argument("--skip-untied", action="store_true")
    ap.add_argument("--dynamic", action="store_true", help="the rows' addresses decided at freeze time (DynWires)")
    a = ap.parse_args()
    t0 = time.time()
    z = np.load(a.point)
    U, W = np.array(z["U"], np.float64), np.array(z["W"], np.float64)
    R, B, C = W.shape
    U /= np.linalg.norm(U, axis=2, keepdims=True)
    tz = np.load(a.tree)
    edges = np.array(tz["tree_edges"] if "tree_edges" in tz.files else tz["edges"], np.int64)
    wires = [tuple(int(v) for v in s.split(":")) for s in a.wires.split(",")]
    r = a.group
    addr = np.load(a.addr)["addr"][r] if a.addr else 76 + 127 * r + np.arange(B)
    tgt = Target()
    I, J = tgt.i, tgt.j
    s = np.sqrt(np.where(I == J, 1.0, 2.0))
    A = np.einsum("rbc,rbp->rcp", W, U[:, :, I] * U[:, :, J])
    S = build_S(r, U, W, A, I, J, s, edges)
    T_ = N - 1
    th, _ = angles_and_mass(U[r], edges)
    x = (th * 64.0).reshape(-1)
    if a.jitter:
        sd, amp = a.jitter.split(":")
        x = x + np.random.default_rng(int(sd)).uniform(-float(amp), float(amp), size=x.shape)
    n = B * T_
    row_of = np.repeat(np.arange(B), T_)
    col_of = np.tile(np.arange(T_), B)
    low_of = np.zeros(n, np.int64)
    sk = np.diag(S).copy()
    Du = float("nan")
    if not a.skip_untied:
        Du, _, _ = staged_coded(S, x, LONG, sk, None, row_of, low_of, col_of)
        print(f"group {r}: D_u (no wires) {Du:.2f} mHa^2 ({time.time() - t0:.0f} s)", flush=True)
    if a.dynamic:
        dyn = DynWires(wires, 76 + 127 * r + np.arange(B), B)
        D, nf, nflip = staged_coded(S, x, LONG, sk, dyn, row_of, low_of, col_of)
        print(f"group {r}: {len(wires)} wires, addresses decided at freeze time: D {D:.2f} ({D - Du:+.2f}, x{D / Du:.4f}); second words left free "
              f"{dyn.free_ok}, forced {dyn.forced_n}, of them moved off their nearest word {nflip} ({time.time() - t0:.0f} s)", flush=True)
        return
    D, nf, nflip = staged_coded(S, x, LONG, sk, Wires(wires, addr), row_of, low_of, col_of)
    print(f"group {r}: {len(wires)} wires, rows {'assigned' if a.addr else 'in the present order'}: D {D:.2f} ({D - Du:+.2f}, x{D / Du:.4f}); "
          f"second words of a relation {nf}, moved off their nearest word {nflip} ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
