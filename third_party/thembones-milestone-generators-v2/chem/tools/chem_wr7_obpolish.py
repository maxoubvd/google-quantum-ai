"""chem_wr7_obpolish (WR7): exact coordinate descent on error_1 for the 76 one-body rows under bit-0 delta wires, with
the moves that keep every relation: one step on an unwired word, two steps on a wired word (bit 0 kept: the
designer's polish_coded of r49_chem_free_obround.py), AND PAIR MOVES: both words of a wire pair (e = 5k-5+j,
f = 5k+j) moved by one step each (four sign combinations): both bits 0 flip, their XOR (the relation) is kept.
Rows stay at their addresses. Same objective and acceptance rule as polish_coded.

    from chem_wr7_obpolish import polish_pairs
    wq, e1 = polish_pairs(wq, edges, H, brot, sweeps, wires)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_onebody import decode_pi, error1  # noqa: E402


def polish_pairs(wq, edges, H, brot, sweeps, wires, inner=25):
    scale, mod = float(1 << (brot - 1)), 1 << brot
    N = edges.shape[0] + 1
    ncol = wq.shape[1]
    wq = np.mod(np.array(wq, np.int64), mod)
    wired = set()
    parent = {}

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for (k_, j_, _, _) in wires:
        e_, f_ = 5 * k_ - 5 + j_, 5 * k_ + j_
        wired |= {e_, f_}
        parent.setdefault(e_, e_)
        parent.setdefault(f_, f_)
        parent[root(e_)] = root(f_)
    # WR7 Nov 17 (seat 2): CHAIN moves. Wires that share a column form a chain (k8-k9-k10 with k9 wired: four words,
    # three wires); a move keeps every relation iff it flips bit 0 of EVERY word of the chain: each word +-1 (2^n
    # sign patterns). For a chain of two words this is exactly the old pair move.
    comps = {}
    for c in wired:
        comps.setdefault(root(c), []).append(c)
    chains = [sorted(v) for v in comps.values()]
    assert max(len(v) for v in chains) <= 6, "chain longer than six words"
    step = np.array([2 if t in wired else 1 for t in range(ncol)], np.int64)
    # candidate move table: (n_moves, ncol) integer deltas
    moves = []
    for t in range(ncol):
        d = np.zeros(ncol, np.int64)
        d[t] = step[t]
        moves.append(d.copy())
        moves.append(-d)
    for ch in chains:
        for sgn in range(1 << len(ch)):
            d = np.zeros(ncol, np.int64)
            for i_, c_ in enumerate(ch):
                d[c_] = 1 if (sgn >> i_) & 1 else -1
            moves.append(d)
    M = np.array(moves, np.int64)
    V = decode_pi(wq / scale, edges)
    e1, perm, D, lam = error1(V, H)
    for sweep in range(sweeps):
        moved = 0
        for j in range(N):
            lj = lam[perm[j]]
            base = D - lj * np.outer(V[j], V[j])
            for _ in range(inner):
                v0 = V[j]
                cand = np.mod(wq[j][None, :] + M, mod)
                Vc = decode_pi(cand / scale, edges)
                quad = np.einsum("ip,pq,iq->i", Vc, base, Vc)
                n4 = np.sum(Vc * Vc, axis=1) ** 2
                val = 2 * lj * quad + lj * lj * n4
                cur = 2 * lj * float(v0 @ base @ v0) + lj * lj * float(v0 @ v0) ** 2
                i = int(np.argmin(val))
                if val[i] >= cur - 1e-13:
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
