"""r50_chem_lbl_wtable: DELTA WIRES ON ADDRESS BITS IMPOSED ON A FINISHED TABLE (R50-CHEM-LBL, construction 9): the
two-body rows. For a wire k:j:0:m every row must satisfy  bit 0 of word 5k+j XOR bit 0 of word 5k-5+j = parity(address
AND m). On a rounded table the rows of each group are (1) ASSIGNED to the group's 127 addresses (the verifier lets
the vectors of a group permute with their weights) and (2) where the assigned address still asks for the other
parity, ONE word of the pair is moved by one step: the move (word e or f, up or down) of least predicted change of
error_2^2, from the exact gradient and the Gauss-Newton diagonal of error_2^2 at the table's own weights:
    dE(move) = grad * step + GNdiag * step^2 .
The assignment minimises the sum of these predicted changes. Output: a `rounded.npz` for the finishing chain
(r50_chem_lbl_table7w.py: words, W, U, tree_edges), with the rows of each group in address order.

    python r50_chem_lbl_wtable.py TABLE.npz OUT.npz --wires "14:0:0:1,..." [--no-assign]
"""
import argparse
import os
import sys
import time

import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_round import N, Target, decode_np  # noqa: E402

STEP = np.pi / 64.0


def lin(addr, mask):
    x = np.asarray(addr, np.int64) & int(mask)
    out = np.zeros_like(x)
    while np.any(x):
        out ^= x & 1
        x = x >> 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("out")
    ap.add_argument("--wires", required=True)
    ap.add_argument("--no-assign", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    z = np.load(a.table)
    edges = np.array(z["tree_edges"], np.int64)
    W = np.array(z["W"], np.float64)
    R, B, C = W.shape
    T_ = N - 1
    w2 = np.array(z["angle_words"], np.int64)[N:]
    zi = np.mod(w2 + 64, 128) - 64                       # signed words
    assert zi.min() >= -32 and zi.max() <= 32, "the two-body words are not in the signed range"
    wires = [tuple(int(v) for v in s.split(":")) for s in a.wires.split(",")]
    assert all(b == 0 for (_, _, b, _) in wires)
    tgt = Target()
    I, J = tgt.i, tgt.j
    s = np.sqrt(np.where(I == J, 1.0, 2.0))
    Tt = np.sqrt(tgt.M) * tgt.T
    th = zi * STEP
    u = decode_np(zi / 64.0, edges)
    u = u / np.linalg.norm(u, axis=1, keepdims=True)
    U = u.reshape(R, B, N)
    e2_0 = tgt.error2_np(U, W)
    P = (U[:, :, I] * U[:, :, J]) * s                    # (R, B, pairs)
    At = np.einsum("rbc,rbp->rcp", W, P).reshape(R * C, -1)
    Rm = At.T @ At - Tt
    AR = (At @ Rm).reshape(R, C, -1)
    At = At.reshape(R, C, -1)
    print(f"table {a.table}: {(R, B, C)}; error_2 at its own weights {1e3 * e2_0:.4f} mHa = {1e6 * e2_0 ** 2:.2f} mHa^2 ({time.time() - t0:.0f} s)", flush=True)
    cols = sorted({5 * k - 5 + j for (k, j, _, _) in wires} | {5 * k + j for (k, j, _, _) in wires})
    assert len(cols) == 2 * len(wires), "a word in two relations: use a chain-free layout"
    new_z = zi.copy().reshape(R, B, T_)
    order = np.zeros((R, B), np.int64)
    tot_pred, nfl, nfl0 = 0.0, 0, 0
    h = 1e-6
    for r in range(R):
        thr = th[r * B:(r + 1) * B]
        zr = zi[r * B:(r + 1) * B]
        G = W[r] @ AR[r]                                 # (B, pairs)
        g = W[r] @ At[r]
        dE = {}                                          # col -> (best dE (B,), dz (B,))
        for c in cols:
            tp, tm = thr.copy(), thr.copy()
            tp[:, c] += h
            tm[:, c] -= h
            v = (decode_np(tp / np.pi, edges) - decode_np(tm / np.pi, edges)) / (2 * h)
            sg = np.sign(np.sum(decode_np(thr / np.pi, edges) * U[r], axis=1))
            v = v * sg[:, None]
            q = (v[:, I] * U[r][:, J] + U[r][:, I] * v[:, J]) * s
            grad = 4.0 * np.sum(q * G, axis=1)
            gn = 2.0 * (np.sum(g * g, axis=1) * np.sum(q * q, axis=1) + np.sum(q * g, axis=1) ** 2)
            up = (grad * STEP + gn * STEP * STEP) * 1e6
            dn = (-grad * STEP + gn * STEP * STEP) * 1e6
            up = np.where(zr[:, c] + 1 <= 31, up, np.inf)
            dn = np.where(zr[:, c] - 1 >= -32, dn, np.inf)
            dE[c] = (np.minimum(up, dn), np.where(up <= dn, 1, -1))
        A = 76 + 127 * r + np.arange(B)
        nat = np.zeros((B, len(wires)), np.int64)
        fc = np.zeros((B, len(wires)))
        for wi, (k, j, _, m) in enumerate(wires):
            e, f = 5 * k - 5 + j, 5 * k + j
            nat[:, wi] = (zr[:, e] ^ zr[:, f]) & 1
            fc[:, wi] = np.minimum(dE[e][0], dE[f][0])
        want = np.stack([lin(A, m) for (_, _, _, m) in wires], axis=1)
        cost = np.einsum("bw,baw->ba", fc, (nat[:, None, :] != want[None, :, :]).astype(np.float64))
        nfl0 += int((nat != want).sum())
        if a.no_assign:
            ri, ci = np.arange(B), np.arange(B)
        else:
            ri, ci = linear_sum_assignment(cost)
        order[r, ci] = ri                                # order[r, position] = old row
        for b_, pos in zip(ri, ci):
            for wi, (k, j, _, m) in enumerate(wires):
                if nat[b_, wi] != want[pos, wi]:
                    e, f = 5 * k - 5 + j, 5 * k + j
                    c = e if dE[e][0][b_] <= dE[f][0][b_] else f
                    new_z[r, b_, c] += dE[c][1][b_]
                    tot_pred += float(dE[c][0][b_])
                    nfl += 1
        print(f"  group {r}: relations broken in the present order {int((nat != want).sum())}, after the assignment {int((nat[ri] != want[ci]).sum())}; "
              f"predicted change of error_2^2 {cost[ri, ci].sum():+.4f} mHa^2 ({time.time() - t0:.0f} s)", flush=True)
    # permute, check, evaluate
    Wn = np.stack([W[r][order[r]] for r in range(R)])
    zn = np.stack([new_z[r][order[r]] for r in range(R)]).reshape(R * B, T_)
    addr = np.arange(R * B) + 76
    for (k, j, _, m) in wires:
        assert np.array_equal((zn[:, 5 * k + j] ^ zn[:, 5 * k - 5 + j]) & 1, lin(addr, m))
    un = decode_np(zn / 64.0, edges)
    un = un / np.linalg.norm(un, axis=1, keepdims=True)
    e2_1 = tgt.error2_np(un.reshape(R, B, N), Wn)
    print(f"{len(wires)} wires on {R * B} two-body rows: words moved by one step {nfl} (present order would need {nfl0}); predicted change {tot_pred:+.3f} mHa^2; "
          f"EXACT error_2 at the same weights {1e3 * e2_1:.4f} mHa = {1e6 * e2_1 ** 2:.2f} mHa^2 ({1e6 * (e2_1 ** 2 - e2_0 ** 2):+.3f}) ({time.time() - t0:.0f} s)", flush=True)
    np.savez(a.out, words=zn, W=Wn, U=un.reshape(R, B, N), tree_edges=edges, order=order, wires=np.array(wires))
    print("wrote", a.out, flush=True)


if __name__ == "__main__":
    main()
