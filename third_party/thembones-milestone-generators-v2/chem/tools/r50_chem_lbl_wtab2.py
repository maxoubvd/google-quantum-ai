"""r50_chem_lbl_wtab2: DELTA WIRES IN CHAINS IMPOSED ON A FINISHED TABLE (R50-CHEM-LBL, construction 9b): the two-body
rows. r50_chem_lbl_wtable.py for any wire list (chains allowed): for every row and every address of its group the
cheapest set of one-step moves that makes every chain's parities t XOR parity(address AND M_i) (t = the chain's free
bit, taken as the cheaper of the two); the rows are assigned to addresses on that cost (predicted change of error_2^2
from the exact gradient and the Gauss-Newton diagonal at the table's own weights), the moves applied, and the EXACT
error_2 at the same weights printed. Output: a `rounded.npz` for the finishing chain (r50_chem_lbl_table7c.py).

    python r50_chem_lbl_wtab2.py TABLE.npz OUT.npz --wires "14:0:0:1,..."|@file [--no-assign] [--polish N]
--polish N: N exact sweeps, one move per row per sweep: a word of a wire column by +-2, any other word by +-1,
the best one taken when the exact error_2^2 falls.
"""
import argparse
import os
import sys
import time

import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_round import N, Target, decode_np  # noqa: E402
from r50_chem_lbl_chain import chains_of, lin, parse  # noqa: E402

STEP = np.pi / 64.0
INNER = int(os.environ.get("R50_INNER", "3"))        # polish: moves per row per sweep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("out")
    ap.add_argument("--wires", required=True)
    ap.add_argument("--no-assign", action="store_true")
    ap.add_argument("--polish", type=int, default=0)
    ap.add_argument("--floor", type=float, default=None, help="least predicted price of a move, mHa^2 (e.g. 1e-4: never move a word for a predicted gain)")
    ap.add_argument("--top", type=int, default=0, help="with --per-block: the N strongest rows of a group (sum |w|) ..")
    ap.add_argument("--per-block", type=int, default=2, help=".. sit at most this many to a block of 16 addresses (the deal of the coded rounding)")
    ap.add_argument("--order-out", default="", help="write perm.npy: perm[address - 76] = the table's present two-body row (for a rounding's --order)")
    ap.add_argument("--step2-from", type=int, default=99, help="polish: the words of batches >= this move by 2 steps (a bit-0 code lives there)")
    ap.add_argument("--polish-only", action="store_true", help="the table already satisfies the wires: no moves, no assignment, only the polish")
    a = ap.parse_args()
    t0 = time.time()
    z = np.load(a.table)
    edges = np.array(z["tree_edges"], np.int64)
    W = np.array(z["W"], np.float64)
    R, B, C = W.shape
    T_ = N - 1
    w2 = np.array(z["angle_words"], np.int64)[N:]
    zi = np.mod(w2 + 64, 128) - 64
    assert zi.min() >= -32 and zi.max() <= 32, "the two-body words are not in the signed range"
    wires = parse(a.wires)
    colinfo, chains = chains_of(wires)
    cols = sorted(colinfo)
    tgt = Target()
    I, J = tgt.i, tgt.j
    s = np.sqrt(np.where(I == J, 1.0, 2.0))
    Tt = np.sqrt(tgt.M) * tgt.T
    th = zi * STEP
    u = decode_np(zi / 64.0, edges)
    u = u / np.linalg.norm(u, axis=1, keepdims=True)
    U = u.reshape(R, B, N)
    e2_0 = tgt.error2_np(U, W)
    P = (U[:, :, I] * U[:, :, J]) * s
    At = np.einsum("rbc,rbp->rcp", W, P).reshape(R * C, -1)
    Rm = At.T @ At - Tt
    AR = (At @ Rm).reshape(R, C, -1)
    At = At.reshape(R, C, -1)
    print(f"table {a.table}: {(R, B, C)}; error_2 at its own weights {1e3 * e2_0:.4f} mHa = {1e6 * e2_0 ** 2:.2f} mHa^2; "
          f"{len(wires)} wires in {len(chains)} chains on {len(cols)} columns ({time.time() - t0:.0f} s)", flush=True)
    new_z = zi.copy().reshape(R, B, T_)
    order = np.zeros((R, B), np.int64)
    tot_pred, nfl, nfl0 = 0.0, 0, 0
    h = 1e-6
    percol = np.zeros(T_)
    for r in range(R):
        if a.polish_only:
            order[r] = np.arange(B)
            continue
        thr = th[r * B:(r + 1) * B]
        zr = zi[r * B:(r + 1) * B]
        G = W[r] @ AR[r]
        g = W[r] @ At[r]
        dE = {}
        sg = np.sign(np.sum(decode_np(thr / np.pi, edges) * U[r], axis=1))
        for c in cols:
            tp, tm = thr.copy(), thr.copy()
            tp[:, c] += h
            tm[:, c] -= h
            v = (decode_np(tp / np.pi, edges) - decode_np(tm / np.pi, edges)) / (2 * h)
            v = v * sg[:, None]
            q = (v[:, I] * U[r][:, J] + U[r][:, I] * v[:, J]) * s
            grad = 4.0 * np.sum(q * G, axis=1)
            gn = 2.0 * (np.sum(g * g, axis=1) * np.sum(q * q, axis=1) + np.sum(q * g, axis=1) ** 2)
            up = (grad * STEP + gn * STEP * STEP) * 1e6
            dn = (-grad * STEP + gn * STEP * STEP) * 1e6
            up = np.where(zr[:, c] + 1 <= 31, up, np.inf)
            dn = np.where(zr[:, c] - 1 >= -32, dn, np.inf)
            dE[c] = (np.minimum(up, dn) if a.floor is None else np.maximum(np.minimum(up, dn), a.floor), np.where(up <= dn, 1, -1))
            percol[c] += float(np.mean(gn * STEP * STEP) * 1e6) / R
        A = 76 + 127 * r + np.arange(B)
        cost = np.zeros((B, B))
        tsel = []                                        # per chain: (B, B) bool = t of the cheaper side
        for ch in chains:
            cc = [c for (c, _) in ch]
            nat = np.stack([zr[:, c] & 1 for c in cc], axis=1)                # (B, n)
            pat = np.stack([lin(A, M) for (_, M) in ch], axis=1)              # (addresses, n)
            fc = np.stack([dE[c][0] for c in cc], axis=1)                     # (B, n)
            X = (nat[:, None, :] ^ pat[None, :, :]).astype(np.float64)        # mismatch at t = 0
            c1 = np.einsum("bn,ban->ba", fc, X)                               # moves if t = 0
            c0 = fc.sum(axis=1)[:, None] - c1                                 # moves if t = 1
            cost += np.minimum(c0, c1)
            tsel.append(c0 < c1)
        if a.no_assign:
            ri, ci = np.arange(B), np.arange(B)
        else:
            ri, ci = linear_sum_assignment(cost)
            if a.top:                                    # the strongest rows spread over the blocks, as the deal does
                strong = np.argsort(-np.abs(W[r]).sum(axis=1))[:a.top]
                blk = A // 16
                cost2 = cost.copy()
                big = 1e3 + float(np.abs(cost).max())
                for _ in range(200):
                    pos_of = np.zeros(B, np.int64)
                    pos_of[ri] = ci
                    cnt = {}
                    for b_ in strong:
                        cnt.setdefault(int(blk[pos_of[b_]]), []).append(int(b_))
                    over = [(k_, v_) for k_, v_ in cnt.items() if len(v_) > a.per_block]
                    if not over:
                        break
                    k_, v_ = over[0]
                    b_ = min(v_, key=lambda x: np.min(cost[x][blk != k_]) - cost[x, pos_of[x]])
                    cost2[b_, blk == k_] += big
                    ri, ci = linear_sum_assignment(cost2)
                else:
                    raise SystemExit("the strongest rows could not be spread")
        order[r, ci] = ri
        g_pred, g_n = 0.0, 0
        for b_, pos in zip(ri, ci):
            for ch, ts in zip(chains, tsel):
                t = int(ts[b_, pos])
                for (c, M) in ch:
                    if (int(zr[b_, c]) & 1) != (t ^ int(lin(A[pos], M))):
                        new_z[r, b_, c] += dE[c][1][b_]
                        g_pred += float(dE[c][0][b_])
                        g_n += 1
        tot_pred += g_pred
        nfl += g_n
        nfl0 += 0
        print(f"  group {r}: words moved {g_n}; predicted change of error_2^2 {g_pred:+.4f} mHa^2 ({time.time() - t0:.0f} s)", flush=True)
    Wn = np.stack([W[r][order[r]] for r in range(R)])
    zn = np.stack([new_z[r][order[r]] for r in range(R)]).reshape(R * B, T_)
    addr = np.arange(R * B) + 76
    for (k, j, _, m) in wires:
        assert np.array_equal((zn[:, 5 * k + j] ^ zn[:, 5 * k - 5 + j]) & 1, lin(addr, m))
    un = decode_np(zn / 64.0, edges)
    un = un / np.linalg.norm(un, axis=1, keepdims=True)
    e2_1 = tgt.error2_np(un.reshape(R, B, N), Wn)
    print(f"{len(wires)} wires on {R * B} two-body rows: words moved by one step {nfl} ({nfl / (R * B):.2f} a row); predicted change {tot_pred:+.3f} mHa^2; "
          f"EXACT error_2 at the same weights {1e3 * e2_1:.4f} mHa = {1e6 * e2_1 ** 2:.2f} mHa^2 ({1e6 * (e2_1 ** 2 - e2_0 ** 2):+.3f}) ({time.time() - t0:.0f} s)", flush=True)
    print("  mean Gauss-Newton price of one step by batch (mHa^2 per word): " +
          " ".join(f"b{k}:{percol[5 * k:5 * k + 5].mean():.5f}" for k in range(15) if percol[5 * k:5 * k + 5].any()), flush=True)
    if a.polish:
        zn, e2_1 = polish(tgt, zn, Wn, edges, sorted(set(cols) | set(range(5 * min(a.step2_from, 15), zn.shape[1]))), a.polish, t0, R, B)
        un = decode_np(zn / 64.0, edges)
        un = un / np.linalg.norm(un, axis=1, keepdims=True)
        for (k, j, _, m) in wires:
            assert np.array_equal((zn[:, 5 * k + j] ^ zn[:, 5 * k - 5 + j]) & 1, lin(addr, m))
        print(f"after the polish: EXACT error_2 {1e3 * e2_1:.4f} mHa = {1e6 * e2_1 ** 2:.2f} mHa^2 ({1e6 * (e2_1 ** 2 - e2_0 ** 2):+.3f}) ({time.time() - t0:.0f} s)", flush=True)
    np.savez(a.out, words=zn, W=Wn, U=un.reshape(R, B, N), tree_edges=edges, order=order, wires=np.array(wires))
    if a.order_out:
        np.save(a.order_out, np.concatenate([r * B + order[r] for r in range(R)]))
        print("wrote", a.order_out, "(perm[address - 76] = present two-body row)", flush=True)
    print("wrote", a.out, flush=True)


def polish(tgt, zn, Wn, edges, cols, sweeps, t0, R, B):
    """Exact coordinate descent at fixed weights: a word of a wire column moves by +-2, any other word of the row
    by +-1. The change of error_2^2 of one row's move is exact (rank-two update of the residual)."""
    I, J = tgt.i, tgt.j
    s = np.sqrt(np.where(I == J, 1.0, 2.0))
    Tt = np.sqrt(tgt.M) * tgt.T
    C = Wn.shape[2]
    allc = list(range(zn.shape[1]))                      # every word of the row may compensate
    step = {c: (2 if c in set(cols) else 1) for c in allc}
    zn = zn.copy()

    def vec(zrow):
        u = decode_np(zrow / 64.0, edges)
        return u / np.linalg.norm(u, axis=1, keepdims=True)

    U = vec(zn)
    P = (U[:, I] * U[:, J]) * s                                        # (rows, pairs)
    At = np.einsum("rbc,rbp->rcp", Wn, P.reshape(R, B, -1)).reshape(R * C, -1)
    Rm = At.T @ At - Tt
    e = float(np.sum(Rm * Rm))
    for sw in range(sweeps):
        nacc = 0
        for row in range(R * B):
            r, b = divmod(row, B)
            w = Wn[r, b]
            al = float(w @ w)
            for _ in range(INNER):
                gA = w @ At[r * C:(r + 1) * C]                         # (pairs,)
                cand, zs = [], []
                for c in allc:
                    for d in (step[c], -step[c]):
                        zc = zn[row].copy()
                        zc[c] += d
                        if -32 <= zc[c] <= 31:
                            zs.append(zc)
                            cand.append((c, d))
                Uc = vec(np.array(zs))
                D = (Uc[:, I] * Uc[:, J]) * s - P[row]                 # (cand, pairs)
                RD = D @ Rm
                Rg = Rm @ gA
                dd, dg, gg = np.sum(D * D, axis=1), D @ gA, float(gA @ gA)
                dE = 4.0 * (D @ Rg) + 2.0 * al * np.sum(RD * D, axis=1) + 2.0 * dd * gg + 2.0 * dg ** 2 + al * al * dd ** 2 + 4.0 * al * dd * dg
                i = int(np.argmin(dE))
                if dE[i] >= -1e-18:
                    break
                c, d = cand[i]
                dv = D[i]
                zn[row, c] += d
                P[row] += dv
                At[r * C:(r + 1) * C] += np.outer(w, dv)
                Rm += np.outer(dv, gA) + np.outer(gA, dv) + al * np.outer(dv, dv)
                e += float(dE[i])
                nacc += 1
        Rm = At.T @ At - Tt
        e = float(np.sum(Rm * Rm))
        e2 = tgt.error2_np(vec(zn).reshape(R, B, N), Wn)
        print(f"  polish sweep {sw + 1}: accepted {nacc}; error_2^2 {1e6 * e2 ** 2:.2f} mHa^2 ({time.time() - t0:.0f} s)", flush=True)
        if nacc == 0:
            break
    return zn, e2


if __name__ == "__main__":
    main()
