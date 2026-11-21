"""r50_chem_lbl_wcost: THE PRICE OF DELTA WIRES ON ADDRESS BITS (R50-CHEM-LBL, construction 9), conditional model.

A wire k:j:b:m asks, in every row:  bit b of word 5k+j  XOR  bit b of word 5k-5+j  =  parity(address AND m)
(m = an address mask; 1, 2, 4, .. = one address bit).
The rows of a group may sit at any of the group's 127 addresses (76 + 127 r .. 76 + 127 r + 126), so each group is a
linear assignment of rows to addresses on
    cost[row, a] = sum over the wires whose natural relation bit differs from bit i of a of the row's FLIP PRICE,
flip price = the cheaper of moving one of the two words to its nearest integer with bit b changed, priced as
k * ((x - z')^2 - (x - z)^2) with k = the conditional one-step price C of the word (spectrum file), x = theta / step,
z = round(x): an independent-rounding reading (the judge is the staged model: r50_chem_lbl_wround.py).

    python r50_chem_lbl_wcost.py SPEC_a.npz,SPEC_b.npz,SPEC_c.npz --wires "14:0:0:0,14:1:0:1,.." [--out ADDR.npz]
            [--table TABLE.npz]
Prints per group: the price with the rows in their present order, after the assignment, the number of flips; with
--table also the 76 one-body rows (addresses 0..75; flips weighted by the squared amplitude entering the flipped word).
--out: addr (R, B) = the address of each row of the point, and the wires.
"""
import argparse
import os
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment

STEP = np.pi / 64.0


def lin(addr, mask):
    """parity(address AND mask): any GF(2)-linear function of the address bits (mask 1, 2, 4, .. = one bit)."""
    x = np.asarray(addr, np.int64) & int(mask)
    out = np.zeros_like(x)
    while np.any(x):
        out ^= x & 1
        x = x >> 1
    return out


def flip_price(x, k, b):
    """(price, z') of the cheapest change of bit b of z = round(x): arrays over rows."""
    z = np.round(x).astype(np.int64)
    best = np.full(x.shape, np.inf)
    for dz in range(-(1 << (b + 1)), (1 << (b + 1)) + 1):
        if dz == 0:
            continue
        z2 = z + dz
        ok = (((z2 >> b) & 1) != ((z >> b) & 1)) & (z2 >= -32) & (z2 <= 32)
        p = k * ((x - z2) ** 2 - (x - z) ** 2)
        best = np.where(ok & (p < best), p, best)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("spec")
    ap.add_argument("--wires", required=True)
    ap.add_argument("--out", default="")
    ap.add_argument("--table", default="")
    a = ap.parse_args()
    D = {}
    for f in a.spec.split(","):
        z = np.load(f)
        for key in ("C", "theta"):
            if key not in D:
                D[key] = np.zeros_like(z[key])
            for g in z["groups"]:
                D[key][g] = z[key][g]
    R, B, T_ = D["C"].shape
    wires = [tuple(int(v) for v in s.split(":")) for s in a.wires.split(",")]
    X = D["theta"] / STEP
    Z = np.round(X).astype(np.int64)
    tot0 = tot1 = 0.0
    nfl = 0
    addr = np.zeros((R, B), np.int64)
    print(f"{len(wires)} wires {a.wires}; conditional price of ONE step of every word, summed: {D['C'].sum():.0f} (the staged damage is 0.21 x that / 12)")
    for r in range(R):
        A = 76 + 127 * r + np.arange(B)
        nat = np.zeros((B, len(wires)), np.int64)
        fp = np.zeros((B, len(wires)))
        for w, (k, j, b, i) in enumerate(wires):
            e, f = 5 * k - 5 + j, 5 * k + j
            nat[:, w] = ((Z[r][:, e] >> b) & 1) ^ ((Z[r][:, f] >> b) & 1)
            fp[:, w] = np.minimum(flip_price(X[r][:, e], D["C"][r][:, e], b), flip_price(X[r][:, f], D["C"][r][:, f], b))
        want = np.stack([lin(A, i) for (_, _, _, i) in wires], axis=1)            # (addresses, wires)
        cost = np.einsum("bw,baw->ba", fp, (nat[:, None, :] != want[None, :, :]).astype(np.float64))
        c0 = float(cost[np.arange(B), np.arange(B)].sum())
        ri, ci = linear_sum_assignment(cost)
        c1 = float(cost[ri, ci].sum())
        addr[r, ri] = A[ci]
        fl = int((nat[ri] != want[ci]).sum())
        nrow = int((nat[ri] != want[ci]).any(axis=1).sum())
        tot0, tot1, nfl = tot0 + c0, tot1 + c1, nfl + fl
        print(f"  group {r}: present order {c0:8.3f} -> assigned {c1:8.4f} mHa^2; flips {fl} in {nrow} rows (of {B * len(wires)} relations)", flush=True)
    print(f"ALL TWO-BODY ROWS: present order {tot0:.2f} -> assigned {tot1:.4f} mHa^2 (conditional, independent rounding); flips {nfl} of {R * B * len(wires)}")
    if a.table:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from r50_chem_lbl_dead import decode_track
        tz = np.load(a.table)
        w1 = np.array(tz["angle_words"], np.int64)[:76]
        edges = np.array(tz["tree_edges"], np.int64)
        _, ent = decode_track(w1, edges)
        m = ent ** 2
        A = np.arange(76)
        nat = np.zeros((76, len(wires)), np.int64)
        fp = np.zeros((76, len(wires)))
        for w, (k, j, b, i) in enumerate(wires):
            e, f = 5 * k - 5 + j, 5 * k + j
            nat[:, w] = ((w1[:, e] >> b) & 1) ^ ((w1[:, f] >> b) & 1)
            fp[:, w] = np.minimum(m[:, e], m[:, f])
        want = np.stack([lin(A, i) for (_, _, _, i) in wires], axis=1)
        cost = np.einsum("bw,baw->ba", fp, (nat[:, None, :] != want[None, :, :]).astype(np.float64))
        ri, ci = linear_sum_assignment(cost)
        mis = nat[ri] != want[ci]
        print(f"ONE-BODY ROWS (table {os.path.basename(a.table)}; proxy = the squared amplitude entering the softer word of a flipped relation): present order "
              f"{cost[np.arange(76), np.arange(76)].sum():.4f} -> assigned {cost[ri, ci].sum():.6f}; flips {int(mis.sum())} in {int(mis.any(axis=1).sum())} rows; "
              f"largest single flipped mass {float((fp[ri] * mis).max()):.2e}")
    if a.out:
        np.savez(a.out, addr=addr, wires=np.array(wires))
        print("wrote", a.out)


if __name__ == "__main__":
    main()
