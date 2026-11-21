"""r50_chem_lbl_dead: the EXACT don't-care words of a rounded table (R50-CHEM-LBL).

In a rounded table a word 0 sends nothing to the child's subtree and a quarter turn leaves nothing with the parent:
every word below such a split acts on an amplitude that is exactly zero, so ANY value of it decodes to the same row.
This tool decodes the two-body rows of a table word by word and reports the share of dead words (entering amplitude
exactly zero), by batch and for every pair of edges; and for tie candidates (pairs of word positions) the rows in
which BOTH words are alive (only those rows pay for a tie).

    python r50_chem_lbl_dead.py TABLE.npz [--pairs "70-73,65-68,58-72,52-69,60-63"] [--out dead.npz]
"""
import argparse

import numpy as np

N = 76


def decode_track(words, edges, bits=7):
    """rows (M, 75) of words -> decoded vectors (M, 76) and the entering amplitude (M, 75) of each word."""
    th = np.pi * words.astype(np.float64) / float(1 << (bits - 1))
    M = len(words)
    v = np.zeros((M, N))
    v[:, 0] = 1.0
    ent = np.zeros((M, len(edges)))
    cs, sn = np.cos(th), np.sin(th)
    q = float(1 << (bits - 2))
    wm = np.mod(words, 1 << bits)
    wm = np.where(wm >= (1 << (bits - 1)), wm - (1 << bits), wm)          # signed value
    for t, (p, c) in enumerate(edges):
        a = v[:, p].copy()
        ent[:, t] = a
        c_, s_ = cs[:, t].copy(), sn[:, t].copy()
        c_[np.abs(np.abs(wm[:, t]) - q) < 0.5] = 0.0                       # a quarter turn: exactly nothing stays
        s_[np.mod(wm[:, t], 2 * q) == 0] = 0.0                             # a multiple of a half turn: nothing goes
        v[:, p] = c_ * a
        v[:, c] = s_ * a
    return v, ent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--pairs", default="")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    z = np.load(a.table)
    edges = np.array(z["tree_edges"], np.int64)
    words = np.array(z["angle_words"], np.int64)[N:]
    W = np.array(z["W"], np.float64)
    wn2 = np.sum(W ** 2, axis=2).ravel()
    v, ent = decode_track(words, edges)
    dead = ent == 0.0
    wm = np.where(words >= 64, words - 128, words)
    triv = (wm == 0) | (np.abs(wm) == 32)
    print(f"{len(words)} two-body rows x {words.shape[1]} words: dead (entering amplitude exactly 0) {100 * dead.mean():.1f}%; "
          f"trivial words (0 or a quarter turn) {100 * triv.mean():.1f}%, of the alive words {100 * triv[~dead].mean():.1f}%; "
          f"decoded norm min {np.linalg.norm(v, axis=1).min():.6f}")
    print("  dead share by batch (%): " + " ".join(f"{100 * dead[:, 5 * k:5 * k + 5].mean():.0f}" for k in range(15)))
    print("  dead share by row weight |w|^2 rank (heaviest 10% .. lightest 10%): " + " ".join(
        f"{100 * dead[np.argsort(-wn2)[i * len(wn2) // 10:(i + 1) * len(wn2) // 10]].mean():.0f}" for i in range(10)))
    alive = ~dead
    both = alive.T.astype(np.float64) @ alive.astype(np.float64)          # (75, 75) rows with both alive
    iu = np.triu_indices(75, 1)
    print(f"  rows with BOTH words alive, over all 2,775 word pairs: min {int(both[iu].min())}, 1% {int(np.quantile(both[iu], 0.01))}, "
          f"10% {int(np.quantile(both[iu], 0.1))}, median {int(np.median(both[iu]))}; per word alive rows: min {int(alive.sum(0).min())} "
          f"median {int(np.median(alive.sum(0)))} max {int(alive.sum(0).max())}")
    al = alive.sum(0)
    print("  alive rows per word (position:count), 20 fewest: " + " ".join(f"{t}:{int(al[t])}" for t in np.argsort(al)[:20]))
    for pr in filter(None, a.pairs.split(",")):
        e, f = (int(x) for x in pr.split("-"))
        b = alive[:, e] & alive[:, f]
        me, mf = np.abs(wm[:, e]), np.abs(wm[:, f])
        sameA = int((b & (me == mf)).sum())
        sameB = int((b & (me == 32 - mf)).sum())
        print(f"  pair ({e},{f}): alive rows {int(alive[:, e].sum())} / {int(alive[:, f].sum())}, both alive {int(b.sum())}; of those already "
              f"equal in magnitude {sameA}, already mirror {sameB}; |w|^2 share of the both-alive rows {100 * wn2[b].sum() / wn2.sum():.1f}%")
    if a.out:
        np.savez(a.out, dead=dead, ent=ent, words=words)


if __name__ == "__main__":
    main()
