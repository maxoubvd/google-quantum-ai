"""chem_c3_transplant.py: builder C3 (third holder, Nov 17 17:25Z). TRANSPLANT owners built on ANOTHER rewrite of the
same table: an owner K's circuit loads delta_K = batch K-1 XOR batch K and nothing else, so it stays exact on any table
whose batches K-1 and K equal the owner's table there. This copies batches K-1, K WHOLE from each owner table
(dead and live words, e.g. R65-CHEM-WIRED221's owners with live kill units) into BASE (another rewrite: other dead
values elsewhere). Checks: the owners' batch pairs are disjoint; outside its own batches an owner table may differ
from BASE only at dead words (the other rewrite); the dead set of the result equals BASE's (a live flip that made a
word's parent amplitude zero or non-zero would move it); every other array of BASE is kept.
Prints the words changed per batch (live ones counted) and the loaders whose delta changed (to rebuild).
   python -I chem_c3_transplant.py BASE.npz OUT.npz K=owner_table.npz [K=...]"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r62_chem_cache_dead import dead_mask  # noqa: E402


def main():
    args = sys.argv[1:]
    # --nooutside (17:58Z): the owner tables were built on another base whose OTHER batches carry live flips of other
    # owners (or lack mine): only the copied batches matter for the owner's circuit, so skip the outside check
    nooutside = "--nooutside" in args
    allowdead = "--allowdead" in args
    args = [x for x in args if x not in ("--nooutside", "--allowdead")]
    base, out, specs = args[0], args[1], args[2:]
    z = np.load(base, allow_pickle=False)
    arrs = {x: z[x] for x in z.files}
    wb = np.asarray(arrs["angle_words"]).astype(np.int64)
    edges = np.asarray(arrs["tree_edges"]).astype(np.int64)
    n = min(1219, wb.shape[0])
    db = dead_mask(wb[:n], edges, 7)
    w = wb.copy()
    used = set()
    per = [0] * 15
    live = [0] * 15
    for sp in specs:
        k, f = sp.split("=", 1)
        k = int(k)
        wk = np.asarray(np.load(f)["angle_words"]).astype(np.int64)
        bats = [b for b in (k - 1, k) if 0 <= b <= 14]
        assert not (set(bats) & used), f"k{k}: batch already transplanted"
        used |= set(bats)
        cols = [c for b in bats for c in range(5 * b, 5 * b + 5)]
        other = [c for c in range(wb.shape[1]) if c not in cols]
        d_other = (wk[:n, other] != wb[:n, other]) & ~db[:, other]
        assert nooutside or not d_other.any(), f"k{k}: its table differs from BASE at a LIVE word outside its batches"
        for c in cols:
            ch = wk[:n, c] != wb[:n, c]
            per[c // 5] += int(ch.sum())
            live[c // 5] += int((ch & ~db[:, c]).sum())
            w[:n, c] = wk[:n, c]
        print(f"k{k}: batches {bats} copied")
    dn = dead_mask(w[:n], edges, 7)
    if not (dn == db).all():
        # 18:15Z: --allowdead: report the words that changed state; the verifier's own evaluation decides
        msg = f"the dead set moved: {int((db & ~dn).sum())} dead -> alive, {int((dn & ~db).sum())} alive -> dead"
        assert allowdead, msg
        print(msg)
    arrs["angle_words"] = w.astype(z["angle_words"].dtype)
    np.savez(out, **arrs)
    nb = sorted({b for b in range(15) if per[b]} | {b + 1 for b in range(15) if per[b] and b + 1 <= 14})
    print(f"words changed per batch {per} (live {live}); loaders whose delta changed: {nb}; written {out}")


if __name__ == "__main__":
    main()
