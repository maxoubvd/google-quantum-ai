"""C2 (22nd holder): merge the dead-word rewrites of several OWNER loaders (R63-CHEM-WIRED's tools/r63_chem_wired_dc.py
writes <emit>/table_kK.npz = the whole table with only K's owned dead words changed) into one table.
  python -I tools/chem_ld2_e22merge.py BASE.npz OUT.npz K=table_kK.npz [K=table_kK.npz ...]
Checks [P]: every owner changed words only in its own batches K-1 and K (columns 5(K-1) .. 5K+4) and only at dead words
(r62_chem_cache_dead.dead_mask of the BASE table: rewriting a dead word keeps its whole subtree at zero amplitude, so
the dead set cannot move); the owners' batch sets are disjoint. Writes angle_words + tree_edges (and every other array
of BASE unchanged) and prints the words changed per batch and the loaders whose delta changed (they must be rebuilt)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r62_chem_cache_dead import dead_mask  # noqa: E402


def main():
    base, out, pairs = sys.argv[1], sys.argv[2], sys.argv[3:]
    z = np.load(base)
    arrs = {k: z[k] for k in z.files}
    w0 = np.asarray(arrs["angle_words"], np.int64)
    edges = np.asarray(arrs["tree_edges"], np.int64)
    dead = dead_mask(w0, edges, 7)
    w = w0.copy()
    owned = {}
    for p in pairs:
        k, f = p.split("=", 1)
        k = int(k)
        wk = np.asarray(np.load(f)["angle_words"], np.int64)
        assert wk.shape == w0.shape, (f, wk.shape)
        d = wk != w0
        assert not (d & ~dead).any(), f"owner k{k}: a live word changed"
        cols = sorted(set(np.nonzero(d)[1].tolist()))
        bats = sorted(set(c // 5 for c in cols))
        assert set(bats) <= {k - 1, k}, f"owner k{k} changed batches {bats}"
        for b in bats:
            assert b not in owned, f"batch {b} owned by k{owned[b]} and k{k}"
            owned[b] = k
        w[d] = wk[d]
        print(f"owner k{k}: {int(d.sum())} dead words rewritten in batches {bats}")
    diff = w != w0
    per = [int(diff[:, 5 * b:5 * b + 5].sum()) for b in range(15)]
    print(f"words changed per batch {per}; total {int(diff.sum())} of {int(dead.sum())} dead")
    ch = [b for b in range(15) if per[b]]
    rebuild = sorted(set(k for b in ch for k in (b, b + 1) if k <= 14))
    print(f"batches changed {ch}; loaders whose delta changed (must be rebuilt on this table): {rebuild}")
    arrs["angle_words"] = w.astype(z["angle_words"].dtype)
    np.savez(out, **arrs)
    print(f"written {out}")


if __name__ == "__main__":
    main()
