"""r49_chem_free_costmap: the spectrum of a point (r49_chem_free_spec.py) carried to the ROW ORDER OF A TABLE.

The table's two-body rows are matched to the point's rows inside each group by the overlap of the decoded vectors
(the rounding deals the rows over the address blocks). Output: COST_<key>.npy (table rows, 75), one-body rows 0.

    python r49_chem_free_costmap.py TABLE.npz POINT.npz OUTPREFIX spec_012.npz spec_345.npz spec_678.npz
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_round import decode_np  # noqa: E402

tab, pt, pref = sys.argv[1:4]
z = np.load(tab)
words, edges = np.array(z["angle_words"], np.int64), np.array(z["tree_edges"], np.int64)
p = np.load(pt)
U = np.array(p["U"], np.float64)
R, B, N = U.shape
U /= np.linalg.norm(U, axis=2, keepdims=True)
z0 = np.load(sys.argv[4])
D = {k: np.zeros_like(z0[k]) for k in ("H", "S", "C", "Rw")}
for f in sys.argv[4:]:
    zz = np.load(f)
    for g in zz["groups"]:
        for k in D:
            D[k][g] = zz[k][g]
assert np.array_equal(z0["edges"], edges), "spectrum and table use different trees"
u_tab = decode_np(words[76:] / 64.0, edges).reshape(R, B, N)
perm = np.zeros((R, B), np.int64)
worst = 1.0
for r in range(R):
    ov = np.abs(u_tab[r] @ U[r].T)                 # (table row, point row)
    from scipy.optimize import linear_sum_assignment
    ri, ci = linear_sum_assignment(-ov)
    perm[r, ri] = ci
    worst = min(worst, float(ov[ri, ci].min()))
print(f"rows matched; smallest overlap of a matched pair {worst:.4f}")
for k in D:
    out = np.zeros((words.shape[0], 75))
    for r in range(R):
        out[76 + r * B:76 + (r + 1) * B] = D[k][r][perm[r]]
    np.save(f"{pref}_{k}.npy", out)
np.save(f"{pref}_perm.npy", perm)
print("wrote", pref + "_{H,S,C,Rw}.npy")
