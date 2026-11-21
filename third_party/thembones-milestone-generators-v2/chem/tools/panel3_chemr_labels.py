"""panel3_chemr_labels: orbital labels for the log-depth tree chosen for the rounding (CHEM-R).

On the log-depth tree an orbital j sits, for a row whose bulk is at orbital h, in the sibling clade of depth
d = v2(label_j - label_h) (2-adic valuation). The amplitude of every such clade is set by ONE angle of full mass
(step 0.049 of the bulk); inside a clade the angles are relative (fine). A clade that holds no other large
orbital of the row is pure tail: its coarse angle cannot be re-expressed by the weights. Proxy for the damage
made there:  cost = sum_rows |w|^2 sum_{pure-tail clades} min(A^2 / bulk^2, step^2 / 12),  A = clade amplitude.
Simulated annealing over the labels of orbitals 1..75 (orbital 0 keeps label 0). Writes perm (perm[old] = new).
    python panel3_chemr_labels.py INIT.npz OUT.npy [--iters 200000]
"""
import argparse, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("init"); ap.add_argument("out"); ap.add_argument("--iters", type=int, default=200000)
ap.add_argument("--seed", type=int, default=0); ap.add_argument("--core", type=float, default=0.04)
ap.add_argument("--t0", type=float, default=0.004)
a = ap.parse_args()
z = np.load(a.init)
U, W = np.array(z["U"], float), np.array(z["W"], float)
U /= np.linalg.norm(U, axis=2, keepdims=True)
R, B, N = U.shape
u2 = (U ** 2).reshape(R * B, N)
wn2 = np.sum(W ** 2, axis=2).ravel()
keep = np.argsort(-wn2)[: int(0.6 * R * B)]            # the rows that carry the damage
u2, wn2 = u2[keep], wn2[keep] / wn2.sum()
home = np.argmax(u2, axis=1)
core = u2 >= a.core                                     # large orbitals of the row (never "tail")
S12 = (2 * np.pi / 128) ** 2 / 12
V2 = np.zeros((N + 64, N + 64), int)
lab = np.arange(N + 64)
d = np.abs(lab[:, None] - lab[None, :])
for k in range(7):
    V2[(d % (2 ** (k + 1)) != 0) & (d % (2 ** k) == 0)] = k
NDEP = 7


def cost(perm):
    """perm[orbital] = label."""
    lh = perm[home]                                      # (rows,)
    dep = V2[lh[:, None], perm[None, :]]                 # (rows, N) depth of each orbital's clade for the row
    tot = 0.0
    for k in range(NDEP):
        m = dep == k
        m[np.arange(len(home)), home] = False
        A2 = np.sum(u2 * m, axis=1)
        pure = ~np.any(core & m, axis=1)
        bulk2 = np.maximum(u2[np.arange(len(home)), home], 1e-3)
        tot += np.sum(wn2 * pure * np.minimum(A2 / bulk2, S12))
    return tot / S12


rng = np.random.default_rng(a.seed)
perm = np.arange(N)
c0 = cost(perm)
rnd = [cost(np.concatenate([[0], 1 + rng.permutation(N - 1)])) for _ in range(5)]
print(f"proxy (in units of one full coarse step per unit |w|^2): today's labels {c0:.4f}; random labels {np.mean(rnd):.4f}", flush=True)
best, cb, cur, cc = perm.copy(), c0, perm.copy(), c0
T0 = a.t0 * c0
for it in range(a.iters):
    i, j = rng.integers(1, N, 2)
    if i == j:
        continue
    cur[i], cur[j] = cur[j], cur[i]
    cn = cost(cur)
    Tk = T0 * (1 - it / a.iters) + 1e-9
    if cn <= cc or rng.random() < np.exp((cc - cn) / Tk):
        cc = cn
        if cc < cb:
            cb, best = cc, cur.copy()
    else:
        cur[i], cur[j] = cur[j], cur[i]
    if it % 20000 == 0:
        print(f"  it {it}: current {cc:.4f} best {cb:.4f}", flush=True)
print(f"best proxy {cb:.4f} ({cb / c0:.3f} of today's)", flush=True)
np.save(a.out, best)
