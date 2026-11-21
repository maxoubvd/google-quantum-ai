#!/usr/bin/env python
"""chem_lens: "orbit tables" - DF-THC fits in which the vectors of an orbit share their fine Givens angles.

Structure tested.  The two-body rows are split into orbits of 2**k rows (inside a group).  All rows of an orbit
share one stored row of 75 fine angles phi (the low brot-2 bits of every word); a row may differ from it only by
quarter turns (the top two bits of each word) and by the reflected branch:

        theta[row, t] = s[row, t] * phi[orbit(row), t] + q[row, t] * pi/2 ,   s = +-1, q in {0,1,2,3}.

s = -1 has the same |cos|, |sin| at that rotation, so the row has the same amplitude magnitudes as a quarter-turn
variant of the stored row and differs from it by a sign pattern, which bit-7 masks reach (all sign patterns but
one parity).  The invariant of the allowed set is cos(4 theta).  (Version 1 of this tool allowed s = +1 only,
which excludes every sign-related pair of vectors in the signed gauge: its k = 4 run ended at 523 mHa.)
A lookup then loads one stored row per orbit plus the top two bits (see research/ws_chem_lens.md, section 2).

Method (no hand-made masks): start from a trained unstructured table, cluster each group's vectors into orbits
(balanced k-means on cos 4 theta, each rotation weighted by the amplitude mass it carries, so rotations a row
does not use conflict with nothing), then continue with a penalty

        F = error_2 + kappa D(w) + mu * sum_rows,t  ( cos 4 theta - c[orbit, t] )^2

for growing mu (branches and quadrants are chosen by the dynamics), snap, and polish in the exactly constrained
parametrisation (phi, w).  Fits run in float32 with the float64-difference objective of chem_factor_fit (dlbfgs).

  python chem_lens_orbit.py INIT.npz OUT_DIR --k 2 --mus 1e-6,1e-5,1e-4,1e-3,1e-2 --iters 5000 --final 30000
  python chem_lens_orbit.py INIT.npz OUT_DIR --k 0 --final 4000          # control: unstructured, angle space
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_factor_lib as L  # noqa: E402

N = 76


def angles_and_mass(vectors: np.ndarray, edges: np.ndarray):
    """Signed-gauge theta/pi (as chem_factor_lib.signed_angles) and the mass entering each rotation."""
    a = np.array(vectors, dtype=np.float64)
    a /= np.linalg.norm(a, axis=1, keepdims=True)
    theta = np.zeros((len(a), len(edges)))
    mass = np.zeros((len(a), len(edges)))
    for t in range(len(edges) - 1, -1, -1):
        p, q = edges[t]
        s = np.where(a[:, p] > 0, 1.0, np.where(a[:, p] < 0, -1.0, np.where(a[:, q] > 0, -1.0, 1.0)))
        theta[:, t] = np.arctan2(s * a[:, q], s * a[:, p])
        h = np.hypot(a[:, p], a[:, q])
        mass[:, t] = h * h
        a[:, p] = s * h
        a[:, q] = 0.0
    return theta / np.pi, mass


def balanced_orbits(feat: np.ndarray, size: int, rng, sweeps: int = 12, wts: np.ndarray | None = None) -> np.ndarray:
    """Balanced k-means: orbit index for every row of feat; orbit sizes `size` (the last one may be smaller).
    With `wts` (same shape as feat) the cost of a row in an orbit is sum_t wts[row, t] (feat[row, t] - centre[t])^2
    and centres are weighted means: a rotation that carries no amplitude of a row conflicts with nothing."""
    from scipy.optimize import linear_sum_assignment
    n = feat.shape[0]
    k = -(-n // size)
    caps = [size] * (n // size) + ([n % size] if n % size else [])
    slots = np.repeat(np.arange(k), caps)
    w = np.ones_like(feat) if wts is None else wts

    def dist_to(seed):          # symmetric conflict with one seed row
        return (np.sqrt(w * w[seed][None]) * (feat - feat[seed][None]) ** 2).sum(-1)

    seeds = [int(rng.integers(n))]
    d = dist_to(seeds[0])
    for _ in range(k - 1):
        pr = np.maximum(d, 0) + 1e-12
        seeds.append(int(rng.choice(n, p=pr / pr.sum())))
        d = np.minimum(d, dist_to(seeds[-1]))
    cent = feat[seeds].copy()
    assign = None
    for it in range(sweeps):
        if it == 0:
            cost = np.stack([dist_to(sd) for sd in seeds], axis=1)
        else:
            cost = (w[:, None, :] * (feat[:, None, :] - cent[None]) ** 2).sum(-1)      # (n, k)
        r, c = linear_sum_assignment(cost[:, slots])
        new = np.empty(n, np.int64)
        new[r] = slots[c]
        if assign is not None and np.array_equal(new, assign):
            break
        assign = new
        for j in range(k):
            mm = assign == j
            cent[j] = (w[mm] * feat[mm]).sum(0) / np.maximum(w[mm].sum(0), 1e-30)
    return assign


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("init")
    ap.add_argument("out")
    ap.add_argument("--k", type=int, default=2, help="orbit size 2**k (0: unstructured control)")
    ap.add_argument("--mus", default="1e-6,1e-5,1e-4,1e-3,1e-2")
    ap.add_argument("--iters", type=int, default=5000, help="L-BFGS iterations per penalty stage")
    ap.add_argument("--final", type=int, default=30000, help="L-BFGS iterations of the constrained polish")
    ap.add_argument("--kappa", type=float, default=5e-5)
    ap.add_argument("--chunk", type=int, default=200)
    ap.add_argument("--recenter", type=int, default=1000)
    ap.add_argument("--memory", type=int, default=30)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--assign", default="fine", choices=["fine", "abs", "random", "order"])
    ap.add_argument("--stall", type=float, default=1e-5)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    logf = open(os.path.join(a.out, "log.txt"), "a")

    def log(*x):
        s = " ".join(str(v) for v in x)
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    jax, jnp = L.jax_setup(x64=False)
    import optax
    refc = L.Reference()
    edges = L.log_depth_tree(N)
    decode = L.make_decode(edges)
    i_, j_ = jnp.asarray(refc.i), jnp.asarray(refc.j)
    sM32 = jnp.asarray(np.sqrt(refc.M), jnp.float32)
    sM64 = np.sqrt(refc.M)
    rng = np.random.default_rng(a.seed)

    z = np.load(a.init)
    U, w = np.array(z["U"], np.float64), np.array(z["W"], np.float64)
    R, B, C = w.shape
    U = U / np.linalg.norm(U, axis=2, keepdims=True)
    theta, mass = angles_and_mass(U.reshape(-1, N), edges)
    sw = float(np.sqrt(np.mean(w ** 2)))
    nb = np.linalg.norm(w, axis=2).reshape(-1)
    nb = np.maximum(nb / np.sqrt(np.mean(nb ** 2)), 3e-2)
    srow = 1.0 / nb                                           # per-row variable scale (equalises curvature)

    # ---------------------------------------------------------------- float64 reference pieces (numpy)
    def A_np(th, w_):
        u = L.decode_np(th, edges).reshape(R, B, N)
        return np.einsum("rbc,rbp->rcp", w_, u[:, :, refc.i] * u[:, :, refc.j]).reshape(R * C, -1)

    def reference(th, w_):
        A0 = A_np(th, w_)
        r0 = sM64 * (A0.T @ A0 - refc.T)
        return A0, r0, float(np.sqrt(np.sum(r0 * r0)))

    def dpen_np(w_):
        return 0.5 * float(np.sum(np.abs(w_).sum(1) ** 2 - w_.sum(1) ** 2))

    # ---------------------------------------------------------------- orbits
    size = 1 << a.k
    if a.k > 0:
        orb = np.zeros(R * B, np.int64)
        base = 0
        for r in range(R):
            sl = slice(r * B, (r + 1) * B)
            fw = None
            if a.assign == "fine":
                f, fw = np.cos(4 * np.pi * theta[sl]), mass[sl]
            elif a.assign == "abs":
                f = np.abs(U[r])
            else:
                f = None
            if f is None:
                o = np.arange(B) // size if a.assign == "order" else rng.permutation(B) // size
            else:
                o = balanced_orbits(f, size, rng, wts=fw)
            orb[sl] = base + o
            base += int(o.max()) + 1
        n_orb = base
    else:
        orb, n_orb = np.arange(R * B), R * B
    orb_j = jnp.asarray(orb)
    log(f"init {a.init} shape {(R, B, C)} k {a.k}: {n_orb} stored rows for {R * B} vectors (assign {a.assign})")

    def centres(th):
        """Mass-weighted mean of cos(4 theta) over each orbit (the invariant of the allowed set +-phi + k pi/2)."""
        num = np.zeros((n_orb, N - 1))
        den = np.zeros((n_orb, N - 1))
        np.add.at(num, orb, mass * np.cos(4 * np.pi * th))
        np.add.at(den, orb, mass)
        return np.clip(num / np.maximum(den, 1e-30), -1.0, 1.0)

    # ---------------------------------------------------------------- objective (float32, difference form)
    srow_j = jnp.asarray(srow, jnp.float32)[:, None]

    def make_run(mode, sgn=None):
        opt = optax.lbfgs(memory_size=a.memory)
        sg = None if sgn is None else jnp.asarray(sgn, jnp.float32)

        def build_theta(p, q):
            if mode == "hard":
                return sg * p["phi"][orb_j] + 0.5 * q
            return p["tau"] * srow_j

        @jax.jit
        def run(params, st, A0, r0, eps0, off, mu, q):
            def f(p):
                th = build_theta(p, q)
                u = decode(th).reshape(R, B, N)
                A = jnp.einsum("rbc,rbp->rcp", p["v"] * sw, u[:, :, i_] * u[:, :, j_]).reshape(R * C, -1)
                dA = A - A0
                S = A0.T @ dA
                dr = sM32 * (S + S.T + dA.T @ dA)
                delta = jnp.sum(dr * (2.0 * r0 + dr))
                deps = delta / (eps0 + jnp.sqrt(jnp.maximum(eps0 * eps0 + delta, 0.0)))
                val = off + deps + a.kappa * L.d_penalty(p["v"] * sw)
                if mode == "pen":
                    val = val + mu * jnp.sum((jnp.cos(4 * jnp.pi * th) - p["cen"][orb_j]) ** 2)
                return val
            vg = optax.value_and_grad_from_state(f)

            def body(_, carry):
                p, s = carry
                value, grad = vg(p, state=s)
                upd, s = opt.update(grad, s, p, value=value, grad=grad, value_fn=f)
                return optax.apply_updates(p, upd), s
            return jax.lax.fori_loop(0, a.chunk, body, (params, st))
        return opt, run, build_theta

    t0 = time.time()

    def stage(runner, mode, params, q, mu, iters, label):
        opt, run, build_theta = runner
        st = opt.init(params)
        done, e_first, last, stalls = 0, None, None, 0
        qj = jnp.asarray(q, jnp.float32)
        eps = None
        while done < iters:
            th = np.asarray(build_theta(params, qj), np.float64)
            w_ = np.asarray(params["v"], np.float64) * sw
            A0, r0, eps0 = reference(th, w_)
            if e_first is None:
                e_first = last = eps0
            ref = (jnp.asarray(A0, jnp.float32), jnp.asarray(r0, jnp.float32), jnp.float32(eps0), jnp.float32(eps0 - e_first))
            for _ in range(max(1, a.recenter // a.chunk)):
                params, st = run(params, st, *ref, jnp.float32(mu), qj)
                done += a.chunk
            th = np.asarray(build_theta(params, qj), np.float64)
            w_ = np.asarray(params["v"], np.float64) * sw
            eps = reference(th, w_)[2]
            msg = f"  {label} {done:6d}: error_2 {1e3 * eps:8.4f} mHa  D {dpen_np(w_):6.2f}"
            if mode == "pen":
                dev = (np.cos(4 * np.pi * th) - np.asarray(params["cen"], np.float64)[orb]) ** 2
                msg += f"  mu {mu:.1e}  pen {dev.sum():9.3f}  rms dev of cos4theta {np.sqrt(dev.mean()):6.4f}"
            log(msg + f"  ({time.time() - t0:.0f} s)")
            if not np.isfinite(eps):
                raise SystemExit("diverged")
            if mode != "pen":
                stalls = stalls + 1 if last - eps < a.stall * last else 0
                if stalls >= 3:
                    log(f"  {label}: no progress over three recenterings, ended at {done}")
                    break
            last = eps
        return params, eps

    e0 = reference(theta, w)[2]
    log(f"start: error_2 {1e3 * e0:.4f} mHa  D {dpen_np(w):.2f}")
    zero_q = np.zeros_like(theta)

    if a.k == 0:
        params = {"tau": jnp.asarray(theta / srow[:, None], jnp.float32), "v": jnp.asarray(w / sw, jnp.float32)}
        params, eps = stage(make_run("free"), "free", params, zero_q, 0.0, a.final, "free")
        th = np.asarray(params["tau"], np.float64) * srow[:, None]
        np.savez(os.path.join(a.out, "last.npz"), theta=th, W=np.asarray(params["v"], np.float64) * sw, orb=orb)
        log(f"done control: error_2 {1e3 * eps:.4f} mHa ({time.time() - t0:.0f} s)")
        return

    cen = centres(theta)
    params = {"tau": jnp.asarray(theta / srow[:, None], jnp.float32), "cen": jnp.asarray(cen, jnp.float32),
              "v": jnp.asarray(w / sw, jnp.float32)}
    pen_runner = make_run("pen")
    for mu in [float(x) for x in a.mus.split(",") if x]:
        params, eps = stage(pen_runner, "pen", params, zero_q, mu, a.iters, "pen")
        theta = np.asarray(params["tau"], np.float64) * srow[:, None]
        w = np.asarray(params["v"], np.float64) * sw
        np.savez(os.path.join(a.out, "pen.npz"), theta=theta, cen=np.asarray(params["cen"], np.float64), orb=orb, W=w)
    cen = np.asarray(params["cen"], np.float64)

    # ---------------------------------------------------------------- snap and polish
    ang = np.arccos(np.clip(cen, -1.0, 1.0)) / (4 * np.pi)                    # stored fine angle, in [0, 1/4]
    best = None
    for sg in (1.0, -1.0):
        qq = np.round((theta - sg * ang[orb]) / 0.5)
        dd = np.abs(theta - sg * ang[orb] - 0.5 * qq)
        if best is None:
            best = (dd, np.full(theta.shape, sg), qq)
        else:
            take = dd < best[0]
            best = (np.where(take, dd, best[0]), np.where(take, sg, best[1]), np.where(take, qq, best[2]))
    _, sgn, q = best
    snapped = sgn * ang[orb] + 0.5 * q
    e_snap = reference(snapped, w)[2]
    log(f"snap: error_2 {1e3 * e_snap:.4f} mHa (before the constrained polish); reflected branch share {np.mean(sgn < 0):.3f}; "
        f"quadrant histogram {np.bincount(np.mod(q, 4).astype(int).ravel(), minlength=4).tolist()}")
    params = {"phi": jnp.asarray(ang, jnp.float32), "v": jnp.asarray(w / sw, jnp.float32)}
    params, eps = stage(make_run("hard", sgn), "hard", params, q, 0.0, a.final, "hard")
    ang = np.asarray(params["phi"], np.float64)
    w = np.asarray(params["v"], np.float64) * sw
    th = sgn * ang[orb] + 0.5 * q
    np.savez(os.path.join(a.out, "hard.npz"), phi=ang, q=q, sgn=sgn, orb=orb, W=w, theta=th)
    u = L.decode_np(th, edges).reshape(R, B, N)
    np.savez(os.path.join(a.out, "hard_UW.npz"), U=u, W=w)
    e_fin = refc.error2_np(u, w)
    log(f"done k {a.k}: {n_orb} stored rows, error_2 {1e3 * e_fin:.4f} mHa  D {dpen_np(w):.2f} ({time.time() - t0:.0f} s)")
    json.dump({"k": a.k, "stored": int(n_orb), "error_2_mHa": 1e3 * e_fin}, open(os.path.join(a.out, "result.json"), "w"))


if __name__ == "__main__":
    main()
