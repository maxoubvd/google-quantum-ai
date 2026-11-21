#!/usr/bin/env python
"""r49_chem_dmg_lmgpu (design seat R49-CHEM-DMG): the angle-space Levenberg-Marquardt of r49_chem_dmg_lmhold.py in
JAX (float64, one GPU): the refit of a HELD state of the 7-bit rounding at card speed.

    python r49_chem_dmg_lmgpu.py STATE.npz ROUNDED.npz OUT_DIR [--minutes 10] [--iters 100000] [--cg 1500] [--cgtol 1e-2]
                                 [--lam 1e-2] [--free-all] [--check]

STATE.npz: theta (units of pi), w (R, B, C), frozen (bool) of tools/r48_chem_q7c_round.py; ROUNDED.npz: any file with
`tree_edges` of the same run. Frozen words are held exactly; free angles and hard-sign weights (w = s v, v >= 0, active
set) move. Model, Gauss-Newton product and exact diagonal as in r49_chem_dmg_lmhold.py (checked there against finite
differences); here the push and the pull are jax.linearize / its transpose of the same map (th, v) -> At.
OUT_DIR/last.npz: theta (units of pi), w, frozen, hist; OUT_DIR/log.txt.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r49_chem_dmg_fitlib as L  # noqa: E402

N = 76


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("state")
    ap.add_argument("rounded")
    ap.add_argument("out")
    ap.add_argument("--minutes", type=float, default=10.0)
    ap.add_argument("--iters", type=int, default=100000)
    ap.add_argument("--cg", type=int, default=1500)
    ap.add_argument("--cgtol", type=float, default=1e-2)
    ap.add_argument("--lam", type=float, default=1e-2)
    ap.add_argument("--lam-min", type=float, default=1e-7)
    ap.add_argument("--free-all", action="store_true")
    ap.add_argument("--check", action="store_true", help="compare value, gradient, GN product and diagonal with the numpy tool")
    ap.add_argument("--start-ref", type=float, default=144.49)
    ap.add_argument("--save-every", type=int, default=20)
    ap.add_argument("--cache", default="out/target_cache.npz")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    logf = open(os.path.join(a.out, "log.txt"), "a")

    def log(*x):
        s = " ".join(str(y) for y in x)
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    t0 = time.time()
    jax, jnp = L.jax_setup(True)
    try:                                    # a persistent compilation cache: the second launch skips most of the 165 s
        jax.config.update("jax_compilation_cache_dir", os.path.abspath("out/jax_cache"))
        jax.config.update("jax_persistent_cache_min_compile_time_secs", 5.0)
    except Exception:                       # noqa: BLE001 - older JAX: no cache, nothing else changes
        pass
    if a.cache and os.path.exists(a.cache):
        zc = np.load(a.cache)
        T_np, I, J = zc["T"], zc["i"], zc["j"]
        refc = None
    else:
        refc = L.Reference()
        T_np, I, J = refc.T, refc.i, refc.j
        if a.cache:
            os.makedirs(os.path.dirname(a.cache) or ".", exist_ok=True)
            np.savez(a.cache, T=T_np, i=I, j=J)
    m = np.where(I == J, 1.0, 2.0)
    sm_np = np.sqrt(m)
    Tt = jnp.asarray(np.outer(sm_np, sm_np) * T_np)
    sm = jnp.asarray(sm_np)
    Ij, Jj = jnp.asarray(I), jnp.asarray(J)
    z = np.load(a.state)
    edges = np.array(np.load(a.rounded)["tree_edges"], np.int64)
    th0 = np.pi * np.array(z["theta"], np.float64)
    w0 = np.array(z["w"], np.float64)
    frozen = np.array(z["frozen"], bool)
    if a.free_all:
        frozen = np.zeros_like(frozen)
    R, B, C = w0.shape
    P = len(I)
    sign_np = np.sign(w0.sum(axis=1, keepdims=True))
    sign_np[sign_np == 0] = 1.0
    v0 = np.maximum(w0 * sign_np, 0.0)
    sign = jnp.asarray(sign_np)
    decode_pi = L.make_decode(edges)

    def decode(th):
        return decode_pi(th / jnp.pi)

    def amap(th, v):
        u = decode(th)
        pt = (sm * u[:, Ij] * u[:, Jj]).reshape(R, B, P)
        return jnp.einsum("rbc,rbp->rcp", sign * v, pt).reshape(R * C, P)

    def fval(th, v):
        At = amap(th, v)
        r = At.T @ At - Tt
        return 0.5 * jnp.sum(r * r)

    fval_j = jax.jit(fval)

    def diag_fn(th, v, G):
        u = decode(th).reshape(R, B, N)
        dv = jax.vmap(jax.jacfwd(lambda t: decode(t[None])[0]))(th)          # (n, 76, 75)
        dv = jnp.transpose(dv, (0, 2, 1)).reshape(R, B, N - 1, N)
        w = sign * v
        dth, dw = [], []
        for r in range(R):
            Gr = G[r * C:(r + 1) * C, r * C:(r + 1) * C]
            O = u[r] @ u[r].T
            UV = jnp.einsum("aj,bkj->abk", u[r], dv[r])
            T1 = 2.0 * jnp.einsum("ac,abk->cbk", w[r], O[:, :, None] * UV)
            wGw = jnp.einsum("bc,cd,bd->b", w[r], Gr, w[r])
            q2 = 2.0 * jnp.sum(dv[r] ** 2, axis=2)
            dth.append(2.0 * wGw[:, None] * q2 + 2.0 * jnp.einsum("bc,cbk->bk", w[r], T1) ** 2)
            AP = (w[r].T @ (O * O)).T
            dw.append(2.0 * jnp.diag(Gr)[None, :] + 2.0 * AP ** 2)
        return jnp.stack(dth).reshape(R * B, N - 1), jnp.stack(dw)

    mth_j = jnp.asarray((~frozen).astype(np.float64))     # an ARGUMENT of the compiled step: the cache then serves any frozen set
    CG, TOL = a.cg, a.cgtol

    @jax.jit
    def lm_iter(th, v, lam, mth):
        At, jvp_fn = jax.linearize(amap, th, v)
        vjp_fn = jax.linear_transpose(jvp_fn, th, v)
        G = At @ At.T
        r = At.T @ At - Tt
        f = 0.5 * jnp.sum(r * r)
        gth, gv = vjp_fn(2.0 * (At @ r))
        dgt, dgw = diag_fn(th, v, G)
        dgt, dgw = jnp.maximum(dgt, 1e-30), jnp.maximum(dgw, 1e-30)
        mv = ((v > 0) | (gv < 0)).astype(th.dtype)
        bt, bv = -gth * mth, -gv * mv

        def gn(d1, d2):
            dA = jvp_fn(d1, d2)
            return vjp_fn(2.0 * (G @ dA + (At @ dA.T) @ At))

        def op(d1, d2):
            o1, o2 = gn(d1 * mth, d2 * mv)
            return mth * (o1 + lam * dgt * d1), mv * (o2 + lam * dgw * d2)

        i1, i2 = mth / ((1.0 + lam) * dgt), mv / ((1.0 + lam) * dgw)
        z1, z2 = i1 * bt, i2 * bv
        rz = jnp.sum(bt * z1) + jnp.sum(bv * z2)
        b2 = rz

        def cond(s):
            return (s[0] < CG) & (s[9] > TOL * TOL * b2)

        def body(s):
            k, x1, x2, r1, r2, p1, p2, rz_, _, _ = s
            a1, a2 = op(p1, p2)
            pAp = jnp.sum(p1 * a1) + jnp.sum(p2 * a2)
            al = rz_ / jnp.maximum(pAp, 1e-300)
            x1, x2 = x1 + al * p1, x2 + al * p2
            r1, r2 = r1 - al * a1, r2 - al * a2
            y1, y2 = i1 * r1, i2 * r2
            rz2 = jnp.sum(r1 * y1) + jnp.sum(r2 * y2)
            be = rz2 / jnp.maximum(rz_, 1e-300)
            return (k + 1, x1, x2, r1, r2, y1 + be * p1, y2 + be * p2, rz2, pAp, rz2)

        s0 = (0, jnp.zeros_like(th), jnp.zeros_like(v), bt, bv, z1, z2, rz, rz, rz)
        s = jax.lax.while_loop(cond, body, s0)
        d1, d2 = s[1] * mth, s[2] * mv
        o1, o2 = gn(d1, d2)
        pred = -(jnp.sum(gth * d1) + jnp.sum(gv * d2) + 0.5 * (jnp.sum(o1 * d1) + jnp.sum(o2 * d2)))
        th2, v2 = th + d1, jnp.maximum(v + d2, 0.0)
        f2 = fval(th2, v2)
        gnorm = jnp.sqrt(jnp.sum(bt * bt) + jnp.sum(bv * bv))
        return f, f2, th2, v2, pred, s[0], gnorm, jnp.sum(mv)

    th, v = jnp.asarray(th0), jnp.asarray(v0)
    f_in = float(fval_j(th, v))
    log(f"[lmgpu] {a.state}: shape {(R, B, C)}, frozen {frozen.mean():.4f}, error_2^2 {2e6 * f_in:.3f} mHa^2 (damage "
        f"{2e6 * f_in - a.start_ref:.2f}); devices {jax.devices()}  ({time.time() - t0:.0f} s)")

    if a.check:
        from r49_chem_dmg_lmhold import Model
        refc = refc or L.Reference()
        mod = Model(refc, edges, (R, B, C), sign_np)
        mod.setup(th0, v0)
        f_np, (g1, g2) = mod.grad()
        d1n, d2n = mod.diag()
        At, jvp_fn = jax.linearize(amap, th, v)
        vjp_fn = jax.linear_transpose(jvp_fn, th, v)
        G = At @ At.T
        r = At.T @ At - Tt
        gth, gv = vjp_fn(2.0 * (At @ r))
        dgt, dgw = diag_fn(th, v, G)
        log(f"[check] f numpy {f_np:.12e} jax {f_in:.12e}")
        log(f"[check] gradient: max |diff| theta {np.abs(np.asarray(gth) - g1).max():.2e} (max {np.abs(g1).max():.2e}), "
            f"weights {np.abs(np.asarray(gv) - g2).max():.2e} (max {np.abs(g2).max():.2e})")
        log(f"[check] diagonal: max rel diff theta {np.abs(np.asarray(dgt) - d1n).max() / d1n.max():.2e}, weights "
            f"{np.abs(np.asarray(dgw) - d2n).max() / d2n.max():.2e}")
        rng = np.random.default_rng(0)
        e1, e2 = rng.standard_normal(th0.shape), rng.standard_normal(v0.shape)
        o1, o2 = mod.gn(e1, e2)
        dA = jvp_fn(jnp.asarray(e1), jnp.asarray(e2))
        p1, p2 = vjp_fn(2.0 * (G @ dA + (At @ dA.T) @ At))
        log(f"[check] GN product: max rel diff {np.abs(np.asarray(p1) - o1).max() / np.abs(o1).max():.2e} / "
            f"{np.abs(np.asarray(p2) - o2).max() / np.abs(o2).max():.2e}")
        return

    lam, it, hist = a.lam, 0, []
    t_lm = None                             # the clock of --minutes starts after the first (compiling) step
    while it < a.iters and (t_lm is None or (time.time() - t_lm) < 60 * a.minutes):
        tries, accepted = 0, False
        while not accepted and tries < 8:
            tries += 1
            f, f2, th2, v2, pred, ncg, gnorm, nfree = lm_iter(th, v, lam, mth_j)
            f, f2, pred = float(f), float(f2), float(pred)
            rho = (f - f2) / max(pred, 1e-300)
            if f2 < f:
                accepted = True
                th, v = th2, v2
                if rho > 0.6:
                    lam = max(lam / 3.0, a.lam_min)
                elif rho < 0.25:
                    lam *= 2.0
            else:
                lam *= 4.0
        it += 1
        if t_lm is None:
            t_lm = time.time()
        e2 = 2e6 * (f2 if accepted else f)
        hist.append((time.time() - t0, e2))
        if it <= 10 or it % 10 == 0:
            log(f"[lmgpu] it {it:5d}  error_2^2 {e2:9.3f} mHa^2  damage {e2 - a.start_ref:8.3f}  rho {rho:6.3f}  lam {lam:.1e}  "
                f"cg {int(ncg):4d}  tries {tries}  |g| {float(gnorm):.2e}  free v {int(nfree)}  ({time.time() - t0:.0f} s)")
        if not accepted:
            log("[lmgpu] no step accepted; stop")
            break
        if it % a.save_every == 0:
            np.savez(os.path.join(a.out, "last.npz"), theta=np.asarray(th) / np.pi, w=sign_np * np.asarray(v), frozen=frozen,
                     hist=np.array(hist))
    np.savez(os.path.join(a.out, "last.npz"), theta=np.asarray(th) / np.pi, w=sign_np * np.asarray(v), frozen=frozen,
             hist=np.array(hist))
    f_out = float(fval_j(th, v))
    log(f"[lmgpu] DONE {it} steps: error_2^2 {2e6 * f_out:.3f} mHa^2, damage {2e6 * f_out - a.start_ref:.2f} (state: "
        f"{2e6 * f_in - a.start_ref:.2f})  ({time.time() - t0:.0f} s)")
    json.dump(dict(state=a.state, steps=it, start=2e6 * f_in, end=2e6 * f_out, seconds=time.time() - t0),
              open(os.path.join(a.out, "result.json"), "w"))


if __name__ == "__main__":
    main()
