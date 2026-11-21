#!/usr/bin/env python
"""chem_own_fit: continuous DF-THC fits on the GPU: our tools/chem_factor_fit.py with one more term.

This file is a copy of tools/chem_factor_fit.py (ours, chem_factor workstream) with the option --flat MU:
    objective = error_2 + kappa * D(w) + MU * sum_rows max(0, s_rb / mean(s) - CAP)^2,   s_rb = |w_rb|^2,
a penalty on rows that are much stronger than the rest (--flat-cap CAP, in units of the mean row strength).
Why: our fits end with three or four rows two to three times stronger than any other; the first entries the
rounding freezes are the top angles of those rows, and long rounding stages repair them far worse than on a
point with a flat top (research/ws_chem_own.md, section 4). The penalty moves their weight onto other rows
through the fit itself.
Second change (dlbfgs), the three float32 repairs of tools/chem_own_round.py (found 2026-11-11, where float64
reached in 12,000 iterations what float32 missed in 40,000, and the repaired float32 equals float64):
  1. the value handed to L-BFGS is measured from the CURRENT float64 reference, not from the error at the start
     of the phase (a fit that started at 0.8 Ha had its value near -0.8: one float32 step is 6e-8 Ha, the
     decrease per iteration 2e-9 Ha, and the line search ran on rounding);
  2. the change of the table since the reference is computed as a change (unit vectors, pairs, table), not as
     A(p) - A0 with A(p) evaluated in float32;
  3. the L-BFGS variables are the changes since the reference, zero at every recentering, with float64 master
     copies on the CPU (a late step is smaller than the last float32 digit of a vector entry or a weight).
Everything else below is unchanged.

chem_factor: continuous DF-THC fits on the GPU (JAX, float64 by default) with shrink steps and checkpoints.

Objective (the public pipeline's): error_2 + kappa * D(w), error_2 in packed pair space (chem_factor_lib).

  python chem_factor_fit.py INIT OUT --phases "lbfgs:4000;dropb:8;lbfgs:20000;adam:200000:1e-4:1e-6"

INIT is an npz with U (R,B,76) and W (R,B,C), or a submission directory (its two-body words are decoded).
Phases, run in order (a finished phase is skipped on restart; Adam and L-BFGS resume from the last checkpoint):
  lbfgs:N              optax L-BFGS, N iterations
  adam:N:lr0:lr1       Adam, geometric learning-rate decay
  dropb:k / dropc:k    fade (x0.5, x0.2, `--fade` L-BFGS iterations each) then remove each group's k weakest
                       vectors / copies, ranked by the exact error_2^2 increase of removing each alone
                       (the ranking is `impacts` of the pipeline's experiments/shrink.py, Apache-2.0)
  addc:k               add k copies per group with 1e-4-scale random weights
OUT holds last.npz, best.npz (lowest error_2 seen at a checkpoint), state.json and log.txt.
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


def impacts(u: np.ndarray, w: np.ndarray, refc) -> tuple[np.ndarray, np.ndarray]:
    """Exact error_2^2 increase from removing each copy (R, C) and each basis vector (R, B) alone."""
    T, M, i, j = refc.T, refc.M, refc.i, refc.j
    m = np.sqrt(np.diag(M))
    pi = u[:, :, i] * u[:, :, j]
    A = np.einsum('rbc,rbp->rcp', w, pi)
    H = np.einsum('rcp,rcq->pq', A, A) - T
    At = A * m
    copies = np.einsum('rcp,rcp->rc', At, At @ H) * -2 + np.einsum('rcp->rc', A * A * m) ** 2
    ut = np.einsum('rbc,rcp->rbp', w, A) * m
    pt = pi * m
    Hp = pt @ H
    beta = np.sum(w ** 2, axis=2)
    inner = -2 * np.sum(ut * Hp, axis=2) + beta * np.sum(pt * Hp, axis=2)
    us, ps = ut / np.sqrt(m), pt / np.sqrt(m)
    up, pp, uu = np.sum(us * ps, axis=2), np.sum(ps * ps, axis=2), np.sum(us * us, axis=2)
    bases = 2 * inner + 2 * uu * pp + 2 * up ** 2 + beta ** 2 * pp ** 2 - 4 * beta * up * pp
    return copies, bases


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("init")
    ap.add_argument("out")
    ap.add_argument("--phases", required=True)
    ap.add_argument("--kappa", type=float, default=5e-5)
    ap.add_argument("--f32", action="store_true")
    ap.add_argument("--fade", type=int, default=600)
    ap.add_argument("--chunk", type=int, default=200)
    ap.add_argument("--report", type=int, default=2000)
    ap.add_argument("--save", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--precond", action="store_true", help="dlbfgs: scaled variables (see run_dlbfgs)")
    ap.add_argument("--memory", type=int, default=30, help="dlbfgs: L-BFGS memory")
    ap.add_argument("--recenter", type=int, default=1000, help="dlbfgs: iterations between float64 reference refreshes")
    ap.add_argument("--flat", type=float, default=0.0, help="weight of the outlier-row penalty (Ha)")
    ap.add_argument("--flat-cap", type=float, default=9.0, help="row strength, in units of the mean, above which a row is penalised")
    ap.add_argument("--stall", type=float, default=2e-5, help="relative decrease per report interval below which L-BFGS counts as stalled")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    logf = open(os.path.join(a.out, "log.txt"), "a")

    def log(*x):
        s = " ".join(str(v) for v in x)
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    jax, jnp = L.jax_setup(x64=not a.f32)
    import optax
    DT = jnp.float32 if a.f32 else jnp.float64
    refc = L.Reference()
    error2, _ = L.make_error2(refc, DT)

    state_path = os.path.join(a.out, "state.json")
    state = {"phase": 0, "step": 0, "best": None}
    if os.path.exists(state_path):
        state = json.load(open(state_path))
        z = np.load(os.path.join(a.out, "last.npz"))
        U, w = np.array(z["U"], np.float64), np.array(z["W"], np.float64)
        log(f"resume at phase {state['phase']} step {state['step']} shape {w.shape}")
    elif os.path.isdir(a.init):
        z = np.load(os.path.join(a.init, "submission.npz"))
        w = np.array(z["W"], np.float64)
        U = L.decode_np(z["angle_words"][N:] / 128.0, z["tree_edges"]).reshape(w.shape[0], w.shape[1], N)
    else:
        z = np.load(a.init)
        U, w = np.array(z["U"], np.float64), np.array(z["W"], np.float64)

    def flat_pen(w_):
        s_ = jnp.sum(w_ * w_, axis=2)
        return jnp.sum(jnp.maximum(s_ / jnp.mean(s_) - a.flat_cap, 0.0) ** 2)

    def objective(p, sb, sc):
        U, w = p
        u = U / jnp.linalg.norm(U, axis=2, keepdims=True)
        we = w * sb[:, :, None] * sc[:, None, :]
        eps = error2(u, we)
        return eps + a.kappa * L.d_penalty(we) + a.flat * flat_pen(we), eps

    def metrics(U, w):
        u = U / np.linalg.norm(U, axis=2, keepdims=True)
        return refc.error2_np(u, w), float(L.d_penalty(jnp.asarray(w)))

    def top(w_):          # the strongest row in units of the median row strength (the number quoted in the log)
        s_ = np.sum(np.asarray(w_, np.float64) ** 2, axis=2)
        return float(s_.max() / np.median(s_))

    def save(U, w, phase, step, eps):
        np.savez(os.path.join(a.out, "last.npz"), U=U, W=w)
        if state["best"] is None or eps < state["best"] - 1e-9 or tuple(state.get("best_shape", [])) != w.shape:
            if state["best"] is None or eps < state["best"] or tuple(state.get("best_shape", [])) != w.shape:
                state["best"], state["best_shape"] = eps, list(w.shape)
                np.savez(os.path.join(a.out, f"best_{w.shape[0]}x{w.shape[1]}x{w.shape[2]}.npz"), U=U, W=w)
        state["phase"], state["step"] = phase, step
        json.dump(state, open(state_path, "w"))

    t0 = time.time()

    def run_lbfgs(U, w, iters, sb=None, sc=None, phase=None, step0=0, label="lbfgs"):
        sb = jnp.ones(w.shape[:2], DT) if sb is None else jnp.asarray(sb, DT)
        sc = jnp.ones((w.shape[0], w.shape[2]), DT) if sc is None else jnp.asarray(sc, DT)
        f = lambda p: objective(p, sb, sc)[0]  # noqa: E731
        opt = optax.lbfgs(memory_size=30)
        vg = optax.value_and_grad_from_state(f)

        @jax.jit
        def run(params, st):
            def body(_, carry):
                p, s = carry
                value, grad = vg(p, state=s)
                upd, s = opt.update(grad, s, p, value=value, grad=grad, value_fn=f)
                return optax.apply_updates(p, upd), s
            return jax.lax.fori_loop(0, a.chunk, body, (params, st))

        params = (jnp.asarray(U, DT), jnp.asarray(w, DT))
        st = opt.init(params)
        e0 = float(objective(params, sb, sc)[1])
        done = step0
        last_eps, stalls = e0, 0
        while done < iters:
            params, st = run(params, st)
            done += a.chunk
            if done % a.report == 0 or done >= iters:
                eps = float(objective(params, sb, sc)[1])
                log(f"  {label} {done:7d}: error_2 {1e3 * eps:8.4f} mHa  D {float(L.d_penalty(params[1])):6.2f}  ({time.time() - t0:.0f} s)")
                if not np.isfinite(eps):
                    raise SystemExit("diverged")
                # float32 L-BFGS freezes when its line search can no longer resolve the decrease: restart the
                # memory once, and give the phase up if the restart does not move it either.
                if last_eps - eps < a.stall * last_eps:
                    stalls += 1
                    if stalls >= 2:
                        log(f"  {label}: stalled twice in a row, phase ended at {done}")
                        break
                    st = opt.init(params)
                else:
                    stalls = 0
                last_eps = eps
            if phase is not None and (done % a.save == 0) and done < iters:
                save(np.asarray(params[0], np.float64), np.asarray(params[1], np.float64), phase, done, float(objective(params, sb, sc)[1]))
        U, w = np.asarray(params[0], np.float64), np.asarray(params[1], np.float64)
        e1 = float(objective(params, sb, sc)[1])
        log(f"{label} x{iters - step0}: error_2 {1e3 * e0:8.4f} -> {1e3 * e1:8.4f} mHa  shape {w.shape}  ({time.time() - t0:.0f} s)")
        return U, w, e1

    def run_dlbfgs(U, w, iters, phase, step0=0):
        """L-BFGS at float32 speed with float64 accuracy: the objective is written as a DIFFERENCE from a reference
        point (U0, w0) whose residual r0 is computed in float64 on the CPU,

            r = r0 + sqrt(M) * (S + S^T + dA^T dA),  S = A0^T dA,  dA = A - A0,
            eps^2 - eps0^2 = sum dr (2 r0 + dr),      F = (eps - eps0) + kappa (D - D0),

        so nothing subtracts two large nearly equal numbers (G - T) in float32. The reference is refreshed every
        --recenter iterations (which also restarts the L-BFGS memory)."""
        i_, j_ = jnp.asarray(refc.i), jnp.asarray(refc.j)
        sM32 = jnp.asarray(np.sqrt(refc.M), jnp.float32)
        R, B, C = w.shape

        # --precond: the variables are (U~, v) with w = sw v (sw = rms weight) and |U~_b| = the row's weight norm
        # relative to the average. The objective does not depend on |U~_b|, but the curvature along a row's
        # direction scales as |w_b|^2 / |U~_b|^2, so this choice equalises it across rows and with the v block.
        sw = float(np.sqrt(np.mean(w ** 2))) if a.precond else 1.0
        if a.precond:
            nb = np.linalg.norm(w, axis=2)
            nb = np.maximum(nb / np.sqrt(np.mean(nb ** 2)), 1e-2)
            U = U / np.linalg.norm(U, axis=2, keepdims=True) * nb[:, :, None]

        opt = optax.lbfgs(memory_size=a.memory)

        @jax.jit
        def run(params, st, A0, r0, eps0, off, U0, n0, u0, p0, v0):   # the reference is an argument: one compilation
            def f(p):
                # p = (dU, dv): CHANGES of the scaled vectors and weights since the float64 reference (zero at
                # every recentering; the master copies are float64 on the CPU). The change of the unit vectors,
                # of their pairs and of the table is computed as a change, never as a difference of two float32
                # evaluations (see tools/chem_own_round.py and research/ws_chem_own.md, 2026-11-11 15:30Z).
                dU, dv = p
                s2 = jnp.sum(dU * (2.0 * U0 + dU), axis=2, keepdims=True)        # |U0 + dU|^2 - |U0|^2
                n1 = jnp.sqrt(n0 * n0 + s2)
                du = (dU - u0 * (s2 / (n1 + n0))) / n1                           # unit vector minus reference unit vector
                dui, duj, ui, uj = du[:, :, i_], du[:, :, j_], u0[:, :, i_], u0[:, :, j_]
                dp = dui * uj + (ui + dui) * duj
                wv = (v0 + dv) * sw
                dA = (jnp.einsum("rbc,rbp->rcp", dv * sw, p0) + jnp.einsum("rbc,rbp->rcp", wv, dp)).reshape(R * C, -1)
                S = A0.T @ dA
                dr = sM32 * (S + S.T + dA.T @ dA)
                delta = jnp.sum(dr * (2.0 * r0 + dr))
                deps = delta / (eps0 + jnp.sqrt(jnp.maximum(eps0 * eps0 + delta, 0.0)))
                # `off` = (eps0 - eps_base) - kappa D_base keeps the value continuous across reference refreshes,
                # so the L-BFGS memory and line-search state carry over
                return off + deps + a.kappa * L.d_penalty(wv) + a.flat * flat_pen(wv)
            vg = optax.value_and_grad_from_state(f)

            def body(_, carry):
                p, s = carry
                value, grad = vg(p, state=s)
                upd, s = opt.update(grad, s, p, value=value, grad=grad, value_fn=f)
                return optax.apply_updates(p, upd), s
            return jax.lax.fori_loop(0, a.chunk, body, (params, st))

        def reference(U_, w_):
            n0 = np.linalg.norm(U_, axis=2, keepdims=True)
            u = U_ / n0
            p0 = u[:, :, refc.i] * u[:, :, refc.j]
            A0 = np.einsum("rbc,rbp->rcp", w_, p0).reshape(R * C, -1)
            r0 = np.sqrt(refc.M) * (A0.T @ A0 - refc.T)
            s_ = np.sum(w_ * w_, axis=2)
            pen = float(np.sum(np.maximum(s_ / np.mean(s_) - a.flat_cap, 0.0) ** 2))
            D0 = 0.5 * float(np.sum(np.abs(w_).sum(1) ** 2 - w_.sum(1) ** 2))
            return A0, r0, float(np.sqrt(np.sum(r0 * r0))), a.kappa * D0 + a.flat * pen, (n0, u, p0)

        f32 = jnp.float32
        done, e_first, last_eps, stalls = step0, None, None, 0
        v = w / sw                                  # float64 master copies: U (scaled rows) and v
        zero = (jnp.zeros(U.shape, f32), jnp.zeros(w.shape, f32))
        params = zero
        st = opt.init(params)
        while done < iters:
            A0, r0, eps0, pen0, (n0, u0, p0) = reference(U, v * sw)
            if e_first is None:
                e_first = last_eps = eps0
            else:
                # new origin and new zero: the point kept by L-BFGS moves with the origin, and the value kept for
                # the line search is the value at the current point, which is 0 in the new zero
                prev = optax.tree_utils.tree_get(st, "params")
                st = optax.tree_utils.tree_set(st, params=(prev[0] - params[0], prev[1] - params[1]),
                                               value=jnp.asarray(0.0, f32))
                params = zero
            ref = (jnp.asarray(A0, f32), jnp.asarray(r0, f32), f32(eps0), f32(-pen0), jnp.asarray(U, f32),
                   jnp.asarray(n0, f32), jnp.asarray(u0, f32), jnp.asarray(p0, f32), jnp.asarray(v, f32))
            for _ in range(max(1, a.recenter // a.chunk)):
                params, st = run(params, st, *ref)
                done += a.chunk
            U = U + np.asarray(params[0], np.float64)          # rows are NOT renormalized: the state carries over
            v = v + np.asarray(params[1], np.float64)
            w = v * sw
            eps, D = metrics(U, w)
            if done % a.report < a.recenter:
                log(f"  dlbfgs {done:7d}: error_2 {1e3 * eps:8.4f} mHa  D {D:6.2f}  top {top(w):5.1f}  ({time.time() - t0:.0f} s)")
            if not np.isfinite(eps):
                raise SystemExit("diverged")
            stalls = stalls + 1 if last_eps - eps < a.stall * last_eps * a.recenter / a.report else 0
            last_eps = eps
            if stalls >= 3:
                log(f"  dlbfgs: no progress over three recenterings, phase ended at {done}")
                break
            if phase is not None and done % a.save < a.recenter and done < iters:
                save(U, w, phase, done, eps)
        eps, D = metrics(U, w)
        log(f"dlbfgs x{done - step0}: error_2 {1e3 * e_first:8.4f} -> {1e3 * eps:8.4f} mHa  shape {w.shape}  ({time.time() - t0:.0f} s)")
        return U, w, eps

    def run_adam(U, w, steps, lr0, lr1, phase, step0=0):
        sb, sc = jnp.ones(w.shape[:2], DT), jnp.ones((w.shape[0], w.shape[2]), DT)
        opt = optax.inject_hyperparams(optax.adam)(learning_rate=lr0)
        vg = jax.value_and_grad(lambda p: objective(p, sb, sc), has_aux=True)

        @jax.jit
        def run(params, st, lr):
            st.hyperparams['learning_rate'] = lr

            def body(_, carry):
                p, s = carry
                (_, _), g = vg(p)
                upd, s = opt.update(g, s, p)
                return optax.apply_updates(p, upd), s
            params, st = jax.lax.fori_loop(0, a.chunk, body, (params, st))
            return params, st, vg(params)[0][1]

        params = (jnp.asarray(U, DT), jnp.asarray(w, DT))
        st = opt.init(params)
        e0 = float(vg(params)[0][1])
        best = (e0, U, w)
        step = step0
        ta = time.time()
        while step < steps:
            lr = lr0 * (lr1 / lr0) ** (step / steps)
            params, st, eps = run(params, st, lr)
            step += a.chunk
            eps = float(eps)
            if eps < best[0]:
                best = (eps, np.asarray(params[0], np.float64), np.asarray(params[1], np.float64))
            if step % a.report == 0:
                log(f"  adam {step:8d}: error_2 {1e3 * eps:8.4f} mHa  D {float(L.d_penalty(params[1])):6.2f}  lr {lr:.2e}  best {1e3 * best[0]:8.4f}  "
                    f"({time.time() - t0:.0f} s, {(step - step0) / (time.time() - ta):.0f} steps/s)")
            if step % a.save == 0 and step < steps:
                save(np.asarray(params[0], np.float64), np.asarray(params[1], np.float64), phase, step, eps)
        log(f"adam x{steps - step0}: error_2 {1e3 * e0:8.4f} -> best {1e3 * best[0]:8.4f} mHa  ({time.time() - t0:.0f} s)")
        return best[1], best[2], best[0]

    e, D = metrics(U, w)
    log(f"start shape {w.shape}: error_2 {1e3 * e:.4f} mHa  D {D:.2f}  phases {a.phases}")
    phases = [p for p in a.phases.split(";") if p]
    rng = np.random.default_rng(a.seed)
    for k, ph in enumerate(phases):
        if k < state["phase"]:
            continue
        step0 = state["step"] if k == state["phase"] else 0
        parts = ph.split(":")
        kind = parts[0]
        if kind == "lbfgs":
            U, w, e = run_lbfgs(U, w, int(parts[1]), phase=k, step0=step0)
        elif kind == "dlbfgs":
            U, w, e = run_dlbfgs(U, w, int(parts[1]), k, step0)
        elif kind == "adam":
            U, w, e = run_adam(U, w, int(parts[1]), float(parts[2]), float(parts[3]), k, step0)
        elif kind in ("dropb", "dropc"):
            n = int(parts[1])
            u = U / np.linalg.norm(U, axis=2, keepdims=True)
            cop, bas = impacts(u, w, refc)
            R, B, C = w.shape
            sb, sc = np.ones((R, B)), np.ones((R, C))
            if kind == "dropb":
                kill = np.argsort(bas, axis=1)[:, :n]
                log(f"dropb {n}: removal impacts (mHa^2) mean {1e6 * np.take_along_axis(bas, kill, 1).mean():.3f} max {1e6 * np.take_along_axis(bas, kill, 1).max():.3f}")
            else:
                kill = np.argsort(cop, axis=1)[:, :n]
                log(f"dropc {n}: removal impacts (mHa^2) mean {1e6 * np.take_along_axis(cop, kill, 1).mean():.3f} max {1e6 * np.take_along_axis(cop, kill, 1).max():.3f}")
            for s in (0.5, 0.2):
                if kind == "dropb":
                    np.put_along_axis(sb, kill, s, axis=1)
                else:
                    np.put_along_axis(sc, kill, s, axis=1)
                # fold the scale into w so the fade continues from the scaled weights
                U, w, e = run_lbfgs(U, w, a.fade, sb=sb, sc=sc, label=f"fade x{s}")
            keep = np.ones((R, B if kind == "dropb" else C), bool)
            np.put_along_axis(keep, kill, False, axis=1)
            if kind == "dropb":
                U = np.stack([U[r][keep[r]] for r in range(R)])
                w = np.stack([w[r][keep[r]] for r in range(R)])
            else:
                w = np.stack([w[r][:, keep[r]] for r in range(R)])
            e, D = metrics(U, w)
            log(f"{ph}: now {w.shape}, error_2 {1e3 * e:.4f} mHa before the refit")
        elif kind == "addc":
            n = int(parts[1])
            w = np.concatenate([w, 1e-4 * rng.normal(size=(w.shape[0], w.shape[1], n))], axis=2)
            e, D = metrics(U, w)
            log(f"{ph}: now {w.shape}, error_2 {1e3 * e:.4f} mHa")
        else:
            raise SystemExit(f"unknown phase {ph}")
        save(U, w, k + 1, 0, e)
    e, D = metrics(U, w)
    log(f"done shape {w.shape}: error_2 {1e3 * e:.4f} mHa  D {D:.2f}  top {top(w):.1f}  ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
