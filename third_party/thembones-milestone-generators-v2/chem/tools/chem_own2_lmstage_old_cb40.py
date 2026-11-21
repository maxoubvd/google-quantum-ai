#!/usr/bin/env python
"""chem_own2_lmstage (R8 holder 6, Nov 16 19:15Z): the IN-STAGE second-order refit of the staged 7-bit rounding.

A Levenberg-Marquardt step (Gauss-Newton + lam x exact GN diagonal, PCG preconditioned by that diagonal) on every
FREE angle and on the hard-sign weights, the frozen words held exactly. The arithmetic of one step is R49-CHEM-DMG's
`tools/r49_chem_dmg_lmgpu.py` (lm_iter, diag_fn), copied here unchanged in substance; what is new:
  * it is a callable that the staged tool runs INSIDE a stage (after its L-BFGS refit, or instead of it), in the same
    process: no stop / hold / hand-back / restart, the compiled step is reused at every stage;
  * it runs in float64 under `jax.enable_x64(True)` (a context), so the staged tool's float32 L-BFGS outside it is
    byte-for-byte the recipe it was;
  * walled free angles (entries whose crossing of +-pi/2 the gauge rules forbid) are projected into the signed range
    after every step (the staged tool would clip them afterwards);
  * the damping lam carries over from stage to stage;
  * a stage's LM stops at --lm-steps, at --lm-seconds, or when the last three accepted steps returned less than
    --lm-stop mHa^2 together.
Carriers must be off (K == 0): the decode has no half turns here.

Model (as in lmgpu): f = 1/2 || At^T At - Tt ||_F^2, At[(r,c), p] = sqrt(m_p) sum_b s_rc v_rbc u_rb,i_p u_rb,j_p,
Tt = sqrt(m m^T) T; error_2^2 = 2 f. Variables: free angles (radians) and v >= 0 with the channel sign s_rc fixed.
"""
from __future__ import annotations

import os
import time

import numpy as np

N = 76


class StageLM:
    def __init__(self, T, I, J, edges, shape, lo, hi, cg=1500, cgtol=1e-2, lam=1e-2, lam_min=1e-7, cache_dir=None):
        import jax
        self.jax = jax
        if cache_dir:
            try:
                jax.config.update("jax_compilation_cache_dir", os.path.abspath(cache_dir))
                jax.config.update("jax_persistent_cache_min_compile_time_secs", 5.0)
            except Exception:       # noqa: BLE001 - no cache: only the first compile is slower
                pass
        self.R, self.B, self.C = shape
        self.lo_r, self.hi_r = np.pi * lo, np.pi * (hi - 1e-9)
        self.lam, self.lam_min = lam, lam_min
        self.CG, self.TOL = cg, cgtol
        self.edges = np.array(edges, np.int64)
        self.levels = _levels(self.edges)
        m = np.where(I == J, 1.0, 2.0)
        self.sm_np = np.sqrt(m)
        self.Tt_np = np.outer(self.sm_np, self.sm_np) * T
        self.I, self.J = np.asarray(I), np.asarray(J)
        self._built = False

    # ------------------------------------------------------------------ the compiled step (built inside x64)
    def _build(self):
        jax = self.jax
        import jax.numpy as jnp
        R, B, C = self.R, self.B, self.C
        P = len(self.I)
        edges, lv = self.edges, self.levels
        Tt = jnp.asarray(self.Tt_np)
        sm = jnp.asarray(self.sm_np)
        Ij, Jj = jnp.asarray(self.I), jnp.asarray(self.J)
        lo_r, hi_r = self.lo_r, self.hi_r
        CG, TOL = self.CG, self.TOL

        def decode(th):             # th in radians, (n, 75) -> (n, 76): the Givens tree from e_0
            c, s = jnp.cos(th), jnp.sin(th)
            u = jnp.zeros((th.shape[0], N), th.dtype).at[:, 0].set(1.0)
            for ts in lv:
                p, q = edges[ts, 0], edges[ts, 1]
                up, uq = u[:, p], u[:, q]
                u = u.at[:, p].set(c[:, ts] * up - s[:, ts] * uq).at[:, q].set(s[:, ts] * up + c[:, ts] * uq)
            return u

        def amap(th, v, sign):
            u = decode(th)
            pt = (sm * u[:, Ij] * u[:, Jj]).reshape(R, B, P)
            return jnp.einsum("rbc,rbp->rcp", sign * v, pt).reshape(R * C, P)

        def fval(th, v, sign):
            At = amap(th, v, sign)
            r = At.T @ At - Tt
            return 0.5 * jnp.sum(r * r)

        def diag_fn(th, v, sign, G):
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

        @jax.jit
        def lm_iter(th, v, sign, lam, mth, wl):
            amap_s = lambda t, x: amap(t, x, sign)          # noqa: E731
            At, jvp_fn = jax.linearize(amap_s, th, v)
            vjp_fn = jax.linear_transpose(jvp_fn, th, v)
            G = At @ At.T
            r = At.T @ At - Tt
            f = 0.5 * jnp.sum(r * r)
            gth, gv = vjp_fn(2.0 * (At @ r))
            dgt, dgw = diag_fn(th, v, sign, G)
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
            th2 = th + d1
            th2 = jnp.where(wl > 0, jnp.clip(th2, lo_r, hi_r), th2)     # walled free angles stay in the range
            v2 = jnp.maximum(v + d2, 0.0)
            f2 = fval(th2, v2, sign)
            gnorm = jnp.sqrt(jnp.sum(bt * bt) + jnp.sum(bv * bv))
            return f, f2, th2, v2, pred, s[0], gnorm

        self._lm_iter = lm_iter
        self._fval = jax.jit(fval)
        self._built = True

    # ------------------------------------------------------------------ one stage's refit
    def refit(self, theta, w, frozen, wall, steps=25, seconds=25.0, stop=0.01, lam0=0.0, lam_start=0.0, stop_lam=0.0):
        """theta (R*B, 75) in units of pi (frozen entries on their words), w (R, B, C), frozen / wall (bool).
        lam0 > 0: the stage starts with lam = max(carried lam, lam0) (fewer rejected first steps after a freeze).
        lam_start > 0: the stage starts with lam = lam_start whatever was carried (a reset).
        stop_lam > 0: the stop rule (last 3 steps < stop) is only applied while lam < stop_lam: an over-damped LM is
        not taken for a converged one (holder 6's run went quiet from stage 72 that way: lam 3-26, 3 tiny steps).
        Returns theta (units of pi), w, info dict (error_2^2 in mHa^2 before / after, steps, seconds, cg)."""
        jax = self.jax
        t0 = time.time()
        with jax.enable_x64(True):
            import jax.numpy as jnp
            if not self._built:
                self._build()
            sign_np = np.sign(w.sum(axis=1, keepdims=True))
            sign_np[sign_np == 0] = 1.0
            v_np = np.maximum(w * sign_np, 0.0)
            th = jnp.asarray(np.pi * np.asarray(theta, np.float64))
            v = jnp.asarray(v_np)
            sign = jnp.asarray(sign_np)
            mth = jnp.asarray((~frozen).astype(np.float64))
            wl = jnp.asarray((wall & ~frozen).astype(np.float64))
            f_in = float(self._fval(th, v, sign))
            f_cur, gains, it, cgs, rejected = f_in, [], 0, [], 0
            lam = lam_start if lam_start > 0 else max(self.lam, lam0)
            while it < steps and (time.time() - t0) < seconds:
                tries, accepted = 0, False
                while not accepted and tries < 8:
                    tries += 1
                    f, f2, th2, v2, pred, ncg, gnorm = self._lm_iter(th, v, sign, lam, mth, wl)
                    f, f2, pred = float(f), float(f2), float(pred)
                    cgs.append(int(ncg))
                    rho = (f - f2) / max(pred, 1e-300)
                    if f2 < f:
                        accepted = True
                        th, v = th2, v2
                        if rho > 0.6:
                            lam = max(lam / 3.0, self.lam_min)
                        elif rho < 0.25:
                            lam *= 2.0
                    else:
                        rejected += 1
                        lam *= 4.0
                it += 1
                if not accepted:
                    break
                gains.append(2e6 * (f_cur - f2))
                f_cur = f2
                if len(gains) >= 3 and sum(gains[-3:]) < stop and (stop_lam <= 0 or lam < stop_lam):
                    break
            self.lam = max(lam, self.lam_min)
            th_out = np.asarray(th, np.float64) / np.pi
            w_out = sign_np * np.asarray(v, np.float64)
        th_out = np.where(frozen, theta, th_out)            # frozen words exactly as they were
        return th_out, w_out, {"e2_in": 2e6 * f_in, "e2_out": 2e6 * f_cur, "steps": it, "rejected": rejected,
                               "seconds": time.time() - t0, "cg_mean": float(np.mean(cgs)) if cgs else 0.0,
                               "lam": self.lam}


class SkkGPU:
    """The weight-compensated stiffness S_kk of every entry (tools/r48_chem_q7c_sdiag.py group_sdiag, same formulas)
    on the GPU in float64: per group, H_aa's diagonal minus the part the group's own weights absorb (Schur complement
    over the B x C weights, Cholesky with the same 1e-10 trace ridge). The tangents d u / d theta come from jacfwd of
    the decode (exact) instead of central differences. Rows in the caller's order; theta in units of pi (+ K = 0)."""

    def __init__(self, I, J, edges, shape):
        import jax
        self.jax = jax
        self.R, self.B, self.C = shape
        self.edges = np.array(edges, np.int64)
        self.levels = _levels(self.edges)
        self.I, self.J = np.asarray(I), np.asarray(J)
        self._fn = None

    def _build(self):
        jax = self.jax
        import jax.numpy as jnp
        import jax.scipy.linalg as jsl
        B, C = self.B, self.C
        edges, lv = self.edges, self.levels
        Ij, Jj = jnp.asarray(self.I), jnp.asarray(self.J)
        s = jnp.asarray(np.sqrt(np.where(self.I == self.J, 1.0, 2.0)))
        bk = np.repeat(np.arange(B), N - 1)
        wb, wc = np.repeat(np.arange(B), C), np.tile(np.arange(C), B)

        def decode(th):
            c, sn = jnp.cos(th), jnp.sin(th)
            u = jnp.zeros((th.shape[0], N), th.dtype).at[:, 0].set(1.0)
            for ts in lv:
                p, q = edges[ts, 0], edges[ts, 1]
                up, uq = u[:, p], u[:, q]
                u = u.at[:, p].set(c[:, ts] * up - sn[:, ts] * uq).at[:, q].set(sn[:, ts] * up + c[:, ts] * uq)
            return u

        @jax.jit
        def group(th_r, W_r, A_r):          # th_r (B, 75) radians, W_r (B, C), A_r (C, P)
            u = decode(th_r)
            V = jnp.transpose(jax.vmap(jax.jacfwd(lambda t: decode(t[None])[0]))(th_r), (0, 2, 1))   # (B, 75, 76)
            p = (u[:, Ij] * u[:, Jj]) * s
            q = ((V[:, :, Ij] * u[:, None, Jj] + u[:, None, Ij] * V[:, :, Jj]) * s).reshape(B * (N - 1), -1)
            g = (W_r @ A_r) * s
            ar = A_r * s
            qq = jnp.sum(q * q, axis=1)
            gg = jnp.sum(g * g, axis=1)
            qg = jnp.sum(q * g[bk], axis=1)
            Haa_d = 2.0 * (gg[bk] * qq + qg * qg)
            Gqp, Gqa = q @ p.T, q @ ar.T
            Ggp, Gga = g @ p.T, g @ ar.T
            Gpp, Gaa, Gap = p @ p.T, ar @ ar.T, ar @ p.T
            Haw = 2.0 * (Gga[bk][:, wc] * Gqp[:, wb] + Ggp[bk][:, wb] * Gqa[:, wc])
            Hww = 2.0 * (Gaa[wc][:, wc] * Gpp[wb][:, wb] + Gap[wc][:, wb] * Gap[wc][:, wb].T)
            nw = B * C
            lam = 1e-10 * jnp.trace(Hww) / nw
            L = jnp.linalg.cholesky(Hww + lam * jnp.eye(nw))
            Y = jsl.solve_triangular(L, Haw.T, lower=True)          # (nw, na): Haw H^-1 Haw^T = |Y_k|^2
            sd = Haa_d - jnp.sum(Y * Y, axis=0)
            return jnp.maximum(sd, 0.0)

        self._fn = group

    def table(self, theta, w):
        """theta (R*B, 75) units of pi, w (R, B, C) -> S_kk (R*B, 75)."""
        jax = self.jax
        with jax.enable_x64(True):
            import jax.numpy as jnp
            if self._fn is None:
                self._build()
            R, B, C = self.R, self.B, self.C
            th = jnp.asarray(np.pi * np.asarray(theta, np.float64)).reshape(R, B, N - 1)
            W = jnp.asarray(np.asarray(w, np.float64))
            out = []
            # A of every group (C, P) from the decoded rows: the same unit rows the CPU table uses
            for r in range(R):
                ur = np.asarray(self._decode_rows(np.asarray(theta, np.float64).reshape(R, B, N - 1)[r]))
                Ar = np.einsum("bc,bp->cp", np.asarray(w, np.float64)[r], ur[:, self.I] * ur[:, self.J])
                out.append(np.asarray(self._fn(th[r], W[r], jnp.asarray(Ar))))
            return np.stack(out).reshape(R * B, N - 1)

    def _decode_rows(self, theta_pi):
        th = np.pi * theta_pi
        v = np.zeros((th.shape[0], N))
        v[:, 0] = 1.0
        for t, (p, q) in enumerate(self.edges):
            a = v[:, p].copy()
            v[:, p] = np.cos(th[:, t]) * a
            v[:, q] = np.sin(th[:, t]) * a
        return v / np.linalg.norm(v, axis=1, keepdims=True)


def _levels(edges):
    """Runs of consecutive edges touching disjoint orbitals (one vectorised update each)."""
    runs, current, used = [], [], set()
    for t, (p, q) in enumerate(edges):
        if p in used or q in used:
            runs.append(np.array(current))
            current, used = [], set()
        current.append(t)
        used.update((int(p), int(q)))
    return runs + [np.array(current)]
