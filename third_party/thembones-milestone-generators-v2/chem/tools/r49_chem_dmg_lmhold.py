#!/usr/bin/env python
"""r49_chem_dmg_lmhold (design seat R49-CHEM-DMG): Levenberg-Marquardt refit of a HELD state of the 7-bit rounding,
in ANGLE space, on the CPU (numpy, float64).

    python r49_chem_dmg_lmhold.py STATE.npz ROUNDED.npz OUT_DIR [--minutes 20] [--iters 1000] [--cg 300] [--cgtol 1e-2]
                                  [--lam 1e-2] [--free-all] [--check] [--start-ref MHA2]

STATE.npz   a state of tools/r48_chem_q7c_round.py (theta in units of pi (rows x 75), w (R, B, C), frozen bool);
ROUNDED.npz any file with `tree_edges` of the same run (the labelled tree).
The frozen words stay exactly where they are; every free angle and every weight moves. Hard sign: w = s_rc v with
v >= 0 (s_rc = the sign of the channel's sum in the state), handled as an active set (projected LM).

The model is the builders' (tools/r48_chem_q7c_fitlib.py): pt_b = sqrt(m) u_b u_b^T packed (P = 2926),
At[(r,c)] = sum_b w_rbc pt_rb, residual r = At^T At - Tt, error_2 = ||r||_F, f = 1/2 ||r||^2.
Gauss-Newton product without P x P work: J d = At^T dA + dA^T At, J^T y = pull(2 At y), so
GN d = pull(2 (At At^T) dA + 2 (At dA^T) At) with dA the push of d. Inner solves: PCG with the exact GN diagonal.
Why: the builders' stage refit is float32 L-BFGS; R48-CHEM-LAT measured that a held 1% head was still falling 2 mHa^2
per 10k iterations after 100k. This tool measures how much of the damage booked at a held state a second-order
refit returns.
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


def decode_tan(th, edges):
    """th in radians (n, 75) -> u (n, 76) and dv (n, 75, 76) = d u / d th_k."""
    n = th.shape[0]
    c, s = np.cos(th), np.sin(th)
    v = np.zeros((n, N))
    v[:, 0] = 1.0
    dv = np.zeros((n, N - 1, N))
    for t, (p, q) in enumerate(edges):
        a = v[:, p].copy()
        da = dv[:, :, p].copy()
        v[:, p] = c[:, t] * a
        v[:, q] = s[:, t] * a
        dv[:, :, p] = c[:, t, None] * da
        dv[:, :, q] = s[:, t, None] * da
        dv[:, t, p] -= s[:, t] * a
        dv[:, t, q] += c[:, t] * a
    return v, dv


def decode(th, edges):
    n = th.shape[0]
    c, s = np.cos(th), np.sin(th)
    v = np.zeros((n, N))
    v[:, 0] = 1.0
    for t, (p, q) in enumerate(edges):
        a = v[:, p].copy()
        v[:, p] = c[:, t] * a
        v[:, q] = s[:, t] * a
    return v


class Model:
    def __init__(self, refc, edges, shape, sign):
        self.i, self.j = refc.i, refc.j
        self.sm = np.sqrt(np.where(self.i == self.j, 1.0, 2.0))
        self.Tt = np.sqrt(refc.M) * refc.T
        self.edges = edges
        self.R, self.B, self.C = shape
        self.sign = sign                              # (R, 1, C)
        self.P = len(self.i)

    # -- the map (th, v) -> At
    def amap(self, th, v):
        u = decode(th, self.edges)
        pt = self.sm * u[:, self.i] * u[:, self.j]
        w = self.sign * v
        R, B, C = self.R, self.B, self.C
        ptr = pt.reshape(R, B, self.P)
        return np.concatenate([w[r].T @ ptr[r] for r in range(R)], axis=0)

    def value(self, th, v):
        At = self.amap(th, v)
        r = At.T @ At - self.Tt
        return 0.5 * float(np.sum(r * r))

    # -- linearisation at a point
    def setup(self, th, v):
        R, B, C = self.R, self.B, self.C
        self.u, self.dv = decode_tan(th, self.edges)
        self.w = self.sign * v
        self.pt = (self.sm * self.u[:, self.i] * self.u[:, self.j]).reshape(R, B, self.P)
        self.At = np.concatenate([self.w[r].T @ self.pt[r] for r in range(R)], axis=0)
        self.G = self.At @ self.At.T

    def grad(self):
        r = self.At.T @ self.At - self.Tt
        f = 0.5 * float(np.sum(r * r))
        return f, self.pull(2.0 * (self.At @ r))

    def pull(self, GA):
        R, B, C = self.R, self.B, self.C
        GAr = GA.reshape(R, C, self.P)
        gw = np.stack([self.pt[r] @ GAr[r].T for r in range(R)])                 # (R, B, C)
        Z = np.concatenate([self.w[r] @ GAr[r] for r in range(R)], axis=0)       # (n, P)
        S = np.zeros((Z.shape[0], N, N))
        S[:, self.i, self.j] = Z * self.sm
        S = S + S.transpose(0, 2, 1)
        gu = np.einsum("nij,nj->ni", S, self.u)
        gth = np.einsum("nkj,nj->nk", self.dv, gu)
        return gth, gw * self.sign

    def push(self, dth, dvv):
        R, B, C = self.R, self.B, self.C
        du = np.einsum("nkj,nk->nj", self.dv, dth)
        dpt = (self.sm * (du[:, self.i] * self.u[:, self.j] + self.u[:, self.i] * du[:, self.j])).reshape(R, B, self.P)
        dw = self.sign * dvv
        return np.concatenate([dw[r].T @ self.pt[r] + self.w[r].T @ dpt[r] for r in range(R)], axis=0)

    def gn(self, dth, dvv):
        dA = self.push(dth, dvv)
        return self.pull(2.0 * (self.G @ dA + (self.At @ dA.T) @ self.At))

    def diag(self):
        R, B, C = self.R, self.B, self.C
        u = self.u.reshape(R, B, N)
        dv = self.dv.reshape(R, B, N - 1, N)
        dth = np.zeros((R, B, N - 1))
        dw = np.zeros((R, B, C))
        for r in range(R):
            Gr = self.G[r * C:(r + 1) * C, r * C:(r + 1) * C]
            w = self.w[r]
            O = u[r] @ u[r].T                                           # (B, B) u_a . u_b
            UV = np.einsum("aj,bkj->abk", u[r], dv[r])                  # u_a . dv_bk
            M = (O[:, :, None] * UV).reshape(B, B * (N - 1))
            T1 = 2.0 * (w.T @ M).reshape(C, B, N - 1)                   # (At q_bk)_c
            wGw = np.einsum("bc,cd,bd->b", w, Gr, w)
            q2 = 2.0 * np.sum(dv[r] ** 2, axis=2)                       # ||q_bk||^2
            dth[r] = 2.0 * wGw[:, None] * q2 + 2.0 * np.einsum("bc,cbk->bk", w, T1) ** 2
            AP = (w.T @ (O * O)).T                                      # (B, C): (At pt_b)_(rc)
            dw[r] = 2.0 * np.diag(Gr)[None, :] + 2.0 * AP ** 2
        return dth.reshape(R * B, N - 1), dw


def pcg(op, b, minv, maxit, tol):
    x = [np.zeros_like(b[0]), np.zeros_like(b[1])]
    r = [b[0].copy(), b[1].copy()]
    z = [minv[0] * r[0], minv[1] * r[1]]
    p = [z[0].copy(), z[1].copy()]
    rz = float(np.sum(r[0] * z[0]) + np.sum(r[1] * z[1]))
    b2 = float(np.sum(b[0] * z[0]) + np.sum(b[1] * z[1]))
    it = 0
    for it in range(1, maxit + 1):
        Ap = op(p)
        pAp = float(np.sum(p[0] * Ap[0]) + np.sum(p[1] * Ap[1]))
        if pAp <= 0:
            break
        al = rz / pAp
        x[0] += al * p[0]
        x[1] += al * p[1]
        r[0] -= al * Ap[0]
        r[1] -= al * Ap[1]
        z = [minv[0] * r[0], minv[1] * r[1]]
        rz2 = float(np.sum(r[0] * z[0]) + np.sum(r[1] * z[1]))
        if rz2 <= tol * tol * b2:
            break
        be = rz2 / rz
        rz = rz2
        p = [z[0] + be * p[0], z[1] + be * p[1]]
    return x, it


def lm_run(mod, th, v, frozen, log, minutes=20.0, iters=100000, cg=300, cgtol=1e-2, lam=1e-2, lam_min=1e-7, free_sign=False,
           start_ref=144.49, save=None, save_every=5, tag="lmhold", quiet=False):
    """Projected Levenberg-Marquardt on the free angles (th, radians; frozen = bool mask held exactly) and the weights
    v (>= 0 unless free_sign). Returns (th, v, history [(seconds, error_2^2 in mHa^2)], steps)."""
    t0 = time.time()
    it = 0
    mth = (~frozen).astype(np.float64)
    hist = []
    while it < iters and (time.time() - t0) < 60 * minutes:
        mod.setup(th, v)
        f, (gth, gv) = mod.grad()
        dgt, dgw = mod.diag()
        mv = np.ones_like(v) if free_sign else ((v > 0) | (gv < 0)).astype(np.float64)
        gth_m, gv_m = gth * mth, gv * mv
        gnorm = float(np.sqrt(np.sum(gth_m ** 2) + np.sum(gv_m ** 2)))
        dgt_s, dgw_s = np.maximum(dgt, 1e-30), np.maximum(dgw, 1e-30)
        accepted = False
        tries = 0
        while not accepted and tries < 8:
            tries += 1

            def op(d, lam_=lam):
                o1, o2 = mod.gn(d[0] * mth, d[1] * mv)
                return [mth * (o1 + lam_ * dgt_s * d[0]), mv * (o2 + lam_ * dgw_s * d[1])]
            minv = [mth / ((1.0 + lam) * dgt_s), mv / ((1.0 + lam) * dgw_s)]
            d, ncg = pcg(op, [-gth_m, -gv_m], minv, cg, cgtol)
            d[0] *= mth
            d[1] *= mv
            o1, o2 = mod.gn(d[0], d[1])
            pred = -(float(np.sum(gth * d[0]) + np.sum(gv * d[1])) + 0.5 * float(np.sum(o1 * d[0]) + np.sum(o2 * d[1])))
            th2 = th + d[0]
            v2 = v + d[1] if free_sign else np.maximum(v + d[1], 0.0)
            f2 = mod.value(th2, v2)
            rho = (f - f2) / max(pred, 1e-300)
            if f2 < f:
                accepted = True
                th, v = th2, v2
                if rho > 0.6:
                    lam = max(lam / 3.0, lam_min)
                elif rho < 0.25:
                    lam *= 2.0
            else:
                lam *= 4.0
        it += 1
        e2 = 2e6 * (f2 if accepted else f)
        hist.append((time.time() - t0, e2))
        if not quiet:
            log(f"[{tag}] it {it:4d}  error_2^2 {e2:9.3f} mHa^2  damage {e2 - start_ref:8.3f}  (was {2e6 * f:9.3f})  rho {rho:6.3f}  "
            f"lam {lam:.1e}  cg {ncg:3d}  tries {tries}  |g| {gnorm:.2e}  free v {int(mv.sum())}  ({time.time() - t0:.0f} s)")
        if not accepted:
            log(f"[{tag}] no step accepted; stop")
            break
        if it % save_every == 0:
            if save is not None:
                save(th, v, hist)
    return th, v, hist, it


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("state")
    ap.add_argument("rounded")
    ap.add_argument("out")
    ap.add_argument("--minutes", type=float, default=20.0)
    ap.add_argument("--iters", type=int, default=100000)
    ap.add_argument("--cg", type=int, default=300)
    ap.add_argument("--cgtol", type=float, default=1e-2)
    ap.add_argument("--lam", type=float, default=1e-2)
    ap.add_argument("--lam-min", type=float, default=1e-7)
    ap.add_argument("--free-all", action="store_true", help="control: nothing frozen (the fit itself in angle space)")
    ap.add_argument("--free-sign", action="store_true", help="no bound on the weights (the sign of a weight is free)")
    ap.add_argument("--check", action="store_true", help="finite-difference checks of the gradient, the GN product and its diagonal")
    ap.add_argument("--start-ref", type=float, default=144.49, help="fit^2 of the start point in mHa^2 (damage = error_2^2 - this)")
    ap.add_argument("--save-every", type=int, default=5)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    logf = open(os.path.join(a.out, "log.txt"), "a")

    def log(*x):
        s = " ".join(str(y) for y in x)
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    t0 = time.time()
    refc = L.Reference()
    z = np.load(a.state)
    edges = np.load(a.rounded)["tree_edges"]
    th = np.pi * np.array(z["theta"], np.float64)
    w0 = np.array(z["w"], np.float64)
    frozen = np.array(z["frozen"], bool)
    if a.free_all:
        frozen = np.zeros_like(frozen)
    R, B, C = w0.shape
    sign = np.sign(w0.sum(axis=1, keepdims=True))
    sign[sign == 0] = 1.0
    v = w0 * sign
    nneg = int((v < 0).sum())
    if not a.free_sign:
        v = np.maximum(v, 0.0)
    mod = Model(refc, edges, (R, B, C), sign)
    e_in = refc.error2_np(decode(th, edges), w0)
    f = mod.value(th, v)
    log(f"[lmhold] {a.state}: shape {(R, B, C)}, frozen {frozen.mean():.4f}, state error_2 {1e3 * e_in:.4f} mHa = "
        f"{1e6 * e_in ** 2:.2f} mHa^2 (damage {1e6 * e_in ** 2 - a.start_ref:.2f}); {nneg} minority weights clipped -> "
        f"{2e6 * f:.2f} mHa^2  ({time.time() - t0:.0f} s)")

    if a.check:
        rng = np.random.default_rng(0)
        mod.setup(th, v)
        f0, (gth, gv) = mod.grad()
        for k in range(3):
            dth = rng.standard_normal(th.shape) * 1e-3
            dvv = rng.standard_normal(v.shape) * 1e-5
            h = 1e-4
            fp = mod.value(th + h * dth, v + h * dvv)
            fm = mod.value(th - h * dth, v - h * dvv)
            an = float(np.sum(gth * dth) + np.sum(gv * dvv))
            log(f"[check] gradient: finite difference {(fp - fm) / (2 * h):.8e}  analytic {an:.8e}")
            # GN: d' GN d = ||J d||^2 with J d = At^T dA + dA^T At
            dA = mod.push(dth, dvv)
            Mx = mod.At.T @ dA
            jd2 = float(np.sum((Mx + Mx.T) ** 2))
            g1, g2 = mod.gn(dth, dvv)
            log(f"[check] GN: ||J d||^2 {jd2:.8e}  d' GN d {float(np.sum(g1 * dth) + np.sum(g2 * dvv)):.8e}")
            # push against finite differences of amap
            Ap = mod.amap(th + h * dth, v + h * dvv)
            Am = mod.amap(th - h * dth, v - h * dvv)
            log(f"[check] push: |fd - dA| / |dA| {np.linalg.norm((Ap - Am) / (2 * h) - dA) / np.linalg.norm(dA):.2e}")
        dgt, dgw = mod.diag()
        for k in range(4):
            e1, e2 = np.zeros_like(th), np.zeros_like(v)
            if k < 2:
                idx = (int(rng.integers(th.shape[0])), int(rng.integers(N - 1)))
                e1[idx] = 1.0
                ex = dgt[idx]
            else:
                idx = (int(rng.integers(R)), int(rng.integers(B)), int(rng.integers(C)))
                e2[idx] = 1.0
                ex = dgw[idx]
            g1, g2 = mod.gn(e1, e2)
            log(f"[check] diagonal {idx}: formula {ex:.8e}  e' GN e {float(np.sum(g1 * e1) + np.sum(g2 * e2)):.8e}")
        ts = time.time()
        for _ in range(5):
            mod.gn(gth, gv)
        log(f"[check] one GN product {(time.time() - ts) / 5:.3f} s")
        return

    th, v, hist, it = lm_run(mod, th, v, frozen, log, minutes=a.minutes - (time.time() - t0) / 60.0, iters=a.iters, cg=a.cg,
                             cgtol=a.cgtol, lam=a.lam, lam_min=a.lam_min, free_sign=a.free_sign, start_ref=a.start_ref,
                             save=lambda th_, v_, h_: np.savez(os.path.join(a.out, "last.npz"), theta=th_ / np.pi, w=sign * v_,
                                                               frozen=frozen, hist=np.array(h_)), save_every=a.save_every)
    np.savez(os.path.join(a.out, "last.npz"), theta=th / np.pi, w=sign * v, frozen=frozen, hist=np.array(hist))
    e_out = refc.error2_np(decode(th, edges), sign * v)
    log(f"[lmhold] DONE {it} steps: error_2 {1e3 * e_out:.4f} mHa = {1e6 * e_out ** 2:.2f} mHa^2, damage "
        f"{1e6 * e_out ** 2 - a.start_ref:.2f} (state: {1e6 * e_in ** 2 - a.start_ref:.2f})  ({time.time() - t0:.0f} s)")
    json.dump(dict(state=a.state, steps=it, start=1e6 * e_in ** 2, end=1e6 * e_out ** 2, start_ref=a.start_ref,
                   seconds=time.time() - t0), open(os.path.join(a.out, "result.json"), "w"))


if __name__ == "__main__":
    main()
