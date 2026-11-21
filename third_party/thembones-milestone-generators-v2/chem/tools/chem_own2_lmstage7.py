#!/usr/bin/env python
"""chem_own2_lmstage7 (R8 holder 7, Nov 17 01:40Z): the in-stage LM of tools/chem_own2_lmstage.py with a cheaper
handling of a step that does not descend.

Measured defect of the present step [M r8w_d / r8s_r1 / r8l_r1 logs, head stages 5-40]: 7 to 13 of about 23 tries per
stage are REJECTED, and every rejected try costs a full PCG solve (about 1 s) for nothing: 15 accepted steps in 25 s.
Two changes, nothing else (the arithmetic of a step is StageLM's compiled lm_iter, unchanged):
  * a rejected step is first BACKTRACKED along its own direction (the convex combination x + t (x_new - x),
    t = 0.5 then 0.25: inside the walls and v >= 0 by convexity; one function value each, no PCG); the first t that
    descends is taken;
  * the damping follows Nielsen's rule: on a full step lam *= max(1/3, 1 - (2 rho - 1)^3), nu = 2; on a rejected
    (or backtracked) step lam *= nu, nu *= 2. (StageLM: /3 above rho 0.6, x2 below 0.25, x4 on a rejection.)
Same call and return as StageLM.refit; info gains "backtracked".
"""
from __future__ import annotations

import time

import numpy as np

from chem_own2_lmstage import StageLM


class StageLM7(StageLM):
    def refit(self, theta, w, frozen, wall, steps=25, seconds=25.0, stop=0.01, lam0=0.0, lam_start=0.0, stop_lam=0.0,
              ts=(0.5, 0.25)):
        jax = self.jax
        t0 = time.time()
        with jax.enable_x64(True):
            import jax.numpy as jnp
            if not self._built:
                self._build()
            sign_np = np.sign(w.sum(axis=1, keepdims=True))
            sign_np[sign_np == 0] = 1.0
            v_np = np.maximum(w * sign_np, 0.0)
            if self.wmask is not None:      # t10 seat 6 (the lead's 16:20Z): masked weights exactly 0 on the way in
                v_np = np.where(self.wmask, v_np, 0.0)
            th = jnp.asarray(np.pi * np.asarray(theta, np.float64))
            v = jnp.asarray(v_np)
            sign = jnp.asarray(sign_np)
            mth = jnp.asarray((~frozen).astype(np.float64))
            wl = jnp.asarray((wall & ~frozen).astype(np.float64))
            f_in = float(self._fval(th, v, sign))
            f_cur, gains, it, cgs, rejected, back = f_in, [], 0, [], 0, 0
            lam = lam_start if lam_start > 0 else max(self.lam, lam0)
            nu = 2.0
            while it < steps and (time.time() - t0) < seconds:
                tries, accepted = 0, False
                while not accepted and tries < 8:
                    tries += 1
                    f, f2, th2, v2, pred, ncg, gnorm = self._lm_iter(th, v, sign, lam, mth, wl)
                    f, f2, pred = float(f), float(f2), float(pred)
                    cgs.append(int(ncg))
                    if f2 < f:
                        rho = (f - f2) / max(pred, 1e-300)
                        accepted = True
                        th, v = th2, v2
                        lam = lam * max(1.0 / 3.0, 1.0 - (2.0 * rho - 1.0) ** 3)
                        nu = 2.0
                    else:
                        for t_ in ts:           # backtrack along the same step: function values only
                            th_t, v_t = th + t_ * (th2 - th), v + t_ * (v2 - v)
                            f_t = float(self._fval(th_t, v_t, sign))
                            if f_t < f:
                                th, v, f2, accepted = th_t, v_t, f_t, True
                                back += 1
                                break
                        if not accepted:
                            rejected += 1
                        lam *= nu
                        nu *= 2.0
                    lam = max(lam, self.lam_min)
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
            if self.wmask is not None:      # t10 seat 6: and exactly 0 on the way out
                w_out = np.where(self.wmask, w_out, 0.0)
        th_out = np.where(frozen, theta, th_out)            # frozen words exactly as they were
        return th_out, w_out, {"e2_in": 2e6 * f_in, "e2_out": 2e6 * f_cur, "steps": it, "rejected": rejected,
                               "backtracked": back, "seconds": time.time() - t0,
                               "cg_mean": float(np.mean(cgs)) if cgs else 0.0, "lam": self.lam}
