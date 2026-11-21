"""r50_chem_lbl_obchain: THE ONE-BODY ROWS UNDER DELTA WIRES IN CHAINS (R50-CHEM-LBL, construction 9b).
progressive_dyn of r50_chem_lbl_obwires.py with DynChain (r50_chem_lbl_chain.py): FREE's progressive rounding, the
address of a row decided at freeze time, one free bit per chain; a word is forced to a parity only when its relation
is already decided by the row's earlier words.
"""
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r48_chem_q7c_onebody import STAGES, angles_and_mass, decode_jac, decode_pi  # noqa: E402
from r49_chem_free_obround import near_par  # noqa: E402
from r50_chem_lbl_chain import DynChain  # noqa: E402


def progressive_dyn(th0, edges, lam_at, H, scale, iters, wires):
    nrow = th0.shape[0]
    dyn = DynChain(wires, np.arange(nrow), nrow)
    theta = th0.copy()
    frozen = np.zeros(theta.shape, bool)
    qv = np.zeros(theta.shape)
    lam2 = lam_at[:, None] ** 2
    ncol = theta.shape[1]

    def fun(x, mask, base):
        t = base.copy()
        t[mask] = x
        V, J = decode_jac(t, edges)
        D = (V.T * lam_at) @ V - H
        gV = 4.0 * lam_at[:, None] * (V @ D)
        return float(np.sum(D * D)), np.einsum("jn,jkn->jk", gV, J)[mask]

    for frac in STAGES:
        need = int(round(frac * theta.size)) - int(frozen.sum())
        if need > 0:
            now = np.where(frozen, qv, theta)
            mass = angles_and_mass(decode_pi(now, edges), edges)[1]
            xs = now * scale
            q = np.round(xs)
            impact = lam2 * mass * (q / scale - now) ** 2
            impact[frozen] = -np.inf
            pick = np.argsort(-impact, axis=None)[:need]
            for flat in pick:
                j, c = divmod(int(flat), ncol)
                zc = int(q[j, c])
                p = dyn.ask(j, c)
                if p is not None:
                    zc = int(near_par(xs[j, c], p))
                dyn.commit(j, c, zc)
                frozen[j, c] = True
                qv[j, c] = zc / scale
        free = ~frozen
        if free.any():
            base = np.where(frozen, qv, theta)
            r = minimize(fun, theta[free], args=(free, base), jac=True, method="L-BFGS-B",
                         options={"maxiter": iters, "maxcor": 20, "ftol": 1e-14, "gtol": 1e-10})
            theta = base.copy()
            theta[free] = r.x
        else:
            theta = qv.copy()
    return np.round(qv * scale).astype(np.int64), dyn
