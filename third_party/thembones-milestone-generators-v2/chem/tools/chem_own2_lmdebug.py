"""chem_own2_lmdebug (R8 holder 6, Nov 16 22:35Z): why does the in-stage LM reject every try on some states? Runs the
LM step of chem_own2_lmstage on one kept state at several lam and prints f, f2, pred, the step sizes and NaN flags.

  python chem_own2_lmdebug.py STATE.npz EDGES.npz
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chem_own2_lmround import Target  # noqa: E402
from chem_own2_lmstage import StageLM  # noqa: E402


def main():
    z = np.load(sys.argv[1])
    th, w, fr = np.array(z["theta"], np.float64), np.array(z["w"], np.float64), np.array(z["frozen"], bool)
    edges = np.array(np.load(sys.argv[2])["edges"], np.int64)
    refc = Target()
    lm = StageLM(refc.T, refc.i, refc.j, edges, w.shape, -0.5 - 0.5 / 64, 0.5 - 0.5 / 64, cache_dir="jax_cache")
    jax = lm.jax
    with jax.enable_x64(True):
        import jax.numpy as jnp
        lm._build()
        sign_np = np.sign(w.sum(axis=1, keepdims=True))
        sign_np[sign_np == 0] = 1.0
        v_np = np.maximum(w * sign_np, 0.0)
        print(f"weights: minority-sign mass {np.abs(w[w * sign_np < 0]).sum():.3e} of {np.abs(w).sum():.3e}; channels with "
              f"sum 0: {int((w.sum(axis=1) == 0).sum())}; v==0 share {(v_np == 0).mean():.4f}", flush=True)
        thj, vj, sj = jnp.asarray(np.pi * th), jnp.asarray(v_np), jnp.asarray(sign_np)
        mth = jnp.asarray((~fr).astype(np.float64)) * (0.0 if os.environ.get("R8_DBG_WONLY") else 1.0)
        wl = jnp.zeros_like(mth)
        for lam in (1e-3, 1e-1, 10.0, 1e3):
            f, f2, th2, v2, pred, ncg, gnorm = lm._lm_iter(thj, vj, sj, lam, mth, wl)
            d1 = np.asarray(th2) - np.asarray(thj)
            d2 = np.asarray(v2) - np.asarray(vj)
            print(f"lam {lam:.0e}: f {2e6 * float(f):.4f} f2 {2e6 * float(f2):.4f} pred {2e6 * float(pred):.4e} cg {int(ncg)} "
                  f"|g| {float(gnorm):.3e} max|d_theta| {np.nanmax(np.abs(d1)):.3e} max|d_v| {np.nanmax(np.abs(d2)):.3e} "
                  f"nan {bool(np.isnan(d1).any() or np.isnan(d2).any())}", flush=True)
            if os.environ.get("R8_DBG_LINE"):      # f along the step: th + t d1, max(v + t d2, 0)
                row = []
                for t in (1.0, 0.3, 0.1, 0.03, 0.01, 0.003):
                    ft = float(lm._fval(thj + t * jnp.asarray(d1), jnp.maximum(vj + t * jnp.asarray(d2), 0.0), sj))
                    row.append(f"t {t}: {2e6 * (ft - float(f)):+.4f}")
                print("   line: " + "  ".join(row), flush=True)


if __name__ == "__main__":
    main()
