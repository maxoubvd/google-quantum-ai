"""chem_rk7_whold (RK7 holder 2): the WEIGHTS-ONLY LM of a finished table (rounded.npz form), on the GPU.
Every word frozen; R8's LM step (tools/chem_own2_lmstage7.StageLM7, imported unchanged: hard sign, Nielsen damping,
backtracking) re-solves the weights; prints error_2^2 after every chunk; writes the same rounded.npz with the new W.
Measured motive [M R8 row 42-45]: an END's weights are not converged (r8w_d: -1.27 mHa^2 after 220 steps).
    python chem_rk7_whold.py IN.npz OUT.npz [--steps 300] [--chunk 20] [--seconds 1800] [--stop 0.002]
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp")
    ap.add_argument("out")
    ap.add_argument("--steps", type=int, default=300)
    ap.add_argument("--chunk", type=int, default=20)
    ap.add_argument("--seconds", type=float, default=1800.0)
    ap.add_argument("--stop", type=float, default=0.002)
    a = ap.parse_args()
    t0 = time.time()
    from chem_rk7_round import Target
    from chem_own2_lmstage7 import StageLM7
    z = np.load(a.inp, allow_pickle=True)
    keep = {k: z[k] for k in z.files}
    words = np.asarray(keep["words"], np.int64)
    signed = np.where(words >= 64, words - 128, words)
    assert signed.min() >= -32 and signed.max() <= 31, "words outside the signed range"
    th = signed.astype(np.float64) / 64.0                      # units of pi (theta = pi * word / 64)
    W = np.asarray(keep["W"], np.float64)
    R, B, C = W.shape
    edges = np.asarray(keep["tree_edges"], np.int64)
    refc = Target()
    early = 0.5 / 64.0
    lm = StageLM7(refc.T, refc.i, refc.j, edges, (R, B, C), -0.5 - early, 0.5 - early, cache_dir="jax_cache")
    if "allowed" in keep:       # t10 seat 6 (the lead's 16:20Z): the table's zero mask (True = a free weight)
        lm.wmask = np.asarray(keep["allowed"], bool)
        print(f"[whold] zero mask: {int((~lm.wmask).sum())} masked weights, max |W| on them "
              f"{float(np.abs(W[~lm.wmask]).max()):.1e} at the start", flush=True)
    fr = np.ones(th.shape, bool)
    wall = np.zeros(th.shape, bool)
    e0 = 1e6 * refc.error2_np(keep["U"], W) ** 2
    print(f"[whold] {a.inp}: error_2^2 {e0:.4f} mHa^2 at its weights", flush=True)
    done, info = 0, None
    while done < a.steps and time.time() - t0 < a.seconds:
        th2, W, info = lm.refit(th, W, fr, wall, steps=a.chunk, seconds=1e9, stop=-1.0)
        done += info["steps"]
        print(f"[whold] steps {done}: {info['e2_out']:.4f} mHa^2 (chunk {info['e2_out'] - info['e2_in']:+.4f}, "
              f"{info['rejected']} rejected, {info['seconds']:.0f} s)", flush=True)
        if info["steps"] < a.chunk or (info["e2_in"] - info["e2_out"]) < a.stop * a.chunk / 3.0:
            break
    e1 = 1e6 * refc.error2_np(keep["U"], W) ** 2
    assert np.all(np.sign(W.sum(axis=1))[:, None, :] * W >= -1e-15), "a weight crossed its channel's sign"
    if lm.wmask is not None:    # t10 seat 6: the mask holds exactly
        assert np.all(W[~lm.wmask] == 0.0), "a masked weight moved"
        print(f"[whold] zero mask held: max |W| on the masked weights {float(np.abs(W[~lm.wmask]).max()):.1e}", flush=True)
    keep["W"] = W
    tmp = a.out + ".tmp.npz"
    np.savez(tmp, **keep)
    os.replace(tmp, a.out)
    print(f"[whold] wrote {a.out}: error_2^2 {e0:.4f} -> {e1:.4f} mHa^2 ({e1 - e0:+.4f}) in {done} steps, "
          f"{time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
