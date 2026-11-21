#!/usr/bin/env python
"""r48_chem_q7c_vtab (design seat R48-CHEM-Q7C): independent check of a finished table with the hosts' own functions,
in a fresh process: verify.load_inputs(TABLE.npz) -> verify.evaluate_words(inputs, angle_words, tree_edges, brot) ->
verify.passes_accuracy; then the formula's price lambda_eff x (fixed(R, B, C, brot) + 2 L) x Q for a loader count L
per pass (the loaders themselves are C2's: this is the price IF the loaders cost L).

    python r48_chem_q7c_vtab.py TABLE.npz [--brot 7] [--loader 4045] [--json OUT.json]
Also prints: the range of the words, the share of trivial two-body words (0 or a quarter turn), D = the hard-sign
defect (share of sum |w| carried by the minority sign of each channel, summed in the units of the fit tools).
"""
import argparse
import json
import os
import sys

import numpy as np

DATA = os.environ.get("QCC_DATA", os.path.expanduser("~/qcc/data"))
sys.path[:0] = [os.path.join(DATA, "kickmix"), os.path.join(DATA, "verifiers")]
os.environ.setdefault("DFTHC_DATA", DATA)
from dfthcsos import verify as V  # noqa: E402

ENTRY = 289164662.05


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--brot", type=int, default=7)
    ap.add_argument("--loader", type=int, default=4045)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    with np.load(a.table, allow_pickle=False) as z:
        d = {k: np.array(z[k]) for k in z.files}
    words, edges = np.asarray(d["angle_words"], np.int64), np.asarray(d["tree_edges"], np.int64)
    R, B, C = d["W"].shape
    N = edges.shape[0] + 1
    assert words.shape == (N + R * B, N - 1), words.shape
    assert words.min() >= 0 and words.max() < (1 << a.brot), (words.min(), words.max())
    seen = {0}
    for p, c in edges:          # a tree rooted at 0, parent before child
        assert int(p) in seen and int(c) not in seen, (p, c)
        seen.add(int(c))
    inputs = V.load_inputs({k: d[k] for k in ("U", "W", "H2_I", "H1_I", "Hsym_lower")})
    m = V.evaluate_words(inputs, words, edges, a.brot)
    ok = V.passes_accuracy(m)
    out = {k: float(m[k]) for k in ("error_1", "error_2", "Lambda", "E_gap", "lambda_eff")}
    h = V.heuristic_costs(N, R, B, C, a.brot)
    T = h["toffoli"] + 2 * a.loader
    Q = h["persistent_qubits"] + 5 * a.brot + max(h["temporary_qubits"], 77)
    cost = out["lambda_eff"] * T * Q
    w = np.asarray(d["W"], np.float64)
    pos, neg = np.clip(w, 0, None).sum(axis=1), np.clip(-w, 0, None).sum(axis=1)
    minority = float(np.minimum(pos, neg).sum() / np.abs(w).sum())
    q = 1 << (a.brot - 2)
    triv = float(np.mean(np.mod(words[N:], q) == 0))
    out.update(shape=[R, B, C], brot=a.brot, passes_accuracy=bool(ok), fixed_toffoli=int(h["toffoli"]),
               loader_per_pass=a.loader, T=int(T), Q=int(Q), cost=cost, vs_entry_pct=100 * (cost / ENTRY - 1),
               minority_sign_share=minority, trivial_two_body_words=triv, table=os.path.abspath(a.table))
    print(f"{a.table}: shape {(R, B, C)} brot {a.brot}  error_1 {1e3 * out['error_1']:.4f} mHa  error_2 "
          f"{1e3 * out['error_2']:.4f} mHa  E_gap {out['E_gap']:.6f}  Lambda {out['Lambda']:.4f}  lambda_eff "
          f"{out['lambda_eff']:.6f}  passes_accuracy {ok}")
    print(f"  price at {a.loader} loader Toffolis per pass: T {T} (fixed {h['toffoli']}) x Q {Q} x lambda_eff = "
          f"{cost:,.0f} ({out['vs_entry_pct']:+.2f}% against the entry {ENTRY:,.0f}); minority-sign share of sum|w| "
          f"{minority:.2e}; trivial two-body words {100 * triv:.1f}%")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
