#!/usr/bin/env python
"""chem_loaders_compile: ANGLE_k lookup circuits for a dfthcsos angle-word table (word width 8, J = 5).

Usage (inside WSL, ~/qcc/.venv/bin/python):

  price a table (no circuits written; exact Toffoli counts of what `build` would emit):
      chem_loaders_compile.py price  <table> [--lb auto|3|4|5] [--effort N] [--workers 3] [--json out.json]
  build the circuits:
      chem_loaders_compile.py build  <table> <out_dir> [--lb ...] [--effort N] [--workers 3] [--verify]
                                     [--npz submission.npz]
  build and check by simulation only (tables without a valid factorization, e.g. synthetic shapes):
      chem_loaders_compile.py check  <table> <out_dir> [--lb ...] [--effort N] [--shots 4]
  search exact form changes of the words (theta, theta+pi, -theta, pi-theta) for a cheaper table:
      chem_loaders_compile.py gauge  <table> <out_words.npy> [--rounds 1] [--effort 2]
  options: --engine auto|blocks|groups (default auto: both, keep the cheaper per circuit), --seed S

<table> is a submission.npz holding `angle_words` (rows x 75, values 0..255) and `tree_edges`, or a .npy file
of the words. The shape is read from the table: rows = 76 + R*B, selection bits = ceil(log2(rows)), blocks of
2**lb addresses. `build` writes angle_k_0..14.kmb into <out_dir> and copies the submission.npz next to them
(`--npz` when the table is a .npy); `--verify` runs the official verifier on salts "" and "alt-salt-1".

From Python (the rounding loop of the factorization workstream):

    sys.path.insert(0, ".../tools"); import chem_loaders_compile as clc
    res = clc.price_table(words, lb="auto", effort=4, workers=3)   # words: (rows, 75) ints
    res["per_pass"]      # lookup Toffolis per pass = sum over the 15 circuits
    res["circuits"]      # per circuit: lb, infra, products, toffoli, phases, carriers, levels, log
    clc.price_table(words, effort=1, workers=1, engine="blocks")   # fastest exact price: about 1 s per table
    clc.model_table(words)                                         # instant model estimate, not exact

The price is exact because pricing runs the same deterministic decomposition that `build` emits (same seed,
same effort, same lb, same engine); `build` asserts that the emitted Toffoli count equals it and `check` /
`--verify` assert that the simulated count equals it. T_total of the verifier is its formula part (6,662 at
(8, 143, 44)) plus 2 x per_pass.

Measured on the public (8, 143, 44) table (1,220 rows), 3 workers unless noted:
    effort 1, engine blocks, 1 worker   1.1 s   4,187 per pass
    effort 2 (default engines)        10-19 s   4,161
    effort 8                          15-19 s   4,159   (shipped circuits: 4,133; older public compiler: 4,613)
    model_table                       < 1 s     4,176   (estimate)
Compare tables at the same effort; effort changes the count by up to about 0.6%.

What the compiler reads from the table (nothing is passed in):
  * shape: rows -> selection bits; address split lb (blocks of 2**lb addresses) = the smallest whose caches fit
    under 128 qubits per circuit (lb 3 for 648 rows, 4 for 791..1,220, see R = 9 below);
  * carriers: any word wire, or bit7 xor bit6 of a word, that holds the same function of the low address bits
    in the batch before and after a circuit; a carrier of a new function saves one AND, any other carrier
    becomes a free workspace qubit;
  * coded bits: output combinations that are affine or quadratic in the low address bits inside every block
    (single bit columns or XORs of columns), found by elimination;
  * rows past the table are don't-cares.
R = 9 (1,363 rows) needs 78 cached high monomials for 77 workspace qubits: it compiles with lb 4 only when the
bit-7 wires are free (words in the signed gauge, bit 7 = bit 6: about 4,550 per pass on the synthetic test);
with carriers in bit 7 it falls back to lb 5 (about 6,900). Keep R <= 8 or the signed gauge at R = 9.

Method and credits: tools/chem_loaders_core.py and research/ws_chem_loaders.md. The scheme is the one of the
shipped public loaders (mrmeows / marculera / Syed lineage), re-derived by dissection; all code here is ours.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402

DATA = os.environ.get("QCC_DATA", os.path.expanduser("~/qcc/data"))
K = 15


def load_words(path):
    if path.endswith(".npy"):
        return np.load(path).astype(np.int64), None
    with np.load(path, allow_pickle=False) as z:
        return np.array(z["angle_words"]).astype(np.int64), np.array(z["tree_edges"]).astype(np.int64)


def _job(args):
    """One circuit; args[2] is a list of address splits tried in order (the first that fits is used)."""
    errs = []
    for lb in args[2]:
        try:
            return core.compile_circuit(args[:2] + (lb,) + args[3:])
        except ValueError as exc:      # this address split does not fit under the qubit cap
            errs.append(str(exc))
    return dict(k=args[1], error="; ".join(errs))


def _run(words, lbs, seed, effort, workers, want_text, engine="auto"):
    jobs = [(words, k, tuple(lbs), seed, effort, want_text, engine) for k in range(K)]
    if workers <= 1:
        res = [_job(j) for j in jobs]
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(min(3, workers)) as ex:
            res = list(ex.map(_job, jobs, chunksize=1))
    return res


def feasible_lbs(words):
    """Address splits whose high cache can fit at all, smallest block first: more rows per product is
    cheaper, so the first split that fits a circuit is the one to use."""
    A = words.shape[0]
    out = []
    for lb in (3, 4, 5):
        H = -(-A // (1 << lb))
        hi = sum(1 for M in range(1, H) if M & (M - 1))
        if hi <= core.QUBIT_CAP - max(1, (A - 1).bit_length()) - 40 + 5:
            out.append(lb)
    return out


def price_table(words, lb="auto", effort=4, workers=3, seed=0, want_text=False, engine="auto"):
    """Exact lookup Toffolis per pass and per circuit. lb="auto" uses, per circuit, the smallest address
    block (2**lb addresses) whose caches fit under the 128-qubit cap; circuits share only the public
    registers, so different circuits may use different splits."""
    words = np.asarray(words, np.int64)
    t0 = time.time()
    lbs = feasible_lbs(words) if lb == "auto" else [int(lb)]
    res = _run(words, lbs, seed, effort, workers, want_text, engine)
    best = []
    for r in res:
        if "error" in r:
            raise ValueError(f"circuit {r['k']}: no address split fits: {r['error']}")
        best.append(r)
    return dict(per_pass=sum(r["toffoli"] for r in best), circuits=best, seconds=round(time.time() - t0, 2),
                rows=int(words.shape[0]), effort=effort, seed=seed)


def model_table(words, lb="auto"):
    """Instant estimate (no pairing search): infrastructure plus half the feature count of the block plan at
    its default span codimensions. Within about 1% of the exact price on the tables tried; not exact."""
    words = np.asarray(words, np.int64)
    lbs = feasible_lbs(words) if lb == "auto" else [int(lb)]
    tot = 0
    for k in range(K):
        for x in lbs:
            try:
                an = core.analyze(words, k, x)
            except ValueError:
                continue
            bd = core.BlockDecomp(an, 0, 1)
            bd.prepare()
            feats = sum(bi["amb"] - e for bi, e in zip(bd.binfo, bd.default_e()))
            tot += an.infra + (feats + 1) // 2
            break
        else:
            raise ValueError(f"circuit {k}: no address split fits")
    return tot


def build(words, out_dir, lb="auto", effort=4, workers=3, seed=0, npz=None, verify=False, engine="auto",
          shots=0):
    sys.path[:0] = [os.path.join(DATA, "kickmix")]
    import kickmix as km
    res = price_table(words, lb, effort, workers, seed, want_text=True, engine=engine)
    os.makedirs(out_dir, exist_ok=True)
    for r in res["circuits"]:
        c = km.Circuit(r.pop("kmx"))
        assert c.num_qubits <= core.QUBIT_CAP, (r["k"], c.num_qubits)
        st = c.operation_statistics().num_toffoli_like
        assert st == r["toffoli"], f"circuit {r['k']}: emitted {st} Toffolis, priced {r['toffoli']}"
        c.to_file(os.path.join(out_dir, f"angle_k_{r['k']}.kmb"), format="kmb")
    if npz:
        shutil.copyfile(npz, os.path.join(out_dir, "submission.npz"))
    if shots:
        # the verifier's own simulation and register checks, without the chemistry gates
        sys.path[:0] = [os.path.join(DATA, "verifiers")]
        from dfthcsos.scoreinfra import circuits as vc
        circ = vc.load_circuits(out_dir)
        got, layout, tof = vc.recover_words_from_circuits(circ, num_vectors=words.shape[0], num_angles=words.shape[1],
                                                          shots_per_address=shots, seed=12345)
        assert np.array_equal(got, words), "circuits load different words"
        assert abs(tof - 2 * res["per_pass"]) < 1e-9, (tof, res["per_pass"])
        res["check"] = dict(words_match=True, toffoli_forward_and_reverse=tof, selection_bits=len(layout.selection),
                            max_qubits=max(c.num_qubits for c in circ), shots_per_address=shots)
    if verify:
        sys.path[:0] = [os.path.join(DATA, "verifiers")]
        from dfthcsos import verify as v
        res["verify"] = {}
        for salt in ("", "alt-salt-1"):
            out = v.verify(out_dir, salt=salt)
            res["verify"][salt] = {x: out[x] for x in ("status", "score", "lambda_eff", "toffoli_count", "num_qubits",
                                                        "verified_shots", "submission_sha256")}
            res["verify"][salt]["error_2"] = out["metrics"]["error_2"]
    return res


def _show(res):
    print(f"rows {res['rows']}  per pass {res['per_pass']}  ({res['seconds']} s, effort {res['effort']})")
    print(" k  lb infra products total phases  outputs(affine,quadratic,generic)")
    for r in res["circuits"]:
        print(f"{r['k']:2d}  {r['lb']}  {r['infra']:4d}  {r['products']:5d}  {r['toffoli']:5d}  {r['phases']:3d}     {r['levels']}")
    print("per-circuit Toffolis:", [r["toffoli"] for r in res["circuits"]])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["price", "build", "check", "gauge", "model"])
    ap.add_argument("table")
    ap.add_argument("out", nargs="?")
    ap.add_argument("--lb", default="auto")
    ap.add_argument("--effort", type=int, default=4)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--npz", default=None)
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--json", default=None)
    ap.add_argument("--log", action="store_true")
    ap.add_argument("--rounds", type=int, default=1)
    ap.add_argument("--engine", default="auto", choices=["auto", "blocks", "groups"])
    ap.add_argument("--shots", type=int, default=4)
    args = ap.parse_args()
    words, edges = load_words(args.table)
    if args.mode == "model":
        print("model estimate per pass:", model_table(words, args.lb))
        return
    if args.mode == "price":
        res = price_table(words, args.lb, args.effort, args.workers, args.seed, engine=args.engine)
    elif args.mode == "build":
        npz = args.npz or (args.table if args.table.endswith(".npz") else None)
        res = build(words, args.out, args.lb, args.effort, args.workers, args.seed, npz, args.verify, args.engine)
    elif args.mode == "check":
        res = build(words, args.out, args.lb, args.effort, args.workers, args.seed, None, False, args.engine,
                    shots=args.shots)
        print("simulation check:", res["check"])
    else:
        import chem_loaders_gauge as gauge
        if edges is None:
            raise SystemExit("gauge needs tree_edges: pass a submission.npz")
        gauge.search(words, edges, args.out, rounds=args.rounds, effort=args.effort, workers=args.workers,
                     lb=4 if args.lb == "auto" else int(args.lb))
        return
    _show(res)
    if args.log:
        for r in res["circuits"]:
            print(f"-- circuit {r['k']}")
            for line in r["log"]:
                print("   ", line)
    if "verify" in res:
        for salt, out in res["verify"].items():
            print(f"verifier salt {salt!r}: {json.dumps(out)}")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(res, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
