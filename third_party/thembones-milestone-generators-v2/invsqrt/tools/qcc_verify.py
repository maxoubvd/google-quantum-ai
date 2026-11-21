#!/usr/bin/env python
"""Local scorer: run the official verifiers on a submission directory and print costs,
per-task ratios against the hosts' reference and the leaderboard-style geometric mean.

Submission layout: <dir>/invsqrt.kmb, <dir>/dqi.kmb, <dir>/dfthcsos/{angle_k_*.kmb, submission.npz}.

Usage (inside WSL, with the 3.13 venv):
  python qcc_verify.py <dir> [--tasks invsqrt,dqi,chem] [--salt S] [--quick] [--workers N]
--quick uses 4 shots per x for invsqrt and 10,240 shots for dqi (not comparable to the board).
"""
import argparse
import json
import math
import os
import sys
import time

DATA = os.environ.get("QCC_DATA", os.path.expanduser("~/qcc/data"))
sys.path[:0] = [os.path.join(DATA, "kickmix"), os.path.join(DATA, "verifiers")]

REF = {"invsqrt": 2076774.233715282, "dqi": 24439463583672.93, "chem": 3355599588.3155365}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("subdir")
    ap.add_argument("--tasks", default="invsqrt,dqi,chem")
    ap.add_argument("--salt", default="")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()
    import kickmix as km

    out = {}
    tasks = args.tasks.split(",")
    if "invsqrt" in tasks:
        from invsqrt import verify_submission_invsqrt as v
        t0 = time.time()
        c = km.Circuit.from_file(os.path.join(args.subdir, "invsqrt.kmb"), format="kmb")
        spec = v.ProblemSpecInvSqrt.for_n_p(10, shots_per_x=4 if args.quick else 64)
        r = v.verify(c, spec, salt=args.salt, batch_size=16 if args.quick else 256, max_workers=args.workers)
        r["seconds"] = round(time.time() - t0, 1)
        out["invsqrt"] = r
        print(f"invsqrt  cost {r['score']:.2f}  Q {r['num_qubits']}  Tof {r['toffoli_count']:.2f}  "
              f"T {r['t_count']:.1f}  rot {r['rotation_count']:.1f}  worst {r['worst_error_over_tolerance']:.5f} "
              f"@x={r['worst_x']}  [{r['seconds']}s]", flush=True)
    if "dqi" in tasks:
        from dqi import verify_submission_dqi as v
        t0 = time.time()
        c = km.Circuit.from_file(os.path.join(args.subdir, "dqi.kmb"), format="kmb")
        spec = v.ProblemSpecDQI.benchmark_instance(num_shots=10240 if args.quick else v.DEFAULT_NUM_SHOTS)
        r = v.verify(c, spec, salt=args.salt, max_workers=args.workers or v.DEFAULT_MAX_WORKERS)
        r.pop("shot_operation_counts", None)
        r["seconds"] = round(time.time() - t0, 1)
        out["dqi"] = r
        print(f"dqi      cost {r['score']:.6g}  Q {r['num_qubits']}  Tof {r['toffoli_count']:.1f}  "
              f"errors {r['value_errors']}/{r['phase_errors']}  R {r['repetition_factor']:.5f}  [{r['seconds']}s]", flush=True)
        del c
    if "chem" in tasks:
        from dfthcsos import verify as v
        t0 = time.time()
        r = v.verify(os.path.join(args.subdir, "dfthcsos"), salt=args.salt)
        r["seconds"] = round(time.time() - t0, 1)
        out["chem"] = r
        m = r["metrics"]
        print(f"chem     cost {r['score']:.6g}  lambda_eff {r['lambda_eff']:.5f}  T {r['toffoli_count']:.1f}  "
              f"Q {r['num_qubits']}  err1 {1e3 * m['error_1']:.3f}  err2 {1e3 * m['error_2']:.3f} mHa  [{r['seconds']}s]", flush=True)
    ratios = {k: REF[k] / out[k]["score"] for k in out}
    for k, x in ratios.items():
        print(f"  ratio {k}: {x:.5f}")
    if len(ratios) == 3:
        print(f"GEOMETRIC MEAN: {math.prod(ratios.values()) ** (1 / 3):.5f}" + ("  (quick mode: not board-comparable)" if args.quick else ""))
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
