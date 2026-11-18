"""Targeted dqi test: run decoders on shots whose error weight is forced, with the official sampler and simulator.

    python tools/dqi_targeted.py CIRCUIT.kmb [CIRCUIT2.kmb ...] --weights 14 15 16 17 18 19 --shots 2048 [--seed 1]

For each weight k, the shots are drawn by the official `sample_shots` with the weight distribution replaced by a
one-hot at k (positions, error values and the classical MBUC record are drawn exactly as in the official test),
then passed to the official `verify(..., shots=...)`, which simulates them with kickmix and counts value and phase
errors. The same shots (same seed per weight) are given to every circuit, so the circuits are compared on
identical inputs. Purpose (idea DQI-A03): confirm which error weights a decoder with a given extractor cutoff no
longer decodes, and that it fails on nothing else.

Run inside WSL with the venv that has kickmix (`~/venvs/gqai`). Cost: one simulation of `--shots` shots per
weight and circuit (about 2–4 s for 2,048 shots of the 455 MB decoder with 2 workers, plus ~5 s to load it).
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time

import numpy as np

VERIFIER_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "google-verifiable-quantum-circuit-compilation", "verifiers", "dqi")
sys.path.insert(0, VERIFIER_DIR)
import kickmix as km  # noqa: E402
import verify_submission_dqi as V  # noqa: E402


def shots_for_weight(k: int, n_shots: int, seed: int):
    spec = V.ProblemSpecDQI.benchmark_instance(num_shots=n_shots)
    tables = V.FieldTables(km.GF2Field(spec.degree))
    probs = np.zeros(spec.max_errors + 1)
    probs[k] = 1.0
    rng = np.random.default_rng([seed, k])
    return spec, V.sample_shots(spec, tables, rng, n_shots, probs)


def run(circuit, spec, shots, workers: int) -> dict:
    try:
        r = V.verify(circuit, spec, shots=shots, max_workers=workers)
        return {"value_errors": r["value_errors"], "phase_errors": r["phase_errors"], "toffolis": r.get("toffoli_count")}
    except AssertionError as e:  # the official check fails when the error rate is >= 1/16; read the counts back
        m = re.search(r"value_errors=(\d+)/\d+, phase_errors=(\d+)/\d+", str(e))
        if not m:
            raise
        return {"value_errors": int(m.group(1)), "phase_errors": int(m.group(2)), "toffolis": None}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("circuits", nargs="+")
    ap.add_argument("--weights", type=int, nargs="+", default=[14, 15, 16, 17, 18, 19, 20])
    ap.add_argument("--shots", type=int, default=2048, help="shots per weight (multiple of 256)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--json")
    a = ap.parse_args()
    t0 = time.perf_counter()
    data = {k: shots_for_weight(k, a.shots, a.seed) for k in a.weights}
    results = {}
    for path in a.circuits:
        c = km.Circuit.from_file(path)
        row = {}
        for k, (spec, shots) in data.items():
            row[k] = run(c, spec, shots, a.workers)
        results[path] = row
        del c
    print(f"shots per weight: {a.shots} (seed {a.seed}); errors = value + phase")
    print("weight | " + " | ".join(os.path.basename(os.path.dirname(p)) or p for p in a.circuits))
    for k in a.weights:
        cells = []
        for p in a.circuits:
            r = results[p][k]
            cells.append(f"{r['value_errors'] + r['phase_errors']:5d} ({r['value_errors']}v/{r['phase_errors']}p)")
        print(f"{k:6d} | " + " | ".join(cells))
    print(f"done in {time.perf_counter() - t0:.1f} s")
    if a.json:
        with open(a.json, "w") as f:
            json.dump({p: {str(k): v for k, v in r.items()} for p, r in results.items()}, f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
