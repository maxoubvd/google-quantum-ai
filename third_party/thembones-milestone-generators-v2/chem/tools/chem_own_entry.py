#!/usr/bin/env python
"""chem_own_entry: one command from a table of angle words to a complete, verified dfthcsos entry.

    chem_own_entry.py <table.npz> <label> [--weights refit|scratch|keep] [--effort N] [--seed S]
                      [--rounds N] [--workers 3] [--work DIR] [--publish DIR] [--quick-refit] [--no-publish]
                      [--only weights] [--spec plan_spec.json]

Steps (every one is our code; inputs are the table and the hosts' verifier data, nothing else):
  1. weights   tools/chem_own_refit.py: W, H2_I for the fixed words.
               refit   start from the W / H2_I stored in the table and refit them (default);
               scratch ignore the stored weights: least squares from the decoded vectors, then the same refit;
               keep    take the stored weights unchanged (only for a table whose weights are already ours).
  2. plans     tools/chem_own_compile.py search: relation spans for the 15 circuits (deterministic in
               (words, seed, effort, rounds)); plans are stored in <work>/plans and reused if they belong to
               the same words, seed and effort.
               --spec gives seed and rounds per circuit (the recipe written by tools/chem_own_pick.py when
               the best plan of each circuit is taken from several seeds).
  3. circuits  tools/chem_own_compile.py build: angle_k_0..14.kmb, checked on every address.
  4. verify    the official verifier, full shots, salt "" and salt "alt-salt-1".
  5. publish   <publish>/best_<label>/: the 16 submission files, both verifier JSON files, MANIFEST.json with
               the sha256 of every file, the submission hash, the per-circuit counts and the exact command.

<table.npz>: `angle_words`, `tree_edges`, and `W`, `H2_I`, `U` (shape (R, B, C) and, unless --weights scratch,
the starting weights). See tools/chem_own_compile.py and tools/chem_own_refit.py for credits.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chem_own_compile as coc  # noqa: E402

WS = os.path.expanduser("~/qcc/ws/chem_own")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table")
    ap.add_argument("label")
    ap.add_argument("--weights", default="refit", choices=["refit", "scratch", "keep"])
    ap.add_argument("--effort", type=int, default=1)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--rounds", type=int, default=1)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--work", default=None)
    ap.add_argument("--publish", default=WS)
    ap.add_argument("--quick-refit", action="store_true")
    ap.add_argument("--no-publish", action="store_true")
    ap.add_argument("--old-effort", type=int, default=2)
    ap.add_argument("--only", default="all", choices=["all", "weights"],
                    help="weights: run step 1 and stop (the full command reuses its result)")
    ap.add_argument("--relax", action="store_true", help="plan search with the last-block don't-cares")
    ap.add_argument("--spec", default=None,
                    help="plan_spec.json of chem_own_pick.py: seed and rounds per circuit (overrides --seed/--rounds)")
    a = ap.parse_args()
    t0 = time.time()
    work = a.work or os.path.join(WS, "build", a.label)
    os.makedirs(work, exist_ok=True)
    table = os.path.abspath(os.path.expanduser(a.table))
    words, edges = coc.load_words(table)
    wsha = coc.words_sha(words)
    print(f"[entry] table {table}: {words.shape[0]} rows, words sha256 {wsha[:16]}...", flush=True)

    # 1. weights (reused when this work directory already holds the refit of the same table and mode)
    npz = os.path.join(work, "submission.npz")
    meta_path = os.path.join(work, "weights_meta.json")
    meta = dict(table_sha256=sha256_file(table), weights=a.weights, quick=bool(a.quick_refit))
    if a.weights == "keep":
        shutil.copyfile(table, npz)
        print("[entry] weights: kept from the table", flush=True)
    elif os.path.exists(npz) and os.path.exists(meta_path) and json.load(open(meta_path)) == meta:
        print("[entry] weights: refit of this table found in the work directory, reused", flush=True)
    else:
        cmd = [sys.executable, os.path.join(HERE, "chem_own_refit.py"), table, npz,
               "--start", "table" if a.weights == "refit" else "scratch"]
        if a.quick_refit:
            cmd.append("--quick")
        print("[entry] weights:", " ".join(cmd), flush=True)
        with open(os.path.join(work, "refit.log"), "w") as fh:
            rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
        tail = open(os.path.join(work, "refit.log")).read().strip().splitlines()[-2:]
        print("[entry]   " + "\n[entry]   ".join(tail), flush=True)
        if rc:
            raise SystemExit(f"refit failed (exit {rc}); see {work}/refit.log")
        with open(meta_path, "w") as fh:
            json.dump(meta, fh)
    if a.only == "weights":
        print(f"[entry] --only weights: done ({time.time() - t0:.0f} s)", flush=True)
        return

    # 2. plans
    plan_dir = os.path.join(work, "plans")
    os.makedirs(plan_dir, exist_ok=True)

    # the recipe: (seed, rounds) per circuit; --spec (written by chem_own_pick.py) or --seed / --rounds for all
    recipe = {k: (a.seed, a.rounds, int(a.relax)) for k in range(coc.K)}
    spec = None
    if a.spec:
        spec = json.load(open(a.spec))
        a.effort = spec.get("effort", a.effort)
        if spec.get("words_sha256") not in (None, wsha):
            raise SystemExit("the plan spec belongs to another table")
        for k, c in spec["circuits"].items():
            recipe[int(k)] = (int(c["seed"]), int(c["rounds"]), int(c.get("relax", 0)))

    def plan_state(k):
        """Rounds of the stored plan of circuit k if it is a prefix of the recipe, else 0."""
        path = coc.plan_path(plan_dir, k)
        if not os.path.exists(path):
            return 0
        p = json.load(open(path))
        pr, rr = int(p.get("relax", 0)), recipe[k][2]       # first relaxed round of the plan / of the recipe
        same = (p.get("words_sha256") == wsha and p.get("seed") == recipe[k][0] and p.get("effort") == a.effort
                and p.get("version") == coc.VERSION and p["rounds"] <= recipe[k][1]
                and (pr == rr or (pr == 0 and rr > p["rounds"])))
        return p["rounds"] if same else 0

    for seed in sorted({s for s, _, _ in recipe.values()}):
        ks = [k for k in range(coc.K) if recipe[k][0] == seed]
        for rnd in range(max(recipe[k][1] for k in ks)):
            for relax in (False, True):      # a circuit is relaxed from the round its recipe says
                todo = [k for k in ks if plan_state(k) <= rnd < recipe[k][1]
                        and (0 < recipe[k][2] <= rnd + 1) == relax]
                if not todo:
                    continue
                print(f"[entry] plan search seed {seed} round {rnd + 1}{' (relaxed)' if relax else ''}: circuits "
                      f"{todo} (effort {a.effort})", flush=True)
                for k in todo:               # a stale plan (other table / seed / effort / longer run) is not resumed
                    if plan_state(k) == 0 and os.path.exists(coc.plan_path(plan_dir, k)):
                        os.remove(coc.plan_path(plan_dir, k))
                res = coc.search(words, plan_dir, a.effort, seed, a.workers, todo, resume=rnd > 0, relax=relax)
                print("[entry]   " + " ".join(f"k{r['k']}:{r['price'] + r['infra']}" for r in res), flush=True)
    coc.show_plans(plan_dir)

    # 3 + 4. circuits, simulation check, official verifier on both salts
    out = coc.build(words, plan_dir, work, a.workers, None, 4, True, a.old_effort)
    per = [r["toffoli"] for r in out["circuits"]]
    print(f"[entry] circuits: per pass {out['per_pass']}  {per}", flush=True)
    print(f"[entry] engines: {sorted({r['source'] for r in out['circuits']})}; check {out['check']}", flush=True)
    for salt, v in out["verify"].items():
        print(f"[entry] verifier salt {salt!r}: {json.dumps(v)}", flush=True)
    ok = all(v["status"] == "PASS" for v in out["verify"].values())

    # 5. publish
    manifest = dict(
        label=a.label, created=time.strftime("%Y-%m-%d %H:%M:%S"), table=table, words_sha256=wsha,
        command=" ".join([os.path.basename(sys.argv[0])] + sys.argv[1:]),
        weights=a.weights, effort=a.effort, compiler_version=coc.VERSION,
        plan_recipe={str(k): dict(seed=recipe[k][0], rounds=recipe[k][1], relax=recipe[k][2])
                     for k in range(coc.K)},
        lookups_per_pass=out["per_pass"], per_circuit_toffoli=per,
        per_circuit_qubits=[r["qubits"] for r in out["circuits"]],
        per_circuit_engine=[r["source"] for r in out["circuits"]],
        verify=out["verify"], all_pass=ok, seconds=round(time.time() - t0, 1))
    files = [f"angle_k_{k}.kmb" for k in range(coc.K)] + ["submission.npz", "verify_salt0.json",
                                                         "verify_alt-salt-1.json"]
    manifest["sha256"] = {f: sha256_file(os.path.join(work, f)) for f in files}
    with open(os.path.join(work, "MANIFEST.json"), "w") as fh:
        json.dump(manifest, fh, indent=1)
    if ok and not a.no_publish:
        dest = os.path.join(a.publish, f"best_{a.label}")
        os.makedirs(dest, exist_ok=True)
        for f in files + ["MANIFEST.json"]:
            shutil.copyfile(os.path.join(work, f), os.path.join(dest, f))
        pl = os.path.join(dest, "plans")
        os.makedirs(pl, exist_ok=True)
        for k in range(coc.K):
            if os.path.exists(coc.plan_path(plan_dir, k)):
                shutil.copyfile(coc.plan_path(plan_dir, k), coc.plan_path(pl, k))
        if os.path.exists(os.path.join(work, "refit.log")) and a.weights != "keep":
            shutil.copyfile(os.path.join(work, "refit.log"), os.path.join(dest, "refit.log"))
        if a.spec:
            shutil.copyfile(a.spec, os.path.join(dest, "plan_spec.json"))
        print(f"[entry] published {dest}", flush=True)
    s0 = out["verify"][""]
    print(f"[entry] {'PASS' if ok else 'FAIL'}: cost {s0['score']:.6e} = lambda_eff {s0['lambda_eff']:.6f} x "
          f"T {s0['toffoli_count']:.0f} x Q {s0['num_qubits']}; lookups {out['per_pass']} per pass; "
          f"submission sha256 {s0['submission_sha256']}  ({time.time() - t0:.0f} s)", flush=True)
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
