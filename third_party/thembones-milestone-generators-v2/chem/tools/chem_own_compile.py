#!/usr/bin/env python
"""chem_own_compile: the fifteen ANGLE_k lookup circuits of a dfthcsos entry, from the angle-word table alone.

    search  <table> <plan_dir> [--effort N] [--seed S] [--workers 3] [--circuits 0,3,7] [--resume]
        Plan every circuit (span search of tools/chem_own_engine.py) and store one JSON per circuit in
        <plan_dir> (relation spans, price, history). Deterministic for given (table, seed, effort): a rerun
        writes the same plans. --resume continues from stored plans with a new round of the same search.
    build   <table> <plan_dir> <out_dir> [--workers 3] [--shots 4] [--verify] [--npz submission.npz]
        Re-derive the products from the table and the stored spans, emit angle_k_0..14.kmb, and check them:
        every wire on every address classically (at emission), then the verifier's own simulation
        (`--shots` per address), then with --verify the official verifier on salts "" and "alt-salt-1".
        A circuit whose stored plan is missing or worse than the old engine is compiled by the old engine
        (tools/chem_loaders_core.py, `decompose`), so the result is never worse than chem_loaders_compile.
    price   <table> <plan_dir>
        Table of per-circuit prices of the stored plans.

<table> is a submission.npz holding `angle_words` (rows x 75, values 0..255) and `tree_edges`, or a .npy of the
words. Nothing else is read: no circuit file, no plan from anybody else.

Credits. Circuit scheme and emitter: tools/chem_loaders_core.py (ours; the scheme is that of the public loaders
of the mrmeows / marculera / Syed lineage, re-derived by dissection in ws_chem_loaders). Pairing engine:
tools/chem_own_engine.py, which implements two ideas described (without code) in marculera's notebook
"DF-THC 11.2x: Codes on the XOR of Five Words" (plans priced exactly by linear matroid parity; shrink the
spare classes and swap at fixed class count until the pairing is perfect again) and takes maoxin123's
"five rank-86 blocks" (notebook "LB 25.303 | Combined invsqrt + refit chem") as its yardstick. The notebook
ideas that concern the TABLE (codes on the XOR of five words, the four word forms) are not implemented here:
this compiler takes the words as given. "Five dirty words in 128 qubits" (two low-cache phases for circuit 0)
was already in chem_loaders_core and is used unchanged.
"""
from __future__ import annotations

import argparse
import hashlib
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
VERSION = 5


def load_words(path):
    if path.endswith(".npy"):
        return np.load(path).astype(np.int64), None
    with np.load(path, allow_pickle=False) as z:
        return np.array(z["angle_words"]).astype(np.int64), np.array(z["tree_edges"]).astype(np.int64)


def words_sha(words):
    return hashlib.sha256(np.ascontiguousarray(words.astype(np.int64)).tobytes()).hexdigest()


def feasible_lbs(words):
    """Address splits whose high cache can fit at all, smallest block first (as in chem_loaders_compile)."""
    A = words.shape[0]
    out = []
    for lb in (3, 4, 5):
        H = -(-A // (1 << lb))
        hi = sum(1 for M in range(1, H) if M & (M - 1))
        if hi <= core.QUBIT_CAP - max(1, (A - 1).bit_length()) - 40 + 5:
            out.append(lb)
    return out


def analyze_first(words, k, lbs, relax=False):
    """Analysis of circuit k at the first address split that fits. relax: the don't-care addresses of the last
    partial block are used to shrink the block row spaces (chem_own_engine.relax_last_row); a plan records
    whether it was searched that way, and is always rebuilt the same way."""
    errs = []
    for lb in lbs:
        try:
            an = core.analyze(words, k, lb)
        except ValueError as exc:
            errs.append(str(exc))
            continue
        if relax:
            import chem_own_engine as eng
            an.relaxed = eng.relax_last_row(an)
        return an, lb
    raise ValueError(f"circuit {k}: no address split fits: {'; '.join(errs)}")


def plan_path(plan_dir, k):
    return os.path.join(plan_dir, f"plan_k{k}.json")


def _search_job(args):
    words, k, lbs, seed, effort, plan_dir, resume, relax = args
    import chem_own_engine as eng
    t0 = time.time()
    path = plan_path(plan_dir, k)
    sha = words_sha(words)
    rounds = 0
    prev = None
    relax_round = 1 if relax else 0          # the first round that used the last-block don't-cares (0: never)
    if resume and os.path.exists(path):
        prev = json.load(open(path))
        if prev.get("words_sha256") == sha and prev.get("version") == VERSION:
            rounds = prev["rounds"]
            if int(prev.get("relax", 0)):    # once relaxed, every later round is relaxed
                relax_round, relax = int(prev["relax"]), True
            elif relax:                      # an unrelaxed plan continued on the relaxed instance: its
                relax_round = rounds + 1     # relations annihilate the smaller row spaces too
        else:
            prev = None
    an, lb = analyze_first(words, k, lbs, relax)
    if prev is not None and prev.get("lb") != lb:
        prev, rounds, relax_round = None, 0, (1 if relax else 0)
    engine = None
    if prev is not None:
        engine = eng.Engine(an, seed + 1000 * rounds)
        engine.restore([[int(x, 16) for x in R] for R in prev["spans"]])
    if an.X0.shape[1] == 0:
        out = dict(k=k, lb=lb, price=0, spans=[], codim=[], rounds=1, hist=[], seconds=0.0)
    else:
        E, price, hist = eng.search(an, seed + 1000 * rounds, effort, engine=engine)
        if prev is not None and prev["price"] < price:      # a new round never replaces a better plan
            E.restore([[int(x, 16) for x in R] for R in prev["spans"]])
            price = prev["price"]
        out = dict(k=k, lb=lb, price=int(price), codim=E.codim, spans=[[hex(v) for v in R] for R in E.R],
                   rounds=rounds + 1, hist=(prev["hist"] if prev else []) + [list(h) for h in hist],
                   evals=E.evals, trials=E.trials, exacts=E.exacts, mismatch=E.mismatch,
                   seconds=round(time.time() - t0, 1))
    out.update(seed=seed, effort=effort, infra=an.infra, words_sha256=sha, version=VERSION, relax=relax_round)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(out, fh)
    os.replace(tmp, path)
    return dict(k=k, price=out["price"], infra=an.infra, codim=out.get("codim"), seconds=out["seconds"],
                rounds=out["rounds"])


def search(words, plan_dir, effort=1, seed=0, workers=3, circuits=None, resume=False, relax=False):
    os.makedirs(plan_dir, exist_ok=True)
    lbs = feasible_lbs(words)
    ks = list(range(K)) if circuits is None else circuits
    jobs = [(words, k, tuple(lbs), seed, effort, plan_dir, resume, relax) for k in ks]
    if workers <= 1:
        res = [_search_job(j) for j in jobs]
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(min(3, workers)) as ex:
            res = list(ex.map(_search_job, jobs, chunksize=1))
    return res


def _build_job(args):
    words, k, lbs, plan_dir, old_effort = args
    import chem_own_engine as eng
    t0 = time.time()
    path = plan_path(plan_dir, k) if plan_dir else None
    plan = json.load(open(path)) if path and os.path.exists(path) else None
    an, lb = analyze_first(words, k, lbs, bool(plan and plan.get("relax", 0)))
    log = []
    prods = None
    source = None
    if plan is not None and an.X0.shape[1]:
        if plan.get("words_sha256") == words_sha(words) and plan.get("lb") == lb:
            E = eng.Engine(an, plan.get("seed", 0))
            E.restore([[int(x, 16) for x in R] for R in plan["spans"]])
            try:
                prods = E.plan()
                source = "chem_own_engine"
                log += E.log
            except RuntimeError as exc:
                log.append(f"stored plan not extracted: {exc}")
        else:
            log.append("stored plan belongs to another table or split: ignored")
    if an.X0.shape[1] == 0:
        prods, source = [], "none needed"
    if old_effort or prods is None:
        old, olog = core.decompose(an, 0, max(1, old_effort), "auto")
        if prods is None or len(old) < len(prods):
            prods, source = old, "chem_loaders_core"
            log += olog
    text, st = core.emit(an, prods)
    st.update(k=k, lb=lb, source=source, log=log, seconds=round(time.time() - t0, 1), kmx=text)
    return st


def build(words, plan_dir, out_dir, workers=3, npz=None, shots=4, verify=False, old_effort=2):
    sys.path[:0] = [os.path.join(DATA, "kickmix")]
    import kickmix as km
    t0 = time.time()
    lbs = feasible_lbs(words)
    jobs = [(words, k, tuple(lbs), plan_dir, old_effort) for k in range(K)]
    if workers <= 1:
        res = [_build_job(j) for j in jobs]
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(min(3, workers)) as ex:
            res = list(ex.map(_build_job, jobs, chunksize=1))
    os.makedirs(out_dir, exist_ok=True)
    for r in res:
        c = km.Circuit(r.pop("kmx"))
        assert c.num_qubits <= core.QUBIT_CAP, (r["k"], c.num_qubits)
        st = c.operation_statistics().num_toffoli_like
        assert st == r["toffoli"], f"circuit {r['k']}: emitted {st} Toffolis, planned {r['toffoli']}"
        c.to_file(os.path.join(out_dir, f"angle_k_{r['k']}.kmb"), format="kmb")
    out = dict(per_pass=sum(r["toffoli"] for r in res), circuits=res, rows=int(words.shape[0]),
               words_sha256=words_sha(words))
    if npz:
        shutil.copyfile(npz, os.path.join(out_dir, "submission.npz"))
    if shots:
        sys.path[:0] = [os.path.join(DATA, "verifiers")]
        from dfthcsos.scoreinfra import circuits as vc
        circ = vc.load_circuits(out_dir)
        got, layout, tof = vc.recover_words_from_circuits(circ, num_vectors=words.shape[0],
                                                          num_angles=words.shape[1], shots_per_address=shots,
                                                          seed=12345)
        assert np.array_equal(got, words), "circuits load different words"
        assert abs(tof - 2 * out["per_pass"]) < 1e-9, (tof, out["per_pass"])
        out["check"] = dict(words_match=True, toffoli_forward_and_reverse=float(tof),
                            selection_bits=len(layout.selection), max_qubits=max(c.num_qubits for c in circ),
                            shots_per_address=shots)
    if verify:
        sys.path[:0] = [os.path.join(DATA, "verifiers")]
        from dfthcsos import verify as v
        out["verify"] = {}
        for salt in ("", "alt-salt-1"):
            full = v.verify(out_dir, salt=salt)
            with open(os.path.join(out_dir, f"verify_{'salt0' if salt == '' else salt}.json"), "w") as fh:
                json.dump(full, fh, indent=1, default=str)
            out["verify"][salt] = {x: full[x] for x in ("status", "score", "lambda_eff", "toffoli_count",
                                                         "num_qubits", "verified_shots", "submission_sha256")}
            out["verify"][salt]["error_1"] = full["metrics"]["error_1"]
            out["verify"][salt]["error_2"] = full["metrics"]["error_2"]
    out["seconds"] = round(time.time() - t0, 1)
    return out


def show_plans(plan_dir):
    tot = 0
    print(" k  lb infra products total  codim              rounds  seconds")
    for k in range(K):
        path = plan_path(plan_dir, k)
        if not os.path.exists(path):
            print(f"{k:2d}  (no plan)")
            continue
        p = json.load(open(path))
        tot += p["infra"] + p["price"]
        print(f"{k:2d}  {p['lb']}  {p['infra']:4d}  {p['price']:5d}   {p['infra'] + p['price']:4d}  "
              f"{str(p.get('codim')):18s} {p['rounds']:3d}   {p['seconds']}"
              f"{'  relaxed from round ' + str(p['relax']) if p.get('relax') else ''}")
    print("per pass (stored plans):", tot)
    return tot


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["search", "build", "price"])
    ap.add_argument("table")
    ap.add_argument("plan_dir")
    ap.add_argument("out", nargs="?")
    ap.add_argument("--effort", type=int, default=1)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--circuits", default=None)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--relax", action="store_true",
                    help="search: use the don't-care addresses of the last partial block; the plan records the "
                         "first round that did (with --resume an unrelaxed plan is continued relaxed)")
    ap.add_argument("--npz", default=None)
    ap.add_argument("--shots", type=int, default=4)
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--old-effort", type=int, default=2)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()
    words, _ = load_words(args.table)
    if args.mode == "search":
        ks = [int(x) for x in args.circuits.split(",")] if args.circuits else None
        res = search(words, args.plan_dir, args.effort, args.seed, args.workers, ks, args.resume, args.relax)
        for r in res:
            print(f"k={r['k']:2d} products {r['price']} total {r['price'] + r['infra']} codim {r['codim']} "
                  f"round {r['rounds']} ({r['seconds']} s)")
        show_plans(args.plan_dir)
    elif args.mode == "price":
        show_plans(args.plan_dir)
    else:
        npz = args.npz or (args.table if args.table.endswith(".npz") else None)
        out = build(words, args.plan_dir, args.out, args.workers, npz, args.shots, args.verify, args.old_effort)
        print(f"rows {out['rows']}  per pass {out['per_pass']}  ({out['seconds']} s)")
        print("per-circuit Toffolis:", [r["toffoli"] for r in out["circuits"]])
        print("per-circuit qubits:  ", [r["qubits"] for r in out["circuits"]])
        print("engines:", sorted({r["source"] for r in out["circuits"]}))
        if "check" in out:
            print("simulation check:", out["check"])
        for salt, v in out.get("verify", {}).items():
            print(f"verifier salt {salt!r}: {json.dumps(v)}")
        if args.json:
            with open(args.json, "w") as fh:
                json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
