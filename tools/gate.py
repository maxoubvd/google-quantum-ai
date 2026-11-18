"""Gate checks per task, with a result cache. Run inside WSL with the venv that has kickmix + kickmix-engine.

    python tools/gate.py invsqrt FILE.kmb [--salts N] [--official] [--proof-tree TREE]
    python tools/gate.py dqi FILE.kmb [--salts N] [--official]
    python tools/gate.py dfthcsos DIR
    python tools/gate.py all DIR [--salts N] [--official] [--proof-tree TREE]   # DIR as for tools/package.py

What runs (fast first; the official verifiers stay the authority at a gate, CLAUDE.md):
- invsqrt: kickmix-engine full verification (64 shots x 3,133,452 inputs, 15 s) at the default salt and
  `--salts N - 1` extra salts; `--proof-tree` runs Syed Asad Ali's exact full-domain proof (`prove_any.py`, every x
  and every HMR transcript, ~13 s; valid for circuits of that generator family, whose integer phase model it
  checks); `--official` adds the official verifier at the default salt (3 min).
- dqi: kickmix-engine at the default salt and N - 1 extra salts (16-18 s and 0.55 GB each); reports value and
  phase errors per salt, the cost per salt, the mean (expected) cost and the worst; `--official` adds the official
  verifier at the default salt (2 min, 4.7 GB with 2 workers).
- dfthcsos: the official verifier (45 s; the engine brings no speed-up there).
- all: the three above, `tools/package.py --check-only`, and the projected score (geometric mean of the
  reference/cost ratios, dqi at its mean over the salts and at its worst).

Extra salts are named `private-1`, `private-2`, ... (the names used since S1.5). Results are cached in
`~/work/gate-cache/` keyed by (task, file SHA-256, salt, tool, tool version), so a repeated gate costs nothing.
Exit code 0 when every check passes, 1 otherwise. Run time: the sum of the runs above (minus cache hits).
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import math
import os
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VERIFIERS = os.path.join(REPO, "google-verifiable-quantum-circuit-compilation", "verifiers")
CACHE = os.path.expanduser("~/work/gate-cache")
PY = sys.executable
ENGINE = os.path.join(os.path.dirname(PY), "kickmix-engine")
PY313 = os.path.expanduser("~/venvs/gqai313/bin/python")
REFERENCE = {"invsqrt": 2076774.233715, "dqi": 24439463583672.93, "dfthcsos": 3355599588.3155375}


def sha256_path(path: str) -> str:
    """SHA-256 of a file, or of a directory as the hash of its sorted (name, file hash) list."""
    h = hashlib.sha256()
    if os.path.isdir(path):
        for f in sorted(os.listdir(path)):
            h.update(f.encode() + b"\0" + sha256_path(os.path.join(path, f)).encode())
        return h.hexdigest()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def engine_version() -> str:
    try:
        from importlib.metadata import version
        return version("kickmix-engine")
    except Exception:
        return "unknown"


def parse_json(text: str) -> dict:
    return json.loads(text[text.index("{"):text.rindex("}") + 1])


def cached_run(key: dict, cmd: list[str], env: dict | None = None, cwd: str | None = None) -> dict:
    """Run cmd (prints a JSON summary) unless (key) is cached. Adds 'seconds' and 'cached'."""
    os.makedirs(CACHE, exist_ok=True)
    name = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()[:24]
    path = os.path.join(CACHE, name + ".json")
    if os.path.exists(path):
        with open(path) as f:
            out = json.load(f)
        out["cached"] = True
        return out
    t0 = time.perf_counter()
    p = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=cwd)
    try:
        out = parse_json(p.stdout)
    except ValueError:
        out = {"status": "ERROR", "error": (p.stderr or p.stdout)[-800:]}
    out["seconds"] = round(time.perf_counter() - t0, 2)
    out["key"] = key
    if out.get("status") in ("PASS", "FAIL") or str(out.get("status", "")).startswith("PASS"):
        with open(path, "w") as f:  # only definite verdicts are cached, never tool errors
            json.dump(out, f)
    out["cached"] = False
    return out


def salts(n: int) -> list[str]:
    return [""] + [f"private-{i}" for i in range(1, max(n, 1))]


def engine_env() -> dict:
    return dict(os.environ, KICKMIX_ENGINE_VERIFIERS=VERIFIERS,
                KICKMIX_ENGINE_CACHE=os.path.expanduser("~/work/kecache"))


def gate_invsqrt(path: str, n_salts: int, official: bool, proof_tree: str | None) -> dict:
    sha = sha256_path(path)
    runs = []
    for s in salts(n_salts):
        key = {"task": "invsqrt", "sha": sha, "salt": s, "tool": "kickmix-engine", "v": engine_version()}
        r = cached_run(key, [ENGINE, "verify", "invsqrt", path, "--salt", s], env=engine_env())
        runs.append(r)
    res = {"task": "invsqrt", "file": path, "sha256": sha, "runs": runs}
    if proof_tree:
        out_json = os.path.join(CACHE, f"proof-{sha[:16]}.json")
        key = {"task": "invsqrt-proof", "sha": sha, "tool": "prove_any", "tree": os.path.basename(os.path.dirname(proof_tree))}
        cmd = ["bash", "-c", f'cd "{proof_tree}" && source tools/env.sh >/dev/null && "{PY313}" tools/macrun.py '
               f'tools/prove_any.py --circuit "{os.path.abspath(path)}" --out "{out_json}" >/dev/null 2>&1; cat "{out_json}"']
        p = cached_run(key, cmd)
        res["proof"] = {k: p.get(k) for k in ("status", "proved_inputs", "every_HMR_transcript_proved", "scratch_zero",
                                              "phase_exactly_matches_frozen_integer_model", "seconds", "cached", "error")}
    if official:
        key = {"task": "invsqrt", "sha": sha, "salt": "", "tool": "official"}
        res["official"] = cached_run(key, [PY, os.path.join(VERIFIERS, "invsqrt", "verify_submission_invsqrt.py"), path],
                                     cwd=VERIFIERS)
    statuses = [r.get("status") for r in runs] + ([res["official"].get("status")] if official else [])
    proof_ok = (not proof_tree) or str(res["proof"].get("status", "")).startswith("PASS")
    res["pass"] = all(s == "PASS" for s in statuses) and proof_ok
    res["cost"] = max((r.get("score") or math.inf) for r in runs)
    return res


def gate_dqi(path: str, n_salts: int, official: bool) -> dict:
    sha = sha256_path(path)
    runs = []
    for s in salts(n_salts):
        key = {"task": "dqi", "sha": sha, "salt": s, "tool": "kickmix-engine", "v": engine_version()}
        runs.append(cached_run(key, [ENGINE, "verify", "dqi", path, "--salt", s], env=engine_env()))
    res = {"task": "dqi", "file": path, "sha256": sha, "runs": runs}
    if official:
        key = {"task": "dqi", "sha": sha, "salt": "", "tool": "official"}
        res["official"] = cached_run(key, [PY, os.path.join(VERIFIERS, "dqi", "verify_submission_dqi.py"), path,
                                           "--max-workers", "2"], cwd=VERIFIERS)
    costs = [r.get("score") for r in runs if r.get("status") == "PASS"]
    statuses = [r.get("status") for r in runs] + ([res["official"].get("status")] if official else [])
    res["pass"] = all(s == "PASS" for s in statuses)
    res["errors_per_salt"] = [(r.get("value_errors"), r.get("phase_errors")) for r in runs]
    res["salts_with_errors"] = sum(1 for r in runs if (r.get("value_errors") or 0) + (r.get("phase_errors") or 0))
    res["cost_mean"] = sum(costs) / len(costs) if costs and res["pass"] else math.inf
    res["cost_worst"] = max(costs) if costs and res["pass"] else math.inf
    res["cost"] = res["cost_mean"]
    return res


def gate_dfthcsos(path: str) -> dict:
    sha = sha256_path(path)
    key = {"task": "dfthcsos", "sha": sha, "salt": "", "tool": "official"}
    r = cached_run(key, [PY, os.path.join(VERIFIERS, "dfthcsos", "verify.py"), os.path.abspath(path)], cwd=VERIFIERS)
    return {"task": "dfthcsos", "file": path, "sha256": sha, "runs": [r], "pass": r.get("status") == "PASS",
            "cost": r.get("score") or math.inf}


def score(costs: dict) -> float:
    return math.prod(REFERENCE[t] / costs[t] for t in REFERENCE) ** (1 / 3)


def report(res: dict) -> None:
    t = res["task"]
    print(f"== {t}: {'PASS' if res['pass'] else 'FAIL'}  sha {res['sha256'][:12]}  cost {res['cost']:.6g}")
    for r in res["runs"]:
        salt = (r.get("key") or {}).get("salt", "")
        extra = f" errors {r.get('value_errors')}/{r.get('phase_errors')}" if t == "dqi" else ""
        print(f"   salt {salt or '(default)':10s} {r.get('status')} cost {r.get('score')}{extra} "
              f"[{r.get('seconds')} s{', cached' if r.get('cached') else ''}]" + (f" {r.get('error')}" if r.get("error") else ""))
    if "proof" in res:
        print(f"   exact proof: {res['proof']}")
    if "official" in res:
        o = res["official"]
        print(f"   official verifier: {o.get('status')} cost {o.get('score')} [{o.get('seconds')} s{', cached' if o.get('cached') else ''}]")
    if t == "dqi":
        print(f"   mean cost {res['cost_mean']:.6g}, worst {res['cost_worst']:.6g}, salts with errors {res['salts_with_errors']}/{len(res['runs'])}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("task", choices=["invsqrt", "dqi", "dfthcsos", "all"])
    ap.add_argument("path")
    ap.add_argument("--salts", type=int, default=4, help="number of salts (default salt + N-1 extra); default 4")
    ap.add_argument("--official", action="store_true", help="also run the official verifier at the default salt")
    ap.add_argument("--proof-tree", help="invsqrt: relocated build tree of Syed Asad Ali's generator (exact proof)")
    ap.add_argument("--json", help="write the full results to this file")
    a = ap.parse_args()
    t0 = time.perf_counter()
    results = []
    if a.task in ("invsqrt", "all"):
        results.append(gate_invsqrt(a.path if a.task == "invsqrt" else os.path.join(a.path, "invsqrt.kmb"),
                                    a.salts, a.official, a.proof_tree))
    if a.task in ("dqi", "all"):
        results.append(gate_dqi(a.path if a.task == "dqi" else os.path.join(a.path, "dqi.kmb"), a.salts, a.official))
    if a.task in ("dfthcsos", "all"):
        results.append(gate_dfthcsos(a.path if a.task == "dfthcsos" else os.path.join(a.path, "dfthcsos")))
    for r in results:
        report(r)
    ok = all(r["pass"] for r in results)
    if a.task == "all":
        p = subprocess.run([PY, os.path.join(REPO, "tools", "package.py"), "--check-only", a.path], capture_output=True, text=True)
        pkg_ok = p.returncode == 0
        print(f"== package checks: {'PASS' if pkg_ok else 'FAIL'}" + ("" if pkg_ok else f" {p.stdout[-300:]}{p.stderr[-300:]}"))
        ok = ok and pkg_ok
        if ok:
            c = {r["task"]: r["cost"] for r in results}
            worst = dict(c, dqi=next(r for r in results if r["task"] == "dqi")["cost_worst"])
            print(f"== projected score: {score(c):.5f} (dqi at its mean over the salts), {score(worst):.5f} (dqi worst salt)")
    print(f"== {'ALL PASS' if ok else 'FAILED'} in {time.perf_counter() - t0:.1f} s")
    if a.json:
        with open(a.json, "w") as f:
            json.dump(results, f, indent=1, default=str)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
