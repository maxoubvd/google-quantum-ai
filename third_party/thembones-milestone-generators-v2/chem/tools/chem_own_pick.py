#!/usr/bin/env python
"""chem_own_pick: choose, for every circuit, the cheapest stored plan among several search runs.

    chem_own_pick.py <out_plan_dir> <plan_dir> [<plan_dir> ...]

Every <plan_dir> holds plan_k*.json files written by `chem_own_compile.py search` (one seed per directory,
any number of rounds). For each circuit the plan with the fewest products wins (ties: fewer rounds, then the
directory listed first). The winners are copied to <out_plan_dir> together with `plan_spec.json`:

    {"effort": E, "version": V, "words_sha256": ..., "circuits": {"k": {"seed": s, "rounds": r, ...}}}

which is the recipe of the plans: `chem_own_entry.py ... --spec plan_spec.json` reuses the copied plans, and on a
clean work directory recomputes exactly them (search is deterministic in (words, seed, effort, rounds)).
"""
from __future__ import annotations

import json
import os
import shutil
import sys

K = 15


def main():
    out = sys.argv[1]
    dirs = sys.argv[2:]
    os.makedirs(out, exist_ok=True)
    best = {}
    for order, d in enumerate(dirs):
        for k in range(K):
            path = os.path.join(d, f"plan_k{k}.json")
            if not os.path.exists(path):
                continue
            p = json.load(open(path))
            key = (p["price"], p["rounds"], order)
            if k not in best or key < best[k][0]:
                best[k] = (key, path, p)
    spec = dict(circuits={})
    tot = 0
    print(" k  total products seed rounds  codim               source")
    for k in range(K):
        if k not in best:
            print(f"{k:2d}  (no plan in any directory)")
            continue
        _, path, p = best[k]
        for name in ("effort", "version", "words_sha256"):
            if name in spec and spec[name] != p[name]:
                raise SystemExit(f"circuit {k}: plans disagree on {name}")
            spec[name] = p[name]
        shutil.copyfile(path, os.path.join(out, f"plan_k{k}.json"))
        spec["circuits"][str(k)] = dict(seed=p["seed"], rounds=p["rounds"], relax=int(p.get("relax", 0)),
                                        products=p["price"], total=p["price"] + p["infra"])
        tot += p["price"] + p["infra"]
        print(f"{k:2d}  {p['price'] + p['infra']:5d}  {p['price']:5d}   {p['seed']:3d}  {p['rounds']:3d}    "
              f"{str(p.get('codim')):19s} {path}")
    spec["lookups_per_pass"] = tot
    with open(os.path.join(out, "plan_spec.json"), "w") as fh:
        json.dump(spec, fh, indent=1)
    print(f"per pass: {tot}  ({len(best)} of {K} circuits)  spec: {os.path.join(out, 'plan_spec.json')}")


if __name__ == "__main__":
    main()
