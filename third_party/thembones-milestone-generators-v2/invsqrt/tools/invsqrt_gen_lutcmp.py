"""invsqrt_gen: compare lookup programs with the block specs of the last local build
(~/qcc/ws/invsqrt_gen/lut/) and, on request, install the cheapest fitting one per block (own tool).

Sources of programs (all read here, none written):
  ~/qcc/ws/lut_synth/out/            programs delivered by the lookup workstream
  ~/qcc/ws/invsqrt_gen/lut/ext/      programs made with its generator by tools/sh/invsqrt_gen_ext.sh, programs
                                     of table sets (ties_<set>/, adapted ones in ties_<set>/trim/)
  ~/qcc/ws/invsqrt_gen/lut/cache/    programs of OUR OWN solver from searches at any depth and fill count (the
                                     same format). Installed as a program file, such a program is replayed
                                     whatever `deep` the build line says, so every block can have its own
                                     best depth and this PC (no C kernel) never has to search.
usage: invsqrt_gen_lutcmp.py            list
       invsqrt_gen_lutcmp.py install    copy, for every block, the cheapest program that fits the current spec
                                        to lut/<block>.prog.json, unless the program the last build made itself
                                        (lut/<block>.own.json) is at least as cheap; an installed program that
                                        is worse than that is removed
       invsqrt_gen_lutcmp.py clear      remove lut/<block>.prog.json files (builds use our own solver again)
       add the word "quiet" to print only the result lines.
A program is replayed through the generator's own truth-table bookkeeping at build time (invsqrt_gen_pool.py),
so a wrong program stops the build; nothing is trusted here beyond the spec hash.
"""
from __future__ import annotations

import glob
import json
import os
import sys

LUT = os.path.expanduser("~/qcc/ws/invsqrt_gen/lut")
SOURCES = [os.path.expanduser("~/qcc/ws/lut_synth/out"), os.path.join(LUT, "ext"), os.path.join(LUT, "cache")]


def n_and(prog):
    return sum(1 for s in prog["steps"] if s[0] == "and")


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "list"
    quiet = "quiet" in sys.argv
    if mode == "clear":
        for path in glob.glob(os.path.join(LUT, "*.prog.json")):
            os.remove(path)
            print("removed", path)
        return
    specs = {}
    for path in sorted(glob.glob(os.path.join(LUT, "*.spec.json"))):
        d = json.load(open(path))
        own = None
        try:
            o = json.load(open(path.replace(".spec.json", ".own.json")))
            if o.get("spec_sha") == d["spec_sha"]:
                own = n_and(o)
        except OSError:
            pass
        specs[d["block"]] = (d["spec_sha"], own)
        print(f"spec {d['block']:6s} sha {d['spec_sha']}  outputs {len(d['outputs'])}  scratch {d['max_scratch']}  "
              f"program of the last build's own solver: {own} ANDs")
    best = {}
    for src in SOURCES:
        for path in sorted(glob.glob(os.path.join(src, "*.json")) + glob.glob(os.path.join(src, "*", "*.json"))
                           + glob.glob(os.path.join(src, "*", "*", "*.json"))):      # e.g. ext/ties_<set>/trim/
            try:
                d = json.load(open(path))
            except (OSError, ValueError):
                continue
            if not isinstance(d, dict) or "steps" not in d:
                continue
            blk, sha = d.get("block"), d.get("spec_sha")
            fit = blk in specs and specs[blk][0] == sha
            if not quiet and (fit or src != SOURCES[2]):
                print(f"prog {path.replace(os.path.expanduser('~'), '~')}: block {blk} sha {sha} ANDs {n_and(d)}  "
                      f"fits current spec: {fit}")
            if fit and (blk not in best or n_and(d) < best[blk][0]):
                best[blk] = (n_and(d), path, src == SOURCES[2])
    for blk, (sha, own) in sorted(specs.items()):
        ip = os.path.join(LUT, blk + ".prog.json")
        try:
            cur = json.load(open(ip))
            cur_n = n_and(cur) if cur.get("spec_sha") == sha else None
        except (OSError, ValueError):
            cur_n = None
        if blk not in best:
            continue
        n, path, ours = best[blk]
        short = path.replace(os.path.expanduser("~"), "~")
        if own is not None and n >= own:
            print(f"{blk}: best program {n} is not below the last build's own {own}")
            if mode == "install" and cur_n is not None and cur_n > own:
                os.remove(ip)
                print(f"{blk}: installed program ({cur_n} ANDs) removed")
        elif mode == "install":
            if cur_n is not None and cur_n <= n:
                print(f"{blk}: installed program has {cur_n} ANDs already (best available {n})")
                continue
            prog = json.load(open(path))
            prog["origin"] = ("own solver: " + os.path.basename(path)) if ours else short
            with open(ip, "w") as fh:
                json.dump(prog, fh)
            print(f"{blk}: installed {n} ANDs (last build's own {own}) from {short}")
        else:
            print(f"{blk}: {n} ANDs available (last build's own {own}; installed {cur_n}) in {short}")


if __name__ == "__main__":
    main()
