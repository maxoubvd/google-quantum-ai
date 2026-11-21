#!/usr/bin/env python3
"""r59_chem_groups_g42.py: designer R59-CHEM-GROUPS (Nov 17): a COPY of builder C2's tools/chem_ld2_e19g42.py (its
driver for the 7-output group at 42) that runs the builder's group tool over the TIE engine
(tools/r59_chem_groups_grp.py: groups 1-2 planned with column ties, group 2 at 38) and first scans engine seeds
(--seeds) for one whose dump keeps the later kills (groups 39 | 38 | 43 | 36): wired loader 155 products = 229.
Everything below the seed scan is the builder's driver, unchanged.

ORIGINAL DOCSTRING: BOX side driver: ONE placed-kill loader with its
7-output groups planned at 42 where the engine reads 43.

Stage 1 (tools/chem_ld2_e19grp.py): dump the groups of the placed-kill engine at the given seed; annealers on generic
group 2's relation space (one core each, --secs each, rounds of len(cores) runs); whenever a run finds a space pricing
at 42, every explicit plan of every found space is offered to the engine with the closure (one plan per attempt of
the group); the first closed loader is kept (<out>/angle_k_<K>.kmx, emit_c4 checked) with its closing plan and the
groups as built. Stage 2 (relation loaders: group 3 is a 7-group with a column kill too): the same on group 3's
residual as left by stage 1's closing plan; both groups injected.

  python3 r59_chem_groups_g42.py --table sub.npz --k 4 --seeds 0-15 --cores 0,1 --rounds 4 --secs 600 --out g42_k4
"""
import argparse
import glob
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--seed", type=int, default=-1)
    ap.add_argument("--seeds", default="0-15", help="engine seeds scanned for a dump that keeps the later kills")
    ap.add_argument("--modes", default="reg,tie3,-,-")
    ap.add_argument("--wantg", default="39,38,43,36", help="generic groups' products the seed scan accepts")
    ap.add_argument("--plan", default="auto", help="placed-kill plan; auto = the plan of the loader's generic count")
    ap.add_argument("--P", type=int, default=42)
    ap.add_argument("--cores", required=True)
    ap.add_argument("--rounds", type=int, default=4)
    ap.add_argument("--secs", type=float, default=600)
    ap.add_argument("--rseed0", type=int, default=100)
    ap.add_argument("--nostage2", action="store_true")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    log = os.path.join(a.out, "driver.log")
    grp = os.path.join(HERE, "r59_chem_groups_grp.py")
    env = dict(os.environ, CHEM_BROT="7", CHEM_HB="77", OMP_NUM_THREADS="1", R59_MODES=a.modes)
    cores = a.cores.split(",")
    kmx = os.path.join(a.out, f"angle_k_{a.k}.kmx")

    def say(msg):
        with open(log, "a") as fh:
            fh.write(time.strftime("%H:%M:%SZ ", time.gmtime()) + msg + "\n")

    def stage(sn, gi, dump, extra):
        """annealer rounds on group gi of `dump` until a closed build; extra = build args fixing earlier groups.
        Returns True when closed (kmx written, closed<sn>.json and dumpg<sn>.json saved)."""
        d = json.load(open(dump))
        g = d["groups"][gi]
        say(f"stage {sn}: group {gi}: b {g['b']} a {g['a']} dimW {g['dim_w']} -> {g['made']}")
        if g["b"] != 7 or g["a"] != 14 or g["made"] <= a.P:
            say(f"stage {sn}: group {gi} is not a 7-group with its column kill above the target: skip")
            return False
        sdir = os.path.join(a.out, f"s{sn}")
        os.makedirs(sdir, exist_ok=True)
        if sn == 1:
            for f in glob.glob(os.path.join(a.out, "an_*.json")):     # spaces handed in (seed spaces)
                os.replace(f, os.path.join(sdir, os.path.basename(f)))
        nsp = 0
        for rd in range(a.rounds):
            procs = []
            for i, c in enumerate(cores):
                rs = a.rseed0 + 1000 * sn + rd * len(cores) + i
                procs.append(subprocess.Popen(
                    ["nice", "-n", "10", "taskset", "-c", c, PY, "-I", grp, "anneal", dump, "--g", str(gi), "--P",
                     str(a.P), "--secs", str(a.secs), "--every", "20000", "--price", str(a.P - 2), "--rseed", str(rs),
                     "--save", os.path.join(sdir, f"an_{rs}.json")], env=env,
                    stdout=open(os.path.join(sdir, f"an_{rs}.log"), "w"), stderr=subprocess.STDOUT))
            say(f"stage {sn} round {rd}: {len(procs)} annealers")
            while True:
                time.sleep(20)
                spaces = []
                for f in sorted(glob.glob(os.path.join(sdir, "an_*.json"))):
                    try:
                        spaces += json.load(open(f))
                    except Exception:  # noqa: BLE001
                        pass
                if len(spaces) > nsp:
                    nsp = len(spaces)
                    sp = os.path.join(sdir, "spaces.json")
                    json.dump(spaces, open(sp, "w"))
                    pl = os.path.join(sdir, "plans.json")
                    subprocess.run(["taskset", "-c", cores[0], PY, "-I", grp, "plan", dump, sp, "--g", str(gi),
                                    "--tries", "20000", "--out", pl], env=env, stdout=open(log, "a"),
                                   stderr=subprocess.STDOUT)
                    n = len(json.load(open(pl)))
                    bl = os.path.join(sdir, "build.log")
                    if sn == 1:
                        args = [a.table, pl, "--g", str(gi)]
                    else:
                        args = [a.table, extra["plans1"], "--g", str(extra["g1"]), "--fix1", extra["closed1"],
                                "--inj2", pl, "--g2", str(gi)]
                    gd = d["groups"][2]
                    subprocess.run(["taskset", "-c", cores[0], PY, "-I", grp, "build"] + args +
                                   ["--k", str(a.k), "--plan", a.plan, "--a", str(gd["a"]), "--b", str(gd["b"]),
                                    "--seed", str(a.seed), "--maxplans", "5000", "--tries", str(max(400, min(n, 5000))),
                                    "--saveclosed", os.path.join(a.out, f"closed{sn}.json"),
                                    "--dumpg", os.path.join(a.out, f"dumpg{sn}.json"), "--emit", sdir],
                                   env=env, stdout=open(bl, "a"), stderr=subprocess.STDOUT)
                    txt = open(bl).read().strip().splitlines()
                    say(f"stage {sn}: {len(spaces)} spaces, {n} plans: " + " ".join(txt[-2:])[:300])
                    sk = os.path.join(sdir, f"angle_k_{a.k}.kmx")
                    if os.path.exists(sk):
                        for p in procs:
                            p.kill()
                        os.replace(sk, kmx)
                        with open(os.path.join(a.out, "RESULT.txt"), "a") as fh:
                            fh.write(f"stage {sn}: " + txt[-1] + "\n")
                        say(f"CLOSED stage {sn}: " + txt[-1])
                        return True
                if all(p.poll() is not None for p in procs):
                    break
        say(f"stage {sn}: no closed loader in the rounds given")
        return False

    dump = os.path.join(a.out, f"g_k{a.k}.json")
    if not os.path.exists(dump):
        seeds = []
        for part in (a.seeds if a.seed < 0 else str(a.seed)).split(","):
            lo, _, hi = part.partition("-")
            seeds += list(range(int(lo), int(hi or lo) + 1))
        for sd in seeds:                       # R59: the first engine seed whose dump keeps every planned kill
            subprocess.run(["taskset", "-c", cores[0], PY, "-I", grp, "dump", a.table, "--k", str(a.k), "--plan",
                            a.plan, "--avoid", "0,1", "--pad", "--seed", str(sd), "--out", dump + ".try"], env=env,
                           stdout=open(log, "a"), stderr=subprocess.STDOUT)
            try:
                dd = json.load(open(dump + ".try"))
            except Exception:  # noqa: BLE001
                continue
            say(f"seed {sd}: {dd['products']} products: " + " | ".join(str(g["made"]) for g in dd["groups"]))
            if ",".join(str(g["made"]) for g in dd["groups"]) == a.wantg:
                os.replace(dump + ".try", dump)
                a.seed = sd
                with open(os.path.join(a.out, "SEED.txt"), "w") as fh:
                    fh.write(str(sd))
                break
        if not os.path.exists(dump):
            say("no engine seed keeps the kills: done (no closed loader)")
            return
    elif a.seed < 0:
        a.seed = int(open(os.path.join(a.out, "SEED.txt")).read().strip())
    d = json.load(open(dump))
    say(f"k{a.k} seed {a.seed}: {d['products']} products; groups " +
        " | ".join(f"b{g['b']} a{g['a']} -> {g['made']}" for g in d["groups"]))
    ok1 = stage(1, 2, dump, None)
    if ok1 and not a.nostage2 and len(d["groups"]) > 3:
        stage(2, 3, os.path.join(a.out, "dumpg1.json"),
              dict(plans1=os.path.join(a.out, "s1", "plans.json"), g1=2, closed1=os.path.join(a.out, "closed1.json")))
    say("done" + (" (closed)" if os.path.exists(kmx) else " (no closed loader)"))


if __name__ == "__main__":
    main()
