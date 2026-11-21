"""invsqrt_gen: the routine hand jobs as single commands (own tool, workstream invsqrt_gen). Runbook with the
checks and the traps: top of research/ws_invsqrt_gen.md, section "How to pin a module or take a table set".

Run inside WSL through tools/sh/invsqrt_gen_run.sh (written here as `hand`):
  hand stop                      ask the unattended loop to stop, wait until it and every build of ours has
                                 ended, print the state (tools/sh/invsqrt_gen_state.sh)
  hand start                     start the loop again (it takes its line from the top of best/)
  hand line <name> [k=v ...]     write hand/<name>.line = the entry line of the build at the top of best/ with
                                 these keys changed (k= removes a key); prints it
  hand pin <n> <module.py> <folder with invsqrt_ph_net.py, _stream.py, _sx.py> [k=v ...]
                                 new pinned copy of a phase module: tools/invsqrt_gen_pin/ph<n>/ (n >= 6) with
                                 PIN.json (sources, sha256 of the four files, constructor options k=v). The
                                 builder option ph=<n> imports it and passes the options to PhaseDelivery.
  hand rpin <k> <rows2 module.py>
                                 new pinned copy of the second rows module: tools/invsqrt_gen_rows2p<k>.py
                                 (k >= 2), builder option rpin=<k>
  hand set <source folder> <set> <name>
                                 take the table set <source folder>/<set>.tables.json for the line
                                 hand/<name>.line (which must say tabs=<set> and the set's tie_mask): import and
                                 check on all inputs, write the specs of the line into lut/, take the set's
                                 programs, install the cheapest fitting program per block, price on a sample
  hand try <name> [N] [:: k=v ... :: k=v ...]
                                 price the line (and option sets laid over it) with the exact simulator on N
                                 random + 4,096 edge inputs (default 20000)
  hand entry <tag> <name> "<what the build adds>"
                                 the entry script (two builds, exact check on every input, official verifier
                                 with both salts on the verification box, publication to best/<tag>/); writes
                                 the sample price to the log before the run, and after a publication that is
                                 the new top of best/ one line to the log and the facts of the state block
`set`, `try` and `entry` work in lut/ and need the loop stopped (`hand stop` first; they refuse otherwise).
Nothing is published except by the entry script, so a mistake here cannot put a wrong file into best/.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time

T = "<HOME>/Quantum/tools"
sys.path.insert(0, T)
import invsqrt_gen_auto as AU                       # run, with_args, local_args, price, best, quiet_time, say
import invsqrt_gen_model as M

W, PY, LOG = AU.W, AU.PY, AU.LOG
H = os.path.join(W, "hand")
BLOCKS = ("LUT1", "LUT1b", "LUT3", "LUT2r", "ADDR")
SUPPORT = ("invsqrt_ph_net.py", "invsqrt_ph_stream.py", "invsqrt_ph_sx.py")


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def script(name):
    """CR-stripped copy of tools/sh/invsqrt_gen_<name>.sh in hand/ (the files in tools/sh have CRLF ends)."""
    os.makedirs(H, exist_ok=True)
    dst = os.path.join(H, f"{name}.sh")
    open(dst, "w", newline="\n").write(open(f"{T}/sh/invsqrt_gen_{name}.sh").read().replace("\r", ""))
    return dst


def procs(pattern):
    r = subprocess.run(["pgrep", "-af", pattern], capture_output=True, text=True)
    return [l for l in r.stdout.splitlines() if "invsqrt_gen_hand" not in l and "pgrep" not in l]


def loop_running():
    return bool(procs("invsqrt_gen_auto[.]py"))


def builds_running():
    # only OUR builds (other workstreams use box_verify.sh too; their runs are not our business here)
    return procs("invsqrt_gen_build|invsqrt_gen/build/box_verify_run|invsqrt_gen_sweep")


def need_quiet():
    if loop_running() or not os.path.exists(os.path.join(AU.A, "STOP")):
        sys.exit("the unattended loop is running: `hand stop` first (it rewrites lut/ every few minutes)")
    if builds_running():
        sys.exit("a build, sweep or verification of ours is running: wait for it\n" + "\n".join(builds_running()))


def read_line(name):
    p = name if os.path.exists(name) else os.path.join(H, f"{name}.line")
    return open(p).read().split()


def blocks_now():
    out = []
    for b in BLOCKS:
        try:
            spec = json.load(open(os.path.join(W, "lut", f"{b}.spec.json")))["spec_sha"]
            pr = json.load(open(os.path.join(W, "lut", f"{b}.prog.json")))
            n = sum(1 for s in pr["steps"] if s[0] == "and")
            fit = "fits" if pr.get("spec_sha") == spec else "DOES NOT FIT the spec"
            out.append(f"  {b}: spec {spec}  program {n} ANDs, {fit}  ({pr.get('origin', 'table set / lookup workstream')})")
        except Exception as ex:
            out.append(f"  {b}: no spec or no program in lut/ ({type(ex).__name__}): our own solver would run "
                       f"at the line's depth")
    return "\n".join(out)


def log(text):
    with open(LOG, "a", encoding="utf-8", newline="\n") as fh:
        fh.write(text.rstrip("\n") + "\n")


def cmd_stop():
    os.makedirs(AU.A, exist_ok=True)
    open(os.path.join(AU.A, "STOP"), "w").close()
    print("stop requested; waiting for the loop and for builds of ours to end (an entry run in progress finishes)")
    t0 = time.time()
    while (loop_running() or builds_running()) and time.time() - t0 < 2400:
        time.sleep(5)
    if loop_running() or builds_running():
        sys.exit("still running after 40 minutes:\n" + "\n".join(procs("invsqrt_gen_auto[.]py") + builds_running()))
    print(f"the loop has ended and nothing of ours is building ({time.time() - t0:.0f} s)")
    print(AU.run(f"bash {script('state')} 3"))


def cmd_start():
    if builds_running():
        sys.exit("a build of ours is still running: wait for it, then start the loop")
    print(AU.run(f"bash {script('auto')} start").strip())


def cmd_line(name, *changes):
    args = AU.with_args(AU.best()[2], **{k: (v or None) for k, v in (c.split("=", 1) for c in changes)})
    os.makedirs(H, exist_ok=True)
    open(os.path.join(H, f"{name}.line"), "w", newline="\n").write(" ".join(args) + "\n")
    print(f"hand/{name}.line (from the top of best/ = {AU.best()[1]}):\n" + " ".join(args))


def cmd_pin(n, module, support, *opts):
    n = int(n)
    assert n >= 6, "ph1 .. ph5 are fixed in the builder; new pins start at 6"
    dst = os.path.join(T, "invsqrt_gen_pin", f"ph{n}")
    if os.path.exists(dst):
        sys.exit(f"{dst} exists: a pinned copy is never changed, take the next number")
    for f in SUPPORT:
        assert os.path.exists(os.path.join(support, f)), f"{support} has no {f}"
    os.makedirs(dst)
    shutil.copy(module, os.path.join(dst, "invsqrt_ph_mod.py"))
    for f in SUPPORT:
        shutil.copy(os.path.join(support, f), os.path.join(dst, f))
    options = {k: (int(v) if v.lstrip("-").isdigit() else v) for k, v in (o.split("=", 1) for o in opts)}
    rec = dict(pinned=time.strftime("%Y-%m-%d %H:%MZ", time.gmtime()), module=os.path.abspath(module),
               support=os.path.abspath(support), options=options,
               sha256={f: sha(os.path.join(dst, f)) for f in ("invsqrt_ph_mod.py",) + SUPPORT})
    json.dump(rec, open(os.path.join(dst, "PIN.json"), "w", newline="\n"), indent=1)
    print(f"pinned as ph={n}: {dst}")
    for f, h in rec["sha256"].items():
        print(f"  {f}  {h[:16]}")
    print(f"  options passed to PhaseDelivery: {options}")
    print(f"next: hand line <name> ph={n}; hand try <name> :: ph=<old number> :: ph={n}")


def cmd_rpin(k, module):
    k = int(k)
    assert k >= 2, "rows2p1 is the copy in use; new pins start at 2"
    dst = os.path.join(T, f"invsqrt_gen_rows2p{k}.py")
    if os.path.exists(dst):
        sys.exit(f"{dst} exists: a pinned copy is never changed, take the next number")
    shutil.copy(module, dst)
    rec = dict(pinned=time.strftime("%Y-%m-%d %H:%MZ", time.gmtime()), module=os.path.abspath(module), sha256=sha(dst))
    json.dump(rec, open(os.path.join(T, "invsqrt_gen_pin", f"rows2p{k}.json"), "w", newline="\n"), indent=1)
    print(f"pinned as rpin={k}: {dst}  {rec['sha256'][:16]}")
    print(f"next: hand line <name> rpin={k}; hand try <name> :: rpin=<old number> :: rpin={k}")


def cmd_set(src, name, line):
    need_quiet()
    args = read_line(line)
    d = dict(a.split("=", 1) for a in args)
    assert d.get("tabs") == name, f"the line must say tabs={name} (hand line <name> tabs={name} tie_mask=...)"
    params = [a for a in args if a.split("=")[0] in M.Params.__dataclass_fields__]
    tfile = os.path.join(os.path.expanduser(src), f"{name}.tables.json")
    out = AU.run(f"cd {W} && {PY} {T}/invsqrt_gen_tabs.py import {tfile} {name} {' '.join(params)}")
    print("\n".join(l[:400] for l in out.strip().splitlines()[-3:]))
    bad = [k for k in ("all inputs: PASS", "'tie_t0': True", "'tie_c1': True", "'c0_even': True", "'window': True",
                       "'u_in_range': True") if k not in out]
    if d.get("zone") == "1" and "'zone_copies': True" not in out:
        bad.append("'zone_copies': True")
    if bad:
        sys.exit(f"IMPORT CHECK FAILED ({bad}): the set is not used; send this output back to its owner")
    print("== specs of the line into lut/ (build at depth 4, no check)")
    out = AU.run(f"cd {W} && nice -n 5 {PY} {T}/invsqrt_gen_build.py g6 {W}/build/hand_tmp "
                 f"{' '.join(AU.with_args(args, deep=4))} --no-check")
    if "kmb sha256" not in out:
        sys.exit("the line does not build with these tables:\n" + out[-1500:])
    print("== programs of the set")
    print(AU.run(f"bash {script('takeset')} {os.path.expanduser(src)} {name}").strip())
    print("== install (cheapest fitting program per block: the set's, the lookup workstream's, our own cache)")
    AU.run(f"cd {W} && {PY} {T}/invsqrt_gen_lutcmp.py install quiet")
    print(blocks_now())
    cmd_try(line)


def cmd_try(line, *rest):
    need_quiet()
    args = AU.local_args(read_line(line))            # writes the specs, installs, keeps the depth if replayable
    n = "20000"
    if rest and rest[0].isdigit():
        n, rest = rest[0], rest[1:]
    out = AU.run(f"cd {W} && nice -n 5 {PY} {T}/invsqrt_gen_sweep.py g6 {' '.join(args)} sample={n} {' '.join(rest)}")
    print("\n".join(l[:440] for l in out.strip().splitlines() if "Q " in l or "Error" in l or "assert" in l.lower()))
    print(f"(best published: {AU.best()[1]} {AU.best()[0]:,.2f}; line built at deep="
          f"{dict(a.split('=', 1) for a in args).get('deep')})")
    return args


def price(args, n=20000):
    """[(set, ok, score, Tof, Q)] of the line on a sample of n random + 4,096 edge inputs (as AU.price, which
    uses 3,000: about 15 points of noise; 20,000 gives about 5)."""
    out = []
    for l in AU.run(f"cd {W} && nice -n 5 {PY} {T}/invsqrt_gen_sweep.py g6 {' '.join(args)} sample={n}").splitlines():
        m = re.match(r"(.*?): Q (\d+) .*\| sample (ok|FAIL)\s+Tof ([0-9.]+).*?score ([0-9,]+)", l)
        if m:
            out.append((m.group(1).strip(), m.group(3) == "ok", float(m.group(5).replace(",", "")),
                        float(m.group(4)), int(m.group(2))))
    return out


def cmd_entry(tag, line, what):
    need_quiet()
    if os.path.exists(os.path.join(W, "best", tag)):
        sys.exit(f"best/{tag} exists: take a new tag")
    if AU.quiet_time():
        sys.exit("23:12Z to 23:40Z: the coordinator takes the top of best/ at 23:30Z; no entry run now")
    args = AU.local_args(read_line(line))
    res = price(args)
    bscore, btag, _ = AU.best()
    if not res or not res[0][1]:
        sys.exit(f"the line does not pass the sample check: {res}")
    _, _, score, tof, q = res[0]
    print(f"sample price {score:,.0f} ({tof:.2f} Toffolis, Q {q}); best published {btag} {bscore:,.2f}")
    if score >= bscore and os.environ.get("FORCE") != "1":
        sys.exit("not below the best published build: no entry run (FORCE=1 to run it anyway)")
    now = time.strftime("%H:%MZ", time.gmtime())
    log(f"- HAND ({now}, `tools/invsqrt_gen_hand.py entry`): entry run {tag} started ({what}); price on the "
        f"sample before the run [M, exact simulator on 20,000 + 4,096 inputs]: {score:,.0f}, {tof:.2f} Toffolis, "
        f"Q {q}; best published {btag} {bscore:,.2f}. Line: `{' '.join(args)}`")
    out = AU.run(f"bash {script('build')} {tag} g6 {' '.join(args)}", timeout=3 * 3600)
    open(os.path.join(W, "logs", f"hand_{tag}.log"), "w").write(out)
    tail = [l.strip() for l in out.splitlines() if any(k in l for k in (
        "PUBLISHED", "scores:", "top of best", "not published", "FAILED", "NOT DETERMINISTIC", "sha256 build"))]
    print("\n".join(tail))
    now = time.strftime("%H:%MZ", time.gmtime())
    if any("top of best/ updated" in t for t in tail):
        st = open(os.path.join(W, "best", "STATE.txt")).read().strip()
        log(f"- HAND ({now}): **published and now at the top of `best/` [M, entry script: exact check on all x, "
            f"two bit-identical builds, official verifier with both salts]: {st}**. Entry line: "
            f"`invsqrt_gen_build.sh {tag} g6 {' '.join(args)}`")
        counts = []
        for b in BLOCKS:
            pr = json.load(open(os.path.join(W, "best", tag, "lut", f"{b}.prog.json")))
            counts.append(f"{b} {sum(1 for s in pr['steps'] if s[0] == 'and')}")
        text = f"{what}; lookup programs " + ", ".join(counts)
        print("state block: " + AU.run(f"{PY} {T}/invsqrt_gen_newbest.py '{text}'").strip()[-200:])
    else:
        log(f"- HAND ({now}): entry run {tag} did NOT change the top of `best/`: {' | '.join(tail)[:400]} "
            f"(full output: `logs/hand_{tag}.log`)")
        sys.exit(1)


def main():
    cmds = dict(stop=cmd_stop, start=cmd_start, line=cmd_line, pin=cmd_pin, rpin=cmd_rpin, set=cmd_set,
                entry=cmd_entry)
    if len(sys.argv) < 2 or sys.argv[1] not in tuple(cmds) + ("try",):
        sys.exit(__doc__)
    if sys.argv[1] == "try":
        rest = " ".join(sys.argv[3:]).split()
        cmd_try(sys.argv[2], *rest)
    else:
        cmds[sys.argv[1]](*sys.argv[2:])


if __name__ == "__main__":
    main()
