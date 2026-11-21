#!/usr/bin/env python3
"""r63_inv_ghost_state.py -- R63-INV-GHOST (designer): END's entry state on ALL 3,133,452 inputs and the questions
about ghosting root wires there.  The state builder is R58-INV-38's (tools/r58_inv_38_redund.py state(), read-only
copy of its text): k9f's exact model through R57's model (R54's main_trace): root wires q15..q21 (v), word wires
w1..w25 (w26 is the constant one), MAIN's two leftover wires r0, r1 (R58's reading), flags f0..f3.

  dump   SRC WD HTAB.npz OUT.npz            pack the state (wire i at bit i of an int64), save it with the names
  progs  STATE.npz PROGDIR                  replay FC's rebuild programs sc8_q19 / q20 / q21 on every input
  fn     STATE.npz TARGET POOL              is TARGET a function of the wires POOL (comma list, ranges w20-w25)?
  supp   STATE.npz TARGET POOL [RMAX]       every inclusion-minimal support inside POOL up to RMAX wires (exhaustive
                                            by size, smallest first; stops at the first size with a function)
  greedy STATE.npz TARGET POOL [ORDERS]     inclusion-minimal supports by greedy elimination
  prob   STATE.npz TARGET POOL OUT.json     the minimum-AND problem file (F1's format: rows, inputs tt, target)
  check  STATE.npz PROG.json                replay one program file (F1's format) on every input
Safe to run: reads the model inputs in my folder; writes only the files named on the command line (my folder)."""
import itertools
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def build(hpath):
    import r57_inv_arch2_model as M       # noqa: E402  (read-only; needs SRC WD on argv like R58's tool)
    DC = M.DC
    n, k, d, delta, y = M.chain_y({})
    h, seen = DC.htab(hpath)
    U = y + h[(y >> 9) + 1024]
    Y, Rs, sg, fin = DC.main_trace(U, d)
    R1 = Rs[1]
    neg = R1 < 0
    Mw = R1 >> 8
    A = 1 + ((d >> 7) & 1)
    L = (Mw >> 1) & 3
    res = np.where(neg, L + A, L - A)
    names, cols = [], []
    v = (d >> 1) & 127
    for i in range(7):
        names.append("q%d" % (15 + i))
        cols.append(((v >> i) & 1).astype(np.int64))
    for j in range(1, 26):
        names.append("w%d" % j)
        cols.append(((Y >> j) & 1).astype(np.int64))
    names += ["r0", "r1"]
    cols += [(res & 1).astype(np.int64), (Mw & 1).astype(np.int64)]
    for i in range(4):
        names.append("f%d" % i)
        cols.append(((k >> i) & 1).astype(np.int64))
    code = np.zeros(len(d), dtype=np.int64)
    for i, c in enumerate(cols):
        code |= c << i
    return names, code, dict(Y=Y.astype(np.int64), d=d.astype(np.int64), k=k.astype(np.int64),
                             U=np.asarray(U).astype(np.int64))


def load(path):
    z = np.load(path, allow_pickle=False)
    names = [str(s) for s in z["names"]]
    return names, z["code"]


def parse_pool(names, s):
    out = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-")
            pre = a[0]
            lo, hi = int(a[1:]), int(b[1:])
            out += ["%s%d" % (pre, i) for i in range(lo, hi + 1)]
        else:
            out.append(part)
    for nm in out:
        assert nm in names, nm
    return out


def col(names, code, nm):
    return (code >> names.index(nm)) & 1


def mask_of(names, pool):
    m = 0
    for nm in pool:
        m |= 1 << names.index(nm)
    return np.int64(m)


def classes(names, code, pool, tname):
    """distinct patterns of the pool with the target's value; None if the target is not a function of the pool."""
    key = code & mask_of(names, pool)
    t = col(names, code, tname)
    u = np.unique(key * 2 + t)
    kk = u >> 1
    if np.any(kk[1:] == kk[:-1]):
        return None
    return kk, (u & 1)


def is_fn_sampled(samp_code, samp_t, mask):
    u = np.unique((samp_code & mask) * 2 + samp_t)
    kk = u >> 1
    return not np.any(kk[1:] == kk[:-1])


def supports(names, code, tname, pool, rmax, nsamp=60000):
    rng = np.random.default_rng(11)
    pick = rng.choice(len(code), size=min(nsamp, len(code)), replace=False)
    sc = code[pick]
    st = col(names, code, tname)[pick]
    full_t = col(names, code, tname)
    found = []
    tested = 0
    for r in range(0, rmax + 1):
        for sup in itertools.combinations(pool, r):
            tested += 1
            m = mask_of(names, sup)
            if not is_fn_sampled(sc, st, m):
                continue
            u = np.unique((code & m) * 2 + full_t)
            kk = u >> 1
            if np.any(kk[1:] == kk[:-1]):
                continue
            found.append((list(sup), int(len(u))))
        if found:
            break
    return found, tested


def greedy(names, code, tname, pool, norders=6):
    full_t = col(names, code, tname)

    def ok(S):
        u = np.unique((code & mask_of(names, S)) * 2 + full_t)
        kk = u >> 1
        return not np.any(kk[1:] == kk[:-1]), len(u)

    good, npat = ok(pool)
    if not good:
        return None
    rng = np.random.default_rng(7)
    orders = [list(pool), list(reversed(pool))]
    for _ in range(max(0, norders - 2)):
        o = list(pool)
        rng.shuffle(o)
        orders.append(o)
    res = []
    for order in orders:
        S = list(pool)
        for w in order:
            S2 = [i for i in S if i != w]
            if ok(S2)[0]:
                S = S2
        S = [nm for nm in pool if nm in S]
        if S not in [r_[0] for r_ in res]:
            res.append((S, ok(S)[1]))
    return res


def problem(names, code, tname, pool, outp):
    cl = classes(names, code, pool, tname)
    assert cl is not None, "%s is not a function of %s" % (tname, pool)
    kk, tv = cl
    N = len(kk)
    ins = []
    for nm in pool:
        b = (kk >> names.index(nm)) & 1
        tt = 0
        for r in np.flatnonzero(b):
            tt |= 1 << int(r)
        ins.append(dict(name=nm, tt="%x" % tt))
    tg = 0
    for r in np.flatnonzero(tv):
        tg |= 1 << int(r)
    json.dump(dict(rows=N, inputs=ins, target="%x" % tg, target_name=tname), open(outp, "w"))
    return N


def replay(names, code, prog, tname=None):
    """F1's program format: gates = AND(XOR left, XOR right) over inputs then earlier gates; -1 = constant one."""
    if "stages" in prog:                               # my staged format: the target is the XOR of all stages
        nin = len(prog["inputs"])
        tot = np.zeros(len(code), dtype=np.int64)
        for st in prog["stages"]:
            sub = dict(inputs=prog["inputs"], gates=st["gates"],
                       outputs=[list(st["out"]) + [nin + g for g in range(int(st["internal"]), len(st["gates"]))]])
            tot ^= replay(names, code, sub)[0]
        return [tot]
    vals = [col(names, code, nm) for nm in prog["inputs"]]
    for g in prog["gates"]:
        a = np.zeros(len(code), dtype=np.int64)
        b = np.zeros(len(code), dtype=np.int64)
        for i in g["left"]:
            a ^= 1 if i == -1 else vals[i]
        for i in g["right"]:
            b ^= 1 if i == -1 else vals[i]
        vals.append(a & b)
    outs = []
    for out in prog["outputs"]:
        o = np.zeros(len(code), dtype=np.int64)
        for i in out:
            o ^= 1 if i == -1 else vals[i]
        outs.append(o)
    return outs


def main():
    cmd = sys.argv[1]
    t0 = time.time()
    if cmd == "dump":
        hpath, outp = os.path.expanduser(sys.argv[4]), os.path.expanduser(sys.argv[5])
        sys.argv = [sys.argv[0]] + sys.argv[2:]          # R57's model reads SRC WD from argv[1:3]
        names, code, extra = build(hpath)
        np.savez_compressed(outp, names=np.array(names), code=code, **extra)
        print("state: %d inputs, %d wires, %d distinct codes, %.0f s -> %s" % (
            len(code), len(names), len(np.unique(code)), time.time() - t0, outp), flush=True)
        return
    names, code = load(os.path.expanduser(sys.argv[2]))
    if cmd == "progs":
        pdir = os.path.expanduser(sys.argv[3])
        for q in sys.argv[4:] or ["q19", "q20", "q21"]:
            prog = json.load(open(os.path.join(pdir, "sc8_%s.prog.json" % q)))
            (o,) = replay(names, code, prog)
            bad = int(np.count_nonzero(o != col(names, code, q)))
            print("sc8_%s: inputs %s, %d ANDs: %d of %d inputs wrong" % (q, prog["inputs"], len(prog["gates"]), bad,
                                                                         len(code)), flush=True)
        return
    if cmd == "check":
        prog = json.load(open(os.path.expanduser(sys.argv[3])))
        outs = replay(names, code, prog)
        tn = prog.get("targets") or [prog.get("target_name") or sys.argv[4]]
        nand = sum(len(st["gates"]) for st in prog["stages"]) if "stages" in prog else len(prog["gates"])
        for o, t in zip(outs, tn):
            bad = int(np.count_nonzero(o != col(names, code, t)))
            print("%s: inputs %s, %d ANDs: %d of %d inputs wrong" % (t, prog["inputs"], nand, bad, len(code)),
                  flush=True)
        return
    tname = sys.argv[3]
    for a in [a_ for a_ in sys.argv if a_.startswith("where=")]:      # where=w25:0,w24:1  restricts the inputs
        sys.argv.remove(a)
        for term in a[6:].split(","):
            nm, val = term.split(":")
            code = code[col(names, code, nm) == int(val)]
        print("(restricted to %s: %d inputs)" % (a[6:], len(code)), flush=True)
    if "&" in tname:                                   # a derived target: q18&!w25 = q18 AND NOT w25 (a range piece)
        a_, b_ = tname.split("&")
        neg = b_.startswith("!")
        piece = col(names, code, a_) & (col(names, code, b_.lstrip("!")) ^ (1 if neg else 0))
        code = code | (piece << 62)
        names = names + [None] * (62 - len(names)) + [tname]
    pool = [p for p in parse_pool([n_ for n_ in names if n_], sys.argv[4]) if p != tname]
    # only the distinct joint patterns of pool + target matter for every question below (exact reduction)
    code = np.unique(code & mask_of(names, pool + [tname]))
    print("(%d distinct joint patterns of the %d pool wires + %s)" % (len(code), len(pool), tname), flush=True)
    if cmd == "fn":
        cl = classes(names, code, pool, tname)
        print("%s a function of %d wires %s: %s%s (%.0f s)" % (
            tname, len(pool), ",".join(pool), "NO" if cl is None else "YES",
            "" if cl is None else ", %d care classes" % len(cl[0]), time.time() - t0), flush=True)
    elif cmd == "supp":
        rmax = int(sys.argv[5]) if len(sys.argv) > 5 else 8
        found, tested = supports(names, code, tname, pool, rmax)
        print("%s in pool of %d: %d supports tested up to %d wires, %d functions (%.0f s)" % (
            tname, len(pool), tested, rmax, len(found), time.time() - t0), flush=True)
        for sup, npat in found:
            print("   %d wires, %d classes: %s" % (len(sup), npat, ",".join(sup)), flush=True)
    elif cmd == "greedy":
        res = greedy(names, code, tname, pool, int(sys.argv[5]) if len(sys.argv) > 5 else 6)
        if res is None:
            print("%s is NOT a function of the pool %s" % (tname, ",".join(pool)), flush=True)
        else:
            for S, npat in sorted(res, key=lambda r_: len(r_[0])):
                print("%s: support %d wires, %d classes: %s" % (tname, len(S), npat, ",".join(S)), flush=True)
        print("(%.0f s)" % (time.time() - t0), flush=True)
    elif cmd == "prob":
        N = problem(names, code, tname, pool, os.path.expanduser(sys.argv[5]))
        print("%s from %s: %d care classes -> %s" % (tname, ",".join(pool), N, sys.argv[5]), flush=True)


if __name__ == "__main__":
    main()
