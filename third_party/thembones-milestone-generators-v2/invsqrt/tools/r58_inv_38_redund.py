#!/usr/bin/env python3
"""r58_inv_38_redund.py -- R58-INV-38 (designer): WHICH WIRE OF END'S ENTRY STATE IS A FUNCTION OF THE OTHERS, and
of how few?  The state after MAIN on k9f's exact model (R54's main_trace through R57's model, read-only), ALL
3,133,452 inputs: root wires q15..q21 (v), word wires w1..w25 (w26 is the constant), MAIN's two leftover wires
(r0 = the low bit of the last row's 2-bit result, r1 = the skipped lowest target of the last row), flags f0..f3 (k).
For every target wire: is it a function of ALL the others (else it is information: no retire at any price); if so an
inclusion-minimal support found by greedy elimination in several orders, and the ANF size of one completion when the
support has at most 18 wires.  Scenario "after": q19, q20, q21 are not available (END has retired them).
  r58_inv_38_redund.py SRC WD HTAB.npz OUT.jsonl NPROC [orders] [targets, comma separated] [small=RMAX]
      small=RMAX: for each target (END's state after its three retires) EVERY support of at most RMAX wires
Env R58_BAN=q18,q19,q20,q21: the wires not available as inputs in scenario "after" (the target is always excluded).
Env R58_SLOT=par: the same scan on the parity slot's chain at m 7 (ideal parity table), where c0 does not exist.
Safe to run: reads only; writes OUT.jsonl in my folder."""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r57_inv_arch2_model as M       # noqa: E402

DC = M.DC
G = {}


def state(hpath):
    if "cols" in G:
        return G
    if os.environ.get("R58_SLOT"):                      # the parity slot's chain at m 7 (R56's cheap chain), ideal
        par = int(os.environ["R58_SLOT"])               # parity table of parity T_0 xor par: U = 2 y19 + hp(T)
        n, k, d, delta, y = M.chain_y(dict(m=7, F1=4, sched1=[0, 0, 0, 0, 1, 2, 4, 5], sched2=[0, 0, 0, 0, 1, 3, 4, 5]))
        T = y >> 8
        yt = (T.astype(np.float64) + 0.5) * 512.0
        gi = M.RM.f_ideal(yt) - yt + 0.5
        t = (T & 1) ^ par
        hp = (2.0 * np.round((gi - t) / 2.0) + t).astype(np.int64)
        U = 2 * y + hp
        print("slot chain (m 7, F1 4, ideal parity table, par %d): U %d..%d, distinct y %d, distinct U %d" % (
            par, int(U.min()), int(U.max()), len(np.unique(y)), len(np.unique(U))), flush=True)
    else:
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
    G.update(names=names, cols=cols, n=len(d))
    return G


def ban_set(scen):
    """wires not available as inputs: END's three retires in scenario "after"; env R58_BAN = another list."""
    if scen != "after":
        return set()
    if os.environ.get("R58_BAN"):
        return set(os.environ["R58_BAN"].split(","))
    return {"q19", "q20", "q21"}


def code_of(idx):
    c = np.zeros(G["n"], dtype=np.int64)
    for t, i in enumerate(idx):
        c |= G["cols"][i] << t
    return c


def full_code():
    if "full" not in G:
        c = np.zeros(G["n"], dtype=np.int64)
        for i, col in enumerate(G["cols"]):
            c |= col << i                              # wire i at bit i (38 wires: fits an int64)
        G["full"] = c
    return G["full"]


def is_fn(idx, tcol):
    """tcol a function of the wires idx?  one sort of the masked state code, then neighbours compared."""
    mask = 0
    for i in idx:
        mask |= 1 << i
    c = full_code() & np.int64(mask)
    order = np.argsort(c, kind="stable")
    cs, ts = c[order], tcol[order]
    same = cs[1:] == cs[:-1]
    bad = bool(np.any(same & (ts[1:] != ts[:-1])))
    return (not bad), int(len(cs) - int(same.sum()))


def anf(idx, tcol):
    r = len(idx)
    c = code_of(idx)
    f = np.zeros(1 << r, dtype=np.uint8)
    f[c[tcol == 1]] = 1
    for i in range(r):
        step = 1 << i
        fv = f.reshape(-1, 2 * step)
        fv[:, step:] ^= fv[:, :step]
    idxs = np.flatnonzero(f)
    pc = np.array([bin(int(i)).count("1") for i in idxs]) if len(idxs) else np.array([0])
    return int(pc.max()), int((pc >= 2).sum()), int(len(idxs))


def work(args):
    tname, scen, hpath, norders = args
    g = state(hpath)
    names = g["names"]
    ti = names.index(tname)
    tcol = g["cols"][ti]
    banned = ban_set(scen)
    others = [i for i, nm in enumerate(names) if i != ti and nm not in banned]
    t0 = time.time()
    ok, npat = is_fn(others, tcol)
    out = dict(target=tname, scenario=scen, function_of_all_others=bool(ok), patterns=int(npat))
    if not ok:
        out["t"] = round(time.time() - t0, 1)
        return out
    best = None
    rng = np.random.default_rng(7)
    orders = [list(others), list(reversed(others))]
    for _ in range(max(0, norders - 2)):
        o = list(others)
        rng.shuffle(o)
        orders.append(o)
    for order in orders:
        S = list(others)
        for w in order:
            S2 = [i for i in S if i != w]
            if is_fn(S2, tcol)[0]:
                S = S2
        if best is None or len(S) < len(best):
            best = S
    out["support"] = [names[i] for i in best]
    out["support_size"] = len(best)
    if len(best) <= 18:
        deg, ands, mons = anf(best, tcol)
        out.update(anf_degree=deg, anf_monomials_deg2plus=ands, anf_monomials=mons)
    out["t"] = round(time.time() - t0, 1)
    return out


def small(args):
    """every support of at most `rmax` wires for one target: refuted on a random sample of inputs (a conflict in the
    sample is a proof that the target is not a function of that support), survivors re-tested on all inputs."""
    tname, scen, hpath, rmax, nsamp = args
    import itertools
    g = state(hpath)
    names = g["names"]
    ti = names.index(tname)
    banned = ban_set(scen)
    others = [i for i, nm in enumerate(names) if i != ti and nm not in banned]
    rng = np.random.default_rng(11)
    pick = rng.choice(g["n"], size=nsamp, replace=False)
    fc = full_code()[pick]
    tc = g["cols"][ti][pick]
    t0 = time.time()
    found, tested = [], 0
    for r in range(0, rmax + 1):
        for sup in itertools.combinations(others, r):
            tested += 1
            mask = 0
            for i in sup:
                mask |= 1 << i
            c = (fc & np.int64(mask)) * 2 + tc
            u = np.unique(c)
            if np.any((u[1:] >> 1) == (u[:-1] >> 1)):
                continue                                 # both values of the target under one pattern: refuted
            ok_, npat = is_fn(list(sup), g["cols"][ti])
            if ok_:
                deg, ands, mons = anf(list(sup), g["cols"][ti])
                found.append(dict(support=[names[i] for i in sup], reachable_patterns=npat, anf_degree=deg,
                                  anf_monomials_deg2plus=ands, anf_monomials=mons))
        if found:
            break
    return dict(target=tname, scenario=scen, mode="every support up to %d wires" % rmax, supports_tested=tested,
                sample=nsamp, functions_found=found, t=round(time.time() - t0, 1))


def main():
    hpath, outp, nproc = os.path.expanduser(sys.argv[3]), os.path.expanduser(sys.argv[4]), int(sys.argv[5])
    norders = int(sys.argv[6]) if len(sys.argv) > 6 else 4
    g = state(hpath)
    names = g["names"]
    print("END entry state: %d wires, %d inputs, distinct states %d" % (len(names), g["n"], len(np.unique(code_of(range(len(names)))))), flush=True)
    jobs = []
    full_code()
    tlist = sys.argv[7].split(",") if len(sys.argv) > 7 else (
        ["r0", "r1"] + ["q%d" % i for i in range(15, 22)] + ["w%d" % j for j in range(25, 0, -1)] + ["f%d" % i for i in range(4)])
    for scen in ("after", "before"):
        for nm in tlist:
            if scen == "after" and nm in ("q19", "q20", "q21") and not os.environ.get("R58_BAN"):
                continue
            jobs.append((nm, scen, hpath, norders))
    import multiprocessing as mp
    if len(sys.argv) > 8 and sys.argv[8].startswith("small="):       # small=RMAX: the exhaustive small-support search
        rmax = int(sys.argv[8].split("=")[1])
        jobs = [(nm, scen, hpath, rmax, 3000) for (nm, scen, _, _) in jobs if scen == "after"]
        with mp.get_context("fork").Pool(nproc) as pool, open(outp, "a") as f:
            for o in pool.imap(small, jobs):
                f.write(json.dumps(o) + chr(10))
                f.flush()
                print(json.dumps(o), flush=True)
        print("DONE", flush=True)
        return
    with mp.get_context("fork").Pool(nproc) as pool, open(outp, "a") as f:
        for o in pool.imap(work, jobs):
            f.write(json.dumps(o) + "\n")
            f.flush()
            print(json.dumps(o), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
