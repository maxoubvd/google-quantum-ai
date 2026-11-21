#!/usr/bin/env python3
"""r54_inv_tail_alt.py -- R54-INV-TAIL (designer): THE PROGRAM AND THE ANGLES CHOSEN TOGETHER for the last j rows of
a divider tail at K ANDs (fix rounds).  Round: (1) the don't-care table of the current angles (tools/r54_inv_tail_dc.py:
minslack per class and word, all 3,133,452 inputs); (2) relational SAT (tools/r54_inv_tail_sat.py, as subprocesses
under coreutils `timeout`) on a grid of margins M0, in parallel: the highest M0 with a K-AND program; (3) that program
on all inputs with the LP RE-SOLVED: its true margin and new angles; next round from those angles.  Stops when no
higher M0 is SAT.  Every line is appended to OUTDIR/TAG.jsonl; programs OUTDIR/TAG_r<round>.prog.json.

  r54_inv_tail_alt.py SRC WD HTAB.npz OUTDIR TAG DIV1|DIV2|MAIN j K [nproc=8] [tlim=600] [rounds=8] [extra=0] [keep0=0]
                      [solver=cadical153] [multi=1] [free=1] [dv0=free_dv.npz]
                      [base=DIV2:prog.json,MAIN:prog.json] [fam=1 [digits=1] [fixops=1] [taper=r]]
Safe to run: reads only; writes under OUTDIR (my own folder)."""
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r54_inv_tail_dc as DC          # noqa: E402  (imports the R37 model, which reads sys.argv[1] = SRC)


def apply_prog(prog, H, vtop, nv, nH, j):
    ins = [((vtop >> i) & 1) for i in range(nv)] + [((H >> i) & 1) for i in range(nH)]
    vals = list(ins)
    for g in prog["gates"]:
        a = np.zeros(len(H), dtype=np.int64)
        for i in g["left"]:
            a ^= 1 if i == -1 else vals[i]
        b = np.zeros(len(H), dtype=np.int64)
        for i in g["right"]:
            b ^= 1 if i == -1 else vals[i]
        vals.append(a & b)
    wnew = np.zeros(len(H), dtype=np.int64)
    for lst in prog["outputs"]:
        o = np.zeros(len(H), dtype=np.int64)
        for i in lst:
            o ^= 1 if i == -1 else vals[i]
        wnew = (wnew << 1) | o
    return wnew


def wval(w, j):
    v = np.zeros(len(w), dtype=np.int64)
    for i in range(j + 1):
        v += (1 - 2 * ((w >> (j - i)) & 1)) << (j - i)
    return v


def main():
    wd, hpath, outdir, tag = (os.path.expanduser(sys.argv[i]) for i in (2, 3, 4, 5))
    which, j, K = sys.argv[6], int(sys.argv[7]), int(sys.argv[8])
    opt = dict(a.split("=", 1) for a in sys.argv[9:] if "=" in a)
    nproc, tlim, rounds = int(opt.get("nproc", 8)), int(opt.get("tlim", 600)), int(opt.get("rounds", 8))
    extra, keep0, solver = int(opt.get("extra", 0)), int(opt.get("keep0", 0)), opt.get("solver", "cadical153")
    logf = open(os.path.join(outdir, tag + ".jsonl"), "a")

    def log(**kw):
        kw["t"] = round(time.time() - t0, 1)
        logf.write(json.dumps(kw) + "\n")
        logf.flush()
        print(json.dumps(kw), flush=True)
    t0 = time.time()
    progs = {}
    for item in [x for x in opt.get("base", "").split(",") if x]:        # base=DIV2:file.prog.json,MAIN:file
        w_, f_ = item.split(":", 1)
        progs[w_] = json.load(open(os.path.expanduser(f_)))
    assert which not in progs, "the target tail is opened from the circuit's own rows"
    DC.PROGS.update(progs)
    B = DC.base(wd, hpath, progs)
    key = DC.KEY[which]
    err0, fire0, bad0 = DC.downstream("Y", B["Y"], B["d"], B["k"], hpath, B["tgt"])
    assert not bad0.any() and not B.get("T_outside_table")
    mg0, conv, bound, dv = DC.lp_dv(err0, B["tol"], fire0)
    H, vtop, nv, wt, shift, kkmin = DC.classes(B, which, j, extra)
    hmin, hmax = int(H.min()), int(H.max())
    nH = max(int(hmax).bit_length(), int(-hmin - 1).bit_length()) + 1
    code = ((H - hmin) << nv) | vtop
    ucode, cid = np.unique(code, return_inverse=True)
    ncls = len(ucode)
    wmin = np.full(ncls, 1 << 30)
    np.minimum.at(wmin, cid, wt)
    cnt = np.bincount(cid, minlength=ncls)
    vt = DC.wvals(which, wt, j)
    log(ev="base", which=which, j=j, K=K, margin=mg0, conv=conv, classes=ncls, inputs=nv + nH, extra=extra,
        y_differs_from_circuit=B.get("y_differs_from_circuit"), base_progs=sorted(progs))
    DC._B[0], DC._B[1], DC._B[2], DC._B[3] = B, wt, cid, ncls
    import multiprocessing as mp
    cur = None                                   # the current program's margin (None: no program yet)
    mstart = None
    if "dv0" in opt:                             # start from the angles of a free-table run (its margin = mstart)
        z0 = np.load(os.path.expanduser(opt["dv0"]))
        dv, mstart = z0["dv"], float(z0["margin"])
        log(ev="dv0", file=opt["dv0"], margin_free=mstart)
    for rnd in range(rounds):
        jobs = [(which, j, w, wd, hpath, dv) for w in range(1 << (j + 1))]
        ms = np.zeros((ncls, 1 << (j + 1)))
        with mp.get_context("fork").Pool(nproc) as pool:
            for w, col in pool.imap_unordered(DC.one_word, jobs):
                ms[:, w] = col
        tabp = os.path.join(outdir, "%s_r%d.npz" % (tag, rnd))
        np.savez_compressed(tabp, ms=ms, Hval=(ucode >> nv) + hmin, vtop=ucode & ((1 << nv) - 1), nH=nH, nv=nv, j=j,
                            wtoday=wmin, cnt=cnt, margin=mg0, dv=dv, which=which, shift=shift, kkmin=kkmin)
        if int(opt.get("free", 0)):
            # THE FREE TABLE (no AND limit): every class takes the word with the largest minslack under the current
            # angles, then the LP is re-solved: the upper bound of ANY program on these inputs, by fix rounds
            if rnd == 0:
                wcur = wmin.copy()
            if int(opt.get("inj", 1)):
                # one-to-one on the window: inside a group (same v bits, same lowest window bit) the classes take
                # DIFFERENT words: per group the assignment with the largest smallest minslack (bottleneck), then
                # the largest sum; today's map is such an assignment (checked)
                from scipy.optimize import linear_sum_assignment
                Hc = (ucode >> nv) + hmin
                gkey = ((ucode & ((1 << nv) - 1)) << 1) | (Hc & 1)
                wnext = wcur.copy()
                nw = ms.shape[1]
                for gk in np.unique(gkey):
                    idx = np.flatnonzero(gkey == gk)
                    capw = 2 if which == "MAIN" else 1              # MAIN ends with two leftover wires
                    cnt_ = np.bincount(wmin[idx], minlength=nw)
                    assert cnt_.max() <= capw, "today's map takes a word %d times in a group" % cnt_.max()
                    sub = np.repeat(ms[idx], capw, axis=1)          # classes x (words, each capw times)
                    vals_ = np.unique(sub[np.isfinite(sub)])
                    lo_, hi_ = 0, len(vals_) - 1
                    cur_floor = float(ms[idx][np.arange(len(idx)), wcur[idx]].min())
                    best_assign = wcur[idx]

                    def try_thr(thr):
                        cost = np.where(sub >= thr, -np.minimum(sub, max(thr, 1.0) * 4) / (abs(thr) + 1.0), 1e9)
                        r_, c_ = linear_sum_assignment(cost)
                        if (cost[r_, c_] >= 1e8).any():
                            return None
                        out_ = np.zeros(len(idx), dtype=np.int64)
                        out_[r_] = c_ // capw
                        return out_
                    while lo_ <= hi_:                               # binary search of the bottleneck value
                        mid = (lo_ + hi_) // 2
                        a_ = try_thr(vals_[mid])
                        if a_ is not None:
                            if vals_[mid] > cur_floor + 1.0:
                                best_assign, cur_floor = a_, float(vals_[mid])
                            lo_ = mid + 1
                        else:
                            hi_ = mid - 1
                    wnext[idx] = best_assign
            else:
                best_w = np.argmax(ms, axis=1)
                keep = ms[np.arange(ncls), wcur] >= ms[np.arange(ncls), best_w] - 1.0
                wnext = np.where(keep, wcur, best_w)
            floor_ = float(ms[np.arange(ncls), wnext].min())
            wnew = wnext[cid]
            word = B[key] - vt + DC.wvals(which, wnew, j)
            err, fire, bad = DC.downstream(key, word, B["d"], B["k"], hpath, B["tgt"])
            if bad.any():
                log(ev="free", round=rnd, bad=int(bad.sum()))
                break
            m2, c2, b2, dv2 = DC.lp_dv(err, B["tol"], fire)
            full = (B["d"].astype(np.int64) << 23) | (word + (1 << 21))
            log(ev="free", round=rnd, classes_moved=int((wnext != wmin).sum()), inputs_changed=int((wnew != wt).sum()),
                floor_old_angles=floor_, margin=m2, conv=c2, ratio=round(m2 / mg0, 4),
                distinct_q_new=int(len(np.unique(full))))
            np.save(os.path.join(outdir, "%s_free_r%d.npy" % (tag, rnd)), wnext)
            np.savez(os.path.join(outdir, "%s_free_dv.npz" % tag), dv=dv2, margin=m2, words=wnext)
            if cur is not None and m2 <= cur * 1.002:
                cur = max(cur, m2)
                break
            cur, dv, wcur = m2, dv2, wnext
            continue
        if cur is None and mstart is not None:
            grid = [mstart * f for f in (0.999, 0.97, 0.94, 0.9, 0.86, 0.82, 0.78, 0.74, 0.7, 0.6, 0.5, 0.25, 0.0)]
        elif cur is None:
            grid = [mg0 * f for f in (0.999, 0.74, 0.5, 0.25, 0.0, -0.25, -0.5, -0.75, -1.0, -1.5, -2.0, -3.0, -5.0)]
        else:
            grid = [cur + mg0 * f for f in (2.0, 1.5, 1.0, 0.75, 0.5, 0.35, 0.25, 0.15, 0.1, 0.05, 0.02, 0.005)]
        procs = []
        for gi, M0 in enumerate(grid):
            pp = os.path.join(outdir, "%s_r%d_g%d.prog.json" % (tag, rnd, gi))
            if os.path.exists(pp):
                os.remove(pp)
            cmd = ["timeout", str(tlim), sys.executable, os.path.join(HERE, "r54_inv_tail_sat.py"), tabp, str(K),
                   "%.1f" % M0, "out=" + pp, "solver=" + solver, "keep0=%d" % keep0]
            cmd += ["%s=%s" % (k_, opt[k_]) for k_ in ("fam", "digits", "fixops", "taper", "pol", "inj") if k_ in opt]
            procs.append((M0, pp, subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)))
        res = []
        for M0, pp, p in procs:
            out = p.communicate()[0]
            verdict = "SAT" if ": SAT" in out else ("UNSAT" if ": UNSAT" in out else "timeout")
            res.append((M0, verdict, pp))
        log(ev="sat", round=rnd, K=K, grid=[[round(M0), v] for M0, v, _ in res])
        hits = [(M0, pp) for M0, v, pp in res if v == "SAT"]
        if not hits:
            log(ev="stop", round=rnd, note="no M0 of the grid is SAT")
            break
        # every hit of the grid goes through the re-solved LP (a lower M0 can re-solve higher); the best is kept
        hits.sort(reverse=True)
        if not int(opt.get("multi", 1)):
            hits = hits[:1]
        best = None
        seen_w = set()
        for M0, pp in hits:
            prog = json.load(open(pp))
            wnew = apply_prog(prog, H, vtop, nv, nH, j)
            sig = hash(wnew.tobytes())
            if sig in seen_w:
                continue
            seen_w.add(sig)
            word = B[key] - vt + DC.wvals(which, wnew, j)
            err, fire, bad = DC.downstream(key, word, B["d"], B["k"], hpath, B["tgt"])
            if bad.any():
                log(ev="eval", round=rnd, M0=M0, bad=int(bad.sum()), note="T or the word leaves its range")
                continue
            sl = DC.slack_of(err, B["tol"], fire, dv)
            m_old = float(sl.min())
            m2, c2, b2, dv2 = DC.lp_dv(err, B["tol"], fire)
            full = (B["d"].astype(np.int64) << 23) | (word + (1 << 21))
            log(ev="eval", round=rnd, M0=M0, prog=os.path.basename(pp), changed=int((wnew != wt).sum()),
                margin_old_angles=m_old, margin=m2, conv=c2, lp_bound=b2, distinct_q_new=int(len(np.unique(full))))
            if best is None or m2 > best[0]:
                best = (m2, dv2, pp)
        if best is None:
            log(ev="stop", round=rnd, note="no usable hit")
            break
        m2, dv2, pp = best
        pj = json.load(open(pp))
        pj["meta"] = dict(which=which, j=j, extra=extra, nH=nH, nv=nv, margin=m2, base_progs=sorted(progs))
        json.dump(pj, open(os.path.join(outdir, "%s_r%d.prog.json" % (tag, rnd)), "w"), indent=1)
        json.dump(pj, open(os.path.join(outdir, "%s_best.prog.json" % tag), "w"), indent=1)
        os.remove(pp)
        if cur is not None and m2 <= cur + 1:
            log(ev="stop", round=rnd, note="no gain over the previous round")
            break
        cur, dv = m2, dv2
        log(ev="round", round=rnd, margin=cur, prog="%s_r%d.prog.json" % (tag, rnd))
    log(ev="done", best=cur, base=mg0)


if __name__ == "__main__":
    main()
