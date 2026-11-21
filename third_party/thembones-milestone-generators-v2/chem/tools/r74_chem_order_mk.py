#!/usr/bin/env python
"""r74_chem_order_mk (design seat R74-CHEM-ORDER): generates tools/r74_chem_order_fround.py from T10's
tools/chem_t10_froundgjc.py (never edited) by ASSERTED replacements. Added options (all off by default: without them
the generated tool does what T10's tool does):

  --seq-m M          BABAI INSIDE THE STAGE: the words of one stage are decided one after the other (rank order): each
                     takes the allowed word nearest to its present value SHIFTED by what the earlier picks of the
                     stage do to it in the group's linear model (the block Sigma_MM of the conditional table's own
                     inverse on the M stiffest free words of the group, tools/r74_chem_order_cond.py; conditioned on
                     the block's words frozen since the table). The refit, the placer and the order are untouched; a
                     pick outside the block is rounded as before.
  --schedule fine8 | fine16   every stage of fine4 split in 2 / 4 (smaller batches; give --iters accordingly). A state
                     kept under another schedule is resumed at the stage its frozen share stands at.
  --cond-from S      the conditional order only from the frozen share S on (the plain S_kk order before): for runs
                     from zero with a switch
  --fp-cost cond     the FreezePlacer's costs are the CONDITIONAL stiffness (the rank value of --rank-cond) instead of
                     S_kk
  --cond-sync S      A TABLE OF AGE ZERO IN THE LEVER PHASE: while the frozen share is below S every stage WAITS for a
                     conditional table computed from ITS OWN present state before it picks (the workers' 110 to 120 s
                     a stage; today a stage is ordered by a table 2.3 to 3.6 stages old). After S the tables come
                     as before (asynchronously).
  --fresh-m M        A FRESH RANK EVERY STAGE: the conditional table in use is 2 to 8 stages old (its workers take 100 to
                     200 s); with the block of its inverse on the M stiffest free words a group the conditional
                     stiffness of those words is CONDITIONED EXACTLY on every block word frozen since the table's state
                     (d' = 1 / diag(Sigma_PP - Sigma_PQ Sigma_QQ^-1 Sigma_QP)) before every pick: a rank as fresh as
                     the last freeze, in the linear model of the table's state. Nothing else changes.
  --opp L --opp-wait W [--opp-until S]   OPPORTUNISTIC FILL ON THE CONDITIONAL RANK: the candidates of a stage are the
                     first L x need words of the rank; a candidate that has waited W stages is taken; the rest of the
                     stage takes the candidates whose rounding costs least NOW (conditional stiffness x offset^2)
  --rank-walpha A    wire-aware rank: a free wired word whose parity is forced OFF its nearest word NOW has its rank
                     value x (1 + A) (decided earlier); --rank-wdef A: x 1 / (1 + A) (deferred)

    python tools/r74_chem_order_mk.py tools/chem_t10_froundgjc.py tools/r74_chem_order_fround.py
"""
import sys

src, dst = sys.argv[1], sys.argv[2]
s = open(src, encoding="utf-8").read()


def rep(a, b, n=1):
    global s
    assert s.count(a) == n, (s.count(a), a[:90])
    s = s.replace(a, b)


rep('"""T10 copy (tools/chem_t10_mkgjc.py: + R67\'s --rank-cond)',
    '"""R74-CHEM-ORDER copy (tools/r74_chem_order_mk.py: + --seq-m, --schedule fine8 / fine16, --cond-from, --fp-cost,\n'
    '--rank-walpha, --rank-wdef) of\nT10 copy (tools/chem_t10_mkgjc.py: + R67\'s --rank-cond)')
rep('    ap.add_argument("--lm-until", type=float, default=0.0,',
    '    ap.add_argument("--seq-m", type=int, default=0, help="r74: Babai inside the stage with the block of the conditional table on its M stiffest free words a group (0 = off)")\n'
    '    ap.add_argument("--decor", type=float, default=0.0, help="r74: decorrelated stages: a candidate whose conditional correlation with a word already taken in the stage exceeds this waits (0 = off; needs --decor-m)")\n'
    '    ap.add_argument("--decor-m", type=int, default=0, help="r74: block size a group for --decor (the M stiffest free words by the conditional table)")\n'
    '    ap.add_argument("--decor-l", type=float, default=4.0, help="r74: --decor looks at the first L x need words of the rank")\n'
    '    ap.add_argument("--fine4-from", type=float, default=2.0, help="r74: with --schedule fine8 / fine16: the stages beyond this frozen share are fine4 stages")\n'
    '    ap.add_argument("--cond-sync", type=float, default=0.0, help="r74: below this frozen share every stage waits for a conditional table of its own present state (0 = off)")\n'
    '    ap.add_argument("--fresh-m", type=int, default=0, help="r74: the conditional rank conditioned every stage on the block words frozen since the table (block of the M stiffest free words a group; 0 = off)")\n'
    '    ap.add_argument("--harv", type=float, default=0.0, help="r74: harvest: this share of every stage is taken from the next L x need words of the conditional rank by the smallest cost now; the rest is the plain top of the rank (0 = off)")\n'
    '    ap.add_argument("--harv-l", type=float, default=3.0, help="r74: --harv looks at the next L x need words of the rank")\n'
    '    ap.add_argument("--harv-until", type=float, default=1.0, help="r74: --harv only below this frozen share")\n'
    '    ap.add_argument("--opp", type=float, default=0.0, help="r74: opportunistic fill: candidates = the first L x need words of the conditional rank (0 = off)")\n'
    '    ap.add_argument("--opp-wait", type=int, default=8, help="r74: a candidate that has waited this many stages is taken")\n'
    '    ap.add_argument("--opp-until", type=float, default=1.0, help="r74: the opportunistic fill only below this frozen share")\n'
    '    ap.add_argument("--fp-cost", default="sk", choices=("sk", "cond"), help="r74: weights of the FreezePlacer\'s costs: S_kk (the tool) or the conditional stiffness of --rank-cond")\n'
    '    ap.add_argument("--cond-from", type=float, default=0.0, help="r74: --rank-cond orders the freeze only from this frozen share on (S_kk before)")\n'
    '    ap.add_argument("--rank-walpha", type=float, default=0.0, help="r74: rank value x (1 + A) for a free wired word forced off its nearest word now")\n'
    '    ap.add_argument("--rank-wdef", type=float, default=0.0, help="r74: rank value / (1 + A) for a free wired word forced off its nearest word now")\n'
    '    ap.add_argument("--lm-until", type=float, default=0.0,')
# finer schedules; a state of another schedule resumes at the stage of its frozen share
rep('                    if cond_job is None:\n                        cond_job, cond_from = cond_start(now.copy(), w.copy(), frozen.copy()), k\n',
    '                    if cond_job is None and not r74sy:      # r74: no second table of the same state under --cond-sync\n                        cond_job, cond_from = cond_start(now.copy(), w.copy(), frozen.copy()), k\n')
rep('choices=["full", "long", "fine2", "fine4"],', 'choices=["full", "long", "fine2", "fine4", "fine8", "fine16"],')
rep('            out_ += [round(prev_ + (f_ - prev_) * j_ / kf, 7) for j_ in range(1, kf + 1)]\n',
    '            r74kf = 4 if (kf > 4 and f_ > a.fine4_from + 1e-9) else kf      # r74: --fine4-from\n'
    '            out_ += [round(prev_ + (f_ - prev_) * j_ / r74kf, 7) for j_ in range(1, r74kf + 1)]\n')
rep('        log(f"resume at stage {stage0} (frozen {frozen.mean():.4f}, {spent / 60:.1f} minutes spent so far)")\n',
    '        if stage0 > 0 and abs(float(sched[min(stage0, len(sched)) - 1]) - float(frozen.mean())) > 2e-4:      # r74: a state of another schedule\n'
    '            r74s0 = int(np.searchsorted(np.array(sched), float(frozen.mean()) + 1e-9))\n'
    '            log(f"r74: the state stands at frozen {frozen.mean():.4f}, stage {stage0} of another schedule: resumed at stage {r74s0} of {a.schedule}")\n'
    '            stage0 = r74s0\n'
    '        log(f"resume at stage {stage0} (frozen {frozen.mean():.4f}, {spent / 60:.1f} minutes spent so far)")\n')
# the worker of my module (R67's + the block)
rep('        from r67_chem_path_cond import group_cond as _group_cond\n',
    '        from r74_chem_order_cond import group_cond as _group_cond      # r74: R67\'s worker + the block (M > 0)\n'
    '        r74_blk = [None] * R\n')
rep("                                                         cond_cpus[g_ % len(cond_cpus)] if cond_cpus else None)) for g_ in range(R)], time.time()\n",
    "                                                         cond_cpus[g_ % len(cond_cpus)] if cond_cpus else None, max(a.seq_m, a.decor_m, a.fresh_m))) for g_ in range(R)], time.time()\n")
rep('            return np.concatenate([j_.get()[0] for j_ in job_[0]], axis=0)\n',
    '            res_ = [j_.get() for j_ in job_[0]]\n'
    '            if max(a.seq_m, a.decor_m, a.fresh_m) > 0:      # r74: the blocks of this table (kept only when every group came back)\n'
    '                r74_blk[:] = [(r_[2], r_[3]) for r_ in res_]\n'
    '            return np.concatenate([r_[0] for r_ in res_], axis=0)\n')
# the conditional order from a share on
rep('                if COND:        # r67: the order by the conditional stiffness (the latest finished table)',
    '                if COND:        # r67: the order by the conditional stiffness (the latest finished table)\n'
    '                    r74_on = float(frozen.mean()) >= a.cond_from - 1e-9      # r74: --cond-from\n'
    '                    r74sy = a.cond_sync > 0 and r74_on and float(frozen.mean()) < a.cond_sync - 1e-9      # r74: --cond-sync\n'
    '                    if r74sy:      # every stage of the lever phase waits for the table of ITS OWN state\n'
    '                        if cond_job is not None:      # a table of an older state is under way: its workers finish, its result is dropped\n'
    '                            try:\n'
    '                                cond_get(cond_job)\n'
    '                            except Exception:      # noqa: BLE001\n'
    '                                pass\n'
    '                            cond_job = None\n'
    '                        r74t0 = time.time()\n'
    '                        try:\n'
    '                            rank_cd = cond_get(cond_start(now.copy(), w.copy(), frozen.copy()))\n'
    '                            log(f"r74 sync: the conditional table of the present state (start of stage {k}) computed and in use ({time.time() - r74t0:.0f} s)")\n'
    '                        except Exception as ex_:      # noqa: BLE001\n'
    '                            log(f"r74 sync: WARNING the table failed ({ex_!r}); the last one stays in use")')
rep('                    val_c = rank_cd if a.rank_geo <= 0 else (np.maximum(rank_sk, 1e-300) ** a.rank_geo) * (np.maximum(rank_cd, 1e-300) ** (1.0 - a.rank_geo))\n',
    '                    r74cd = rank_cd\n'
    '                    if a.fresh_m > 0 and r74_blk[0] is not None and r74_on:      # r74: the rank conditioned on the block words frozen since the table\n'
    '                        r74cd = rank_cd.copy()\n'
    '                        r74nq = 0\n'
    '                        for r74gg in range(R):\n'
    '                            r74i, r74K = r74_blk[r74gg]\n'
    '                            if len(r74i) == 0:\n'
    '                                continue\n'
    '                            r74fl = r74gg * B * (N - 1) + r74i\n'
    '                            r74fz = frozen.reshape(-1)[r74fl]\n'
    '                            r74Q, r74P = np.flatnonzero(r74fz), np.flatnonzero(~r74fz)\n'
    '                            if len(r74Q) == 0 or len(r74P) == 0:\n'
    '                                continue\n'
    '                            r74Kq = r74K[np.ix_(r74Q, r74Q)].copy()\n'
    '                            r74Kq[np.diag_indices_from(r74Kq)] += 1e-12 * float(np.trace(r74Kq)) / len(r74Q)\n'
    '                            r74Kpq = r74K[np.ix_(r74P, r74Q)]\n'
    '                            r74dg = np.diag(r74K)[r74P] - np.sum(r74Kpq * np.linalg.solve(r74Kq, r74Kpq.T).T, axis=1)\n'
    '                            r74cd.reshape(-1)[r74fl[r74P]] = np.maximum(1.0 / np.maximum(r74dg, 1e-300), rank_cd.reshape(-1)[r74fl[r74P]])\n'
    '                            r74nq += len(r74Q)\n'
    '                        r74st = np.argsort(-np.where(frozen, -np.inf, rank_cd), axis=None)[:need]\n'
    '                        r74fr = np.argsort(-np.where(frozen, -np.inf, r74cd), axis=None)[:need]\n'
    '                        log(f"r74 fresh: the rank conditioned on {r74nq} block words frozen since the table; {int(np.isin(r74fr, r74st).sum())} of the {need} picks are the stale table\'s picks")\n'
    '                    val_c = r74cd if a.rank_geo <= 0 else (np.maximum(rank_sk, 1e-300) ** a.rank_geo) * (np.maximum(r74cd, 1e-300) ** (1.0 - a.rank_geo))\n')
rep('                    impact = val_c + 1e-30 * c\n',
    '                    if r74_on:      # r74\n'
    '                        impact = val_c + 1e-30 * c\n')
# the opportunistic fill
rep('            else:\n                pick = np.argsort(-impact, axis=None)[:need]\n',
    '            else:\n                pick = np.argsort(-impact, axis=None)[:need]\n'
    '                if a.opp > 1 and COND and rank_cd is not None and a.cond_from - 1e-9 <= float(frozen.mean()) < a.opp_until:      # r74: opportunistic fill\n'
    '                    if "r74_wait" not in globals():\n'
    '                        globals()["r74_wait"] = np.zeros(theta.size, np.int64)\n'
    '                    r74w = globals()["r74_wait"]\n'
    '                    r74c = np.argsort(-impact, axis=None)[:min(max(need, int(round(a.opp * need))), int((~frozen).sum()))]\n'
    '                    r74must = r74c[r74w[r74c] >= a.opp_wait][:need]\n'
    '                    r74rest = np.setdiff1d(r74c, r74must, assume_unique=True)\n'
    '                    r74sc = (((q - now) * scale).reshape(-1)[r74rest] ** 2) * np.maximum(rank_cd.reshape(-1)[r74rest], 0.0)\n'
    '                    pick = np.concatenate([r74must, r74rest[np.argsort(r74sc, kind="stable")[:need - len(r74must)]]]).astype(np.int64)\n'
    '                    r74w[r74c] += 1\n'
    '                    log(f"r74 opp: {len(r74c)} candidates, {len(r74must)} taken after waiting {a.opp_wait} stages, the rest by the smallest cost now; mean squared offset of the stage words {float(np.mean((((q - now) * scale).reshape(-1)[pick]) ** 2)):.4f} steps^2")\n'
    '                if a.harv > 0 and COND and rank_cd is not None and a.cond_from - 1e-9 <= float(frozen.mean()) < a.harv_until:      # r74: harvest\n'
    '                    r74n1 = int(round((1.0 - a.harv) * need))\n'
    '                    r74o = np.argsort(-impact, axis=None)\n'
    '                    r74pool = r74o[r74n1:min(r74n1 + int(round(a.harv_l * need)), int((~frozen).sum()))]\n'
    '                    r74sc = (((q - now) * scale).reshape(-1)[r74pool] ** 2) * np.maximum(rank_cd.reshape(-1)[r74pool], 0.0)\n'
    '                    pick = np.concatenate([r74o[:r74n1], r74pool[np.argsort(r74sc, kind="stable")[:need - r74n1]]]).astype(np.int64)\n'
    '                    log(f"r74 harv: {r74n1} words by rank, {need - r74n1} harvested from the next {len(r74pool)}; mean squared offset of the harvested words {float(np.mean((((q - now) * scale).reshape(-1)[pick[r74n1:]]) ** 2)) if need > r74n1 else 0.0:.4f}, of the ranked words {float(np.mean((((q - now) * scale).reshape(-1)[pick[:r74n1]]) ** 2)) if r74n1 else 0.0:.4f} steps^2")\n'
    '                if a.decor > 0 and a.opp <= 1 and COND and r74_blk[0] is not None and float(frozen.mean()) >= a.cond_from - 1e-9:      # r74: decorrelated stages\n'
    '                    r74c = np.argsort(-impact, axis=None)[:min(max(need, int(round(a.decor_l * need))), int((~frozen).sum()))]\n'
    '                    r74tk, r74bp, r74ps, r74dg, r74nd = [], [[] for _ in range(R)], [None] * R, [None] * R, 0\n'
    '                    for r74ff in r74c:\n'
    '                        if len(r74tk) >= need:\n'
    '                            break\n'
    '                        r74gg = (int(r74ff) // (N - 1)) // B\n'
    '                        if r74ps[r74gg] is None:\n'
    '                            r74ps[r74gg] = {int(r74vv): r74nn for r74nn, r74vv in enumerate(r74_blk[r74gg][0])}\n'
    '                            r74dg[r74gg] = np.sqrt(np.maximum(np.diag(r74_blk[r74gg][1]), 1e-300))\n'
    '                        r74ll = r74ps[r74gg].get(int(r74ff) - r74gg * B * (N - 1), -1)\n'
    '                        if r74ll >= 0 and r74bp[r74gg]:\n'
    '                            r74tt = np.array(r74bp[r74gg])\n'
    '                            if float(np.max(np.abs(r74_blk[r74gg][1][r74ll, r74tt]) / (r74dg[r74gg][r74ll] * r74dg[r74gg][r74tt]))) > a.decor:\n'
    '                                r74nd += 1\n'
    '                                continue\n'
    '                        r74tk.append(int(r74ff))\n'
    '                        if r74ll >= 0:\n'
    '                            r74bp[r74gg].append(r74ll)\n'
    '                    if len(r74tk) < need:\n'
    '                        r74in = set(r74tk)\n'
    '                        r74tk += [int(r74ff) for r74ff in r74c if int(r74ff) not in r74in][:need - len(r74tk)]\n'
    '                    r74p0 = pick\n'
    '                    pick = np.array(r74tk[:need], np.int64)\n'
    '                    log(f"r74 decor: {r74nd} candidates wait (conditional correlation over {a.decor} with a word taken in this stage); {int(np.isin(pick, r74p0).sum())} of the {need} picks are the plain top of the rank")\n')
# the placer's costs; the wire-aware rank (the pick is taken again after the placer's residues)
rep('                fpc = (rank_sk + 1e-30 * c) if (a.rank in ("sdiag", "sdyn") and rank_sk is not None) else c\n'
    '                fpr, fpm = FP.residues(fpt, fpc, frozen, np.round(qv * scale).astype(np.int64))\n',
    '                fpc = (rank_sk + 1e-30 * c) if (a.rank in ("sdiag", "sdyn") and rank_sk is not None) else c\n'
    '                if a.fp_cost == "cond" and COND and rank_cd is not None and float(frozen.mean()) >= a.cond_from - 1e-9:      # r74\n'
    '                    fpc = np.maximum(rank_cd, 0.0) + 1e-30 * c\n'
    '                fpr, fpm = FP.residues(fpt, fpc, frozen, np.round(qv * scale).astype(np.int64))\n'
    '                if a.rank_walpha > 0 or a.rank_wdef > 0:      # r74: the wire-aware rank: the pick is taken again\n'
    '                    r74z = fpr + 2 * np.round((fpt - fpr) / 2)\n'
    '                    r74z = np.where(r74z > HALF - 1, r74z - 2, np.where(r74z < -HALF, r74z + 2, r74z))\n'
    '                    r74off = (fpm == 2)[None, :] & ~frozen & (r74z != np.round(q * scale))\n'
    '                    r74f = (1.0 + a.rank_walpha) if a.rank_walpha > 0 else 1.0 / (1.0 + a.rank_wdef)\n'
    '                    impact = np.where(r74off & (impact > 0), impact * r74f, impact)\n'
    '                    pick = np.argsort(-impact, axis=None)[:need]\n'
    '                    log(f"r74 wire-aware rank: {int(r74off.sum())} free wired words forced off their nearest word now, {int(r74off.flat[pick].sum())} of them in this stage")\n')
# Babai inside the stage: after the placer's words, before the freeze
rep('            frozen.flat[pick] = True\n            qv.flat[pick] = q.flat[pick]\n',
    '            if a.seq_m > 0 and COND and r74_blk[0] is not None and float(frozen.mean()) >= a.cond_from - 1e-9:      # r74: Babai inside the stage\n'
    '                r74t = now * scale\n'
    '                r74n = r74m = r74o = 0\n'
    '                r74g = (pick // (N - 1)) // B\n'
    '                for r74gg in range(R):\n'
    '                    r74i, r74K = r74_blk[r74gg]\n'
    '                    r74sel = pick[r74g == r74gg]\n'
    '                    if len(r74sel) < 2 or len(r74i) == 0:\n'
    '                        continue\n'
    '                    r74pos = {int(r74vv): r74ii for r74ii, r74vv in enumerate(r74i)}\n'
    '                    r74loc = np.array([r74pos.get(int(r74ff) - r74gg * B * (N - 1), -1) for r74ff in r74sel], np.int64)\n'
    '                    r74o += int((r74loc < 0).sum())\n'
    '                    r74f, r74P = r74sel[r74loc >= 0], r74loc[r74loc >= 0]\n'
    '                    if len(r74P) < 2:\n'
    '                        continue\n'
    '                    r74Kp = r74K[np.ix_(r74P, r74P)].copy()\n'
    '                    r74Q = np.flatnonzero(frozen.reshape(-1)[r74gg * B * (N - 1) + r74i])      # block words frozen since the table\n'
    '                    if len(r74Q):\n'
    '                        r74Kq = r74K[np.ix_(r74Q, r74Q)].copy()\n'
    '                        r74Kq[np.diag_indices_from(r74Kq)] += 1e-12 * float(np.trace(r74Kq)) / len(r74Q)\n'
    '                        r74Kpq = r74K[np.ix_(r74P, r74Q)]\n'
    '                        r74Kp -= r74Kpq @ np.linalg.solve(r74Kq, r74Kpq.T)\n'
    '                    r74y = r74t.reshape(-1)[r74f].copy()\n'
    '                    for r74ii in range(len(r74P)):\n'
    '                        r74r, r74c = divmod(int(r74f[r74ii]), N - 1)\n'
    '                        if FP is not None and fpm[r74c] == 2:\n'
    '                            r74z = float(fpr[r74r, r74c] + 2 * np.round((r74y[r74ii] - fpr[r74r, r74c]) / 2))\n'
    '                            r74z = r74z - 2 if r74z > HALF - 1 else (r74z + 2 if r74z < -HALF else r74z)\n'
    '                        else:\n'
    '                            r74z = float(min(max(np.round(r74y[r74ii]), -HALF), HALF - 1))\n'
    '                        r74n += 1\n'
    '                        r74m += int(r74z != np.round(q[r74r, r74c] * scale))\n'
    '                        q[r74r, r74c] = r74z / scale\n'
    '                        r74d = float(r74Kp[r74ii, r74ii])\n'
    '                        if r74d > 1e-300:\n'
    '                            r74col = r74Kp[:, r74ii] / r74d\n'
    '                            r74y = r74y + r74col * (r74z - r74y[r74ii])\n'
    '                            r74Kp = r74Kp - np.outer(r74col, r74Kp[r74ii, :])\n'
    '                        r74y[r74ii] = r74z\n'
    '                log(f"r74 seq: {r74n} words decided one after the other inside their groups, {r74m} took another word than their own nearest; {r74o} picks outside the blocks")\n'
    '            frozen.flat[pick] = True\n            qv.flat[pick] = q.flat[pick]\n')
open(dst, "w", encoding="utf-8", newline="\n").write(s)
print("written", dst, len(s))
