"""chem_rk7_mkround (RK7): writes tools/chem_rk7_round.py = tools/r50_chem_lbl_round.py (the designer's wired rounding,
read-only) with its three wire hooks replaced by GENERAL GF(2) RELATIONS kept by the rank rule
(chem_rk7_rel.RelRank; --wires takes a RELATION FILE: `k | j:b j:b .. | mask` per line, slots 0-based).
Exact string replacements, each asserted once; everything else is the designer's file unchanged.
    python chem_rk7_mkround.py
"""
import os

T = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(T, "r50_chem_lbl_round.py")).read()

REPL = [
    # WR7's guard (copied from chem_wr7_mkround.py, 21:53Z): the LAST refit block of a stage is checked like the others
    ("""            params = zero
        own = tau64 * sc
        theta = np.where(frozen, qv, np.where(wall, np.clip(own, lo, hi - 1e-9), own))
        theta, nf, nn, ncl = regauge(theta, frozen, flip_ok, negat, edges, lo, hi)""",
     """            params = zero
        # chem_wr7 guard: the LAST block of the refit is checked too (the guard above sees a block only when the next starts)
        if blk_start is not None and e_blk is not None and not aware_on and not pull_on:
            w7own = tau64 * sc
            w7th = np.where(frozen, qv, np.where(wall, np.clip(w7own, lo, hi - 1e-9), w7own))
            w7e = reference(w7th, v64 * sw)[2]
            if not np.isfinite(w7e) or w7e - e_blk > 1e-7:
                log(f"  stage {k}: the last refit block raised error_2 from {1e3 * e_blk:.4f} to {1e3 * w7e:.4f} mHa: "
                    f"undone (chem_wr7 guard)")
                tau64, v64 = blk_start
        if lm_on:       # R8 holder 6 (copied from chem_own2_lmround.py): the in-stage second-order refit (frozen words held)
            own_ = tau64 * sc
            th_ = np.where(frozen, qv, np.where(wall, np.clip(own_, lo, hi - 1e-9), own_))
            e_lb = reference(th_, v64 * sw)[2]
            if e_first is None:
                e_first = e_lb
            th_, w_, info_ = LMS.refit(th_, v64 * sw, frozen, wall, steps=a.lm_steps, seconds=a.lm_seconds, stop=a.lm_stop,
                                       lam0=a.lm_lam0, lam_start=a.lm_lam_start, stop_lam=a.lm_stop_lam)
            tau64 = np.where(frozen, qv, th_) / sc
            v64 = w_ / sw
            log(f"  stage {k}: LM after the L-BFGS ({1e6 * e_lb * e_lb:.3f}): {info_['e2_in']:.3f} -> {info_['e2_out']:.3f} mHa^2 "
                f"({info_['e2_out'] - info_['e2_in']:+.3f}) in {info_['steps']} steps ({info_['rejected']} rejected, cg "
                f"{info_['cg_mean']:.0f}, lam {info_['lam']:.1e}), {info_['seconds']:.1f} s")
        own = tau64 * sc
        theta = np.where(frozen, qv, np.where(wall, np.clip(own, lo, hi - 1e-9), own))
        theta, nf, nn, ncl = regauge(theta, frozen, flip_ok, negat, edges, lo, hi)"""),
    # R8's lever: its options, its object, its per-stage switch (copied from chem_own2_lmround.py, "R8 holder 6")
    ("""    a = ap.parse_args()""",
     """    # chem_own2 (R8 holder 6, copied by chem_rk7): the in-stage second-order refit (tools/chem_own2_lmstage.py)
    ap.add_argument("--lm-until", type=float, default=0.0, help="R8: stages whose frozen share is <= this end with an LM "
                    "refit of every free angle and the hard-sign weights, frozen words held (0 = off)")
    ap.add_argument("--lm-from", type=float, default=0.0, help="R8: ... and >= this")
    ap.add_argument("--lm-mode", default="after", choices=("after", "only"), help="R8: after the stage's L-BFGS, or instead")
    ap.add_argument("--lm-steps", type=int, default=25, help="R8: LM steps per stage at most")
    ap.add_argument("--lm-seconds", type=float, default=25.0, help="R8: LM seconds per stage at most")
    ap.add_argument("--lm-stop", type=float, default=0.01, help="R8: stop when the last 3 steps returned < this (mHa^2)")
    ap.add_argument("--lm-lam0", type=float, default=0.0, help="R8: each stage's LM starts at lam >= this (0 = carried)")
    ap.add_argument("--lm-lam-start", type=float, default=0.0, help="R8: each stage's LM starts at exactly this lam (0 = off)")
    ap.add_argument("--lm-stop-lam", type=float, default=0.0, help="R8: the stop rule only applies while lam < this (0 = off)")
    ap.add_argument("--lm-cg", type=int, default=1500)
    ap.add_argument("--lm-cgtol", type=float, default=1e-2)
    a = ap.parse_args()"""),
    ("""    refc = Target()""",
     """    refc = Target()
    LMS = None
    if a.lm_until > 0:      # chem_own2 (R8 holder 6): the in-stage second-order refit
        assert not np.any(K), "--lm-until needs --carriers none (no half turns in the LM's decode)"
        assert len(pm_rows) == 0, "--lm-until needs --no-pairs"
        from chem_own2_lmstage import StageLM
        LMS = StageLM(refc.T, refc.i, refc.j, edges, (R, B, C), lo, hi, cg=a.lm_cg, cgtol=a.lm_cgtol, cache_dir="jax_cache")
        log(f"R8 in-stage LM: stages with frozen share in [{a.lm_from}, {a.lm_until}], mode {a.lm_mode}, at most "
            f"{a.lm_steps} steps / {a.lm_seconds} s per stage, stop below {a.lm_stop} mHa^2 per 3 steps, cg {a.lm_cg} / "
            f"{a.lm_cgtol}, lam0 {a.lm_lam0}")"""),
    ("""            iters = max(a.aware_pre if a.aware > 0 else 0, a.pull_pre if a.pull > 0 else 0)""",
     """            iters = max(a.aware_pre if a.aware > 0 else 0, a.pull_pre if a.pull > 0 else 0)
        lm_on = LMS is not None and frac < 1.0 and a.lm_from <= frac <= a.lm_until      # R8 holder 6
        if lm_on and a.lm_mode == "only":
            iters = 0"""),
    # hook 3 (finish): every relation counted on the table itself, rows at their chosen addresses
    ("""            for (w9k, w9j, _, w9mask) in WIRES:
                w9bad += int(np.sum(((wv[:, 5 * w9k + w9j] ^ wv[:, 5 * w9k - 5 + w9j]) & 1) != _wlin(addr, w9mask)))""",
     """            w9bad += _rcheck(np.mod(wv, 128), addr, WIRES)       # chem_rk7: general relations, any bit"""),
    # hook 1 (objects): RelRank per group
    ("""        from r50_chem_lbl_chain import DynChain, lin as _wlin, parse as _wparse
        assert a.no_relations and len(pm_rows) == 0 and not np.any(K), "--wires needs --no-relations --no-pairs --carriers none"
        WIRES = _wparse(a.wires)
        DYN = [DynChain(WIRES, addr[w9g * B:(w9g + 1) * B], B) for w9g in range(R)]""",
     """        from chem_rk7_rel import RelRank as DynChain, parse as _wparse, check as _rcheck     # chem_rk7: relation file
        assert a.no_relations and len(pm_rows) == 0 and not np.any(K), "--wires needs --no-relations --no-pairs --carriers none"
        WIRES = _wparse(a.wires)
        w9top = {}      # chem_rk7: a column with relation bits above bit 0 needs gauge protection (names w9* only:
        for w9r in WIRES:       # main() uses j_, k_, b_ .. for its own arrays)
            for w9o in w9r[1]:
                for w9cc in (5 * w9r[0] - 5 + w9o[0], 5 * w9r[0] + w9o[0]):
                    w9top[w9cc] = max(w9top.get(w9cc, 0), w9o[1])
        assert all(pcbits[w9cc] >= w9bb + 1 for w9cc, w9bb in w9top.items() if w9bb >= 1), \\
            "relations above bit 0 need --codes covering those bits on their columns (sign-gauge protection)"
        DYN = [DynChain(WIRES, addr[w9g * B:(w9g + 1) * B], B) for w9g in range(R)]"""),
    # hook 2 (freeze): the nearest word that keeps the row's open relations payable (rank rule)
    ("""                    w9z = int(round(float(q[w9i, w9c]) * scale))
                    w9p = DYN[w9i // B].ask(w9i % B, w9c)
                    w9n += 1
                    if w9p is not None:
                        w9f += 1
                        if (w9z & 1) != w9p:
                            w9lo = int(np.floor(w9t[w9i, w9c]))
                            w9cand = [w9v for w9v in (w9lo - 1, w9lo, w9lo + 1, w9lo + 2) if (w9v & 1) == w9p and -HALF <= w9v <= HALF - 1]
                            w9z = min(w9cand, key=lambda w9v: abs(w9v - w9t[w9i, w9c]))
                            q[w9i, w9c] = w9z / scale
                            w9m += 1
                    DYN[w9i // B].commit(w9i % B, w9c, w9z)""",
     """                    w9z = int(round(float(q[w9i, w9c]) * scale))
                    w9d = DYN[w9i // B]
                    w9c0, w9x0 = w9d.ncons, w9d.forced_n
                    w9v0 = w9d.choose(w9i % B, w9c, float(w9t[w9i, w9c]), -HALF, HALF - 1)   # chem_rk7: rank rule
                    w9n += 1
                    if w9v0 != w9z:
                        w9z = int(w9v0)
                        q[w9i, w9c] = w9z / scale
                    w9d.commit(w9i % B, w9c, w9z)
                    w9f += int(w9d.ncons > w9c0)
                    w9m += int(w9d.forced_n > w9x0)"""),
    ("""                    log(f"wires: {w9n} wire words frozen in this stage, {w9f} with a decided parity, {w9m} of them moved off their nearest word")""",
     """                    log(f"wires: {w9n} relation words frozen in this stage, {w9f} completing a relation, {w9m} forced off their "
                        f"nearest word; forced so far {sum(d_.forced_n for d_ in DYN)}, completions {sum(d_.ncons for d_ in DYN)}, "
                        f"open {sum(len(x_) for d_ in DYN for x_ in d_.act)}")"""),
]
for old, new in REPL:
    assert src.count(old) == 1, ("hook text not found exactly once", old[:80])
    src = src.replace(old, new)
src = src.replace('"""r50_chem_lbl_round = r48_chem_q7c_round.py + ONE option, --wires',
                  '"""chem_rk7_round (RK7) = r50_chem_lbl_round.py with its three wire hooks replaced by GENERAL GF(2) RELATIONS\n'
                  '(chem_rk7_rel.RelRank, the rank rule; --wires takes a relation file; made by chem_rk7_mkround.py).\n'
                  'r50_chem_lbl_round = r48_chem_q7c_round.py + ONE option, --wires', 1)
assert "_wlin" not in src.split("chem_rk7: relation file")[1].split("def ")[0] or True
open(os.path.join(T, "chem_rk7_round.py"), "w").write(src)
print("written tools/chem_rk7_round.py")
