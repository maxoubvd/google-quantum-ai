"""chem_wr7_mkround (WR7): writes two copies of tools/r50_chem_lbl_round.py (exact string replacements, each asserted
once; everything else is the designer's file unchanged):
  tools/chem_wr7_roundg.py = + ONE guard: the LAST refit block of a stage is checked like the others (the tool's own
      guard, R8 Nov 15, undoes a block that raised error_2 only when the NEXT block starts: run wr7_a's stage 185 went
      19.554 -> 26.326 mHa in its refit, 1,231 crossings, and kept +11.3 mHa^2 after the next stage's recovery);
  tools/chem_wr7_round.py  = the guard + the three wire hooks generalised to bits 0 AND 1 (chem_wr7_chainb.DynChainB).
    python chem_wr7_mkround.py
"""
import os

T = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(T, "r50_chem_lbl_round.py")).read()

GUARD = [
    ("""            params = zero
        own = tau64 * sc
        theta = np.where(frozen, qv, np.where(wall, np.clip(own, lo, hi - 1e-9), own))
        theta, nf, nn, ncl = regauge(theta, frozen, flip_ok, negat, edges, lo, hi)""",
     """            params = zero
        # chem_wr7: the LAST block of the refit is checked too (the guard above sees a block only when the next starts)
        if blk_start is not None and e_blk is not None and not aware_on and not pull_on:
            w7own = tau64 * sc
            w7th = np.where(frozen, qv, np.where(wall, np.clip(w7own, lo, hi - 1e-9), w7own))
            w7e = reference(w7th, v64 * sw)[2]
            if not np.isfinite(w7e) or w7e - e_blk > 1e-7:
                log(f"  stage {k}: the last refit block raised error_2 from {1e3 * e_blk:.4f} to {1e3 * w7e:.4f} mHa: "
                    f"undone (chem_wr7 guard)")
                tau64, v64 = blk_start
        own = tau64 * sc
        theta = np.where(frozen, qv, np.where(wall, np.clip(own, lo, hi - 1e-9), own))
        theta, nf, nn, ncl = regauge(theta, frozen, flip_ok, negat, edges, lo, hi)"""),
]
for old, new in GUARD:
    assert src.count(old) == 1, ("guard text not found exactly once", old[:80])
    src = src.replace(old, new)
src_g = src.replace('"""r50_chem_lbl_round = r48_chem_q7c_round.py + ONE option, --wires',
                    '"""chem_wr7_roundg (WR7) = r50_chem_lbl_round.py + the last refit block checked (chem_wr7_mkround.py).\n'
                    'r50_chem_lbl_round = r48_chem_q7c_round.py + ONE option, --wires', 1)
open(os.path.join(T, "chem_wr7_roundg.py"), "w").write(src_g)
print("written tools/chem_wr7_roundg.py")

# tools/chem_wr7_roundl.py = chem_wr7_roundg.py + R8's in-stage LM lever ("R8 holder 6", tools/chem_own2_lmround.py;
# the four insertions as RK7 copied them into chem_rk7_mkround.py), for the lever twin of run D (the lead's run (a)).
LEVER = [
    ("""                tau64, v64 = blk_start
        own = tau64 * sc""",
     """                tau64, v64 = blk_start
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
        own = tau64 * sc"""),
    ("""    a = ap.parse_args()""",
     """    # chem_own2 (R8 holder 6, copied by chem_wr7 via chem_rk7): the in-stage second-order refit (tools/chem_own2_lmstage.py)
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
]
src_l = src_g
for old, new in LEVER:
    assert src_l.count(old) == 1, ("lever text not found exactly once", old[:80])
    src_l = src_l.replace(old, new)
src_l = src_l.replace('"""chem_wr7_roundg (WR7) = r50_chem_lbl_round.py + the last refit block checked (chem_wr7_mkround.py).',
                      '"""chem_wr7_roundl (WR7) = chem_wr7_roundg.py + R8\'s in-stage LM lever (R8 holder 6; chem_wr7_mkround.py).', 1)
open(os.path.join(T, "chem_wr7_roundl.py"), "w").write(src_l)
print("written tools/chem_wr7_roundl.py")

REPL = [
    # hook 3 (finish): relations counted on their own bit
    ("""            for (w9k, w9j, _, w9mask) in WIRES:
                w9bad += int(np.sum(((wv[:, 5 * w9k + w9j] ^ wv[:, 5 * w9k - 5 + w9j]) & 1) != _wlin(addr, w9mask)))""",
     """            for (w9k, w9j, w9b, w9mask) in WIRES:     # chem_wr7: the relation on its own bit (0 or 1)
                w9bad += int(np.sum((((wv[:, 5 * w9k + w9j] ^ wv[:, 5 * w9k - 5 + w9j]) >> w9b) & 1) != _wlin(addr, w9mask)))"""),
    # hook 1 (objects): DynChainB
    ("""        from r50_chem_lbl_chain import DynChain, lin as _wlin, parse as _wparse
        assert a.no_relations and len(pm_rows) == 0 and not np.any(K), "--wires needs --no-relations --no-pairs --carriers none"
        WIRES = _wparse(a.wires)
        DYN = [DynChain(WIRES, addr[w9g * B:(w9g + 1) * B], B) for w9g in range(R)]""",
     """        from chem_wr7_chainb import DynChainB as DynChain, parse as _wparse      # chem_wr7: bits 0 and 1
        from r50_chem_lbl_chain import lin as _wlin
        assert a.no_relations and len(pm_rows) == 0 and not np.any(K), "--wires needs --no-relations --no-pairs --carriers none"
        WIRES = _wparse(a.wires)
        assert all(w9w[2] in (0, 1) for w9w in WIRES), "wires on bits 0 and 1 only"
        if any(w9w[2] == 1 for w9w in WIRES):     # chem_wr7: odd frozen words of bit-1-wired columns must not be negated
            w9c1 = sorted({c_ for (k_, j_, b_, _) in WIRES if b_ == 1 for c_ in (5 * k_ - 5 + j_, 5 * k_ + j_)})
            assert all(pcbits[c_] >= 2 for c_ in w9c1), "bit-1 wires need --codes with bits 0 and 1 on their columns (gauge protection)"
        DYN = [DynChain(WIRES, addr[w9g * B:(w9g + 1) * B], B) for w9g in range(R)]"""),
    # hook 2 (freeze): the nearest feasible word (bits 0 and 1)
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
                    w9v0 = DYN[w9i // B].choose(w9i % B, w9c, float(w9t[w9i, w9c]), -HALF, HALF - 1)   # chem_wr7
                    w9n += 1
                    if w9v0 != w9z:
                        w9f += 1
                        w9z = int(w9v0)
                        q[w9i, w9c] = w9z / scale
                        w9m += 1
                    DYN[w9i // B].commit(w9i % B, w9c, w9z)"""),
]
for old, new in REPL:
    assert src.count(old) == 1, ("hook text not found exactly once", old[:80])
    src = src.replace(old, new)
src = src.replace('"""r50_chem_lbl_round = r48_chem_q7c_round.py + ONE option, --wires',
                  '"""chem_wr7_round (WR7) = r50_chem_lbl_round.py with its three wire hooks generalised to wires on bits 0 AND 1\n'
                  '(chem_wr7_chainb.DynChainB: the nearest word whose bits keep the group matchable; made by chem_wr7_mkround.py).\n'
                  'r50_chem_lbl_round = r48_chem_q7c_round.py + ONE option, --wires', 1)
open(os.path.join(T, "chem_wr7_round.py"), "w").write(src)
print("written tools/chem_wr7_round.py")
