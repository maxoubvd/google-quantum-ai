#!/usr/bin/env python
"""chem_lens: progressive 8-bit rounding that keeps every word in the signed range WITHOUT a wall.

Problem (research/ws_chem_lens.md 3g): a table with a code on a bit above bit 0 must have its words in the signed
range (cos >= 0 at every rotation) when the rounding ends, because the gauge fix that normally runs at the end
negates words and would break the code.  Holding every angle inside the range with a clip costs 1.1 mHa at 14%
frozen (controls G, G15): a free angle that wants to cross +-pi/2 cannot.

Fix tried here: let it cross.  Crossing at edge t = (p -> c) is the exact gauge move
        theta_t -> theta_t -+ pi,   theta -> -theta on the edge into p and on the earlier sibling edges out of p,
which leaves the row's vector unchanged up to a global sign.  The half turn keeps the low 7 bits of word t; a
negation keeps bit 0 and changes the bits above the lowest set bit.  So the move is harmless when every negated
entry is "negatable":
  * not frozen yet (its code residues are re-chosen at every stage anyway), or
  * frozen in a column with no code above bit 0 (the stored word becomes -w; -64 is excluded), or
  * frozen in a coded column and, rule "slack": the slack entries of its block and batch are still free, so
    the changed XOR-of-five relation is re-absorbed by them (chem_lens_codes.xor_slack), or rule "strict": w == 0.
Free angles whose crossing is allowed run without any clip; at the end of every stage each row is re-gauged exactly
(vectorised over rows, edges in reverse order); the others keep the clip.  Inside a stage the set of allowed
crossings cannot change (negation keeps w == 0 and w != -64, and no slack entry is frozen inside a stage), so no
re-gauge is needed there; the first version re-gauged every 1,000 iterations and restarted L-BFGS each time, which
alone cost 0.045 mHa at stage 0 (--regauge-at recenter reproduces it).

  python chem_lens_gauge.py INIT.npz OUT --codes 0:1:3-14 --xcodes 1:2:3-7,1:1:8-14,2:2:8-14 --deal [--price]
  --no-relations     keep the structure the codes define (coded columns, slack entries, their late freezing) but
                     round to nearest: the control "signed range enforced, no code" to compare with the free run
  --hcols K          control without code specs: columns >= K count as coded above bit 0
  --neg-rule strict|slack     (default slack)      --no-regauge  the old clip everywhere
  --max-minutes M    stop after the stage that passes M minutes; rerun the same command to resume (OUT/state.npz)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_factor_lib as L  # noqa: E402
from chem_lens_orbit import angles_and_mass  # noqa: E402
import chem_lens_codes as cc  # noqa: E402

rmcode, onebody = cc.rmcode, cc.onebody
N = 76
LONG = tuple(round(1.0 - f, 5) for f in (0.9995, 0.999, 0.998, 0.997, 0.9955, 0.994, 0.992, 0.99, 0.9875, 0.985, 0.98, 0.975, 0.97,
                                         0.965, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.9, 0.885, 0.87, 0.855, 0.84, 0.82, 0.8, 0.78,
                                         0.76, 0.74, 0.72, 0.69, 0.66, 0.63, 0.6, 0.57, 0.54, 0.51, 0.48, 0.45, 0.42, 0.39, 0.36,
                                         0.33, 0.3, 0.27, 0.24, 0.21, 0.18, 0.15, 0.125, 0.1, 0.08, 0.06, 0.045, 0.03, 0.02, 0.012,
                                         0.006, 0.0))


def neighbour_matrix(edges: np.ndarray) -> np.ndarray:
    """NB[t, k] = 1 when a crossing at edge t negates edge k: the edge into t's parent and t's earlier siblings."""
    n = len(edges)
    NB = np.zeros((n, n), np.int32)
    for t, (p, _) in enumerate(edges):
        for k, (pp, qq) in enumerate(edges):
            if (qq == p) or (pp == p and k < t):
                NB[t, k] = 1
    return NB


def negatable(frozen, wv, hcol, NB, protected=None, cascade=True):
    """Which entries may change sign, and which crossings are therefore allowed.
    frozen, wv: (rows, 75) frozen mask and frozen words (signed, -64..63); hcol: columns with a code above bit 0;
    protected: frozen entries of coded columns whose word must not change (rule "slack": the slack entries of
    their block and batch are frozen; rule "strict": every non-zero word).  A frozen word -64 is special: its
    negation is +64 = -64 plus a half turn, so the stored word stays and the sign change passes on to the
    entry's own neighbours; with `cascade` it is negatable exactly when that crossing is allowed (computed as
    the least fixed point; the neighbours of an edge have smaller indices, so it ends)."""
    edge = frozen & (wv == -64)
    negat = ~edge
    if protected is not None:
        negat &= ~(protected & hcol[None, :] & ~edge)
    if cascade:
        for _ in range(NB.shape[0] + 1):
            flip_ok = ((~negat).astype(np.int32) @ NB.T) == 0
            new = np.where(edge, flip_ok, negat)
            if np.array_equal(new, negat):
                break
            negat = new
    flip_ok = ((~negat).astype(np.int32) @ NB.T) == 0
    return negat, flip_ok


def regauge(th, frozen, flip_ok, negat, edges, lo, hi):
    """Exact re-gauging of every row into the signed range [lo, hi), in units of pi.  Free entries whose crossing is
    allowed get +-1 (a half turn) and the entries it negates change sign; free entries whose crossing is not allowed
    are clipped; frozen entries only ever change sign, and only where `negat` says so.  Returns the new angles and
    the counts (crossings, negated frozen entries, clipped entries)."""
    out = np.where(frozen, th, np.mod(th + 1.0, 2.0) - 1.0)
    free = ~frozen
    sign = np.ones((th.shape[0], edges.shape[0] + 1), np.int8)
    nflip = nneg = nclip = 0
    for k in range(edges.shape[0] - 1, -1, -1):
        p, q = edges[k]
        neg = sign[:, p] != sign[:, q]
        assert not np.any(neg & ~negat[:, k]), "a protected entry would be negated"
        a = np.where(neg, -out[:, k], out[:, k])
        nneg += int(np.sum(neg & frozen[:, k]))
        oor = (a >= hi) | (a < lo)
        # a frozen entry leaves the range only as -64 -> +64; its half turn gives the stored word -64 back and the
        # sign change passes on to its own neighbours (allowed only where `negatable` found the whole cascade allowed)
        fl = oor & flip_ok[:, k]
        a = np.where(fl, np.where(a >= hi, a - 1.0, a + 1.0), a)
        cl = oor & free[:, k] & ~flip_ok[:, k]
        a = np.where(cl, np.clip(a, lo, hi - 1e-9), a)
        assert not np.any(oor & frozen[:, k] & ~flip_ok[:, k]), "a frozen entry left the signed range"
        out[:, k] = a
        nflip += int(fl.sum())
        nclip += int(cl.sum())
        sign[fl, p] *= -1
        sign[:, q] = 1
    return out, nflip, nneg, nclip


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("init")
    ap.add_argument("out")
    ap.add_argument("--codes", default="")
    ap.add_argument("--xcodes", default="")
    ap.add_argument("--deal", action="store_true")
    ap.add_argument("--no-relations", action="store_true")
    ap.add_argument("--hcols", type=int, default=-1)
    ap.add_argument("--neg-rule", default="slack", choices=["strict", "slack"])
    ap.add_argument("--no-regauge", action="store_true")
    ap.add_argument("--regauge-at", default="stage", choices=["stage", "recenter"],
                    help="stage: re-gauge once at the end of every stage (the wall mask cannot change inside a stage, "
                         "so nothing else is needed and the L-BFGS memory survives); recenter: the first version, which "
                         "re-gauged and restarted L-BFGS every --recenter iterations (measured worse, log section 6)")
    ap.add_argument("--no-cascade", action="store_true", help="a frozen word -64 is never negated (the first version)")
    ap.add_argument("--absorb", default="last", choices=["last", "complete"],
                    help="last: reserved slack entries absorb the XOR relations and freeze at the very end (nothing but "
                         "other slack entries is free to compensate); complete: no reserved entries, the last free word "
                         "of a (row, batch) unit is frozen as soon as its four partners are")
    ap.add_argument("--schedule", default="full", choices=["full", "long"])
    ap.add_argument("--iters", type=int, default=3000)
    ap.add_argument("--final", type=int, default=6000)
    ap.add_argument("--kappa", type=float, default=5e-5)
    ap.add_argument("--chunk", type=int, default=200)
    ap.add_argument("--recenter", type=int, default=1000)
    ap.add_argument("--max-minutes", type=float, default=30.0)
    ap.add_argument("--price", action="store_true")
    ap.add_argument("--effort", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    logf = open(os.path.join(a.out, "log.txt"), "a")

    def log(*x):
        s = " ".join(str(v) for v in x)
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    z = np.load(a.init)
    U, w = np.array(z["U"], np.float64), np.array(z["W"], np.float64)
    R, B, C = w.shape
    U = U / np.linalg.norm(U, axis=2, keepdims=True)
    if a.deal:
        for r in range(R):
            perm = cc.deal(np.sum(w[r] ** 2, axis=1))
            U[r], w[r] = U[r][perm], w[r][perm]
    edges = L.log_depth_tree(N)
    NB = neighbour_matrix(edges)
    scale = 128.0
    deg = rmcode.parse(a.codes, N - 1, 8)
    xdeg = cc.parse_x(a.xcodes)
    hcol = deg[:, 1:].any(axis=1) | np.repeat(xdeg[:, 1:].any(axis=1), 5)
    if a.hcols >= 0:
        hcol = hcol | (np.arange(N - 1) >= a.hcols)
    ccols = np.flatnonzero(deg.any(axis=1) | np.repeat(xdeg.any(axis=1), 5))
    early = 0.5 / scale
    lo, hi = -0.5 - early, 0.5 - early
    addr = N + np.arange(R * B)
    blk, pos = addr // 16, addr % 16
    nblocks = int(blk.max()) + 1
    sched = cc.FULL if a.schedule == "full" else LONG
    state_path = os.path.join(a.out, "state.npz")

    jax, jnp = L.jax_setup(x64=False)
    import optax
    refc = L.Reference()
    decode = L.make_decode(edges)
    i_, j_ = jnp.asarray(refc.i), jnp.asarray(refc.j)
    sM32 = jnp.asarray(np.sqrt(refc.M), jnp.float32)
    sM64 = np.sqrt(refc.M)
    sw = float(np.sqrt(np.mean(w ** 2)))
    nbw = np.linalg.norm(w, axis=2).reshape(-1)
    srow = 1.0 / np.maximum(nbw / np.sqrt(np.mean(nbw ** 2)), 3e-2)
    srow_j = jnp.asarray(srow, jnp.float32)[:, None]

    def reference(th, w_):
        u = L.decode_np(th, edges).reshape(R, B, N)
        A0 = np.einsum("rbc,rbp->rcp", w_, u[:, :, refc.i] * u[:, :, refc.j]).reshape(R * C, -1)
        r0 = sM64 * (A0.T @ A0 - refc.T)
        return A0, r0, float(np.sqrt(np.sum(r0 * r0)))

    opt = optax.lbfgs(memory_size=30)

    @jax.jit
    def run(params, st, A0, r0, eps0, off, fmask, fq, wall):
        def f(p):
            own = p["tau"] * srow_j
            own = jnp.where(wall > 0, jnp.clip(own, lo, hi - 1e-7), own)
            th = fmask * fq + (1.0 - fmask) * own
            u = decode(th).reshape(R, B, N)
            A = jnp.einsum("rbc,rbp->rcp", p["v"] * sw, u[:, :, i_] * u[:, :, j_]).reshape(R * C, -1)
            dA = A - A0
            S = A0.T @ dA
            dr = sM32 * (S + S.T + dA.T @ dA)
            delta = jnp.sum(dr * (2.0 * r0 + dr))
            deps = delta / (eps0 + jnp.sqrt(jnp.maximum(eps0 * eps0 + delta, 0.0)))
            return off + deps + a.kappa * L.d_penalty(p["v"] * sw)
        vg = optax.value_and_grad_from_state(f)

        def body(_, carry):
            p, s = carry
            value, grad = vg(p, state=s)
            upd, s = opt.update(grad, s, p, value=value, grad=grad, value_fn=f)
            return optax.apply_updates(p, upd), s
        return jax.lax.fori_loop(0, a.chunk, body, (params, st))

    # ------------------------------------------------------------------ state
    if os.path.exists(state_path):
        with np.load(state_path) as s:
            theta, w, frozen, qv, stage0, spent = s["theta"], s["w"], s["frozen"], s["qv"], int(s["stage"]), float(s["spent"])
            RES = s["res"] if "res" in s.files else None
        log(f"resume at stage {stage0} (frozen {frozen.mean():.4f}, {spent / 60:.1f} minutes spent so far)")
    else:
        theta, _ = angles_and_mass(U.reshape(-1, N), edges)
        theta = np.clip(rmcode.gauge_angles(theta, edges, early), lo, hi - 1e-9)
        frozen = np.zeros(theta.shape, bool)
        qv = np.zeros_like(theta)
        stage0, spent, RES = 0, 0.0, None
        log(f"init {a.init} {(R, B, C)} codes [{a.codes}] xcodes [{a.xcodes}] relations {not a.no_relations} deal {a.deal} "
            f"coded-above-bit-0 columns {int(hcol.sum())} rule {a.neg_rule} regauge {not a.no_regauge} schedule {a.schedule} x {a.iters}")
        log(f"start {1e3 * reference(theta, w)[2]:.4f} mHa")

    xbat = np.flatnonzero(xdeg[:, 1:].any(axis=1))                         # batches with an XOR-of-five code

    def rules():
        wv = np.round(qv * scale).astype(np.int64)
        if a.absorb == "complete":
            # a row's five words of a batch are a unit; while one of them is free it absorbs the relation, so the
            # frozen ones may still change sign; a complete unit is protected
            done5 = frozen.reshape(frozen.shape[0], 15, 5).all(axis=2)
            protected = frozen & np.repeat(done5, 5, axis=1)
        elif a.neg_rule == "strict" or RES is None:
            protected = frozen & (wv != 0)
        else:
            fb = np.zeros((nblocks, 15), bool)
            rr, cl = np.nonzero(frozen & RES)
            fb[blk[rr], cl // 5] = True
            protected = frozen & fb[blk][:, np.arange(N - 1) // 5]
        return negatable(frozen, wv, hcol, NB, protected, cascade=not a.no_cascade)

    def coded_round(now, c):
        q = np.round(now * scale) / scale
        if a.no_relations or len(ccols) == 0:
            return q
        t = now * scale
        r_, m_ = cc.unified_residues(t, c, blk, pos, deg, xdeg, frozen, np.round(qv * scale).astype(np.int64), scale, RES)
        v = rmcode.into_gauge(rmcode.residue_round(t, r_, m_), m_, scale)
        q[:, ccols] = (v / scale)[:, ccols]
        return q

    t0 = time.time()
    absorbed = 0
    params = {"tau": jnp.asarray(theta / srow[:, None], jnp.float32), "v": jnp.asarray(w / sw, jnp.float32)}
    for k in range(stage0, len(sched)):
        frac = sched[k]
        need = int(round(frac * theta.size)) - int(frozen.sum())
        if need > 0:
            now = np.where(frozen, qv, theta)
            mass = angles_and_mass(L.decode_np(now, edges), edges)[1]
            c = (nbw[:, None] ** 2) * mass + 1e-30
            if RES is None and xdeg[:, 1:].any() and a.absorb == "last":
                RES = cc.xor_slack(c, blk, pos, xdeg)
                log(f"slack entries: {int(RES.sum())} ({RES.mean():.4f} of all entries)")
            q = coded_round(now, c)
            impact = c * (q - now) ** 2
            if RES is not None:
                impact = rmcode.last(impact, RES)
            impact[frozen] = -np.inf
            pick = np.argsort(-impact, axis=None)[:need]
            frozen.flat[pick] = True
            qv.flat[pick] = q.flat[pick]
            nabs = 0
            if a.absorb == "complete" and len(xbat) and not a.no_relations:
                # completion rule: a unit (row, batch) left with ONE free word is completed now, so that its absorber
                # takes the bits the relation needs while the other units are still free to compensate
                free5 = (~frozen).reshape(frozen.shape[0], 15, 5)
                single = np.zeros((frozen.shape[0], 15), bool)
                single[:, xbat] = free5[:, xbat].sum(axis=2) == 1
                mask = (~frozen) & np.repeat(single, 5, axis=1)
                if mask.any():
                    q = coded_round(np.where(frozen, qv, theta), c)            # residues with this stage's freezes in
                    frozen |= mask
                    qv = np.where(mask, q, qv)
                    nabs = int(mask.sum())
                absorbed += nabs
        negat, flip_ok = rules()
        wall = (~frozen) & (np.ones_like(flip_ok) if a.no_regauge else ~flip_ok)
        st = opt.init(params)
        iters = a.final if frac >= 1.0 else a.iters
        done, e_first, flips, negs, clips, unstuck = 0, None, 0, 0, 0, 0
        while done < iters:
            own = np.asarray(params["tau"], np.float64) * srow[:, None]
            th = np.where(frozen, qv, np.where(wall, np.clip(own, lo, hi - 1e-9), own))
            if not a.no_regauge and a.regauge_at == "recenter":
                th, nf, nn, ncl = regauge(th, frozen, flip_ok, negat, edges, lo, hi)
                flips, negs, clips = flips + nf, negs + nn, clips + ncl
                if nf or nn or ncl:
                    qv = np.where(frozen, th, qv)
                    params = {"tau": jnp.asarray(th / srow[:, None], jnp.float32), "v": params["v"]}
                    negat, flip_ok = rules()
                    wall = (~frozen) & ~flip_ok
                    st = opt.init(params)
            else:
                # a walled angle that stepped past the clip has zero gradient there and would stay stuck: put it back
                # on the boundary (the objective does not change; the L-BFGS memory is kept)
                stuck = wall & ((own < lo) | (own > hi - 1e-9))
                if stuck.any():
                    params = {"tau": jnp.where(jnp.asarray(stuck), jnp.asarray(th / srow[:, None], jnp.float32), params["tau"]),
                              "v": params["v"]}
                    unstuck += int(stuck.sum())
            w_ = np.asarray(params["v"], np.float64) * sw
            A0, r0, eps0 = reference(th, w_)
            if e_first is None:
                e_first = eps0
            ref = (jnp.asarray(A0, jnp.float32), jnp.asarray(r0, jnp.float32), jnp.float32(eps0), jnp.float32(eps0 - e_first))
            fm, fq, wl = jnp.asarray(frozen, jnp.float32), jnp.asarray(qv, jnp.float32), jnp.asarray(wall, jnp.float32)
            for _ in range(max(1, a.recenter // a.chunk)):
                params, st = run(params, st, *ref, fm, fq, wl)
                done += a.chunk
        own = np.asarray(params["tau"], np.float64) * srow[:, None]
        theta = np.where(frozen, qv, np.where(wall, np.clip(own, lo, hi - 1e-9), own))
        if not a.no_regauge:
            theta, nf, nn, ncl = regauge(theta, frozen, flip_ok, negat, edges, lo, hi)
            flips, negs, clips = flips + nf, negs + nn, clips + ncl
            qv = np.where(frozen, theta, qv)
            params = {"tau": jnp.asarray(theta / srow[:, None], jnp.float32), "v": params["v"]}
        w = np.asarray(params["v"], np.float64) * sw
        eps = reference(theta, w)[2]
        now_spent = spent + time.time() - t0
        log(f"stage {k:2d}: frozen {frac:.4f}  error_2 {1e3 * e_first:8.3f} -> {1e3 * eps:8.3f} mHa  crossings {flips} "
            f"negated frozen {negs} walled {int(wall.sum())} clipped {clips} put back {unstuck}"
            f"{f' frozen now {frozen.mean():.4f} absorbers so far {absorbed}' if a.absorb == 'complete' else ''}  ({now_spent:.0f} s)")
        extra = {} if RES is None else {"res": RES}
        np.savez(state_path, theta=theta, w=w, frozen=frozen, qv=qv, stage=k + 1, spent=now_spent, **extra)
        if (time.time() - t0) / 60.0 > a.max_minutes and k + 1 < len(sched):
            log(f"stopped after stage {k} by --max-minutes; rerun to resume")
            return
    raw = np.mod(np.round(qv * scale).astype(np.int64), 256)
    words = onebody.signed_gauge(raw, edges, 8)
    moved = int(np.sum(words != raw))
    if not a.no_relations:
        rmcode.check(words, addr, deg)
        cc.check_x(words, addr, xdeg)
    u = L.decode_np(words / scale, edges).reshape(R, B, N)
    e_fin = refc.error2_np(u, w)
    out = {"codes": a.codes, "xcodes": a.xcodes, "relations": not a.no_relations, "rule": a.neg_rule, "regauge": not a.no_regauge,
           "error_2_mHa": 1e3 * e_fin, "gauge_fix_moved_words": moved, "schedule": a.schedule, "iters": a.iters}
    np.savez(os.path.join(a.out, "rounded.npz"), words=words, W=w, U=u, codes=deg, xcodes=xdeg)
    log(f"all frozen: error_2 {1e3 * e_fin:.4f} mHa; words moved by the final gauge fix: {moved} (must be 0); "
        f"{'codes verified on every block; ' if not a.no_relations else ''}({spent + time.time() - t0:.0f} s)")
    if a.price and not a.no_relations:
        import chem_loaders_compile as clc
        table = cc.full_table(words, edges, deg, xdeg)
        np.save(os.path.join(a.out, "table.npy"), table)
        res = clc.price_table(table, "auto", a.effort, 1, 0)
        out["per_pass"] = int(res["per_pass"])
        log(f"compiler price (effort {a.effort}): {res['per_pass']} per pass; per circuit {[int(r['toffoli']) for r in res['circuits']]}")
    json.dump(out, open(os.path.join(a.out, "result.json"), "w"))


if __name__ == "__main__":
    main()
