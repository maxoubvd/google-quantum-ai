#!/usr/bin/env python
"""chem_lens: the price of a code - progressive 8-bit rounding of a continuous DF-THC table under Reed-Muller
codes on chosen bits, followed by the loader compiler's price of the resulting words.

One run gives one point of the exchange curve: error_2 of the all-frozen coded table (mHa) and lookup Toffolis per
pass.  Compare runs made with the same schedule; the difference of error_2^2 against the uncoded run is the
code's cost in mHa^2 and the difference of the compiler price is what it buys.

Method = the public pipeline's coded rounding (marculera dfthc_refit experiments/quantize_rm.py --open, Apache-2.0;
its rmcode.py is imported unchanged): every stage takes, on each coded column and 16-address block, the codeword
of least weighted squared offset that agrees with what is already frozen, freezes the open angles of largest
impact, and refits the free angles and all weights.  Differences from the pipeline: float32 L-BFGS on the
float64-difference objective (chem_factor_fit's dlbfgs) for every stage, and impact = weight norm^2 x mass entering
the rotation x squared offset instead of sampled Gauss-Newton curvature.

  python chem_lens_codes.py INIT.npz OUT_DIR --codes "0:1:3-14" [--deal] [--iters 3000] [--price]
  codes: 'bit:degree:first-last' over five-word batches, comma separated ('' = uncoded).
  --xcodes: the same syntax for codes on the XOR OF THE FIVE WORDS of a batch at that bit (the codes of the public
      C44 build: "1:2:3-7,1:1:8-14,2:2:8-14" on top of --codes "0:1:3-14"); a violated relation is repaired on the
      cheapest of the five words.
  --pretrain 1e-4,1e-3,1e-2 [--pre-iters 4000]: "codes trained in" - before any angle is frozen, a continuous descent
      with the penalty mu * sum over coded entries of (theta - nearest value on that entry's codeword)^2 (in steps),
      the codewords re-chosen every 1,000 iterations (alternating projection), so the fit itself moves onto the code.
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

sys.path.insert(0, os.path.join(L.PIPE, "experiments"))
import rmcode  # noqa: E402  (public pipeline, read-only)
import onebody  # noqa: E402

N = 76
FULL = tuple(round(1.0 - f, 4) for f in (0.9995, 0.999, 0.9975, 0.995, 0.99, 0.98, 0.965, 0.95, 0.925, 0.9, 0.86, 0.82, 0.77,
                                         0.72, 0.66, 0.6, 0.54, 0.48, 0.42, 0.36, 0.3, 0.25, 0.2, 0.16, 0.12, 0.09, 0.065, 0.045,
                                         0.03, 0.02, 0.012, 0.006, 0.0))
PUB = os.path.expanduser("~/qcc/work/pub/dfthcsos/submission.npz")


def deal(strength: np.ndarray, block: int = 16) -> np.ndarray:
    """deal_rows.deal of the public pipeline: rank k goes to block k % nb."""
    B = len(strength)
    nb = -(-B // block)
    ranks = np.argsort(-strength)
    perm = np.empty(B, dtype=int)
    sizes = [min(block, B - j * block) for j in range(nb)]
    slots = [0] * nb
    j = 0
    for row in ranks:
        while slots[j] >= sizes[j]:
            j = (j + 1) % nb
        perm[j * block + slots[j]] = row
        slots[j] += 1
        j = (j + 1) % nb
    return perm


def parse_x(spec: str, nbatch: int = 15, bits: int = 8) -> np.ndarray:
    xdeg = np.zeros((nbatch, bits), np.int64)
    for item in filter(None, (spec or "").split(",")):
        bit, d, span = item.split(":")
        lo, hi = (int(x) for x in span.split("-")) if "-" in span else (int(span),) * 2
        xdeg[lo:hi + 1, int(bit)] = int(d)
    return xdeg


def check_vec(p: int, d: int) -> int:
    """Column of the parity-check matrix of RM(d, 4) at position p, as a bit mask (monomials of degree <= 3 - d)."""
    b = [(p >> i) & 1 for i in range(4)]
    mons = [1]
    if d <= 2:
        mons += b
    if d <= 1:
        mons += [b[i] & b[j] for i in range(4) for j in range(i + 1, 4)]
    return sum(v << i for i, v in enumerate(mons))


def xor_slack(c, blk, pos, xdeg):
    """Entries kept for last (as rmcode.slack does for per-column codes): for every block and XOR-coded batch, the
    weakest word of each of the weakest rows whose parity checks are independent (5 rows for a degree-2 code, 11
    for an affine one).  They cost nothing in the codeword choice, so every relation is met on an entry that
    hardly matters; without them a repaired entry has a large offset, freezes at once, and repairs pile up."""
    res = np.zeros(c.shape, bool)
    for k in np.flatnonzero(xdeg.any(axis=1)):
        d = int(xdeg[k][xdeg[k] > 0].min())
        need = {1: 11, 2: 5, 3: 1}[d]
        cols = np.arange(5 * k, 5 * k + 5)
        cmin, jmin = c[:, cols].min(axis=1), c[:, cols].argmin(axis=1)
        for b in np.unique(blk):
            idx = np.flatnonzero(blk == b)
            if len(idx) <= 4:
                continue
            piv: dict[int, int] = {}
            for i in idx[np.argsort(cmin[idx])]:
                v = check_vec(int(pos[i]), d)
                while v and v.bit_length() in piv:
                    v ^= piv[v.bit_length()]
                if v:
                    piv[v.bit_length()] = v
                    res[i, cols[jmin[i]]] = True
                    if len(piv) == need:
                        break
    return res


def unified_residues(t, c, blk, pos, deg, xdeg, frozen, value, scale=None, reserved=None):
    """Residue (mod a power of two) of every entry: per-column Reed-Muller codes (deg, as rmcode.residues) and codes
    on the XOR of the five words of a batch (xdeg), bit by bit from bit 0, cheapest codewords that agree with the
    frozen entries.  t in steps (rows, 75), c the weights.  Returns r (rows, 75) and the modulus of every column."""
    rows, ncols = t.shape
    nblocks = int(blk.max()) + 1
    r = np.zeros((rows, ncols), np.int64)
    nbits = np.zeros(ncols, np.int64)
    gauge = (lambda v, m_: rmcode.into_gauge(v, m_, scale)) if scale else (lambda v, m_: v)
    for bit in range(deg.shape[1]):
        m = 1 << bit
        pc = np.flatnonzero(deg[:, bit])
        xb = np.flatnonzero(xdeg[:, bit])
        act = np.concatenate([pc] + [np.arange(5 * k, 5 * k + 5) for k in xb]).astype(np.int64)
        if not len(act):
            continue
        assert np.all(nbits[act] == bit), "coded bits of a column must be contiguous from bit 0"
        cost = [c[:, act] * (gauge(rmcode.residue_round(t[:, act], r[:, act] + v * m, 2 * m), 2 * m) - t[:, act]) ** 2 for v in (0, 1)]
        if reserved is not None and len(xb):
            rs = reserved[:, act]
            rs[:, :len(pc)] = False                     # only the words of XOR-coded batches absorb relations
            cost = [np.where(rs, 1e-9 * x, x) for x in cost]
        have = (value[:, act] >> bit) & 1
        big = 1e4 * max(float(x.max()) for x in cost) + 1e-300
        cost = [np.where(frozen[:, act], np.where(have == v, 0.0, big), cost[v]) for v in (0, 1)]
        table = [np.zeros((nblocks, 16, len(act))) for _ in (0, 1)]
        for v in (0, 1):
            table[v][blk, pos] = cost[v]
        rbit = np.zeros((rows, len(act)), np.int64)
        for d in (np.unique(deg[pc, bit]) if len(pc) else []):
            sidx = np.flatnonzero(deg[act[:len(pc)], bit] == d)
            words = rmcode.code(int(d))
            gain = (table[1] - table[0])[:, :, sidx].transpose(0, 2, 1) @ words.T
            pick = words[gain.argmin(axis=2)]                                   # (blocks, columns, 16)
            rbit[:, sidx] = pick[blk, :, pos].astype(np.int64)
        for q_, k in enumerate(xb):
            sidx = len(pc) + 5 * q_ + np.arange(5)
            t0, t1 = table[0][:, :, sidx], table[1][:, :, sidx]                 # (blocks, 16, 5)
            beta = (t1 < t0).astype(np.int64)
            delta = np.abs(t1 - t0)
            x = beta.sum(-1) % 2
            jmin, dmin = delta.argmin(-1), delta.min(-1)
            words = rmcode.code(int(xdeg[k, bit]))
            costcw = ((words[None, :, :] != x[:, None, :]) * dmin[:, None, :]).sum(-1)
            cw = words[costcw.argmin(1)]
            flip = (cw != x).astype(np.int64)
            bi, pi = np.nonzero(flip)
            beta[bi, pi, jmin[bi, pi]] ^= 1
            rbit[:, sidx] = beta[blk, pos]
        r[:, act] += rbit * m
        nbits[act] += 1
    return r, 1 << nbits


def check_x(words, addr, xdeg):
    blk, pos = addr // 16, addr % 16
    for k, bit in zip(*np.nonzero(xdeg)):
        cw = rmcode.code(int(xdeg[k, bit]))
        x = np.bitwise_xor.reduce((words[:, 5 * k:5 * k + 5] >> bit) & 1, axis=1)
        for b in np.unique(blk):
            rows = np.flatnonzero(blk == b)
            assert np.any(np.all(cw[:, pos[rows]] == x[rows], axis=1)), f"XOR code: bit {bit} of batch {k}, block {b}"


def full_table(words2: np.ndarray, edges: np.ndarray, deg=None, xdeg=None) -> np.ndarray:
    """76 one-body rows + the two-body words.  The one-body rows are the public table's vectors re-encoded on `edges`
    and put on the SAME codes as the two-body rows (the compiler only uses a code that holds on every block,
    including blocks 0-4, which hold the one-body rows and the first four two-body rows).  Their error_1 is not
    optimised here: this table is for pricing the loaders."""
    with np.load(PUB) as z:
        pw, pe = np.array(z["angle_words"]).astype(np.int64), np.array(z["tree_edges"]).astype(np.int64)
    v1 = L.decode_np(pw[:N] / 128.0, pe)
    t1 = L.signed_angles(v1, edges) * 128.0                                  # (76, 75), in steps
    if deg is None or not (np.any(deg) or np.any(xdeg)):
        w1 = np.round(t1).astype(np.int64)
    else:
        w2s = ((words2[:4].astype(np.int64) + 128) % 256) - 128              # rows 76..79 share block 4
        t = np.concatenate([t1, w2s.astype(np.float64)])
        addr = np.arange(N + 4)
        frozen = np.zeros(t.shape, bool)
        frozen[N:] = True
        value = np.concatenate([np.round(t1).astype(np.int64), w2s])
        gauged = bool(np.any(deg[:, 1:]) or np.any(xdeg[:, 1:]))
        r_, m_ = unified_residues(t, np.ones_like(t), addr // 16, addr % 16, deg, xdeg, frozen, value, 128.0 if gauged else None)
        v = rmcode.residue_round(t, r_, m_)
        if gauged:
            v = rmcode.into_gauge(v, m_, 128.0)
        w1 = v[:N].astype(np.int64)
        gz = gauged
    else_gauged = bool(deg is not None and (np.any(deg[:, 1:]) or np.any(xdeg[:, 1:])))
    w1 = np.mod(w1, 256)
    if not else_gauged:
        w1 = onebody.signed_gauge(w1, edges, 8)          # bit 0 is unchanged by the gauge moves, bit 7 = bit 6 afterwards
    table = np.concatenate([w1, words2.astype(np.int64)], axis=0)
    full_addr = np.arange(table.shape[0])
    if deg is not None:
        rmcode.check(table, full_addr, deg)
        check_x(table, full_addr, xdeg)
    return table


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("init")
    ap.add_argument("out")
    ap.add_argument("--codes", default="")
    ap.add_argument("--xcodes", default="")
    ap.add_argument("--deal", action="store_true")
    ap.add_argument("--iters", type=int, default=3000)
    ap.add_argument("--final", type=int, default=6000)
    ap.add_argument("--kappa", type=float, default=5e-5)
    ap.add_argument("--chunk", type=int, default=200)
    ap.add_argument("--recenter", type=int, default=1000)
    ap.add_argument("--price", action="store_true")
    ap.add_argument("--force-gauge", action="store_true", help="control: keep every angle in the signed range (as the runs with coded bits above bit 0 must) without any such code")
    ap.add_argument("--gauge-from", type=int, default=0, help="with a signed range in force, columns below this index stay free (control)")
    ap.add_argument("--max-stage", type=int, default=99, help="stop after this stage (a probe; nothing is written)")
    ap.add_argument("--pretrain", default="")
    ap.add_argument("--pre-iters", type=int, default=4000)
    ap.add_argument("--pen-stage", type=float, default=0.0,
                    help="keep the pretraining penalty with this mu on the not yet frozen entries during the rounding stages "
                         "(without it the refits of the first stages undo the pretraining)")
    ap.add_argument("--pre-all", action="store_true", help="pretraining penalty on every entry (soft rounding), not only on coded columns")
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
            perm = deal(np.sum(w[r] ** 2, axis=1))
            U[r], w[r] = U[r][perm], w[r][perm]
    edges = L.log_depth_tree(N)
    scale = 128.0
    deg = rmcode.parse(a.codes, N - 1, 8)
    xdeg = parse_x(a.xcodes)
    xcols = np.repeat(xdeg.any(axis=1), 5)
    ccols = np.flatnonzero(deg.any(axis=1) | xcols)
    gauged = bool(deg[:, 1:].any() or xdeg[:, 1:].any() or a.force_gauge)
    early = 0.5 / scale if gauged else 0.0
    g_lo, g_hi = -0.5 - early, 0.5 - early - 1e-9
    theta, mass = angles_and_mass(U.reshape(-1, N), edges)
    if gauged:
        theta = np.clip(rmcode.gauge_angles(theta, edges, early), g_lo, g_hi)
    addr = N + np.arange(R * B)
    blk, pos = addr // 16, addr % 16
    log(f"init {a.init} {(R, B, C)} codes [{a.codes}] ({int((deg > 0).sum())} coded column-bits) xcodes [{a.xcodes}] "
        f"({int((xdeg > 0).sum())} XOR-of-five relations per block) deal {a.deal} gauged {gauged}")

    jax, jnp = L.jax_setup(x64=False)
    import optax
    refc = L.Reference()
    decode = L.make_decode(edges)
    i_, j_ = jnp.asarray(refc.i), jnp.asarray(refc.j)
    sM32 = jnp.asarray(np.sqrt(refc.M), jnp.float32)
    sM64 = np.sqrt(refc.M)
    sw = float(np.sqrt(np.mean(w ** 2)))
    nbw = np.linalg.norm(w, axis=2).reshape(-1)
    nb = np.maximum(nbw / np.sqrt(np.mean(nbw ** 2)), 3e-2)
    srow = 1.0 / nb
    srow_j = jnp.asarray(srow, jnp.float32)[:, None]

    def A_np(th, w_):
        u = L.decode_np(th, edges).reshape(R, B, N)
        return np.einsum("rbc,rbp->rcp", w_, u[:, :, refc.i] * u[:, :, refc.j]).reshape(R * C, -1)

    def reference(th, w_):
        A0 = A_np(th, w_)
        r0 = sM64 * (A0.T @ A0 - refc.T)
        return A0, r0, float(np.sqrt(np.sum(r0 * r0)))

    opt = optax.lbfgs(memory_size=30)
    gcol_np = (np.arange(N - 1) >= a.gauge_from)[None, :]
    GCOL = jnp.asarray(gcol_np)

    def clip_np(own):
        return np.where(gcol_np, np.clip(own, g_lo, g_hi), own) if gauged else own

    @jax.jit
    def run(params, st, A0, r0, eps0, off, fmask, fq, mu=0.0, tgt=0.0, cmask=0.0):
        def f(p):
            own = p["tau"] * srow_j
            if gauged:
                own = jnp.where(GCOL, jnp.clip(own, g_lo, g_hi), own)
            th = fmask * fq + (1.0 - fmask) * own
            u = decode(th).reshape(R, B, N)
            A = jnp.einsum("rbc,rbp->rcp", p["v"] * sw, u[:, :, i_] * u[:, :, j_]).reshape(R * C, -1)
            dA = A - A0
            S = A0.T @ dA
            dr = sM32 * (S + S.T + dA.T @ dA)
            delta = jnp.sum(dr * (2.0 * r0 + dr))
            deps = delta / (eps0 + jnp.sqrt(jnp.maximum(eps0 * eps0 + delta, 0.0)))
            pen = mu * jnp.sum(cmask * ((th - tgt) * scale) ** 2)
            return off + deps + a.kappa * L.d_penalty(p["v"] * sw) + pen
        vg = optax.value_and_grad_from_state(f)

        def body(_, carry):
            p, s = carry
            value, grad = vg(p, state=s)
            upd, s = opt.update(grad, s, p, value=value, grad=grad, value_fn=f)
            return optax.apply_updates(p, upd), s
        return jax.lax.fori_loop(0, a.chunk, body, (params, st))

    def coded_round(now, c):
        q = np.round(now * scale) / scale
        if len(ccols) == 0:
            return q
        t = now * scale
        r_, m_ = unified_residues(t, c, blk, pos, deg, xdeg, frozen, np.round(qv * scale).astype(np.int64),
                                  scale if gauged else None, RES[0])
        v = rmcode.residue_round(t, r_, m_)
        q[:, ccols] = ((rmcode.into_gauge(v, m_, scale) if gauged else v) / scale)[:, ccols]
        return q

    RES = [None]
    t0 = time.time()
    frozen = np.zeros(theta.shape, bool)
    qv = np.zeros_like(theta)
    e0 = reference(theta, w)[2]
    log(f"start {1e3 * e0:.4f} mHa")
    params = {"tau": jnp.asarray(theta / srow[:, None], jnp.float32), "v": jnp.asarray(w / sw, jnp.float32)}
    trace = []
    if a.pretrain and (len(ccols) or a.pre_all):
        cmask = np.ones(theta.shape) if a.pre_all else np.zeros(theta.shape)
        cmask[:, ccols] = 1.0
        zero = jnp.zeros(theta.shape, jnp.float32)
        for mu in [float(x) for x in a.pretrain.split(",") if x]:
            st = opt.init(params)
            done, e_first = 0, None
            while done < a.pre_iters:
                own = np.asarray(params["tau"], np.float64) * srow[:, None]
                theta = clip_np(own)
                w_ = np.asarray(params["v"], np.float64) * sw
                A0, r0, eps0 = reference(theta, w_)
                if e_first is None:
                    e_first = eps0
                mass = angles_and_mass(L.decode_np(theta, edges), edges)[1]
                tgt = coded_round(theta, (nbw[:, None] ** 2) * mass + 1e-30)          # nearest point of the coded lattice
                ref = (jnp.asarray(A0, jnp.float32), jnp.asarray(r0, jnp.float32), jnp.float32(eps0), jnp.float32(eps0 - e_first))
                for _ in range(max(1, a.recenter // a.chunk)):
                    params, st = run(params, st, *ref, zero, zero, jnp.float32(mu), jnp.asarray(tgt, jnp.float32), jnp.asarray(cmask, jnp.float32))
                    done += a.chunk
            own = np.asarray(params["tau"], np.float64) * srow[:, None]
            theta = clip_np(own)
            w = np.asarray(params["v"], np.float64) * sw
            mass = angles_and_mass(L.decode_np(theta, edges), edges)[1]
            tgt = coded_round(theta, (nbw[:, None] ** 2) * mass + 1e-30)
            dev = np.sqrt(np.sum(cmask * ((theta - tgt) * scale) ** 2) / cmask.sum())
            log(f"pretrain mu {mu:.0e}: error_2 {1e3 * e_first:8.3f} -> {1e3 * reference(theta, w)[2]:8.3f} mHa; rms distance of the coded "
                f"entries to their codeword {dev:.3f} steps  ({time.time() - t0:.0f} s)")
    for k, frac in enumerate(FULL):
        need = int(round(frac * theta.size)) - int(frozen.sum())
        if need > 0:
            now = np.where(frozen, qv, theta)
            mass = angles_and_mass(L.decode_np(now, edges), edges)[1]
            c = (nbw[:, None] ** 2) * mass + 1e-30
            if RES[0] is None and xdeg.any():
                RES[0] = xor_slack(c, blk, pos, xdeg)
                log(f"slack entries for the XOR relations: {int(RES[0].sum())} ({RES[0].mean():.4f} of all entries)")
            q = coded_round(now, c)
            impact = c * (q - now) ** 2
            if RES[0] is not None:
                impact = rmcode.last(impact, RES[0])
            impact[frozen] = -np.inf
            pick = np.argsort(-impact, axis=None)[:need]
            frozen.flat[pick] = True
            qv.flat[pick] = q.flat[pick]
        fm, fq = jnp.asarray(frozen, jnp.float32), jnp.asarray(qv, jnp.float32)
        st = opt.init(params)
        iters = a.final if frac >= 1.0 else a.iters
        done, e_first = 0, None
        while done < iters:
            own = np.asarray(params["tau"], np.float64) * srow[:, None]
            th = np.where(frozen, qv, clip_np(own))
            w_ = np.asarray(params["v"], np.float64) * sw
            A0, r0, eps0 = reference(th, w_)
            if e_first is None:
                e_first = eps0
            ref = (jnp.asarray(A0, jnp.float32), jnp.asarray(r0, jnp.float32), jnp.float32(eps0), jnp.float32(eps0 - e_first))
            if a.pen_stage > 0 and not frozen.all():
                mass_s = angles_and_mass(L.decode_np(th, edges), edges)[1]
                tgt_s = coded_round(th, (nbw[:, None] ** 2) * mass_s + 1e-30)
                cm_s = (np.ones(theta.shape) if a.pre_all else np.isin(np.arange(N - 1), ccols)[None, :] * np.ones(theta.shape)) * (~frozen)
                pen_args = (jnp.float32(a.pen_stage), jnp.asarray(tgt_s, jnp.float32), jnp.asarray(cm_s, jnp.float32))
            else:
                pen_args = (jnp.float32(0.0), fq, fm)
            for _ in range(max(1, a.recenter // a.chunk)):
                params, st = run(params, st, *ref, fm, fq, *pen_args)
                done += a.chunk
        own = np.asarray(params["tau"], np.float64) * srow[:, None]
        theta = np.where(frozen, qv, clip_np(own))
        w = np.asarray(params["v"], np.float64) * sw
        eps = reference(theta, w)[2]
        trace.append((frac, 1e3 * e_first, 1e3 * eps))
        log(f"stage {k:2d}: frozen {frac:.4f}  error_2 {1e3 * e_first:8.3f} -> {1e3 * eps:8.3f} mHa  ({time.time() - t0:.0f} s)")
        if k >= a.max_stage:
            log("stopped by --max-stage (probe)")
            return
    raw = np.mod(np.round(qv * scale).astype(np.int64), 256)
    words = onebody.signed_gauge(raw, edges, 8)
    if gauged:
        assert np.array_equal(words, raw), "the gauge fix moved a word"
    rmcode.check(words, addr, deg)
    check_x(words, addr, xdeg)
    u = L.decode_np(words / scale, edges).reshape(R, B, N)
    e_fin = refc.error2_np(u, w)
    out = {"codes": a.codes, "xcodes": a.xcodes, "deal": a.deal, "error_2_mHa": 1e3 * e_fin, "iters": a.iters}
    np.savez(os.path.join(a.out, "rounded.npz"), words=words, W=w, codes=deg)
    log(f"all frozen, codes verified on every block: error_2 {1e3 * e_fin:.4f} mHa (continuous {1e3 * e0:.4f})  ({time.time() - t0:.0f} s)")
    if a.price:
        import chem_loaders_compile as clc
        table = full_table(words, edges, deg, xdeg)
        np.save(os.path.join(a.out, "table.npy"), table)
        res = clc.price_table(table, "auto", a.effort, 1, 0)
        out["per_pass"] = int(res["per_pass"])
        out["per_circuit"] = [int(r["toffoli"]) for r in res["circuits"]]
        log(f"compiler price (effort {a.effort}): {res['per_pass']} per pass; per circuit {out['per_circuit']}")
    json.dump(out, open(os.path.join(a.out, "result.json"), "w"))


if __name__ == "__main__":
    main()
