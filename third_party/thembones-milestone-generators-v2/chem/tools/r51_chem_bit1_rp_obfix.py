"""r51_chem_bit1_rp_obfix (R51-CHEM-BIT1 re-pricer, 2026-11-16): can the one-body rows of a bit-1 layout be made
SHIPPABLE (every chain link affine on blocks 0-4 with rows 76-79 fixed, every word in -32..31) near the designer's
full-circle error_1? Two routes, both at fixed h1', on a finished table:

  R2  'bounded'  CD7's native step (tools/chem_cd7_oblay.progressive_layout_fixed) with ONE change: the L-BFGS-B refits
                 of the free angles are bounded to the signed range (theta / pi in [-1/2, 31/64]), so no free angle
                 drifts out and is clamped at its freeze; then CD7's polish_range. Placements light4 / heavy5, best of N.
  R1  'repair'   the designer's full-circle rounding (tools/r51_chem_bit1_obround.py, tool decoder, heavy5 or light4),
                 then the chain's exact gauge move (obsigned: half turns and negations), then a REPAIR of every
                 (bit, chain, block) whose bits the move broke or that disagrees with rows 76-79: the cheapest common
                 coset of RM(1,4) over that block (the rounding tool's shared_coset_bits, weights lam^2 x mass, a word's
                 bit b changed by +-2^b, lower bits kept, rows 76-79 frozen), then CD7's polish_range (codes and range
                 kept).

    python r51_chem_bit1_rp_obfix.py TABLE.npz --item SPEC_OR_FILE.spec ... [--route both|bounded|repair] [--tries 8]
"""
import argparse
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r51_chem_bit1_round as RD  # noqa: E402
import r51_chem_bit1_obround as DOB  # noqa: E402
import chem_cd7_oblay as OL  # noqa: E402
import chem_cd7_table7 as CT  # noqa: E402
from r48_chem_q7c_onebody import STAGES, angles_and_mass, decode_jac, decode_pi, error1  # noqa: E402
from r51_chem_bit1_codemodel import parse_spec  # noqa: E402
from r51_chem_bit1_rp_obn import link_violations, obsigned_move, repair76, signed, table_H  # noqa: E402

LO, HI = -32, 31


def progressive_bounded(th0, edges, lam_at, H, scale, iters, deg, w76):
    """chem_cd7_oblay.progressive_layout_fixed with the refits bounded to the signed range (the only change)."""
    n, nc = th0.shape
    theta = np.clip(th0.copy(), LO / scale, HI / scale)
    frozen = np.zeros(theta.shape, bool)
    qv = np.zeros(theta.shape)
    val = np.zeros(theta.shape, np.int64)
    lam2 = lam_at[:, None] ** 2
    blk_e = np.concatenate([np.arange(n) // 16, np.full(4, 4)])
    pos_e = np.concatenate([np.arange(n) % 16, np.arange(12, 16)])
    xdeg = np.zeros((15, 8), np.int64)
    w76s = signed(w76)

    def fun(x, mask, base):
        t = base.copy()
        t[mask] = x
        V, J = decode_jac(t, edges)
        D = (V.T * lam_at) @ V - H
        gV = 4.0 * lam_at[:, None] * (V @ D)
        return float(np.sum(D * D)), np.einsum("jn,jkn->jk", gV, J)[mask]

    for frac in STAGES:
        need = int(round(frac * theta.size)) - int(frozen.sum())
        if need > 0:
            now = np.where(frozen, qv, theta)
            mass = angles_and_mass(decode_pi(now, edges), edges)[1]
            q = np.round(now * scale) / scale
            impact = lam2 * mass * (q - now) ** 2
            impact[frozen] = -np.inf
            pick = np.argsort(-impact, axis=None)[:need]
            t_e = np.vstack([now * scale, w76s.astype(np.float64)])
            c_e = np.vstack([lam2 * mass + 1e-30, np.ones((4, nc))])
            f_e = np.vstack([frozen, np.ones((4, nc), bool)])
            v_e = np.vstack([val, w76s])
            r_, m_ = RD.residues(t_e, c_e, blk_e, pos_e, deg, xdeg, f_e, v_e)
            for fl in pick:
                j, t = divmod(int(fl), nc)
                xs = now[j, t] * scale
                z = int(RD.into_range(RD.residue_round(xs, r_[j, t], m_[t]), m_[t]))
                frozen[j, t] = True
                qv[j, t] = z / scale
                val[j, t] = z
        free = ~frozen
        if free.any():
            base = np.where(frozen, qv, theta)
            r = minimize(fun, theta[free], args=(free, base), jac=True, method="L-BFGS-B",
                         bounds=[(LO / scale, HI / scale)] * int(free.sum()),
                         options={"maxiter": iters, "maxcor": 20, "ftol": 1e-14, "gtol": 1e-10})
            theta = base.copy()
            theta[free] = r.x
        else:
            theta = qv.copy()
    return np.round(qv * scale).astype(np.int64)


def _try_bounded(args):
    H, edges, deg, step, codes, seed, scheme, w76 = args
    lam_e, vec = np.linalg.eigh(H)
    th0, _ = angles_and_mass(vec.T, edges)
    row_at = CT._placement(lam_e, len(lam_e), seed, scheme)
    th0, lam_at = th0[row_at], lam_e[row_at]
    wq = progressive_bounded(th0, edges, lam_at, H, 64.0, 300, deg, w76)
    e_p = error1(decode_pi(np.mod(wq, 128) / 64.0, edges), H)[0]
    wq, e = OL.polish_range(wq, edges, H, 7, 40, step)
    bad = OL.violations_full(wq, w76, codes)
    out = int(np.sum((signed(wq) < LO) | (signed(wq) > HI)))
    return wq, e, dict(route="bounded", scheme=scheme, seed=int(seed), progressive_mHa=1e3 * e_p, polished_mHa=1e3 * e,
                       violations=int(bad), outside=out)


def repair_codes(w1, codes, w76, lam_at, edges):
    """Re-solve every (bit, chain, block) that is not one coset (rows 76-79 frozen in block 4); returns words, counts."""
    w = signed(w1).copy()
    wf = np.vstack([w, signed(w76)])
    mass = angles_and_mass(decode_pi(np.mod(w, 128) / 64.0, edges), edges)[1]
    cw = np.vstack([lam_at[:, None] ** 2 * mass + 1e-30, np.full((4, w.shape[1]), 1e6)])
    fixed = np.zeros(wf.shape, bool)
    fixed[76:] = True
    nfix = nflip = 0
    for bit, _, chains in sorted(codes, key=lambda x: x[0]):
        for h, ch in enumerate(chains):
            ch = np.asarray(ch)
            for blk in range(5):
                rows = np.arange(16 * blk, 16 * blk + 16)
                d_ok = True
                for c0, c1 in zip(ch[:-1], ch[1:]):
                    d = ((np.mod(wf[rows, c0], 128) ^ np.mod(wf[rows, c1], 128)) >> bit) & 1
                    d_ok &= bool((DOB.CW == d[None, :]).all(axis=1).any())
                if d_ok:
                    continue
                nfix += 1
                cur = (np.mod(wf[np.ix_(rows, ch)], 128) >> bit) & 1               # (16, m)
                cost_flip = cw[np.ix_(rows, ch)] * float(1 << bit) ** 2
                cost_flip = np.where(fixed[np.ix_(rows, ch)], 1e12, cost_flip)
                # a flip must stay in range: +-2^b, one of the two is always inside unless the word is pinned
                t0 = np.where(cur == 0, 0.0, cost_flip)[None]
                t1 = np.where(cur == 1, 0.0, cost_flip)[None]
                bits = RD.shared_coset_bits(t0, t1, ch, 11, xor0=False, chain=np.zeros(len(ch), np.int64))[0]
                for i, r in enumerate(rows):
                    for k, c in enumerate(ch):
                        if bits[i, k] != cur[i, k]:
                            assert not fixed[r, c], "a frozen row 76-79 would have to move"
                            x = int(wf[r, c])
                            cands = [x + (1 << bit), x - (1 << bit)]
                            cands = [y for y in cands if LO <= y <= HI and ((y >> bit) & 1) == bits[i, k]
                                     and (y & ((1 << bit) - 1)) == (x & ((1 << bit) - 1))]
                            y = min(cands, key=lambda y: abs(y))
                            wf[r, c] = y
                            nflip += 1
    return wf[:76], nfix, nflip


def placement(lam_e, N, seed, scheme):
    """CD7's placements (chem_cd7_table7._placement) plus 'l4p2' (re-pricer): rows 76-79 fixed leave block 4 two free
    dimensions per chain, at low addresses 0 and 4 (outside the affine span of 12..15): the 21st and 22nd heaviest
    rows sit there and keep their own values; the 10 lightest fill the rest of block 4; blocks 0-3 as light4."""
    if scheme != "l4p2":
        return CT._placement(lam_e, N, seed, scheme)
    rank = np.argsort(-np.abs(lam_e), kind="stable")
    ind = [16 * b_ + l_ for l_ in (0, 1, 2, 4, 8) for b_ in range(4)] + [64, 68]
    blk4 = [a_ for a_ in range(64, N) if a_ not in (64, 68)]
    rest = [a_ for a_ in range(64) if a_ not in ind]
    rest.sort(key=lambda a_: (a_ % 16, a_ // 16))
    heavy, mid, light = rank[:22].copy(), rank[22:N - len(blk4)].copy(), rank[N - len(blk4):].copy()
    if seed > 0:
        rng = np.random.default_rng(seed)
        for t_ in range(5):
            heavy[4 * t_:4 * t_ + 4] = rng.permutation(heavy[4 * t_:4 * t_ + 4])
        rng.shuffle(mid)
        rng.shuffle(light)
    row_at = np.empty(N, np.int64)
    row_at[np.array(ind + rest + blk4)] = np.concatenate([heavy, mid, light])
    return row_at


def route_repair(H, edges, words, spec_text, scheme, seed):
    _, codes = parse_spec(spec_text, np.arange(75))
    RD.HALF = 64
    deg, chain = RD.layout_arrays([[b_, 1, ch] for b_, d_, chains in codes for ch in chains], (75, 8))
    RD.CHAIN = chain
    tool = (deg, np.zeros((15, 8), np.int64))
    N = edges.shape[0] + 1
    lam_e, vec = np.linalg.eigh(H)
    th0, _ = angles_and_mass(vec.T, edges)
    row_at = placement(lam_e, N, seed, scheme)
    th0, lam_at = th0[row_at], lam_e[row_at]
    wq, nf, nm = DOB.progressive_layout(th0, edges, lam_at, H, 64.0, 300, None, tool)
    step = np.ones(75, np.int64)
    for bit, _, chains in codes:
        for ch in chains:
            for c in ch:
                step[c] = max(step[c], 1 << (bit + 1))
    wq2, e_fc = DOB.polish(wq, edges, H, 7, 40, step)
    w1, off, err = obsigned_move(wq2, edges)
    w76, flips76 = repair76(words[76:80], codes)
    viol_moved = link_violations(w1, codes, w76)
    e_moved = error1(decode_pi(np.mod(w1, 128) / 64.0, edges), H)[0]
    # the exact move leaves words at +32 when a parent amplitude is exactly 0: none expected; clamp count reported
    w1s = signed(w1)
    at32 = int(np.sum(w1s > HI))
    w1s = np.minimum(w1s, HI)
    wr, nfix, nflip = repair_codes(w1s, codes, w76, lam_at, edges)
    e_rep = error1(decode_pi(np.mod(wr, 128) / 64.0, edges), H)[0]
    OL.RD.HALF = 32
    wp, e_pol = OL.polish_range(wr, edges, H, 7, 40, step)
    bad = OL.violations_full(wp, w76, codes)
    out = int(np.sum((signed(wp) < LO) | (signed(wp) > HI)))
    return dict(route="repair", scheme=scheme, seed=int(seed), e1_full_circle_mHa=1e3 * e_fc, move_exact=bool(off < 1e-6 and err < 1e-9),
                e1_after_move_mHa=1e3 * e_moved, words_at_plus32=at32, broken_after_move_incl_76_79=viol_moved,
                blocks_repaired=nfix, words_flipped=nflip, e1_after_repair_mHa=1e3 * e_rep,
                e1_after_polish_mHa=1e3 * e_pol, violations=int(bad), outside=out), wp


def save_table(src, wq, out):
    """The table with its 76 one-body rows replaced by wq (7-bit words, mod 128); rows 76-79 and the rest kept."""
    z = dict(np.load(src))
    w = np.asarray(z["angle_words"], np.int64).copy()
    w[:76] = np.mod(np.asarray(wq, np.int64), 128)
    z["angle_words"] = w.astype(z["angle_words"].dtype)
    np.savez(out, **z)
    print("written", out, flush=True)


def _try_repair(args):
    H, edges, words, spec_text, scheme, seed = args
    info, wp = route_repair(H, edges, words, spec_text, scheme, seed)
    return wp, info["e1_after_polish_mHa"] / 1e3, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--item", action="append", default=[])
    ap.add_argument("--route", default="both", choices=["both", "bounded", "repair"])
    ap.add_argument("--tries", type=int, default=8)
    ap.add_argument("--json", default="")
    ap.add_argument("--save-table", default="", help="prefix: write the best rows of each route into the table (npz)")
    ap.add_argument("--scheme", default="both", choices=["both", "light4", "heavy5", "l4p2"])
    a = ap.parse_args()
    import multiprocessing as mp
    mp.set_start_method("fork", force=True)
    from concurrent.futures import ProcessPoolExecutor
    t0 = time.time()
    words, edges, H, e_own = table_H(a.table)
    print(f"table {a.table}: its own one-body rows error_1 {1e3 * e_own:.4f} mHa ({time.time() - t0:.0f} s)", flush=True)
    res = {}
    for x in a.item:
        name, text = (os.path.basename(x), open(x).read().strip()) if x.endswith(".spec") else (x, x)
        res[name] = {}
        jobs_s = ([(s // 2, ("light4", "heavy5")[s % 2]) for s in range(a.tries)] if a.scheme == "both"
                  else [(s, a.scheme) for s in range(a.tries)])
        if a.route in ("both", "bounded"):
            codes, deg, step, coded = OL.setup(text)
            w76, _ = repair76(words[76:80], codes)
            with ProcessPoolExecutor(a.tries) as ex:
                out = list(ex.map(_try_bounded, [(H, edges, deg, step, codes, sd, sc, w76) for sd, sc in jobs_s]))
            ok = [o for o in out if o[2]["violations"] == 0 and o[2]["outside"] == 0]
            best = min(ok or out, key=lambda o: o[1])
            res[name]["bounded"] = dict(best=best[2], tries=[round(o[2]["polished_mHa"], 3) for o in out],
                                        prog=[round(o[2]["progressive_mHa"], 2) for o in out])
            if a.save_table:
                save_table(a.table, best[0], f"{a.save_table}_{name[:24]}_bounded.npz")
            print(f"[{name}] R2 BOUNDED best of {a.tries}: {json.dumps(res[name]['bounded'])} ({time.time() - t0:.0f} s)", flush=True)
        if a.route in ("both", "repair"):
            with ProcessPoolExecutor(a.tries) as ex:
                out = list(ex.map(_try_repair, [(H, edges, words, text, sc, sd) for sd, sc in jobs_s]))
            ok = [o for o in out if o[2]["violations"] == 0 and o[2]["outside"] == 0]
            best = min(ok or out, key=lambda o: o[1])
            res[name]["repair"] = dict(best=best[2], tries=[round(o[2]["e1_after_polish_mHa"], 3) for o in out],
                                       full_circle=[round(o[2]["e1_full_circle_mHa"], 3) for o in out])
            if a.save_table:
                save_table(a.table, best[0], f"{a.save_table}_{name[:24]}_repair.npz")
            print(f"[{name}] R1 REPAIR best of {a.tries}: {json.dumps(res[name]['repair'])} ({time.time() - t0:.0f} s)", flush=True)
        if a.json:
            with open(a.json, "w") as f:
                json.dump(res, f, indent=1)
    print(f"done ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
