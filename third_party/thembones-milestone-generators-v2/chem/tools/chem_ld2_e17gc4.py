#!/usr/bin/env python3
"""chem_ld2_e17gc4.py: builder C2 (seventeenth holder, Nov 17). ROUTE (b): THE C-4 CLOSURE BUILT INTO THE PLAIN
GROUP ENGINE (the 68-AND high cache on loaders the plain compiler takes: wired 235 -> 234, coded 245 -> 244).

Why it can work [M + P, 01:15Z]: on a 1,219-row table every loader admits a block set V (chem_ld2_c4.vset: the 13
padding rows of block 76 are free), and each loader's padded delta depends on its own V only (pad_table), so every
loader can be closed alone; an affine bit-0 wire is extended affinely by the padding (want76 = D(16*76 + l) [P]), so
the wires survive. With the padded table the solved high parts lie in H automatically (independent features, V-sums
of the table 0 [P]). What is left is the plan condition tau == kappa on every output wire. The group engine's
products share no low planes (no sector moves), so the closure is made INSIDE it: every group is solved again (its
random relations and orders redrawn) until tau == kappa on the wires whose basis support that group completes (the
products of a finished group never change for those wires: kills only extend earlier products to LATER outputs).
Expected tries per group 2^(wires it completes), about 64.

    python chem_ld2_e17gc4.py <table.npz (final words, 1,219 rows)> --ks 4,6,8 [--tries 400] [--effort 3]
                              [--seeds 0-3] [--emit DIR] [--json out.json]

Per loader: the plain group engine's count on the same padded table (69-AND cache) and the closed count (68-AND
cache, emitted by chem_ld2_c4.emit_c4, which checks every wire on every address and tau == kappa per wire).
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_ld2_wwidth  # noqa: E402,F401  (env C2_NOCAR=1: no carrier read, as everywhere in my chain)
import chem_ld2_c4 as c4  # noqa: E402


def rng_list(s):
    out = []
    for part in s.split(","):
        if "-" in part:
            a, b = part.split("-")
            out += list(range(int(a), int(b) + 1))
        elif part:
            out.append(int(part))
    return out


def vset_rel(words, k, an):
    """02:20Z: a C-4 block set V that the closure can meet on a loader with RELATIONS among its outputs. For every
    combination u of output wires whose mixed parts cancel (sum_q u_q G[q] = 0: a relation f(l) + g(block), or a wire),
    tau on that combination is 0 for any plan, so its kappa must be 0: the combined delta's V-sum at low row 0 vanishes
    [P]. chem_ld2_c4pub.vset's system plus one row per such combination (RK7's F-table relations have non-affine g: the
    plain vset's V fails them [M k4 of rk7_3: no closed solve in 400 tries]). None if no V."""
    import chem_ld2_c4pub as cp
    from chem_ld2_hb import HB, LAST
    F, valid = cp.ct.delta_bits(words, k, core.BROT - 1)
    rows = [np.ones(HB, np.uint8)]
    pidx = np.arange(HB)
    for j in range(7):
        rows.append(((pidx >> j) & 1).astype(np.uint8))
    for b in range(F.shape[0]):
        for l in range(1, 16):
            if not valid[LAST, l]:
                continue
            rows.append(F[b, :, l] ^ F[b, :, 0])
    nw, nb = len(an.out_wires), an.X0.shape[1]
    rows_t = [sum(1 << q for q in range(nw) if (an.G[q] >> o) & 1) for o in range(nb)]
    A = words.shape[0]
    J = core.J
    prev = words[:, J * (k - 1):J * k] if k else np.zeros((A, J), np.int64)
    d = prev ^ words[:, J * k:J * k + J]
    for u in core.nullspace(rows_t, nw):
        col = np.zeros(A, np.uint8)
        for q in core.bits(u):
            j, b = an.out_wires[q]
            col ^= ((d[:, j] >> b) & 1).astype(np.uint8)
        rows.append(np.array([col[16 * p] if 16 * p < A else 0 for p in range(HB)], np.uint8))
    N = cp.ct.nullspace_gf2(np.array(rows, np.uint8))
    sols = [v for v in N if v[LAST]]
    if not sols:
        return None
    return sum(1 << p for p in range(HB) if sols[0][p])


class TauDecomp(core.Decomp):
    """core.Decomp whose every group is retried until tau == kappa on the wires that the group completes."""

    def __init__(self, an, seed, effort, V, kap, tries):
        super().__init__(an, seed, effort)
        self.V, self.kap, self.tries = V, kap, tries
        self.done = 0                        # basis outputs whose groups are solved
        self.retries = []

    def _tau_bad(self, wmask):
        P = c4.normalize([list(p) for p in self.prods])
        return (c4.tau_state(self.an, P, self.V) ^ self.kap) & wmask

    def solve_group(self, outs, coords, dl=0, d=0, tag=""):
        omask = sum(1 << o for o in outs)
        newly = self.done | omask
        wmask = 0
        for q, g in enumerate(self.an.G):
            if g and (g & omask) and not (g & ~newly):
                wmask |= 1 << q
        snap = (self.X.copy(), [list(p) for p in self.prods], self.PM.copy())
        made = None
        for t in range(self.tries):
            try:
                made = super().solve_group(outs, coords, dl, d, tag)
            except (RuntimeError, AssertionError):
                made = None
            if made is not None and not self._tau_bad(wmask):
                self.retries.append(t + 1)
                break
            self.X, self.prods, self.PM = snap[0].copy(), [list(p) for p in snap[1]], snap[2].copy()
            made = None
        if made is None:
            raise RuntimeError(f"no closed solve for group {outs} in {self.tries} tries")
        self.done = newly
        return made


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--ks", required=True)
    ap.add_argument("--tries", type=int, default=400)
    ap.add_argument("--effort", type=int, default=3)
    ap.add_argument("--seeds", default="0-3")
    ap.add_argument("--emit", default=None)
    ap.add_argument("--json", default=None)
    # 01:32Z: PLACED KILLS (tools/chem_ld2_e17kill.py) with the closure: plan 6 / 6 / 7 + kill t0 / 6 + kill t1 + 2 row
    # kills, groups 1-2 off the planes silent on t0, t1: 158 products = 232 closed on k14 of wr7_df2 [M, 19 of 60 seeds]
    ap.add_argument("--kplan", default="", help='e.g. "6,6,7:1:0:0,6:1:2:1" (empty = off)')
    ap.add_argument("--kavoid", default="0,1")
    ap.add_argument("--kseeds", default="0-15")
    a = ap.parse_args()
    core.J, core.BROT = 5, 7
    words = np.asarray(np.load(a.table)["angle_words"], np.int64)
    res = {}
    for k in rng_list(a.ks):
        t0 = time.time()
        V = c4.vset(words, k)
        if V is None:
            print(f"k={k}: no block set V", flush=True)
            continue
        Wp = c4.pad_table(words, {k: V})
        an = core.analyze(Wp, k, 4)
        # 02:20Z: relations among the outputs constrain V too (vset_rel); the padded table is made again with it
        try:
            V2 = vset_rel(words, k, an)
        except Exception as exc:  # noqa: BLE001
            print(f"k={k}: vset_rel failed ({exc}); plain V kept", flush=True)
            V2 = None
        if V2 is not None and V2 != V:
            V = V2
            Wp = c4.pad_table(words, {k: V})
            an = core.analyze(Wp, k, 4)
        bad = c4.check_padded(Wp, {k: V})
        kap = c4.kappa_wires(an, V)
        # the plain group engine on the same padded table (69-AND cache), for the comparison
        plain = an.infra + len(core.Decomp(an, 0, max(2, a.effort)).run())
        best = None
        for s in rng_list(a.seeds):
            core.MARGIN, core.MAX_GROUP = 2, 8
            td = TauDecomp(an, s, max(2, a.effort), V, kap, a.tries)
            try:
                prods = td.run()
            except RuntimeError as exc:
                print(f"k={k} seed={s}: {exc}", flush=True)
                continue
            P = c4.normalize([list(p) for p in prods])
            mu = c4.tau_state(an, P, V) ^ kap
            try:
                text, st = c4.emit_c4(an, prods, V)
            except (AssertionError, ValueError) as exc:
                print(f"k={k} seed={s}: {len(prods)} products, mismatch {bin(mu).count('1')}, emit_c4 failed: {exc}",
                      flush=True)
                continue
            tof = st["toffoli"]
            print(f"k={k} seed={s}: closed {tof} (plain on the padded table {plain}; {len(prods)} products; "
                  f"retries per group {td.retries}; {time.time() - t0:.1f} s)", flush=True)
            if best is None or tof < best[0]:
                best = (tof, s, text, td.retries, len(prods))
        if a.kplan:
            import chem_ld2_e17kill as kl                        # lazy: chem_ld2_e17kill imports this module
            avoid = [int(x) for x in a.kavoid.split(",") if x.strip()]
            # 02:14Z: the given plan, then the model-optimal plan for this loader's number of generic outputs
            # (chem_ld2_e17kill.PLANS_BY_N); each stops at its model floor
            n_gen = sum(1 for x in an.out_level if x > 2)
            ad = kl.PLANS_BY_N.get(n_gen)
            plans = ([ad[0]] if ad else []) + ([a.kplan] if not ad or ad[0] != a.kplan else [])
            reached = False
            for pstr in plans:                                  # 02:21Z: the adaptive plan first
                if reached:
                    break
                plan = kl.parse_plan(pstr)
                floor = ad[1] if ad and ad[0] == pstr else None
                for s in rng_list(a.kseeds):
                    kd = kl.KillDecomp(an, s, max(2, a.effort), plan, avoid, V, kap, a.tries)
                    try:
                        prods = kd.run()
                        text, st = c4.emit_c4(an, prods, V)
                    except (RuntimeError, AssertionError, ValueError):
                        continue
                    tof = st["toffoli"]
                    if best is None or tof < best[0]:
                        best = (tof, f"kplan {pstr} seed {s}", text, kd.retries, len(prods))
                        print(f"k={k} ({n_gen} generic) kplan {pstr} seed={s}: closed {tof} ({len(prods)} products)",
                              flush=True)
                    if floor is not None and len(prods) <= getattr(kd, "n_low", 0) + floor:
                        reached = True
                        break
        # 02:45Z: the placed kills WITHOUT the closure (69-AND cache, chem_loaders_core.emit, classical check) on the
        # unpadded words, tried only when its model floor can beat the closed best: a loader whose relations leave no
        # consistent V (k12 of the six-relation test table: no closed solve) still gains the kills over the plain engine
        if a.kplan:
            import chem_ld2_e17kill as kl
            an_u = core.analyze(words, k, 4)
            n_gen_u = sum(1 for x in an_u.out_level if x > 2)
            n_low_u = sum(1 for x in an_u.out_level if x <= 2)
            ad_u = kl.PLANS_BY_N.get(n_gen_u)
            floor_u = an_u.infra + ad_u[1] + 2 * n_low_u if ad_u else None
            if best is None or floor_u is None or floor_u < best[0]:
                avoid = [int(x) for x in a.kavoid.split(",") if x.strip()]
                done_u = False
                for pstr in ([ad_u[0]] if ad_u else []) + ([a.kplan] if not ad_u or ad_u[0] != a.kplan else []):
                    if done_u:
                        break
                    plan = kl.parse_plan(pstr)
                    for s in rng_list(a.kseeds):
                        kd = kl.KillDecomp(an_u, s, max(2, a.effort), plan, avoid)
                        try:
                            prods = kd.run()
                            text, st = core.emit(an_u, prods)
                        except (RuntimeError, AssertionError, ValueError):
                            continue
                        tof = st["toffoli"]
                        if best is None or tof < best[0]:
                            best = (tof, f"unclosed kplan {pstr} seed {s}", text, [], len(prods))
                            print(f"k={k} ({n_gen_u} generic) UNCLOSED kplan {pstr} seed={s}: {tof} ({len(prods)} "
                                  f"products, 69-AND cache)", flush=True)
                        if ad_u and pstr == ad_u[0] and len(prods) <= getattr(kd, "n_low", 0) + ad_u[1]:
                            done_u = True
                            break
        if best is None:
            res[k] = dict(plain=plain, closed=None)
            continue
        res[k] = dict(plain=plain, closed=best[0], seed=best[1], retries=best[3], products=best[4],
                      padded_bad=bad.get(k), V=hex(V))
        print(f"BEST k={k}: {best[0]} closed (plain {plain})", flush=True)
        if a.emit:
            os.makedirs(a.emit, exist_ok=True)
            with open(os.path.join(a.emit, f"angle_k_{k}.kmx"), "w") as fh:
                fh.write(best[2])
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
