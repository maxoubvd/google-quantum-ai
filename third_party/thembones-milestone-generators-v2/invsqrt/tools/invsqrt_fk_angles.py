"""invsqrt_fk (F1): rotation angles of a finished invsqrt circuit tuned by an LP for the largest minimum margin.

The verifier prices every generic rotation at 0.56 log2(1/eps) + 5.3 T with eps = pi * m / 2^64 / n_rot, where
m = min_x (tol(x) - |err(x)|) in units of 2^-64 turn. The phase of input x is sum_k angle_k * bit_k(x) plus exact
Clifford/T parts, so err(x) is LINEAR in the generic angles. We give each generic Z_POW an integer offset d_k
(units of 2^-64 turn; the angles stay dyadic) chosen to maximise m (constraint generation, scipy HiGHS).
Same principle as G2's tools/invsqrt_int_g2angles.py; this file works on any kmx text (no generator needed).

  extract KMX OUT.npz          simulate ALL 3,133,452 inputs (HMR outcomes 0; the circuits are exact for every
                               transcript), save the bits of every generic Z_POW, err0 and tol
  lp OUT.npz [start=200000] [rounds=14]       -> OUT.angles.json {line_index: offset_units, ...} + margin report
  apply KMX OUT.angles.json NEW.kmx [NEW.kmb] write the circuit with the new angles (and compile it)
Runs on box 1 (python ~/qcc/.venv_work/bin/python) or in WSL; uses our own parser and simulator only.
"""
import json
import math
import os
import sys
import time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_lib as L          # noqa: E402  (our kmx parser)
import invsqrt_gen_circ as C     # noqa: E402  (our bitset simulator)

MAX_X = 3133452
M64 = (1 << 64) - 1


def reference(cache):
    """target and tol (uint64, units of 2^-64 turn) for x = 1..MAX_X, exactly as the official verifier."""
    if os.path.exists(cache):
        d = np.load(cache)
        return d["target"], d["tol"]
    tg = np.empty(MAX_X, dtype=np.uint64)
    tl = np.empty(MAX_X, dtype=np.uint64)
    for x in range(1, MAX_X + 1):
        r = math.isqrt((1 << 130) // x)
        tg[x - 1] = ((r + 1) >> 1) & M64
        tl[x - 1] = min(r >> 25, 1 << 63)
    np.savez(cache, target=tg, tol=tl)
    return tg, tl


def is_generic(angle):
    return angle is not None and (angle * 4) % 1 != 0


def extract(kmx, out):
    header, ops = L.parse_kmx(kmx)
    ops = [op for op in ops if op[0] not in ("INFO_TAG",)]
    nq = 1 + max(q for op in ops for q in op[1])
    xs = np.arange(1, MAX_X + 1, dtype=np.uint64)
    target, tol = reference(os.path.join(os.path.dirname(os.path.abspath(out)), "ref_all.npz"))
    t0 = time.time()
    sim = C.Sim(xs, nq, hmr="zero")
    sim.run(ops)
    th = sim.phase_units()
    err0 = (th - target).view(np.int64)
    margin = tol.astype(np.int64) - np.abs(err0)
    zl = [i for i, op in enumerate(ops) if op[0] == "Z_POW"]
    assert len(zl) == len(sim.rot)
    gen = [r for r, i in enumerate(zl) if is_generic(ops[i][4])]
    bits = np.zeros((len(gen), MAX_X), dtype=np.uint8)
    for j, r in enumerate(gen):
        b = sim.rot[r][1]
        if b:
            bits[j] = C.unpack(b, MAX_X)
    angles = [str(ops[zl[r]][4]) for r in gen]
    np.savez_compressed(out, bits=np.packbits(bits, axis=1), err0=err0, tol=tol.astype(np.int64),
                        zpos=np.array([zl[r] for r in gen]), angles=np.array(angles))
    print(f"{kmx}: {MAX_X} inputs in {time.time() - t0:.0f} s; Z_POW {len(zl)}, generic {len(gen)}; "
          f"min margin {int(margin.min())} at x = {int(np.argmin(margin)) + 1}; failing {int((margin < 0).sum())}; "
          f"worst ratio {float((np.abs(err0) / tol).max()):.5f}", flush=True)


def tper(m, nrot):
    return 0.56 * math.log2(nrot * 2 ** 64 / (math.pi * m)) + 5.3


def lp(npz, start=200000, rounds=14):
    from scipy.optimize import linprog
    d = np.load(npz)
    err0 = d["err0"].astype(np.float64)
    tol = d["tol"].astype(np.float64)
    bits = np.unpackbits(d["bits"], axis=1)[:, :len(err0)]
    K = bits.shape[0]
    m0 = tol - np.abs(err0)
    print(f"K {K}, nominal min margin {m0.min():.0f}", flush=True)
    act = set(np.argsort(m0)[:start].tolist())
    best = (m0.min(), np.zeros(K))
    for rnd in range(rounds):
        idx = np.array(sorted(act))
        Bm = bits[:, idx].T.astype(np.float64)
        A = np.vstack([np.hstack([Bm, np.ones((len(idx), 1))]), np.hstack([-Bm, np.ones((len(idx), 1))])])
        b = np.concatenate([tol[idx] - err0[idx], tol[idx] + err0[idx]])
        cvec = np.zeros(K + 1)
        cvec[-1] = -1.0
        t0 = time.time()
        res = linprog(cvec, A_ub=A, b_ub=b, bounds=[(-2.0 ** 40, 2.0 ** 40)] * K + [(None, None)], method="highs")
        assert res.status == 0, res.message
        dv = np.round(res.x[:K])
        err = err0.copy()
        for k in range(K):
            if dv[k]:
                err += dv[k] * bits[k]
        marg = tol - np.abs(err)
        mmin = marg.min()
        if mmin > best[0]:
            best = (mmin, dv.copy())
        viol = np.argsort(marg)[:20000]
        new = [int(i) for i in viol if int(i) not in act]
        print(f"round {rnd}: {len(idx)} constraints, LP m {res.x[-1]:.0f}, all-input m {mmin:.0f}, "
              f"{time.time() - t0:.0f} s, adding {len(new)}", flush=True)
        if not new or marg[viol[0]] >= res.x[-1] - 2 * K:
            break
        act.update(new)
    g0, (mmin, dv) = m0.min(), best
    dT = K * (tper(g0, K) - tper(mmin, K)) / 4
    print(f"min margin {g0:.0f} -> {mmin:.0f} ({mmin / g0:.4f}x); T per rotation {tper(g0, K):.4f} -> "
          f"{tper(mmin, K):.4f}; -{dT:.3f} Toffoli-eq = -{46 * dT:.1f} points at Q 46", flush=True)
    outj = {"offsets_units": {str(int(z)): int(v) for z, v in zip(d["zpos"], dv) if v},
            "min_margin": int(mmin), "nominal_min_margin": int(g0), "K": K}
    json.dump(outj, open(npz.replace(".npz", ".angles.json"), "w"), indent=0)


def solve(err0, tol, bits, start=150000, rounds=10):
    """max-min margin LP over offsets of the rows of bits; returns (all-input min margin, offsets)."""
    from scipy.optimize import linprog
    K = bits.shape[0]
    m0 = tol - np.abs(err0)
    act = set(np.argsort(m0)[:start].tolist())
    best = (m0.min(), np.zeros(K))
    for rnd in range(rounds):
        idx = np.array(sorted(act))
        Bm = bits[:, idx].T.astype(np.float64)
        A = np.vstack([np.hstack([Bm, np.ones((len(idx), 1))]), np.hstack([-Bm, np.ones((len(idx), 1))])])
        b = np.concatenate([tol[idx] - err0[idx], tol[idx] + err0[idx]])
        cvec = np.zeros(K + 1)
        cvec[-1] = -1.0
        res = linprog(cvec, A_ub=A, b_ub=b, bounds=[(-2.0 ** 44, 2.0 ** 44)] * K + [(None, None)], method="highs")
        if res.status != 0:
            return best
        dv = np.round(res.x[:K])
        err = err0.copy()
        for k in range(K):
            if dv[k]:
                err += dv[k] * bits[k]
        marg = tol - np.abs(err)
        if marg.min() > best[0]:
            best = (marg.min(), dv.copy())
        viol = np.argsort(marg)[:20000]
        new = [int(i) for i in viol if int(i) not in act]
        if not new or marg[viol[0]] >= res.x[-1] - 2 * K:
            break
        act.update(new)
    return best


def drops(npz):
    """For each generic rotation k: remove it (angle 0), re-solve the others, report margin and cost change."""
    d = np.load(npz)
    err0 = d["err0"].astype(np.float64)
    tol = d["tol"].astype(np.float64)
    bits = np.unpackbits(d["bits"], axis=1)[:, :len(err0)]
    K = bits.shape[0]
    g0 = (tol - np.abs(err0)).min()
    base = K * tper(g0, K) / 4
    for k in range(K):
        a = Fraction(str(d["angles"][k]))
        units = int(a * (1 << 63)) % (1 << 64)
        if units >= 1 << 63:
            units -= 1 << 64
        e = err0 - units * bits[k].astype(np.float64)
        e = (e + 2.0 ** 63) % 2.0 ** 64 - 2.0 ** 63
        keep = [j for j in range(K) if j != k]
        t0 = time.time()
        m, _ = solve(e, tol, bits[keep])
        cost = (K - 1) * tper(m, K - 1) / 4 if m > 0 else float("nan")
        print(f"drop {k:2d} (line {int(d['zpos'][k])}, angle {float(a):.3e} half turns, ones {int(bits[k].sum())}): "
              f"min margin {m:.0f}; rotation charge {cost:.2f} vs {base:.2f} -> {46 * (cost - base):+.1f} points "
              f"({time.time() - t0:.0f} s)", flush=True)


def apply(kmx, angles_json, new_kmx, new_kmb=None):
    header, ops = L.parse_kmx(kmx)
    ops = [op for op in ops if op[0] not in ("INFO_TAG",)]
    off = {int(k): v for k, v in json.load(open(angles_json))["offsets_units"].items()}
    out = []
    for i, op in enumerate(ops):
        if i in off:
            name, qs, bit, cond, ang = op
            assert name == "Z_POW" and is_generic(ang)
            ang = (ang + Fraction(off[i], 1 << 63)) % 2          # units of 2^-64 turn = 2^-63 half turn
            op = (name, qs, bit, cond, ang)
        out.append(op)
    L.write_kmx(new_kmx, header, out)
    print(f"wrote {new_kmx}: {len(off)} angles changed")
    if new_kmb:
        kdir = os.path.expanduser("~/qcc/data/kickmix")
        sys.path.insert(0, kdir)
        import kickmix as km
        c = km.Circuit.from_file(new_kmx, format="kmx")
        c.to_file(new_kmb, format="kmb")
        print(f"compiled {new_kmb}: Q {c.num_qubits}")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "extract":
        extract(sys.argv[2], sys.argv[3])
    elif cmd == "lp":
        kw = dict(a.split("=") for a in sys.argv[3:])
        lp(sys.argv[2], **{k: int(v) for k, v in kw.items()})
    elif cmd == "apply":
        apply(*sys.argv[2:])
    elif cmd == "drops":
        drops(sys.argv[2])
