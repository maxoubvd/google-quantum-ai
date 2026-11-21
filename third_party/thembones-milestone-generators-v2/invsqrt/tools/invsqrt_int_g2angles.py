"""invsqrt_int (G2, 2026-11-12): rotation angles tuned by a linear programme for the largest minimum margin.

The phase of every input is sum_k units(angle_k) * bit_k(x) (mod 2^64, units of 2^-64 turn) over the circuit's
Z_POW gates (tools/invsqrt_gen_circ.py: Sim.rot, Sim.phase_units), so the error is LINEAR in the angles. Today's
phase module uses the nominal angles -2^-q half turns. Here the 40 rotations (Z_POW gates that are no multiple of
1/4 half turn: Z, S, T stay exact) get integer offsets d_k (units of 2^-64 turn, so the angles stay dyadic and
need no rounding) chosen to maximise m subject to |err0(x) + sum_k d_k bit_k(x)| <= tol(x) - m for all inputs.

  extract <tag> <lut folder> <key=value line ...>   build the line, simulate ALL inputs once, save
        ~/qcc/ws/invsqrt_gen/int/g2ang/<tag>.npz (bits of the 40 rotations, err0, tol, positions, stages)
  lp <tag> [start=200000] [rounds=12]               constraint generation with scipy's HiGHS: prints m and
        the score change; writes <tag>.angles.json = {"lo": {pos: [num, den]}, "hi": {...}} (half turns)
"""
import json
import math
import os
import sys
import time
from fractions import Fraction

import numpy as np

sys.path.insert(0, "<HOME>/Quantum/tools")
OUT = os.path.expanduser("~/qcc/ws/invsqrt_gen/int/g2ang")
os.makedirs(OUT, exist_ok=True)


def extract(tag, lutdir, args):
    os.environ["INVSQRT_GEN_LUT"] = lutdir
    import invsqrt_gen_build as B
    import invsqrt_gen_circ as C
    import invsqrt_gen_model as M
    kw = dict(a.split("=") for a in args)
    tabname = kw.pop("tabs", None)
    kw = {k: int(v) for k, v in kw.items()}
    pf = set(M.Params.__dataclass_fields__)
    P = M.Params(**{k: v for k, v in kw.items() if k in pf})
    tabs = M.load_tables(tabname)
    c = B.VARIANTS["g6"](P, tabs, **{k: v for k, v in kw.items() if k not in pf})
    # stage of every Z_POW op, in op order (= order of Sim.rot)
    stage, info = None, []
    for i, op in enumerate(c.ops):
        if op[0] == "@stage":
            stage = op[1]
        elif op[0] == "Z_POW":
            info.append((i, stage, op[4]))
    xs = np.arange(1, M.MAX_X + 1, dtype=np.int64)
    target, tol = M.reference_all()
    t0 = time.time()
    sim = C.Sim(xs.astype(np.uint64), c.num_qubits(), hmr="rand", seed=1)
    sim.run(c.ops)
    assert len(sim.rot) == len(info)
    th = sim.phase_units()
    err0 = (th - target).view(np.int64)
    used = np.abs(err0)
    margin = tol.astype(np.int64) - used
    print(f"simulated {len(xs)} inputs in {time.time() - t0:.0f} s; Z_POW gates {len(info)}; min margin "
          f"{int(margin.min())} at x = {int(xs[int(np.argmin(margin))])}; failing {int((used > tol).sum())}",
          flush=True)
    tun = [k for k, (_, _, a) in enumerate(info) if (a * 4) % 1 != 0]
    bits = np.zeros((len(tun), len(xs)), dtype=np.uint8)
    for r, k in enumerate(tun):
        b = sim.rot[k][1]
        if b:
            bits[r] = C.unpack(b, len(xs))
    pos = []
    for k in tun:
        a = info[k][2]
        q = round(-math.log2(abs(float(a))))
        assert Fraction(-1, 1 << q) == a, ("not a nominal angle", a)
        pos.append(q)
    stages = [info[k][1] for k in tun]
    np.savez_compressed(os.path.join(OUT, f"{tag}.npz"), bits=np.packbits(bits, axis=1), n=len(xs), err0=err0,
                        tol=tol.astype(np.int64), pos=np.array(pos), stages=np.array(stages))
    print(f"tunable rotations {len(tun)}: " + " ".join(f"{s}:{q}" for s, q in zip(stages, pos)))


def lp(tag, start=200000, rounds=12):
    from scipy.optimize import linprog
    d = np.load(os.path.join(OUT, f"{tag}.npz"))
    n = int(d["n"])
    bits = np.unpackbits(d["bits"], axis=1)[:, :n]                            # K x n, uint8
    err0 = d["err0"].astype(np.float64)

    def err_of(dv):
        e = err0.copy()
        for k in range(len(dv)):
            if dv[k]:
                e += dv[k] * bits[k]
        return e
    tol = d["tol"].astype(np.float64)
    K = bits.shape[0]
    m0 = tol - np.abs(err0)
    print(f"K {K}, nominal min margin {m0.min():.0f}", flush=True)
    act = set(np.argsort(m0)[:start].tolist())
    dvec = np.zeros(K)
    for rnd in range(rounds):
        idx = np.array(sorted(act))
        Bm = bits[:, idx].T.astype(np.float64)                                  # rows: inputs
        # variables: d_1..d_K, m ; maximise m  <=>  minimise -m
        # err0 + B d + m <= tol   and   -err0 - B d + m <= tol
        A = np.vstack([np.hstack([Bm, np.ones((len(idx), 1))]), np.hstack([-Bm, np.ones((len(idx), 1))])])
        b = np.concatenate([tol[idx] - err0[idx], tol[idx] + err0[idx]])
        cvec = np.zeros(K + 1)
        cvec[-1] = -1.0
        bounds = [(-2.0 ** 40, 2.0 ** 40)] * K + [(None, None)]
        t0 = time.time()
        res = linprog(cvec, A_ub=A, b_ub=b, bounds=bounds, method="highs")
        assert res.status == 0, res.message
        dvec = np.round(res.x[:K])
        err = err_of(dvec)
        marg = tol - np.abs(err)
        mmin = marg.min()
        viol = np.argsort(marg)[:20000]
        new = [int(i) for i in viol if int(i) not in act]
        print(f"round {rnd}: {len(idx)} constraints, LP m {res.x[-1]:.0f}, all-input m (rounded) {mmin:.0f}, "
              f"{time.time() - t0:.0f} s, adding {len(new)}", flush=True)
        if not new or marg[viol[0]] >= res.x[-1] - 2 * K:
            break
        act.update(new)
    g0 = m0.min()
    nrot = K
    tper = lambda mm: 0.56 * math.log2(nrot * 2 ** 64 / (math.pi * mm)) + 5.3
    dT = nrot * (tper(g0) - tper(mmin)) / 4
    print(f"min margin {g0:.0f} -> {mmin:.0f} ({mmin / g0:.4f}x); T per rotation {tper(g0):.4f} -> {tper(mmin):.4f}; "
          f"-{dT:.3f} Toffoli-eq = -{59 * dT:.1f} points at Q 59")
    pos, stages = d["pos"], d["stages"]
    outj = {"lo": {}, "hi": {}, "min_margin": int(mmin), "nominal_min_margin": int(g0), "tag": tag}
    for k in range(K):
        q = int(pos[k])
        units = -(1 << (64 - q - 1)) + int(dvec[k])      # -2^-q half turn = -2^-(q+1) turn = -2^(63-q) units
        fr = Fraction(units, 1 << 63)                    # half turns
        key = "lo" if stages[k] == "PH2" else "hi"
        assert str(q) not in outj[key], ("two rotations at one position", key, q)
        outj[key][str(q)] = [fr.numerator, fr.denominator]
    json.dump(outj, open(os.path.join(OUT, f"{tag}.angles.json"), "w"), indent=0)
    print("stages of the tuned rotations:", sorted(set(stages.tolist())))


if __name__ == "__main__":
    if sys.argv[1] == "extract":
        extract(sys.argv[2], sys.argv[3], sys.argv[4:])
    else:
        kw = dict(a.split("=") for a in sys.argv[3:])
        lp(sys.argv[2], **{k: int(v) for k, v in kw.items()})
