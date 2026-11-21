"""F1 (invsqrt fork): truth table of one wire as a function of other wires at a given op of a .kmb, written as a
minimum-AND problem for builder F2's `tools/invsqrt_fl_minand.py json` (rows = distinct input patterns).
Simulates all 3,133,452 inputs to op T (F3's simulator, HMR outcomes 0; data wires do not depend on outcomes).
    python invsqrt_fk_signsat.py <file.kmb> <T> <target=qN> <name=qN,name=qN,...> <out.json>
Prints the number of patterns and whether the target is a function of the inputs.
"""
import json
import sys

import numpy as np

sys.path.insert(0, "<HOME>/Quantum/tools")
import invsqrt_fw_sim as F  # noqa: E402


def wire(c, name):
    return int(name[1:]) if name.startswith("q") else int(name)


def main(path, T, tgt, ins, out):
    c = F.Circ(path)
    s = F.Sim(c, hmr="zero", phase=False, track=False)
    s.run_to(int(T))
    n = F.MAXX
    tname, tw = tgt.split("=")
    pairs = [x.split("=") for x in ins.split(",")]
    cols = [F.unpack(s.V[wire(c, w)], n).astype(np.int64) for _, w in pairs]
    if tw.startswith("ref:v"):               # target = root bit v_i of the integer model (wire may be retired here)
        y = (ref_v(n) >> int(tw[5:])) & 1
    elif "@" in tw:                          # target = a wire's value at ANOTHER op (e.g. a sign at its rebuild)
        w2, op2 = tw.split("@")
        s2 = F.Sim(c, hmr="zero", phase=False, track=False)
        s2.run_to(int(op2))
        y = F.unpack(s2.V[wire(c, w2)], n).astype(np.int64)
    else:
        y = F.unpack(s.V[wire(c, tw)], n).astype(np.int64)
    key = np.zeros(n, dtype=np.int64)
    for i, col in enumerate(cols):
        key |= col << i
    pats, inv = np.unique(key, return_inverse=True)
    yv = np.full(len(pats), -1)
    func = True
    for p_i in range(len(pats)):
        vals = np.unique(y[inv == p_i])
        if len(vals) > 1:
            func = False
        yv[p_i] = vals[0]
    print(f"{tname} at op {T}: {len(pats)} patterns over {len(pairs)} inputs; function: {func}; "
          f"ones {int(yv.sum())}")
    if not func:
        return
    R = len(pats)
    tt = []
    for i in range(len(pairs)):
        v = 0
        for r in range(R):
            if (pats[r] >> i) & 1:
                v |= 1 << r
        tt.append(v)
    t = 0
    for r in range(R):
        if yv[r]:
            t |= 1 << r
    json.dump(dict(rows=R, inputs=[dict(name=nm, tt=hex(v)) for (nm, _), v in zip(pairs, tt)], target=hex(t)),
              open(out, "w"))
    print("wrote", out)


def ref_v(n):
    """v = isqrt(16 n') - 4096 with n' = x 4^k in [2^20, 2^22), for x = 1..n (integer model)."""
    import math
    xs = np.arange(1, n + 1, dtype=np.int64)
    nn = xs.copy()
    while True:
        m = nn < (1 << 20)
        if not m.any():
            break
        nn[m] *= 4
    return np.array([math.isqrt(int(16 * t)) for t in nn], dtype=np.int64) - 4096


def label(path, T):
    """Which wires hold root bit v_i (v = isqrt(16 n) - 4096, n = x 4^k in [2^20, 2^22)) or a flag bit at op T."""
    import math
    c = F.Circ(path)
    s = F.Sim(c, hmr="zero", phase=False, track=False)
    s.run_to(int(T))
    n = F.MAXX
    xs = np.arange(1, n + 1, dtype=np.int64)
    k = np.zeros(n, dtype=np.int64)
    nn = xs.copy()
    while True:
        m = nn < (1 << 20)
        if not m.any():
            break
        nn[m] *= 4
        k[m] += 1
    v = np.array([math.isqrt(int(16 * t)) for t in nn], dtype=np.int64) - 4096
    refs = {f"v{i}": (v >> i) & 1 for i in range(12)}
    kk = np.where(k == 10, 0b1010, k)
    refs.update({f"f{i}": (kk >> i) & 1 for i in range(4)})
    for w in range(c.nq):
        b = F.unpack(s.V[w], n).astype(np.int64)
        if not b.any():
            continue
        hit = [nm for nm, r in refs.items() if np.array_equal(b, r)] + \
              [nm + "'" for nm, r in refs.items() if np.array_equal(b, 1 - r)]
        if hit:
            print(f"q{w}: {' '.join(hit)}")


if __name__ == "__main__":
    if sys.argv[1] == "label":
        label(sys.argv[2], sys.argv[3])
    else:
        main(*sys.argv[1:6])
