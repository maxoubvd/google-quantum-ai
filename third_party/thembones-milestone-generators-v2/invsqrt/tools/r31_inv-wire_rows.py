"""R31 seat INV-WIRE (2026-11-13): per-row Toffolis and wires of the framed adder rows of a .kmb stage.

Our own code; circuit file read as data. A framed row is [X t; CX s t] over the target wires, the adder, the same
frame again. Rows are found by the two frame runs (same sign wire s, >= 3 targets). For each row: op range,
sign wire, target width, exact expected Toffolis from the frame start to the frame end, the expected Toffolis
between the previous row's end and this row's start (feed carry, rebuilds), width at the frame (wires live
there = allocated) and the peak width inside.

    python r31_inv-wire_rows.py <file.kmb> <width.npy> <a> <b>
"""
import importlib.util
import sys

import numpy as np

sys.path.insert(0, "<HOME>/Quantum/tools")
import invsqrt_fw_sim as F  # noqa: E402

spec = importlib.util.spec_from_file_location("r31tof", "<HOME>/Quantum/tools/r31_inv-wire_tof.py")
T = importlib.util.module_from_spec(spec)
spec.loader.exec_module(T)


def frames(c, a, b):
    """Maximal runs X t; CX s t; X t'; CX s t' ... with one sign wire s. Returns (start, end, s, targets)."""
    out = []
    i = a
    ops = c.ops
    while i < b - 1:
        if ops[i][0] == 0x11 and ops[i + 1][0] == 0x21 and ops[i + 1][2][1] == ops[i][2][0]:
            s = ops[i + 1][2][0]
            j = i
            tg = []
            while j < b - 1 and ops[j][0] == 0x11 and ops[j + 1][0] == 0x21 and ops[j + 1][2] == (s, ops[j][2][0]):
                tg.append(ops[j][2][0])
                j += 2
            if len(tg) >= 3:
                out.append((i, j, s, tg))
                i = j
                continue
        i += 1
    return out


def frames_cx(c, a, b):
    """Maximal runs CX s t; CX s t' ... (one control, distinct targets, >= 4): the product rows' sign frames."""
    out = []
    i = a
    ops = c.ops
    while i < b:
        if ops[i][0] == 0x21:
            s = ops[i][2][0]
            j = i
            tg = []
            while j < b and ops[j][0] == 0x21 and ops[j][2][0] == s and ops[j][2][1] not in tg:
                tg.append(ops[j][2][1])
                j += 1
            if len(tg) >= 4:
                out.append((i, j, s, tg))
                i = j
                continue
        i += 1
    return out


def rows(path, wpath, a, b, mode="x"):
    c, cum = T.cumulative(path)
    w = np.load(wpath)
    fr = frames(c, a, b) if mode == "x" else frames_cx(c, a, b)
    res = []
    k = 0
    prev_end = a
    while k + 1 < len(fr):
        f0, f1 = fr[k], fr[k + 1]
        if f0[2] != f1[2] or set(f0[3]) != set(f1[3]):
            k += 1
            continue
        res.append(dict(start=f0[0], end=f1[1], sign=f0[2], n=len(f0[3]), tof=cum[f1[1]] - cum[f0[0]],
                        pre=cum[f0[0]] - cum[prev_end], alloc=int(w[f0[0]:f0[1]].min()),
                        peak=int(w[f0[0]:f1[1]].max())))
        prev_end = f1[1]
        k += 2
    return res, cum, prev_end


if __name__ == "__main__":
    a, b = int(sys.argv[3]), int(sys.argv[4])
    res, cum, pe = rows(sys.argv[1], sys.argv[2], a, b, sys.argv[5] if len(sys.argv) > 5 else "x")
    if len(sys.argv) > 6:
        import json
        json.dump([dict(n=r["n"], alloc=r["alloc"], tof=r["tof"], feed=bool(r["pre"] >= 1.0 and r["n"] == 13))
                   for r in res], open(sys.argv[6], "w"))
    print(f"{'row':>3s} {'ops':>15s} {'sign':>4s} {'n':>2s} {'tof':>7s} {'pre':>6s} {'alloc':>5s} {'peak':>4s}")
    for i, r in enumerate(res):
        print(f"{i:3d} {r['start']:7d}-{r['end']:7d} q{r['sign']:<3d} {r['n']:2d} {r['tof']:7.2f} {r['pre']:6.2f} "
              f"{r['alloc']:5d} {r['peak']:4d}")
    print("rows", len(res), "row Toffolis", sum(r["tof"] for r in res), "between rows", sum(r["pre"] for r in res),
          "tail", cum[b] - cum[pe], "stage", cum[b] - cum[a])
