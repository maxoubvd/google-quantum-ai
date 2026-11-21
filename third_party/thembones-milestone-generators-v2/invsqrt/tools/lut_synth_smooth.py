"""Find the smooth elements of a target space (few non-constant 16-address blocks) and cache them.
Usage: lut_synth_smooth.py [iters=1200]  -> ~/qcc/ws/lut_synth/smooth16.json (incumbent's 35-wire table)
Library: find_smooth(T, iters) -> list of {"bw", "tr", "blocks", "f"} sorted by block weight (independent set)."""
import json
import os
import random
import sys

sys.path.insert(0, "<HOME>/Quantum/tools")
from lut_synth_core import BASE, Span

BL16 = [sum(1 << a for a in range(b * 16, (b + 1) * 16)) for b in range(16)]


def nonconst(f):
    return [b for b in range(16) if (f & BL16[b]) not in (0, BL16[b])]


def transitions(tt):
    n = 0
    for half in (0, 1):
        for a in range(127):
            n += ((tt >> (half * 128 + a)) ^ (tt >> (half * 128 + a + 1))) & 1
    return n


def find_smooth(T, iters=1200, seed=7, keep=14):
    gens = list(T) + BASE[1:]
    nt = len(T)
    rng = random.Random(seed)
    best = {}
    for it in range(iters):
        zero = rng.sample(range(16), 2)
        mask = BL16[zero[0]] | BL16[zero[1]]
        piv = {}
        ker = []
        for i, t in enumerate(gens):
            v = t & mask
            for z in zero:
                if v & (BL16[z] & -BL16[z]):
                    v ^= BL16[z]
            c = 1 << i
            while v:
                h = v.bit_length() - 1
                if h in piv:
                    v ^= piv[h][0]
                    c ^= piv[h][1]
                else:
                    piv[h] = (v, c)
                    break
            if v == 0:
                ker.append(c)
        ker = [c for c in ker if c & ((1 << nt) - 1)]
        if not ker or len(ker) > 13:
            continue
        for m in range(1, 1 << len(ker)):
            c = 0
            for j in range(len(ker)):
                if (m >> j) & 1:
                    c ^= ker[j]
            key = c & ((1 << nt) - 1)
            if not key:
                continue
            f = 0
            for i, t in enumerate(gens):
                if (c >> i) & 1:
                    f ^= t
            w = len(nonconst(f))
            if key not in best or best[key][0] > w:
                best[key] = (w, f)
    sp = Span()
    for v in BASE:
        sp.add(v)
    out = []
    for w, f in sorted(best.values()):
        if sp.add(f)[0]:
            out.append({"bw": w, "tr": transitions(f), "blocks": nonconst(f), "f": hex(f)})
        if len(out) >= keep:
            break
    return out


if __name__ == "__main__":
    from lut_synth_inv_a import load_targets
    T, wires = load_targets()
    iters = int(sys.argv[1]) if len(sys.argv) > 1 else 1200
    out = find_smooth(T, iters)
    json.dump(out, open(os.path.expanduser("~/qcc/ws/lut_synth/smooth16.json"), "w"))
    for o in out:
        print(o["bw"], o["tr"], o["blocks"])
