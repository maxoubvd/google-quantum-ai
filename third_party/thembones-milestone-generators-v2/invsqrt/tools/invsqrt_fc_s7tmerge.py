#!/usr/bin/env python3
"""invsqrt_fc_s7tmerge.py -- FC (seventh holder): FOLD T-TYPE Z_POWs INTO GENERIC ROTATIONS WITH THE SAME FUNCTION.

The invsqrt circuits are classical reversible + diagonal phases (measurement kicks are repaired), so the final phase is
sum_k theta_k f_k(x) over the Z_POW ops, f_k = the value of the op's wire there as a function of the input x. A Z_POW
whose angle is an odd multiple of 1/8 turn costs 1 T in the official verifier; a generic rotation costs
0.56 log2(1/eps) + 5.3 T whatever its angle (eps from the margin). So a T-type Z_POW whose f equals (or complements)
the f of a generic Z_POW is deleted and its angle added to (or subtracted from) that rotation: -1 T each, the phase
on every input unchanged (a complement only adds a global phase). Checked on all 3,133,452 inputs with P1's checker.
  python invsqrt_fc_s7tmerge.py IN.kmb OUT.kmb
"""
import hashlib
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_int_kmbw as kw  # noqa: E402
import invsqrt_pp_sim as sim  # noqa: E402

M64 = (1 << 64) - 1


def rot_phase(path, outcome=None):
    """Per input: sum over executed Z_POWs of angle (2^-64 turns) * wire value, mod 2^64; plus final S, P."""
    k = kmb.read(path)
    zidx = np.flatnonzero((np.asarray(k.ops) & 0x7F) == 0x13).tolist()
    w1 = {i: np.uint64(int(k.angles[n][1])) for n, i in enumerate(zidx)}
    C = sim.load(path)
    acc = np.zeros(sim.NX, np.uint64)

    def hook(i, S, P, bits, act):
        if i in w1 and act:
            v = np.unpackbits(S[C.qs[i][0]].view(np.uint8), bitorder='little')[:sim.NX].astype(np.uint64)
            acc[:] = acc + v * w1[i]            # uint64 arithmetic wraps mod 2^64
    S, P, _, _ = sim.run(C, hook=hook, outcome=outcome)
    return acc, S, P, C


def check(parent, child):
    """Exact on all inputs: final state and sign phase equal (zero + 3 random transcripts), and the total rotation
    phase equal on every input up to one global constant."""
    ok = True
    a0, S0, P0, C = rot_phase(parent)
    for seed in (0, 1, 2, 3):
        oc = None
        if seed:
            D = sim.load(child)
            dr = np.random.default_rng(seed).integers(0, 2, size=D.n)
            oc = (lambda j, d=dr: int(d[j]))
        a1, S1, P1, _ = rot_phase(child, oc)
        d = a1 - a0
        e = bool(np.array_equal(S1, S0)) and bool(np.array_equal(P1, P0)) and bool(np.all(d == d[0]))
        print('  check seed %d: state, sign phase and rotation phase (global offset %d) equal: %s' % (seed, int(d[0]), e),
              flush=True)
        ok &= e
    return ok


def main():
    src, dst = sys.argv[1], sys.argv[2]
    k = kmb.read(src)
    zidx = np.flatnonzero((np.asarray(k.ops) & 0x7F) == 0x13).tolist()
    assert len(zidx) == k.angles.shape[0], (len(zidx), k.angles.shape)
    row = {i: n for n, i in enumerate(zidx)}
    w0 = {i: int(k.angles[row[i]][0]) for i in zidx}
    w1 = {i: int(k.angles[row[i]][1]) for i in zidx}
    assert all(v == 0 for v in w0.values()), 'w0 != 0: angle layout not handled'
    C = sim.load(src)
    nx = sim.NX
    H = {}

    def hook(i, S, P, bits, act):
        if i in row and act and C.level[i] == 0 and C.cnd[i] < 0:
            v = np.unpackbits(S[C.qs[i][0]].view(np.uint8), bitorder='little')[:nx]
            H[i] = (hashlib.sha1(np.packbits(v).tobytes()).hexdigest(),
                    hashlib.sha1(np.packbits(1 - v).tobytes()).hexdigest())
    sim.run(C, hook=hook)

    def kind(i):
        a = w1[i]
        if a % (1 << 62) == 0:
            return 'clifford'
        if a % (1 << 61) == 0:
            return 't'
        return 'rot'
    gen = {}
    for i in zidx:
        if i in H and kind(i) == 'rot':
            gen.setdefault(H[i][0], i)
    delete, new_w1 = [], dict(w1)
    for i in zidx:
        if i not in H or kind(i) != 't':
            continue
        h, hc = H[i]
        if h in gen:
            g = gen[h]
            new_w1[g] = (new_w1[g] + w1[i]) & M64
        elif hc in gen:                       # f_T = 1 - f_g: theta (1 - f_g) = global theta - theta f_g
            g = gen[hc]
            new_w1[g] = (new_w1[g] - w1[i]) & M64
        else:
            continue
        delete.append(i)
        print('fold T-type Z_POW op %d (wire %d) into rotation op %d (wire %d)%s' % (
            i, C.qs[i][0], g, C.qs[g][0], '' if h in gen else ' (complement)'), flush=True)
    if not delete:
        print('no T-type Z_POW shares its function with a generic rotation: no child')
        return
    for g in zidx:                            # a fold must not turn a rotation into a T / Clifford angle
        if new_w1[g] != w1[g]:
            assert new_w1[g] % (1 << 61) != 0, 'fold made a T / Clifford angle'
    k2 = kw.edit(k, delete, {})
    keep_rows = [row[i] for i in zidx if i not in delete]
    ang = np.asarray(k.angles).copy()
    for i in zidx:
        ang[row[i]][1] = np.uint64(new_w1[i])
    k2.angles = ang[keep_rows]
    tmp = dst + '.tmp'
    kw.write(tmp, k2)
    ok = check(src, tmp)
    print('T-type Z_POWs folded: %d (T -%d per shot); CHECK against the input %s' % (
        len(delete), len(delete), 'PASS' if ok else 'FAIL'), flush=True)
    if ok:
        os.replace(tmp, dst)
    else:
        os.remove(tmp)


if __name__ == "__main__":
    main()
