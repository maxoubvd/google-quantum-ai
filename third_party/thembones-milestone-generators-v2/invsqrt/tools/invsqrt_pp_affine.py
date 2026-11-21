#!/usr/bin/env python3
"""invsqrt_pp_affine.py -- workstream invsqrt_pp (P1): EXACT, block-aware affine audit of an invsqrt .kmb.

For every Toffoli-class op, on all 3,133,452 inputs: is its product (CCX a*b, CSWAP c*(x+y)) in the affine span of
the 46 wire values just before it (then the op costs 0 Toffolis: CX network, or HMR + CZ when the target is in the
relation); for CCZ, is a*b*c a quadratic form of the wires (then CZ / Z network). Ops inside conditioned blocks are
evaluated with the block EXECUTED (the only case in which they run) for every inner-outcome pattern of 0 and 1
(both must agree). Prints each site with its expected weight.
Usage: python invsqrt_pp_affine.py <circuit.kmb> [out.json]
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import invsqrt_pp_sim as sim  # noqa: E402
from invsqrt_pp_blocks import gf2_solve, quad_basis, sample_bits, vec_sample  # noqa: E402
from invsqrt_pp_exlinks import exact_reduce  # noqa: E402

TOF = (0x30, 0x31, 0x32)


def product(C, i, S):
    q = C.qs[i]
    b = C.base[i]
    if b == 0x32:
        return S[q[0]] & (S[q[1]] ^ S[q[2]])
    if b == 0x30:
        return S[q[0]] & S[q[1]] & S[q[2]]
    return S[q[0]] & S[q[1]]


def main():
    t0 = time.time()
    C = sim.load(sys.argv[1])
    outp = sys.argv[2] if len(sys.argv) > 2 else None
    rng = np.random.default_rng(11)
    lanes = np.sort(rng.choice(sim.NX, size=2048, replace=False)).astype(np.int64)
    sites = []

    def test(i, S, ctx):
        b = C.base[i]
        p = product(C, i, S)
        w = 0.5 ** (C.level[i] + (1 if C.cnd[i] >= 0 else 0))
        if b in (0x31, 0x32):
            R = exact_reduce(S, p[None, :])
            if not R[0].any():
                sites.append({'op': i, 'kind': kmb.NAMES[b], 'class': 'affine', 'w': w, 'ctx': ctx})
        else:
            # CCZ: quadratic span on samples, then exact check of the solution
            A, cols = quad_basis(sample_bits(S, lanes).T)
            ok, sol = gf2_solve(A, vec_sample(p, lanes)[:, None])
            if ok[0]:
                s = sol[0]
                nq = C.nq
                iu, ju = np.triu_indices(nq, 1)
                q = np.zeros(sim.W, np.uint64)
                if s[0]:
                    q ^= sim.ONES
                for ww in range(nq):
                    if s[1 + ww]:
                        q ^= S[ww]
                for t in np.flatnonzero(s[1 + nq:]).tolist():
                    q ^= S[iu[t]] & S[ju[t]]
                d = q ^ p
                d[-1] &= np.uint64((1 << (sim.NX - 64 * (sim.W - 1))) - 1)
                if not d.any():
                    sites.append({'op': i, 'kind': 'CCZ', 'class': 'quadratic', 'w': w, 'ctx': ctx,
                                  'terms': int(s.sum())})

    def hook(i, S, P, bits, act):
        b = C.base[i]
        if C.level[i] == 0 and b in TOF and C.cnd[i] < 0:
            test(i, S, 'top')
        if C.level[i] == 0 and b == kmb.PUSH:
            hi = C.match[i]
            if not any(C.base[k] in TOF for k in range(i + 1, hi)):
                return
            # run the block as executed, testing every Toffoli-class op inside; inner outcomes 0 then 1
            res = {}
            for oc in (0, 1):
                S2 = S.copy()
                b2 = list(bits)
                found = []

                def inner(j, SS, PP, bb, aa, found=found):
                    if aa and C.base[j] in TOF and (C.cnd[j] < 0 or bb[C.cnd[j]]):
                        n0 = len(sites)
                        test(j, SS, 'block %d oc %d' % (i, oc))
                        if len(sites) > n0:
                            found.append(sites.pop())
                sim.run(C, hook=inner, outcome=(lambda j, oc=oc: oc), S=S2, P=np.zeros(sim.W, np.uint64),
                        bits=b2, lo=i + 1, hi=hi, active_stack=[True])
                res[oc] = found
            ops0 = set(s['op'] for s in res[0])
            ops1 = set(s['op'] for s in res[1])
            for s in res[0]:
                s['both_outcomes'] = s['op'] in ops1
                sites.append(s)
            for s in res[1]:
                if s['op'] not in ops0:
                    s['both_outcomes'] = False
                    sites.append(s)
    sim.run(C, hook=hook)
    tot = sum(s['w'] for s in sites if s.get('both_outcomes', True))
    print('sites', len(sites), 'expected Toffolis (sites valid for every inner outcome tested)', tot,
          round(time.time() - t0, 1), 's')
    for s in sites:
        print('  ', s)
    if outp:
        json.dump(sites, open(outp, 'w'))


if __name__ == '__main__':
    main()
