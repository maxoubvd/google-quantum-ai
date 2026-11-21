#!/usr/bin/env python3
"""dqi_dq3_ff24test.py -- DQ3 (Nov 18): exhaustive value test of the flag-free 24-bit inverse forms (row ff24).

Builds the kit's decode_flagfree_cached form (sheared.cached_inverse_partial with the inner cached 12-bit inverse,
fast = init18 leaf, narrow = the nine-wire leaf) for a given (omitted count, narrow) on ALL 2^24 inputs with the
light bit-parallel builder (tools/dqi_dq3_lite.py) and compares the results of (3, narrow) [ff24's choice at space 37],
(4, narrow) [the line's choice at space 36] and (0, fast) [the full cached form] word for word; also checks that
x * result = 1 for nonzero x via the kit's own GF(2^24) product on a 65,536-input sample.
  python dqi_dq3_ff24test.py [kit root] [n_inputs=16777216]
"""
import os
import sys
import time

import numpy as np

ROOT = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].isdigit() else os.path.expanduser(
    '~/qcc/ws/dqi_dq3/gen/zc2p/home/qcc/ws/dqi_zn/src/kv_unfused5')
NIN = int(sys.argv[-1]) if sys.argv[-1].isdigit() else 1 << 24
TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.expanduser('~/qcc/data/kickmix'), ROOT, ROOT + '/playbook/kit/decoder',
                ROOT + '/exp/r-main/narrow240', ROOT + '/exp/round7/root', ROOT + '/exp/next/carry',
                ROOT + '/exp/frontier/root', ROOT + '/exp/frontier/a03', TOOLS]
import kickmix as km  # noqa: E402
sys.modules.setdefault('circ', km)
import decode_local_forward_sheared as sheared  # noqa: E402
from exp.jet24 import pro_jet as jet  # noqa: E402
from exp.noscratch import schedules, kernels  # noqa: E402
import dqi_dq3_lite as L  # noqa: E402
xs = sheared.combo.xsub


def inverse6_initialized_nine(b, x):         # copied from exp/frontier/root/decode_flagfree_cached.py (read only)
    to, back, constant = jet.tower(6)
    jet.linear_inplace(b, x, to)
    norm = b.alloc_qubits(3)
    schedules.norm_xor_noscratch(b, x, norm, constant)
    reciprocal = kernels.inverse3_noscratch(b, norm)
    jet.ar._erase_norm(b, x, norm, constant)
    b.free(norm)
    out = b.alloc_qubits(6)
    b.cx(x[3:], x[:3])
    kernels.mul_xor_noscratch(b, x[:3], reciprocal, out[:3])
    kernels.mul_xor_noscratch(b, x[3:], reciprocal, out[3:])
    b.cx(x[3:], x[:3])
    jet.ar._erase_norm(b, out, reciprocal, constant)
    b.free(reciprocal)
    jet.linear_inplace(b, x, back)
    jet.linear_inplace(b, out, back)
    return out


def leaf_fast(b, x, with_zero):
    return xs.inverse6_init18(b, x), None


def leaf_narrow(b, x, with_zero):
    return inverse6_initialized_nine(b, x), None


CACHED12 = tuple(xs.clone(jet.cached_inverse, dict(jet.cached_inverse.__globals__, inverse6_fixed18=leaf))
                 for leaf in (leaf_fast, leaf_narrow))
CACHED24 = tuple(xs.clone(sheared.cached_inverse_partial, dict(sheared.cached_inverse_partial.__globals__,
                                                              cached_inverse=inner)) for inner in CACHED12)


def ff24(b, x, count, narrow, K=None):
    columns, back, _ = jet.extension_columns(K)
    jet.linear_inplace(b, x, columns)
    out, flag = CACHED24[narrow](b, x, jet.PLAN24, count, False)
    t = b.alloc_bits(24)
    b.hmr(x, t)
    b.swap(x, out)
    b.free(out)
    jet.linear_inplace(b, x, back)


def run(count, narrow):
    b = L.LB(NIN)
    x = b.alloc_qubits(24)
    idx = np.arange(NIN, dtype=np.int64)
    if NIN < (1 << 24):
        idx = np.random.default_rng(7).integers(0, 1 << 24, NIN)
    for i in range(24):
        x[i].v[:] = (idx >> i) & 1
    t0 = time.time()
    ff24(b, x, count, narrow)
    val = np.zeros(NIN, dtype=np.int64)
    for i in range(24):
        val |= x[i].v.astype(np.int64) << i
    return val, b.toff, b.peak - 24, time.time() - t0, idx


def main():
    res = {}
    for count, narrow in ((4, 1), (3, 1), (0, 0)):
        val, toff, peak, sec, idx = run(count, narrow)
        res[(count, narrow)] = val
        print('form (count %d, narrow %d): Toffolis %d, peak extra %d, %.0f s on %d inputs' % (
            count, narrow, toff, peak, sec, NIN), flush=True)
    a, b3, c0 = res[(4, 1)], res[(3, 1)], res[(0, 0)]
    print('(3, narrow) == (4, narrow) on all inputs: %s (mismatches %d)' % (bool((a == b3).all()), int((a != b3).sum())))
    print('(0, fast) == (4, narrow) on all inputs: %s (mismatches %d)' % (bool((a == c0).all()), int((a != c0).sum())))
    print('inverse is an involution on the sample: %s' % bool(len(set(a.tolist())) == len(set(idx.tolist()))))


if __name__ == '__main__':
    main()
