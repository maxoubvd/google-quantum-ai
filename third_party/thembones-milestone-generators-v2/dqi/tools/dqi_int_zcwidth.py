#!/usr/bin/env python3
"""dqi_int_zcwidth.py -- B5 (Nov 18): Toffolis and peak extra width (qubits allocated beyond the entry) of the line's
12-bit pivot inverse exact_subs.inverse_zero against the kit's cached form cached_inverse(x, PLAN12, True), both
traced on all 4,096 inputs with the bit-parallel builder; checks both return x^-1 and F exactly."""
import os, sys
os.environ.setdefault('DQI_INT_TRACE_N', '4096')
import numpy as np
import dqi_int_inv6trace as T
SRC = os.path.expanduser(os.environ.get('B5_SRC', '~/qcc/ws/dqi_int/src/source_cap901'))
import importlib
import kickmix as km
gf = km.GF2Field(12)
want = np.array([0] + [gf.pow(v, 4094) for v in range(1, 4096)])


class WB(T.TB):
    def __init__(self):
        super().__init__()
        self.peak = 0
    def alloc_qubits(self, n):
        r = super().alloc_qubits(n)
        self.peak = max(self.peak, sum(1 for w in self.live if w.kind == 'q'))
        return r


def run(fn, label):
    b = WB()
    x = b.alloc_qubits(12)
    for i in range(12):
        x[i].v[:] = (np.arange(4096) >> i) & 1
    entry = sum(1 for w in b.live if w.kind == 'q')
    b.peak = entry
    out, F = fn(b, x)
    val = sum(out[i].v.astype(np.int64) << i for i in range(12))
    nt = sum(1 for o in b.ops if o[4])
    after = sum(1 for w in b.live if w.kind == 'q')
    print('%-34s Toffolis %d, peak extra %d, live after (x + out + F) - entry = %d, exact %s, F exact %s' % (
        label, nt, b.peak - entry, after - entry, bool((val == want).all()), bool((F.v == (np.arange(4096) == 0)).all())))


from exp.jet24 import pro_jet as PJ
run(lambda b, x: PJ.cached_inverse(b, x, PJ.PLAN12, True), 'kit cached_inverse PLAN12 (58)')
sys.path.insert(0, SRC)
ES = importlib.import_module('exp.r-main.narrow240.exact_subs'.replace('r-main', 'r-main')) if False else None
import importlib.util
spec = importlib.util.spec_from_file_location('b5_exact_subs', os.path.join(SRC, 'exp/r-main/narrow240/exact_subs.py'))
try:
    ES = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ES)
    run(ES.inverse_zero, 'line exact_subs.inverse_zero (63)')
except Exception as e:
    print('exact_subs import failed:', type(e).__name__, str(e)[:200])
