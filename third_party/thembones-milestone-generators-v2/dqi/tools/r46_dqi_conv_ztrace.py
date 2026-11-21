"""r46_dqi_conv_ztrace.py -- R46-DQI-CONV: allocation trace of a decoder build (no gate change, no marker).

EXT module (append to DQI_INT_EXT): wraps kickmix.CircuitBuilder.alloc_qubits / free and the Toffoli-class emitters
(ccx, ccz, cswap, init_and with quantum operands: the same width rule as tools/dqi_build_patch_r7prof.py, so K = the
static Toffoli index of the built .kmb). For every builder it records (K, allocated qubits) after each alloc and
free; at exit it writes R46_ZTRACE_OUT (npz): for each builder with more than 100,000 static Toffolis the arrays
K_<i>, A_<i> and total_<i>. The real decode is the one whose total equals T_static of dqi.build.json.
"""
import atexit
import os
import sys

import numpy as np

OUT = os.environ.get('R46_ZTRACE_OUT', 'ztrace.npz')
PAID = ('ccx', 'ccz', 'cswap', 'init_and')


def install(km, ns):
    km = sys.modules['kickmix']
    B = km.CircuitBuilder
    if getattr(B, '_r46_ztrace', False):
        return {'ztrace': 'already'}
    B._r46_ztrace = True
    count = {}
    ev = {}

    def width(a, kw):
        n = 1
        for x in list(a) + list(kw.values()):
            try:
                k = len(x)
            except TypeError:
                if type(x).__name__ == 'b':
                    return 0
                continue
            if k and type(x[0]).__name__ == 'b':
                return 0
            n = max(n, k)
        return n

    def wrap(name):
        original = getattr(B, name)

        def paid(self, *a, **kw):
            n = width(a, kw)
            if n:
                count[id(self)] = count.get(id(self), 0) + n
            return original(self, *a, **kw)
        setattr(B, name, paid)

    for name in PAID:
        wrap(name)
    # a builder that dies may hand its id to a later one: reset the counters at construction; if the class does not
    # take a wrapped __init__, keep every builder alive instead (ids are then never reused)
    keep = []
    mode = 'init'
    try:
        o_init = B.__init__

        def init(self, *a, **kw):
            o_init(self, *a, **kw)
            count[id(self)] = 0
            ev[id(self)] = []
        B.__init__ = init
    except Exception:
        mode = 'keep'
    o_alloc, o_free = B.alloc_qubits, B.free

    def alloc_qubits(self, *a, **kw):
        out = o_alloc(self, *a, **kw)
        if mode == 'keep' and id(self) not in ev:
            keep.append(self)
        ev.setdefault(id(self), []).append((count.get(id(self), 0), self.num_allocated_qubits))
        return out

    def free(self, *a, **kw):
        out = o_free(self, *a, **kw)
        ev.setdefault(id(self), []).append((count.get(id(self), 0), self.num_allocated_qubits))
        return out

    B.alloc_qubits = alloc_qubits
    B.free = free

    def report():
        d = {}
        i = 0
        for bid, tot in count.items():
            if tot > 100000 and bid in ev:
                e = np.array(ev[bid], dtype=np.int32)
                d['K_%d' % i] = e[:, 0]
                d['A_%d' % i] = e[:, 1].astype(np.int16)
                d['total_%d' % i] = tot
                i += 1
        np.savez(OUT, **d)
        print('r46_ztrace mode', mode, i, 'builders ->', OUT, [int(d['total_%d' % j]) for j in range(i)], file=sys.stderr)

    atexit.register(report)
    return {'ztrace': OUT}
