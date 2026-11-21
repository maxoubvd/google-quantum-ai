#!/usr/bin/env python3
"""dqi_int_inv6trace.py -- B5 (22nd holder, Nov 18): classical bit-parallel trace of the kit's GF(64) inversion leaves.

A tracing builder (the subset of the kickmix builder API the kit's leaves call) holds every wire as a bool vector over
ALL 64 inputs at once.  Each Toffoli-class op (ccx, init_and, del_and, ccz / cswap on quantum controls) is recorded
with its two control vectors.  Per AND: does a OR b = 1 on every NONZERO input (the AND is then 1 + a + b there:
'redundant except at 0'), is it constant, does it equal an affine form of the live wires.
Usage (WSL): ~/qcc/.venv/bin/python dqi_int_inv6trace.py [leaf=inverse6_joint18_12extra] [module=exp.jet24.pro_jet]
Reads the kit copy (read-only); writes nothing.
"""
import importlib
import os
import sys
from contextlib import contextmanager

import numpy as np

KIT = os.environ.get('DQI_JUNC_KIT', os.path.expanduser('~/qcc/ws/dqi_junc/kit'))
sys.path[:0] = [os.path.expanduser('~/qcc/data/kickmix'), KIT, os.path.dirname(KIT)]
import kickmix as km  # noqa: E402

sys.modules['circ'] = km
from exp.decoder import common as _C  # noqa: E402,F401  (sets up the kit's 'field' module as the decoder does)

N = int(os.environ.get('DQI_INT_TRACE_N', '64'))


class W:
    """one wire (qubit or classical bit)"""
    __slots__ = ('v', 'name', 'kind')

    def __init__(self, name, kind):
        self.v = np.zeros(N, dtype=bool)
        self.name = name
        self.kind = kind


class TB:
    """tracing builder: values over the 64 inputs; records Toffoli-class ops"""

    def __init__(self):
        self.nw = 0
        self.ops = []          # (kind, a_vec, b_vec, target name, live snapshot index)
        self.live = []
        self.rng = np.random.default_rng(1)

    def _w(self, kind):
        w = W('%s%d' % (kind, self.nw), kind)
        self.nw += 1
        self.live.append(w)
        return w

    def alloc_qubits(self, n):
        return [self._w('q') for _ in range(n)]

    def alloc_bits(self, n):
        return [self._w('c') for _ in range(n)]

    def free(self, ws):
        if isinstance(ws, W):
            ws = [ws]
        for w in ws:
            assert not w.v.any() or w.kind == 'c', 'freeing a dirty qubit %s' % w.name
            if w in self.live:
                self.live.remove(w)

    def bit_store(self, ws, val):
        for w in (ws if isinstance(ws, (list, tuple)) else [ws]):
            w.v[:] = bool(val)

    @contextmanager
    def mark(self, *a, **k):
        yield

    # Clifford-free classical gates
    def x(self, t):
        for w in (t if isinstance(t, (list, tuple)) else [t]):
            w.v ^= True

    def cx(self, a, t):
        if isinstance(a, (list, tuple)):
            for p, q in zip(a, t):
                q.v ^= p.v
        else:
            t.v ^= a.v

    def swap(self, a, t):
        if isinstance(a, (list, tuple)):
            for p, q in zip(a, t):
                p.v, q.v = q.v.copy(), p.v.copy()
        else:
            a.v, t.v = t.v.copy(), a.v.copy()

    def _rec(self, kind, a, b, t):
        quantum = (a.kind == 'q') and (b.kind == 'q')
        self.ops.append((kind, a.v.copy(), b.v.copy(), t.name if t is not None else None, quantum,
                         [w for w in self.live if w.kind == 'q'], [w.v.copy() for w in self.live if w.kind == 'q']))

    def ccx(self, a, b, t):
        self._rec('ccx', a, b, t)
        t.v ^= a.v & b.v

    def init_and(self, a, b, t):
        assert not t.v.any()
        self._rec('init_and', a, b, t)
        t.v ^= a.v & b.v

    def del_and(self, a, b, t):
        t.v ^= a.v & b.v
        assert not t.v.any()

    # phase gates: values unchanged
    def z(self, *a):
        pass

    def s(self, *a):
        pass

    def cz(self, *a):
        pass

    def ccz(self, a, b, c):
        if a.kind == 'q' and b.kind == 'q' and c.kind == 'q':
            self._rec('ccz', a, b, c)

    def hmr(self, q, m):
        if isinstance(q, (list, tuple)):
            for p, r in zip(q, m):
                self.hmr(p, r)
            return
        q.v[:] = False
        m.v[:] = self.rng.integers(0, 2, size=N).astype(bool)


def main():
    opts = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
    modname = opts.get('module', 'exp.jet24.pro_jet')
    leaf = opts.get('leaf', 'inverse6_joint18_12extra')
    mod = importlib.import_module(modname)
    fn = getattr(mod, leaf)
    gf = km.GF2Field(6)
    b = TB()
    x = b.alloc_qubits(6)
    for i in range(6):
        x[i].v[:] = (np.arange(N) >> i) & 1
    xin = [w.v.copy() for w in x]
    res = fn(b, x)
    if isinstance(res, tuple):
        out, F = res
    else:
        out, F = res, None
    val = sum(out[i].v.astype(np.int64) << i for i in range(6))
    want = np.array([0] + [gf.pow(v, 62) for v in range(1, 64)])
    print('leaf %s.%s: inverse exact on %d / 64 inputs; x restored %s' % (
        modname, leaf, int((val == want).sum()), all((w.v == v).all() for w, v in zip(x, xin))))
    if F is not None:
        print('  F = [x == 0] exact: %s' % bool((F.v == (np.arange(N) == 0)).all()))
    nz = np.arange(N) != 0
    tof = [o for o in b.ops if o[4]]
    print('  Toffoli-class ops with quantum controls: %d (all ops recorded %d)' % (len(tof), len(b.ops)))
    if 'dump' in opts:            # dump=<json>,<AND index>: live qubit values just before that AND + the targets
        import json
        path, at = opts['dump'].split(',')
        lv = [v for v in tof[int(at)][6] if v.any()]
        tobits = lambda v: int(sum(int(bool(v[i])) << i for i in range(N)))
        tgt = [tobits(((want >> i) & 1).astype(bool)) for i in range(6)] + [tobits(np.arange(N) == 0)]
        prods = [tobits(o[1] & o[2]) for o in tof]
        lefts = [tobits(o[1]) for o in tof]
        rights = [tobits(o[2]) for o in tof]
        json.dump(dict(leaf=leaf, before_and=int(at), ands_before=int(at), basis=[tobits(v) for v in lv], products=prods, lefts=lefts, rights=rights,
                       targets=tgt, target_names=['inv%d' % i for i in range(6)] + ['F']), open(path, 'w'))
        print('  dumped %s: %d live nonzero wires before AND %s' % (path, len(lv), at))
    for k, (kind, a, c, t, q, lw, lv) in enumerate(tof):
        prod = a & c
        orall = bool((a | c)[nz].all())
        const = 'const%d' % int(prod[0]) if (prod == prod[0]).all() else ''
        # affine in the live qubits (on all inputs / on nonzero inputs)?
        M = np.array([np.ones(N, bool)] + lv, dtype=np.uint8)
        def in_span(vec, mask):
            A = M[:, mask].T.astype(np.uint8)
            y = vec[mask].astype(np.uint8)
            return rank2(np.concatenate([A, y[:, None]], 1)) == rank2(A)
        aff_all = in_span(prod, np.ones(N, bool))
        aff_nz = in_span(prod, nz)
        print('  AND %2d %-8s a1 %2d b1 %2d ab %2d | a|b=1 on nonzero: %-5s %s affine(all) %s affine(nonzero) %s' % (
            k, kind, int(a.sum()), int(c.sum()), int(prod.sum()), orall, const, aff_all, aff_nz))


def rank2(M):
    M = (np.array(M, dtype=np.uint8) & 1).copy()
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        nz = np.flatnonzero(M[r:, c])
        if len(nz) == 0:
            continue
        p = r + int(nz[0])
        if p != r:
            M[[r, p]] = M[[p, r]]
        m = M[:, c].astype(bool)
        m[r] = False
        M[m] ^= M[r]
        r += 1
        if r == rows:
            break
    return r


if __name__ == '__main__':
    main()
