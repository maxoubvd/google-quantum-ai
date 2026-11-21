#!/usr/bin/env python3
"""dqi_dq3_kpn.py -- DQ3 (Nov 18): ROW 'kpn' of the census: the reverse partial path's 12-bit inverse that keeps its
norm (exp/frontier/root/decode_local_norm_partial.inverse12_keep_partial_norm, read only) in the kit's CACHED form.

At keep = 6 (the line's 3 calls a decode, unextract, entry 852, +24 [M census cen3]) it computes D = Norm(x) (15),
u = 1/D by its GF(64) leaf (inverse6_init19 as bound in its scope: 17 on the line), t = conj(x) u by two 6-bit
products (30): 62, and returns (t, u, D) for erase_norm_with_partial.  Here: tools/dqi_dq3_npair.cached_pair (the kit's
PLAN12 raw cache, the SAME leaf on D, D copied out, the 25-product final stage, u kept): 57, peak <= +24.  keep != 6:
the original.  Exactness: `python dqi_dq3_kpn.py test` (t, u, D equal on all 4,096 inputs).
Install: token dqi_dq3_kpn LAST in DQI_INT_EXT (needs tools/dqi_dq3_npair.py beside it).  Switches: DQI_DQ3_KPN_CEIL
(901), DQI_DQ3_KPN_EXTRA (24), DQI_DQ3_KPN_OFF=1, DQI_DQ3_KPN_OUT=<json>.
"""
import atexit
import gc
import importlib
import json
import os
import sys
import types

CEIL = int(os.environ.get('DQI_DQ3_KPN_CEIL', '901'))
EXTRA = int(os.environ.get('DQI_DQ3_KPN_EXTRA', '24'))
STATS = dict(cached=0, orig=0, patched=0)
CALLS = []
sys.modules.setdefault('dqi_dq3_kpn', sys.modules[__name__])


def cached(b, x, keep, G):
    npair = sys.modules.get('dqi_dq3_npair') or importlib.import_module('dqi_dq3_npair')
    jet = sys.modules.get('exp.jet24.pro_jet') or importlib.import_module('exp.jet24.pro_jet')
    leaf = G['inverse6_init19']
    t, D, u = npair.cached_pair(b, x, jet, lambda b_, D_: leaf(b_, D_))
    return t, u, D


_STUB = '''
def inverse12_keep_partial_norm(b, x, keep):
    _m = __import__('sys').modules['dqi_dq3_kpn']
    if keep == 6 and len(x) == 12 and b.num_allocated_qubits + _m.EXTRA <= _m.CEIL:
        _m.STATS['cached'] += 1
        _m.CALLS.append((b.num_allocated_qubits, 1))
        return _m.cached(b, x, keep, globals())
    _m.STATS['orig'] += 1
    _m.CALLS.append((b.num_allocated_qubits, 0))
    return _dq3_kpn_orig(b, x, keep)
'''


def _target(fn):
    co = getattr(fn, '__code__', None)
    return (isinstance(fn, types.FunctionType) and co is not None and co.co_name == 'inverse12_keep_partial_norm'
            and co.co_filename.endswith('root/decode_local_norm_partial.py') and 'dqi_dq3_kpn' not in co.co_consts)


def _patch(fn):
    gl = fn.__globals__
    if '_dq3_kpn_orig' not in gl:
        gl['_dq3_kpn_orig'] = types.FunctionType(fn.__code__, gl, 'inverse12_keep_partial_norm_line', fn.__defaults__,
                                                 None)
    ns = {}
    pad = '\n' * (fn.__code__.co_firstlineno - 1)
    exec(compile(pad + _STUB, fn.__code__.co_filename, 'exec'), gl, ns)
    fn.__code__ = ns['inverse12_keep_partial_norm'].__code__
    STATS['patched'] += 1


def _dump():
    path = os.environ.get('DQI_DQ3_KPN_OUT')
    if path:
        try:
            json.dump(dict(STATS, ceil=CEIL, extra=EXTRA, calls=CALLS), open(path, 'w'))
        except Exception:
            pass
    sys.stderr.write('dqi_dq3_kpn: %s\n' % json.dumps(dict(STATS, ceil=CEIL, extra=EXTRA)))


def install(km, ns):
    if os.environ.get('DQI_DQ3_KPN_OFF', '0') == '1':
        return dict(off=True)
    for obj in gc.get_objects():
        if _target(obj):
            _patch(obj)
    atexit.register(_dump)
    return dict(STATS)


def test(root):
    import importlib.util
    import numpy as np
    sys.path[:0] = [os.path.expanduser('~/qcc/data/kickmix'), root, os.path.join(root, 'playbook/kit/decoder'),
                    os.path.join(root, 'exp/r-main/narrow240'), os.path.join(root, 'exp/round7/root'),
                    os.path.join(root, 'exp/next/carry'), os.path.join(root, 'exp/frontier/root'),
                    os.path.join(root, 'exp/frontier/a03'), os.path.dirname(os.path.abspath(__file__))]
    import kickmix as km
    sys.modules.setdefault('circ', km)
    import decode_local_forward_sheared  # noqa: F401  (the kit's import paths)
    spec = importlib.util.spec_from_file_location('kpn_root', os.path.join(root, 'exp/frontier/root/decode_local_norm_partial.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    import dqi_dq3_lite as L

    def run(fn):
        b = L.LB(4096)
        x = b.alloc_qubits(12)
        for i in range(12):
            x[i].v[:] = (np.arange(4096) >> i) & 1
        t, u, D = fn(b, x, 6)
        val = lambda ws: sum(w.v.astype(np.int64) << i for i, w in enumerate(ws))
        ok = all((x[i].v == ((np.arange(4096) >> i) & 1).astype(bool)).all() for i in range(12))
        return val(t), val(u), val(D), ok, b.toff, b.peak - 12
    a = run(mod.inverse12_keep_partial_norm)
    c = run(lambda b, x, keep: cached(b, x, keep, mod.__dict__))
    print('line: Toffolis %d, peak extra %d; cached: Toffolis %d, peak extra %d' % (a[4], a[5], c[4], c[5]))
    print('t equal %s, u equal %s, D equal %s, x restored %s %s (all 4096 inputs)' % (
        bool((a[0] == c[0]).all()), bool((a[1] == c[1]).all()), bool((a[2] == c[2]).all()), a[3], c[3]))


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'test':
    test(sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser(
        '~/qcc/ws/dqi_dq3/gen/zc2p/home/qcc/ws/dqi_zn/src/kv_unfused5'))
