#!/usr/bin/env python3
"""dqi_int_zcache.py -- B5 (22nd holder, Nov 18): WIDTH-ADAPTIVE 12-bit pivot inverse.

The line's 12-bit pivot inverses (convert, unconvert and the other callers of exact_subs.inverse_zero) run the tower
form Norm 15 + GF(64) leaf 18 + two 6-bit products 30 = 63 Toffolis at a peak of 19 qubits beyond the entry.  The kit's
cached 12-bit inverse cached_inverse(x, PLAN12, with_zero=True) computes the same (x^-1, F = [x == 0]) in 58 Toffolis
(norm products kept as raw wires, final stage 25 instead of 30) at a peak of 22 [M tools/dqi_int_zcwidth.py: both exact
on all 4,096 inputs, x restored, out + F the only live wires after].  Where a call starts with live + 22 <= CEIL the
cached form runs, elsewhere the original: -5 a call where it fits, nothing else changes in the callers (same
interface (out, F); the callers' later erase / consume read only out and F).
Install: token dqi_int_zcache LAST in DQI_INT_EXT (B5_FILES_ADD=dqi_int_zcache.py).  Switches: DQI_INT_ZC_CEIL (899),
DQI_INT_ZC_EXTRA (22), DQI_INT_ZC_OFF=1.  DQI_INT_ZC_OUT=<json>: per call (caller, live at entry, form) at exit.
"""
import atexit
import importlib
import json
import os
import sys
import types

CEIL = int(os.environ.get('DQI_INT_ZC_CEIL', '899'))
EXTRA = int(os.environ.get('DQI_INT_ZC_EXTRA', '22'))
STATS = dict(cached=0, orig=0, patched=0, found=0)
CALLS = []

_SRC = '''
def inverse_zero(b, x):
    live = b.num_allocated_qubits
    if live + _b5_zc_extra <= _b5_zc_ceil:
        _b5_zc_stats['cached'] += 1
        _b5_zc_calls.append((live, 1))
        return _b5_zc_cached(b, x, _b5_zc_plan, True)
    _b5_zc_stats['orig'] += 1
    _b5_zc_calls.append((live, 0))
    return _b5_zc_orig(b, x)
'''


def _is_target(fn):
    co = getattr(fn, '__code__', None)
    return (isinstance(fn, types.FunctionType) and co is not None and co.co_name == 'inverse_zero'
            and co.co_filename.endswith('exact_subs.py') and not getattr(fn, '_b5_zc', False)
            and not fn.__closure__ and fn.__code__.co_argcount == 2)


def _patch(fn, pj):
    STATS['found'] += 1
    gl = fn.__globals__
    gl['_b5_zc_orig'] = types.FunctionType(fn.__code__, gl, 'inverse_zero_line', fn.__defaults__, None)
    gl['_b5_zc_cached'] = pj.cached_inverse
    gl['_b5_zc_plan'] = pj.PLAN12
    gl['_b5_zc_ceil'] = CEIL
    gl['_b5_zc_extra'] = EXTRA
    gl['_b5_zc_stats'] = STATS
    gl['_b5_zc_calls'] = CALLS
    ns = {}
    pad = '\n' * (fn.__code__.co_firstlineno - 1)
    exec(compile(pad + _SRC, fn.__code__.co_filename, 'exec'), gl, ns)
    fn.__code__ = ns['inverse_zero'].__code__
    fn._b5_zc = True
    STATS['patched'] += 1


def _cells(fn):
    out = []
    if not isinstance(fn, types.FunctionType):
        return out
    for cl in (fn.__closure__ or ()):
        try:
            out.append(cl.cell_contents)
        except Exception:
            pass
    return out


def _dump():
    path = os.environ.get('DQI_INT_ZC_OUT')
    rec = dict(STATS, ceil=CEIL, extra=EXTRA, calls=CALLS)
    if path:
        try:
            json.dump(rec, open(path, 'w'))
        except Exception:
            pass
    sys.stderr.write('dqi_int_zcache: %s\n' % json.dumps(dict(STATS, ceil=CEIL, extra=EXTRA)))


def install(km, ns):
    if os.environ.get('DQI_INT_ZC_OFF', '0') == '1':
        return dict(STATS, off=True)
    pj = sys.modules.get('exp.jet24.pro_jet') or importlib.import_module('exp.jet24.pro_jet')
    walk = sys.modules['walk_bound']
    seen = set()

    def visit(v):
        if id(v) in seen:
            return
        seen.add(id(v))
        if _is_target(v):
            _patch(v, pj)
    for fn in walk.functions(ns.decode):
        if not isinstance(fn, types.FunctionType):
            continue
        visit(fn)
        cells = _cells(fn)
        ds = [fn.__globals__] + [c for c in cells if isinstance(c, dict)]
        for c in cells:
            if isinstance(c, types.FunctionType):
                visit(c)
        for d in ds:
            try:
                vals = list(d.values())
            except Exception:
                continue
            for v in vals:
                if isinstance(v, types.FunctionType):
                    visit(v)
                    for c in _cells(v):
                        if isinstance(c, types.FunctionType):
                            visit(c)
    atexit.register(_dump)
    return dict(found=STATS['found'], patched=STATS['patched'], ceil=CEIL, extra=EXTRA)
