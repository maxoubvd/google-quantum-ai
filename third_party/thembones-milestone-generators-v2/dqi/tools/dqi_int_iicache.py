#!/usr/bin/env python3
"""dqi_int_iicache.py -- B5 (22nd holder, Nov 18): WIDTH-ADAPTIVE initialized 12-bit inverse.

The kit's inverse_search/inverse20.inverse_initialized(b, x, target=None) at 12 bits (the converter's restore_root /
post sites: Norm 15 + the 6-bit initialized leaf 19 + two 6-bit products 30 = 64) is replaced, where x has 12 qubits,
no target is given and the call starts with live + EXTRA <= CEIL, by the kit's flag-free cached 12-bit inverse
(decode_flagfree_cached.CACHED12[0]: 58, the same zero-correct 1/x returned in fresh qubits, x kept) -- the same choice
the kit's own inverse_initialized_flagfree12 makes below its width line.  -6 a call.  The recursive 6-bit and 3-bit
calls and every call with a target keep the original.
Install: token dqi_int_iicache LAST in DQI_INT_EXT (B5_FILES_ADD=dqi_int_iicache.py).  Switches: DQI_INT_II_CEIL (901),
DQI_INT_II_EXTRA (22), DQI_INT_II_OFF=1, DQI_INT_II_OUT=<json> (per call: live at entry, form).
"""
import atexit
import importlib
import json
import os
import sys
import types

CEIL = int(os.environ.get('DQI_INT_II_CEIL', '901'))
EXTRA = int(os.environ.get('DQI_INT_II_EXTRA', '22'))
STATS = dict(cached=0, orig=0, patched=0, found=0)
CALLS = []
_KIT = {}

_SRC = '''
def inverse_initialized(b, x, target=None):
    if target is None and len(x) == 12:
        live = b.num_allocated_qubits
        if live + _b5_ii_extra <= _b5_ii_ceil:
            _b5_ii_stats['cached'] += 1
            _b5_ii_calls.append((live, 1))
            out, fl = _b5_ii_cached()(b, x, _b5_ii_plan(), False)
            return out
        _b5_ii_stats['orig'] += 1
        _b5_ii_calls.append((live, 0))
    return _b5_ii_orig(b, x, target)
'''


def _cached():
    if 'c' not in _KIT:
        ff = sys.modules.get('decode_flagfree_cached') or importlib.import_module('decode_flagfree_cached')
        _KIT['c'] = ff.CACHED12[0]
    return _KIT['c']


def _plan():
    if 'p' not in _KIT:
        pj = sys.modules.get('exp.jet24.pro_jet') or importlib.import_module('exp.jet24.pro_jet')
        _KIT['p'] = pj.PLAN12
    return _KIT['p']


def _is_target(fn):
    co = getattr(fn, '__code__', None)
    return (isinstance(fn, types.FunctionType) and co is not None and co.co_name == 'inverse_initialized'
            and co.co_filename.endswith('inverse20.py') and not getattr(fn, '_b5_ii', False)
            and not fn.__closure__ and co.co_argcount == 3)


def _patch(fn):
    STATS['found'] += 1
    gl = fn.__globals__
    gl['_b5_ii_orig'] = types.FunctionType(fn.__code__, gl, 'inverse_initialized_kit', fn.__defaults__, None)
    gl['_b5_ii_cached'] = _cached
    gl['_b5_ii_plan'] = _plan
    gl['_b5_ii_ceil'] = CEIL
    gl['_b5_ii_extra'] = EXTRA
    gl['_b5_ii_stats'] = STATS
    gl['_b5_ii_calls'] = CALLS
    ns = {}
    pad = '\n' * (fn.__code__.co_firstlineno - 1)
    exec(compile(pad + _SRC, fn.__code__.co_filename, 'exec'), gl, ns)
    fn.__code__ = ns['inverse_initialized'].__code__
    fn._b5_ii = True
    STATS['patched'] += 1


def _cells(fn):
    out = []
    for cl in (getattr(fn, '__closure__', None) or ()):
        try:
            out.append(cl.cell_contents)
        except Exception:
            pass
    return out


def _dump():
    path = os.environ.get('DQI_INT_II_OUT')
    if path:
        try:
            json.dump(dict(STATS, ceil=CEIL, extra=EXTRA, calls=CALLS), open(path, 'w'))
        except Exception:
            pass
    sys.stderr.write('dqi_int_iicache: %s\n' % json.dumps(dict(STATS, ceil=CEIL, extra=EXTRA)))


def install(km, ns):
    if os.environ.get('DQI_INT_II_OFF', '0') == '1':
        return dict(STATS, off=True)
    walk = sys.modules['walk_bound']
    seen = set()

    def visit(v):
        if id(v) in seen:
            return
        seen.add(id(v))
        if _is_target(v):
            _patch(v)
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
