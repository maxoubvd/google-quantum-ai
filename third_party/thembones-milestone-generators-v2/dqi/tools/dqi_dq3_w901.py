#!/usr/bin/env python3
"""dqi_dq3_w901.py -- DQ3 (Nov 18): ROW 'w901' of the width-slack census: the kit's width gates re-dispatched to OUR
ceiling.

The kit chooses kernel forms by comparing the live count with a ceiling literal: 889 (the public decoder's width) or
900; its own dispatcher (exp/next/a01/dispatch.py, called with ceiling 900 by the line's decode root) rewrites both
literals to 900 in every function it reaches at import.  The line's file is 901 wide, so every gate binds ONE WIRE
EARLY, and a gate the dispatcher did not reach still binds at 889.  Example [M census cen2 + the kit's own selector]:
the extraction's flag-free 24-bit inverse (decode_flagfree_cached.choose_outer(ceiling - live)) enters at 864 live:
space 36 -> 183 Toffolis; space 37 -> 179.
This module does what the dispatcher does, at CEIL (default 901): for every function object (gc scan) whose code (or
nested code) reads `num_allocated_qubits` and holds an int constant 889 or 900, the constant becomes CEIL; closure
cells and defaults named 'ceiling' holding 889 / 900 likewise.  Only code from the kit tree is touched (file path
contains /qcc/ws/dqi_zn/src/ or /research/agents/ or /playbook/ or /exp/), never our tools.  Each form a gate selects
is exact on all inputs on its own and is already used in the line's file at other widths: only the choice moves.
Install: token dqi_dq3_w901 LAST in DQI_INT_EXT.  Switches: DQI_DQ3_W901_CEIL (901), DQI_DQ3_W901_OFF=1,
DQI_DQ3_W901_OUT=<json> (every code object rewritten: name, file tail, first line).
"""
import atexit
import gc
import json
import os
import sys
import types

CEIL = int(os.environ.get('DQI_DQ3_W901_CEIL', '901'))
OLD = (889, 900)
# DQI_DQ3_W901_NAMES=a,b: only code objects with these names (default: every gate). Measured 06:00Z: every gate at
# once reads +13.8 at the marker (w901a), so the row is applied per gate.
NAMES = set(x for x in os.environ.get('DQI_DQ3_W901_NAMES', '').split(',') if x)
KIT = ('/qcc/ws/dqi_zn/src/', '/research/agents/', '/playbook/', '/exp/')
STATS = dict(codes=0, functions=0, cells=0, defaults=0)
PATCHED = []
_NEW = {}


def _kit(code):
    fn = code.co_filename
    return any(k in fn for k in KIT) and '/Quantum/tools/' not in fn


def _gate(code):
    return ('num_allocated_qubits' in code.co_names and _kit(code) and (not NAMES or code.co_name in NAMES)
            and any(type(c) is int and c in OLD for c in code.co_consts))


def _has(code):
    return _gate(code) or any(isinstance(c, types.CodeType) and _has(c) for c in code.co_consts)


def _rebuild(code):
    key = id(code)
    if key in _NEW:
        return _NEW[key][1]
    gate = _gate(code)
    consts = []
    for c in code.co_consts:
        if isinstance(c, types.CodeType) and _has(c):
            consts.append(_rebuild(c))
        elif gate and type(c) is int and c in OLD:
            consts.append(CEIL)
        else:
            consts.append(c)
    new = code.replace(co_consts=tuple(consts))
    _NEW[key] = (code, new)          # keep the old code alive so its id is not reused
    if gate:
        STATS['codes'] += 1
        PATCHED.append(dict(name=code.co_name, file=code.co_filename[-70:], line=code.co_firstlineno))
    return new


def _dump():
    path = os.environ.get('DQI_DQ3_W901_OUT')
    rec = dict(STATS, ceil=CEIL, patched=PATCHED)
    if path:
        try:
            json.dump(rec, open(path, 'w'))
        except Exception:
            pass
    sys.stderr.write('dqi_dq3_w901: %s\n' % json.dumps(dict(STATS, ceil=CEIL)))


def install(km, ns):
    if os.environ.get('DQI_DQ3_W901_OFF', '0') == '1':
        return dict(off=True)
    for obj in gc.get_objects():
        if not isinstance(obj, types.FunctionType):
            continue
        co = obj.__code__
        if _has(co):
            new = _rebuild(co)
            if new is not co:
                obj.__code__ = new
                STATS['functions'] += 1
        if _kit(co) and (not NAMES or co.co_name in NAMES):
            for name, cell in zip(co.co_freevars, obj.__closure__ or ()):
                try:
                    if name == 'ceiling' and type(cell.cell_contents) is int and cell.cell_contents in OLD:
                        cell.cell_contents = CEIL
                        STATS['cells'] += 1
                except ValueError:
                    pass
            if obj.__defaults__ and 'ceiling' in co.co_varnames[:co.co_argcount]:
                d = tuple(CEIL if type(c) is int and c in OLD else c for c in obj.__defaults__)
                if d != obj.__defaults__:
                    obj.__defaults__ = d
                    STATS['defaults'] += 1
    atexit.register(_dump)
    return dict(STATS, ceil=CEIL)
