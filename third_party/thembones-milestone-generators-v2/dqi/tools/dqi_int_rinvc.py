#!/usr/bin/env python3
"""dqi_int_rinvc.py -- B5 (22nd holder, Nov 18): the reverse 24-bit in-place inverse of the line (B7's
inverse24_reverse_i27, tools/dqi_var_rinv27.py, read only) with its INNER 12-bit inverse t = 1/N taken by the kit's
flag-free cached 12-bit inverse (decode_flagfree_cached.CACHED12[0]: the norm's products kept as raw wires, the init18
leaf, a 25-product final stage: 58) in place of Norm 15 + init19 19 + two 6-bit products 30 = 64, wherever the inner
inverse starts with live + EXTRA <= CEIL: -6 a call.  Both forms are zero-correct (1/0 = 0) and leave N in place and
t = 1/N in the same basis; nothing else in the kernel changes.
How: the entry installs the DQI_INT_EXT modules (token 'ext') BEFORE B7's var_rinv27, so this module imports
dqi_var_rinv27 (the entry's later import by the same name returns this object) and rewrites its source string SRC
before B7's make() compiles it: the five inner lines become one call into _inner() below, which keeps those lines as
its fallback path.
Install: token dqi_int_rinvc in DQI_INT_EXT (B5_FILES_ADD=dqi_int_rinvc.py).  Switches: DQI_INT_RC_CEIL (901),
DQI_INT_RC_EXTRA (22), DQI_INT_RC_NZ=1 (the inner leaf at 17 by tools/dqi_int_leafnz.py: these inputs are INV17
list A's junction norms, proved nonzero by R46-DQI-EXT; the file-level INV17 map no longer finds them), DQI_INT_RC_OFF=1, DQI_INT_RC_OUT=<json> (per call: live at the inner entry, form).
"""
import atexit
import importlib
import json
import os
import sys

CEIL = int(os.environ.get('DQI_INT_RC_CEIL', '901'))
EXTRA = int(os.environ.get('DQI_INT_RC_EXTRA', '22'))
NZ = os.environ.get('DQI_INT_RC_NZ', '0') == '1'   # the 14 junction sites of INV17's list A: N = Norm(x), x nonzero (R46-DQI-EXT)
STATS = dict(cached=0, orig=0, rewritten=0)
CALLS = []
_KIT = {}
sys.modules.setdefault('dqi_int_rinvc', sys.modules[__name__])

OLD = ('''    to12, frm12, K12 = tower(12); linear_inplace(b, N, to12)
    D = b.alloc_qubits(6); norm_xor_noscratch(b, N, D, K12); u = inverse6_init19(b, D); ar._erase_norm(b, N, D, K12); b.free(D)
    t = b.alloc_qubits(12); b.cx(N[6:], N[:6]); mul_xor_noscratch(b, N[:6], u, t[:6]); mul_xor_noscratch(b, N[6:], u, t[6:]); b.cx(N[6:], N[:6])
    ar._erase_norm(b, t, u, K12); b.free(u)
    linear_inplace(b, N, frm12); linear_inplace(b, t, frm12)
''')
NEW = ('''    t = __import__('sys').modules['dqi_int_rinvc']._inner(b, N, tower, linear_inplace, norm_xor_noscratch, inverse6_init19, ar, mul_xor_noscratch)
''')


def _inner(b, N, tower, linear_inplace, norm_xor_noscratch, inverse6_init19, ar, mul_xor_noscratch):
    live = b.num_allocated_qubits
    if live + EXTRA <= CEIL:
        if not _KIT:
            pj = sys.modules.get('exp.jet24.pro_jet') or importlib.import_module('exp.jet24.pro_jet')
            if NZ:      # the leaf with INV17's redundant AND written affine at generation time (input proved nonzero)
                _KIT['cached12'] = importlib.import_module('dqi_int_leafnz').make(pj)[1]
            else:
                ff = sys.modules.get('decode_flagfree_cached') or importlib.import_module('decode_flagfree_cached')
                _KIT['cached12'] = ff.CACHED12[0]
            _KIT['plan12'] = pj.PLAN12
        STATS['cached'] += 1
        CALLS.append((live, 1))
        t, fl = _KIT['cached12'](b, N, _KIT['plan12'], False)
        assert fl is None
        return t
    STATS['orig'] += 1
    CALLS.append((live, 0))
    to12, frm12, K12 = tower(12); linear_inplace(b, N, to12)
    D = b.alloc_qubits(6); norm_xor_noscratch(b, N, D, K12); u = inverse6_init19(b, D); ar._erase_norm(b, N, D, K12); b.free(D)
    t = b.alloc_qubits(12); b.cx(N[6:], N[:6]); mul_xor_noscratch(b, N[:6], u, t[:6]); mul_xor_noscratch(b, N[6:], u, t[6:]); b.cx(N[6:], N[:6])
    ar._erase_norm(b, t, u, K12); b.free(u)
    linear_inplace(b, N, frm12); linear_inplace(b, t, frm12)
    return t


def _dump():
    path = os.environ.get('DQI_INT_RC_OUT')
    if path:
        try:
            json.dump(dict(STATS, ceil=CEIL, extra=EXTRA, calls=CALLS), open(path, 'w'))
        except Exception:
            pass
    sys.stderr.write('dqi_int_rinvc: %s\n' % json.dumps(dict(STATS, ceil=CEIL, extra=EXTRA)))


def install(km, ns):
    if os.environ.get('DQI_INT_RC_OFF', '0') == '1':
        return dict(STATS, off=True)
    b7 = importlib.import_module('dqi_var_rinv27')
    src = b7.SRC
    assert src.count(OLD) == 1, 'dqi_int_rinvc: the inner inverse lines of inverse24_reverse_i27 not found once'
    b7.SRC = src.replace(OLD, NEW)
    STATS['rewritten'] = 1
    atexit.register(_dump)
    return dict(rewritten=1, ceil=CEIL, extra=EXTRA, nz=NZ)
