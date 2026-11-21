#!/usr/bin/env python3
"""dqi_dq3_hflag.py -- DQ3 (Nov 18): ROW 'hflag' of the census: a zero flag F = [a == 0] erased by the kit's
measured form instead of a second AND ladder.

The line's converter erases some flags with exp/safe_scalar/converter.clear_flag = zero_toggle(a, F) + free: the
11-AND ladder again, on every shot.  The kit's own final_decoder/conversion.clear_flag erases the same flag by HMR
(outcome m, uniform), then ONLY when m = 1: X(a), a 10-AND tree, CCZ(tree, a[-1], m) (a classical control: a free
CZ_IF), the tree uncomputed by del_and (free), X(a): 10 Toffolis executed half the time = 5 expected [P].  Both
forms need F = [a == 0] on entry with a live, and both leave F erased and a unchanged; the measured form is
exact for that precondition (the HMR phase (-1)^(m [a == 0]) is repaid exactly), and it is narrower (10 nodes, not
11).  The ladder form is the one in the callers' scopes (clip_close, drop, merged_step, ...); this module gives every
function object whose code is safe_scalar's clear_flag the measured form's code (no globals are read by it).
Install: token dqi_dq3_hflag LAST in DQI_INT_EXT.  Switches: DQI_DQ3_HFLAG_OFF=1, DQI_DQ3_HFLAG_OUT=<json>.
"""
import atexit
import gc
import json
import os
import sys
import types

STATS = dict(functions=0, codes=0)

_SRC = '''
def clear_flag(b, a, F):
    m = b.alloc_bits(1); b.hmr(F, m[0]); b.free(F)
    b.push_condition(m[0]); b.x(a)
    nodes = b.alloc_qubits(len(a) - 2)
    b.init_and(a[0], a[1], nodes[0])
    for i in range(1, len(nodes)): b.init_and(nodes[i - 1], a[i + 1], nodes[i])
    b.ccz(nodes[-1], a[-1], m[0])
    for i in reversed(range(len(nodes))):
        if i: b.del_and(nodes[i - 1], a[i + 1], nodes[i])
        else: b.del_and(a[0], a[1], nodes[0])
    b.free(nodes); b.x(a); b.pop_condition(); b.free(m)
'''


def _target(fn):
    co = getattr(fn, '__code__', None)
    return (isinstance(fn, types.FunctionType) and co is not None and co.co_name == 'clear_flag'
            and co.co_filename.endswith('safe_scalar/converter.py') and co.co_argcount == 3
            and 'zero_toggle' in co.co_names)


def _dump():
    path = os.environ.get('DQI_DQ3_HFLAG_OUT')
    if path:
        try:
            json.dump(STATS, open(path, 'w'))
        except Exception:
            pass
    sys.stderr.write('dqi_dq3_hflag: %s\n' % json.dumps(STATS))


def install(km, ns):
    if os.environ.get('DQI_DQ3_HFLAG_OFF', '0') == '1':
        return dict(off=True)
    new = {}
    for obj in gc.get_objects():
        if _target(obj):
            co = obj.__code__
            if id(co) not in new:
                g = {'__builtins__': __builtins__, 'len': len, 'range': range, 'reversed': reversed}
                pad = '\n' * (co.co_firstlineno - 1)
                exec(compile(pad + _SRC, co.co_filename, 'exec'), g)
                new[id(co)] = (co, g['clear_flag'].__code__)
                STATS['codes'] += 1
            obj.__code__ = new[id(co)][1]
            STATS['functions'] += 1
    atexit.register(_dump)
    return dict(STATS)
