#!/usr/bin/env python3
"""dqi_int_livedump.py -- B5 (17th holder, Nov 16): diagnostics only (no gate change).

A DQI_INT_EXT module (install it FIRST) that records the allocation stack of every live qubit and, each time an
allocation brings the live count to >= B5_LD_AT (default 902) while one of the trigger functions is on the call
stack (B5_LD_TRIG, comma list of function names; default the synthesis tail's kernels), writes the live set to
B5_LD_OUT (append): every live wire labelled by the innermost local variable that holds it (frame.var[index]) on the
current call stack, and the unlabelled ones grouped by their allocation stack.  At most B5_LD_MAX dumps (default 30).
B5_LD_OUT unset: installs nothing.
"""
import os
import sys

OUT = os.environ.get('B5_LD_OUT')
AT = int(os.environ.get('B5_LD_AT', '902'))
MAXN = int(os.environ.get('B5_LD_MAX', '30'))
TRIG = set(os.environ.get('B5_LD_TRIG', 'tail_inplace,tail_finish,t4_forward,t4_reverse,reverse_masked_multiply,'
                          'reverse_deferred_multiply,masked_update').split(','))
LOG = dict(dumps=0, allocs=0)
LIVE = {}            # qubit id -> allocation stack (function:line, innermost first)
QT = []              # the qubit type (kickmix.q), learnt at install
SEQ = (list, tuple)  # containers walked for qubits; kickmix.array added at install


def _frames(skip):
    fr = sys._getframe(skip)
    while fr is not None:
        yield fr
        fr = fr.f_back


def _stack(skip=3, n=9):
    out = []
    for fr in _frames(skip):
        out.append('%s:%d' % (fr.f_code.co_name, fr.f_lineno))
        if len(out) >= n:
            break
    return ' < '.join(out)


def _ids(obj, depth=0):
    """(index path, id) of every qubit inside obj (a qubit or nested lists / tuples of them)."""
    if QT and isinstance(obj, QT[0]):
        yield '', obj.id
    elif isinstance(obj, SEQ) and depth < 3 and len(obj) <= 4096:
        for i, x in enumerate(obj):
            for p, q in _ids(x, depth + 1):
                yield '[%d]%s' % (i, p), q


def _ranges(ids):
    ids = sorted(ids)
    out, s, e = [], None, None
    for i in ids:
        if s is None:
            s = e = i
        elif i == e + 1:
            e = i
        else:
            out.append('%d' % s if s == e else '%d-%d' % (s, e))
            s = e = i
    if s is not None:
        out.append('%d' % s if s == e else '%d-%d' % (s, e))
    return ','.join(out)


def _dump(builder):
    label = {}
    names = []
    for fr in _frames(3):
        names.append(fr.f_code.co_name)
        for var, val in list(fr.f_locals.items()):
            if var in ('self', 'b', 'builder'):
                continue
            for path, q in _ids(val):
                if q in LIVE and q not in label:
                    label[q] = '%s.%s%s' % (fr.f_code.co_name, var, path)
        if len(names) > 60:
            break
    with open(OUT, 'a') as f:
        f.write('=== dump %d: live %d (tracked %d) ops %s\n    at %s\n' % (
            LOG['dumps'], builder.num_allocated_qubits, len(LIVE), getattr(builder, "num_ops", "?"), _stack(4, 14)))
        groups = {}
        for q, lab in label.items():
            base = lab.split('[')[0]
            groups.setdefault(base, []).append(q)
        for base in sorted(groups, key=lambda k: min(groups[k])):
            f.write('  L %-48s %4d  %s\n' % (base, len(groups[base]), _ranges(groups[base])))
        un = {}
        for q, st in LIVE.items():
            if q not in label:
                key = ' < '.join(st.split(' < ')[:6])
                un.setdefault(key, []).append(q)
        for key in sorted(un, key=lambda k: -len(un[k])):
            f.write('  U %4d  %s\n        ids %s\n' % (len(un[key]), key, _ranges(un[key])))


def install(km, ns):
    if not OUT:
        return dict(on=False)
    global SEQ
    kmod = sys.modules['kickmix']
    CB = kmod.CircuitBuilder
    if getattr(CB.alloc_qubits, '_b5_livedump', False):
        return dict(on=True, again=True)
    QT.append(kmod.q)
    SEQ = (list, tuple, kmod.array)
    orig_alloc, orig_free = CB.alloc_qubits, CB.free

    def alloc_qubits(self, n, *a, **kw):
        r = orig_alloc(self, n, *a, **kw)
        st = _stack(3, 9)
        for _, q in _ids(r):
            LIVE[q] = st
        LOG['allocs'] += 1
        if LOG['dumps'] < MAXN and self.num_allocated_qubits >= AT:
            if any(fr.f_code.co_name in TRIG for fr in _frames(1)):
                _dump(self)
                LOG['dumps'] += 1
        return r

    def free(self, x, *a, **kw):
        if QT:
            for _, q in _ids(x):
                LIVE.pop(q, None)
        return orig_free(self, x, *a, **kw)

    alloc_qubits._b5_livedump = True
    CB.alloc_qubits = alloc_qubits
    CB.free = free
    import atexit
    atexit.register(lambda: print('dqi_int_livedump', LOG, file=sys.stderr))
    return dict(on=True, at=AT)
