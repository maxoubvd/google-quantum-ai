"""Round 7 diagnostic patch: executed-Toffoli profile by call stack (workstream dqi_build).

Every Toffoli-like builder call (ccx, ccz, cswap, init_and with quantum operands) is keyed by the Python call
stack that emits it and by the stack of `b.mark` tags. Whenever the key changes, one DEBUG_PRINT marker is put
into the circuit (no gate, no qubit; the official verifier ignores it) and the key is appended to the segment
list. `tools/dqi_build_r7_prof_run.py` then reads the simulator's executed Toffoli counter at every marker, so
each executed Toffoli is charged to the stack that emitted it (static counts are recorded here as well).
Writes DQI_R7_PROF_OUT (default ./segments.json): {"stacks": [[frame, ...], ...], "marks": [[tag, ...], ...],
"segments": [[stack id, marks id, static Toffolis], ...]} in circuit order.
"""
import atexit
import json
import os
import sys

OUT = os.environ.get('DQI_R7_PROF_OUT', 'segments.json')
PAID = ('ccx', 'ccz', 'cswap', 'init_and')
SKIP = {'<lambda>', 'paid', '<module>', '<listcomp>', '<genexpr>', 'counted'}


def install(km, module):
    km = sys.modules['kickmix']
    B = km.CircuitBuilder
    o_dbg, o_mark = B.debug_print, B.mark
    probe = km.q(0)
    stacks, stack_ids = [], {}
    marks, mark_ids = [], {}
    segments = []
    tags = []
    state = dict(key=None, row=None, code_cache={})

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

    def frame_name(code):
        s = state['code_cache'].get(code)
        if s is None:
            s = code.co_name + '@' + os.path.basename(code.co_filename)[:-3]
            state['code_cache'][code] = s
        return s

    def wrap(name):
        original = getattr(B, name)

        def paid(self, *a, **kw):
            n = width(a, kw)
            if n:
                names = []
                f = sys._getframe(1)
                while f is not None:
                    co = f.f_code
                    if co.co_name not in SKIP:
                        names.append(frame_name(co))
                    f = f.f_back
                key = (tuple(names), tuple(tags))
                if key != state['key']:
                    si = stack_ids.get(key[0])
                    if si is None:
                        si = stack_ids[key[0]] = len(stacks)
                        stacks.append(list(reversed(key[0])))
                    mi = mark_ids.get(key[1])
                    if mi is None:
                        mi = mark_ids[key[1]] = len(marks)
                        marks.append(list(key[1]))
                    o_dbg(self, probe)
                    state['key'] = key
                    state['row'] = [si, mi, 0]
                    segments.append(state['row'])
                state['row'][2] += n
            return original(self, *a, **kw)
        setattr(B, name, paid)

    for name in PAID:
        wrap(name)

    class Marked:
        def __init__(self, inner, tag):
            self.inner, self.tag = inner, tag

        def __enter__(self):
            tags.append(self.tag)
            return self.inner.__enter__()

        def __exit__(self, *a):
            tags.pop()
            return self.inner.__exit__(*a)

    B.mark = lambda self, tag, *a, **kw: Marked(o_mark(self, tag, *a, **kw), tag)

    def report():
        with open(OUT, 'w') as f:
            json.dump(dict(stacks=stacks, marks=marks, segments=segments), f)
        print('dqb_r7prof', len(segments), 'segments,', len(stacks), 'stacks,', sum(s[2] for s in segments),
              'static Toffolis ->', OUT, file=sys.stderr)

    atexit.register(report)
    return dict(out=OUT)
