"""r46_dqi_conv_zpark.py -- R46-DQI-CONV: keep a gap flag instead of testing its word again (convert / unconvert).

Measured on best_q23j5's recipe (research/panel/R46-DQI-CONV.md): the converter makes 131 twelve-bit zero tests per
decode on 16 word slots; all tests of one slot give the same answer (the level's gap flag), and between two lives of a
flag a wire is free on 59 of the 94 gaps. This EXT module (append LAST to DQI_INT_EXT) bridges those gaps with no
knowledge of the call sites: it wraps, in every scope and module that has them, the kit's four flag primitives.
  hold_flag(b, F)          the kit X-measures F and returns the mask. Here: F stays allocated ('parked') and the
                           caller gets a mask bit stored 0; if the park is evicted later, the HMR outcome is written
                           into that very bit, so the caller's own consume_flag repays it as before.
  fresh_zero_flag(b, a)    if a flag of the same 12 wires is parked: hand it over (no zero test). Else the kit's test.
  clear_flag(b, a, F)      (R46_ZP_FOLD=1) if a mask of the same wires and the same epoch is still to be consumed:
                           park F; an eviction then XORs its HMR outcome into that mask. Else the kit's clear.
  consume_flag(b, F, m)    a park that still writes into m is evicted first; then the kit's CZ.
Evictions (always exact: an HMR whose phase is repaid by the owner's consume): (1) an allocation that would pass
R46_ZP_CAP (default 902) qubits, oldest park first; (2) a role swap (cswap emitted by a frame named aligned,
canonical_permute or rotate) that touches the slot's wires: the slot's word changes there, the epoch of the slot is
advanced; (3) leaving the stage (parks exist only inside convert and unconvert); (4) entering a classically
conditioned block with fewer than R46_ZP_MARGIN (40) free wires (no eviction can be emitted inside a condition).
The hand-over is exact where the kit's own second test would have returned the parked value: measured on 262,144
samples for every test of every slot (tools/r46_dqi_conv_zsame.py: 0 disagreements).
"""
import atexit
import os
import sys

CAP = int(os.environ.get('R46_ZP_CAP', '902'))
MARGIN = int(os.environ.get('R46_ZP_MARGIN', '40'))
# R46_ZP_INCOND=1 (v2): no eviction at a block's entry; when an allocation inside a conditioned block would pass the
# cap, the open conditions are popped, the park is evicted unconditionally and the same conditions are pushed again
# (measured on the line: 21,131 outermost conditioned blocks in convert + unconvert, 27 of them allocate, at most 11
# wires, tools/r46_dqi_conv_condneed.py). Exact as long as a block does not rewrite its own condition bit.
INCOND = os.environ.get('R46_ZP_INCOND', '0') == '1'
FOLD = os.environ.get('R46_ZP_FOLD', '1') == '1'
STAGES = ('convert', 'unconvert')
SWAPS = ('aligned', 'canonical_permute', 'rotate')
NAMES = ('fresh_zero_flag', 'hold_flag', 'consume_flag', 'clear_flag')
LOG = dict(scopes=0, real=0, handed=0, parked=0, clear_parked=0, evict_alloc=0, evict_swap=0, evict_consume=0,
           evict_stage=0, evict_cond=0, evict_incond=0, over_cap_in_cond=0, most_parked=0)
S = {}
RAW = {}


def new_state():
    return dict(parked={}, live={}, masks={}, epoch={}, wire={}, cond=0, stage=[], tick=0, condargs=[])


def state(b):
    st = S.get(id(b))
    if st is None:
        st = S[id(b)] = new_state()
    return st


def active(st):
    return bool(st['stage']) and st['stage'][0] in STAGES


def key_of(a):
    try:
        return tuple(sorted(q.id for q in a))
    except Exception:
        return None


def evict(b, st, key, why):
    if st['cond'] and INCOND and RAW:
        for _ in st['condargs']:
            RAW['pop'](b)
        st['cond'], held = 0, list(st['condargs'])
        try:
            evict(b, st, key, 'incond' if why == 'alloc' else why)
        finally:
            for a, kw in held:
                RAW['push'](b, *a, **kw)
            st['cond'] = len(held)
        return
    e = st['parked'].pop(key)
    if e['own']:
        b.hmr(e['flag'], e['mask'][0])
        b.free(e['flag'])
    else:
        t = b.alloc_bits(1)
        b.hmr(e['flag'], t[0])
        b.free(e['flag'])
        b.cx(t[0], e['mask'][0])
        b.free(t)
    LOG['evict_' + why] += 1


def make(name, orig):
    if name == 'fresh_zero_flag':
        def fresh_zero_flag(b, a, *r, **kw):
            st = state(b)
            key = key_of(a) if (active(st) and not st['cond']) else None
            if key is not None and len(key) == 12:
                e = st['parked'].pop(key, None)
                if e is not None:
                    LOG['handed'] += 1
                    st['live'][e['flag'].id] = key
                    return e['flag']
            F = orig(b, a, *r, **kw)
            if key is not None and len(key) == 12:
                LOG['real'] += 1
                st['live'][F.id] = key
                for w in key:
                    st['wire'][w] = key
            return F
        w = fresh_zero_flag
    elif name == 'hold_flag':
        def hold_flag(b, F, *r, **kw):
            st = state(b)
            key = st['live'].pop(getattr(F, 'id', None), None) if active(st) else None
            if key is None or key in st['parked'] or st['cond']:
                m = orig(b, F, *r, **kw)
                if key is not None:
                    st['masks'][m[0].id] = (key, st['epoch'].get(key, 0), m)
                return m
            m = b.alloc_bits(1)
            b.bit_store(m, False)
            st['tick'] += 1
            st['parked'][key] = dict(flag=F, mask=m, own=True, t=st['tick'])
            st['masks'][m[0].id] = (key, st['epoch'].get(key, 0), m)
            LOG['parked'] += 1
            LOG['most_parked'] = max(LOG['most_parked'], len(st['parked']))
            return m
        w = hold_flag
    elif name == 'consume_flag':
        def consume_flag(b, F, m, *r, **kw):
            st = state(b)
            try:
                mid = m[0].id
            except Exception:
                mid = None
            if mid is not None:
                for key, e in list(st['parked'].items()):
                    if e['mask'][0].id == mid:
                        evict(b, st, key, 'consume')
                st['masks'].pop(mid, None)
            return orig(b, F, m, *r, **kw)
        w = consume_flag
    else:
        def clear_flag(b, a, F, *r, **kw):
            st = state(b)
            key = st['live'].pop(getattr(F, 'id', None), None) if active(st) else None
            if FOLD and key is not None and key == key_of(a) and key not in st['parked'] and not st['cond']:
                ep = st['epoch'].get(key, 0)
                for mid, (k2, e2, m) in reversed(list(st['masks'].items())):
                    if k2 == key and e2 == ep:
                        st['tick'] += 1
                        st['parked'][key] = dict(flag=F, mask=m, own=False, t=st['tick'])
                        LOG['clear_parked'] += 1
                        return None
            return orig(b, a, F, *r, **kw)
        w = clear_flag
    w._r46zp = True
    w.__name__ = name
    return w


def install(km, ns):
    km = sys.modules['kickmix']
    B = km.CircuitBuilder
    walk = sys.modules['walk_bound']
    dicts = {}
    for fn in walk.functions(ns.decode):
        dicts[id(fn.__globals__)] = fn.__globals__
    for mod in list(sys.modules.values()):
        d = getattr(mod, '__dict__', None)
        if isinstance(d, dict) and any(n in d for n in NAMES) and mod.__name__ != __name__:
            dicts[id(d)] = d
    for d in dicts.values():
        hit = False
        for n in NAMES:
            f = d.get(n)
            if callable(f) and not getattr(f, '_r46zp', False):
                d[n] = make(n, f)
                hit = True
        LOG['scopes'] += 1 if hit else 0
    if getattr(B, '_r46zp', False):
        return dict(LOG)
    B._r46zp = True
    o_init, o_alloc, o_free, o_cswap = B.__init__, B.alloc_qubits, B.free, B.cswap
    o_push, o_pop, o_mark = B.push_condition, B.pop_condition, B.mark
    RAW.update(push=o_push, pop=o_pop)

    def init(self, *a, **kw):
        o_init(self, *a, **kw)
        S[id(self)] = new_state()

    def alloc_qubits(self, n, *a, **kw):
        st = S.get(id(self))
        if st is not None and st['parked']:
            need = self.num_allocated_qubits + n - CAP
            if need > 0:
                if st['cond'] and not INCOND:
                    LOG['over_cap_in_cond'] += 1
                else:
                    for key, e in sorted(st['parked'].items(), key=lambda kv: kv[1]['t']):
                        if need <= 0:
                            break
                        evict(self, st, key, 'alloc')
                        need -= 1
        return o_alloc(self, n, *a, **kw)

    def gone(st, q):
        tn = type(q).__name__
        if tn == 'q':
            st['live'].pop(q.id, None)
        elif tn == 'b' and q.id in st['masks']:
            # a mask freed without consume_flag: no park may still write into it
            for e in st['parked'].values():
                if e['mask'][0].id == q.id:
                    raise RuntimeError('r46_zpark: mask freed outside consume_flag while a park depends on it')
            del st['masks'][q.id]

    def free(self, x, *a, **kw):
        st = S.get(id(self))
        if st is not None and (st['live'] or st['masks']):
            for q in (x if isinstance(x, (list, tuple)) else [x]):
                if type(q).__name__ in ('q', 'b'):
                    gone(st, q)
                else:
                    try:
                        for qq in q:
                            gone(st, qq)
                    except TypeError:
                        pass
        return o_free(self, x, *a, **kw)

    def cswap(self, c, x, y, *a, **kw):
        st = S.get(id(self))
        if st is not None and st['wire'] and sys._getframe(1).f_code.co_name in SWAPS:
            keys = set()
            for w in (x, y):
                try:
                    for q in w:
                        k = st['wire'].get(q.id)
                        if k is not None:
                            keys.add(k)
                except TypeError:
                    k = st['wire'].get(getattr(w, 'id', None))
                    if k is not None:
                        keys.add(k)
            for k in keys:
                if k in st['parked']:
                    evict(self, st, k, 'swap')
                st['epoch'][k] = st['epoch'].get(k, 0) + 1
        return o_cswap(self, c, x, y, *a, **kw)

    def push_condition(self, *a, **kw):
        st = S.get(id(self))
        if st is not None:
            if not INCOND and st['cond'] == 0 and st['parked'] and self.num_allocated_qubits + MARGIN > CAP:
                for key in list(st['parked']):
                    evict(self, st, key, 'cond')
            st['cond'] += 1
            st['condargs'].append((a, kw))
        return o_push(self, *a, **kw)

    def pop_condition(self, *a, **kw):
        st = S.get(id(self))
        if st is not None and st['cond']:
            st['cond'] -= 1
            st['condargs'].pop()
        return o_pop(self, *a, **kw)

    class Marked:
        def __init__(self, inner, tag, b):
            self.inner, self.tag, self.b = inner, tag, b

        def __enter__(self):
            state(self.b)['stage'].append(self.tag)
            return self.inner.__enter__()

        def __exit__(self, *a):
            st = state(self.b)
            if len(st['stage']) == 1 and st['stage'][0] in STAGES and a[0] is None:
                for key in list(st['parked']):
                    evict(self.b, st, key, 'stage')
                st['live'].clear()
            st['stage'].pop()
            return self.inner.__exit__(*a)

    B.__init__ = init
    B.alloc_qubits = alloc_qubits
    B.free = free
    B.cswap = cswap
    B.push_condition = push_condition
    B.pop_condition = pop_condition
    B.mark = lambda self, tag, *a, **kw: Marked(o_mark(self, tag, *a, **kw), tag, self)
    atexit.register(lambda: print('r46_zpark', dict(LOG), file=sys.stderr, flush=True))
    return dict(LOG)
