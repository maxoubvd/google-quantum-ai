"""dqi_int_zpk1.py -- B5 (16th holder, Nov 16): COPY of r46_dqi_conv_zpark.py with a call-clock TRACE and an ORACLE
eviction schedule. Without B5_ZP_TRACE / B5_ZP_ORACLE it emits exactly the original's circuit.
  B5_ZP_TRACE=<file.json>   pass 1: record, on a clock that counts only kit-driven calls (fresh / hold / consume /
                            clear wrappers and the cswaps of the role-swap frames), every 12-wire test (clock, key
                            class) and every key touched by a role swap; written at exit. No behaviour change.
  B5_ZP_ORACLE=<file.json>  pass 2: a pass-1 trace of the SAME recipe. A park is DOOMED when the trace shows its
                            key's next test only after a role swap touching it, or never: B5_ZP_DOOM=1 evicts a doomed
                            park at once (frees its wire; the eviction is the same HMR into the same mask, only
                            earlier); B5_ZP_BELADY=1 evicts, at the cap, doomed parks first, then the park whose next
                            test is furthest (zpark: oldest first). Every eviction is exact at any time while the flag
                            is parked (R46: an HMR whose phase the owner's consume repays); only the order changes.
  B5_ZP_DEFER=1 (with an oracle): a kit clear (zero_toggle again: 11 ANDs) whose key's next 12-wire test comes before
                            any role swap touching it and before the stage's end is replaced by the kit's own hold
                            (HMR into a new bit, 0 ANDs); the next test of that key repays it with the kit's consume
                            (CZ of its flag with the bit). Exact where the two tests return the same predicate (R46's
                            one-predicate-per-slot measurement, the same property zpark's hand-overs rest on).
                            Fallback (an oracle miss): before a swap / a free of the word / the stage's end, a real test,
                            the consume and a clear (22).
  B5_ZP_CPARK=1 (with an oracle): the same kit clears are PARKED instead (the flag kept, a bit stored 0): the key's
                            next test takes the flag (no test, no clear: -22); a cap eviction turns it into the
                            deferred clear above (HMR into the bit, repaid at the next test: -11); a role swap or the
                            stage's end clears it the kit's way first (+0 against today).
  B5_ZP_VICTIMS=1 (with an oracle that holds 'evict' records): a park whose creation clock the oracle run saw end in a
                            CAP eviction is evicted at its creation (it never reached a hand-over there; its wire
                            time goes to the rest and to P1's links).
  The trace also records stage exits (clock of the last tick), clears (clock, key class, fold / kit) and evictions
  (clock, why, the park's creation clock, key class).
"""
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
LOG = dict(b5_doom=0, b5_reorder=0, b5_k1_miss=0, b5_clock=0, b5_kit_clear=0, b5_defer=0, b5_defer_paid=0,
           b5_defer_fallback=0, b5_victim=0, b5_cpark=0, b5_cpark_handed=0, b5_cpark_deferred=0, b5_cpark_cleared=0, scopes=0, real=0, handed=0, parked=0, clear_parked=0, evict_alloc=0, evict_swap=0, evict_consume=0,
           evict_stage=0, evict_cond=0, evict_incond=0, over_cap_in_cond=0, most_parked=0)
S = {}
RAW = {}
TRACE = os.environ.get('B5_ZP_TRACE') or None
ORACLE = os.environ.get('B5_ZP_ORACLE') or None
DOOM = os.environ.get('B5_ZP_DOOM', '0') == '1'
BELADY = os.environ.get('B5_ZP_BELADY', '0') == '1'
DEFER = os.environ.get('B5_ZP_DEFER', '0') == '1'
CPARK = os.environ.get('B5_ZP_CPARK', '0') == '1'
VICTIMS = os.environ.get('B5_ZP_VICTIMS', '0') == '1'
CLK = [0]
KC = {}                      # key tuple -> key class (order of first test) in THIS run
TR = dict(fresh=[], swap=[], stage=[], clear=[], evict=[])
OR = None                    # oracle: fresh clock -> key class; key class -> sorted fresh clocks / swap clocks
if ORACLE:
    import bisect
    import json as _json
    _o = _json.load(open(ORACLE))
    OR = dict(k1={c: k for c, k in _o['fresh']}, fc={}, sc={}, stage=sorted(_o.get('stage', [])),
              vict=set(c0 for _c, why, c0, _k in _o.get('evict', []) if why in ('alloc', 'incond')))
    for c, k in _o['fresh']:
        OR['fc'].setdefault(k, []).append(c)
    for c, k in _o['swap']:
        OR['sc'].setdefault(k, []).append(c)
    for d in (OR['fc'], OR['sc']):
        for k in d:
            d[k].sort()


def tick():
    CLK[0] += 1
    LOG['b5_clock'] = CLK[0]
    return CLK[0]


def kclass(key):
    if key not in KC:
        KC[key] = len(KC)
    return KC[key]


def next_use(e, now):
    """(doomed, next test clock) of a park from the oracle; (False, None) when unknown."""
    k1 = e.get('k1')
    if OR is None or k1 is None:
        return False, None
    fc = OR['fc'].get(k1, [])
    i = bisect.bisect_right(fc, now)
    if i >= len(fc):
        return True, None
    nf = fc[i]
    sc = OR['sc'].get(k1, [])
    j = bisect.bisect_right(sc, now)
    if j < len(sc) and sc[j] < nf:
        return True, None
    xs = OR['stage']
    x = bisect.bisect_left(xs, now)
    if x < len(xs) and xs[x] < nf:
        return True, None
    return False, nf


def pay_defers(b, st, key, F):
    """repay the deferred clears of `key` with the flag F of its new test (the kit's consume)."""
    for m, _a in st.get('defer', {}).pop(key, []):
        b.cz(F, m[0])
        b.free(m)
        LOG['b5_defer_paid'] += 1


def fallback_defers(b, st, keys, why='swap'):
    """an oracle miss: repay now by a real test of the word, the consume and a clear (raw primitives). Only before a
    role swap (the word still holds its predicate); a free of the word or the stage's end with a deferred clear open
    is an error (the predicate may be gone): the build stops."""
    if why != 'swap' or st['cond']:
        raise RuntimeError('dqi_int_zpk1: deferred clear still open at a %s%s' % (why, ' (in a condition)' if st['cond'] else ''))
    for key in list(keys):
        for m, a in st.get('defer', {}).pop(key, []):
            F2 = RAW['fresh'](b, a)
            b.cz(F2, m[0])
            b.free(m)
            RAW['clear'](b, a, F2)
            LOG['b5_defer_fallback'] += 1


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
    if TRACE:
        TR['evict'].append((CLK[0], why, e.get('c0'), kclass(key)))
    if e['own'] == 'clr':
        if why in ('swap', 'stage', 'cond'):
            RAW['clear'](b, e['word'], e['flag'])          # the kit's clear while the word still holds it
            b.free(e['mask'])
            LOG['b5_cpark_cleared'] += 1
        else:
            b.hmr(e['flag'], e['mask'][0])
            b.free(e['flag'])
            st.setdefault('defer', {}).setdefault(key, []).append((e['mask'], e['word']))
            for w in key:
                st.setdefault('dwire', {})[w] = key
            LOG['b5_cpark_deferred'] += 1
    elif e['own']:
        b.hmr(e['flag'], e['mask'][0])
        b.free(e['flag'])
    else:
        t = b.alloc_bits(1)
        b.hmr(e['flag'], t[0])
        b.free(e['flag'])
        b.cx(t[0], e['mask'][0])
        b.free(t)
    LOG['evict_' + why] = LOG.get('evict_' + why, 0) + 1


def make(name, orig):
    if name == 'fresh_zero_flag':
        RAW.setdefault('fresh', orig)
        def fresh_zero_flag(b, a, *r, **kw):
            st = state(b)
            now = tick()
            key = key_of(a) if (active(st) and not st['cond']) else None
            k1 = None
            if key is not None and len(key) == 12:
                if TRACE:
                    TR['fresh'].append((now, kclass(key)))
                if OR is not None:
                    k1 = OR['k1'].get(now)
                    if k1 is None:
                        LOG['b5_k1_miss'] += 1
                e = st['parked'].pop(key, None)
                if e is not None:
                    if e['own'] == 'clr':
                        b.free(e['mask'])
                        LOG['b5_cpark_handed'] += 1
                    LOG['handed'] += 1
                    st['live'][e['flag'].id] = key
                    st.setdefault('k1', {})[e['flag'].id] = k1
                    pay_defers(b, st, key, e['flag'])
                    return e['flag']
            F = orig(b, a, *r, **kw)
            if key is not None and len(key) == 12:
                LOG['real'] += 1
                st['live'][F.id] = key
                st.setdefault('k1', {})[F.id] = k1
                for w in key:
                    st['wire'][w] = key
                pay_defers(b, st, key, F)
            return F
        w = fresh_zero_flag
    elif name == 'hold_flag':
        def hold_flag(b, F, *r, **kw):
            st = state(b)
            now = tick()
            key = st['live'].pop(getattr(F, 'id', None), None) if active(st) else None
            if key is None or key in st['parked'] or st['cond']:
                m = orig(b, F, *r, **kw)
                if key is not None:
                    st['masks'][m[0].id] = (key, st['epoch'].get(key, 0), m)
                return m
            m = b.alloc_bits(1)
            b.bit_store(m, False)
            st['tick'] += 1
            st['parked'][key] = dict(flag=F, mask=m, own=True, t=st['tick'], k1=st.get('k1', {}).get(F.id), c0=now)
            st['masks'][m[0].id] = (key, st['epoch'].get(key, 0), m)
            LOG['parked'] += 1
            LOG['most_parked'] = max(LOG['most_parked'], len(st['parked']))
            if DOOM and next_use(st['parked'][key], now)[0]:
                evict(b, st, key, 'doom')
            elif VICTIMS and OR is not None and now in OR['vict']:
                LOG['b5_victim'] += 1
                evict(b, st, key, 'victim')
            return m
        w = hold_flag
    elif name == 'consume_flag':
        def consume_flag(b, F, m, *r, **kw):
            st = state(b)
            tick()
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
        RAW.setdefault('clear', orig)

        def clear_flag(b, a, F, *r, **kw):
            st = state(b)
            now = tick()
            key = st['live'].pop(getattr(F, 'id', None), None) if active(st) else None
            if FOLD and key is not None and key == key_of(a) and key not in st['parked'] and not st['cond']:
                ep = st['epoch'].get(key, 0)
                for mid, (k2, e2, m) in reversed(list(st['masks'].items())):
                    if k2 == key and e2 == ep:
                        st['tick'] += 1
                        st['parked'][key] = dict(flag=F, mask=m, own=False, t=st['tick'], k1=st.get('k1', {}).get(F.id),
                                                 c0=now)
                        LOG['clear_parked'] += 1
                        if TRACE:
                            TR['clear'].append((now, kclass(key), 'fold'))
                        if DOOM and next_use(st['parked'][key], now)[0]:
                            evict(b, st, key, 'doom')
                        elif VICTIMS and OR is not None and now in OR['vict']:
                            LOG['b5_victim'] += 1
                            evict(b, st, key, 'victim')
                        return None
            if key is not None and len(key) == 12 and key == key_of(a) and not st['cond']:
                LOG['b5_kit_clear'] += 1
                if TRACE:
                    TR['clear'].append((now, kclass(key), 'kit'))
                if (DEFER or CPARK) and OR is not None and key not in st['parked']:
                    dm, nf = next_use(dict(k1=st.get('k1', {}).get(F.id)), now)
                    if not dm and nf is not None and CPARK:
                        m = b.alloc_bits(1)
                        b.bit_store(m, False)
                        st['tick'] += 1
                        st['parked'][key] = dict(flag=F, mask=m, own='clr', t=st['tick'], k1=st.get('k1', {}).get(F.id),
                                                 word=a)
                        LOG['b5_cpark'] += 1
                        LOG['most_parked'] = max(LOG['most_parked'], len(st['parked']))
                        return None
                    if not dm and nf is not None:
                        m = b.alloc_bits(1)
                        b.hmr(F, m[0])
                        b.free(F)
                        st.setdefault('defer', {}).setdefault(key, []).append((m, a))
                        for w in key:
                            st.setdefault('dwire', {})[w] = key
                        LOG['b5_defer'] += 1
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
                    order = sorted(st['parked'].items(), key=lambda kv: kv[1]['t'])
                    if BELADY and OR is not None:
                        now = CLK[0]

                        def prio(kv):
                            dm, nf = next_use(kv[1], now)
                            if dm:
                                return (0, 0, kv[1]['t'])
                            if nf is None:
                                return (1, 0, kv[1]['t'])          # unknown: zpark's own rank among the rest
                            return (1, -nf, kv[1]['t'])
                        o2 = sorted(st['parked'].items(), key=prio)
                        if [k for k, _ in o2[:need]] != [k for k, _ in order[:need]]:
                            LOG['b5_reorder'] += 1
                        order = o2
                    for key, e in order:
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
        if st is not None and st.get('defer'):
            hit = set()
            for q in (x if isinstance(x, (list, tuple)) else [x]):
                try:
                    qs = list(q) if type(q).__name__ not in ('q', 'b') else [q]
                except TypeError:
                    qs = [q]
                for qq in qs:
                    if type(qq).__name__ != 'q':
                        continue
                    k = st.get('dwire', {}).get(getattr(qq, 'id', None))
                    if k is not None and k in st['defer']:
                        hit.add(k)
            if hit:
                fallback_defers(self, st, hit, 'free')
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
        now = CLK[0]
        if st is not None and sys._getframe(1).f_code.co_name in SWAPS:
            now = tick()
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
            if st.get('defer'):
                fallback_defers(self, st, [k for k in keys if k in st['defer']])
            for k in keys:
                if TRACE:
                    TR['swap'].append((now, kclass(k)))
                if k in st['parked']:
                    evict(self, st, k, 'swap')
                st['epoch'][k] = st['epoch'].get(k, 0) + 1
            if st.get('defer') and any(k in st['defer'] for k in keys):
                fallback_defers(self, st, [k for k in keys if k in st['defer']])
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
                if TRACE:
                    TR['stage'].append(CLK[0])
                if st.get('defer'):
                    fallback_defers(self.b, st, list(st['defer']), 'stage end')
                for key in list(st['parked']):
                    if key in st['parked']:
                        evict(self.b, st, key, 'stage')
                if st.get('defer'):
                    fallback_defers(self.b, st, list(st['defer']))
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
    if TRACE:
        def _dump():
            import json as _j
            with open(TRACE, 'w') as fh:
                _j.dump(dict(fresh=TR['fresh'], swap=TR['swap'], stage=TR['stage'], clear=TR['clear'],
                             evict=TR['evict'], keys=len(KC), clock=CLK[0]), fh)
        atexit.register(_dump)
    return dict(LOG)
