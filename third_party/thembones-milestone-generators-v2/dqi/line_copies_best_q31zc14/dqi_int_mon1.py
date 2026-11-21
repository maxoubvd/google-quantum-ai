#!/usr/bin/env python3
"""dqi_int_mon1.py -- B5 (25th holder, Nov 19; ledger O-7 / rank 7): the converter's base level (m = 1) without its
`first` context (the 'monic' order of R62-DQI-CONVN idea 6 at level 1, priced by DQ5 on the line's code), with the
gap routing repaired by an exact F-controlled scaling of the parent's tops that costs 120 for both tops at +12 wires
(DQ5 had 2 x 96 at +12, or 2 x 66 at +24, which does not fit: research/ws_dqi_dq5.md L89, L119-124).

WHY THE VALUES ARE THE SAME [P; level 1 has single-word payloads v = [v0], so v[-1] is v[0]]:
 today (q7 base_with_fused_zero, forward): first: v0 /= a, q /= a; q_post: cswap(F, q, p), cswap(F, v0, t),
   v0 /= (q-slot), t ^= F c v0, cswap(F, q, p); last: v0 /= a, q *= a.
   F = 0: v0 -> v0 / (q a), q unchanged.  F = 1: t -> v0 / a + c t_old / p, v0-slot -> t_old / (p a).
 here: the line's own q_post with its masked products held back (captured); then ONE context of a: (1) every parent
   top t <- t s, s = F ? a^-1 : 1 (with F = 1 the child's v0 reaches t undivided); (2) the captured masked products,
   replayed unchanged; (3) v0 /= a.  Same values for F = 0 and F = 1. q is never scaled, so q_post's token is on the
   inverse of q itself, and the reverse pays it with ITS inverse of the same q: the reverse is the mirror:
 reverse: one context of a with the forward's token: v0 *= a, the masked products, t <- t s^-1 (s^-1 = F ? a : 1);
   then the line's reverse q_post with its masked products held back (its F is the one used above); then a08's tail
   WITHOUT its three scale_before_inverse (prefix, finish, consume of the carry, erase_normpair with the level-2
   hand-off of F AND [a = 0]); for a q7 reverse base, q7's first(undo) without its scalings.
THE CONTROLLED SCALING `cscale` (forward x = a^-1, y = a; reverse x = a, y = a^-1); exact on every input wherever
 x y = 1, which holds inside a context of a (the kit's a = 0 trick sets a = ai = 1 there): per top t: u = F t
 (12 ANDs), t ^= u (x + 1) (one K-product, 42), HMR u, whose phase is (-1)^{m . F (t y)}; then ONE shared w = F y
 (12 ANDs) repays every top's phase as the bilinear phase (-1)^{m . (t w)} (ar.phase_product: 0 Toffolis) and w's own
 HMR phase F (m2 . y) by conditioned CZs (0 Toffolis). Two tops: 2 x (12 + 42) + 12 = 120, at most +12 wires (the
 width of one kit scaling).
COUNT [P]: forward - first's context (58 cached / 67 kit) - its 3 scalings (126) - last's q scaling (42) + 120 =
 about -106; reverse - the reverse last's q scaling (42) - the tail's 3 scale_before_inverse (126) + 120 = -48;
 about -154 a decode.
GUARDS: every condition is checked before a gate is emitted (level 1, two single-word payloads, a parent, 2 in
 KEEP_CHILD, F handed in given_gap[2], a known reverse base); a refused forward refuses its reverse (the shared phase
 token carries '_mon1'); a reverse that does not match its forward raises (no silent mixed circuit).
Env: DQI_INT_MON1_OFF=1 (nothing installed), DQI_INT_MON1_OUT=<json> (LOG after every call).
Usage: python -I dqi_int_mon1.py test <kit copy root>   (cscale on every a != 0 x F in {0, 1} by kickmix, both forms)
"""
import json
import os
import sys
import types

NAME = 'dqi_int_mon1'
OFF = os.environ.get('DQI_INT_MON1_OFF', '0') == '1'
OUT = os.environ.get('DQI_INT_MON1_OUT')
NOTAIL = os.environ.get('DQI_INT_MON1_NOTAIL', '0') == '1'   # COUNT-ONLY diagnostic (harness): drops the reverse tail
#                                                               (not exact; never set in a build)
# MERGE (v3, fix round 2, default off): the reverse's ONE context of a also makes the level-2 hand-off F AND [a = 0]
# from its own zero flag and pays its own phase in place (the context function's erase with defer=False), so a08's
# tail (a second inverse of a for the carry and the hand-off: about 69 executed [M]) is not made: about -32 [E].
# The context is the line's own inverse_apply_safe step by step (its inverse_zero / consume / erase, read from its
# globals); the hand-off bit is live through the reverse q_post (+1 wire there).
MERGE = os.environ.get('DQI_INT_MON1_MERGE', '0') == '1'


_REF_IAS = ("def inverse_apply_safe(b,a,action,incoming=None,defer=False):\n"
            "    ai,F=inverse_zero(b,a);consume(b,ai,F,incoming)\n"
            "    b.cx(F,a[0]);b.cx(F,ai[0]);action(ai);b.cx(F,ai[0]);b.cx(F,a[0])\n"
            "    return erase(b,a,ai,F,defer)\n")


_REF_IAS_CLOSURE = ("def _outer():\n"
                    "    inverse_zero = scope = None\n"
                    "    def inverse_apply_safe(b, a, action, incoming=None, defer=False):\n"
                    "        ai, F = inverse_zero(b, a)\n"
                    "        scope['consume'](b, ai, F, incoming)\n"
                    "        b.cx(F, a[0])\n"
                    "        b.cx(F, ai[0])\n"
                    "        action(ai)\n"
                    "        b.cx(F, ai[0])\n"
                    "        b.cx(F, a[0])\n"
                    "        return scope['erase'](b, a, ai, F, defer)\n"
                    "    return inverse_apply_safe\n")


def _unwrap_ia(ia):
    """Follow the line's wrappers of inverse_apply down to an inverse_apply_safe: Z1's pdl2 wrapper (a closure over
    `orig` that runs late paired work queued in its module's LATE after the context returns: none is queued at level 1,
    which the merged path asserts) and the kit's one-line inverse_apply (return inverse_apply_safe(...)).
    Returns (the inverse_apply_safe or None, the chain followed)."""
    chain, f = [], ia
    for _ in range(6):
        if not isinstance(f, types.FunctionType):
            return None, chain
        chain.append(f)
        co = f.__code__
        if co.co_name == 'inverse_apply_safe':
            return f, chain
        cells = dict(zip(co.co_freevars, f.__closure__ or ()))
        if getattr(f, '_zn_pdl2', False) and 'orig' in cells:
            f = cells['orig'].cell_contents
            continue
        if co.co_name == 'inverse_apply' and co.co_names == ('inverse_apply_safe',) and not f.__closure__:
            f = f.__globals__.get('inverse_apply_safe')
            continue
        return None, chain
    return None, chain


def _late_guards(d):
    """The LATE queues of the pdl2 wrappers on the way to the context (their lengths must not change in it)."""
    _, chain = _unwrap_ia(d['inverse_apply'])
    return [f.__globals__['LATE'] for f in chain if getattr(f, '_zn_pdl2', False) and 'LATE' in f.__globals__]


def _code_of(src, *path):
    co = compile(src, '<ref>', 'exec')
    for name in path:
        co = [c for c in co.co_consts if isinstance(c, types.CodeType) and c.co_name == name][0]
    return co


def _context_fns(d):
    """(inverse_zero, consume, erase) of the line's context function behind d['inverse_apply'] when that function is
    the kit's inverse_apply_safe step for step, else None (a note in LOG['merge_refused']). Two forms are accepted, each
    by bytecode equality with a reference body: the kit's module function (inverse_zero, consume, erase as globals) and
    the cache owners' closures (a12 forward_norm_cache, a14 reverse_caches: inverse_zero a free variable, consume and
    erase read from their scope cell); local variable names may differ (a14 writes Z for F)."""
    ias, chain = _unwrap_ia(d['inverse_apply'])
    if ias is None:
        LOG['merge_refused'] = 'no inverse_apply_safe behind inverse_apply: ' + ' > '.join(
            '%s:%s:%d' % (os.path.basename(f.__code__.co_filename), f.__code__.co_name, f.__code__.co_firstlineno)
            for f in chain)
        return None
    co = ias.__code__
    if not ias.__closure__:
        ref = _code_of(_REF_IAS, 'inverse_apply_safe')
        if co.co_code == ref.co_code and co.co_names == ref.co_names:
            g2 = ias.__globals__
            return g2['inverse_zero'], g2['consume'], g2['erase']
    else:
        ref = _code_of(_REF_IAS_CLOSURE, '_outer', 'inverse_apply_safe')
        if co.co_code == ref.co_code and co.co_names == ref.co_names and co.co_freevars == ref.co_freevars:
            cells = dict(zip(co.co_freevars, ias.__closure__))
            sc = cells['scope'].cell_contents
            return cells['inverse_zero'].cell_contents, sc['consume'], sc['erase']
    LOG['merge_refused'] = '%s:%d' % (os.path.basename(co.co_filename), co.co_firstlineno)
    return None
LOG = dict(installed=False, off=OFF, scopes=0, mp_dicts=0, originals=[], calls={'fwd': 0, 'rev': 0}, refused={},
           captured={'fwd': 0, 'rev': 0}, tails={})
_MP_DICTS = []


def _dump():
    if OUT:
        try:
            with open(OUT, 'w') as fh:
                json.dump(LOG, fh)
        except OSError:
            pass


def _refuse(why):
    LOG['refused'][why] = LOG['refused'].get(why, 0) + 1


def cscale(b, ar, F, tops, x, y):
    """t <- t (1 + F (x + 1)) for every 12-wire word t in tops; exact wherever x y = 1. 2 x (12 + 42) + 12 for two."""
    masks = []
    b.x(x[0])
    for t in tops:
        u = b.alloc_qubits(12)
        for j in range(12):
            b.init_and(F, t[j], u[j])
        ar.mul_xor0(b, u, x, t)
        m = b.alloc_bits(12)
        b.hmr(u, m)
        b.free(u)
        masks.append(m)
    b.x(x[0])
    w = b.alloc_qubits(12)
    for j in range(12):
        b.init_and(F, y[j], w[j])
    for t, m in zip(tops, masks):
        ar.phase_product(b, t, w, m)
        b.free(m)
    m2 = b.alloc_bits(12)
    b.hmr(w, m2)
    b.free(w)
    for j in range(12):
        b.ccz(F, y[j], m2[j])
    b.free(m2)


def _kind(orig):
    co = getattr(orig, '__code__', None)
    if co is None:
        return None
    fn = co.co_filename.replace('\\', '/')
    if co.co_name == 'base_case' and fn.endswith('a08/scalar_compact_v3.py'):
        return 'a08'
    if co.co_name == 'base_with_fused_zero':
        return 'q7'
    return None


def _cells(fn):
    out = []
    for c in (fn.__closure__ or ()):
        try:
            out.append(c.cell_contents)
        except ValueError:
            pass
    return out


def _call_q_post(d, b, q, payloads, H, parent, reverse, phase):
    """The line's q_post at m = 1 with every masked_products / masked_product call held back; returns the calls."""
    calls, saved = [], []

    def noop_for(fn):
        def held(*args, **kw):
            calls.append((fn, args, kw))
        held._mon1_held = True
        return held
    for dd in _MP_DICTS:
        for key in ('masked_products', 'masked_product'):
            fn = dd.get(key)
            if isinstance(fn, types.FunctionType) and not getattr(fn, '_mon1_held', False):
                saved.append((dd, key, fn))
                dd[key] = noop_for(fn)
    try:
        d['q_post'](b, q, payloads, H, 1, parent, reverse, phase)
    finally:
        for dd, key, fn in saved:
            dd[key] = fn
    return calls


def _check(d, kind, f, g, payloads, parent, reverse, phase):
    if kind is None:
        return 'unknown base'
    if parent is None or len(parent) != 3:
        return 'no parent'
    if len(f) != 2 or len(g) != 1:
        return 'not level 1'
    if len(payloads) != 2 or any(len(v) != 1 for v in payloads) or len(parent[0]) != 2:
        return 'payload shape'
    words = [g[0], f[-1], parent[1], parent[2]] + [v[0] for v in payloads] + list(parent[0])
    if any(len(w) != 12 for w in words):
        return 'word width'
    kc = d.get('KEEP_CHILD')                     # where bound; the real condition is the phase check below: level 2
    if kc is not None and 2 not in kc:           # hands F in given_gap[2] exactly when it keeps the child's flag
        return 'KEEP_CHILD'
    if not isinstance(phase, dict) or 2 not in phase.get('given_gap', {}) or 'a' not in phase:
        return 'no given gap'
    if reverse and 1 not in phase['a']:
        return 'no forward token'
    return None


def _tail(kind, gl, d, b, a, carry, phase):
    """The reverse base's last part without its three scalings."""
    if kind == 'a08':
        Z, saved = gl['prefix'](b, a)
        ai = gl['finish'](b, a, saved)
        d['consume'](b, ai, Z, carry)
        gl['erase_normpair'](b, a, ai, Z, saved, (phase['returned_gap'][2], phase))
        return
    ai, Z = d['inverse_zero'](b, a)
    d['consume'](b, ai, Z, carry)
    F = phase['returned_gap'][2]
    Z0 = b.alloc_qubits(1)[0]
    b.init_and(F, Z, Z0)
    d['erase_product_zero'](b, a, ai, Z)
    phase.setdefault('_undo_parent_zeros', {})[2] = Z0
    mask = b.alloc_bits(12)
    b.hmr(ai, mask)
    b.free(ai)
    d['phase_inverse_for_gap'](b, a, mask)
    b.free(mask)


def _make(d, orig):
    kind = _kind(orig)
    gl = orig.__globals__

    def base_case(b, f, g, payloads, H, parent, reverse=False, phase=None):
        tag = 'rev' if reverse else 'fwd'
        why = _check(d, kind, f, g, payloads, parent, reverse, phase)
        if reverse:
            fwd = isinstance(phase, dict) and phase.get('_mon1', False)
            if fwd and why:
                raise AssertionError('dqi_int_mon1: forward monic, reverse refused (%s)' % why)
            if not fwd:
                _refuse('rev: forward not monic' if why is None else 'rev: ' + why)
                return orig(b, f, g, payloads, H, parent, reverse, phase)
        elif why:
            _refuse('fwd: ' + why)
            return orig(b, f, g, payloads, H, parent, reverse, phase)
        a = g[0]
        q = f[-1]
        tops, p, c = parent
        runops = d['runops']
        if not reverse:
            calls = _call_q_post(d, b, q, payloads, H, parent, False, phase)
            assert len(calls) == 1, ('dqi_int_mon1: forward q_post held %d masked calls' % len(calls))
            F = phase['returned_gap'][2]
            fn, args, kw = calls[0]
            assert args[1] is F and list(args[4]) == list(tops), 'dqi_int_mon1: unexpected masked call'
            ar = d.get('ar') or fn.__globals__['ar']

            def action(ai):
                cscale(b, ar, F, tops, ai, a)
                fn(*args, **kw)
                runops(b, [('s', v[0], True) for v in payloads], a, ai, False)      # v0 /= a (the kit's own form)
            phase['a'][1] = d['inverse_apply'](b, a, action, None, True)
            phase['_mon1'] = (fn, ar)
        else:
            F = phase['given_gap'][2]
            mp, ar = phase['_mon1']

            def action(ai):
                runops(b, [('s', v[0], True) for v in payloads], a, ai, True)       # v0 *= a (undo of v0 /= a)
                mp(b, F, payloads, c, tops)
                cscale(b, ar, F, tops, a, ai)
            fns = _context_fns(d) if (MERGE and kind == 'a08') else None
            if fns is not None:
                iz, cons, er = fns
                guards = _late_guards(d)
                n0 = [len(x) for x in guards]
                ai, Zf = iz(b, a)
                cons(b, ai, Zf, phase['a'].pop(1))
                b.cx(Zf, a[0])
                b.cx(Zf, ai[0])
                action(ai)
                b.cx(Zf, ai[0])
                b.cx(Zf, a[0])
                Zp = b.alloc_qubits(1)[0]
                b.init_and(F, Zf, Zp)                       # the level-2 hand-off F AND [a = 0], as a08's tail makes it
                phase.setdefault('_undo_parent_zeros', {})[2] = Zp
                er(b, a, ai, Zf, False)                     # the context's own phase paid in place: no carry
                assert [len(x) for x in guards] == n0, 'dqi_int_mon1: late paired work queued inside the merged context'
                carry = None
                LOG['merged'] = LOG.get('merged', 0) + 1
            else:
                carry = d['inverse_apply'](b, a, action, phase['a'].pop(1), True)
            calls = _call_q_post(d, b, q, payloads, H, parent, True, phase)
            assert len(calls) == 1, ('dqi_int_mon1: reverse q_post held %d masked calls' % len(calls))
            assert calls[0][1][1] is F, 'dqi_int_mon1: reverse q_post used another gap flag'
            assert calls[0][0].__code__ is mp.__code__, 'dqi_int_mon1: reverse masked products differ from the forward'
            if fns is not None:
                pass
            elif NOTAIL:          # count only: a fresh zero qubit stands in for the hand-off (NOT exact)
                LOG['notail'] = LOG.get('notail', 0) + 1
                phase.setdefault('_undo_parent_zeros', {})[2] = b.alloc_qubits(1)[0]
            else:
                _tail(kind, gl, d, b, a, carry, phase)
            phase.pop('_mon1', None)
        LOG['calls'][tag] += 1
        LOG['captured'][tag] += 1
        LOG['tails'][kind] = LOG['tails'].get(kind, 0) + (1 if reverse else 0)
        _dump()
    base_case._dqi_int_mon1 = True
    return base_case


_A08 = {}          # (id(scope), id(original)) -> the monic base made by _make(scope, rebuilt original a08 base_case)
_TRAMP_SRC = ("def _outer():\n"
              "    original = scope = None\n"
              "    def base_case(b, f, g, payloads, H, parent, reverse=False, phase=None):\n"
              "        return __import__('dqi_int_mon1').dispatch(scope, original, b, f, g, payloads, H, parent,\n"
              "                                                   reverse, phase)\n"
              "    return base_case\n")


def dispatch(scope, original, b, f, g, payloads, H, parent, reverse=False, phase=None):
    """Called by every patched a08 base_case (v3): the monic base of that closure's own (scope, original)."""
    return _A08[(id(scope), id(original))](b, f, g, payloads, H, parent, reverse, phase)


def _is_a08(fn):
    co = getattr(fn, '__code__', None)
    return (isinstance(fn, types.FunctionType) and co is not None and co.co_name == 'base_case'
            and co.co_filename.replace('\\', '/').endswith('a08/scalar_compact_v3.py')
            and set(co.co_freevars) == {'scope', 'original'} and fn.__closure__ is not None)


def install(km, ns):
    """v3 (fix round 1, Nov 19 00:40Z): patch by CODE IDENTITY, as tools/dqi_int_zcache.py and dqi_int_fdots.py do.
    The line reaches level 1 through a08's base_case CLOSURE (free variables scope, original) held in a dict that no
    walk by name reaches (v1 wrapped 160 copies of dqi_zn_pdl2's base_case and was never called; v2 found 0 scopes).
    Every live function object with a08's base_case code gets a trampoline code with the same free variables; it calls
    dispatch(scope, original, ...), i.e. _make(scope, a rebuilt copy of the original), and the copy serves refusals."""
    if OFF:
        return dict(LOG)
    import gc
    fns = [o for o in gc.get_objects() if isinstance(o, types.FunctionType)]
    walk = sys.modules.get('walk_bound')
    if walk is not None:
        fns += [fn for fn in walk.functions(ns.decode) if isinstance(fn, types.FunctionType)]
    seen_mp, seen_fn = set(), set()
    for fn in fns:
        for dd in [fn.__globals__] + [c for c in _cells(fn) if isinstance(c, dict)]:
            if ('masked_products' in dd or 'masked_product' in dd) and id(dd) not in seen_mp:
                seen_mp.add(id(dd))
                _MP_DICTS.append(dd)
    tramp = [c for c in compile(_TRAMP_SRC, __file__, 'exec').co_consts if isinstance(c, types.CodeType)][0]
    tramp = [c for c in tramp.co_consts if isinstance(c, types.CodeType)][0]
    patched = 0
    for fn in fns:
        if id(fn) in seen_fn or not _is_a08(fn):
            continue
        seen_fn.add(id(fn))
        co = fn.__code__
        if tramp.co_freevars != co.co_freevars or tramp.co_argcount != co.co_argcount:
            LOG['refused']['trampoline layout'] = LOG['refused'].get('trampoline layout', 0) + 1
            continue
        cells = dict(zip(co.co_freevars, fn.__closure__))
        scope, original = cells['scope'].cell_contents, cells['original'].cell_contents
        key = (id(scope), id(original))
        if key not in _A08:
            copy = types.FunctionType(co, fn.__globals__, 'base_case', fn.__defaults__, fn.__closure__)
            _A08[key] = _make(scope, copy)
            LOG['originals'].append('%s:%s:%d' % (os.path.basename(co.co_filename), co.co_name, co.co_firstlineno))
        fn.__code__ = tramp
        patched += 1
    LOG.update(installed=True, scopes=len(_A08), patched=patched, mp_dicts=len(_MP_DICTS))
    _dump()
    return dict(installed=True, scopes=len(_A08), patched=patched, mp_dicts=len(_MP_DICTS), originals=LOG['originals'])


def _test(root):
    """cscale on every a != 0 and F in {0, 1} (8,190 shots, random tops), forward (x = a^-1, y = a) and reverse form
    (x = a, y = a^-1): values, clean ancillas, no phase flip, static Toffoli-like count 120."""
    import importlib.util
    import numpy as np
    sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
    import kickmix as km
    sys.modules['circ'] = km
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, here)
    from dqi_x12c_gf import mul, inv0
    sys.path.insert(0, os.path.join(root, 'playbook/kit/decoder'))
    spec = importlib.util.spec_from_file_location('mon1_kit_field',
                                                  os.path.join(root, 'playbook/kit/decoder/exp/arithmetic/field.py'))
    ar = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ar)
    rng = np.random.default_rng(25)
    av = np.repeat(np.arange(1, 4096, dtype=np.int64), 2)
    Fv = np.tile(np.array([0, 1], dtype=np.int64), 4095)
    n = av.shape[0]
    t1 = rng.integers(0, 4096, size=n, dtype=np.int64)
    t2 = rng.integers(0, 4096, size=n, dtype=np.int64)
    aiv = inv0(av)
    assert (mul(av, aiv) == 1).all()
    ok = True
    for form in ('forward', 'reverse'):
        b = km.CircuitBuilder()
        F = b.alloc_qubits(1)[0]
        a = b.alloc_qubits(12)
        ai = b.alloc_qubits(12)
        T1 = b.alloc_qubits(12)
        T2 = b.alloc_qubits(12)
        if form == 'forward':
            cscale(b, ar, F, [T1, T2], ai, a)
        else:
            cscale(b, ar, F, [T1, T2], a, ai)
        c = b.finish_circuit()
        sim = km.Simulator(batch_size=n, debug_print=False, seed=7)
        sim.use_same_registers_as(c)
        words = [a, ai, T1, T2]
        vals = np.stack([av, aiv, t1, t2], axis=1)
        qa = km.array([q for w in words for q in w])
        bits = ((vals[:, :, None] >> np.arange(12)) & 1 == 1).transpose(1, 2, 0).reshape(-1, n)
        sim.write_across_shots(qa, bits)
        sim.write_across_shots(km.array([F]), (Fv == 1)[None, :])
        sim.do(c)
        out = np.asarray(sim.read_across_shots(qa)).reshape(4, 12, n).astype(np.int64)
        outw = (out * (1 << np.arange(12))[None, :, None]).sum(axis=1).T
        fout = np.asarray(sim.read_across_shots(km.array([F]))).reshape(-1).astype(np.int64)
        s = np.where(Fv == 1, aiv if form == 'forward' else av, 1)
        want = np.stack([av, aiv, mul(t1, s), mul(t2, s)], axis=1)
        wrong = int((outw != want).any(axis=1).sum()) + int((fout != Fv).sum())
        ph = int(np.asarray(sim.read_phase_flipped_across_shots(), dtype=bool).sum())
        anc = km.array([km.q(i) for i in range(49, c.num_qubits)])
        dirty = int(np.asarray(sim.read_across_shots(anc)).any(axis=0).sum()) if c.num_qubits > 49 else 0
        tof = c.operation_statistics().num_toffoli_like
        good = wrong == 0 and ph == 0 and dirty == 0 and tof == 120 and c.num_qubits <= 49 + 12
        ok &= good
        print('%s form: %d shots (every a != 0, F = 0 / 1), wrong %d, phase flips %d, dirty ancillas %d, static '
              'Toffoli-like %d, qubits %d (+%d): %s' % (form, n, wrong, ph, dirty, tof, c.num_qubits,
                                                        c.num_qubits - 49, 'PASS' if good else 'FAIL'), flush=True)
    print('TEST %s' % ('PASS' if ok else 'FAIL'), flush=True)
    return ok


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == 'test':
        sys.exit(0 if _test(sys.argv[2]) else 1)
    print(__doc__)
