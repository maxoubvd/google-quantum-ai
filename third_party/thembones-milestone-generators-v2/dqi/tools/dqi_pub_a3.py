"""dqi_pub_a3.py -- builder D1: install the projective finish (transplant a3) into our decode.

Mechanism credit: the projective finish of the public DQI decoder (Syed Asad Ali, Apache-2.0); our own code on
our own kit kernels (tools/dqi_pub_finish.py).

What changes in the decode (both directions of the converter, and the scan's classical constant):
  converter core (dot_cache_core, square-root domain), forward, after the Euclid and the swaps:
      top ^= sum_{j in J} B'_{j+1} X_j        (17 products: the Laurent top X~_34; measured rule, ws_dqi_pub.md)
      A[0..34] := K' = M' C' mod B'           (finish: 17 x multiply_forward99, free maps)
    instead of: delta (16 products, kept live through the scan), a11's Hankel (4,641; one word kept through
    the scan), top ^= delta, rotation.  Reverse mirrors it (17 x divide, the same 17 products).
  scan: node_polynomial P -> P * M (M = coefficient-wise square of M'), so the classical mask constant
    1 / (x P(x)) becomes 1 / (x P(x) M(x)): at every root x of B, K(x) / M(x) = C(x).
install(km, ns) with ns.decode = the decode function; returns a log dict.

Environment (defaults = the first verified form, v2 with divide_lowwidth369, p1a3b):
  DQI_PUB_A3_CORE  v2 (placeholders, top word both ways) | v3 (delta live) | v4 (node 0 held; converter test only) |
                   v5 (delta live, taken out of K' before the divisions: the best, -2,997 executed)
  DQI_PUB_A3_DIV   369 | 372 (kit divide_lowwidth*) | dr (B6's tools/dqi_kern_div.divide_direct: 315, width 36)
  DQI_PUB_A3_SYNC  1 = free-stack synchronisation: a scratch build of the unpatched decode inside install() records
                   the free stack at the forward core's exit and at the reverse Euclid's entry, and the patched
                   build reproduces both (kickmix reuses freed qubits last-in first-out; without it v5 fails the
                   decode's final routing). Needed for v5.
  DQI_PUB_A3_EXACT 1 = placeholder words on exactly the wires the unpatched decode held (with SYNC; untested).
Verified: v5 + dr + SYNC on q902m = p1m5b (362,383.0, sha 66e5f1f0), research/ws_dqi_pub.md.
"""
import os
import sys

TOOLS = os.path.dirname(os.path.abspath(__file__))
if TOOLS not in sys.path:
    sys.path.insert(0, TOOLS)
import dqi_pub_finish as FIN                                          # noqa: E402

J = tuple(range(0, 34, 2))   # measured: X~_34 = top + sum_i B'_{2i+1} X_{2i} (dqi_pub_topfit.py, 05:46Z)

CORE_SRC_DELTA = '''
def dot_cache_core(b,A,B,reverse=False,history_scratch12=None,phase_token=None,*,monic_root=False):
    assert len(A)==len(B)==35
    phase={'a':{},'q':{},'_short_monic':bool(monic_root)} if phase_token is None else phase_token
    assert not reverse or phase_token is not None
    old.frob(b,B,2048);old.frob(b,A,2048)
    E=B[::2];O=B[1::2];u=A[1::2];w=A[:34:2];top=A[-1]
    def bez(rev):
        if phase['_short_monic']:_bez_monic(b,O,E,[u,w],rev,history_scratch12,phase)
        else:old.bez_lookahead(b,O,E,[u,w],rev,history_scratch12,phase)
    if not reverse:
        for i in range(17):old.ar.mul_xor0(b,top,O[i],w[i])
        bez(False)
        # delta on the old delta's wires (allocated where the old core allocated it), live through the scan
        delta=b.alloc_qubits(12)
        for i in range(17):b.swap(A[2*i],A[2*i+1])
        for j in _D1_J:old.ar.mul_xor0(b,B[j+1],A[j],delta)
        b.cx(delta,top)
        phase['_d1_finish']=_D1_FIN.forward(b,B,A)
        # one clean word on the wires of a11's saved word (layout parity for later allocations)
        phase['_d1_delta']=delta
        phase['_d1_dummy']=(_D1_DUMMY(b,12),)
    else:
        for _w in phase.pop('_d1_dummy'):b.free(_w)
        _D1_FIN.reverse(b,B,A,phase.pop('_d1_finish'))
        delta=phase.pop('_d1_delta');b.cx(delta,top)
        mask=b.alloc_bits(12);b.hmr(delta,mask);b.free(delta)
        for j in _D1_J:old.ar.phase_product(b,B[j+1],A[j],mask)
        b.free(mask)
        for i in range(17):b.swap(A[2*i],A[2*i+1])
        bez(True)
        for i in range(17):old.ar.mul_xor0(b,top,O[i],w[i])
    old.frob(b,A,2);old.frob(b,B,2)
    if reverse:assert not phase['a'] and not phase['q']
    return phase if not reverse else None
'''

CORE_SRC = '''
def dot_cache_core(b,A,B,reverse=False,history_scratch12=None,phase_token=None,*,monic_root=False):
    assert len(A)==len(B)==35
    phase={'a':{},'q':{},'_short_monic':bool(monic_root)} if phase_token is None else phase_token
    assert not reverse or phase_token is not None
    old.frob(b,B,2048);old.frob(b,A,2048)
    E=B[::2];O=B[1::2];u=A[1::2];w=A[:34:2];top=A[-1]
    def bez(rev):
        if phase['_short_monic']:_bez_monic(b,O,E,[u,w],rev,history_scratch12,phase)
        else:old.bez_lookahead(b,O,E,[u,w],rev,history_scratch12,phase)
    if not reverse:
        for i in range(17):old.ar.mul_xor0(b,top,O[i],w[i])
        bez(False)
        for i in range(17):b.swap(A[2*i],A[2*i+1])
        for j in _D1_J:old.ar.mul_xor0(b,B[j+1],A[j],top)
        phase['_d1_finish']=_D1_FIN.forward(b,B,A)
        # Two clean words on the wires the old core's delta and a11's saved word held through the scan, so every
        # later allocation (the parked flags' rebuilds, the final routing) meets the same free wires as before.
        _dd=_D1_DUMMY(b,24);phase['_d1_dummy']=(_dd[:12],_dd[12:])
    else:
        for _w in phase.pop('_d1_dummy'):b.free(_w)
        _D1_FIN.reverse(b,B,A,phase.pop('_d1_finish'))
        for j in _D1_J:old.ar.mul_xor0(b,B[j+1],A[j],top)
        for i in range(17):b.swap(A[2*i],A[2*i+1])
        bez(True)
        for i in range(17):old.ar.mul_xor0(b,top,O[i],w[i])
    old.frob(b,A,2);old.frob(b,B,2)
    if reverse:assert not phase['a'] and not phase['q']
    return phase if not reverse else None
'''

CORE_SRC_HOLD = '''
def dot_cache_core(b,A,B,reverse=False,history_scratch12=None,phase_token=None,*,monic_root=False):
    assert len(A)==len(B)==35
    phase={'a':{},'q':{},'_short_monic':bool(monic_root)} if phase_token is None else phase_token
    assert not reverse or phase_token is not None
    old.frob(b,B,2048);old.frob(b,A,2048)
    E=B[::2];O=B[1::2];u=A[1::2];w=A[:34:2];top=A[-1]
    def bez(rev):
        if phase['_short_monic']:_bez_monic(b,O,E,[u,w],rev,history_scratch12,phase)
        else:old.bez_lookahead(b,O,E,[u,w],rev,history_scratch12,phase)
    if not reverse:
        for i in range(17):old.ar.mul_xor0(b,top,O[i],w[i])
        bez(False)
        for i in range(17):b.swap(A[2*i],A[2*i+1])
        for j in _D1_J:old.ar.mul_xor0(b,B[j+1],A[j],top)
        # node 0's old residue is kept in the first 24 wires allocated after the Euclid: the wires of the old
        # core's delta and a11's saved word (layout parity), so no placeholder words are needed
        phase['_d1_finish']=_D1_FIN.forward(b,B,A,hold=True)
    else:
        _masks,_held=phase.pop('_d1_finish')
        _D1_FIN.reverse_hold(b,B,A,_masks,_held)
        for j in _D1_J:old.ar.mul_xor0(b,B[j+1],A[j],top)
        for i in range(17):b.swap(A[2*i],A[2*i+1])
        bez(True)
        for i in range(17):old.ar.mul_xor0(b,top,O[i],w[i])
    old.frob(b,A,2);old.frob(b,B,2)
    if reverse:assert not phase['a'] and not phase['q']
    return phase if not reverse else None
'''

CORE_SRC_V5 = '''
def dot_cache_core(b,A,B,reverse=False,history_scratch12=None,phase_token=None,*,monic_root=False):
    assert len(A)==len(B)==35
    phase={'a':{},'q':{},'_short_monic':bool(monic_root)} if phase_token is None else phase_token
    assert not reverse or phase_token is not None
    old.frob(b,B,2048);old.frob(b,A,2048)
    E=B[::2];O=B[1::2];u=A[1::2];w=A[:34:2];top=A[-1]
    def bez(rev):
        if phase['_short_monic']:_bez_monic(b,O,E,[u,w],rev,history_scratch12,phase)
        else:old.bez_lookahead(b,O,E,[u,w],rev,history_scratch12,phase)
    if not reverse:
        for i in range(17):old.ar.mul_xor0(b,top,O[i],w[i])
        bez(False)
        delta=b.alloc_qubits(12)
        for i in range(17):b.swap(A[2*i],A[2*i+1])
        for j in _D1_J:old.ar.mul_xor0(b,B[j+1],A[j],delta)
        b.cx(delta,top)
        phase['_d1_finish']=_D1_FIN.forward(b,B,A)
        phase['_d1_delta']=delta
        phase['_d1_dummy']=(_D1_DUMMY(b,12),)
    else:
        for _w in phase.pop('_d1_dummy'):b.free(_w)
        # delta enters K' only as delta * M' (the Laurent top's column): take it out of K' for free, X-measure it
        # before the divisions (they then have the width of B6's dr), repay its phase once X is back
        delta=phase.pop('_d1_delta')
        _D1_ADDM(b,delta,A)
        mask=b.alloc_bits(12);b.hmr(delta,mask);b.free(delta)
        _D1_FIN.reverse(b,B,A,phase.pop('_d1_finish'))
        for j in _D1_J:old.ar.phase_product(b,B[j+1],A[j],mask)
        b.free(mask)
        for i in range(17):b.swap(A[2*i],A[2*i+1])
        bez(True)
        for i in range(17):old.ar.mul_xor0(b,top,O[i],w[i])
    old.frob(b,A,2);old.frob(b,B,2)
    if reverse:assert not phase['a'] and not phase['q']
    return phase if not reverse else None
'''

LOG = {}


def core_source():
    """DQI_PUB_A3_CORE: 'v2' (placeholders, top word both ways; default), 'v3' (delta live), 'v4' (node 0 held),
    'v5' (delta live, taken out of K' in the reverse before the divisions)."""
    v = os.environ.get('DQI_PUB_A3_CORE') or ('v3' if os.environ.get('DQI_PUB_A3_DELTA', '0') == '1' else 'v2')
    return v, {'v2': CORE_SRC, 'v3': CORE_SRC_DELTA, 'v4': CORE_SRC_HOLD, 'v5': CORE_SRC_V5}[v]


def make_addm(fin, cx):
    """words[k] ^= M'_k * delta for k = 0..34 (M' monic of degree 34): free (constant multiplications)."""
    M = fin.geo['M']

    def addm(b, delta, words):
        assert len(words) == 35
        for k in range(35):
            if M[k]:
                cx(b, delta, words[k], M[k])
    return addm


def kit_modules():
    found = {}
    for m in list(sys.modules.values()):
        d = getattr(m, '__dict__', {})
        if 'multiply_forward99' in d and 'divide_lowwidth369' in d and 'divide_lowwidth372' in d:
            found['sched'] = m
        if 'linear_words' in d and 'matrix_ops' in d:
            found['common'] = m
        if 'extension_columns' in d and 'mul_xor0' in d:
            found['ar'] = m
    assert set(found) == {'sched', 'common', 'ar'}, found.keys()
    return found


def divider_globals(decode):
    """The un-synthesis fold divider's globals (as B6's tools/dqi_kern_divtest.py finds them)."""
    walk = sys.modules['walk_bound']
    owner = next(f for f in walk.functions(decode)
                 if f.__code__.co_filename.endswith('decode_inverse_retained.py') and 'produce_retained' in f.__globals__)
    return owner.__globals__


def make_finish(km, decode=None):
    kit = kit_modules()
    gf = km.GF2Field(12)
    KX = kit['ar'].extension_columns()[2]
    name = os.environ.get('DQI_PUB_A3_DIV', '369')
    if name == 'dr':                     # B6's direct division (tools/dqi_kern_div.py): 315 executed, width 36
        import dqi_kern_div as KD
        G = divider_globals(decode)

        def div(b, x, y, K):
            return KD.divide_direct(b, x, y, K, G)
        div.__name__ = 'dqi_kern_div.divide_direct'
    elif name in ('dt', 'dt42', 'tw2', 'tw3'):  # B6's other unflagged forms: dt 294.93 / 48, dt42 309.89 / 42
        import dqi_kern_div as KD
        G = divider_globals(decode)
        fn = KD.VARIANTS[name]

        def div(b, x, y, K):
            return fn(b, x, y, K, G)
        div.__name__ = 'dqi_kern_div.' + name
    else:
        div = getattr(kit['sched'], 'divide_lowwidth' + name)
    return FIN.Finish(km, gf, KX, kit['common'].linear_words, kit['sched'].multiply_forward99, div,
                      getattr(kit['sched'], 'mul24_xor99')), gf


def converter_cells(decode):
    """The converter closures that hold the core (cell 'core' whose function is named dot_cache_core)."""
    walk = sys.modules['walk_bound']
    out = []
    for fn in walk.functions(decode):
        if fn.__closure__ and 'core' in fn.__code__.co_freevars:
            cells = dict(zip(fn.__code__.co_freevars, fn.__closure__))
            core = cells['core'].cell_contents
            if getattr(core, '__name__', '') == 'dot_cache_core':
                out.append(cells['core'])
    return out


FINISH = {}

# ---- free-stack synchronisation (layout parity with the unpatched decode) -------------------------------------
# kickmix hands freed qubits back last-in first-out, so the ORDER of frees decides where every later allocation
# lands, and the decode's flag rebuilds (alloc_pref fallbacks) and its final routing depend on it. A scratch build
# of the unpatched decode records the free stack at two points (S1: exit of the forward converter core; S2: entry
# of the reverse Euclid); the patched build reproduces exactly those stacks there (all spare wires taken and given
# back in the recorded order: no gate). Everything after S2 is then the unpatched decode's own allocation sequence.
STACK = dict(mode=None, s1=[], s2=[], i1=0, i2=0, fixes=0, misses=0)


def _qidx(q):
    return int(str(q)[1:])


def _spare(b):
    n = b.finish_circuit().num_qubits - b.num_allocated_qubits
    return list(b.alloc_qubits(n)) if n > 0 else []


def _snapshot(b):
    got = _spare(b)
    order = [_qidx(q) for q in got]                  # stack top first
    for q in reversed(got):
        b.free(q)
    return order


def _restore(b, order):
    got = _spare(b)
    have = {_qidx(q): q for q in got}
    if set(have) != set(order):
        STACK['misses'] += 1
        for q in reversed(got):
            b.free(q)
        return
    for i in reversed(order):
        b.free(have[i])
    STACK['fixes'] += 1


def _sync(b, key):
    if STACK['mode'] == 'capture':
        STACK[key].append(_snapshot(b))
    elif STACK['mode'] == 'replay':
        i = 'i1' if key == 's1' else 'i2'
        if STACK[i] < len(STACK[key]):
            _restore(b, STACK[key][STACK[i]])
        STACK[i] += 1


def _dummy(b, n):
    """Placeholder of n clean qubits. In replay (DQI_PUB_A3_EXACT=1): exactly the wires the unpatched decode held
    at its forward core's exit beyond what is live here (its delta / saved word), with the recorded free stack
    restored (this is then the S1 synchronisation)."""
    km = sys.modules['circ']
    if (STACK['mode'] == 'replay' and os.environ.get('DQI_PUB_A3_EXACT', '0') == '1'
            and STACK['i1'] < len(STACK['s1'])):
        rec = STACK['s1'][STACK['i1']]
        recset = set(rec)
        got = _spare(b)
        extra = [q for q in got if _qidx(q) not in recset]
        have = {_qidx(q): q for q in got if _qidx(q) in recset}
        if len(extra) == n and set(have) == recset:
            for i in reversed(rec):
                b.free(have[i])
            STACK['i1'] += 1
            STACK['fixes'] += 1
            STACK['exact'] = STACK.get('exact', 0) + 1
            return km.array(extra)
        for q in reversed(got):
            b.free(q)
        STACK['misses'] += 1
    return b.alloc_qubits(n)


def _wrap_core(core):
    def dot_cache_core(b, A, B, reverse=False, history_scratch12=None, phase_token=None, *, monic_root=False):
        i1 = STACK['i1']
        out = core(b, A, B, reverse, history_scratch12, phase_token, monic_root=monic_root)
        if not reverse and (STACK['mode'] == 'capture' or STACK['i1'] == i1):
            _sync(b, 's1')            # capture; or replay when the core's placeholder did not sync already
        return out
    dot_cache_core.__wrapped__ = core
    return dot_cache_core


def _wrap_bez(fn):
    if getattr(fn, '_d1_wrapped', False):
        return fn

    def _bez_monic(b, O, E, payloads, rev, *args, **kw):
        if rev:
            _sync(b, 's2')
        return fn(b, O, E, payloads, rev, *args, **kw)
    _bez_monic._d1_wrapped = True
    return _bez_monic


def _capture(km, decode, cells):
    """Scratch build of the unpatched decode, recording S1 / S2."""
    scopes = {}
    for cell in cells:
        core = cell.cell_contents
        scopes[id(core.__globals__)] = core.__globals__
    for scope in scopes.values():
        scope['_bez_monic'] = _wrap_bez(scope['_bez_monic'])
        old = scope.get('old')                    # the general (non-monic) root path: same S2 point
        if old is not None and hasattr(old, 'bez_lookahead'):
            old.bez_lookahead = _wrap_bez(old.bez_lookahead)
    originals = [cell.cell_contents for cell in cells]
    for cell in cells:
        cell.cell_contents = _wrap_core(cell.cell_contents)
    STACK.update(mode='capture', s1=[], s2=[])
    # a once-per-process INFO_TAG wrapper (the re-roll patch r6tag) must not be spent on the scratch build
    CB = km.CircuitBuilder
    tagged = CB.mark
    fv = getattr(getattr(tagged, '__code__', None), 'co_freevars', ())
    if 'done' in fv and 'original' in fv:
        CB.mark = dict(zip(fv, tagged.__closure__))['original'].cell_contents
    try:
        b0 = CB()
        S0 = [b0.alloc_qubits(12) for _ in range(70)]
        E0 = [b0.alloc_bits(12) for _ in range(4095)]
        decode(b0, S0, E0)
    finally:
        CB.mark = tagged
    del b0, S0, E0
    import gc
    gc.collect()
    for cell, core in zip(cells, originals):
        cell.cell_contents = core
    STACK.update(mode='replay', i1=0, i2=0)


def install(km, ns):
    if km is None:
        km = sys.modules['circ']
    fin, gf = make_finish(km, ns.decode)
    FINISH.update(fin=fin, gf=gf, M=fin.scan_factor(gf))
    cells = converter_cells(ns.decode)
    assert cells, 'no converter core found'
    sync = os.environ.get('DQI_PUB_A3_SYNC', '0') == '1'
    if sync:
        _capture(km, ns.decode, cells)
    done = set()
    for cell in cells:
        core = cell.cell_contents
        if id(core) in done:
            continue
        scope = core.__globals__
        scope['_D1_J'] = J
        scope['_D1_FIN'] = fin
        scope['_D1_ADDM'] = make_addm(fin, kit_modules()['common'].cx)
        scope['_D1_DUMMY'] = _dummy
        _, src = core_source()
        exec(compile(src, __file__, 'exec'), scope)
        new = scope['dot_cache_core']
        if sync:
            new = _wrap_core(new)
        for other in cells:
            if other.cell_contents is core:
                other.cell_contents = new
        done.add(id(core))
        done.add(id(new))
    g = ns.decode.__globals__['g']
    P = list(g['Pscan'])
    M = fin.scan_factor(gf)
    g['Pscan'] = FIN.pmul(FIN.Field(gf), P, M)
    LOG.update(cores=len(done) // 2, cells=len(cells), Pscan_degree=len(g['Pscan']) - 1, KX=fin.KX,
               division=fin.divide.__name__, core=core_source()[0], sync=sync,
               captured=(len(STACK['s1']), len(STACK['s2'])))
    return dict(LOG)
