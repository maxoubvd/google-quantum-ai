#!/usr/bin/env python3
"""dqi_zn_kpark.py -- builder Z1 (Nov 14, fourth holder): the 879th qubit of Z2's un-fused root (kv_unfused4) back.
Measured (tools/dqi_zn_flagidle.py, box_out/fi_u4.rep): every 879 moment of kv_unfused4 outside a14's three-frame
blocks has the level-15 held flag K = phase['held_gap'][15] (= G at the level-16 post) live and IDLE. K is parked
(X-measured, outcome kept) over six intervals and rebuilt = [w == 0] on its own wire (a zero test, 11 Toffolis):
  P1 forward, inside root_child_post's quotient_inverse after its last use (B5's flagpark part c kernel: the CCZ
     (G, K, m0) deferred under a kept bit) -> rebuilt at hook f1 (level-16 post, before b.ccx(F,G,a[0]));
  P2 forward, hook f2 (after consume_clipped) -> f3 (before b.ccx(F,G,g[-2][0]));
  P3 forward, hook f4 (end of the level-16 post) -> r1 (reverse level-16 post, before b.ccx(F,G,a[0])): covers the
     root's return, dot_cache_core and the reverse's division / spare moments;
  R1 reverse, hook r2 (after consume_clipped) -> r3 (before b.ccx(F,G,a[0]); out=retire_clipped);
  R2 reverse, inside quotient_inverse (part c kernel) -> close_quotient_phase (part c's rebuild).
The a14 three-frame moments (reverse levels 12 and 15) are flagpark part a's: run with ZN_FP=a.
Words: env ZN_KP_W_<hook> (python expression in f, g, CH), default ZN_KP_W or "CH['B'][28]".
Diagnostic: ZN_KP_DIAG=1 parks nothing; at every hook K ^ [w == 0] is computed for candidate words into kept
qubits; `python dqi_zn_kpark.py diag <seed>` builds the round trip and prints which candidates are 0 on every shot.
Install: ZN_LCT_P=kpark (after the entry, like every dqi_zn_ module).
"""
import os
import sys

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)

LOG = dict(posts=0, qi=0, parks=0, rebuilds=0, moved=0, diag=0)
TOK = {}
DIAG = []                      # (hook, candidate label, qubit)
DIAG_MODE = os.environ.get('ZN_KP_DIAG') == '1'
FWD_REBUILD_ONLY = os.environ.get('ZN_KP_PARTS', 'P1,FR,R2,RF').split(',')
FR = []
RF = []


def _fp():
    return sys.modules['dqi_int_flagpark']


def _word(hook, f, g):
    expr = os.environ.get('ZN_KP_W_' + hook, os.environ.get('ZN_KP_W', "CH['B'][28]"))
    return eval(expr, {}, dict(f=f, g=g, CH=_fp().CH))


def _active(m, phase, G):
    return m == 16 and phase is not None and phase.get('held_gap', {}).get(15) is G


def _candidates(f, g):
    CH = _fp().CH
    out = [('f[%d]' % i, f[i]) for i in range(-1, -6, -1)] + [('g[%d]' % i, g[i]) for i in range(-1, -5, -1)]
    out += [("CH['B'][%d]" % i, CH['B'][i]) for i in range(len(CH['B']))]
    return out


def _diag(b, hook, f, g, G):
    FP = _fp()
    for lab, w in _candidates(f, g):
        D = b.alloc_qubits(1)[0]
        FP.zero_toggle(b, w, D)
        b.cx(G, D)
        DIAG.append((hook, lab, D))
    LOG['diag'] += 1


def _rebuild(b, token, word, phase):
    """K = [word == 0] again. On its original wire if that wire is free; otherwise on a new wire, with the line's
    references re-pointed (held_gap[15], the root child's scalar context). Returns the flag object to use."""
    FP = _fp()
    k, orig = token
    w = FP._idx(orig)
    F = b.alloc_qubits(1)[0]
    FP.zero_toggle(b, word, F)
    b.cz(F, k[0])
    b.free(k)
    LOG['rebuilds'] += 1
    if FP._idx(F) == w:
        return orig
    held = []
    sp = b.alloc_qubits(1)[0]
    while FP._idx(sp) != w and len(held) < 64:
        held.append(sp)
        sp = b.alloc_qubits(1)[0]
    if FP._idx(sp) == w:
        b.swap(F, sp)
        b.free(F)
        for h in held:
            b.free(h)
        return orig
    held.append(sp)
    for h in held:
        b.free(h)
    LOG['moved'] += 1                                    # original wire taken: K lives on F from now on
    if phase.get('held_gap', {}).get(15) is orig:
        phase['held_gap'][15] = F
    ctx = phase.get('_child_scalar_context')
    if ctx is not None and ctx[2] is orig:
        phase['_child_scalar_context'] = (ctx[0], ctx[1], F)
    return F


def hook(b, name, m, phase, f, g, G):
    """Called from the rewritten level post; returns the flag object G to use from here on."""
    if not _active(m, phase, G):
        return G
    if DIAG_MODE:
        if name == 'f1':
            _diag(b, name, f, g, G)
        if name == 'f0':
            _diag(b, 'f0', f, g, phase['returned_gap'][17])
        if name == 'rq':
            _diag(b, 'rq', f, g, phase['given_gap'][17])
        return G
    FP = _fp()
    if name == 'f0' and 'FR' in FWD_REBUILD_ONLY:          # Froot waits through the level-16 post
        FR.append(FP.park(b, phase['returned_gap'][17]))
        LOG['parks'] += 1
        return G
    if name == 'r0' and 'RF' in FWD_REBUILD_ONLY:          # Froot (given_gap[17]) waits through the reverse post
        RF.append(FP.park(b, phase['given_gap'][17]))
        LOG['parks'] += 1
        return G
    if name == 'rq':
        if RF and RF[-1][1] is phase['given_gap'][17]:
            Fr = _rebuild(b, RF.pop(), eval(os.environ.get('ZN_KP_W_rq', 'f[-1]'), {}, dict(f=f, g=g, CH=_fp().CH)),
                          phase)
            if Fr is not phase['given_gap'][17]:
                phase['given_gap'][17] = Fr
        return G
    if name == 'f1' and 'P1' in FWD_REBUILD_ONLY:          # P1's rebuild (the part c kernel parked K in qi)
        if FP.QI_DEBT and FP.QI_DEBT[-1][0][1] is G:
            tK, d, F0, H0 = FP.QI_DEBT.pop()
            K = _rebuild(b, tK, _word('f1', f, g), phase)
            b.x(F0)
            b.push_condition(d[0])
            b.ccz(H0, F0, K)
            b.pop_condition()
            b.x(F0)
            b.free(d)
            return K
        return G
    if name in ('f2', 'r2') and ('P2' if name == 'f2' else 'R1') in FWD_REBUILD_ONLY:
        TOK.setdefault(name, {})[id(phase)] = FP.park(b, G)
        LOG['parks'] += 1
        return G
    if name in ('f3', 'r3'):
        src = 'f2' if name == 'f3' else 'r2'
        t = TOK.get(src, {}).pop(id(phase), None)
        return G if t is None else _rebuild(b, t, _word(name, f, g), phase)
    if name == 'f4' and 'P3' in FWD_REBUILD_ONLY:
        P3.append(FP.park(b, G))                          # a stack: forward / reverse converter calls nest
        LOG['parks'] += 1
        return G
    if name == 'r1':
        if not P3 or P3[-1][1] is not G:
            return G
        return _rebuild(b, P3.pop(), _word('r1', f, g), phase)
    return G


P3 = []


EDITS = [   # (line of decode_backhalf_q7.post, where, hook): forward branch indented 8, reverse branch 4
    ("        q_post(b,q1,payloads,H,m,parent,phase=phase)", 'after', 'f0'),
    ("        b.ccx(F,G,a[0]);ai=inverse_initialized(b,a)", 'before', 'f1'),
    ("    pivot(b,g,F);a=g[-1]", 'before', 'r0'),
    ("    pivot(b,g,F);q_post(b,q1,payloads,H,m,parent,True,phase)", 'before', 'rq'),
]


def rewrite_post(orig):
    import inspect
    import textwrap
    src = textwrap.dedent(inspect.getsource(orig))
    lines = src.split('\n')
    out = []
    used = set()
    for ln in lines:
        hit = [e for e in EDITS if ln == e[0]]
        if hit:
            text, where, name = hit[0]
            assert name not in used, ('edit twice', name)
            used.add(name)
            ind = text[:len(text) - len(text.lstrip())]
            call = ind + "G=_KP.hook(b,'%s',m,phase,f,g,G)" % name
            out += [call, ln] if where == 'before' else [ln, call]
        else:
            out.append(ln)
    assert used == {e[2] for e in EDITS}, ('edits missing', {e[2] for e in EDITS} - used)
    import types
    assert not orig.__code__.co_freevars, 'post has a closure'
    tmp = {}
    exec(compile('\n'.join(out), orig.__code__.co_filename, 'exec'), dict(orig.__globals__), tmp)
    orig.__globals__['_KP'] = sys.modules[__name__]
    new = types.FunctionType(tmp[orig.__name__].__code__, orig.__globals__, orig.__name__, orig.__defaults__)
    new.__kwdefaults__ = orig.__kwdefaults__
    new._zn_kpark = True                 # same live globals as the original: later installers stay visible
    return new


def install_qi(decode):
    """Part c's quotient_inverse (K parked inside, deferred CCZ) in every root_child_post scope; reverse rebuild in
    close_quotient_phase (part c's k_rebuild); forward rebuild at hook f1 (not part c's a22 root_pivot)."""
    FP = _fp()
    walk = sys.modules['walk_bound']
    n = 0
    for fn in walk.functions(decode):
        if fn.__name__ != 'root_child_post' or not fn.__code__.co_filename.endswith('ephemeral_root_block.py'):
            continue
        gl = fn.__globals__
        if gl.get('_zn_kpark_qi'):
            continue
        gl['_zn_kpark_qi'] = True
        if DIAG_MODE or 'P1' not in FWD_REBUILD_ONLY:
            continue
        gl['quotient_inverse'] = FP.make_quotient_inverse(gl['quotient_inverse'].__globals__)
        real = gl['close_quotient_phase']

        def close_quotient_phase(b, p, g, F, H, K, mask, _real=real):
            K = FP.k_rebuild(b, K, eval(os.environ.get('ZN_KP_W_cqp', os.environ.get('ZN_KP_W', "CH['B'][28]")),
                                        {}, dict(p=p, g=g, CH=FP.CH)))
            return _real(b, p, g, F, H, K, mask)

        gl['close_quotient_phase'] = close_quotient_phase
        n += 1
    LOG['qi'] = n
    return n


def install(km, ns):
    assert 'dqi_int_flagpark' in sys.modules, 'needs B5 flagpark (its chart record CH and kernels)'
    walk = sys.modules['walk_bound']
    fns = list(walk.functions(ns.decode))
    origs = [f for f in fns if f.__name__ == 'post' and f.__code__.co_filename.endswith('decode_backhalf_q7.py')
             and "G=phase['held_gap'][m-1]" in _src(f)]
    origs = list({id(f): f for f in origs}.values())
    assert origs, 'no decode_backhalf_q7 post'
    for orig in origs:
        new = rewrite_post(orig)
        n = 0
        for f in fns:
            for cell in (f.__closure__ or ()):
                try:
                    if cell.cell_contents is orig:
                        cell.cell_contents = new
                        n += 1
                except ValueError:
                    pass
            for k, v in list(f.__globals__.items()):
                if v is orig:
                    f.__globals__[k] = new
                    n += 1
        LOG['posts'] += n
    install_qi(ns.decode)
    install_fr(ns.decode)
    if 'RD' in FWD_REBUILD_ONLY and not DIAG_MODE:      # part c's reverse park of given_gap[17] (B5's install_d)
        FP = _fp()
        FP.PARTS.add('c')
        LOG['rd'] = FP.install_d(ns.decode)
    return dict(LOG)


def install_fr(decode):
    """Rebuild Froot (= [w == 0], default w = f[-1], the restored top) right after the original level-16 post
    returns inside a07z's direct_post wrapper (before the wrapper's root-pivot debt, which reads Froot)."""
    walk = sys.modules['walk_bound']
    n = 0
    for fn in walk.functions(decode):
        if fn.__name__ != 'direct_post' or not fn.__code__.co_filename.endswith('ephemeral_root_block.py'):
            continue
        cells = dict(zip(fn.__code__.co_freevars, fn.__closure__ or ()))
        if 'original' not in cells or getattr(cells['original'].cell_contents, '_zn_kpark', False):
            continue
        inner = cells['original'].cell_contents

        def original(b, f, g, payloads, F, H, parent, reverse=False, incoming=None, defer=False, phase=None,
                     _inner=inner):
            out = _inner(b, f, g, payloads, F, H, parent, reverse, incoming, defer, phase)
            if not reverse and len(f) - 1 == 16 and phase is not None:
                if DIAG_MODE:
                    _diag(b, 'fr', f, g, phase['returned_gap'][17])
                elif FR and FR[-1][1] is phase['returned_gap'][17]:
                    expr = os.environ.get('ZN_KP_W_fr', 'f[-1]')
                    Fr = _rebuild(b, FR.pop(), eval(expr, {}, dict(f=f, g=g, CH=_fp().CH)), phase)
                    if Fr is not phase['returned_gap'][17]:
                        phase['returned_gap'][17] = Fr
            return out
        original._zn_kpark = True
        cells['original'].cell_contents = original
        n += 1
    LOG['fr'] = n
    return n


def _src(f):
    import inspect
    try:
        return inspect.getsource(f)
    except (OSError, TypeError):
        return ''


def diag_main(seed):
    """Build the round trip in diagnostic mode and report which candidate words equal K at every hook."""
    import importlib
    import importlib.util
    import types
    import numpy as np
    import dqi_int_convtest as CT
    km = CT.km
    spec = importlib.util.spec_from_file_location('dqi_sw_entry', os.path.join(TOOLS, 'dqi_sw_entry.py'))
    mod = importlib.util.module_from_spec(spec)
    sys.modules['dqi_sw_entry'] = mod
    spec.loader.exec_module(mod)
    for env, pre in (('ZN_LCT_P', 'dqi_zn_'), ('Z2_P', 'dqi_z2_')):
        for name in [p for p in os.environ.get(env, '').split(',') if p]:
            m = importlib.import_module(pre + name)
            print('installed', name, m.install(km, types.SimpleNamespace(decode=mod.decode)), flush=True)
    conv = mod.decode.__globals__['g']['convert_safe']
    saved = CT.hankel_saved(mod.decode)
    b = km.CircuitBuilder()
    A = [b.alloc_qubits(12) for _ in range(35)]
    B = [b.alloc_qubits(12) for _ in range(35)]
    ct = conv(b, A, B, monic_root=True)
    conv(b, A, B, True, phase_token=ct)
    saved.clear()
    c = b.finish_circuit()
    DIAG = sys.modules['dqi_zn_kpark'].DIAG          # the installed copy's list (this file also runs as __main__)
    print('diag qubits', len(DIAG), 'circuit qubits', c.num_qubits, sys.modules['dqi_zn_kpark'].LOG, flush=True)
    rng = np.random.default_rng(seed)
    res = {}
    import dqi_zn_lct as L
    for label, An, Bn in L.chart_list(['0', '1', '2', '3', '4', 'g'], seed)[:12]:
        sim = km.Simulator(batch_size=An.shape[0], debug_print=False, counting_operations=False, seed=seed)
        sim.use_same_registers_as(c)
        sim.write_across_shots(km.array([q for w in A for q in w]), CT.to_bits(An))
        sim.write_across_shots(km.array([q for w in B for q in w]), CT.to_bits(Bn))
        sim.do(c)
        bits = np.asarray(sim.read_across_shots(km.array([d for _, _, d in DIAG])))
        for i, (hk, lab, _) in enumerate(DIAG):
            res.setdefault((hk, lab), []).append(int(bits[i].sum()))
    del rng
    for (hk, lab), v in res.items():
        if sum(v) == 0:
            print('K == [%s == 0] at hook %s on every shot (%d batches)' % (lab, hk, len(v)), flush=True)
    for hk in sorted({h for h, _ in res}):
        best = sorted((sum(v), lab) for (h, lab), v in res.items() if h == hk)[:3]
        print('hook %s: best candidates (shots differing, word): %s' % (hk, best), flush=True)


if __name__ == '__main__':
    if sys.argv[1] == 'diag':
        os.environ['ZN_KP_DIAG'] = '1'
        DIAG_MODE = True
        diag_main(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
