"""dqi_int_a1.py -- workstream dqi_int (B5): builder D1's item a1 (tools/dqi_pub_a1.py + tools/dqi_pub_a1model.py,
read-only): the un-synthesis body in CRT form, one rank-one product fewer per body, as a DQI_INT_P patch of
tools/dqi_int_entry.py. Mechanism credit: the synthesis closing of the public DQI decoder (Syed Asad Ali, Apache-2.0),
re-implemented by D1 on our kit's kernels.
Install it after div6f (it copies the sweep's globals, so the division bound there at install time is the one it
uses) and after a3: DQI_INT_P=retj,polepark2,flagpark,retl,div6f,a3,a1.
DQI_INT_A1_SNAP=<dir>: load a frozen copy of D1's two files from <dir> (provenance of a build); default: live files.
DQI_INT_A1_FIX=1 (B5's fix, 08:45Z): r6n's UNFLAGGED division (policy W's NOFLAG folds, latem / div6f keep it) never
calls `before(F)`, and D1's premap sits at the end of `before`; so in those folds the residue term was never added.
With the fix the sweep (D1's, copied) runs the premap right before such a division, where U is still U_pre.
"""
import importlib
import importlib.util
import os
import sys

SWEEP_SRC = '''
def sweep(b, C, B, folds, divide):
    previous.qadic(b, C, B, k=previous.K)
    beta = b.alloc_qubits(12)
    b.x(beta[0])
    leaders = [C[-2], C[-1], beta]
    H = b.alloc_qubits(24)
    h0 = [H[:12], H[12:]]
    b.free(H)
    A = h0 + list(C[:-2])
    heads = previous.prior.rr.tail.heads_from_banks(A, B)
    frames, _, _ = previous.frames_for(A, B)
    for j, (a, d, k) in enumerate(frames):
        n = len(a)
        a1, a2, beta = leaders
        U = previous.combine(a[:2])
        W = previous.combine(d[:2])
        ft = list(a[2:]) + [a1, a2]
        fs = list(d) + [beta]
        st = list(d[2:]) + [beta]
        ss = list(a[2:]) + [a1, a2]
        kfirst, ksecond, _ = previous.rank.bodies(b, a, d, leaders, k)
        if j:
            previous.qframe(b, a, d, leaders, frames[j-1][2], True)
            previous.qframe(b, a, d, leaders, k)
            nf = id(folds[2*j][1]) in _B5_NOFLAG
            crt = nf or _B5_MODE != 'noflagonly'
            def correct_first(F):
                target = a + [a1, a2]
                if crt and _B5_MODE == 'fix':
                    co = _D1A1_CO(k, n + 1, len(target), len(ft))
                else:
                    co = previous.prior.constant_qadic(k, n + 1, len(target))
                previous.controlled_linear_xor(b, F, a2, target, co)
                if crt:
                    _D1A1_PRE(b, a[:2], ft, k)
            def first():
                _D1A1_BODY(b, a[:2], fs, ft, k)
            if nf:
                _B5_LOG['pre_outside'] += 1
                _D1A1_PRE(b, a[:2], ft, k)
            divide_from_norms(b, U, W, folds[2*j][1], correct_first, first if crt else kfirst, divide, K=k)
            _B5_LOG['crt_bodies' if crt else 'kit_bodies'] += 1
        source = previous.combine(a[2:4]) if n > 3 else previous.combine([a[-1], a1])
        nf2 = id(folds[2*j+1][1]) in _B5_NOFLAG
        crt2 = nf2 or _B5_MODE != 'noflagonly'
        def correct_second(F):
            target = d + [beta]
            if crt2 and _B5_MODE == 'fix':
                co = _D1A1_CO(k, n, len(target), len(st))
            else:
                co = previous.prior.constant_qadic(k, n, len(target))
            previous.controlled_linear_xor(b, F, beta, target, co)
            if crt2:
                _D1A1_PRE(b, d[:2], st, k)
        def second():
            _D1A1_BODY(b, d[:2], ss, st, k)
        if nf2:
            _B5_LOG['pre_outside'] += 1
            _D1A1_PRE(b, d[:2], st, k)
        divide_from_norms(b, W, source, folds[2*j+1][1], correct_second, second if crt2 else ksecond, divide, K=k)
        _B5_LOG['crt_bodies' if crt2 else 'kit_bodies'] += 1
    return A, heads, leaders
'''


def _install_fixed(mod, ns):
    """D1's install with B5's sweep (premap outside `before` in NOFLAG folds)."""
    walk = sys.modules['walk_bound']
    co, premap, body = mod.helpers()
    w = sys.modules.get('dqi_build_patch_r5w')
    noflag = w.NOFLAG if w is not None else {}
    log = dict(pre_outside=0, crt_bodies=0, kit_bodies=0)
    # DQI_INT_A1_MODE: fix (D1's flagged-fold constants with premap(co_U)), kitco (flagged folds: the kit's completion
    # constants, premap after the completion), noflagonly (CRT bodies in the NOFLAG folds only; flagged folds keep the
    # kit's quotient_add, exact whatever F is)
    mode = os.environ.get('DQI_INT_A1_MODE', 'fix')
    assert mode in ('fix', 'kitco', 'noflagonly'), mode
    log['mode'] = mode
    made = {}
    done = 0
    for fn in walk.functions(ns.decode):
        g = fn.__globals__
        s = g.get('sweep')
        if fn.__name__ != 'consume_synthesis' or s is None or not hasattr(s, '__code__'):
            continue
        if not s.__code__.co_filename.endswith('norm_deferred_flag_return.py'):
            continue
        if id(s) not in made:
            env = dict(s.__globals__)
            env.update(_D1A1_CO=co, _D1A1_PRE=premap, _D1A1_BODY=body, _B5_NOFLAG=noflag, _B5_LOG=log,
                       _B5_MODE=mode)
            exec(compile(SWEEP_SRC, __file__, 'exec'), env)
            made[id(s)] = env['sweep']
        g['sweep'] = made[id(s)]
        done += 1
    assert done, 'no live consume_synthesis (norm_deferred_flag_return.sweep) found'
    import atexit
    atexit.register(lambda: print('dqi_int_a1 fix:', log, mod.LOG, flush=True))
    return dict(installed=done, sweeps=len(made), noflag=len(noflag), fix=True)


def install(km, ns):
    snap = os.environ.get('DQI_INT_A1_SNAP')
    if snap:
        if snap not in sys.path:
            sys.path.insert(0, snap)                  # dqi_pub_a1 imports dqi_pub_a1model from sys.path
        for name in ('dqi_pub_a1model', 'dqi_pub_a1'):
            spec = importlib.util.spec_from_file_location(name, os.path.join(snap, name + '.py'))
            m = importlib.util.module_from_spec(spec)
            sys.modules[name] = m
            spec.loader.exec_module(m)
        mod = sys.modules['dqi_pub_a1']
        assert os.path.dirname(os.path.abspath(sys.modules['dqi_pub_a1model'].__file__)) == os.path.abspath(snap)
    else:
        mod = importlib.import_module('dqi_pub_a1')
    if os.environ.get('DQI_INT_A1_FIX') == '1':
        return {'pub_a1_b5fix': _install_fixed(mod, ns), 'snap': snap or 'live'}
    return {'pub_a1': mod.install(km, ns), 'snap': snap or 'live'}
