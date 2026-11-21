#!/usr/bin/env python3
"""dqi_zn_lct.py -- builder Z1: the LINE's converter alone (as the decode calls it), on fixed charts; one
configuration per run, results saved for a word-by-word comparison with another configuration.
Environment: the line's build env (tools/sh/dqi_zn_lenv.sh <cores> sets r12hz's DQI_* and S3's env), entry
tools/dqi_sw_entry.py (B5's + var_ / pub_ modules), B5's harness tools/dqi_int_convtest.py (read-only, imported).
Extra modules of mine: ZN_LCT_P=<name,...> installs tools/dqi_zn_<name>.py install(km, ns) after the entry.
Usage: python dqi_zn_lct.py run <out.npz> <kinds e.g. 0,1,2,3,4,r,g> <seed>
       python dqi_zn_lct.py cmp <a.npz> <b.npz>
"""
import importlib
import importlib.util
import os
import sys
import types

import numpy as np

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
REAL = os.path.expanduser('~/qcc/ws/dqi_zn/box_out/charts_real.npz')
SITE_LEVELS = set()


def chart_list(kinds, seed):
    import dqi_int_convtest as CT
    rng = np.random.default_rng(seed)
    out = []
    for kd in kinds:
        if kd == 'r':
            z = np.load(REAL)
            for i in range(0, z['An'].shape[0], 256):
                out.append(('real %d' % i, z['An'][i:i + 256], z['Bn'][i:i + 256]))
        elif kd == 'g':
            pats = [{m} for m in range(17, 1, -1)] + [{m, m - d} for d in (2, 3, 4) for m in range(17, 2 + d, -1)]
            CT.N = 64
            A = np.concatenate([CT.gap_charts(rng, 64, p)[0] for p in pats])
            B = np.concatenate([CT.gap_charts(rng, 64, p)[1] for p in pats])
            CT.N = 256
            for i in range(0, A.shape[0], 256):
                out.append(('gaps %d' % i, A[i:i + 256], B[i:i + 256]))
        else:
            CT.N = 256
            A, B, label = CT.charts(int(kd), rng)
            out.append((label, A, B))
    return out


def run(path, kinds, seed):
    import dqi_int_convtest as CT
    spec = importlib.util.spec_from_file_location('dqi_sw_entry', os.path.join(TOOLS, 'dqi_sw_entry.py'))
    mod = importlib.util.module_from_spec(spec)
    sys.modules['dqi_sw_entry'] = mod
    spec.loader.exec_module(mod)
    for name in [p for p in os.environ.get('ZN_LCT_P', '').split(',') if p]:
        m = importlib.import_module('dqi_zn_' + name)
        print('installed', name, m.install(CT.km, types.SimpleNamespace(decode=mod.decode)), flush=True)
    conv = mod.decode.__globals__['g']['convert_safe']
    saved = CT.hankel_saved(mod.decode)
    CT.TR['at'] = int(os.environ.get('ZN_LCT_AT', '10000'))
    lvmax = {}
    sites = {}
    global SITE_LEVELS
    SITE_LEVELS = {int(x) for x in os.environ.get('ZN_LCT_SITES', '').split(',') if x}
    if os.environ.get('ZN_LCT_LEVELS') == '1':          # live count high-water mark per recursion level and phase
        inner = CT.CB.alloc_qubits

        def alloc_lv(self, n, _inner=inner):
            r = _inner(self, n)
            f = sys._getframe(1)
            lvl, fn = -1, '?'
            while f is not None:
                m = f.f_locals.get('m')
                if isinstance(m, int) and 'convert' not in f.f_code.co_name:
                    lvl, fn = m, f.f_code.co_name
                    break
                f = f.f_back
            key = (CT.TR.get('mode', '?'), lvl)
            live = len(CT.TR['live'])
            if live > lvmax.get(key, (0, ''))[0]:
                lvmax[key] = (live, fn)
            if os.environ.get('ZN_LCT_SITES') and lvl in SITE_LEVELS:
                ff = sys._getframe(1)
                names = []
                while ff is not None and len(names) < 4:
                    if ff.f_code.co_name not in ('alloc_lv', '_a'):
                        names.append(ff.f_code.co_name)
                    ff = ff.f_back
                skey = (CT.TR.get('mode', '?'), lvl, '<'.join(names))
                if live > sites.get(skey, 0):
                    sites[skey] = live
            return r
        CT.CB.alloc_qubits = alloc_lv
    calls = {}
    if os.environ.get('ZN_LCT_CALLS'):                  # functions actually called, by recursion level (fwd build)
        want = {int(x) for x in os.environ['ZN_LCT_CALLS'].split(',')}

        def prof(frame, event, arg):
            if event != 'call':
                return
            co = frame.f_code
            loc = frame.f_locals
            m = loc.get('m')
            if not isinstance(m, int):
                f = loc.get('f')
                g = loc.get('g')
                if isinstance(f, list) and isinstance(g, list) and len(f) == len(g) + 1:
                    m = len(f) - 1
            if isinstance(m, int) and m in want:
                key = (m, co.co_name, os.path.basename(co.co_filename), co.co_firstlineno,
                       '%x' % (id(frame.f_globals) & 0xffffff))
                calls[key] = calls.get(key, 0) + 1
    circ = {}
    for mode in ('fwd', 'round'):
        CT.TR['mode'] = mode
        if calls is not None and os.environ.get('ZN_LCT_CALLS') and mode == 'round':
            sys.setprofile(prof)
        circ[mode] = CT.build(conv, mode, saved)
        sys.setprofile(None)
        for name in [p for p in os.environ.get('ZN_LCT_P', '').split(',') if p]:
            lg = getattr(sys.modules.get('dqi_zn_' + name), 'LOG', None)
            if lg is not None:
                print('after %s build: %s' % (mode, {k: v for k, v in lg.items() if not isinstance(v, set)}), flush=True)
                evs = getattr(sys.modules.get('dqi_zn_' + name), 'EVENTS', None)
                if evs and mode == 'round':
                    print('events:', evs[-60:], flush=True)
        c, _, _, pk, _ = circ[mode]
        print('%-5s qubits %d, live peak %d, static Toffoli-like %d' % (
            mode, c.num_qubits, pk, c.operation_statistics().num_toffoli_like), flush=True)
    if calls:
        for key in sorted(calls, key=lambda k: (-k[0], k[1])):
            print('call level %2d %-28s %s:%d scope %s x%d' % (key[0], key[1], key[2], key[3], key[4], calls[key]))
        if os.environ.get('ZN_LCT_LEVELS_ONLY') == '1':
            return
    if sites:
        for key in sorted(sites, key=lambda k: (k[0], -k[1], -sites[k]))[:400]:
            print('site %-5s level %2d live %d  %s' % (key[0], key[1], sites[key], key[2]))
    if lvmax:
        for key in sorted(lvmax, key=lambda k: (k[0], -k[1])):
            print('live high-water %-5s level %3d: %d (%s)' % (key[0], key[1], lvmax[key][0], lvmax[key][1]))
        if os.environ.get('ZN_LCT_LEVELS_ONLY') == '1':
            return
    res = dict(C=[], B=[], clean=[], tf=[], tr=[], labels=[])
    for bi, (label, A, B) in enumerate(chart_list(kinds, seed)):
        f = CT.run(circ['fwd'][0], circ['fwd'][1], circ['fwd'][2], A, B, seed * 1000 + bi)
        r = CT.run(circ['round'][0], circ['round'][1], circ['round'][2], A, B, seed * 1000 + bi)
        vals_ok = (r[0] == A).all(axis=1) & (r[1] == B).all(axis=1)
        clean = vals_ok & ~r[2] & ~r[3]
        res.setdefault('vals', []).append(vals_ok)
        res.setdefault('phase', []).append(r[2])
        res.setdefault('dirty', []).append(r[3])
        res['C'].append(f[0])
        res['B'].append(f[1])
        res['clean'].append(clean)
        res['tf'].append(np.full(A.shape[0], f[4]))
        res['tr'].append(np.full(A.shape[0], r[4]))
        res['labels'] += [label] * A.shape[0]
        print('%-22s round trip clean %d/%d (values back %d, phase flips %d, dirty %d); executed fwd %.1f round %.1f'
              % (label, clean.sum(), len(clean), vals_ok.sum(), r[2].sum(), r[3].sum(), f[4], r[4]), flush=True)
    np.savez(path, **{k: (np.concatenate(v) if k != 'labels' else np.array(v)) for k, v in res.items()},
             qubits=np.array([circ['fwd'][0].num_qubits, circ['round'][0].num_qubits]),
             peak=np.array([circ['fwd'][3], circ['round'][3]]))
    print('saved', path, 'mean executed fwd %.1f round %.1f' % (np.concatenate(res['tf']).mean(),
                                                               np.concatenate(res['tr']).mean()), flush=True)


def _align(za, zb):
    """Rows of a and b with the same chart label, in a's order (runs may use other kind orders or subsets)."""
    la, lb = za['labels'], zb['labels']
    ia, ib = [], []
    for L in dict.fromkeys(la):
        pa, pb = np.flatnonzero(la == L), np.flatnonzero(lb == L)
        if len(pa) == len(pb):
            ia.append(pa)
            ib.append(pb)
    ia, ib = np.concatenate(ia), np.concatenate(ib)
    keys = ('C', 'B', 'clean', 'tf', 'tr', 'labels')
    return ({k: za[k][ia] for k in keys} | {'qubits': za['qubits'], 'peak': za['peak']},
            {k: zb[k][ib] for k in keys} | {'qubits': zb['qubits'], 'peak': zb['peak']})


def cmp(a, b):
    za, zb = _align(np.load(a), np.load(b))
    eq = (za['C'] == zb['C']).all(axis=1) & (za['B'] == zb['B']).all(axis=1)
    lab = za['labels']
    print('forward words equal %d/%d; clean a %d b %d; qubits a %s b %s; peak a %s b %s' % (
        eq.sum(), len(eq), za['clean'].sum(), zb['clean'].sum(), za['qubits'], zb['qubits'], za['peak'], zb['peak']))
    for L in sorted(set(lab)):
        s = lab == L
        if not eq[s].all() or not zb['clean'][s].all():
            print('  %-22s equal %d/%d clean b %d/%d' % (L, eq[s].sum(), s.sum(), zb['clean'][s].sum(), s.sum()))
    print('executed fwd a %.1f b %.1f (%+.1f); round a %.1f b %.1f (%+.1f)' % (
        za['tf'].mean(), zb['tf'].mean(), zb['tf'].mean() - za['tf'].mean(), za['tr'].mean(), zb['tr'].mean(),
        zb['tr'].mean() - za['tr'].mean()))


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run(sys.argv[2], sys.argv[3].split(','), int(sys.argv[4]))
    else:
        cmp(sys.argv[2], sys.argv[3])
