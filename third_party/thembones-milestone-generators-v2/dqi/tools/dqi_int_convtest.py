#!/usr/bin/env python3
"""dqi_int_convtest.py -- workstream dqi_int (B5, round 3): test of converter patches (DQI_INT converter modules) on
the converter alone, as the decode calls it (convert, then unconvert with the forward's token), before and after
installing the patches. Method of tools/dqi_x12d_convtest24.py, plus charts with Euclid gaps forced at chosen levels.

Checks per batch: forward words (C, B) equal to the unpatched converter's; round trip clean (A and B restored, no
phase flip, every ancilla back to 0) for both. Prints the converter's qubits (high-water mark), the number of
allocations that bring the live count to >= the old peak, and executed Toffolis per shot.

Chart kinds (cycled): 0 random monic; 1 forced root gap (D1's generator); 2 top gaps (B's top words zeroed);
3 one gap at a random level 2..16; 4 two or three separated gaps at random levels. A gap at level m: in the Euclid
remainder sequence of (O, E) (B = E^2 + z O^2, O monic of degree 17), the remainder after a polynomial of degree m
has degree m - 2 (built backward from the bottom, every quotient's leader nonzero).
Usage (in the environment of tools/sh/dqi_int_env.sh, DQI_INT_P = the base patches):
  python dqi_int_convtest.py <batches> <seed> <module[,module...]> [kinds=0,1,2,3,4]
modules: tools/dqi_int_<name>.py with install(km, ns) (ns.decode = the decode)
"""
import importlib
import importlib.util
import os
import sys
import types

import numpy as np

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
import kickmix as km                                                      # noqa: E402
sys.modules['circ'] = km
from dqi_x12c_gf import mul, inv0                                         # noqa: E402  (GF(2^12), numpy)

N = 256


def to_bits(vals):
    bits = (vals[:, :, None] >> np.arange(12)) & 1 == 1
    return bits.transpose(1, 2, 0).reshape(-1, vals.shape[0])


def from_bits(bits, nwords):
    b = bits.reshape(nwords, 12, -1).astype(np.int64)
    return (b * (1 << np.arange(12))[None, :, None]).sum(axis=1).T


def pmul(a, b):
    """Polynomial product over GF(2^12), coefficient arrays (n, la) x (n, lb) -> (n, la + lb - 1)."""
    out = np.zeros((a.shape[0], a.shape[1] + b.shape[1] - 1), dtype=np.int64)
    for i in range(a.shape[1]):
        for j in range(b.shape[1]):
            out[:, i + j] ^= mul(a[:, i], b[:, j])
    return out


def rand_poly(rng, n, deg):
    p = rng.integers(0, 4096, size=(n, deg + 1), dtype=np.int64)
    p[:, deg] = rng.integers(1, 4096, size=n, dtype=np.int64)
    return p


def gap_charts(rng, n, gaps):
    """(A, B) with gaps at the levels in `gaps` (same for every shot)."""
    d = [17]
    while d[-1] > 0:
        d.append(d[-1] - 2 if (d[-1] in gaps and d[-1] >= 2) else d[-1] - 1)
    r = [None] * len(d)
    r[-1] = rand_poly(rng, n, 0)
    r[-2] = rand_poly(rng, n, d[-2])
    for k in range(len(d) - 2, 0, -1):
        q = rand_poly(rng, n, d[k - 1] - d[k])
        t = pmul(q, r[k])
        lo = r[k + 1]
        t[:, :lo.shape[1]] ^= lo
        r[k - 1] = t[:, :d[k - 1] + 1]
    O, R0 = r[0], r[1]
    lam = inv0(O[:, 17])
    O = np.stack([mul(lam, O[:, i]) for i in range(18)], axis=1)
    R0 = np.stack([mul(lam, R0[:, i]) for i in range(R0.shape[1])], axis=1)
    c = rng.integers(0, 4096, size=n, dtype=np.int64)
    E = np.stack([mul(c, O[:, i]) for i in range(18)], axis=1)
    E[:, :R0.shape[1]] ^= R0
    Bfull = np.zeros((n, 36), dtype=np.int64)
    Bfull[:, 0::2] = mul(E, E)
    Bfull[:, 1::2] = mul(O, O)
    A = rng.integers(0, 4096, size=(n, 35), dtype=np.int64)
    return A, Bfull[:, :35]


CYC = [-1]


def charts(kind, rng):
    if kind == 0:
        return (rng.integers(0, 4096, size=(N, 35), dtype=np.int64),
                rng.integers(0, 4096, size=(N, 35), dtype=np.int64), 'random monic')
    if kind == 1:
        from dqi_mid_conv import random_inputs
        A, B = random_inputs(N, rng, root_gap=True)
        return A, B, 'forced root gap'
    if kind == 2:
        A = rng.integers(0, 4096, size=(N, 35), dtype=np.int64)
        B = rng.integers(0, 4096, size=(N, 35), dtype=np.int64)
        top = rng.integers(0, 4, size=N)
        for i in range(N):
            if top[i]:
                B[i, 35 - top[i]:] = 0
        return A, B, 'top gaps'
    if kind == 3:
        lv = int(rng.integers(2, 17))
        A, B = gap_charts(rng, N, {lv})
        return A, B, 'gap at %d' % lv
    if kind == 9:                                   # charts from DQI_INT_CT_LOAD (npz with An, Bn), tiled to N
        zf = np.load(os.environ['DQI_INT_CT_LOAD'])
        reps = -(-N // zf['An'].shape[0])
        return np.tile(zf['An'], (reps, 1))[:N], np.tile(zf['Bn'], (reps, 1))[:N], 'loaded charts'
    if kind == 5:                                   # gaps at the levels of DQI_INT_CT_GAPS, one level per batch
        lv = [int(x) for x in os.environ.get('DQI_INT_CT_GAPS', '10,11,12,14,15,16,17').split(',')]
        CYC[0] += 1
        x = lv[CYC[0] % len(lv)]
        A, B = gap_charts(rng, N, {x})
        return A, B, 'gap at %d' % x
    k = int(rng.integers(2, 4))
    lv = set()
    while len(lv) < k:
        x = int(rng.integers(2, 18))
        if all(abs(x - y) >= 3 for y in lv):
            lv.add(x)
    A, B = gap_charts(rng, N, lv)
    return A, B, 'gaps at %s' % sorted(lv)


def hankel_saved(decode):
    walk = sys.modules['walk_bound']
    for fn in walk.functions(decode):
        if fn.__code__.co_filename.endswith('/research/agents/a11/hankel_cache_even.py') and fn.__name__ == 'hankel':
            return dict(zip(fn.__code__.co_freevars, fn.__closure__))['saved'].cell_contents
    raise AssertionError('no a11 hankel')


CB = km.CircuitBuilder
_alloc, _free = CB.alloc_qubits, CB.free
TR = dict(live=set(), at=10 ** 9, hits=0, peak=0)


def _flat(x, out):
    if isinstance(x, (list, tuple)) or type(x).__name__ == 'array':
        for y in x:
            _flat(y, out)
        return
    s = str(x)
    if len(s) > 1 and s[0] == 'q' and s[1:].isdigit():
        out.append(int(s[1:]))


def _where():
    """Call chain of the allocation (kit / patch frames only), innermost first: 'file:line:function < ...'."""
    import traceback
    st = traceback.extract_stack()[:-3]
    keep = [f for f in st if 'convtest' not in f.filename and 'kickmix' not in f.filename]
    return ' < '.join('%s:%d:%s' % (os.path.basename(f.filename), f.lineno, f.name) for f in reversed(keep[-5:]))


def _a(self, n):
    r = _alloc(self, n)
    out = []
    _flat(r, out)
    TR['live'].update(out)
    TR['peak'] = max(TR['peak'], len(TR['live']))
    if len(TR['live']) >= TR['at']:
        TR['hits'] += 1
        if TR.get('where'):
            print('  at %d: %s' % (len(TR['live']), _where()), flush=True)
    return r


def _f(self, regs):
    out = []
    _flat(regs, out)
    TR['live'].difference_update(out)
    return _free(self, regs)


CB.alloc_qubits, CB.free = _a, _f


def build(conv, mode, saved):
    TR.update(live=set(), hits=0, peak=0)
    b = km.CircuitBuilder()
    A = [b.alloc_qubits(12) for _ in range(35)]
    B = [b.alloc_qubits(12) for _ in range(35)]
    ct = conv(b, A, B, monic_root=True)
    if mode == 'round':
        conv(b, A, B, True, phase_token=ct)
    saved.clear()
    return b.finish_circuit(), A, B, TR['peak'], TR['hits']


def run(c, A, B, An, Bn, seed):
    sim = km.Simulator(batch_size=An.shape[0], debug_print=False, counting_operations=True, seed=seed)
    sim.use_same_registers_as(c)
    qa = km.array([q for w in A for q in w])
    qb = km.array([q for w in B for q in w])
    sim.write_across_shots(qa, to_bits(An))
    sim.write_across_shots(qb, to_bits(Bn))
    sim.do(c)
    anc = km.array([km.q(i) for i in range(840, c.num_qubits)])
    return (from_bits(sim.read_across_shots(qa), 35), from_bits(sim.read_across_shots(qb), 35),
            np.asarray(sim.read_phase_flipped_across_shots(), dtype=bool),
            sim.read_across_shots(anc).any(axis=0), sim.operation_statistics().num_toffoli_like / An.shape[0])


def main():
    batches, seed = int(sys.argv[1]), int(sys.argv[2])
    mods = [p for p in sys.argv[3].split(',') if p]
    kinds = [int(x) for x in (sys.argv[4] if len(sys.argv) > 4 else '0,1,2,3,4').split(',')]
    spec = importlib.util.spec_from_file_location('dqi_int_entry', os.path.join(TOOLS, 'dqi_int_entry.py'))
    mod = importlib.util.module_from_spec(spec)
    sys.modules['dqi_int_entry'] = mod
    spec.loader.exec_module(mod)
    conv = mod.decode.__globals__['g']['convert_safe']
    saved = hankel_saved(mod.decode)
    TR['at'] = int(os.environ.get('DQI_INT_CT_AT', '879'))
    circ = {}
    for mode in ('fwd', 'round'):
        circ[('ref', mode)] = build(conv, mode, saved)
    for name in mods:
        log = importlib.import_module('dqi_int_' + name).install(km, types.SimpleNamespace(decode=mod.decode))
        print('installed', name, log, flush=True)
    conv = mod.decode.__globals__['g']['convert_safe']           # a patch may wrap the converter's entry
    for mode in ('fwd', 'round'):
        TR['where'] = mode == 'round' and os.environ.get('DQI_INT_CT_WHERE') == '1'   # call sites of the hits
        circ[('new', mode)] = build(conv, mode, saved)
        TR['where'] = False
    for name in mods:
        print('patch log', name, getattr(sys.modules.get('dqi_int_' + name), 'LOG', None), flush=True)
    for key, (c, _, _, pk, hits) in circ.items():
        print('%-3s %-5s qubits %d, live peak %d, allocations at >= %d: %d, static Toffoli-like %d' % (
            key[0], key[1], c.num_qubits, pk, TR['at'], hits, c.operation_statistics().num_toffoli_like), flush=True)
    rng = np.random.default_rng(seed)
    tot = dict(shots=0, fwd_equal=0, ref_clean=0, new_clean=0)
    tof = {}
    for bi in range(batches):
        An, Bn, label = charts(kinds[bi % len(kinds)], rng)
        out = {}
        for key, (c, A, B, _, _) in circ.items():
            out[key] = run(c, A, B, An, Bn, seed * 1000 + bi)
            tof[key] = tof.get(key, 0) + out[key][4]
        rf, nf = out[('ref', 'fwd')], out[('new', 'fwd')]
        eq = (rf[0] == nf[0]).all(axis=1) & (rf[1] == nf[1]).all(axis=1)
        clean = {}
        for who in ('ref', 'new'):
            r = out[(who, 'round')]
            clean[who] = (r[0] == An).all(axis=1) & (r[1] == Bn).all(axis=1) & ~r[2] & ~r[3]
        nr = out[('new', 'round')]
        tot['shots'] += len(eq)
        tot['fwd_equal'] += int(eq.sum())
        tot['ref_clean'] += int(clean['ref'].sum())
        tot['new_clean'] += int(clean['new'].sum())
        print('batch %d (%s): forward equal %d/%d, round trip clean ref %d new %d (new: values %d, phase %d, '
              'ancilla %d)' % (bi, label, eq.sum(), len(eq), clean['ref'].sum(), clean['new'].sum(),
                               int(((nr[0] != An).any(axis=1) | (nr[1] != Bn).any(axis=1)).sum()),
                               int(nr[2].sum()), int(nr[3].sum())), flush=True)
    print('TOTAL', tot)
    print('executed Toffolis per shot:', {'%s %s' % k: round(v / batches, 1) for k, v in tof.items()})


if __name__ == '__main__':
    main()
