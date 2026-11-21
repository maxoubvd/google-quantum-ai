#!/usr/bin/env python3
"""dqi_int_mon1ct.py -- B5 (25th holder, Nov 19): kill test (b) of tools/dqi_int_mon1.py on the converter alone, as the
decode calls it (convert, then unconvert with the forward's token), before and after installing the module.
Method of tools/dqi_int_fdots.py's convtest (B5's harness tools/dqi_int_convtest.py, imported read-only).

Charts (cycled by `kinds`): 0 random; 3 one gap at a random level 2..16; 11 a gap at level 2 (the base level's gap
flag F = 1 on every shot); 10 the last remainder zero (a = 0 at the base level: E, O share a degree-1 factor); 12 both
(a gap at level 2 and a zero last remainder). Checks per batch: forward words equal to the unpatched converter's on
every shot; round trip clean (A, B restored, no phase flip, ancillas 0) for the new circuit wherever the reference is.
Prints the module's LOG (wrapped originals, calls, refusals) and the static Toffoli-like counts.
Usage (line environment): bash tools/sh/dqi_zn_lenv.sh <cores> ~/qcc/.venv/bin/python -B dqi_int_mon1ct.py
  <batches> <seed> [kinds = 0,11,10,12,3]
"""
import importlib.util
import json
import os
import sys
import types

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dqi_int_convtest as CT                               # noqa: E402  (sets sys.modules['circ'] = kickmix)
from dqi_x12c_gf import mul, inv0                           # noqa: E402


def zero_tail_charts(rng, n, gaps):
    """CT.gap_charts with the last remainder set to 0 (a = 0 at the converter's base level)."""
    d = [17]
    while d[-1] > 0:
        d.append(d[-1] - 2 if (d[-1] in gaps and d[-1] >= 2) else d[-1] - 1)
    r = [None] * len(d)
    r[-1] = np.zeros((n, 1), dtype=np.int64)
    r[-2] = CT.rand_poly(rng, n, d[-2])
    for k in range(len(d) - 2, 0, -1):
        q = CT.rand_poly(rng, n, d[k - 1] - d[k])
        t = CT.pmul(q, r[k])
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


def charts(kind, rng):
    if kind == 11:
        A, B = CT.gap_charts(rng, CT.N, {2})
        return A, B, 'gap at 2 (F = 1 at the base)'
    if kind == 10:
        A, B = zero_tail_charts(rng, CT.N, set())
        return A, B, 'zero last remainder (a = 0)'
    if kind == 12:
        A, B = zero_tail_charts(rng, CT.N, {2})
        return A, B, 'gap at 2 and zero last remainder'
    return CT.charts(kind, rng)


def main():
    batches, seed = int(sys.argv[1]), int(sys.argv[2])
    kinds = [int(x) for x in (sys.argv[3] if len(sys.argv) > 3 else '0,11,10,12,3').split(',')]
    spec = importlib.util.spec_from_file_location('dqi_sw_entry', os.path.join(HERE, 'dqi_sw_entry.py'))
    mod = importlib.util.module_from_spec(spec)
    sys.modules['dqi_sw_entry'] = mod
    spec.loader.exec_module(mod)
    import dqi_int_mon1 as M
    km = CT.km
    conv = mod.decode.__globals__['g']['convert_safe']
    saved = CT.hankel_saved(mod.decode)
    CT.TR['at'] = 10 ** 9
    circ = {}
    for mode in ('fwd', 'round'):
        circ[('ref', mode)] = CT.build(conv, mode, saved)
        print('built ref', mode, flush=True)
    print('install:', json.dumps(M.install(km, types.SimpleNamespace(decode=mod.decode))), flush=True)
    conv = mod.decode.__globals__['g']['convert_safe']
    for mode in ('fwd', 'round'):
        circ[('new', mode)] = CT.build(conv, mode, saved)
        print('built new %s; module LOG %s' % (mode, json.dumps(M.LOG)), flush=True)
    stat = {}
    for key, (c, _, _, pk, _) in circ.items():
        stat[key] = c.operation_statistics().num_toffoli_like
        print('%-3s %-5s qubits %d, live peak %d, static Toffoli-like %d' % (key[0], key[1], c.num_qubits, pk,
                                                                             stat[key]), flush=True)
    print('static Toffoli-like, new - ref: fwd %d, round %d' % (stat[('new', 'fwd')] - stat[('ref', 'fwd')],
                                                               stat[('new', 'round')] - stat[('ref', 'round')]),
          flush=True)
    rng = np.random.default_rng(seed)
    tot = dict(shots=0, fwd_equal=0, ref_clean=0, new_clean=0, new_clean_where_ref=0)
    for bi in range(batches):
        An, Bn, label = charts(kinds[bi % len(kinds)], rng)
        out = {}
        for key, (c, A, B, _, _) in circ.items():
            out[key] = CT.run(c, A, B, An, Bn, seed * 1000 + bi)
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
        tot['new_clean_where_ref'] += int((clean['new'] & clean['ref']).sum())
        print('batch %d (%s): forward equal %d/%d, round trip clean ref %d new %d (new: values %d, phase %d, '
              'ancilla %d); executed Toffolis a shot fwd %.1f -> %.1f, round %.1f -> %.1f' % (
                  bi, label, eq.sum(), len(eq), clean['ref'].sum(), clean['new'].sum(),
                  int(((nr[0] != An).any(axis=1) | (nr[1] != Bn).any(axis=1)).sum()), int(nr[2].sum()),
                  int(nr[3].sum()), rf[4], nf[4], out[('ref', 'round')][4], nr[4]), flush=True)
    print('TOTAL', tot, flush=True)
    ok = (tot['fwd_equal'] == tot['shots'] and tot['new_clean'] == tot['ref_clean'] == tot['new_clean_where_ref']
          and circ[('new', 'round')][0].num_qubits == circ[('ref', 'round')][0].num_qubits
          and M.LOG['calls']['fwd'] >= 1 and M.LOG['calls']['rev'] >= 1)
    print('MON1CT %s (forward equal on every shot; new clean exactly where ref clean; qubits equal; both directions '
          'taken)' % ('PASS' if ok else 'FAIL'), flush=True)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
