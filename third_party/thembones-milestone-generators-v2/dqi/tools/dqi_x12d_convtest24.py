#!/usr/bin/env python3
"""dqi_x12d_convtest24.py -- standalone test of my converter-side patches on builder J's 24-bit flow (workstream
dqi_x12d). Builds the flow's converter alone (convert, then unconvert with the forward's phase token, exactly as
the decode calls them) before and after installing the patches, and runs both on the same charts:
  forward words equal to the unpatched converter's; round trip clean (A and B restored, no phase flip, every
  ancilla back to 0); qubits of the converter alone and executed Toffolis per shot, before and after.
Charts cycle: random monic / forced root gap (B = E^2 + z O^2, D1's generator) / top gaps (B's top words zeroed).
Usage (J's ip4c907 environment, as tools/sh/dqi_x12d_pbbuild.sh sets it):
  python dqi_x12d_convtest24.py <batches> <seed> <patch[,patch...]>      patches = my modules dqi_x12d_<name>
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


def to_bits(vals):
    bits = (vals[:, :, None] >> np.arange(12)) & 1 == 1
    return bits.transpose(1, 2, 0).reshape(-1, vals.shape[0])


def from_bits(bits, nwords):
    b = bits.reshape(nwords, 12, -1).astype(np.int64)
    return (b * (1 << np.arange(12))[None, :, None]).sum(axis=1).T


def charts(kind, rng, n=256):
    if kind == 0:
        return (rng.integers(0, 4096, size=(n, 35), dtype=np.int64),
                rng.integers(0, 4096, size=(n, 35), dtype=np.int64), 'random monic')
    if kind == 1:
        from dqi_mid_conv import random_inputs                            # D1's generator (numpy only)
        A, B = random_inputs(n, rng, root_gap=True)
        return A, B, 'forced root gap'
    A = rng.integers(0, 4096, size=(n, 35), dtype=np.int64)
    B = rng.integers(0, 4096, size=(n, 35), dtype=np.int64)
    top = rng.integers(0, 4, size=n)                                     # 0..3 top words of B zeroed
    for i in range(n):
        if top[i]:
            B[i, 35 - top[i]:] = 0
    return A, B, 'top gaps'


def hankel_saved(decode):
    """The a11 Hankel pair's handoff list (a forward-only build leaves its entry behind; cleared between builds)."""
    walk = sys.modules['walk_bound']
    for fn in walk.functions(decode):
        if fn.__code__.co_filename.endswith('/research/agents/a11/hankel_cache_even.py') and fn.__name__ == 'hankel':
            return dict(zip(fn.__code__.co_freevars, fn.__closure__))['saved'].cell_contents
    raise AssertionError('no a11 hankel')


def build(conv, mode, saved):
    b = km.CircuitBuilder()
    A = [b.alloc_qubits(12) for _ in range(35)]
    B = [b.alloc_qubits(12) for _ in range(35)]
    ct = conv(b, A, B, monic_root=True)
    if mode == 'round':
        conv(b, A, B, True, phase_token=ct)
    saved.clear()
    return b.finish_circuit(), A, B


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
    batches, seed, patches = int(sys.argv[1]), int(sys.argv[2]), [p for p in sys.argv[3].split(',') if p]
    spec = importlib.util.spec_from_file_location('dqi_junc_r8_entry', os.path.join(TOOLS, 'dqi_junc_r8_entry.py'))
    mod = importlib.util.module_from_spec(spec)
    sys.modules['dqi_junc_r8_entry'] = mod
    spec.loader.exec_module(mod)
    conv = mod.decode.__globals__['g']['convert_safe']
    saved = hankel_saved(mod.decode)
    circuits = {}
    for mode in ('fwd', 'round'):
        circuits[('ref', mode)] = build(conv, mode, saved)
    for name in patches:
        log = importlib.import_module('dqi_x12d_' + name).install(None, types.SimpleNamespace(decode=mod.decode))
        print('installed', name, log, flush=True)
    for mode in ('fwd', 'round'):
        circuits[('new', mode)] = build(conv, mode, saved)
    for key, (c, _, _) in circuits.items():
        print('%-3s %-5s qubits %d, static Toffoli-like %d' % (key[0], key[1], c.num_qubits,
              c.operation_statistics().num_toffoli_like), flush=True)
    rng = np.random.default_rng(seed)
    tot = dict(shots=0, fwd_equal=0, ref_clean=0, new_clean=0)
    tof = {}
    for bi in range(batches):
        An, Bn, label = charts(bi % 3, rng)
        out = {}
        for key, (c, A, B) in circuits.items():
            out[key] = run(c, A, B, An, Bn, seed * 1000 + bi)
            tof[key] = tof.get(key, 0) + out[key][4]
        rf, nf = out[('ref', 'fwd')], out[('new', 'fwd')]
        eq = (rf[0] == nf[0]).all(axis=1) & (rf[1] == nf[1]).all(axis=1)
        clean = {}
        for who in ('ref', 'new'):
            r = out[(who, 'round')]
            clean[who] = (r[0] == An).all(axis=1) & (r[1] == Bn).all(axis=1) & ~r[2] & ~r[3]
        tot['shots'] += len(eq)
        tot['fwd_equal'] += int(eq.sum())
        tot['ref_clean'] += int(clean['ref'].sum())
        tot['new_clean'] += int(clean['new'].sum())
        print('batch %d (%s): forward equal %d/%d, round trip clean ref %d new %d' % (
            bi, label, eq.sum(), len(eq), clean['ref'].sum(), clean['new'].sum()), flush=True)
    print('TOTAL', tot)
    print('executed Toffolis per shot:', {'%s %s' % k: round(v / batches, 1) for k, v in tof.items()})


if __name__ == '__main__':
    main()
