#!/usr/bin/env python3
"""dqi_dq5_c12.py -- DQ5 (Nov 18), ledger row 6 (the un-synthesis division's inverse machinery): the 17 UNFLAGGED
un-synthesis divisions invert the 6-bit subnorm M in place with the kit's own a14 consuming plan
(exp/frontier/a14/decode_inverse20_12.inverse6_consume20_12: zero-correct, 20 products, 12 live wires with M)
instead of J's div6 inverse (tools/dqi_junc_r8p_div6.py, read only: 21 products, the same 12 live wires).

How the line reaches the unflagged divisions: R48's late-flag wrapper (tools/r48_dqi_ret2_lateflag.py, R48_DIV6ALL=1)
builds, per divider scope, the cells div6._make(G) and swaps them into the divider's produce_retained /
rebuild_retained for each unflagged call; the flagged calls (R48_N3FLAG=1) use _make_n3, which takes ONLY the rebuild
from div6._make. This module wraps div6._make and replaces ONLY its produce: the consuming plan, returning div6's own
interface (transcript, [zero 3-bit mask, zero 3-bit mask]); div6's rebuild is returned unchanged (it pays the two
masks by CZ: here no-ops). The transcript, a phase mask on the ORIGINAL M, is paid by the kit's own
`b.cz(M, inverse_mask)` right after the rebuild (div6 returns an all-zero mask there). So the flagged calls build
exactly as before, any rebuild that takes div6's interface (div6's, or B5's div6k at 14 ANDs) composes with it, and the
unflagged calls lose one Toffoli each: -1 x 17 = -17 a decode [P count; M bench 338.94 against 339.91].
Width: identical (12 live with M at the inverse; nothing kept through the parts).
Install: token dqi_dq5_c12 in DQI_INT_EXT anywhere (R48 builds its cells lazily, at the first division of the build).
Switch: DQI_DQ5_C12_OFF=1. Record: DQI_DQ5_C12_OUT=<json> (counts).
Tests: `python tools/dqi_dq5_c12.py test6` (the GF(64) kernel alone on all 64 inputs: value, phase, clean),
the division bench through tools/dqi_dq5_divprof.py (forms lateu+div6, lateu+c12).
"""
import atexit
import json
import os
import sys

STATS = dict(produce=0, makes=0)


def _consume_fn():
    """The kit's a14 executor (exp/frontier/a14/decode_inverse20_12.py of the build's own kit tree)."""
    m = sys.modules.get('decode_inverse20_12')
    if m is None:
        import importlib
        for name in ('decode_inverse_retained', 'decode_retained19', 'decode_inverse20_12'):
            f = getattr(sys.modules.get(name), '__file__', None)
            if f:
                a14 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(f))), 'a14')
                if os.path.isfile(os.path.join(a14, 'decode_inverse20_12.py')) and a14 not in sys.path:
                    sys.path.append(a14)
                break
        m = importlib.import_module('decode_inverse20_12')
    return m.inverse6_consume20_12


def _wrap_make(div6):
    old_make = div6._make
    if getattr(old_make, '_dq5_c12', False):
        return

    def make(G):
        p0, r0 = old_make(G)
        consume = _consume_fn()
        STATS['makes'] += 1

        def produce(b, M):
            STATS['produce'] += 1
            transcript = consume(b, M)               # M <- inverse0(M) in place; phase debt on the original M
            masks = []
            for _ in range(2):                         # div6's interface: two 3-bit masks for its rebuild (zero here)
                m = b.alloc_bits(3)
                b.bit_store(m, False)
                masks.append(m)
            return transcript, masks
        produce._dq5_c12 = True
        return produce, r0
    make._dq5_c12 = True
    div6._make = make


def _dump():
    path = os.environ.get('DQI_DQ5_C12_OUT')
    if path:
        try:
            json.dump(STATS, open(path, 'w'))
        except Exception:
            pass
    sys.stderr.write('dqi_dq5_c12: %s\n' % json.dumps(STATS))


def install(km, ns):
    if os.environ.get('DQI_DQ5_C12_OFF', '0') == '1':
        return dict(off=True)
    div6 = sys.modules.get('dqi_junc_r8p_div6')
    if div6 is None:
        import importlib
        div6 = importlib.import_module('dqi_junc_r8p_div6')
    _consume_fn()                                  # fail at install, not mid-build, if the kit file is missing
    _wrap_make(div6)
    atexit.register(_dump)
    return dict(wrapped=True)


def _test6():
    """inverse6_consume20_12 alone on all 64 inputs: x -> inverse0(x), phase paid by cz(x_orig, transcript), clean."""
    import numpy as np
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
    sys.path.insert(0, here)
    import kickmix as km
    sys.modules['circ'] = km
    import dqi_junc_r8_entry as entry                    # noqa: F401  (puts the kit on sys.path, as in the bench)
    consume = _consume_fn()
    ar = sys.modules[consume.__module__].ar
    f = ar.field(6)
    b = km.CircuitBuilder()
    x = b.alloc_qubits(6)
    keep = b.alloc_qubits(6)                             # a copy of the input, to pay the transcript's phase
    b.cx(x, keep)
    tr = consume(b, x)
    b.cz(keep, tr)
    b.free(tr)
    nalloc = b.num_allocated_qubits
    circ = b.finish_circuit()
    st = circ.operation_statistics()
    n = 64
    sim = km.Simulator(batch_size=n, debug_print=False, counting_operations=True, seed=5)
    vals = np.arange(n)
    bitsv = ((vals[None, :] >> np.arange(6)[:, None]) & 1).astype(bool)
    sim.write_across_shots(km.array(list(x)), bitsv)
    sim.do(circ)
    nq = circ.num_qubits
    state = np.asarray(sim.read_across_shots(km.array([km.q(i) for i in range(nq)]))).astype(np.int64)
    out = sum(state[q.id] << k for k, q in enumerate(x))
    kp = sum(state[q.id] << k for k, q in enumerate(keep))
    flips = np.asarray(sim.read_phase_flipped_across_shots())
    want = np.array([0 if v == 0 else f.invert(int(v)) for v in vals])
    others = [i for i in range(nq) if i not in set(q.id for q in x) | set(q.id for q in keep)]
    dirty = state[others].any(axis=0) if others else np.zeros(n, bool)
    tof = int(sim.operation_statistics().num_toffoli_like) / n
    print('consume20_12 on 64 inputs: values %d/64 (basis as given), keep %d/64, phase flips %d, dirty %d, '
          'executed Toffoli-like %.2f, live at end %d' % (int((out == want).sum()), int((kp == vals).sum()),
                                                         int(flips.sum()), int(dirty.sum()), tof, nalloc))


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'test6':
    _test6()
