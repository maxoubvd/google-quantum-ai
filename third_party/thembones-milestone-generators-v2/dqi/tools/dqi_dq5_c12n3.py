#!/usr/bin/env python3
"""dqi_dq5_c12n3.py -- DQ5 (Nov 18), ledger row 6, second piece: the 15 FLAGGED un-synthesis divisions invert the
6-bit subnorm M with the kit's a14 consuming plan (20 products, 12 live with M) instead of R48's _make_n3 (div6's
inverse with the zero flag read from its 3-bit norm: 21 + 2 = 23), with the zero flag read the same way, from the
plan's OWN 3-bit norm: 20 + 2 = 22: -1 a flagged division, -15 a decode [P count].

How: a copy of exp/frontier/a14/decode_inverse20_12.inverse6_consume20_12 (credit: the kit's a14 agent) with one hook
after its prefix has written the GF(8) norm Norm3(M) into raw[:3] (before any raw wire is measured):
  F = zero_flag(raw[:3])  (2 ANDs: Norm3(M) = 0 <=> M = 0);
  M and the raw wires made SAFE, i.e. set to the values they take for the input 1, by CNOTs from F (free): on M = 0
  every raw value is 0, so XOR-ing F times the input-1 constants gives exactly the input-1 state.
The plan then inverts M_safe = M + F (1 -> 1). Interface = R48's _make_n3 produce: (transcript, [zero 3-bit mask, zero
3-bit mask], F); the caller makes N safe from F (R48's PRODUCE3 edit), so the rebuilt M is Norm6(N_safe) = M_safe,
and the kit's `b.cz(M, inverse_mask)` pays the transcript (a phase mask on the plan's input, M_safe); div6's rebuild
pays the two zero masks (no-ops).
Install: token dqi_dq5_c12n3 in DQI_INT_EXT anywhere (it replaces r48_dqi_ret2_lateflag._make_n3, which R48's divider
looks up at the first flagged call). Switch DQI_DQ5_C12N3_OFF=1. Test: `python tools/dqi_dq5_c12n3.py test` (all 64
inputs against _make_n3's contract: M_out = inverse0(M_safe), F = [M = 0], phase 0 after each payment, clean).
"""
import atexit
import json
import os
import sys

STATS = dict(produce=0, makes=0)
_SAFE = {}


def _a14():
    m = sys.modules.get('decode_inverse20_12')
    if m is None:
        import importlib
        for name in ('decode_inverse_retained', 'decode_retained19'):
            f = getattr(sys.modules.get(name), '__file__', None)
            if f:
                a14 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(f))), 'a14')
                if os.path.isfile(os.path.join(a14, 'decode_inverse20_12.py')) and a14 not in sys.path:
                    sys.path.append(a14)
                break
        m = importlib.import_module('decode_inverse20_12')
    return m


def _prefix(b, x, rows, hook):
    """a14's retained_prefix, unchanged except for `hook(b, x, raw)` after the GF(8) norm is in raw[:3]."""
    A = _a14()
    ar, ref, basic, solve_columns, inverse3_inplace = A.ar, A.ref, A.basic, A.solve_columns, A.inverse3_inplace
    to, back, K = ar.tower(6)
    ar.linear_inplace(b, x, to)
    raw = b.alloc_qubits(6)
    forms = basic.raw_forms19()
    for j, (left, right) in enumerate(forms):
        with ref.expose_pair(b, x, left, right) as (a, c):
            b.init_and(a, c, raw[j])
    norm_forms = ref.norm_forms(K)
    norm_rows = [(value >> 7) & 63 for value in norm_forms]
    active = norm_rows[:]
    for row in rows:
        if solve_columns(row, active) is None:
            active.append(row)
    full = active[:]
    for i in range(6):
        if solve_columns(1 << i, full) is None:
            full.append(1 << i)
    columns = [sum(((row >> j) & 1) << i for i, row in enumerate(full)) for j in range(6)]
    ar.linear_inplace(b, raw, columns)
    for j, form in enumerate(norm_forms):
        ref.form_xor(b, x, form & 127, raw[j])
    for q, row in zip(raw[len(active):], full[len(active):]):
        mask = b.alloc_bits(1)[0]
        b.hmr(q, mask)
        for j in ar.bits(row):
            with ref.expose_pair(b, x, *forms[j]) as (a, c):
                b.ccz(a, c, mask)
        b.free(mask)
    hook(b, x, raw, list(raw[len(active):]))         # DQ5: after the extra raw wires are measured (reset to 0)
    transformed = [1, 2, 4]
    zero_positions = list(range(len(active), 6))
    coefficients = []
    for row in rows:
        coefficient = solve_columns(row, active)
        coefficients.append(coefficient)
        physical = coefficient
        if solve_columns(physical, transformed) is not None:
            physical ^= 1 << zero_positions.pop()
        transformed.append(physical)
    columns = [sum(((row >> j) & 1) << i for i, row in enumerate(transformed)) for j in range(6)]
    ar.linear_inplace(b, raw, columns)
    for q, coefficient in zip(raw[3:], coefficients):
        affine = 0
        for j in ar.bits(coefficient & 7):
            affine ^= norm_forms[j] & 127
        ref.form_xor(b, x, affine, q)
    inverse3_inplace(b, raw[:3])
    return raw


def _safe_constants():
    """(x bits, raw bits) at the hook for the input 1 (standard basis), read once from a 1-shot simulation."""
    if _SAFE:
        return _SAFE['x'], _SAFE['raw']
    import numpy as np
    km = sys.modules.get('kickmix') or __import__('kickmix')
    A = _a14()
    b = km.CircuitBuilder()
    x = b.alloc_qubits(6)
    b.x(x[0])                                          # the element 1 in the standard basis
    got = {}

    def hook(bb, xx, raw, zeros):
        got['x'], got['raw'] = list(xx), list(raw)
        raise StopIteration
    try:
        _prefix(b, x, A.PLAN['raw_rows'], hook)
    except StopIteration:
        pass
    circ = b.finish_circuit()
    sim = km.Simulator(batch_size=1, debug_print=False, counting_operations=False, seed=3)
    sim.do(circ)
    st = np.asarray(sim.read_across_shots(km.array([km.q(i) for i in range(circ.num_qubits)]))).astype(int)[:, 0]
    _SAFE['x'] = [int(st[q.id]) for q in got['x']]
    _SAFE['raw'] = [int(st[q.id]) for q in got['raw']]
    return _SAFE['x'], _SAFE['raw']


def make_produce(zero_flag):
    """R48 _make_n3's produce contract on the a14 consuming plan."""
    A = _a14()
    PLAN, ref, retire, change_basis = A.PLAN, A.ref, A.retire, A.change_basis
    sx, sraw = _safe_constants()

    def produce(b, M):
        STATS['produce'] += 1
        box = {}

        def hook(bb, xx, raw, zeros):
            if zeros:                                  # [Norm3(M) = 0] = [M = 0]: 2 ANDs, temp on a reset raw wire
                tmp = zeros[0]
                n3 = raw[:3]
                for q in n3:
                    bb.x(q)
                bb.init_and(n3[0], n3[1], tmp)
                F = bb.alloc_qubits(1)[0]
                bb.init_and(tmp, n3[2], F)
                bb.del_and(n3[0], n3[1], tmp)
                for q in n3:
                    bb.x(q)
                STATS['tmp_on_raw'] = STATS.get('tmp_on_raw', 0) + 1
            else:
                F = zero_flag(bb, raw[:3])
            for q, v in zip(xx, sx):
                if v:
                    bb.cx(F, q)
            for q, v in zip(raw, sraw):
                if v:
                    bb.cx(F, q)
            box['F'] = F
        raw = _prefix(b, M, PLAN['raw_rows'], hook)
        pool = list(b.alloc_qubits(PLAN['max_live'] - 12))
        all_wires = list(M) + list(raw) + pool[:]
        live = list(M) + list(raw)
        transcript = b.alloc_bits(6)
        b.bit_store(transcript, False)
        retire(b, live, pool, transcript, PLAN['begin'])
        for step in PLAN['steps']:
            if step['kind'] == 'new':
                target = pool.pop()
                with ref.expose_pair(b, live, step['left'], step['right']) as (a, c):
                    b.init_and(a, c, target)
                live.append(target)
            else:
                change_basis(b, live, step['pre_rows'])
                with ref.expose_pair(b, live, step['left'], step['right']) as (a, c):
                    b.ccx(a, c, live[step['target']])
            retire(b, live, pool, transcript, step)
        retire(b, live, pool, transcript, PLAN['final'])
        assert len(live) == 6
        for i, target in enumerate(M):
            source = live[i]
            if str(source) == str(target):
                continue
            other = next((j for j, q in enumerate(live) if str(q) == str(target)), None)
            b.swap(source, target)
            live[i] = target
            if other is not None:
                live[other] = source
        kept = set(map(str, M))
        b.free([q for q in all_wires if str(q) not in kept])
        masks = []
        for _ in range(2):
            m = b.alloc_bits(3)
            b.bit_store(m, False)
            masks.append(m)
        return transcript, masks, box['F']
    produce._dq5_c12n3 = True
    return produce


def install(km, ns):
    if os.environ.get('DQI_DQ5_C12N3_OFF', '0') == '1':
        return dict(off=True)
    r48 = sys.modules.get('r48_dqi_ret2_lateflag')
    assert r48 is not None, 'install after r48_dqi_ret2_lateflag'
    old = r48._make_n3
    if getattr(old, '_dq5_c12n3', False):
        return dict(already=True)

    def make_n3(G, div6, zero_flag):
        _, rebuild = old(G, div6, zero_flag)          # R48's rebuild (div6's: pays the two masks by CZ)
        STATS['makes'] += 1
        return make_produce(zero_flag), rebuild
    make_n3._dq5_c12n3 = True
    r48._make_n3 = make_n3
    atexit.register(lambda: sys.stderr.write('dqi_dq5_c12n3: %s\n' % json.dumps(STATS)))
    return dict(wrapped=True)


def _test():
    """All 64 inputs: c12n3's produce against the contract (values, flag, phase after payment, clean) and the same
    harness on R48's _make_n3 produce (div6 + N3 flag) for comparison of Toffolis."""
    import numpy as np
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
    sys.path.insert(0, here)
    import kickmix as km
    sys.modules['circ'] = km
    import dqi_junc_r8_entry as entry                    # noqa: F401  (kit on sys.path)
    import dqi_junc_r8p_div6 as div6
    import r48_dqi_ret2_lateflag as r48
    walk = sys.modules['walk_bound']
    owner = next(f for f in walk.functions(entry.decode)
                 if f.__code__.co_filename.endswith('decode_inverse_retained.py') and 'produce_retained' in f.__globals__)
    G = owner.__globals__
    from exp.arithmetic.division import zero_flag
    A = _a14()
    f6 = A.ar.field(6)
    TO6, FRM6, K6 = A.ar.tower(6)
    forms = {'c12n3': make_produce(zero_flag)}
    p_n3, rebuild_n3 = r48._make_n3(G, div6, zero_flag)
    forms['make_n3'] = p_n3
    for name, produce in forms.items():
        b = km.CircuitBuilder()
        M = b.alloc_qubits(6)
        keep = b.alloc_qubits(6)                         # a copy of the input, made safe below, to pay the phases
        b.cx(M, keep)
        tr, masks, F = produce(b, M)
        b.cx(F, keep[0])                                 # keep = M_safe (the caller's N safe -> rebuilt M = M_safe)
        b.cz(keep, tr)                                   # the kit's cz(M, inverse_mask)
        A.ar.linear_inplace(b, keep, TO6)                # div6's rebuild payment on conj(M_safe)'s halves
        b.cx(keep[3:], keep[:3])
        b.cz(keep[:3], masks[0])
        b.cz(keep[3:], masks[1])
        b.cx(keep[3:], keep[:3])
        A.ar.linear_inplace(b, keep, FRM6)
        circ = b.finish_circuit()
        n = 64
        sim = km.Simulator(batch_size=n, debug_print=False, counting_operations=True, seed=9)
        vals = np.arange(n)
        sim.write_across_shots(km.array(list(M)), ((vals[None, :] >> np.arange(6)[:, None]) & 1).astype(bool))
        sim.do(circ)
        nq = circ.num_qubits
        st = np.asarray(sim.read_across_shots(km.array([km.q(i) for i in range(nq)]))).astype(np.int64)
        out = sum(st[q.id] << k for k, q in enumerate(M))
        fl = st[F.id]
        safe = np.where(vals == 0, 1, vals)
        want = np.array([f6.invert(int(v)) for v in safe])
        flips = np.asarray(sim.read_phase_flipped_across_shots())
        others = [i for i in range(nq) if i not in set(q.id for q in M) | set(q.id for q in keep) | {F.id}]
        dirty = st[others].any(axis=0) if others else np.zeros(n, bool)
        tof = int(sim.operation_statistics().num_toffoli_like) / n
        print('%-8s on 64 inputs: values %d/64, flag %d/64, phase flips %d, dirty %d, Toffoli-like %.2f, qubits %d'
              % (name, int((out == want).sum()), int((fl == (vals == 0)).sum()), int(flips.sum()), int(dirty.sum()),
                 tof, nq), flush=True)


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'test':
    _test()
