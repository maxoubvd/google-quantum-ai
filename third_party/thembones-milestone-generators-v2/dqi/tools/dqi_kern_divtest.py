#!/usr/bin/env python
"""dqi_kern_divtest.py -- the un-synthesis division x -> x / y (24-bit words, y kept) as a kickmix kernel, checked
and counted (workstream dqi_kern, builder B6, item D3).

Variants (argument 1):
  ip9u, ip9f  the decoder's division as in ip9: kit divide_from_norms + r6n + latem + div6 (unflagged / flagged)
  x20u, x20f  the same source with the kit's out-of-place 6-bit inverse (inverse6_xor20: M kept, part masks paid
              at once, no rebuild)
  any name in dqi_kern_div.VARIANTS: my divisions, called as fn(b, x, y, K, g) (g = the fold divider's globals)
Usage (through tools/sh/dqi_kern_jtest.sh, DQI_JUNC_R8P=latem,div6): dqi_kern_divtest.py <variant> [shots] [seed]
Checks on random x, y != 0 (flagged variants also get y = 0 with x = 0): x_new * y == x_old in the kit's own
24-bit field (the 99-term tensor of exp/mul99, evaluated classically), y unchanged, no phase flip, every other
qubit back to 0. Prints executed Toffolis per division and the peak of qubits beyond x and y.
"""
import os
import sys
import types

import numpy as np

sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
sys.dont_write_bytecode = True
import kickmix as km                                                   # noqa: E402

sys.modules['circ'] = km
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dqi_junc_r8_entry as entry                                      # noqa: E402

variant = sys.argv[1] if len(sys.argv) > 1 else 'ip9u'
shots = int(sys.argv[2]) if len(sys.argv) > 2 else 512
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1

walk = sys.modules['walk_bound']
owner = next(f for f in walk.functions(entry.decode)
             if f.__code__.co_filename.endswith('decode_inverse_retained.py') and 'produce_retained' in f.__globals__)
G = owner.__globals__
ar = G['ar']
K = ar.extension_columns()[2]
mod = sys.modules['norm_deferred_flag_return']
r6n = sys.modules['dqi_build_patch_r6n']
latem = sys.modules.get('dqi_junc_r8p_latem')
mul99 = next(m for m in list(sys.modules.values()) if getattr(m, 'tensor99_extension', None) is not None)
TENSOR = mul99.tensor99_extension(K)


def kit_variant(flagged, streamed, flagn=False):
    """The decoder's division source (r6n edits + latem's late M); `streamed` False drops the in-place 6-bit
    inverse (produce_retained) so the kit's out-of-place inverse6_xor20 path runs."""
    text = open(mod.__file__).read()
    start = text.index('def divide_from_norms(')
    src = text[start:text.index('\ndef sweep(', start)].rstrip() + '\n'
    edits = (r6n.FLAGGED if flagged else r6n.UNFLAGGED) + (latem.LATE_M if latem else ())
    if flagn:                                    # B6's flag from the 12-bit norm after the product
        import dqi_kern_flagn as fl
        edits = fl.edits(r6n, latem)
    fn = r6n._variant(mod, src, edits, 'kern_div_' + ('f' if flagged else 'u'))
    own = owner
    if not streamed:
        g2 = {k: v for k, v in G.items() if k not in ('produce_retained', 'norm_keep2')}
        if 'inverse6_xor20' not in g2:
            g2['inverse6_xor20'] = __import__('exp.inverse_search.inverse20',
                                              fromlist=['inverse6_xor20']).inverse6_xor20
        own = types.FunctionType(owner.__code__, g2, owner.__name__)

    def run(b, x, y):
        fm = b.alloc_bits(1)
        b.bit_store(fm, False)
        import dqi_kern_div as D
        fn(b, x, y, fm, D.test_before(b, x), lambda: None, own, K)
    return run


flagged = (variant.endswith('f') and variant[:-1] in ('ip9', 'x20')) or variant in ('ip9g', 'x20g')
if variant in ('ip9g', 'x20g'):
    body = kit_variant(True, variant == 'ip9g', flagn=True)
elif variant in ('ip9u', 'ip9f'):
    body = kit_variant(variant == 'ip9f', True)
elif variant in ('x20u', 'x20f'):
    body = kit_variant(variant == 'x20f', False)
else:
    import dqi_kern_div as D
    fn = D.VARIANTS[variant]
    flagged = getattr(fn, 'flagged', False)

    def body(b, x, y):
        fn(b, x, y, K, G)

real = km.CircuitBuilder()
x = real.alloc_qubits(24)
y = real.alloc_qubits(24)
base = real.num_allocated_qubits
peak = [base]


class Proxy:
    def __getattr__(self, n):
        return getattr(real, n)

    def alloc_qubits(self, n):
        r = real.alloc_qubits(n)
        peak[0] = max(peak[0], real.num_allocated_qubits)
        return r


body(Proxy(), x, y)
circ = real.finish_circuit()
nq = circ.num_qubits

rng = np.random.default_rng(seed)
xv = rng.integers(0, 1 << 24, size=shots, dtype=np.int64)
yv = rng.integers(1, 1 << 24, size=shots, dtype=np.int64)
xv[3::17] = 0
yv[5::19] = 1
yv[7::23] = 1 << 12
YZX = os.environ.get('DQI_KERN_YZERO_X', '0') == '1'      # y = 0 with x kept (the decoder's terminated folds)
if flagged:
    yv[::13] = 0
    if not YZX:
        xv[::13] = 0


def bits(v, n):
    return ((v[None, :] >> np.arange(n)[:, None]) & 1).astype(bool)


def mul(a, c):
    """The kit's 24-bit product (99-term tensor), vectorised over shots."""
    z = np.zeros_like(a)
    for mask, v in TENSOR:
        pa = np.zeros_like(a)
        pc = np.zeros_like(a)
        for i in range(24):
            if (mask >> i) & 1:
                pa ^= (a >> i) & 1
                pc ^= (c >> i) & 1
        z ^= np.where((pa & pc) == 1, v, 0)
    return z


sim = km.Simulator(batch_size=shots, debug_print=False, counting_operations=True, seed=seed + 100)
sim.write_across_shots(km.array(list(x)), bits(xv, 24))
sim.write_across_shots(km.array(list(y)), bits(yv, 24))
sim.do(circ)
state = np.asarray(sim.read_across_shots(km.array([km.q(i) for i in range(nq)]))).astype(np.int64)
flips = np.asarray(sim.read_phase_flipped_across_shots())
tof = int(sim.operation_statistics().num_toffoli_like) / shots


def value(ids):
    return sum(state[q] << k for k, q in enumerate(ids))


xn = value([q.id for q in x])
yn = value([q.id for q in y])
okv = (mul(xn, yv) == xv) & (yn == yv)
if flagged:                                  # y = 0: x must end as x_old ^ the caller's correction (DQI_KERN_BEFORE)
    want0 = (xv if YZX else 0) ^ int(os.environ.get('DQI_KERN_BEFORE', '0'))
    okv = np.where(yv == 0, (xn == want0) & (yn == yv), okv)
others = [q for q in range(nq) if q not in set(q.id for q in x) | set(q.id for q in y)]
dirty = state[others].any(axis=0) if others else np.zeros(shots, bool)
print(f'{variant}: executed Toffolis per division {tof:.2f}; peak qubits beyond x, y {peak[0] - base}; '
      f'values {int(okv.sum())}/{shots}, phase flips {int(flips.sum())}, dirty {int(dirty.sum())}'
      f'{", y=0 cases " + str(int((yv == 0).sum())) if flagged else ""}', flush=True)
