#!/usr/bin/env python
"""r32_dqi-ext_ztkern.py -- construction Z-T of R32-DQI-EXT as kickmix kernels, checked and counted.

Z-T: the node z = 0 run as a telescoped group on a 12-bit unit of the carried scale G = (low, high):
  split     f = [low = 0], CSWAP(f, low, high): the low half is a unit nu on every shot with G != 0
  zj_in     first step with the unit: D_1 = nu (u0 + a_0), L_1 = kap_0 nu^2 (v0 + b_0), N_1 = kap_0 nu;
            the two old words are X-measured, their masks deferred to the return (as T1's j_in does); the
            chain's D then moves into the clean wires of L, so the group runs inside chart + G (12 fresh wires
            are held only during zj_in)
  chain     the kit's own `step` (unchanged), here k - 1 levels
  zj_out    the kit's pole toggle and `normal_form` (unchanged); the D register is NOT measured: it is the new
            unit nu' = nu * (nonzero one of N_f, D_f); split undone, f erased by its zero test
            (after the forward half: u in the unit's old wires, v in the wires of N, nu' in the wires of L;
            a decoder relabels them by free swaps)
  zret_a    return start: 1 / nu' rebuilt, the forward's normal-form mask paid by CZ, u and v scaled by nu', nu'^2
            in place, 1 / nu' erased by the kit's inverse phase; pole swap and pole flag undone
  chain     the kit's `step(..., reverse=True)`
  zret_b    return end: 1 / nu, the two words divided in place, the deferred masks of zj_in paid by CZ, 1 / nu erased
            by the kit's inverse phase; split undone
No gauge is measured and no scalar chain is run. Only kit primitives are used (final_decoder.decode: step,
normal_form, scale_word, square12, inverse_initialized, phase_inverse34, guarded.zero_toggle; ar.mul_xor0).

Checks on random inputs (zero low halves and pole cases forced on a share of the shots): values after the forward
half against the classical model (tools/r32_dqi-ext_zt.py conventions), identity of the whole circuit, no phase
flip, every other qubit back to zero; prints executed Toffolis and peak ancillas per part.
Usage (WSL): tr -d '\r' < tools/sh/dqi_kern_jtest.sh > /tmp/j.sh; CORES=0,1 bash /tmp/j.sh r32_dqi-ext_ztkern.py [shots=256] [levels=4]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
sys.dont_write_bytecode = True
import kickmix as km                                                   # noqa: E402

sys.modules['circ'] = km
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dqi_junc_r8_entry as entry                                      # noqa: E402
import dqi_research_model as M                                         # noqa: E402

shots = int(sys.argv[1]) if len(sys.argv) > 1 else 256
K = int(sys.argv[2]) if len(sys.argv) > 2 else 4                       # levels of the chain (even)
assert K % 2 == 0 and K >= 2
walk = sys.modules['walk_bound']
G = next(f.__globals__ for f in walk.functions(entry.decode)
         if f.__name__ == 'external_node_chain' and 'carry_pole' in f.__code__.co_names)
ref = G['ref']
ar, scale, inv = ref.ar, ref.scale, ref.inv
ks, corr = ref.ks, ref.corrections
phase_zero = sys.modules['exp.arithmetic.division'].phase_zero
assert corr[0] == 0 and (K == len(ks) or corr[K] == 0)

real = km.CircuitBuilder()
names = []
for i in range(K):
    names += ['a%d' % i, 'h%d' % i]
names += ['N', 'L', 'C', 'MU']
reg = {n: real.alloc_qubits(12) for n in names}
heads = [(reg['a%d' % i], reg['h%d' % i]) for i in range(K)]
peak = [0]


class Proxy:
    """Delegates to the builder and records the peak number of allocated qubits."""

    def __getattr__(self, name):
        return getattr(real, name)

    def alloc_qubits(self, n):
        r = real.alloc_qubits(n)
        peak[0] = max(peak[0], real.num_allocated_qubits)
        return r


b = Proxy()
base = b.num_allocated_qubits
parts = []


def cut(name):
    parts.append((name, real.finish_circuit(), peak[0] - base))
    peak[0] = real.num_allocated_qubits


def split(low, high):
    f = ref.zero_flag(b, low)
    b.cswap(f, low, high)
    return f


def unsplit(f, low, high):
    b.cswap(f, low, high)
    bit = b.alloc_bits(1)
    b.hmr(f, bit[0])
    b.free(f)
    phase_zero(b, low, bit[0])
    b.free(bit)


N, L, C, MU = reg['N'], reg['L'], reg['C'], reg['MU']
a0, h0 = heads[0]

# ---- forward ----
f = split(C, MU)
cut('split')
b.cx(a0, N)
b.cx(h0, L)
D = b.alloc_qubits(12)
ar.mul_xor0(b, C, N, D)
mN = b.alloc_bits(12)
b.hmr(N, mN)
ref.square12(b, C)
ar.mul_xor0(b, C, L, N)
ref.square12(b, C, True)
mL = b.alloc_bits(12)
b.hmr(L, mL)
scale(b, C, ks[0])
scale(b, N, ks[0])
b.swap(D, L)                               # the wires of L are clean: the chain's D moves there, the 12 fresh
b.free(D)                                  # wires are returned; the group runs inside chart + G from here on
Nc, Dc, Lc = C, L, N                       # chain registers: N in the unit's wires, D in L's, L in N's
cut('zj_in')
for i in range(1, K):
    ref.step(b, heads[i], Nc, Dc, Lc, ks[i], corr[i])
cut('chain')
pole = b.alloc_qubits(1)[0]
ref.guarded.zero_toggle(b, Dc, pole)
b.cswap(pole, Nc, Dc)
m_norm = ref.normal_form(b, Nc, Dc, Lc)
cut('zj_out')
unsplit(f, Dc, MU)
cut('unsplit')
FWD = len(parts)

# ---- return ----
f = split(Dc, MU)
cut('split (return)')
di = ref.inverse_initialized(b, Dc)
b.cz(di, m_norm)
b.free(m_norm)
ref.scale_word(b, Nc, Dc, di)
ref.square12(b, Dc)
ref.square12(b, di)
ref.scale_word(b, Lc, Dc, di)
ref.square12(b, di, True)
ref.square12(b, Dc, True)
m2 = b.alloc_bits(12)
b.hmr(di, m2)
b.free(di)
ref.phase_inverse34(b, Dc, m2)
b.free(m2)
b.cswap(pole, Nc, Dc)
ref.guarded.zero_toggle(b, Dc, pole)
b.free(pole)
cut('zret_a')
for i in reversed(range(1, K)):
    ref.step(b, heads[i], Nc, Dc, Lc, ks[i], corr[i], True)
cut('chain (return)')
scale(b, Nc, inv(ks[0]))
scale(b, Lc, inv(ks[0]))
di = ref.inverse_initialized(b, Nc)
ref.scale_word(b, Dc, di, Nc)
ref.square12(b, Nc)
ref.square12(b, di)
ref.scale_word(b, Lc, di, Nc)
ref.square12(b, di, True)
ref.square12(b, Nc, True)
b.cz(Dc, mN)
b.free(mN)
b.cz(Lc, mL)
b.free(mL)
m3 = b.alloc_bits(12)
b.hmr(di, m3)
b.free(di)
ref.phase_inverse34(b, Nc, m3)
b.free(m3)
b.cx(a0, Dc)
b.cx(h0, Lc)
b.swap(N, L)                               # u0 back in the wires of N, v0 in the wires of L
cut('zret_b')
unsplit(f, Nc, MU)
cut('unsplit (return)')
nq = max(c.num_qubits for _, c, _ in parts)

# ---- inputs ----
rng = np.random.default_rng(7)
val = {n: rng.integers(0, 4096, size=shots, dtype=np.int64) for n in names}
val['C'][::5] = 0                         # zero low half: the flag swaps the halves
val['MU'][val['MU'] == 0] = 1
val['C'][(val['C'] == 0) & (np.arange(shots) % 5 != 0)] = 1


def chain_c(hd, N0, D0, L0):
    for i in range(K):
        a, h = hd[i]
        L0 = M.mul(L0 ^ M.mul(M.mul(a, corr[i]) ^ h, M.mul(D0, D0)), ks[i])
        N0, D0 = M.mul(D0, ks[i]), N0 ^ M.mul(a, D0)
    return N0, D0, L0


# pole cases: choose u0 with D_f = 0 on every seventh shot
for s in range(3, shots, 7):
    hd = [(int(val['a%d' % i][s]), int(val['h%d' % i][s])) for i in range(K)]
    dc = chain_c(hd, 1, 0, 0)[1]
    dd = chain_c(hd, 0, 1, 0)[1]
    if dc:
        val['N'][s] = M.mul(dd, M.inv(dc))


def bits(v, n):
    return ((v[None, :] >> np.arange(n)[:, None]) & 1).astype(bool)


def run(circuit):
    sim = km.Simulator(batch_size=shots, debug_print=False, counting_operations=True, seed=11)
    for n in names:
        sim.write_across_shots(km.array(list(reg[n])), bits(val[n], 12))
    sim.do(circuit)
    allq = km.array([km.q(i) for i in range(nq)])
    state = np.asarray(sim.read_across_shots(allq)).astype(np.int64)
    return state, np.asarray(sim.read_phase_flipped_across_shots()), int(sim.operation_statistics().num_toffoli_like)


def value(state, ids):
    return sum(state[q] << k for k, q in enumerate(ids))


ids = {n: [q.id for q in reg[n]] for n in names}
prev = 0
tof_by = {}
for idx, (name, circuit, anc) in enumerate(parts):
    state, flipped, tof = run(circuit)
    tof_by[name] = (tof - prev) / shots
    print(f'{name:18s} executed Toffolis {(tof - prev) / shots:8.1f}   cumulative {tof / shots:8.1f}   '
          f'peak ancillas beyond the {len(names)} words {anc}')
    prev = tof
    if idx == FWD - 1:
        bad = npole = nswap = 0
        for s in range(shots):
            hd = [(int(val['a%d' % i][s]), int(val['h%d' % i][s])) for i in range(K)]
            low, high = int(val['C'][s]), int(val['MU'][s])
            fz = int(low == 0)
            nu, other = (high, low) if fz else (low, high)
            u0, v0 = int(val['N'][s]), int(val['L'][s])
            Nf, Df, Lf = chain_c(hd, M.mul(nu, u0), nu, M.mul(M.mul(nu, nu), v0))
            pl = int(Df == 0)
            if pl:
                Nf, Df = Df, Nf
            di_ = M.inv(Df)
            want_u, want_v = M.mul(Nf, di_), M.mul(Lf, M.mul(di_, di_))
            low2, high2 = (other, Df) if fz else (Df, other)
            got = (int(value(state[:, s], ids['C'])), int(value(state[:, s], ids['N'])),
                   int(value(state[:, s], ids['L'])), int(value(state[:, s], ids['MU'])))
            if got != (want_u, want_v, low2, high2):
                bad += 1
            npole += pl
            nswap += fz
        print(f'   forward values against the classical model: {shots - bad} of {shots} correct '
              f'(pole cases {npole}, zero low halves {nswap})')
final, flipped, _ = run(parts[-1][1])
ok = np.ones(shots, bool)
for n in names:
    ok &= value(final, ids[n]) == val[n]
used = set(q for n in names for q in ids[n])
others = [q for q in range(nq) if q not in used]
dirty = final[others].any(axis=0)
print(f'RESULT identity {int(ok.sum())} of {shots}, phase flips {int(flipped.sum())}, shots with dirty ancillas '
      f'{int(dirty.sum())}; circuit qubits {nq}')
fw = tof_by['split'] + tof_by['zj_in'] + tof_by['zj_out'] + tof_by['unsplit']
rt = tof_by['split (return)'] + tof_by['zret_a'] + tof_by['zret_b'] + tof_by['unsplit (return)']
print(f'Z-T outside the chain: forward {fw:.1f} (of which the kit\'s own pole toggle and normal form '
      f'{tof_by["zj_out"]:.1f}), return {rt:.1f}')
