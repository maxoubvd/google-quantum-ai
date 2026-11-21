"""dqi_junc_kit.py -- access to kickmix and to MY COPY of the kit's decoder sources (workstream dqi_junc).

No kit file is edited. The copy lives in ~/qcc/ws/dqi_junc/kit (copied read-only from
~/qcc/ws/dqi_build/src/source/playbook/kit/decoder). Exposes the primitives the junction module uses:
  ar.mul_xor0(b,x,y,z)            z ^= x*y (GF(2^12) or GF(2^6) or GF(2^3)), zero workspace, rank 42/15/6
  mul_xor_noscratch(b,x,y,z)      same, kit's no-scratch module
  norm_xor_noscratch(b,x,z,K)     z ^= Norm(x) for x in tower coordinates, zero workspace
  inverse3_noscratch(b,x)         fresh GF(8) inverse, 3 Toffolis, zero workspace
  ar.tower(n), ar.phase_product, ar._phase_norm, ar.phase_linear, ar.linear_inplace
plus helpers: simulate, to_bits, from_bits, measured.
"""
import os
import sys

sys.dont_write_bytecode = True
DATA = os.path.expanduser('~/qcc/data')
KIT = os.environ.get('DQI_JUNC_KIT', os.path.expanduser('~/qcc/ws/dqi_junc/kit'))
sys.path[:0] = [os.path.join(DATA, 'kickmix')]
import kickmix as km  # noqa: E402

sys.modules['circ'] = km
sys.path[:0] = [KIT, os.path.dirname(KIT)]
from exp.decoder import common as C                                      # noqa: E402
from exp.decoder.common import ar, cx, scale                              # noqa: E402,F401
from exp.noscratch.kernels import mul_xor_noscratch, inverse3_noscratch, rankone_xor  # noqa: E402,F401
from exp.noscratch.schedules import norm_xor_noscratch                    # noqa: E402,F401
import numpy as np                                                        # noqa: E402

gf = km.GF2Field(12)
assert gf.modulus == 0x1053
SQ = tuple(ar.square_columns(12))
SQRT = tuple(gf.pow(1 << j, 2048) for j in range(12))


def square12(b, x, reverse=False):
    """x <- x^2 in place (x <- sqrt(x) with reverse=True): a CNOT network, 0 Toffolis."""
    ar.linear_inplace(b, x, SQRT if reverse else SQ)


def measured(b, reg):
    """X-measure reg (HMR): returns the classical outcome bits; the register is reset to 0."""
    m = b.alloc_bits(len(reg))
    b.hmr(reg, m)
    return m


def to_bits(vals, n=12):
    vals = np.asarray(vals, dtype=np.int64)
    return ((vals[None, :] >> np.arange(n)[:, None]) & 1) == 1


def from_bits(bits):
    bits = np.asarray(bits)
    return (bits.astype(np.int64) * (1 << np.arange(bits.shape[0]))[:, None]).sum(axis=0)


def simulate(circuit, writes, shots, seed=1):
    """writes: list of (register, int values per shot). Returns the simulator after the run."""
    sim = km.Simulator(batch_size=shots, debug_print=False, counting_operations=True, seed=seed)
    for reg, vals in writes:
        sim.write_across_shots(km.array(list(reg)), to_bits(vals, len(reg)))
    sim.do(circuit)
    return sim


def read(sim, reg):
    return from_bits(sim.read_across_shots(km.array(list(reg))))
