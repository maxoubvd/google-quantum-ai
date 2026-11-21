"""chem_ld2_wwidth: word width / words per batch for my loader tools (builder C2, ninth holder).

The builders' compiler library tools/chem_loaders_core.py (read-only) takes the word width and the words per batch
from two module constants (BROT = 8, J = 5); R31's tools/r31_chem-table_loadj.py sets them for its 7-bit prices.
Importing this module after chem_loaders_core sets them from the environment (CHEM_BROT, CHEM_J) in the process
that runs a search, so that my chain (chem_ld2_try / downs / dwalk / bfind / mkbasis) can plan 7-bit tables.
Without the variables nothing changes.
"""
import os

import chem_loaders_core as _core

if os.environ.get("CHEM_BROT"):
    _core.BROT = int(os.environ["CHEM_BROT"])
if os.environ.get("CHEM_J"):
    _core.J = int(os.environ["CHEM_J"])
# fourteenth holder: env C2_QCAP sets the qubit cap of EVERY loader's analysis (engine kill test at 7 bits: 123 =
# the 8-bit engine's 77 workspace wires). Without the variable nothing changes.
if os.environ.get("C2_QCAP"):
    _core.QUBIT_CAP = int(os.environ["C2_QCAP"])

# tenth holder: loader 0's qubit cap for the ANALYSIS (env C2_QCAP0, e.g. 126). k0's 11 low monomials have
# 85 - H free slots under 128 qubits (H = high blocks): 10 at B 139, 9 at B 141, 8 at B 143. chem_loaders_core's
# phase plan then splits them as core 9 + one exclusive coordinate per phase at B 139, and the engine's blocks
# degenerate (a 4-coordinate block with 57 planar relations: k0's protected greedy crashes or crawls). Under a cap
# of 126 the split is B 143's (core 5 + 3 exclusive per phase, the single-phase blocks 9-11 / 12-14). Only k = 0
# is analysed under the cap; every emitted circuit stays under 128 qubits. Without the variable nothing changes.
if os.environ.get("C2_QCAP0") and not getattr(_core, "_c2_qcap0", False):
    _orig_analyze = _core.analyze

    def _analyze_qcap0(words, k, *args, **kwargs):
        if k != 0:
            return _orig_analyze(words, k, *args, **kwargs)
        old = _core.QUBIT_CAP
        _core.QUBIT_CAP = int(os.environ["C2_QCAP0"])
        try:
            return _orig_analyze(words, k, *args, **kwargs)
        finally:
            _core.QUBIT_CAP = old

    _core.analyze = _analyze_qcap0
    _core._c2_qcap0 = True

# fifteenth holder (Nov 16 17:25Z; LEAD_NOTES 17:07Z "keep a carrier-free build ready"): env C2_NOCAR=1 makes EVERY
# loader a pure XOR of its delta, right for every incoming word state: the analysis sees prev and cur XORed with the
# same random pattern (the delta is unchanged), so no bit is "the same pure low function in both directions": no
# carrier is read and no bit-6 wire is taken as clean workspace; all 35 word wires are outputs. k0 (prev = zeros
# inside the engine) is analysed as loader 1 of the table [N | W_0 XOR N] and relabelled 0. Without the variable
# nothing changes.
if os.environ.get("C2_NOCAR") == "1" and not getattr(_core, "_c2_nocar", False):
    import numpy as _np
    _prev_analyze = _core.analyze

    def _analyze_nocar(words, k, *args, **kwargs):
        w = _np.array(words, dtype=_np.int64, copy=True)
        A, ncol = w.shape
        J = _core.J
        rng = _np.random.default_rng(424242 + 7 * k)
        N = rng.integers(0, 1 << _core.BROT, size=(A, J), dtype=_np.int64)
        if k == 0:
            fake = _np.zeros((A, 2 * J), _np.int64)
            fake[:, :J] = N
            fake[:, J:2 * J] = w[:, :J] ^ N
            an = _prev_analyze(fake, 1, *args, **kwargs)
            an.k = 0
            return an
        K = -(-ncol // J)
        pad = _np.zeros((A, K * J), _np.int64)
        pad[:, :ncol] = w
        pad[:, J * (k - 1):J * k] ^= N
        pad[:, J * k:J * k + J] ^= N
        return _prev_analyze(pad[:, :max(ncol, J * k + J)], k, *args, **kwargs)

    _core.analyze = _analyze_nocar
    _core._c2_nocar = True
