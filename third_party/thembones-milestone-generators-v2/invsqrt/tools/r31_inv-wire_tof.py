"""R31 seat INV-WIRE (2026-11-13): exact expected Toffoli count per op of a .kmb, by structure.

Our own code; the circuit file is read as data (tools/kmb.py through F3's decoder invsqrt_fw_sim.Circ).
A Toffoli-type op (CCX, CCZ, CSWAP) under d distinct classical condition bits (PUSH_CONDITION stack plus its own
condition) runs with probability 2^-d when the bits are independent HMR outcomes. The total must equal the
generator's exact expectation (1511.5703125 for the public file): that is the check of the independence reading.

    python r31_inv-wire_tof.py <file.kmb> <out.npy>      writes the cumulative expectation after every op
    python r31_inv-wire_tof.py <file.kmb> <out.npy> a b [a b ...]   also prints the sum over each op range [a, b)
"""
import sys

import numpy as np

sys.path.insert(0, "<HOME>/Quantum/tools")
import invsqrt_fw_sim as F  # noqa: E402


def cumulative(path):
    c = F.Circ(path)
    stack = []
    cum = np.zeros(len(c.ops) + 1)
    t = 0.0
    for i, (base, cnd, q, b, a, tag) in enumerate(c.ops):
        if base == 0x81:
            stack.append(cnd)
        elif base == 0x01:
            stack.pop()
        elif base in (0x30, 0x31, 0x32):
            bits = set(stack)
            if cnd >= 0:
                bits.add(cnd)
            t += 0.5 ** len(bits)
        cum[i + 1] = t
    return c, cum


if __name__ == "__main__":
    c, cum = cumulative(sys.argv[1])
    np.save(sys.argv[2], cum)
    print("ops", len(c.ops), "expected Toffoli-type", cum[-1])
    r = sys.argv[3:]
    for k in range(0, len(r), 2):
        a, b = int(r[k]), int(r[k + 1])
        print(f"[{a:6d}, {b:6d})  {cum[b] - cum[a]:10.4f}")
