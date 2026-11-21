"""Micro-tests of the OFFICIAL invsqrt verifier's counting (audit workstream).

Spec: width=1, max_x=1 (only x=1 is tested, target phase 1/sqrt(1) = 1 turn = 0, tol = 2^-24 turns),
256 shots. Every circuit below has zero net phase on x=1, so the official verify() passes and reports
its own counts. This confirms how the verifier counts, without relying on documentation.

Usage (WSL): python audit_micro_invsqrt.py
"""
import math
import os
import sys
from fractions import Fraction

DATA = os.path.expanduser("~/qcc/data")
sys.path[:0] = [os.path.join(DATA, "kickmix"), os.path.join(DATA, "verifiers")]
import kickmix as km
from invsqrt import verify_submission_invsqrt as v

HEAD = 'APPEND_TO_REGISTER q0 r0\nREGISTER r0 "x"\n'
spec = v.ProblemSpecInvSqrt(width=1, max_x=1, beta=Fraction(1), precision_bits=24, shots_per_x=256,
                            system_qubits=0)

CASES = {
    "A empty": "",
    "B CCX on idle zero qubits": "CCX q1 q2 q3\n",
    "C CCX if b (b stored 0)": "BIT_STORE0 b0\nCCX q1 q2 q3 if b0\n",
    "C2 CCX if b (b stored 1)": "BIT_STORE1 b0\nCCX q1 q2 q3 if b0\n",
    "D CCX if b (b = HMR outcome)": "HMR q1 b0\nCCX q1 q2 q3 if b0\n",
    "E CCX+CCZ+CSWAP inside PUSH on HMR bit": "HMR q1 b0\nPUSH_CONDITION if b0\nCCX q1 q2 q3\nCCZ q1 q2 q3\nCSWAP q1 q2 q3\nPOP_CONDITION\n",
    "E2 nested PUSH on two HMR bits": "HMR q1 b0\nHMR q2 b1\nPUSH_CONDITION if b0\nPUSH_CONDITION if b1\nCCX q1 q2 q3\nPOP_CONDITION\nPOP_CONDITION\n",
    "E3 PUSH on b0 then CCX if b1": "HMR q1 b0\nHMR q2 b1\nPUSH_CONDITION if b0\nCCX q1 q2 q3 if b1\nPOP_CONDITION\n",
    "F T and T-dagger on x": "Z_POW q0 0.25\nZ_POW q0 1.75\n",
    "F2 S and Z on x (twice each, net 0)": "Z_POW q0 0.5\nZ_POW q0 1.5\nZ_POW q0 1.0\nZ_POW q0 1.0\n",
    "G two level-4 rotations (+1/16, -1/16 turn)": "Z_POW q0 0.125\nZ_POW q0 1.875\n",
    "H same but both if HMR bit": "HMR q1 b0\nZ_POW q0 0.125 if b0\nZ_POW q0 1.875 if b0\n",
    "I level-4 rotation on idle zero qubit": "Z_POW q1 0.125\n",
    "I2 two tiny-angle rotations level 60": f"Z_POW q0 {2**-59:.70f}\nZ_POW q0 {2 - 2**-59:.70f}\n",
    "J gap in ids (X q100 twice)": "X q100\nX q100\n",
    "K DEBUG_PRINT q200": "DEBUG_PRINT q200\n",
    "L 3/4 half-turn twice + S (net 2 turns)": "Z_POW q0 0.75\nZ_POW q0 0.75\nZ_POW q0 0.5\n",
    "M SWAP and CX do not count": "SWAP q1 q2\nCX q1 q2\nCZ q1 q2\n",
    "N CCX if bit that was set by BIT_INVERT if HMR": "HMR q1 b0\nBIT_INVERT b1 if b0\nCCX q1 q2 q3 if b1\n",
    "O classical AND of two coins gates a CCX": "HMR q1 b0\nHMR q2 b1\nPUSH_CONDITION if b0\nBIT_INVERT b2 if b1\nPOP_CONDITION\nCCX q1 q2 q3 if b2\n",
}

for name, body in CASES.items():
    c = km.Circuit(HEAD + body)
    try:
        r = v.verify(c, spec, max_workers=1)
        eps = r["rotation_eps"]
        per_rot = v.rotation_t_count(eps) if eps else 0.0
        print(f"{name:48s} Q {r['num_qubits']:4d}  tof {r['toffoli_count']:.4f}  static {r['max_toffoli_count']}  "
              f"t {r['t_count']:.4f}  rot {r['rotation_count']:.4f}  eps {eps}  T/rot {per_rot:.3f}  "
              f"score {r['score']:.4f}")
    except Exception as ex:  # noqa: BLE001
        print(f"{name:48s} FAIL {type(ex).__name__}: {ex}")

# The eps formula check: eps = pi * margin / n_rot with margin = tol(1) = 2^-24 turns when the error is 0.
print("expected eps for n_rot=2:", math.pi * 2**-24 / 2, " n_rot=1:", math.pi * 2**-24)
