"""R94-CHEM-MARKET: every formula step of the buy and sell side, at TODAY's entry, by the hosts' own cost model.

The hosts' module (dfthcsos.verify, unmodified, imported from ~/qcc/data/verifiers as R91-CHEM-FORMULA did) gives the
heuristic part for any (R, B, C, brot). Today's entry r89_pB1q9e4: lambda_eff 46.907957, loader Toffolis 3,490 a pass
(T_angles 6,980), T 13,480, Q 405 [V, LEAD_NOTES 05:58Z]. A step is priced in PASS-EQUIVALENTS: the loader Toffolis a pass
that change the cost lambda x T x Q by the same factor at today's T and Q:
    pe = ((T + dT) (Q + dQ) / (T Q) - 1) / (2 / T)
(at dQ = 0 that is dT / 2). 1 a pass = 2 / 13,480 = 0.01484 % of the task = 0.00154 of score (1 % of a task = 0.104).
Loaders are held at today's count for the C axis (C moves only the heuristic: R91 note (c)); for the B and R axes the
loader side is printed separately (rows, blocks) because no single law holds (R91 note (d)).

Usage (WSL): ~/qcc/.venv/bin/python -I r94_chem_market_formula.py
"""
import math
import sys
from pathlib import Path

DATA = Path.home() / "qcc" / "data"
sys.path[:0] = [str(DATA / "kickmix"), str(DATA / "verifiers")]
from dfthcsos import verify as V  # noqa: E402  (hosts' module, unmodified)

N = 76
LAM = 46.907957
L_PASS = 3490
T0_ANG = 2 * L_PASS
R0, B0, C0, BROT, J = 9, 127, 46, 7, 5
CAP = V.QUBIT_CAP


def heur(R, B, C, brot=BROT):
    return V.heuristic_costs(N, R, B, C, brot)


def total(R, B, C, angles_T=T0_ANG):
    h = heur(R, B, C)
    rows = N + R * B
    sel = math.ceil(math.log2(rows))
    pers_ang = J * BROT
    tmp_ang = CAP - sel - pers_ang
    T = angles_T + h["toffoli"]
    Q = h["persistent_qubits"] + pers_ang + max(h["temporary_qubits"], tmp_ang)
    return T, Q, h, rows, math.ceil(rows / 16)


T0, Q0, H0, ROWS0, BLK0 = total(R0, B0, C0)
assert (T0, Q0) == (13480, 405), (T0, Q0)
SCORE0 = LAM * T0 * Q0


def pe(T, Q):
    return ((T * Q) / (T0 * Q0) - 1.0) / (2.0 / T0)


def main():
    print(f"TODAY [P hosts' code]: T {T0} (heuristic {H0['toffoli']} + angles {T0_ANG}), Q {Q0}, lambda {LAM}, "
          f"cost {SCORE0:,.2f} (the verified 256,089,299.53); rows {ROWS0}, blocks {BLK0}")
    print("1 a pass = %.5f %% of the task = %.5f of score" % (100 * 2 / T0, 0.104 * 100 * 2 / T0))
    print("\nC axis (loaders fixed at 3,490 a pass): step = this C against the C one lower")
    prev = None
    for C in range(42, 55):
        T, Q, h, _, _ = total(R0, B0, C)
        step = "" if prev is None else f"  step from C {C - 1}: {T - prev:+d} T = {(T - prev) / 2:+.1f} a pass"
        print(f"  C {C}: heuristic {h['toffoli']} (outer {h['outer_state_prep_toffoli']}, inner "
              f"{h['inner_state_prep_toffoli']}, givens {h['givens_rotations_toffoli']}), T {T}, Q {Q}, "
              f"against today {pe(T, Q):+.1f} a pass = {100 * ((T * Q) / (T0 * Q0) - 1):+.3f} %{step}")
        prev = T
    print("\nB axis (C 46; loader side NOT included: rows and blocks printed)")
    for B in range(124, 131):
        T, Q, h, rows, blk = total(R0, B, C0)
        print(f"  B {B}: heuristic {h['toffoli']} (inner {h['inner_state_prep_toffoli']}), persistent "
              f"{h['persistent_qubits']}, temporary {h['temporary_qubits']}, T {T}, Q {Q}: {pe(T, Q):+.1f} a pass "
              f"({100 * ((T * Q) / (T0 * Q0) - 1):+.3f} %), rows {rows} (+{rows - ROWS0}), blocks {blk} (+{blk - BLK0})")
    print("\nR axis (C 46; loader side NOT included)")
    for R in (8, 10):
        T, Q, h, rows, blk = total(R, B0, C0)
        print(f"  R {R}: heuristic {h['toffoli']}, T {T}, Q {Q}: {pe(T, Q):+.1f} a pass "
              f"({100 * ((T * Q) / (T0 * Q0) - 1):+.3f} %), rows {rows} ({rows - ROWS0:+d}), blocks {blk} ({blk - BLK0:+d})")
    print("\nCombined C shapes with the five rungs (C3 05:33Z / 05:46Z counts on e4 [M + P]; lambda held at today's):")
    rungs = {"0 rungs": 3490, "rung 10": 3468, "rungs 10-11": 3446, "rungs 10-12": 3424, "rungs 10-13": 3401,
             "rungs 10-14": 3378, "rungs 10-14, 2 relations each": 3333}
    for C in (46, 47, 48, 49):
        for name, L in rungs.items():
            T, Q, h, _, _ = total(R0, B0, C, 2 * L)
            print(f"  C {C} + {name:30s} L {L}: T {T}, Q {Q}, cost at today's lambda {LAM * T * Q:,.0f} "
                  f"({pe(T, Q):+.1f} a pass)")


if __name__ == "__main__":
    main()
