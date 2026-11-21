"""T10 seat 6: THE PRICE OF A TABLE AT ITS OWN SHAPE by the hosts' cost code (the lead's 10:27Z order (3): no constants).

The table step's printed price (r48_chem_q7c_vtab.py) uses fixed T 13,600 / Q 405, i.e. shape (9, 127, 46) only. This tool
reads the table's (R, B, C) from its W and prices it with the hosts' module unmodified (dfthcsos.verify.heuristic_costs, as
tools/r94_chem_market_formula.py does): rows = 76 + R B, T = 2 x LOADERS + the heuristic Toffolis, Q = persistent + J brot +
max(temporary, QUBIT_CAP - ceil(log2 rows) - J brot) (J = 5, brot = 7). Checked: (9, 127, 46) gives T 13,480 / Q 405 at
3,490 loader Toffolis a pass (today's entry [P R94]).

usage: ~/qcc/.venv/bin/python -I chem_t10_price.py TABLE.npz LAMBDA_EFF [--loaders 3490]
prints: PRICE R B C rows T Q lambda_eff cost
"""
import math
import os
import sys
from pathlib import Path

import numpy as np

DATA = Path.home() / "qcc" / "data"
sys.path[:0] = [str(DATA / "kickmix"), str(DATA / "verifiers")]
from dfthcsos import verify as V  # noqa: E402  (hosts' module, unmodified)

N, BROT, J = 76, 7, 5


def price(R, B, C, loaders):
    h = V.heuristic_costs(N, R, B, C, BROT)
    rows = N + R * B
    sel = math.ceil(math.log2(rows))
    pers_ang = J * BROT
    T = 2 * loaders + h["toffoli"]
    Q = h["persistent_qubits"] + pers_ang + max(h["temporary_qubits"], V.QUBIT_CAP - sel - pers_ang)
    return rows, T, Q


def main():
    tab, lam = sys.argv[1], float(sys.argv[2])
    loaders = int(sys.argv[sys.argv.index("--loaders") + 1]) if "--loaders" in sys.argv else 3490
    assert price(9, 127, 46, 3490)[1:] == (13480, 405), "the hosts' code no longer gives today's T / Q"
    R, B, C = np.load(tab)["W"].shape
    rows, T, Q = price(R, B, C, loaders)
    if os.environ.get("T10_PRICE", ""):     # t10 seat 6 (the lead's 16:20Z): a format priced by hand, "T0:Q" (half-space: 5476:395)
        t0_, q_ = os.environ["T10_PRICE"].split(":")
        T, Q = int(t0_) + 2 * loaders, int(q_)
    print(f"PRICE {R} {B} {C} rows {rows} T {T} Q {Q} lambda_eff {lam:.6f} cost {lam * T * Q:,.0f} (loaders {loaders} a pass)")


if __name__ == "__main__":
    main()
