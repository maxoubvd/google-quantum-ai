#!/usr/bin/env python3
"""r33_inv-wire2_report.py -- reads the JSON lines of r33_inv-wire2_cm.py and re-prices every set with the
refined row law (own code; reads only its own result files).

Row law [P: tools/r33_inv-wire2_rowlaw.py, exhaustive for m = 3..7; equal to the line's main divider, where a
13-bit row costs 12 and a feed row one more, M F1]:
  exact row that takes a dividend bit in (feed row): m + 1 ANDs = span - 3
  exact row on an appended zero bit (no feed: the lowest carry is a root wire): m ANDs = span - 4
  tapered row k (divisor rounded to its top m + 2 - k bits): m + 2 - k ANDs = span - 2 in the middle reading
      (the line's truncated tail rows cost one less, my rounded form on a register with free low bits one more:
      PESS = 1 prices span - 1)
  undo: the same + 0.5 per carry not held (n - 2 - c, n = adder bits, c clean wires) [P R31 / R32]
  main divider: the line's rows capped at m + 1 bits, + one feed carry each way per root bit below 12
Higher-order stage: see ho_cost (pairs counted [P], the rest [E]).

  r33_inv-wire2_report.py FILE.jsonl [min_margin [PESS [THR]]]     table of the best sets per root split
     PESS -2 / 0 / 1: tapered rows at the line's tail-row price / middle / pessimistic; THR 26 / 24: pairs kept
"""
import json
import math
import sys

BASE = 75436.96
SHARED = dict(NORM=51.0, NORM_INV=43.25, ENDPOINT=114.0, PASSES=-8.0)
W12 = [13] * 15 + [11, 10, 9, 8, 7, 6, 5, 4, 3]
MAIN_FWD12, MAIN_CLEAN12 = 260.0, 309.0


PESS = 0


def div_cost(spans, n_exact, n_nofeed, c):
    f = b = 0.0
    for i, sp in enumerate(spans):
        if i < n_exact - n_nofeed:
            a = max(0, sp - 3)
        elif i < n_exact:
            a = max(0, sp - 4)
        else:
            a = max(0, sp - 2 + PESS)
        n = a + 1
        f += a
        b += a + 0.5 * max(0, n - 2 - c)
    return f, b


THR = 26
THR_SET = None          # a set run in mode "grpt" carries its own threshold


def ho_cost(m, ho, t, q_max, cubic_max):
    """the higher-order stage as the grouped in-place recurrence (design file, section 2 stage 6), both ways.
    [P] digit pairs (e, f), e < f, among the t digits it fixes, whose product x 3 / 2^28 reaches the unit kept:
        e + f >= THR (26: products of 3/4 unit and more; 24: of 3/16 unit and more)
    [E] two adder cells per pair (the factor 3 = two add / subtracts; the sign-controlled form has no product AND),
        the borrow chain inside each group (t (m - 3) / 2) and at its joint (3 per group), one carry chain through
        U (27 - m), the cubic term (its bits squared).
    Sets run with the accuracy proxy ("top", "ideal") are priced with one digit more: the recurrence needs it [M]."""
    if ho == "none":
        return 0.0
    if ho == "ideal":
        t = 26 - 2 * m
    if ho in ("ideal", "top"):
        t = t + 1
    Lb = 27 - m
    lo = max(0, Lb - t)
    thr = THR if THR_SET is None else THR_SET
    pairs = sum(1 for e in range(lo, Lb) for f in range(e + 1, Lb) if e + f >= thr)
    k = max(2, m - 2)
    ripple = t * (m - 3) / 2.0 + 3.0 * math.ceil(t / k)
    cb = 0.0
    if cubic_max > 0.5:
        cb = float(int(math.ceil(math.log2(cubic_max + 1))) ** 2)
    return 2.0 * (2.0 * pairs + ripple + (27 - m) + cb)


def reprice(r):
    m = r["m"]
    prm = r["prm"]
    st = dict(SHARED)
    root = 96.0 * sum(i + 1 for i in range(1, m + 1)) / 90.0
    st["ROOT"] = st["ROOT_INV"] = root
    if m <= 9:
        st["CENTRE"] = 2.0 * (m + 2)
    else:
        st["DELTA_LO"] = 22.0
    red = sum(w - min(w, m + 1) for w in W12)
    st["DIV_FWD"] = MAIN_FWD12 - red + (12 - m) + (2 if int(prm.get("bt", m - 1)) == m - 2 else 0)
    st["DIV_CLEAN"] = MAIN_CLEAN12 - (MAIN_CLEAN12 / MAIN_FWD12) * red + (12 - m)
    q = r["q"]
    ne1 = r["n1"] - (m + 2 + int(prm.get("r1", 0)))
    ne2 = r["n2"] - (m + 2 + int(prm.get("r2", 0)))
    c1 = max(0, q - (r["ledger"]["div1"] - 2))
    c2 = max(0, q - (r["ledger"]["div2"] - 2))
    st["DIV1"], st["DIV1_INV"] = div_cost(r["spans1"], ne1, min(ne1, r["F1"]), c1)
    st["DIV2"], st["DIV2_INV"] = div_cost(r["spans2"], ne2, 0, c2)
    cubic = {5: 1265.0, 6: 158.1, 7: 19.9, 8: 2.5, 9: 0.31}.get(m, 0.0)
    global THR_SET
    THR_SET = int(prm["thr"]) if "thr" in prm else None
    st["HO"] = ho_cost(m, prm.get("ho", "top"), int(prm.get("t", 0)), r.get("q_max", 0.0), cubic)
    tof = sum(st.values())
    marg = r["margin"]
    if marg is None or marg <= 0:
        return st, tof, float("inf"), float("nan")
    tr = 0.56 * math.log2(33 * 2.0 ** 64 / (math.pi * marg)) + 5.3
    cost = q * (tof + (1 + 33 * tr) / 4.0)
    return st, tof, cost, tr


def main():
    global PESS
    path = sys.argv[1]
    mmin = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    PESS = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    global THR
    THR = int(sys.argv[4]) if len(sys.argv) > 4 else 26
    rows = []
    nerr = 0
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        r = json.loads(line)
        if "error" in r:
            nerr += 1
            continue
        st, tof, cost, tr = reprice(r)
        r["_st"], r["_tof"], r["_cost"], r["_tr"] = st, tof, cost, tr
        rows.append(r)
    print(f"{len(rows)} sets, {nerr} errors; feasible (LP margin > 0): {sum(1 for r in rows if r['margin'] > 0)}")
    print("| m | sets | feasible | best cost | gain | Q | Toffolis | LP margin | T/rot | U error | set |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for m in sorted(set(r["m"] for r in rows)):
        rm = [r for r in rows if r["m"] == m]
        fe = [r for r in rm if r["margin"] > mmin]
        if not fe:
            best = max(rm, key=lambda r: r["margin"])
            print(f"| {m} | {len(rm)} | 0 | - | - | {best['q']} | {best['_tof']:.1f} | {best['margin']:.3g} | - | "
                  f"{best['U_err'][0]:.1f}..{best['U_err'][1]:.1f} | `{json.dumps(best['prm'])}` |")
            continue
        best = min(fe, key=lambda r: r["_cost"])
        print(f"| {m} | {len(rm)} | {len(fe)} | {best['_cost']:,.0f} | {100 * (1 - best['_cost'] / BASE):.1f}% | {best['q']} | "
              f"{best['_tof']:.1f} | {best['margin']:.3g} | {best['_tr']:.2f} | {best['U_err'][0]:.1f}..{best['U_err'][1]:.1f} "
              f"| `{json.dumps(best['prm'])}` |")
    print()
    for m in sorted(set(r["m"] for r in rows)):
        fe = sorted([r for r in rows if r["m"] == m and r["margin"] > mmin], key=lambda r: r["_cost"])[:6]
        for r in fe:
            print(f"m {m} cost {r['_cost']:,.0f} Q {r['q']} tof {r['_tof']:.1f} marg {r['margin']:.3g} n1 {r['n1']} n2 {r['n2']} "
                  f"F1 {r['F1']} left {r['left1']} {r['left2']} fn {r['l1_fn_last']} {r['l2_fn_last']} drop {r['dropped']} "
                  f"g {r.get('g')} led {r['ledger']} prm {json.dumps(r['prm'])}")
        if fe:
            b = fe[0]
            print("   stages: " + ", ".join(f"{k} {v:.1f}" for k, v in b["_st"].items()))
            print(f"   spans1 {b['spans1']}\n   spans2 {b['spans2']}")


if __name__ == "__main__":
    main()
