"""R31 seat INV-WIRE (2026-11-13): exact cost model of the fork's width-limited inverse adder, own code.

Model of what the public generator's inverse adder does with `avail` clean wires (read in its source,
paired_arithmetic.inverse_add / _terminal_partition and fusion.arithmetic.partition; nothing is imported):
  * a retained segment of L bits costs L ANDs and keeps its boundary carry on a clean wire to the end; the
    boundary is then X-measured and repaired by a comparator of L - 1 ANDs on outcome 1: L + 0.5 (L - 1);
  * the ordinary final segment of F bits costs F - 1 and needs F - 2 clean wires;
  * with no clean wire left the rest is a dirty terminal (carries stored in source wires): 2 per source bit,
    plus a dirty increment over the zero-source top bits;
  * feed rows (disposable incoming carry): the input carry is released after the first segment (+0.5, +1 wire).
`python r31_inv-wire_invadd.py check <rows.json>` compares the model with the rows measured in the circuit file
(tools/r31_inv-wire_rows.py); `table` prints cost against clean wires; `sched Q` prices the divider clean-up.
"""
import json
import sys


def dirty_mcx_cost(c):
    return 0 if c < 2 else (1 if c == 2 else 4 * c - 8)


def increment_cost(w):
    return sum(dirty_mcx_cost(p + 1) for p in range(w))


def terminal_cost(w, real):
    prefix = min(real, max(0, w - 1))
    return 2 * prefix + (0 if real == w else increment_cost(w - prefix))


def partition(width, clean):
    if width <= 0:
        return []
    if width <= clean + 2:
        return [width]
    for b in range(1, clean + 1):
        final = clean - b + 2
        cap = sum(clean - j for j in range(b)) + final
        if cap < width:
            continue
        rem = width - final
        sizes = []
        for j in range(b):
            L = min(clean - j, rem - (b - j - 1))
            if L <= 0:
                break
            sizes.append(L)
            rem -= L
        if len(sizes) == b and rem == 0:
            return sizes + [final]
    return None


def ordinary_cost(sizes):
    return sum(L + 0.5 * (L - 1) for L in sizes[:-1]) + (sizes[-1] - 1 if sizes[-1] > 1 else 0)


def terminal_partition(width, real_total, clean, release, ordinary_final, dirty=True):
    states = {(0, 0): (0.0, ())}
    cands = []
    for used in range(clean + 1 + int(release)):
        for (start, count), (cost, sizes) in list(states.items()):
            if count != used:
                continue
            capacity = clean - used + int(release and used > 0)
            rel = 0.5 if release and used > 0 else 0.0
            if ordinary_final and width - start <= capacity + 2:
                cands.append((cost + width - start - 1 + rel, sizes + (width - start,), "ord"))
            real = max(0, real_total - start)
            if dirty and real < width - start or dirty and real == width - start:
                cands.append((cost + terminal_cost(width - start, real) + rel, sizes + (width - start,), "dirty"))
            if used == clean + int(release):
                continue
            for L in range(1, min(capacity, width - start - 1) + 1):
                key = (start + L, used + 1)
                val = (cost + L + 0.5 * (L - 1), sizes + (L,))
                if key not in states or val < states[key]:
                    states[key] = val
    return min(cands, key=lambda t: (t[0], t[1])) if cands else None


def inv_cost(n, real, avail, feed):
    """Expected Toffolis of one inverse add: n target bits, `real` source bits, `avail` clean wires."""
    clean = int(min(avail, max(0, n - 2)))
    if n == 1:
        return 0.0, (1,), "ord"
    if feed:
        return terminal_partition(n, real, clean, True, True)
    p = partition(n, clean)
    if p is not None:
        return ordinary_cost(p), tuple(p), "ord"
    return terminal_partition(n, real, clean, False, False)


def table():
    print("feed row n=13 (11 source bits): clean -> cost, sizes")
    for c in range(0, 12):
        print("  ", c, inv_cost(13, 11, c, True))
    print("tail rows (no feed carry): n, clean -> cost")
    for n in (13, 11, 10, 9, 8, 7, 6, 5, 4):
        print("  n", n, [inv_cost(n, n - 2, c, False)[0] for c in range(0, 9)])


if __name__ == "__main__":
    if sys.argv[1] == "table":
        table()
    elif sys.argv[1] == "check":
        rows = json.load(open(sys.argv[2]))
        Q = 46
        tot = 0.0
        for r in rows:
            feed = r["feed"]
            c, sizes, kind = inv_cost(r["n"], r["n"] - 2, Q - r["alloc"], feed)
            tot += abs(c - r["tof"])
            print(f"n {r['n']:2d} alloc {r['alloc']} avail {Q - r['alloc']} feed {int(feed)}: model {c:5.1f} {sizes} {kind}; "
                  f"measured {r['tof']:5.1f}  {'OK' if abs(c - r['tof']) < 1e-9 else 'DIFF'}")
        print("sum |model - measured| =", tot)


def stage_cost(rows, Q, extra_free=0):
    """(sum of model costs of a stage's inverse rows at width Q, rows that do not fit, wires short at the worst
    row). extra_free = wires retired on top of the file's. A row needs alloc <= Q; with 0 clean wires it runs
    as a dirty terminal. The offset measured - model at 46 (rounding carries of the product rows) is kept."""
    tot, short, worst = 0.0, 0, 0
    for r in rows:
        avail = Q - r["alloc"] + extra_free
        if avail < 0:
            short += 1
            worst = max(worst, -avail)
            avail = 0
        off = r["tof"] - inv_cost(r["n"], r.get("real", r["n"] - 2), 46 - r["alloc"], r["feed"])[0]
        tot += inv_cost(r["n"], r.get("real", r["n"] - 2), avail, r["feed"])[0] + off
    return tot, short, worst


if __name__ == "__main__" and sys.argv[1] == "sched":
    # rows files written by r31_inv-wire_rows.py (alloc = wires live at the row's frame in the public file)
    names = sys.argv[2:]
    for nm in names:
        rows = json.load(open(nm))
        base = stage_cost(rows, 46)[0]
        meas = sum(r["tof"] for r in rows)
        print(f"{nm}: measured rows {meas:.1f}; model at 46: {base:.1f}")
        for extra in (0, 1, 2, 3, 4):
            out = []
            for Q in (48, 47, 46, 45, 44, 43, 42, 41, 40):
                c, short, worst = stage_cost(rows, Q, extra)
                out.append(f"Q{Q}: {c - base:+6.1f}" + (f"[{short}r,-{worst}]" if short else ""))
            print(f"   +{extra} retired: " + "  ".join(out))
