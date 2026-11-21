"""invsqrt_ph: candidate reflection sets for the backward swap network (UNNORM on the half thermometer); own
code, workstream invsqrt_ph. Reflection idea credited in invsqrt_ph_netsearch.py.

Backward direction: for exponent e (s = 21 - e) the data bits sit on wires s..20, the leading one on 21; they
must end on 0..e-1 and e. A layer is a reflection q <-> c - q; an even number of reflections c1, c2, ... moves
a bit by (c2 - c1) + (c4 - c3). Controls are of the family alpha(e >> 1) xor beta (e & 1): the chains of the
exponents 2m and 2m + 1 must differ by one fixed set B of layers.

This script lists ordered sets of L sums for which every exponent >= emin has a chain, with the family
condition, and prices each. Every bit is followed on its own: a bit whose partner under a reflection would be
outside the wires 0..21 stays where it is; a chain is accepted when at most `relax` data bits end on a wrong
wire (those and a leading one that does not ride need patch swaps). Price = swaps needed in the reflections
(a swap is needed if some exponent has a bit on one of its wires while the layer is on) + misplaced bits.
The exact search (invsqrt_ph_netsat.py, mode bws) then finishes the small exponents and the patches for the
best sets; this script only ranks sum sets for it.

usage: invsqrt_ph_bwdsets.py [L=5] [emin=8] [top=25] [lo=7] [hi=34] [relax=1] [cs=19,20,21,22,23]
       invsqrt_ph_bwdsets.py sums=26,22,20,19,11 dump=<json> [emin=6] [relax=1]   chains of one set, for the
                             option chains= of invsqrt_ph_netsat.py
"""
import itertools
import json
import sys

EVEN = {}


def even_masks(n):
    if n not in EVEN:
        EVEN[n] = [m for m in range(1 << n) if not bin(m).count("1") & 1]
    return EVEN[n]


def run_chain(sums, mask, e):
    """Follow every bit. Returns (misplaced data bits, leading one misplaced, set of (layer, low wire) used)."""
    s = 21 - e
    pos = list(range(s, 21)) + [21]                # data bit i on s + i, leading one last
    used = set()
    for i, c in enumerate(sums):
        if (mask >> i) & 1:
            for k, q in enumerate(pos):
                p = c - q
                if 0 <= p <= 21 and p != q:
                    pos[k] = p
                    used.add((i, min(p, q)))
    bad = sum(1 for k in range(e) if pos[k] != k)
    return bad, int(pos[e] != e), used


def options(sums, e, relax):
    out = []
    for mask in even_masks(len(sums)):
        bad, lead, used = run_chain(sums, mask, e)
        if bad <= relax:
            out.append((mask, bad + lead, used))
    return out


def price(sums, emin, relax):
    n = len(sums)
    opts = {}
    for e in range(21, emin - 1, -1):
        o = options(sums, e, relax)
        if not o:
            return None
        opts[e] = o
    best = None
    for B in range(1, 1 << n):
        per_m = []
        for m in range((emin + 1) >> 1, 11):
            lst = [(oa, ob) for oa in opts[2 * m] for ob in opts[2 * m + 1] if (oa[0] ^ ob[0]) == B]
            if not lst:
                per_m = None
                break
            lst.sort(key=lambda t: t[0][1] + t[1][1])
            per_m.append(lst[:3])
        if per_m is None:
            continue
        for combo in itertools.product(*per_m):
            used, bad = set(), 0
            for oa, ob in combo:
                used |= oa[2] | ob[2]
                bad += oa[1] + ob[1]
            cost = len(used) + bad
            if best is None or cost < best[0]:
                best = (cost, len(used), bad, B, [o[0] for pr in combo for o in pr])
    return best


def main():
    opt = dict(a.split("=") for a in sys.argv[1:])
    L, emin, top = int(opt.get("L", 5)), int(opt.get("emin", 8)), int(opt.get("top", 25))
    lo, hi, relax = int(opt.get("lo", 7)), int(opt.get("hi", 34)), int(opt.get("relax", 1))
    cs = [int(v) for v in opt.get("cs", "19,20,21,22,23").split(",")]
    if "sums" in opt:
        sums = [int(v) for v in opt["sums"].split(",")]
        pr = price(sums, emin, relax)
        assert pr is not None, "no chains for this set"
        es = list(range(2 * ((emin + 1) >> 1), 22))
        ch = {str(e): [sums[i] for i in range(len(sums)) if (mk >> i) & 1] for e, mk in zip(es, pr[4])}
        print("cost %d swaps %d misplaced %d chains %s" % (pr[0], pr[1], pr[2], ch))
        with open(opt["dump"], "w") as fh:
            json.dump({"sums": sums, "chains": ch}, fh)
        return
    res = []
    n_try = 0
    # the chain of e = 20 (s = 1) is B itself when the chain of e = 21 is empty: two sums c, c - 1 in this order
    for c in cs:
        for ia, ib in itertools.combinations(range(L), 2):
            others = [i for i in range(L) if i not in (ia, ib)]
            for rest in itertools.product(range(lo, hi + 1), repeat=L - 2):
                sums = [0] * L
                sums[ia], sums[ib] = c, c - 1
                for i, v in zip(others, rest):
                    sums[i] = v
                if len(set(sums)) < L:
                    continue
                n_try += 1
                ok = True
                for e in (19, 18, 17, 16, 13, 10, 8):          # quick filter before the full pricing
                    if e >= emin and not options(sums, e, relax):
                        ok = False
                        break
                if not ok:
                    continue
                pr = price(sums, emin, relax)
                if pr is not None:
                    res.append((pr[0], pr[1], pr[2], sums, pr[3], pr[4]))
                    print("found cost %d swaps %d misplaced %d sums %s" % (pr[0], pr[1], pr[2], sums), flush=True)
    res.sort(key=lambda r: (r[0], r[1]))
    print("tried", n_try, "feasible", len(res), flush=True)
    for cost, used, bad, sums, B, masks in res[:top]:
        es = list(range(2 * ((emin + 1) >> 1), 22))
        names = {e: [sums[i] for i in range(L) if (mk >> i) & 1] for e, mk in zip(es, masks)}
        print("cost %d swaps %d misplaced %d sums %s B %s chains %s" % (
            cost, used, bad, sums, [sums[i] for i in range(L) if (B >> i) & 1], names), flush=True)


if __name__ == "__main__":
    main()
