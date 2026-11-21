#!/usr/bin/env python3
"""invsqrt_fc_s18m9q.py -- FC eighteenth holder: with the three affine forms fixed, solve the quadratic part Q of
  target = Q + L1 L2 L3   on every reached pattern of a spec (GF(2) elimination; exact, no search)
and write the device for the generator in ROLE names (one gated CCZ on L1, L2, L3 and the CZ / Z / sign terms of Q).
  invsqrt_fc_s18m9q.py SPEC.json "L1" "L2" "L3" ROLEMAP OUT.json
  L = wire names joined by ^ (and ^1);  ROLEMAP = q12=t0,q13=t1,...  (spec wire name -> generator role)
Safe to run: reads the spec, writes OUT.json."""
import itertools
import json
import sys


def main():
    spec, f1, f2, f3, rmap, out = sys.argv[1:7]
    d = json.load(open(spec))
    names = d["inputs"]
    n = len(names)
    idx = {nm: j for j, nm in enumerate(names)}
    role = dict(kv.split("=") for kv in rmap.split(","))

    def form(s):
        ws, k = [], 0
        for t in s.split("^"):
            t = t.strip()
            if t == "1":
                k ^= 1
            else:
                ws.append(idx[t])
        return ws, k
    Ls = [form(f) for f in (f1, f2, f3)]
    mons = [()] + [(j,) for j in range(n)] + list(itertools.combinations(range(n), 2))
    rows = []
    for p, y in d["rows"]:
        p = int(p)
        b = [(p >> j) & 1 for j in range(n)]
        prod = 1
        for ws, k in Ls:
            prod &= (sum(b[j] for j in ws) + k) & 1
        rhs = (int(y) ^ prod) & 1
        vec = 0
        for i, m in enumerate(mons):
            if all(b[j] for j in m):
                vec |= 1 << i
        rows.append((vec, rhs))
    # Gaussian elimination over GF(2): rows (vec, rhs); find Q with vec . Q = rhs for all rows
    piv = {}
    for vec, rhs in rows:
        while vec:
            h = vec.bit_length() - 1
            if h not in piv:
                piv[h] = (vec, rhs)
                break
            pv, pr = piv[h]
            vec ^= pv
            rhs ^= pr
        else:
            if rhs:
                print("NOT quadratic after the product: the forms do not serve this target")
                sys.exit(1)
    q = 0
    for h in sorted(piv):                      # back substitution, free variables 0
        pv, pr = piv[h]
        low = pv & ((1 << h) - 1)
        val = pr ^ (bin(low & q).count("1") & 1)
        if val:
            q |= 1 << h
    for vec, rhs in rows:                      # verify on every reached pattern
        assert (bin(vec & q).count("1") & 1) == rhs
    terms = [mons[i] for i in range(len(mons)) if (q >> i) & 1]
    R = lambda j: role[names[j]]
    dev = dict(forms=[[[R(j) for j in ws], k] for ws, k in Ls],
               quad=[[R(m[0]), R(m[1])] for m in terms if len(m) == 2],
               lin=[R(m[0]) for m in terms if len(m) == 1],
               const=int(() in terms), spec=spec, rows=len(rows))
    json.dump(dev, open(out, "w"), indent=1)
    print("Q found and checked on all %d reached patterns: %d CZ terms, %d Z terms, sign %d; one CCZ on %s" % (
        len(rows), len(dev["quad"]), len(dev["lin"]), dev["const"], dev["forms"]))
    print("device", out)


if __name__ == "__main__":
    main()
