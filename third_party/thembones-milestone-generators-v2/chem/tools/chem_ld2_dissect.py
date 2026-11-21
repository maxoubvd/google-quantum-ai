"""chem_ld2_dissect: structure of the products of one ANGLE_k loader circuit (data only; the circuit is parsed and
simulated classically here, nothing from it is executed).

Usage: python -I chem_ld2_dissect.py <dir with submission.npz and angle_k_*.kmb> <k> [--json out.json]

For circuit k, forward direction: selection holds the address a (bit i on selection wire i), the word wires hold
the words of batch k-1 (zero for k = 0); at the end they must hold batch k. Every CCX is classified:
  * infrastructure: both controls are functions of the high bits only, or of the low bits only;
  * product: anything else. Its controls are decomposed as A(h) + a(l) and B(h) + b(l) (checked on every valid
    address), and its fan-out is the set of word bits that hold it at the end (symbolic XOR tracking).
Printed: counts, the low parts used, the blocks of low monomials they live in, and per block the dimension of
the span of its features (y (x) s), i.e. the plan's codimension e_b = 3 * outputs - dim.
"""
import json
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402

NADDR_FULL = 1232  # 77 high values x 16 low values


def mobius(arr, nbits):
    sh = arr.shape[:-1]
    for i in range(nbits):
        v = arr.reshape(sh + (-1, 2, 1 << i))
        v[..., 1, :] ^= v[..., 0, :]
    return arr


def to_int(bits_arr):
    return int.from_bytes(np.packbits(np.asarray(bits_arr, np.uint8), bitorder="little").tobytes(), "little")


def rank_gf2(vecs):
    piv = {}
    r = 0
    for v in vecs:
        while v:
            hb = v.bit_length() - 1
            if hb in piv:
                v ^= piv[hb]
            else:
                piv[hb] = v
                r += 1
                break
    return r


def load(d, k):
    with np.load(os.path.join(d, "submission.npz"), allow_pickle=False) as z:
        words = np.array(z["angle_words"]).astype(np.int64)
    c = kmb.read(os.path.join(d, f"angle_k_{k}.kmb"), mmap=False)
    return words, c


def simulate(words, c, k):
    rows = words.shape[0]
    regs = {name: ids for name, ids, _ in c.registers}
    sel = [int(x) for x in regs["selection"]]
    wid = [[int(x) for x in regs[f"word_{j}"]] for j in range(5)]
    a = np.arange(NADDR_FULL)
    val = np.zeros((c.num_qubits, NADDR_FULL), np.uint8)
    for i, q in enumerate(sel):
        val[q] = (a >> i) & 1
    prev = np.zeros((NADDR_FULL, 5), np.int64)
    if k > 0:
        prev[:rows] = words[:, 5 * (k - 1):5 * k]
    for j in range(5):
        for b in range(8):
            val[wid[j][b]] = (prev[:, j] >> b) & 1
    sym = [0] * c.num_qubits
    valid = a < rows
    H = (rows + 15) // 16

    def parts(v):
        lo = v[:16].astype(np.uint8).copy()
        mobius(lo, 4)
        hi = np.zeros(128, np.uint8)
        hi[:H] = v[::16][:H]
        mobius(hi, 7)
        c0 = int(v[0])
        lo_i = to_int(lo) & ~1
        hi_i = to_int(hi[:H]) & ~1
        # check v == hi(h) + lo(l) + c0 on valid addresses
        h = a >> 4
        l = a & 15
        recon = (v[h * 16] ^ v[l] ^ c0) & 1  # f(h,0) + f(0,l) + f(0,0)
        pure_sum = bool(np.all(recon[valid] == v[valid]))
        return lo_i, hi_i, c0, pure_sum

    ops = c.base
    i3 = i2 = i1 = 0
    infra = []
    prods = []
    nonsep = 0
    for op in ops.tolist():
        n = (op & 0x30) >> 4 if op != kmb.PUSH else 0
        if op == kmb.PUSH or op == kmb.INFO_TAG or op == kmb.POP:
            # PUSH has one bit operand; INFO_TAG / POP none (classical bookkeeping)
            pass
        if n == 3:
            q0, q1, q2 = (int(c.q3[x][i3]) & kmb.ID_MASK for x in range(3))
            i3 += 1
            if op == kmb.CCX:
                l1, h1, c1, s1 = parts(val[q0])
                l2, h2, c2, s2 = parts(val[q1])
                if not (s1 and s2):
                    nonsep += 1
                is_infra = (l1 == 0 and l2 == 0) or (h1 == 0 and h2 == 0)
                if is_infra:
                    infra.append((q0, q1, q2))
                else:
                    prods.append(dict(l1=l1, h1=h1, c1=c1, l2=l2, h2=h2, c2=c2, sep=bool(s1 and s2),
                                      chained=bool(sym[q0] or sym[q1])))
                    sym[q2] ^= 1 << (len(prods) - 1)
                val[q2] ^= val[q0] & val[q1]
            elif op == kmb.CCZ:
                pass
            elif op == kmb.CSWAP:
                raise SystemExit("CSWAP not handled")
        elif n == 2:
            q0, q1 = (int(c.q2[x][i2]) & kmb.ID_MASK for x in range(2))
            i2 += 1
            if op == kmb.CX:
                val[q1] ^= val[q0]
                sym[q1] ^= sym[q0]
            elif op == kmb.SWAP:
                val[[q0, q1]] = val[[q1, q0]]
                sym[q0], sym[q1] = sym[q1], sym[q0]
        elif n == 1:
            q0 = int(c.q1[i1]) & kmb.ID_MASK
            i1 += 1
            if op == kmb.X:
                val[q0] ^= 1
            elif op == kmb.HMR:
                val[q0] = 0
                sym[q0] = 0
    # check the words
    target = np.zeros((NADDR_FULL, 5), np.int64)
    target[:rows] = words[:, 5 * k:5 * k + 5]
    ok = True
    for j in range(5):
        got = np.zeros(NADDR_FULL, np.int64)
        for b in range(8):
            got |= val[wid[j][b]].astype(np.int64) << b
        if not np.array_equal(got[:rows], target[:rows, j]):
            ok = False
    fan = []
    wires = [wid[j][b] for j in range(5) for b in range(8)]
    for p in range(len(prods)):
        f = 0
        for t, q in enumerate(wires):
            if (sym[q] >> p) & 1:
                f |= 1 << t
        fan.append(f)
    return dict(infra=len(infra), prods=prods, fan=fan, ok=ok, nonsep=nonsep, nq=c.num_qubits)


def main():
    d, k = sys.argv[1], int(sys.argv[2])
    out = sys.argv[sys.argv.index("--json") + 1] if "--json" in sys.argv else None
    words, c = load(d, k)
    r = simulate(words, c, k)
    prods, fan = r["prods"], r["fan"]
    print(f"k={k} nq={r['nq']} infra={r['infra']} products={len(prods)} words_ok={r['ok']} "
          f"nonseparable_controls={r['nonsep']} chained={sum(p['chained'] for p in prods)} "
          f"zero_fanout={sum(1 for f in fan if f == 0)}")
    lows = Counter()
    for p in prods:
        lows[p["l1"]] += 1
        lows[p["l2"]] += 1
    print("distinct low parts:", len(lows), "zero low parts:", lows.get(0, 0))
    # blocks: connected components of low monomials appearing together in one low part
    parent = list(range(16))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for lp in lows:
        ms = [m for m in range(1, 16) if (lp >> m) & 1]
        for m in ms[1:]:
            parent[find(m)] = find(ms[0])
    comp = {}
    for m in range(1, 16):
        comp.setdefault(find(m), []).append(m)
    blocks = sorted(comp.values())
    print("blocks (low monomials, as 4-bit masks):", [[format(m, "04b") for m in b] for b in blocks])
    bidx = {}
    for i, b in enumerate(blocks):
        for m in b:
            bidx[m] = i

    def block_of(lp):
        ms = [m for m in range(1, 16) if (lp >> m) & 1]
        bs = {bidx[m] for m in ms}
        return bs.pop() if len(bs) == 1 else (-1 if not bs else -2)
    pair_kind = Counter()
    feats = {i: [] for i in range(len(blocks))}
    for p, f in zip(prods, fan):
        b1, b2 = block_of(p["l1"]), block_of(p["l2"])
        pair_kind[(min(b1, b2), max(b1, b2))] += 1
        # features: (l2, f) with coefficient h1; (l1, f) with coefficient h2
        if p["l2"] and p["h1"]:
            feats[b2].append((p["l2"], f))
        if p["l1"] and p["h2"]:
            feats[b1].append((p["l1"], f))
    print("pairs by blocks (block1, block2): count:", dict(sorted(pair_kind.items())))
    nout = 40
    tot_codim = 0
    for i, b in enumerate(blocks):
        # feature y (x) s as a vector over (monomial position in block, output bit)
        pos = {m: j for j, m in enumerate(b)}
        vecs = []
        for y, s in feats[i]:
            v = 0
            for m in b:
                if (y >> m) & 1:
                    v |= s << (nout * pos[m])
            vecs.append(v)
        rk = rank_gf2(vecs)
        fs = rank_gf2([s for _, s in feats[i]])
        print(f"block {i} {[format(m, '04b') for m in b]}: features {len(feats[i])} span dim {rk} "
              f"(distinct y {len({y for y, _ in feats[i]})}, fan-out rank {fs})")
        tot_codim += len(b) * 35 - rk
    print("sum over blocks of (3 x 35 - span dim):", tot_codim)
    if out:
        with open(out, "w") as fh:
            json.dump(dict(k=k, infra=r["infra"], prods=prods, fan=fan, blocks=blocks), fh)


if __name__ == "__main__":
    main()
