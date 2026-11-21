"""chem_ld2_c4pub: STUDY (data only; builder C2, seventh holder). How a public 68-AND loader meets the C-4 conditions.

Usage: python -I chem_ld2_c4pub.py <dir with submission.npz and angle_k_*.kmb> <k> [<k> ...]

Per circuit k (the circuit is parsed and simulated classically by chem_ld2_dissect; nothing from it is executed):
  * infra ANDs split into high-only and low-only (by the high / low parts of both controls);
  * V = the block set of chem_ld2_c4test (with block 76) for the table's delta of batch k;
  * per product: are both high factors (constant included) in H = {f : sum_V f = 0}?  sigma = sum_V H1 H2;
  * tau_b = XOR of sigma over the products feeding payload bit b (bits 0-6 of the 5 words), against
    kappa_b = the table's V-sum of the delta at low row 0; plus the top bit (bit 7) of each word: products feeding it.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_ld2_dissect as ds  # noqa: E402
import chem_ld2_c4test as ct  # noqa: E402
import kmb  # noqa: E402
from chem_ld2_hb import HB, LAST  # noqa: E402  (tenth holder: 77 / 76 unless CHEM_HB is set)


def truth(anf, H=HB):
    """Truth table over the H high blocks (bit p = value at block p) of a function given by ANF coefficients."""
    t = 0
    for h in range(H):
        s = 0
        m = anf
        while m:
            b = m & -m
            i = b.bit_length() - 1
            if (h & i) == i:
                s ^= 1
            m ^= b
        t |= s << h
    return t


def infra_split(words, c, k):
    """Re-simulate and classify the infra ANDs (high-only / low-only); returns (n_high, n_low)."""
    rows = words.shape[0]
    regs = {name: ids for name, ids, _ in c.registers}
    sel = [int(x) for x in regs["selection"]]
    wid = [[int(x) for x in regs[f"word_{j}"]] for j in range(5)]
    a = np.arange(ds.NADDR_FULL)
    val = np.zeros((c.num_qubits, ds.NADDR_FULL), np.uint8)
    for i, q in enumerate(sel):
        val[q] = (a >> i) & 1
    prev = np.zeros((ds.NADDR_FULL, 5), np.int64)
    if k > 0:
        prev[:rows] = words[:, 5 * (k - 1):5 * k]
    for j in range(5):
        for b in range(8):
            val[wid[j][b]] = (prev[:, j] >> b) & 1
    H = (rows + 15) // 16

    def parts(v):
        lo = v[:16].astype(np.uint8).copy()
        ds.mobius(lo, 4)
        hi = np.zeros(128, np.uint8)
        hi[:H] = v[::16][:H]
        ds.mobius(hi, 7)
        return ds.to_int(lo) & ~1, ds.to_int(hi[:H]) & ~1

    nh = nl = 0
    hi_funcs = []
    i3 = i2 = i1 = 0
    for op in c.base.tolist():
        n = (op & 0x30) >> 4 if op != kmb.PUSH else 0
        if n == 3:
            q0, q1, q2 = (int(c.q3[x][i3]) & kmb.ID_MASK for x in range(3))
            i3 += 1
            if op == kmb.CCX:
                l1, h1 = parts(val[q0])
                l2, h2 = parts(val[q1])
                if l1 == 0 and l2 == 0:
                    nh += 1
                    hi_funcs.append(ds.to_int(val[q2][::16][:H] ^ (val[q0][::16][:H] & val[q1][::16][:H])))
                elif h1 == 0 and h2 == 0:
                    nl += 1
                val[q2] ^= val[q0] & val[q1]
        elif n == 2:
            q0, q1 = (int(c.q2[x][i2]) & kmb.ID_MASK for x in range(2))
            i2 += 1
            if op == kmb.CX:
                val[q1] ^= val[q0]
            elif op == kmb.SWAP:
                val[[q0, q1]] = val[[q1, q0]]
        elif n == 1:
            q0 = int(c.q1[i1]) & kmb.ID_MASK
            i1 += 1
            if op == kmb.X:
                val[q0] ^= 1
            elif op == kmb.HMR:
                val[q0] = 0
    return nh, nl


def vset(words, k):
    import chem_loaders_core as _core                 # tenth holder: payload bits = word width - 1 (7 at 8 bits)
    F, valid = ct.delta_bits(words, k, _core.BROT - 1)
    rows = [np.ones(HB, np.uint8)]
    pidx = np.arange(HB)
    for j in range(7):
        rows.append(((pidx >> j) & 1).astype(np.uint8))
    for b in range(F.shape[0]):
        for l in range(1, 16):
            if not valid[LAST, l]:
                continue
            rows.append(F[b, :, l] ^ F[b, :, 0])
    N = ct.nullspace_gf2(np.array(rows, np.uint8))
    sols = [v for v in N if v[LAST]]
    if not sols:
        return None
    v = sols[0]
    return sum(1 << p for p in range(HB) if v[p])


def main():
    d = sys.argv[1]
    for k in (int(x) for x in sys.argv[2:]):
        words, c = ds.load(d, k)
        r = ds.simulate(words, c, k)
        nh, nl = infra_split(words, c, k)
        V = vset(words, k)
        print(f"k={k} products {len(r['prods'])} infra {r['infra']} (high {nh}, low {nl}) words_ok {r['ok']}")
        if V is None:
            print("  no V with block 76")
            continue
        Vs = [p for p in range(HB) if (V >> p) & 1]
        print(f"  V: {len(Vs)} blocks, xor {np.bitwise_xor.reduce(np.array(Vs))}")
        allH = (1 << HB) - 1
        sig = lambda t: bin(t & V).count("1") & 1  # noqa: E731
        out_of_H = 0
        nsig = 0
        tau = 0
        top_fed = 0
        for p, f in zip(r["prods"], r["fan"]):
            t1 = truth(p["h1"]) ^ (allH if p["c1"] else 0)
            t2 = truth(p["h2"]) ^ (allH if p["c2"] else 0)
            out_of_H += sig(t1) + sig(t2)
            s = sig(t1 & t2)
            nsig += s
            for j in range(5):
                if (f >> (8 * j + 7)) & 1:
                    top_fed += 1
            if s:
                for j in range(5):
                    for b in range(7):
                        if (f >> (8 * j + b)) & 1:
                            tau ^= 1 << (7 * j + b)
        F, valid = ct.delta_bits(words, k, 7)
        kappa = 0
        for b in range(35):
            kappa |= (int(np.bitwise_xor.reduce(F[b, Vs, 0])) & 1) << b
        print(f"  factors outside H: {out_of_H}; products with sigma=1: {nsig} of {len(r['prods'])}; "
              f"product fan-outs into a top bit: {top_fed}")
        print(f"  tau   {tau:035b}\n  kappa {kappa:035b}  equal {tau == kappa}  weight {bin(kappa).count('1')}")


if __name__ == "__main__":
    sys.exit(main())
