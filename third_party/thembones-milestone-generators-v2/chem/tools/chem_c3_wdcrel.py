"""chem_c3_wdcrel.py: builder C3 (third holder, Nov 17 15:05Z). R63-CHEM-WIRED's owner tool (tools/r63_chem_wired_dc.py,
imported read-only, run unchanged through its main()) for RELATION-BEARING wired loaders (WR7's P family: 24 or 23
generic outputs on 30 wires). The tool's free units are single-wire flips mapped through G^-1 (ginv_masks), which
asserts one pivot per wire: with a bit-1 relation the wires are DEPENDENT (30 wires, 24 outputs) and the assert fails.
Here two functions of that module are replaced before main() runs (marked C3):
  free_bits      : at every dead address a (low point >= lmin, below the last block), the wires of the loader's dead
                   words there (bits 1-5 of every word j whose cell (a, j) of an owned batch is dead) span a flip space
                   S_a; a flip f is loadable iff it keeps every relation of the loader (H f = 0, H = left null space of
                   G over the generic wires). Each vector of a basis of S_a cap ker H is one unit: its output mask
                   g = G^-1 f (unique: G has full column rank), filed under the lowest output of g (the tool solves
                   outputs in order and a unit changes only its own and later outputs), with the list of (column, bit)
                   flips it makes in the words.
  modified_words : applies a unit's flip list (bit 6 follows bit 5, as the tool does).
For a loader without relations every unit is a single wire and the result is the tool's own (one unit per wire).
   CHEM_BROT=7 CHEM_HB=77 python -I chem_c3_wdcrel.py <same arguments as r63_chem_wired_dc.py>"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import r63_chem_wired_dc as m  # noqa: E402


def _relations(G):
    """left null space of the wire rows G[q] (ints over outputs): list of wire masks h with XOR_{q in h} G[q] = 0."""
    piv = {}                            # pivot bit -> (row value, wire mask)
    rel = []
    for q, g in enumerate(G):
        v, t = g, 1 << q
        while v:
            b = v.bit_length() - 1
            if b not in piv:
                piv[b] = (v, t)
                break
            pv, pt = piv[b]
            v ^= pv
            t ^= pt
        if not v:
            rel.append(t)
    return rel


def _solver(G):
    """returns f -> x with G x = f over GF(2) (x: output mask), or None when f is not in the image of G."""
    nw = len(G)
    nout = max((g.bit_length() for g in G), default=0)
    rows = [(G[q], 1 << q) for q in range(nw)]        # equation q: parity(G[q] & x) = f_q
    piv = []                                          # (pivot bit, row value, combination of equations)
    for v, t in rows:
        for b, pv, pt in piv:
            if (v >> b) & 1:
                v ^= pv
                t ^= pt
        if v:
            b = v.bit_length() - 1
            for i, (b2, pv2, pt2) in enumerate(piv):     # keep reduced
                if (pv2 >> b) & 1:
                    piv[i] = (b2, pv2 ^ v, pt2 ^ t)
            piv.append((b, v, t))
    assert len(piv) == nout, "G must have full column rank"
    rel = _relations(G)

    def solve(f):
        for h in rel:                                  # f must keep every relation
            if bin(h & f).count("1") & 1:
                return None
        x = 0
        for b, v, t in piv:                            # x_b = parity of f over the equations combined into row b
            if bin(t & f).count("1") & 1:
                x |= 1 << b
        # the reduced rows have one pivot each and no other pivot bits: v = 1 << b plus non-pivot bits; with full
        # column rank every bit is a pivot, so v == 1 << b
        return x
    for b, v, t in piv:
        assert v == (1 << b), "reduction left a non-pivot bit"
    return solve, rel


def neighbour_relations(words, k2):
    """C3 (15:40Z): relations of loader k2 over its wires, as {(word j, bit b): relation bit mask}. A dead word of a
    batch is read by TWO loaders: a flip that keeps the owner's relations but breaks a relation of the other reader
    gives that loader one more generic output (RK7 run 24: k14 owner, k13 rebuilt at 233 instead of 226 [M]). A wire
    whose row is zero is a relation by itself (a flip there would make it generic)."""
    import chem_ld2_c4 as c4
    install_vset_rel()
    if k2 < 1 or k2 > 14:
        return {}, 0
    V2 = c4.vset(words, k2)
    an2 = core.analyze(c4.pad_table(words, {k2: V2}), k2, 4)
    G2 = [int(g) for g in an2.G]
    rel2 = _relations(G2)
    sig = {}
    for q, wjb in enumerate(an2.out_wires):
        s = 0
        for r, hm in enumerate(rel2):
            if (hm >> q) & 1:
                s |= 1 << r
        sig[tuple(wjb)] = s
    return sig, len(rel2)


def free_bits(words, edges, k, an, V, own, lmin, near=None):
    """C3: relation-keeping units (see the module docstring). Same return as the tool's free_bits."""
    assert near is None, "chem_c3_wdcrel: --near is not supported"
    A0, ncol = words.shape
    J = core.J
    nr = an.nr
    dead = m.dead_mask(words, edges, core.BROT)
    qidx = {w: q for q, w in enumerate(an.out_wires)}
    G = [int(g) for g in an.G]
    solve, rel = _solver(G)
    # C3 (15:40Z): the other reader of each owned batch: batch k (cur) -> loader k + 1, batch k - 1 (prev) -> k - 1
    nsig = {}
    off = len(rel)
    for nb_k, bat in ((k + 1, k), (k - 1, k - 1)):
        if 1 <= nb_k <= 14 and 0 <= bat:
            s, nrel = neighbour_relations(words, nb_k)
            nsig[bat] = (s, off, nb_k, nrel)
            off += nrel
    print(f"C3 wdcrel: relations kept: own {len(rel)}; "
          + ", ".join(f"batch {b} -> k{x[2]}: {x[3]}" for b, x in sorted(nsig.items())), flush=True)
    L = an.L
    dl = []
    for l in range(L):
        v = 0
        for mm in range(1, L):
            if (mm & l) == l:
                v ^= an.mono_coords[mm]
        dl.append(v)
    last = an.H - 1
    cols = []
    if "cur" in own:
        cols += [(J * k + j, j) for j in range(J)]
    if "prev" in own and k > 0:
        cols += [(J * (k - 1) + j, j) for j in range(J)]
    free, nfree, nunits = {}, 0, 0
    for a in range(min(A0, dead.shape[0])):
        h, l = a >> 4, a & 15
        if l < lmin or h >= last:
            continue
        dcol = {}                                       # word j -> the dead column used for it (cur first)
        for col, j in cols:
            if dead[a, col] and j not in dcol:
                dcol[j] = col
        if not dcol:
            continue
        cand = []                                       # (wire q, column, bit, word)
        for j, col in sorted(dcol.items()):
            for b in range(1, core.BROT - 1):
                q = qidx.get((j, b))
                if q is not None and G[q]:
                    cand.append((q, col, b, j))
        if not cand:
            continue
        # kernel of the relations restricted to the candidate wires: c with XOR_i c_i [H applied to e_qi] = 0;
        # C3 15:40Z: H also holds the relations of the other reader of the flipped word's batch
        hv = []
        for q, col, b, j in cand:
            s = 0
            for r, hm in enumerate(rel):
                if (hm >> q) & 1:
                    s |= 1 << r
            ns = nsig.get(col // J)
            if ns is not None:
                s |= ns[0].get((j, b), 0) << ns[1]
            hv.append(s)
        basis = []                                      # Gaussian elimination: (value, combination)
        kern = []
        for i, s in enumerate(hv):
            v, t = s, 1 << i
            for bv, bt in basis:
                if v & (1 << (bv.bit_length() - 1)):
                    v ^= bv
                    t ^= bt
            if v:
                basis.append((v, t))
                basis.sort(key=lambda x: -(x[0].bit_length()))      # leading bits in decreasing order
            else:
                kern.append(t)
        vec = np.zeros(nr, np.uint32)
        for M in range(1, nr + 1):
            if (M & h) == h:
                vec[M - 1] = dl[l]
        if (V is not None) and ((V >> h) & 1) and last <= nr:
            vec[last - 1] ^= np.uint32(dl[l])
        if not vec.any():
            continue
        for t in kern:
            f, flips = 0, []
            for i, (q, col, b, _j) in enumerate(cand):
                if (t >> i) & 1:
                    f ^= 1 << q
                    flips.append((col, b))
            x = solve(f)
            assert x is not None and x, "C3: a kernel flip that the outputs cannot carry"
            o = (x & -x).bit_length() - 1               # lowest output of g
            free.setdefault(o, []).append(dict(vec=vec, g=x, a=a, flips=flips, col=flips[0][0], bit=flips[0][1],
                                               h=h, l=l))
            nfree += 1
            nunits += len(flips) > 1
    print(f"C3 wdcrel: relations of the loader {len(rel)}; units {nfree} ({nunits} of them multi-wire)", flush=True)
    owned = [c for c, _ in cols]
    return free, nfree, int(dead[:, owned].sum())


def modified_words(words, used):
    """C3: a unit flips every (column, bit) of its list; bit 6 follows bit 5 (the carrier stays)."""
    w = words.copy()
    tb = core.BROT - 1
    for o, lst in used.items():
        for f in lst:
            for col, b in f.get("flips", [(f["col"], f["bit"])]):
                w[f["a"], col] ^= (1 << b)
                if b == tb - 1:
                    w[f["a"], col] ^= (1 << tb)
    return w


def install_vset_rel():
    """C3 (15:20Z): the block set V as C2's chem_ld2_e19grp.block_set makes it (R59's tie engine uses that): c4.vset,
    then gc.vset_rel when relations among the outputs constrain V. R63's main calls c4.vset alone, which on a relation
    loader gives a V whose closure never solves (RK7 run 24 k14: 'no closed solve for group [0..5]' even with --nofree,
    while R59's tie engine closes it at 152 [M]). Installed by replacing chem_ld2_c4.vset for this process."""
    import chem_ld2_c4 as c4
    import chem_ld2_e17gc4 as gc
    if getattr(c4.vset, "_c3", False):
        return
    orig = c4.vset

    def vset(words, k):
        V = orig(words, k)
        an = core.analyze(c4.pad_table(words, {k: V}), k, 4)
        try:
            V2 = gc.vset_rel(words, k, an)
        except Exception:  # noqa: BLE001  (block_set does the same)
            V2 = None
        return V2 if V2 is not None else V
    vset._c3 = True
    c4.vset = vset


def install(mod):
    """C3 (15:45Z): put the three C3 changes into an owner-tool module with R63's main (R63's own, or builder C2's
    pinned copy chem_ld2_e22dcr, with or without C2's annealed last group chem_ld2_e22dci)."""
    install_vset_rel()
    mod.free_bits = free_bits
    mod.modified_words = modified_words
    return mod


if __name__ == "__main__":
    # C3_BASE: r63 (default: tools/r63_chem_wired_dc.py), e22dcr (C2's pinned copy), e22dci (C2's copy + its annealed
    # last group: DCI_DUMP / DCI_INJECT as in tools/chem_ld2_e22dci.py)
    base = os.environ.get("C3_BASE", "r63")
    if base == "e22dci":
        import chem_ld2_e22dci as dci  # noqa: E402  (installs its make_tie into chem_ld2_e22dcr at import)
        mod = dci.dc
    elif base == "e22dcr":
        import chem_ld2_e22dcr as mod  # noqa: E402
    else:
        mod = m
    install(mod).main()
