"""chem_c3_vrun.py: builder C3 (fourth holder, Nov 17 18:50Z). Run a loader tool's main (e.g. design seat
R65-CHEM-WIRED221's r65_chem_wired221_dc.py / _build.py, read-only) with builder C2's relation-aware block set V
installed in place of chem_ld2_c4.vset (chem_c3_wdcrel.install_vset_rel). WHY [M box 6 gtest, 18:40Z]: on WR7's G
tables (a wire's bit-0 delta = parity(address AND m) XOR g(address >> 4): an 'S' wire) R65's sweep stops in the
emitter with "wire (1, 0): pure-high part outside H (tau != kappa)", while my owner tool (which installs this V)
emits the same loader with every wire equal to its model.
19:05Z: ALSO (C3_NREL=1, the default) R65's units keep the OTHER READER's relations: a unit flips one (word j, bit b)
of batch K-1 or K, which loader K-1 or K+1 reads too; a unit whose (j, b) carries a relation of that reader (its
bit-1 relations on WR7's / T10's P tables, its bit-0 wires) is dropped (chem_c3_wdcrel.neighbour_relations, the test
my own owner tool uses). On tables without such relations nothing is dropped (R65's units never flip bit 0).
   python -I chem_c3_vrun.py r65_chem_wired221_dc TABLE --k 14 ...   (the module name, then its own arguments)"""
import os
import runpy
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_c3_wdcrel as w  # noqa: E402

w.install_vset_rel()


def keep_neighbour_relations(D):
    """wrap D.free_bits (R65's unit builder) so that no unit breaks a relation of the other reader of its batch"""
    import chem_loaders_core as core
    orig = D.free_bits

    def free_bits(words, edges, k, an, V, own, lmin, live=None):
        free, nfree, nlive, ndead = orig(words, edges, k, an, V, own, lmin, live)
        J = core.J
        bad = {}
        for nb_k, bat in ((k + 1, k), (k - 1, k - 1)):
            if 1 <= nb_k <= 14 and bat >= 0:
                sig, _n = w.neighbour_relations(words, nb_k)
                bad[bat] = {tuple(jb) for jb, s in sig.items() if s}
        dropped = {"dead": 0, "live": 0}
        for o in list(free):
            keep = []
            for f in free[o]:
                bat = f["col"] // J
                if bat in bad and (f["col"] % J, f["bit"]) in bad[bat]:
                    dropped["live" if f.get("cost", 0.0) != 0.0 else "dead"] += 1
                    continue
                keep.append(f)
            free[o] = keep
        print(f"C3 vrun: units dropped for the other reader's relations: {dropped} "
              f"(relation wires: {', '.join(f'batch {b}: {len(s)}' for b, s in sorted(bad.items()))})", flush=True)
        return free, nfree - dropped["dead"], nlive - dropped["live"], ndead
    D.free_bits = free_bits


def relation_owner_units(D):
    """C3 (19:25Z): R65's unit builder for a RELATION-BEARING owner (24 / 23 generic outputs: R65's own builder asserts
    one pivot per wire there): my relation-keeping dead units (chem_c3_wdcrel.free_bits: kernels of the loader's and the
    other readers' relations, multi-wire) plus R65's priced live units on the wires that carry no relation of the
    loader nor of the batch's other reader (a single flip keeps every relation there), each mapped to its outputs by
    the loader's own GF(2) solve. A plain loader keeps R65's builder (wrapped by keep_neighbour_relations)."""
    import numpy as np
    import chem_loaders_core as core
    import r63_chem_wired_dc as r63
    from r62_chem_cache_dead import dead_mask
    r63.modified_words = w.modified_words          # handles multi-wire units ('flips') and single ones alike
    plain = D.free_bits

    def free_bits(words, edges, k, an, V, own, lmin, live=None):
        G = [int(g) for g in an.G]
        if not w._relations(G):
            return plain(words, edges, k, an, V, own, lmin, live)
        solve, rel = w._solver(G)
        free, nfree, ndead = w.free_bits(words, edges, k, an, V, own, lmin)
        for o in free:
            for f in free[o]:
                f.setdefault("cost", 0.0)
        nlive = 0
        if live is not None:
            tau, maxbit, price, n1 = live
            J = core.J
            A0 = words.shape[0]
            last = an.H - 1
            dead = dead_mask(words, edges, core.BROT)
            qidx = {tuple(wb): q for q, wb in enumerate(an.out_wires)}
            osig = {q: any((hm >> q) & 1 for hm in rel) for q in range(len(G))}
            nbad = {}
            for nb_k, bat in ((k + 1, k), (k - 1, k - 1)):
                if 1 <= nb_k <= 14 and bat >= 0:
                    sig, _n = w.neighbour_relations(words, nb_k)
                    nbad[bat] = {tuple(jb) for jb, s in sig.items() if s}
            cols = []
            if "cur" in own:
                cols += [(J * k + j, j) for j in range(J)]
            if "prev" in own and k > 0:
                cols += [(J * (k - 1) + j, j) for j in range(J)]
            dl = []
            for l in range(an.L):
                v = 0
                for m in range(1, an.L):
                    if (m & l) == l:
                        v ^= an.mono_coords[m]
                dl.append(v)
            for j in range(J):
                cj = [c for c, jj in cols if jj == j]
                for a in range(n1, A0):
                    h, l = a >> 4, a & 15
                    if l < lmin or h >= last or dead[a, cj].any():
                        continue
                    vec = None
                    for b in range(1, maxbit + 1):
                        q = qidx.get((j, b))
                        if q is None or not G[q] or osig[q]:
                            continue
                        best = None
                        for c in cj:
                            if (c // J) in nbad and (j, b) in nbad[c // J]:
                                continue
                            p = float(price[b][c][a - n1])
                            if best is None or p < best[0]:
                                best = (p, c)
                        if best is None or best[0] > tau:
                            continue
                        if vec is None:
                            vec = np.zeros(an.nr, np.uint32)
                            for M in range(1, an.nr + 1):
                                if (M & h) == h:
                                    vec[M - 1] = dl[l]
                            if (V is not None) and ((V >> h) & 1) and last <= an.nr:
                                vec[last - 1] ^= np.uint32(dl[l])
                            if not vec.any():
                                break
                        x = solve(1 << q)
                        if not x:
                            continue
                        o = (x & -x).bit_length() - 1
                        free.setdefault(o, []).append(dict(vec=vec, g=x, a=a, col=best[1], bit=b, h=h, l=l,
                                                           cost=best[0], flips=[(best[1], b)]))
                        nlive += 1
            for o in free:
                free[o].sort(key=lambda f: f["cost"])
        print(f"C3 vrun: relation-bearing owner k{k}: {len(rel)} relations; dead units {nfree}, live units {nlive}",
              flush=True)
        return free, nfree, nlive, ndead
    D.free_bits = free_bits


def placed_order_all_wires(D):
    """C3 (20:08Z): R65's output placement (--oauto) places sum(plan sizes) outputs; a relation-bearing loader has one
    bit 1-5 wire more than generic outputs (24 / 23 of 25), and R65's order left that wire out of the rebase order.
    The wires not placed go right after the placed ones (the dependent wire then sits after the generic basis)."""
    orig = D.placed_order

    def placed_order(an, sup, sizes, pattern="", explicit=""):
        order, seq = orig(an, sup, sizes, pattern, explicit)
        if len(set(order)) == len(an.out_wires):
            return order, seq
        wires = {w2: q for q, w2 in enumerate(an.out_wires)}
        placed = set(order)
        missing = [q for q, w2 in enumerate(an.out_wires) if q not in placed]
        nfirst = sum(1 for w2 in an.out_wires if w2[1] == 0)
        cut = nfirst + len(seq)
        new = order[:cut] + missing + order[cut:]
        print(f"C3 vrun: placement: {len(missing)} wire(s) not placed appended after the placed ones: "
              f"{[an.out_wires[q] for q in missing]}", flush=True)
        assert sorted(new) == list(range(len(an.out_wires)))
        return new, seq
    D.placed_order = placed_order


mod = sys.argv[1]
args = sys.argv[2:]
if os.environ.get("C3_NREL", "1") == "1" and mod.startswith("r65_chem_wired221_"):
    import r65_chem_wired221_dc as D
    keep_neighbour_relations(D)
    relation_owner_units(D)
    placed_order_all_wires(D)
    if mod == "r65_chem_wired221_dc":
        sys.argv = [D.__file__] + args
        D.main()
        sys.exit(0)
sys.argv = [os.path.join(os.path.dirname(os.path.abspath(__file__)), mod + ".py")] + args
runpy.run_module(mod, run_name="__main__")
