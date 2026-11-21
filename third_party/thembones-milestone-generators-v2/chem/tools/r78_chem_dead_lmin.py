"""r78_chem_dead_lmin.py: design seat R78-CHEM-DEAD. THE DEAD WORDS AT LOW POINTS 0-2 AS KILL UNKNOWNS OF AN OWNER.

R63-CHEM-WIRED's owner tool draws free bits from dead words at low points >= 3 only (--lmin 3: "keeps V and kappa of
the C-4 closure fixed"). The closure (tools/chem_ld2_c4.py, read-only) asks that the V-sum of every payload bit of a
loader's delta be the same on every low row; at low points 3..15 the padding row of block 76 follows, at 0, 1, 2
block 76 holds real rows. A flip at (block h, low point 0-2) leaves every V-sum as it was when h is in the block set
of NEITHER loader that reads the word (the owner's and the neighbour's). This tool
  count: prints, per owner, the dead words at low points 0-2 and how many stand at such a block;
  run:   runs the owner tool itself (imported, nothing copied or edited) with those bits added to its free bits; the
         tool's own checks (the block set V must not move, kappa must not move, the tensor model against the
         rewritten table, the emitter on every address) decide whether the reading holds.

    CHEM_BROT=7 CHEM_HB=77 python tools/r78_chem_dead_lmin.py count TABLE.npz [--owners 4,8,10,12,14]
    CHEM_BROT=7 CHEM_HB=77 python tools/r78_chem_dead_lmin.py run [--base] TABLE.npz --k 8 <the owner tool's arguments>
(--base: the owner tool as it is, same arguments: the control.)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("CHEM_BROT", "7")
os.environ.setdefault("CHEM_HB", "77")
import chem_loaders_core as core  # noqa: E402
import chem_ld2_c4 as c4  # noqa: E402
import r63_chem_wired_dc as dc  # noqa: E402
from r62_chem_cache_dead import dead_mask  # noqa: E402

J = 5
_orig_free_bits = dc.free_bits
STAT = {}


def vsets(words, ks):
    out = {}
    for k in ks:
        if 1 <= k < words.shape[1] // J:
            try:
                out[k] = c4.vset(words, k)
            except Exception as exc:      # a loader without a block set
                out[k] = None
    return out


def free_bits_ext(words, edges, k, an, V, own, lmin, near=None):
    free, nfree, ndead = _orig_free_bits(words, edges, k, an, V, own, lmin, near)
    if V is None:
        return free, nfree, ndead
    dead = dead_mask(words, edges, core.BROT)
    nb = words.shape[1] // J
    Vs = vsets(words, [k - 1, k + 1])
    qidx = {w: q for q, w in enumerate(an.out_wires)}
    ginv = dc.ginv_masks(an)
    L, nr = an.L, an.nr
    dl = []
    for l in range(L):
        v = 0
        for m in range(1, L):
            if (m & l) == l:
                v ^= an.mono_coords[m]
        dl.append(v)
    last = an.H - 1
    cols = []
    if "cur" in own:
        cols += [(J * k + j, j, k + 1) for j in range(J)]            # batch k is read by loaders k and k + 1
    if "prev" in own and k > 0:
        cols += [(J * (k - 1) + j, j, k - 1) for j in range(J)]      # batch k - 1 by loaders k - 1 and k
    seen = set()
    add = cells = barred = 0
    for col, j, other in cols:
        Vo = Vs.get(other) if 1 <= other < nb else 0                 # loader 0 / none: no second block set to keep
        for a in np.nonzero(dead[:, col])[0].tolist():
            h, l = a >> 4, a & 15
            if l >= lmin or h >= last or (a, j) in seen:
                continue
            if (V >> h) & 1 or Vo is None or (Vo >> h) & 1:
                barred += 1
                continue
            seen.add((a, j))
            vec = np.zeros(nr, np.uint32)
            for M in range(1, nr + 1):
                if (M & h) == h:
                    vec[M - 1] = dl[l]
            if not vec.any():
                continue
            cells += 1
            for b in range(1, core.BROT - 1):
                q = qidx.get((j, b))
                if q is None or q not in ginv:
                    continue
                o, gm = ginv[q]
                free.setdefault(o, []).append(dict(vec=vec, g=gm, a=a, col=col, bit=b, h=h, l=l))
                add += 1
    STAT.update(cells=cells, add=add, barred=barred)
    print(f"r78: low points 0-{lmin - 1}: {cells} dead addresses added ({add} free bits) at blocks outside both readers' "
          f"block sets; {barred} left out (inside one of them)", flush=True)
    return free, nfree + add, ndead


_real_vset = c4.vset
_VCACHE = {}


def vset_keep(words, k):
    """c4.vset, but on a rewritten table the block set of the first call is kept WHEN IT STILL HOLDS there (the
    solver returns the first vector of a null space; flips at low points 0-2 change that space, so its pick can move
    although the old set still satisfies every equation). The check is the builders' own check_padded."""
    v = _real_vset(words, k)
    if k not in _VCACHE:
        _VCACHE[k] = (v, words.copy())
        return v
    v0, w0 = _VCACHE[k]
    if v == v0 or v0 is None or np.array_equal(words, w0):
        return v
    bad = c4.check_padded(c4.pad_table(words, {k: v0}), {k: v0})[k]
    nb = words.shape[1] // J
    oth = {}
    for kk in (k - 1, k + 1):
        if 1 <= kk < nb:
            vk = _VCACHE[kk][0] if kk in _VCACHE else _real_vset(w0, kk)
            oth[kk] = None if vk is None else c4.check_padded(c4.pad_table(words, {kk: vk}), {kk: vk})[kk]
    print(f"r78: on the rewritten table the solver picks another block set for k{k}; the FIRST one checked there: "
          f"{bad} payload bits off (check_padded); the neighbours' first block sets: {oth}", flush=True)
    return v0 if bad == 0 else v


def free_bits_ext65(words, edges, k, an, V, own, lmin, live=None):
    """R65's free_bits + the dead words at low points 0..lmin-1 at blocks outside both readers' block sets."""
    free, nfree, nlive, ndead = _orig65(words, edges, k, an, V, own, lmin, live)
    if V is None:
        return free, nfree, nlive, ndead
    dead = dead_mask(words, edges, core.BROT)
    nb = words.shape[1] // J
    Vs = vsets(words, [k - 1, k + 1])
    qidx = {w: q for q, w in enumerate(an.out_wires)}
    ginv = dc.ginv_masks(an)
    L, nr = an.L, an.nr
    dl = []
    for l in range(L):
        v = 0
        for m in range(1, L):
            if (m & l) == l:
                v ^= an.mono_coords[m]
        dl.append(v)
    last = an.H - 1
    cols = []
    if "cur" in own:
        cols += [(J * k + j, j, k + 1) for j in range(J)]
    if "prev" in own and k > 0:
        cols += [(J * (k - 1) + j, j, k - 1) for j in range(J)]
    seen = set()
    add = cells = barred = 0
    for col, j, other in cols:
        Vo = Vs.get(other) if 1 <= other < nb else 0
        for a in np.nonzero(dead[:, col])[0].tolist():
            h, l = a >> 4, a & 15
            if l >= lmin or h >= last or (a, j) in seen:
                continue
            if (V >> h) & 1 or Vo is None or (Vo >> h) & 1:
                barred += 1
                continue
            seen.add((a, j))
            vec = np.zeros(nr, np.uint32)
            for M in range(1, nr + 1):
                if (M & h) == h:
                    vec[M - 1] = dl[l]
            if not vec.any():
                continue
            cells += 1
            for b in range(1, core.BROT - 1):
                q = qidx.get((j, b))
                if q is None or q not in ginv:
                    continue
                o, gm = ginv[q]
                free.setdefault(o, []).append(dict(vec=vec, g=gm, a=a, col=col, bit=b, h=h, l=l, cost=0.0))
                add += 1
    for o in free:
        free[o].sort(key=lambda f: f["cost"])
    print(f"r78: low points 0-{lmin - 1}: {cells} dead addresses added ({add} free bits) at blocks outside both readers' "
          f"block sets; {barred} left out (inside one of them)", flush=True)
    return free, nfree + add, nlive, ndead


_orig65 = None
_origw = None


def free_bits_extw(words, edges, k, an, V, own, lmin, near=None):
    """builder C3's relation-keeping unit builder (chem_c3_wdcrel.free_bits, called unchanged, twice) + its own units
    at low points 0..lmin-1 whose block stands outside the block set of every loader that reads a flipped word."""
    free, nfree, ndead = _origw(words, edges, k, an, V, own, lmin, near)
    if V is None or lmin <= 0:
        return free, nfree, ndead
    nb = words.shape[1] // J
    Vs = vsets(words, [k - 1, k + 1])
    f0, _n0, _d0 = _origw(words, edges, k, an, V, own, 0, near)
    add = cells = barred = 0
    seen = set()
    for o, lst in f0.items():
        for f in lst:
            if f["l"] >= lmin:
                continue
            h = f["h"]
            ok = not ((V >> h) & 1)
            for col, _b in f.get("flips", [(f["col"], f["bit"])]):
                bat = col // J
                other = k + 1 if bat == k else k - 1
                if 1 <= other < nb:
                    Vo = Vs.get(other)
                    if Vo is None or (Vo >> h) & 1:
                        ok = False
            if not ok:
                barred += 1
                continue
            free.setdefault(o, []).append(f)
            add += 1
            seen.add((f["a"], f["col"]))
    cells = len(seen)
    STAT.update(cells=cells, add=add, barred=barred)
    print(f"r78: low points 0-{lmin - 1}: {add} units added at {cells} dead words (blocks outside every reader's block "
          f"set); {barred} units left out (a block inside one of them)", flush=True)
    return free, nfree + add, ndead


def install_emul(seed):
    """COUNT-ONLY EMULATION OF RULE (c) (the placer keeps rows that hold dead words off the low points 0-2): the dead
    mask every tool sees is the table's own with ROWS exchanged inside their group (127 rows a group from address 76):
    every row at a low point 0-2 (or in the last block) that holds a dead word changes place with a row of the same
    group at a usable address that holds none. The words are NOT moved (a moved row would need its wire bits refitted
    by the rounding), so the emitted loader is a count for a table whose dead rows sit there, not a loader of this
    table. Patches dead_mask in the modules the owner tools read it from."""
    import r62_chem_cache_dead as r62
    orig = r62.dead_mask
    cache = {}

    def perm_for(d):
        A = d.shape[0]
        key = (A, d.shape[1])
        if key in cache:
            return cache[key]
        rng = np.random.default_rng(seed)
        a = np.arange(A)
        last = (A + 15) // 16 - 1
        bad = ((a & 15) < 3) | ((a >> 4) >= last)
        cnt = d.sum(axis=1)
        perm = a.copy()
        moved = short = 0
        for g in range(9):
            rows = 76 + 127 * g + np.arange(127)
            rows = rows[rows < A]
            src = [int(r) for r in rows if bad[r] and cnt[r] > 0]
            dst = [int(r) for r in rows if (not bad[r]) and cnt[r] == 0]
            rng.shuffle(dst)
            src.sort(key=lambda r: -cnt[r])
            for r, s in zip(src, dst):
                perm[r], perm[s] = s, r
                moved += 1
            short += max(0, len(src) - len(dst))
        cache[key] = perm
        d2 = d[perm]
        us = ~bad
        print(f"r78 EMULATION (seed {seed}): {moved} rows exchanged inside their groups ({short} left for want of an empty "
              f"row); dead words at usable addresses {int(d[us].sum())} -> {int(d2[us].sum())} of {int(d.sum())}", flush=True)
        return perm

    def dead_mask(w, edges, brot, thresh=0):
        d = orig(w, edges, brot, thresh)
        if d.shape[0] < 1219:
            return d
        dd = d[:1219]
        out = d.copy()
        out[:1219] = dd[perm_for(dd)]
        return out

    r62.dead_mask = dead_mask
    for name in ("r63_chem_wired_dc", "r65_chem_wired221_dc", "chem_c3_wdcrel", "chem_ld2_e22dcr"):
        mod = sys.modules.get(name)
        if mod is not None and hasattr(mod, "dead_mask"):
            mod.dead_mask = dead_mask
    return dead_mask


def count(table, owners):
    z = np.load(table)
    words = np.asarray(z["angle_words"], np.int64)
    edges = np.asarray(z["tree_edges"], np.int64)
    core.J, core.BROT = 5, 7
    A = words.shape[0]
    dead = dead_mask(words, edges, 7)
    nb = words.shape[1] // J
    Vs = vsets(words, range(1, nb))
    print("block sets V (blocks in V of 77): " + ", ".join(f"k{k}: {bin(v).count('1') if v is not None else None}" for k, v in Vs.items()))
    last = (A + 15) // 16 - 1
    a = np.arange(A)
    h, l = a >> 4, a & 15
    tot_low = tot_ok = tot_pair = 0
    for k in owners:
        row = []
        for b, other in ((k - 1, k - 1), (k, k + 1)):
            Vk = Vs.get(k)
            Vo = Vs.get(other) if 1 <= other < nb else 0
            for j in range(J):
                m = dead[:, J * b + j] & (l < 3) & (h < last)
                inV = np.array([((Vk >> int(x)) & 1) if Vk is not None else 1 for x in h]).astype(bool)
                inO = np.array([((Vo >> int(x)) & 1) if Vo is not None else 1 for x in h]).astype(bool)
                ok = m & ~inV & ~inO
                row.append((int(m.sum()), int(ok.sum())))
        nlow = sum(x[0] for x in row)
        nok = sum(x[1] for x in row)
        lastb = int((dead[:, J * (k - 1):J * k + J] & (h >= last)[:, None]).sum())
        tot_low += nlow
        tot_ok += nok
        print(f"owner k{k}: dead words at low points 0-2 under the last block {nlow} (last block {lastb}); at a block "
              f"outside both readers' block sets {nok}; per word slot (both batches) "
              f"{[row[j][1] + row[J + j][1] for j in range(J)]} of {[row[j][0] + row[J + j][0] for j in range(J)]}")
    print(f"ALL OWNERS: {tot_ok} of {tot_low} dead words at low points 0-2 stand at a block outside both block sets")


def main():
    global _orig65, _real_vset
    mode = sys.argv[1]
    if mode == "count":
        owners = [4, 8, 10, 12, 14]
        if "--owners" in sys.argv:
            owners = [int(x) for x in sys.argv[sys.argv.index("--owners") + 1].split(",")]
        count(sys.argv[2], owners)
        return
    assert mode in ("run", "run65", "vrun65", "vcount", "vrun")
    args = sys.argv[2:]
    base = "--base" in args
    args = [x for x in args if x != "--base"]
    module = "r65_chem_wired221_dc"
    if mode == "vrun":          # FOR BUILDER C3: any module of its chain: `vrun <module> <its arguments>` = chem_c3_vrun.py's call
        module, args = args[0], args[1:]
    if mode in ("vrun65", "vcount", "vrun"):
        # THE CHAIN'S OWN CONDITIONS: builder C3's wrapper tools/chem_c3_vrun.py (run, not copied): the relation-aware
        # block set in place of c4.vset, units that keep the other reader's relations; my bits go UNDER its wrappers.
        import runpy
        import chem_c3_wdcrel as w
        w.install_vset_rel()
        _real_vset = c4.vset
        if mode == "vcount":
            owners = [4, 8, 10, 12, 14]
            count(args[0], owners)
            return
        import r65_chem_wired221_dc as r65
        _orig65 = r65.free_bits
        if "--emul" in args:                      # count-only emulation of rule (c); the tools stay as they are
            i = args.index("--emul")
            es = int(args[i + 1])
            del args[i:i + 2]
            install_emul(es)
            w.m.dead_mask = sys.modules["r62_chem_cache_dead"].dead_mask if hasattr(w, "m") else None
            base = True
        if not base:
            global _origw
            r65.free_bits = free_bits_ext65       # a plain loader keeps R65's builder
            _origw = w.free_bits                  # a loader with relations (every wired one on a G table) uses C3's
            w.free_bits = free_bits_extw
            vset_keep._c3 = True
            c4.vset = vset_keep
        here = os.path.dirname(os.path.abspath(__file__))
        sys.argv = [os.path.join(here, "chem_c3_vrun.py"), module] + args
        runpy.run_path(os.path.join(here, "chem_c3_vrun.py"), run_name="__main__")
        return
    sys.argv = [sys.argv[0]] + args
    if mode == "run":
        if not base:
            dc.free_bits = free_bits_ext
            c4.vset = vset_keep
        dc.main()
        return
    # run65: R65-CHEM-WIRED221's owner tool (the route of the entry: dead bits + exact-priced live units)
    import r65_chem_wired221_dc as r65
    _orig65 = r65.free_bits
    if not base:
        r65.free_bits = free_bits_ext65
        c4.vset = vset_keep
    r65.main()


if __name__ == "__main__":
    main()
