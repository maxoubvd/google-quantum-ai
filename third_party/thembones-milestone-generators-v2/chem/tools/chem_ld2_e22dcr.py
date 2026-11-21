#!/usr/bin/env python3
# C2 (22nd holder): PINNED COPY of tools/r63_chem_wired_dc.py (R63-CHEM-WIRED; read-only there), made by
# tools/chem_ld2_e22mkdcr.py; one change: ginv_masks skips wires bound by a relation instead of asserting.
"""r63_chem_wired_dc.py: designer R63-CHEM-WIRED. DEAD WORDS AS A SECOND SUPPLY OF KILL UNKNOWNS in the placed-kill
group engine of a reduced (wired) loader.

The group engine pays every forced relation (column kill, row kill) with re-fan-out bits of EARLIER products: 75
equations a column kill a output, about 14 a row kill; the first group has no pool, the second 39 a output, and the
plan 39 | 38 | 42 | 36 is the floor of that supply (R59-CHEM-GROUPS). R62-CHEM-CACHE measured that 2.6% of the words
of a rounded table are never read by the decoder (their parent amplitude is exactly zero): bits 1-5 of such a word
are FREE TABLE BITS at zero mHa^2 (tools/r62_chem_cache_dead.py, read-only import). A free bit of wire q at address
(block h, low point l) adds d_l to the rows M >= h (bitwise) of q's column of the mixed tensor: it is one more
unknown in the kill systems of q's basis output, independent of the pool. This tool adds those unknowns
(class DcKill, a subclass of the builder's KillDecomp; nothing of the builder's is edited), runs the engine,
rewrites the table at the dead words it used, re-pads, re-analyses, checks the tensor model against the new table
and emits with the builder's emit_c4 (every wire on every address).

    CHEM_BROT=7 CHEM_HB=77 python -I tools/r63_chem_wired_dc.py <submission.npz> --k 14 --plan "6,6:1:0:0,6:1:3:1,7:2:0:0" \
        [--own cur,prev] [--nofree] [--seeds 0-7] [--pad] [--emit DIR] [--lmin 3]

--own: which batches' dead words this loader owns (cur = batch k, prev = batch k - 1); --nofree: the control (the
builder's engine on the same plan); --lmin: lowest low point used (3 keeps V and kappa of the C-4 closure fixed).
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_ld2_wwidth  # noqa: E402,F401
import chem_ld2_e17kill as ek  # noqa: E402
from r62_chem_cache_dead import dead_mask  # noqa: E402

bits = core.bits


def wire_tensor(an, X0):
    """per wire: the XOR of its basis columns (the wire's mixed tensor)."""
    out = []
    for q, g in enumerate(an.G):
        v = np.zeros(an.nr, np.uint32)
        for o in bits(g):
            v ^= X0[:, o]
        out.append(v)
    return out


def ginv_masks(an):
    """for every generic wire q: (pivot basis output, mask over basis outputs that change when wire q alone changes)."""
    gen = [(g.bit_length() - 1, q) for q, g in enumerate(an.G) if g]
    gen.sort()
    out = {}
    for pos, (b, q) in enumerate(gen):
        m = 1 << b
        for b2, q2 in gen[pos + 1:]:
            if (an.G[q2] & m).bit_count() & 1:
                m ^= 1 << b2
        # check: G[q'] . m = [q' == q]; C2's copy (22nd holder): a wire bound by a relation (bit-1 relations of
        # WR7's P family: the wire is a sum of others) has no such mask: it gets no free bits instead of an assert
        if all(((an.G[q2] & m).bit_count() & 1) == (1 if q2 == q else 0) for b2, q2 in gen):
            out[q] = (b, m)
    return out


def rebase(an, order):
    """the output basis of analyze() rebuilt with the wires taken in another order (same procedure, same pivots
    rule): the groups of the engine follow the basis order, so this chooses which wires sit in which group."""
    nr, nl = an.nr, an.nl
    T = wire_tensor(an, an.X0)
    newpos = sorted(range(nl), key=lambda t: (an.level[t], t))
    pos_of = {t: i for i, t in enumerate(newpos)}
    spo = core.Span()
    G = [0] * len(an.G)
    Xcols, levels = [], []
    for q in order:
        vec = 0
        for t in range(nl):
            col = core.to_int((T[q] >> np.uint32(t)) & 1)
            if col:
                vec |= col << (pos_of[t] * nr)
        r, c = spo.reduce(vec)
        if r:
            b = len(Xcols)
            spo.piv[r.bit_length() - 1] = (r, 1 << b)
            G[q] = c | (1 << b)
            col = np.zeros(nr, np.uint32)
            for t in range(nl):
                cbits = (r >> (pos_of[t] * nr)) & ((1 << nr) - 1)
                if cbits:
                    col |= core.int_to_bits(cbits, nr).astype(np.uint32) << np.uint32(t)
            Xcols.append(col)
            levels.append(an.level[newpos[(r.bit_length() - 1) // nr]])
        else:
            G[q] = c
    T2 = None
    an.G, an.X0, an.out_level = G, np.stack(Xcols, axis=1), levels
    T2 = wire_tensor(an, an.X0)
    assert all(np.array_equal(x, y) for x, y in zip(T, T2)), "rebase changed a wire"


def dead_per_word(words, edges, k, own, lmin):
    dead = dead_mask(words, edges, core.BROT)
    J = core.J
    low = np.arange(words.shape[0]) & 15
    out = []
    for j in range(J):
        m = np.zeros(words.shape[0], bool)
        if "cur" in own:
            m |= dead[:, J * k + j]
        if "prev" in own and k > 0:
            m |= dead[:, J * (k - 1) + j]
        out.append(int((m & (low >= lmin)).sum()))
    return out


def near_cost(table, words, edges, brot):
    """proxy price (mHa^2 of error_2^2, no refit) of flipping bit b of every two-body word [E]: the decoded vector
    moves by 2 |sin(delta / 2)| x (the parent amplitude when the edge is applied), and the row's weight in the
    tensor is 8 sum_c W_c^2 ||X_rc||^2. Returns amp (rows, cols) and rowS (rows), one-body rows at infinity."""
    A, ncol = words.shape
    N = ncol + 1
    th = np.pi * (words % (1 << brot)) / float(1 << (brot - 1))
    cur = np.zeros((A, N))
    cur[:, 0] = 1.0
    amp = np.zeros((A, ncol))
    for c, (p, ch) in enumerate(edges):
        par = cur[:, p].copy()
        amp[:, c] = np.abs(par)
        cur[:, ch] = np.sin(th[:, c]) * par
        cur[:, p] = np.cos(th[:, c]) * par
    W = np.asarray(table["W"], float)                    # (R, B, C)
    R, B, C = W.shape
    n1 = A - R * B
    u = cur[n1:].reshape(R, B, N)                        # decoded vectors of the two-body rows
    rowS = np.full(A, np.inf)
    for r in range(R):
        G2 = (u[r] @ u[r].T) ** 2                        # (B, B): (u_b . u_b')^2
        xn = np.einsum("bc,bd,dc->c", W[r], G2, W[r])    # ||X_rc||_F^2 per channel
        rowS[n1 + r * B:n1 + (r + 1) * B] = 8.0 * (W[r] ** 2 @ xn)
    return amp, rowS


def free_bits(words, edges, k, an, V, own, lmin, near=None):
    """dict basis output -> list of free bits (vec uint32[nr], gmask, row, column, bit)."""
    A0, ncol = words.shape
    J = core.J
    nr = an.nr
    dead = dead_mask(words, edges, core.BROT)
    qidx = {w: q for q, w in enumerate(an.out_wires)}
    ginv = ginv_masks(an)
    L = an.L
    dl = []
    for l in range(L):
        v = 0
        for m in range(1, L):
            if (m & l) == l:
                v ^= an.mono_coords[m]
        dl.append(v)
    last = an.H - 1
    free = {}
    seen = set()
    cols = []
    if "cur" in own:
        cols += [(J * k + j, j) for j in range(J)]
    if "prev" in own and k > 0:
        cols += [(J * (k - 1) + j, j) for j in range(J)]
    nfree = 0
    for col, j in cols:
        for a in np.nonzero(dead[:, col])[0].tolist():
            h, l = a >> 4, a & 15
            if l < lmin or h >= last or (a, j) in seen:
                continue
            seen.add((a, j))
            rows = [M for M in range(1, nr + 1) if (M & h) == h]
            vec = np.zeros(nr, np.uint32)
            for M in rows:
                vec[M - 1] = dl[l]
            if (V is not None) and ((V >> h) & 1) and last <= nr:
                vec[last - 1] ^= np.uint32(dl[l])            # the padding row (last block, l) follows the V-sum
            if not vec.any():
                continue
            for b in range(1, core.BROT - 1):
                q = qidx.get((j, b))
                if q is None or q not in ginv:
                    continue
                o, gm = ginv[q]
                free.setdefault(o, []).append(dict(vec=vec, g=gm, a=a, col=col, bit=b, h=h, l=l))
                nfree += 1
    nnear = 0
    if near is not None:
        tau, maxbit, amp, rowS = near
        tb = core.BROT - 1
        for col, j in cols:
            for a in range(words.shape[0]):
                h, l = a >> 4, a & 15
                if dead[a, col] or l < lmin or h >= last or (a, j) in seen or not np.isfinite(rowS[a]):
                    continue
                base = 1e6 * rowS[a] * amp[a, col] ** 2
                if base * (2 * np.sin(np.pi * 2 / (1 << core.BROT))) ** 2 > tau:
                    continue
                vec = None
                for b in range(1, min(maxbit, tb - 1) + 1):
                    du = 2 * abs(np.sin(np.pi * ((1 << b) + ((1 << tb) if b == tb - 1 else 0)) / (1 << core.BROT)))
                    cost = base * du * du
                    if cost > tau:
                        continue
                    q = qidx.get((j, b))
                    if q is None or q not in ginv:
                        continue
                    if vec is None:
                        vec = np.zeros(nr, np.uint32)
                        for M in range(1, nr + 1):
                            if (M & h) == h:
                                vec[M - 1] = dl[l]
                        if (V is not None) and ((V >> h) & 1) and last <= nr:
                            vec[last - 1] ^= np.uint32(dl[l])
                        if not vec.any():
                            break
                    o, gm = ginv[q]
                    free.setdefault(o, []).append(dict(vec=vec, g=gm, a=a, col=col, bit=b, h=h, l=l, cost=cost))
                    nnear += 1
                seen.add((a, j))
        for o in free:
            free[o].sort(key=lambda f: f.get("cost", 0.0))
    return free, nfree + nnear, int(dead[:, [c for c, _ in cols]].sum())


class DcMixin:
    """kill systems that also hold the free table bits of each output (dead words); first in the MRO of a KillDecomp."""

    def _dc_init(self, free):
        self.free = free or {}
        self.used = {}                      # basis output -> list of free bits flipped for it
        self.gcount = []                    # (outputs, products, last kill stat) of every solved group
        self.laststat = {}

    def solve_group(self, outs, coords, dl=0, d=0, tag=""):
        n0 = len(self.prods)
        made = super().solve_group(outs, coords, dl, d, tag)
        self.gcount.append((len(outs), len(self.prods) - n0, self.laststat.get(tuple(outs))))
        return made

    def apply_kills(self, outs, coords, dl, d, shift=0):
        shift = (shift + self.shift0) % 6
        an, nr = self.an, self.nr
        for o in outs:
            self.used.pop(o, None)
        order = core.kill_order(coords, an)
        killed = order[shift:shift + dl] if shift + dl <= len(order) else order[:dl]
        live = coords
        for t in killed:
            live &= ~(1 << t)
        usable = [i for i, p in enumerate(self.prods) if not ((p[0] | p[2]) & ~coords)]
        npool = len(usable)
        if not usable and not any(self.free.get(o) for o in outs):
            return None
        P = self.PM[usable] if usable else np.zeros((0, nr), np.uint32)
        zs = []
        for _ in range(d):
            z = np.array([self.rng.random() < 0.5 for _ in range(nr)])
            if not z.any():
                z[self.rng.randrange(nr)] = True
            zs.append(z)

        def vectors(R):
            """R: (n, nr) uint32 rows -> list of ints over the kill equations."""
            n = R.shape[0]
            cols = [0] * n
            width = 0
            for t in killed:
                for i, v in enumerate(core.rows_to_ints((R >> np.uint32(t)) & 1)):
                    cols[i] |= v << width
                width += nr
            for z in zs:
                red = np.bitwise_xor.reduce(R[:, z], axis=1) & np.uint32(live)
                for i, v in enumerate(red.tolist()):
                    cols[i] |= v << width
                width += an.nl
            return cols

        sp = core.Span()
        if npool:
            for i, c in enumerate(vectors(P)):
                sp.add(c, 1 << i)
        Xw = self.X.copy()
        plan = []
        nf_used = 0
        for o in sorted(outs):
            fl = self.free.get(o, [])
            spo = sp.copy()
            if fl:
                F = np.stack([f["vec"] for f in fl])
                for i, c in enumerate(vectors(F)):
                    if c:
                        spo.add(c, 1 << (npool + i))
            tgt = vectors(Xw[:, [o]].T.copy())[0]
            r, c = spo.reduce(tgt)
            if r:
                return None
            idx = bits(c)
            pidx = [i for i in idx if i < npool]
            fidx = [i - npool for i in idx if i >= npool]
            if pidx:
                Xw[:, o] ^= np.bitwise_xor.reduce(P[pidx], axis=0)
            for i in fidx:
                for o2 in bits(fl[i]["g"]):
                    Xw[:, o2] ^= fl[i]["vec"]
            plan.append((o, pidx, [fl[i] for i in fidx]))
            nf_used += len(fidx)
        self.X = Xw
        for o, pidx, fused in plan:
            for i in pidx:
                self.prods[usable[i]][4] ^= 1 << o
            self.used[o] = fused
        self.laststat[tuple(outs)] = (dl, d, npool, nf_used)
        return live


class DcKill(DcMixin, ek.KillDecomp):
    def __init__(self, an, seed, effort, plan, avoid_pos, V=None, kap=0, tries=1, free=None):
        ek.KillDecomp.__init__(self, an, seed, effort, plan, avoid_pos, V, kap, tries)
        self._dc_init(free)


def make_tie(an, seed, effort, plan, avoid_pos, modes, V, kap, tries, free, gtrials=40):
    """the R59 TieKill (column ties on the second group, regular planes on the first) with the free columns."""
    import r59_chem_groups_tie as tie

    class DcTie(DcMixin, tie.TieKill):
        pass
    kd = DcTie(an, seed, effort, plan, avoid_pos, modes, V, kap, tries, gtrials)
    kd._dc_init(free)
    return kd


def modified_words(words, used):
    w = words.copy()
    tb = core.BROT - 1
    for o, lst in used.items():
        for f in lst:
            w[f["a"], f["col"]] ^= (1 << f["bit"])
            if f["bit"] == tb - 1:
                w[f["a"], f["col"]] ^= (1 << tb)        # bit 6 follows bit 5 (the carrier bit6 ^ bit5 stays)
    return w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--plan", required=True)
    ap.add_argument("--avoid", default="")
    ap.add_argument("--seeds", default="0-3")
    ap.add_argument("--effort", type=int, default=3)
    ap.add_argument("--pad", action="store_true")
    ap.add_argument("--tries", type=int, default=400)
    ap.add_argument("--own", default="cur,prev")
    ap.add_argument("--nofree", action="store_true")
    ap.add_argument("--lmin", type=int, default=3)
    ap.add_argument("--emit", default=None)
    ap.add_argument("--stop", type=int, default=0)
    ap.add_argument("--modes", default="", help="R59 group modes (reg,tie3,-,-): the tie engine under the free columns")
    ap.add_argument("--near", type=float, default=0.0, help="also free: bits of live two-body words whose proxy price is at most this (mHa^2 a flip)")
    ap.add_argument("--maxbit", type=int, default=5)
    ap.add_argument("--worder", default="", help="word order of the output basis: auto (poorest word first) or e.g. 2,3,4,1,0")
    a = ap.parse_args()
    core.J, core.BROT = 5, 7
    z = np.load(a.table)
    words = np.asarray(z["angle_words"], np.int64)
    edges = np.asarray(z["tree_edges"], np.int64)
    V, kap = None, 0
    wp = words
    if a.pad:
        import chem_ld2_c4 as c4
        V = c4.vset(words, a.k)
        wp = c4.pad_table(words, {a.k: V})
    an = core.analyze(wp, a.k, 4)
    if a.pad:
        kap = c4.kappa_wires(an, V)
    if a.worder:
        dpw = dead_per_word(words, edges, a.k, a.own.split(","), a.lmin)
        wo = sorted(range(core.J), key=lambda j: (dpw[j], j)) if a.worder == "auto" else [int(x) for x in a.worder.split(",")]
        order = sorted(range(len(an.out_wires)), key=lambda q: (wo.index(an.out_wires[q][0]), an.out_wires[q][1]))
        rebase(an, order)
        print(f"k{a.k}: dead addresses a word {dpw}; word order of the basis {wo}", flush=True)
    free, nfree, ndead = ({}, 0, 0)
    if not a.nofree:
        near = None
        if a.near > 0:
            amp, rowS = near_cost(z, words, edges, core.BROT)
            near = (a.near, a.maxbit, amp, rowS)
        free, nfree, ndead = free_bits(words, edges, a.k, an, V, a.own.split(","), a.lmin, near)
    per = sorted((o, len(v)) for o, v in free.items())
    print(f"k{a.k}: generic outputs {an.X0.shape[1]}, infra {an.infra}; dead words owned {ndead}; free bits {nfree}; "
          f"per output {[n for _, n in per]}", flush=True)
    plan = ek.parse_plan(a.plan)
    avoid = [int(x) for x in a.avoid.split(",") if x.strip()]
    best = None
    for s in ek.rng_list(a.seeds):
        t0 = time.time()
        if a.modes:
            kd = make_tie(an, s, max(2, a.effort), plan, avoid, a.modes.split(","), V, kap, a.tries if a.pad else 1, free)
        else:
            kd = DcKill(an, s, max(2, a.effort), plan, avoid, V, kap, a.tries if a.pad else 1, free)
        try:
            prods = kd.run()
        except (RuntimeError, AssertionError) as exc:
            print(f"seed {s}: FAILED {exc}", flush=True)
            continue
        n = len(prods)
        nflip = sum(len(v) for v in kd.used.values())
        line = " | ".join(f"{b}->{m}" + (f" [kills {st[0]},{st[1]} pool {st[2]} free {st[3]}]" if st else "")
                          for b, m, st in kd.gcount)
        print(f"seed {s}: {n} products = {an.infra + n}{' (-1 with C-4)' if a.pad else ''}; {line}; free bits "
              f"flipped {nflip} ({time.time() - t0:.1f} s)", flush=True)
        if best is None or n < best[0]:
            best = (n, s, prods, {o: list(v) for o, v in kd.used.items()})
        if a.stop and n <= a.stop:
            break
    if best is None:
        return
    n, s, prods, used = best
    print(f"BEST: {n} products, seed {s}", flush=True)
    # ---- the new table and the check of the tensor model against it
    X0 = an.X0.copy()
    for o, lst in used.items():
        for f in lst:
            for o2 in bits(f["g"]):
                X0[:, o2] ^= f["vec"]
    w2 = modified_words(words, used)
    nch = int((w2 != words).sum())
    dead = dead_mask(words, edges, core.BROT)
    live_changed = int(((w2 != words) & ~dead).sum())
    proxy = sum(f.get("cost", 0.0) for lst in used.values() for f in lst)
    if a.near > 0:
        print(f"LIVE words changed {live_changed}; proxy price of the flips {proxy:.4f} mHa^2 [E]", flush=True)
    else:
        assert live_changed == 0, "a live word was changed"
    w2p = w2
    if a.pad:
        V2 = c4.vset(w2, a.k)
        assert V2 == V, "the block set V moved"
        w2p = c4.pad_table(w2, {a.k: V})
    an2 = core.analyze(w2p, a.k, 4)
    t_new = wire_tensor(an2, an2.X0)
    an2.G, an2.X0, an2.out_level = an.G, X0, an.out_level
    t_mod = wire_tensor(an2, X0)
    assert an2.out_wires == an.out_wires
    bad = sum(1 for x, y in zip(t_new, t_mod) if not np.array_equal(x, y))
    print(f"table rewritten at {nch} dead words; wires whose tensor differs from the model: {bad} of {len(t_new)}",
          flush=True)
    if bad:
        return
    if a.pad:
        assert c4.kappa_wires(an2, V) == kap, "kappa moved"
    if a.emit:
        os.makedirs(a.emit, exist_ok=True)
        if a.pad:
            text, st = c4.emit_c4(an2, prods, V)
        else:
            text, st = core.emit(an2, prods)
        with open(os.path.join(a.emit, f"angle_k_{a.k}.kmx"), "w") as fh:
            fh.write(text)
        np.savez(os.path.join(a.emit, f"table_k{a.k}.npz"), angle_words=w2, tree_edges=edges)
        print(f"emitted {a.emit}/angle_k_{a.k}.kmx: toffoli {st['toffoli']}, qubits {st.get('qubits')} (checked on "
              f"every address by the emitter); new words in {a.emit}/table_k{a.k}.npz", flush=True)


if __name__ == "__main__":
    main()
