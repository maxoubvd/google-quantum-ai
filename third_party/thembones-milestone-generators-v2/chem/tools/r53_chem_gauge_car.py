"""r53_chem_gauge_car (design seat R53-CHEM-GAUGE): CARRIERS THAT KEEP THE RELATIONS: the exact rule and its test.

C2's carriers ask, on every row and column t: bit 6 XOR bit 5 of the stored word = c_t(l) = the slot's low-address
monomial (slot = t mod 5; monomials 3, 5, 6, 9, 10 = l0l1, l0l2, l1l2, l0l3, l1l3; l = address mod 16). By the gauge's
table (r53_chem_gauge_group.py) that FIXES the row's gauge: the half turns are h = c(l) (+ a correction at words that
sit at the pole +-32), the negations n = N h, and
        stored word_t = (-1)^n_t z_t + 64 h_t,      z = the signed-range word the rounding produced.
So a relation on bits 0..5 that must hold in the entry has to be asked, in the rounding, of the bits of
        y_t = eps_t z_t,   eps_t = (-1)^n_t(l)     (16 sign patterns, one per low address; y in [-32, 31]).

  predict   TABLE.npz [--cmp C2OUT.npz]   my recursion -> stored words; with --cmp: equality with the output of C2's own
            tool (chem_ld2_e11car7.py) on the same table, word for word; pole statistics.
  impose    TABLE.npz OUT.npz --code "code:4,6,..:0=5,1=5" [--seed 11] [--naive]
            count-only: the relations imposed on the bits of y (the stored words), written back as signed-range words z
            (pole words are never touched or made); --naive imposes them on z itself (what a rounding without the rule
            does). The table then goes through C2's tool and `check`.
  check     C2OUT.npz --code ".." [--seed 11]   relations on the stored words of a carrier table; carriers exact?
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

MONO = (3, 5, 6, 9, 10)
HS = {5: (1, 2, 4, 8, 16), 4: (0b10001, 0b10010, 0b10100, 0b11000), 3: (0b00011, 0b11000, 0b01101),
      2: (0b00111, 0b11100), 1: (0b11111,)}


def gauge_lists(edges):
    """for every edge t: the edges u whose half turn negates word t (later edges of the same parent, all edges of the child)."""
    T = len(edges)
    return [[u for u in range(T) if (u > t and edges[u][0] == edges[t][0]) or edges[u][0] == edges[t][1]] for t in range(T)]


def carriers(nrow, T, slots=MONO):
    low = np.arange(nrow) % 16
    return np.stack([((low & slots[t % 5]) == slots[t % 5]).astype(np.int64) for t in range(T)], 1)


def carrier_gauge(z, edges, want, brot=7):
    """z: signed-range words (any representation mod 2^brot). -> stored words (0..2^brot-1), n, h."""
    full, half = 1 << brot, 1 << (brot - 1)
    z = np.asarray(z, np.int64) % full
    R, T = z.shape
    deps = gauge_lists(edges)
    h = np.zeros((R, T), np.int64)
    n = np.zeros((R, T), np.int64)
    st = np.zeros((R, T), np.int64)
    for t in range(T - 1, -1, -1):
        if deps[t]:
            n[:, t] = h[:, deps[t]].sum(1) & 1
        w = np.where(n[:, t] == 1, (-z[:, t]) % full, z[:, t])
        car = ((w >> (brot - 1)) ^ (w >> (brot - 2))) & 1
        h[:, t] = car ^ want[:, t]
        st[:, t] = w ^ (half * h[:, t])
    return st, n, h


def eps_rows(z_now, edges, addresses, slots=MONO, brot=7):
    """THE HOOK'S ONE CALL. z_now (rows, 75): the rows' words as they stand, in the signed range (frozen words at their
    values, free words at their rounded means); addresses: the rows' addresses (only address mod 16 is used).
    -> eps (rows, 75) in {+1, -1} and y (rows, 75) in [-32, 31]: y = the bits 0..5 of the word the entry will hold
    (stored word = y + 64 x carrier). A relation is asked of the bits of y; a word about to freeze takes the nearest
    integer z with (eps z) mod 8 in its residue set. Words at the pole (z = -32) have y = -32 whatever eps says."""
    z_now = np.asarray(z_now, np.int64)
    low = np.asarray(addresses, np.int64) % 16
    T = z_now.shape[1]
    want = np.stack([((low & slots[t % 5]) == slots[t % 5]).astype(np.int64) for t in range(T)], 1)
    st, n, h = carrier_gauge(z_now % (1 << brot), edges, want, brot)
    half = 1 << (brot - 1)
    y = ((st + half // 2) % half) - half // 2
    return 1 - 2 * n, y


def tables(code, nrow, seed, hs=HS):
    """[(k, b, mask over slots, T over the rows)] for a layout string code:<loaders>:<bit>=<r>,..[+..]"""
    addr = np.arange(nrow)
    out = []
    for seg in code[5:].split("+"):
        ld, bits = seg.split(":")
        for k in [int(v) for v in ld.split(",")]:
            for t in bits.split(","):
                b, r = int(t.split("=")[0]), int(t.split("=")[1])
                for i in range(r):
                    rng = np.random.default_rng([int(seed), k, b, i])
                    f_, g_ = rng.integers(0, 2, 16), rng.integers(0, 2, 128)
                    out.append((k, b, hs[r][i], (f_[addr & 15] ^ g_[addr >> 4]).astype(np.int64)))
    return out


def syndromes(w, rels):
    """w: words whose bits are tested. -> (rows, relations) 0 = holds."""
    out = []
    for (k, b, hm, T) in rels:
        v = T.copy()
        for j in range(5):
            if (hm >> j) & 1:
                v ^= ((w[:, 5 * k + j] >> b) & 1) ^ ((w[:, 5 * k - 5 + j] >> b) & 1)
        out.append(v)
    return np.stack(out, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("predict", "impose", "impose2", "check", "selftest"))
    ap.add_argument("table")
    ap.add_argument("out", nargs="?")
    ap.add_argument("--cmp", default="")
    ap.add_argument("--code", default="")
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--naive", action="store_true")
    ap.add_argument("--brot", type=int, default=7)
    a = ap.parse_args()
    zf = dict(np.load(a.table))
    w = np.asarray(zf["angle_words"], np.int64) % 128
    edges = np.asarray(zf["tree_edges"], np.int64)
    R, T = w.shape
    want = carriers(R, T)
    if a.mode == "check":
        car = ((w >> 6) ^ (w >> 5)) & 1
        rels = tables(a.code, R, a.seed)
        syn = syndromes(w, rels)
        lost = np.asarray(zf["r53_lost_rows"], bool) if "r53_lost_rows" in zf else np.zeros(R, bool)
        print(f"check {a.table}: carriers exact on {int((car == want).all(0).sum())}/75 columns ({int((car != want).sum())} wrong entries); "
              f"relations {len(rels)} a row: broken {int(syn.sum())} of {syn.size} (rows with any {int(syn.any(1).sum())}); on the {int((~lost).sum())} rows "
              f"where every flip was placed: broken {int(syn[~lost].sum())} of {syn[~lost].size}")
        return
    st, n, h = carrier_gauge(w, edges, want)
    sg = np.where(w >= 64, w - 128, w)
    if a.mode == "selftest":
        # eps_rows on rows at SHUFFLED addresses against C2's own move (chem_ld2_e11car7.flip) applied until the carriers are exact
        rng = np.random.default_rng(5)
        addr = rng.permutation(R)
        eps, y = eps_rows(sg, edges, addr)
        low = addr % 16
        want2 = np.stack([((low & MONO[t % 5]) == MONO[t % 5]).astype(np.int64) for t in range(T)], 1)
        wc = w.copy()
        made_by = {int(c): t for t, (p_, c) in enumerate(edges)}
        for it in range(40):
            bad = (((wc >> 6) ^ (wc >> 5)) & 1) != want2
            if not bad.any():
                break
            for g in np.flatnonzero(bad.any(0)):
                rws = ((((wc[:, g] >> 6) ^ (wc[:, g] >> 5)) & 1) != want2[:, g])
                if not rws.any():
                    continue
                x_ = int(edges[g][0])
                wc[rws, g] ^= 64
                neg = [t for t in range(g) if int(edges[t][0]) == x_]
                if x_ in made_by:
                    neg.append(made_by[x_])
                for t in neg:
                    wc[rws, t] = (-wc[rws, t]) % 128
        yc = ((wc + 32) % 64) - 32
        ok_y = int((yc == y).sum())
        nonpole = sg != -32
        ok_eps = int(((eps * sg) == yc)[nonpole].sum())
        print(f"selftest on {a.table}, rows at shuffled addresses: carriers exact after C2's moves: {bool(not bad.any())} ({it} passes); "
              f"y of eps_rows equal to the low six bits of C2's words on {ok_y} of {y.size}; eps x z = y on {ok_eps} of {int(nonpole.sum())} words off the pole")
        return
    pole = np.abs(sg) == 32
    dev = h != want
    naive_n = np.zeros_like(n)
    deps = gauge_lists(edges)
    for t in range(T):
        if deps[t]:
            naive_n[:, t] = want[:, deps[t]].sum(1) & 1
    sens = (sg % 8 != 0) & (sg % 8 != 4)
    if a.mode == "predict":
        top = ((w >> 6) ^ (w >> 5)) & 1
        print(f"table {a.table}: {R} x {T}; words outside the signed range {int(top.sum())}; poles (+-32): {int(pole.sum())} "
              f"(+32: {int((sg == 32).sum())}, -32: {int((sg == -32).sum())}), rows with any {int(pole.any(1).sum())}")
        print(f"carrier gauge by my recursion: half turns {int(h.sum())} (the carriers ask {int(want.sum())}); words negated {int(n.sum())} = {n.mean():.3f}; "
              f"half turns that differ from the carrier (pole corrections) {int(dev.sum())} on {int(dev.any(1).sum())} rows, all at poles: {bool((dev <= pole).all())}")
        diff = n != naive_n
        print(f"negations that differ from the pole-blind rule n = N c(l): {int(diff.sum())} words on {int(diff.any(1).sum())} rows; of them with bits 1..2 "
              f"that a negation changes (word not 0 or 4 mod 8): {int((diff & sens).sum())}; in batches 3..14: {int((diff & sens)[:, 15:].sum())}")
        car = ((st >> 6) ^ (st >> 5)) & 1
        print(f"my stored words: carriers exact on {int((car == want).all(0).sum())}/75 columns; bit-0 changes against the table {int(((st ^ w) & 1).sum())}")
        if a.cmp:
            c2 = np.asarray(np.load(a.cmp)["angle_words"], np.int64) % 128
            print(f"AGAINST C2's TOOL ({a.cmp}): words equal {int((c2 == st).sum())} of {st.size}; different {int((c2 != st).sum())}")
        return
    if a.mode == "impose2":
        # count-only, by INTEGER moves like a rounding, iterated to a fixed point: every relation is asked of the stored
        # words (the carrier gauge recomputed after every sweep); a word moves by +-2^b inside [-31, 31]; a pair of two
        # poles gives one pole up (it stops being a pole and the gauge is recomputed). With --naive: asked of z, no gauge.
        rels = tables(a.code, R, a.seed)
        z = sg.copy()
        moves = broken_poles = 0
        for sweep in range(40):
            if a.naive:
                y, eps = z.copy(), np.ones_like(z)
            else:
                st_, n_, h_ = carrier_gauge(z % 128, edges, want)
                eps = 1 - 2 * n_
                y = ((st_ + 32) % 64) - 32                      # bits 0..5 of the stored word as a number in [-32, 31]
            syn = syndromes(y, rels)
            if not syn.any():
                break
            done_pairs = set()
            for row, ri in zip(*np.nonzero(syn)):
                k, b, hm, Tt = rels[ri]
                j = [jj for jj in range(5) if (hm >> jj) & 1][-1]
                key = (int(row), k, j)
                if key in done_pairs:
                    continue                                   # one move per pair and sweep (the lowest bit first)
                done_pairs.add(key)
                cols_ = [c for c in (5 * k + j, 5 * k - 5 + j) if y[row, c] != -32] or [5 * k + j]
                col = cols_[0]
                yv = int(y[row, col])
                broken_poles += int(yv == -32)
                cand = [v for v in (yv + (1 << b), yv - (1 << b)) if -31 <= v <= 31]
                yn = min(cand, key=abs)
                z[row, col] = int(eps[row, col]) * yn
                moves += 1
        zn = z % 128
        if a.naive:
            left = int(syndromes(z, rels).sum())
        else:
            st_, n_, h_ = carrier_gauge(zn, edges, want)
            left = int(syndromes(((st_ + 32) % 64) - 32, rels).sum())
        assert (np.abs(z) <= 32).all()
        zf["angle_words"] = zn.astype(zf["angle_words"].dtype)
        np.savez(a.out, **zf)
        print(f"impose2 ({'on z, no gauge' if a.naive else 'on the stored words, gauge recomputed every sweep'}): {len(rels)} relations a row, sweeps {sweep + 1}, "
              f"integer moves {moves}, poles given up {broken_poles}; relations broken at the end {left} of {R * len(rels)}; words outside the signed range "
              f"{int(((((zn >> 6) ^ (zn >> 5)) & 1) & (z != -32)).sum())}; written {a.out}")
        return
    # impose
    rels = tables(a.code, R, a.seed)
    eps_neg = n == 1
    y = np.where(eps_neg, (-w) % 128, w)                  # bits 0..5 of the stored word
    work = w.copy() if a.naive else y.copy()
    flips = skipped = 0
    lost = np.zeros(R, bool)                              # rows where a flip found no word (poles): not a test of the rule
    coded = set(k for (k, b, hm, Tt) in rels)
    by_kb = {}
    for (k, b, hm, Tt) in rels:
        by_kb.setdefault((k, b), []).append((hm, Tt))
    for (k, b), lst in sorted(by_kb.items()):
        r = len(lst)
        fix = {}
        for e in sorted(range(32), key=lambda v: bin(v).count("1")):
            s = 0
            for i in range(r):
                s |= (bin(lst[i][0] & e).count("1") & 1) << i
            fix.setdefault(s, e)
        syn = np.zeros(R, np.int64)
        for i, (hm, Tt) in enumerate(lst):
            v = Tt.copy()
            for j in range(5):
                if (hm >> j) & 1:
                    v ^= ((work[:, 5 * k + j] >> b) & 1) ^ ((work[:, 5 * k - 5 + j] >> b) & 1)
            syn |= v << i
        for row in np.nonzero(syn)[0]:
            e = fix[int(syn[row])]
            for j in range(5):
                if (e >> j) & 1:
                    done = False
                    cols_ = (5 * k + j,) if (k - 1) in coded else (5 * k + j, 5 * k - 5 + j)
                    for col in cols_:                                  # the later word first; never touch or make a pole
                        nv = int(work[row, col]) ^ (1 << b)
                        if pole[row, col] or (nv & 63) == 32:
                            continue
                        work[row, col] = nv
                        flips += 1
                        done = True
                        break
                    skipped += 0 if done else 1
                    lost[row] = lost[row] or not done
    znew = work if a.naive else np.where(eps_neg, (-work) % 128, work)
    sgn = np.where(znew >= 64, znew - 128, znew)
    assert (np.abs(sgn) <= 32).all() and (((znew >> 6) ^ (znew >> 5)) & 1)[~pole].sum() == 0, "left the signed range"
    assert ((np.abs(sgn) == 32) == pole).all(), "pole set changed"
    st2, n2, h2 = carrier_gauge(znew, edges, want)
    syn2 = syndromes(st2, rels)
    zf["angle_words"] = znew.astype(zf["angle_words"].dtype)
    np.savez(a.out, **zf)
    zf["r53_lost_rows"] = lost
    np.savez(a.out, **zf)
    print(f"imposed {len(rels)} relations a row on the bits of {'z (NAIVE: the signed-range words)' if a.naive else 'y = eps z (the stored words)'}: "
          f"bit flips {flips}, not placed {skipped} (on {int(lost.sum())} rows: a pole in the way); predicted after the carrier gauge: relations broken "
          f"{int(syn2.sum())} of {syn2.size}, on the {int((~lost).sum())} rows where every flip was placed {int(syn2[~lost].sum())}; written {a.out}")


if __name__ == "__main__":
    main()
