"""chem_wr7_g2rel (WR7, seat 5): a WIRED table with R65-CHEM-ALL15's fitted block constants -> a RELATION FILE that the
relation chain and the loaders read with the constants (R65's point 5).

INPUT: rounded.npz of `tools/r65_chem_all15_froundg.py --fp --fp-g` (two-body rows at addresses 76.., `wire_g`: one row
per wire of the list, one column per block 4..76, block = address >> 4) and the run's wire list (k:j:0:m,...).

A wire k:j:0:m asks  bit0(word 5k+j) XOR bit0(word 5k-5+j) = parity(address AND m) XOR g(address >> 4), g(0..3) = 0
(the one-body blocks 0..3 are not fitted by the run: the masks hold there, as today). It is written as R65's target
kind (`tools/chem_rk7_rel.py` as patched by R65: ~/qcc/ws/r65_chem_all15/patched/chem_rk7_rel.py; C2's reader takes
any non-integer target as 'derived from the table', f(low 4 bits) XOR g(block)):

    k | j:0 | G/<m'>/<hex'>      m' = m XOR 64 g(4),   bit b of hex' = g(b) XOR g(4) * (bit 2 of b)

the SAME target on every address (parity(a AND 64) = bit 6 of a = bit 2 of a >> 4), with block 4's bit of hex' = 0:
the one-body rows (addresses 0..75) read the integer base m' as in R65's relp-g files (block 4 pinned) [P].

    python chem_wr7_g2rel.py ROUNDED.npz WIRES.txt OUT_REL.txt          writes the file, checks it on ROUNDED.npz
    python chem_wr7_g2rel.py --check TABLE.npz REL.txt                  a table (angle_words, all rows) or a rounded.npz
    python chem_wr7_g2rel.py --selftest
Exit 1 if a relation is broken.
"""
import sys

import numpy as np

N_ONE = 76


def parse_wires(s):
    if s.startswith("@"):
        s = open(s[1:]).read()
    out = []
    for p in s.replace("\n", ",").split(","):
        p = p.strip()
        if p:
            k, j, b, m = (int(v) for v in p.split(":"))
            assert b == 0, "bit-0 wires only"
            out.append((k, j, m))
    return out


def par(a, m):
    a = np.asarray(a, np.int64) & int(m)
    r = np.zeros(a.shape, np.int64)
    while np.any(a):
        r ^= a & 1
        a = a >> 1
    return r


def target(t, addrs):
    """R65's patched ttable for the kinds written here: an int mask, or G/<int base>/<hex>."""
    addrs = np.asarray(addrs, np.int64)
    if isinstance(t, str) and t.startswith("G/"):
        _, base, hx = t.split("/")[:3]
        gv = int(hx, 16)
        gb = np.array([(gv >> i) & 1 for i in range(128)], np.int64)
        return par(addrs, int(base, 0)) ^ gb[addrs >> 4]
    return par(addrs, int(t))


def read_rel(path):
    rels = []
    for ln in open(path):
        ln = ln.split("#")[0].strip()
        if not ln:
            continue
        k, outs, t = (v.strip() for v in ln.split("|"))
        o = [tuple(int(x) for x in q.split(":")) for q in outs.split()]
        rels.append((int(k), o, t if t.startswith("G/") else int(t)))
    return rels


def broken(words, addrs, rels):
    w = np.mod(np.asarray(words, np.int64), 128)
    per = []
    for (k, o, t) in rels:
        acc = np.zeros(len(w), np.int64)
        for (j, b) in o:
            acc ^= ((w[:, 5 * k + j] ^ w[:, 5 * k - 5 + j]) >> b) & 1
        per.append(int(np.sum(acc != target(t, addrs))))
    return np.array(per, np.int64)


def load_rows(path):
    z = np.load(path)
    if "angle_words" in z.files:
        w = np.asarray(z["angle_words"], np.int64)
        return z, w, np.arange(len(w))
    w = np.asarray(z["words"], np.int64)
    return z, w, N_ONE + np.arange(len(w))


def g2rel(wires, wire_g):
    """wires [(k, j, m)], wire_g (n_wires, 73) for blocks 4..76 -> relation lines."""
    wire_g = np.asarray(wire_g, np.int64)
    assert wire_g.shape[0] == len(wires), (wire_g.shape, len(wires))
    lines = []
    for (k, j, m), gw in zip(wires, wire_g):
        g = np.zeros(128, np.int64)
        g[4:4 + len(gw)] = gw
        g4 = int(g[4])
        mp = int(m) ^ (64 * g4)
        hx = 0
        for b in range(128):
            if g[b] ^ (g4 * ((b >> 2) & 1)):
                hx |= 1 << b
        assert not (hx >> 4) & 1
        lines.append(f"{k} | {j}:0 | G/{mp}/{hx:x}")
    return lines


def selftest():
    rng = np.random.default_rng(5)
    wires = [(k, j, int(rng.integers(1, 128))) for k in (2, 3, 13) for j in range(5)]
    A = N_ONE + 1143
    a = np.arange(A)
    G = rng.integers(0, 2, (len(wires), 73))
    w = rng.integers(0, 128, (A, 75))
    for i, (k, j, m) in enumerate(wires):         # impose the wire with its constants on every row
        g = np.zeros(128, np.int64)
        g[4:77] = G[i]
        want = par(a, m) ^ g[a >> 4]
        cur = (w[:, 5 * k + j] ^ w[:, 5 * k - 5 + j]) & 1
        w[:, 5 * k + j] ^= cur ^ want
    rels = [(int(s.split("|")[0]), [(int(s.split("|")[1].split(":")[0]), 0)], s.split("|")[2].strip()) for s in g2rel(wires, G)]
    assert broken(w, a, rels).sum() == 0
    plain = [(k, [(j, 0)], m) for (k, j, m) in wires]
    nb = int(broken(w[N_ONE:], a[N_ONE:], plain).sum())
    one = [(k, [(j, 0)], m ^ 0) for (k, j, m) in wires]
    assert broken(w[:64], a[:64], one).sum() == 0          # blocks 0..3: the masks hold
    print(f"selftest: {len(wires)} wires with random constants on blocks 4..76: relations broken 0 with the G/ file on all "
          f"{A} rows; with the masks alone {nb} on the two-body rows; blocks 0..3 hold the masks; block 4's hex bit 0 in every line")


def main():
    if sys.argv[1] == "--selftest":
        selftest()
        return 0
    if sys.argv[1] == "--check":
        z, w, a = load_rows(sys.argv[2])
        rels = read_rel(sys.argv[3])
        per = broken(w, a, rels)
        one = int(broken(w[a < N_ONE], a[a < N_ONE], rels).sum()) if np.any(a < N_ONE) else 0
        print(f"[g2rel check] {sys.argv[2]}: {len(w)} rows, {len(rels)} relations: relations broken {int(per.sum())} "
              f"(one-body rows {one})")
        return 0 if per.sum() == 0 else 1
    rnd, wl, out = sys.argv[1], sys.argv[2], sys.argv[3]
    z, w, a = load_rows(rnd)
    wires = parse_wires("@" + wl if not wl.startswith("@") else wl)
    wg = np.asarray(z["wire_g"], np.int64)
    lines = g2rel(wires, wg)
    hdr = (f"# chem_wr7_g2rel: {rnd} ({wl}): {len(lines)} wires with R65's fitted block constants as G/<m'>/<hex'> "
           f"(constants set to 1: {int(wg.sum())} of {wg.size})\n")
    open(out, "w").write(hdr + "\n".join(lines) + "\n")
    rels = read_rel(out)
    per = broken(w, a, rels)
    plain = [(k, [(j, 0)], m) for (k, j, m) in wires]
    nplain = int(broken(w, a, plain).sum())
    print(f"[g2rel] wrote {out}: {len(lines)} relations; constants set to 1: {int(wg.sum())} of {wg.size}; on {rnd}: relations "
          f"broken {int(per.sum())} with the constants, {nplain} with the masks alone")
    return 0 if per.sum() == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
