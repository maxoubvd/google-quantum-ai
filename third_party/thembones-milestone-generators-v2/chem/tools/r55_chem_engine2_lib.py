"""R55-CHEM-ENGINE2 library: the staircase's spill system (shared by the plain and the builders'-engine drivers).

A feature of an earlier stage is (product index, slot, H, s): H = its high coefficient (int over the nr high rows),
s = its fan-out (int over the nb outputs). A spill bit (f, j) adds H_f (x) e_j (x) s_f to the table, i.e. it puts the
low coordinate j of a later stage into the low form of feature f. All conditions below are LINEAR in the spill bits
and are solved in one system over all coordinates of the stage:
  * emptied outputs (designed rank-one relations e_j (x) e_o for every coordinate j of the stage): output o has no
    content on the stage's coordinates: nr equations per (o, j);
  * extra designed relations: a list of cells (j, o) whose XOR vanishes in every row: nr equations each. Same output
    on several coordinates = a rank-one relation (direction (x) output); different outputs on two coordinates = a
    PLANAR relation e_j1 (x) e_o1 + e_j2 (x) e_o2. They tie coordinates, hence the joint system;
  * rank kills: a combination z of the high rows vanishes on every live cell of the stage: one equation per
    (z, coordinate, live output).
Relations of a block after the solve: cells - rank = (old relations) + (rank kills); the designed ones are among
them and cost no fan-out room at the points they vanish on.
"""
import numpy as np

import chem_loaders_core as core

bits = core.bits


def par(x):
    return x.bit_count() & 1


def colint(X, o, j):
    """int over the high rows: bit M = bit j of X[M, o]."""
    return core.to_int((X[:, o] >> np.uint32(j)) & 1)


def tensor_of(prods, nr, nb):
    X = np.zeros((nr, nb), np.uint32)
    for y1, a1, y2, a2, s, _ in prods:
        row = core.int_to_bits(a1, nr).astype(np.uint32) * np.uint32(y1)
        if a2:
            row ^= core.int_to_bits(a2, nr).astype(np.uint32) * np.uint32(y2)
        for o in bits(s):
            X[:, o] ^= row
    return X


def feats_of(prods, n0=0):
    out = []
    for pi, (y1, a1, y2, a2, s, _) in enumerate(prods):
        if a1:
            out.append((n0 + pi, 0, a1, s))
        if a2:
            out.append((n0 + pi, 1, a2, s))
    return out


def solve_spill(X, feats, coords, killed, Z, extra, nr, nb):
    """Spill masks (one int over the low coordinates per feature) or None if the system has no solution.
    coords: the stage's low coordinates; killed: emptied outputs; Z: rank-kill row combinations;
    extra: list of relations, each a list of cells (coordinate, output)."""
    live_mask = ((1 << nb) - 1) & ~sum(1 << o for o in killed)
    nc, nk, nz = len(coords), len(killed), len(Z)
    off_cell = {}                      # (ci, ki) -> offset of nr equations
    off_row = {}                       # (ci, zi) -> offset of nb equations
    pos = 0
    for ci in range(nc):
        for ki in range(nk):
            off_cell[(ci, ki)] = pos
            pos += nr
        for zi in range(nz):
            off_row[(ci, zi)] = pos
            pos += nb
    off_extra = []
    for _ in extra:
        off_extra.append(pos)
        pos += nr
    zpar = [[par(z & H) for z in Z] for (_, _, H, _) in feats]
    sp = core.Span()
    nf = len(feats)
    for ci, j in enumerate(coords):
        ex_here = [(e, o2) for e, cells in enumerate(extra) for (jj, o2) in cells if jj == j]
        for f, (_, _, H, s) in enumerate(feats):
            col = 0
            for ki, o in enumerate(killed):
                if s >> o & 1:
                    col |= H << off_cell[(ci, ki)]
            sl = s & live_mask
            if sl:
                for zi in range(nz):
                    if zpar[f][zi]:
                        col |= sl << off_row[(ci, zi)]
            for e, o2 in ex_here:
                if s >> o2 & 1:
                    col ^= H << off_extra[e]
            if col:
                sp.add(col, 1 << (ci * nf + f))
    tgt = 0
    for ci, j in enumerate(coords):
        cj = [colint(X, o, j) for o in range(nb)]
        for ki, o in enumerate(killed):
            tgt |= cj[o] << off_cell[(ci, ki)]
        for zi, z in enumerate(Z):
            m = 0
            for o in range(nb):
                if live_mask >> o & 1 and par(z & cj[o]):
                    m |= 1 << o
            tgt |= m << off_row[(ci, zi)]
        for e, cells in enumerate(extra):
            for jj, o2 in cells:
                if jj == j:
                    tgt ^= cj[o2] << off_extra[e]
    r, c = sp.reduce(tgt)
    if r:
        return None, pos, len(sp)
    spill = [0] * nf
    for q in bits(c):
        ci, f = divmod(q, nf)
        spill[f] |= 1 << coords[ci]
    return spill, pos, len(sp)


def apply_spill(Xres, prods, feats, spill, nr):
    for f, (pi, slot, H, s) in enumerate(feats):
        if spill[f]:
            row = core.int_to_bits(H, nr).astype(np.uint32) * np.uint32(spill[f])
            for o in bits(s):
                Xres[:, o] ^= row
            prods[pi][0 if slot == 0 else 2] |= spill[f]


def check_spill(Xres, coords, killed, Z, extra, nr):
    m = np.uint32(sum(1 << t for t in coords))
    assert not (Xres[:, killed] & m).any(), "emptied outputs are not empty"
    for z in Z:
        zb = core.int_to_bits(z, nr).astype(bool)
        assert not (np.bitwise_xor.reduce(Xres[zb], axis=0) & m).any(), "rank kill failed"
    for cells in extra:
        acc = np.zeros(nr, np.uint32)
        for j, o2 in cells:
            acc ^= (Xres[:, o2] >> np.uint32(j)) & np.uint32(1)
        assert not acc.any(), "designed relation failed"


# ---- designed relations with no rank-one combination: a rank-distance-2 (Gabidulin) code of 3 x 3 matrices
def _f8mul(a, b):
    r = 0
    for i in range(3):
        if b >> i & 1:
            r ^= a << i
    for i in (4, 3):
        if r >> i & 1:
            r ^= 0b1011 << (i - 3)
    return r


def mrd_matrices(m):
    """m (<= 6) independent 3 x 3 binary matrices M[j][t] (j = coordinate of the block, t = output of the shared
    triple) such that every non-zero XOR of them has rank >= 2: the maps v -> a v and v -> b v^2 of F8."""
    maps = [lambda v, a=a: _f8mul(a, v) for a in (1, 2, 4)] + [lambda v, b=b: _f8mul(b, _f8mul(v, v)) for b in (1, 2, 4)]
    mats = []
    for f in maps[:m]:
        mats.append([[f(1 << j) >> t & 1 for t in range(3)] for j in range(3)])
    return mats


def mrd_relations(block, outs3, m):
    """The m designed relations of one 3-coordinate block on the three shared outputs, as lists of cells."""
    rel = []
    for M in mrd_matrices(m):
        rel.append([(block[j], outs3[t]) for j in range(3) for t in range(3) if M[j][t]])
    return rel
