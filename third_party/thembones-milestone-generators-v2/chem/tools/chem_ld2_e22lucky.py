"""C2 (22nd holder): R66-CHEM-PLAN244's LUCKY RANK KILLS ('For C2') in the dead-word group engine.
A mixin over my pinned copy of R63's tool (tools/chem_ld2_e22dcr.py DcMixin): on ONE generic group (env LUCKY_G,
0-based), after the kill systems of its outputs are solved, the null combinations of each output's system (pool
re-fans and free table bits that leave every kill equation unchanged; only free bits whose basis-output mask meets
the group in that output alone) are walked at random (one direction flipped a step) until the group's flattening
(high rows x live cells) loses LUCKY_N (default 1) more ranks, at most LUCKY_S samples (default 200,000). No
equation is added: every kill still holds; the free bits flipped go into `used` (the table rewrite) and the pool
re-fans into the products' fan-outs, exactly as the kill solve records its own.
Use: tools/chem_ld2_e22dci.py with C2_LUCKY=1 (its make_tie puts this mixin first).
"""
import os

import numpy as np

import chem_loaders_core as core
from chem_ld2_e22dcr import DcMixin, bits

LUCKY_G = int(os.environ.get("LUCKY_G", "-1"))
LUCKY_N = int(os.environ.get("LUCKY_N", "1"))
LUCKY_S = int(os.environ.get("LUCKY_S", "200000"))


def rank_rows(rows):
    piv = {}
    for v in rows:
        while v:
            b = v.bit_length() - 1
            e = piv.get(b)
            if e is None:
                piv[b] = v
                break
            v ^= e
    return len(piv)


class LuckyMixin:
    def _lucky_init(self):
        self._lgi = -1
        self.lucky_log = []

    def solve_group(self, outs, coords, dl=0, d=0, tag=""):
        if tag == "generic":
            self._lgi += 1
        return super().solve_group(outs, coords, dl, d, tag)

    def apply_kills(self, outs, coords, dl, d, shift=0):
        if self._lgi != LUCKY_G:
            return super().apply_kills(outs, coords, dl, d, shift)
        # ---- a copy of DcMixin.apply_kills (tools/chem_ld2_e22dcr.py) that also keeps each output's null space
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
        # 18:02Z: the row-kill vectors and the walk from their OWN seeded generators (not the engine's stream): every
        # attempt of this group, in a dump run and in an injection run alike, sees the same kills and the same walk
        import random as _random
        zr = _random.Random(int(os.environ.get("LUCKY_SEED", "1")) * 1000 + self._lgi)
        self._wrng = _random.Random(int(os.environ.get("LUCKY_SEED", "1")) * 1000 + 500 + self._lgi)
        zs = []
        for _ in range(d):
            z = np.array([zr.random() < 0.5 for _ in range(nr)])
            if not z.any():
                z[zr.randrange(nr)] = True
            zs.append(z)

        def vectors(R):
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
        null_pool = []
        if npool:
            for i, c in enumerate(vectors(P)):
                if not sp.add(c, 1 << i):
                    r, cc = sp.reduce(c, 1 << i)
                    null_pool.append(cc)
        Xw = self.X.copy()
        plan, nulls = [], {}
        nf_used = 0
        gmask_outs = sum(1 << o for o in outs)
        for o in sorted(outs):
            fl = self.free.get(o, [])
            spo = sp.copy()
            nl_o = list(null_pool)
            if fl:
                F = np.stack([f["vec"] for f in fl])
                for i, c in enumerate(vectors(F)):
                    if not c:
                        if (fl[i]["g"] & gmask_outs) == (1 << o):
                            nl_o.append(1 << (npool + i))
                        continue
                    if not spo.add(c, 1 << (npool + i)):
                        r, cc = spo.reduce(c, 1 << (npool + i))
                        if not any((fl[j - npool]["g"] & gmask_outs) != (1 << o) for j in bits(cc) if j >= npool):
                            nl_o.append(cc)
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
            plan.append([o, pidx, [fl[i] for i in fidx], fl])
            nulls[o] = nl_o
            nf_used += len(fidx)
        # ---- the walk: directions = (output, null combination); effect on the group's live cells
        lc = bits(live)
        cell = {(o, t): i for i, (o, t) in enumerate((o, t) for o in sorted(outs) for t in lc)}

        def eff_of(o, comb, fl):
            col = np.zeros(nr, np.uint32)
            for i in bits(comb):
                if i < npool:
                    col ^= P[i]
                else:
                    col ^= fl[i - npool]["vec"]
            e = [0] * nr
            for M in range(nr):
                x = int(col[M]) & live
                if x:
                    v = 0
                    for t in bits(x):
                        v |= 1 << cell[(o, t)]
                    e[M] = v
            return e

        dirs = []
        for o, pidx, fused, fl in plan:
            for comb in nulls[o]:
                e = eff_of(o, comb, fl)
                if any(e):
                    dirs.append((o, comb, e))
        rows = [0] * nr
        for M in range(nr):
            v = 0
            for (o, t), i in cell.items():
                if int(Xw[M, o]) >> t & 1:
                    v |= 1 << i
            rows[M] = v
        r0 = rank_rows(rows)
        found = None
        if dirs:
            state = [0] * len(dirs)
            best = r0
            for it in range(LUCKY_S):
                q = self._wrng.randrange(len(dirs))
                state[q] ^= 1
                rows = [x ^ y for x, y in zip(rows, dirs[q][2])]
                rk = rank_rows(rows)
                best = min(best, rk)
                if rk <= r0 - LUCKY_N:
                    found = (it + 1, rk, [dirs[j] for j in range(len(dirs)) if state[j]])
                    break
            self.lucky_log.append((len(dirs), r0, best, found[0] if found else LUCKY_S))
            print(f"   lucky: group {self._lgi}: {len(dirs)} directions over {len(cell)} cells, rank {r0} -> best {best}"
                  f" in {found[0] if found else LUCKY_S} samples; {'FOUND' if found else 'none'}", flush=True)
        else:
            print(f"   lucky: group {self._lgi}: no direction (rank {r0})", flush=True)
        if found:
            # apply the walk's combination: pool re-fans and free bits on top of the particular solutions
            per_o = {}
            for o, comb, _ in found[2]:
                per_o[o] = per_o.get(o, 0) ^ comb
            for item in plan:
                o, pidx, fused, fl = item
                comb = per_o.get(o, 0)
                if not comb:
                    continue
                ps = set(pidx)
                fs = {id(f): f for f in fused}
                for i in bits(comb):
                    if i < npool:
                        ps ^= {i}
                        Xw[:, o] ^= P[i]
                    else:
                        f = fl[i - npool]
                        if id(f) in fs:
                            fs.pop(id(f))
                        else:
                            fs[id(f)] = f
                        for o2 in bits(f["g"]):
                            Xw[:, o2] ^= f["vec"]
                item[1] = sorted(ps)
                item[2] = list(fs.values())
        self.X = Xw
        for o, pidx, fused, fl in plan:
            for i in pidx:
                self.prods[usable[i]][4] ^= 1 << o
            self.used[o] = fused
        self.laststat[tuple(outs)] = (dl, d, npool, nf_used)
        return live
