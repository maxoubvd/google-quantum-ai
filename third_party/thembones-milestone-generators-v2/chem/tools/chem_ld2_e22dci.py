"""C2 (22nd holder): R63-CHEM-WIRED's NEXT STEP (a): the annealed relation space on the LAST group of a dead-word
owner (7 outputs after two column kills: 91 cells; 11 relations taken = 40 products; 13 = 39).
Runs my pinned copy tools/chem_ld2_e22dcr.py (R63's dead-word tool) with two hooks over its tie class (DcTie), ported
from my tools/chem_ld2_e19grp.py (Probe of cmd_dump, PlanKill of cmd_build; read there):
  DCI_DUMP=g.json  : the generic groups of the run (same format as e19grp dump: outs, coords, a, b, dim_w, wb, prods)
                     -> anneal / plan with tools/chem_ld2_e19grp.py on that file (group index --g 3 = the last group);
  DCI_INJECT=G=plans.json : explicit plans for generic group G (e19grp plan output), one plan per attempt of that
                     group (the closure retries cycle through them); everything else of the run unchanged (free bits,
                     table rewrite, tensor check, emit_c4).
Same command line as tools/chem_ld2_e22dcr.py (use one seed: the dump and the build must see the same groups 1-3).
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_ld2_e22dcr as dc  # noqa: E402
import chem_ld2_e19grp as eg  # noqa: E402

DUMP = os.environ.get("DCI_DUMP", "")
INJ = {}
for item in [x for x in os.environ.get("DCI_INJECT", "").split(";") if x.strip()]:   # 15:59Z: several groups
    gi, f = item.split("=", 1)
    INJ[int(gi)] = [[tuple(x) for x in pl["planes"]] for pl in json.load(open(f))[:int(os.environ.get("DCI_MAXPLANS", "400"))]]


class XMixin:
    """dump and injection hooks (first in the MRO, above DcMixin)."""

    def _x_init(self):
        self._pidx, self._cur, self._inject, self._planes = {}, None, None, None
        self._gi, self._dumped, self._wb, self._live = -1, [], None, 0
        self._inj_calls = 0

    def apply_kills(self, outs, coords, dl, d, shift=0):
        live = super().apply_kills(outs, coords, dl, d, shift)
        if live is not None:
            self._live = live
        return live

    def solve_group(self, outs, coords, dl=0, d=0, tag=""):
        self._live, self._cur, self._wb, self._inject = coords, None, None, None
        n0 = len(self.prods)
        if tag == "generic":
            self._gi += 1
            if self._gi in INJ:
                self._cur = self._gi
                self._inject = (outs, coords)
        made = super().solve_group(outs, coords, dl, d, tag)
        if tag == "generic" and self._wb is not None and DUMP:
            lc = eg.bits(self._live)
            loc = {t: i for i, t in enumerate(lc)}
            oi = {o: i for i, o in enumerate(outs)}
            pr = []
            for p in self.prods[n0:]:
                y1 = sum(1 << loc[t] for t in eg.bits(p[0]) if t in loc)
                y2 = sum(1 << loc[t] for t in eg.bits(p[2]) if t in loc) if p[2] else 0
                pr.append([y1, y2, sum(1 << oi[o] for o in eg.bits(p[4]) if o in oi)])
            self._dumped.append(dict(outs=list(outs), coords=lc, a=len(lc), b=len(outs), dim_w=len(self._wb),
                                     wb=[hex(v) for v in self._wb], prods=pr, made=made,
                                     log=self.log[-1] if getattr(self, "log", None) else "", injected=self._cur))
            json.dump(dict(groups=self._dumped), open(DUMP, "w"), indent=0)
        self._cur = None
        return made

    def _select(self, keys, wb, dim_w):
        self._wb = list(wb)
        if self._cur is None:
            return super()._select(keys, wb, dim_w)
        outs, coords = self._inject
        # 16:02Z: one call = one attempt at the injected group (also for a group without kills): its next plan
        pi = self._pidx.get(self._cur, -1) + 1
        self._pidx[self._cur] = pi
        self._planes = INJ[self._cur][pi % len(INJ[self._cur])]
        lc = eg.bits(self._live)
        A = len(lc)
        ylocal = []
        for ym in self.an.ymasks:
            m = sum(1 << ci for ci, t in enumerate(lc) if ym >> t & 1)
            if m and m not in ylocal:
                ylocal.append(m)
        chosen, F, why = [], eg.Span(), None
        for (y1, y2, s) in self._planes:
            ps = [p for p, ym in enumerate(ylocal) if not ((y1 | y2) & ~ym)]
            if not ps:
                why = f"no phase holds the plane (phases {len(ylocal)})"
                break
            key = (s, ps[0])
            for y in (y1, y2):
                flat = sum(y << (gi * A) for gi in eg.bits(s))
                if not F.add(flat):
                    why = "dependent features"
                chosen.append((key, y, flat))
        if why is None and any(F.reduce(w) for w in wb):
            why = f"W not inside the plan's span (dim W {dim_w})"
        self._inj_calls += 1
        if why is not None:
            if self._inj_calls <= 3:
                print(f"  injection refused: {why} (a {A}, b {len(outs)})", flush=True)
            return None
        if self._inj_calls == 1 or self._inj_calls % 50 == 0:
            print(f"  injection accepted: {len(self._planes)} products (call {self._inj_calls})", flush=True)
        return len(self._planes), chosen


def make_tie(an, seed, effort, plan, avoid_pos, modes, V, kap, tries, free, gtrials=40):
    import importlib   # C2_TIEMOD=1: my pinned copy tools/chem_ld2_e22tie.py (C2_TIE_EXTRA more free relations)
    tie = importlib.import_module("chem_ld2_e22tie" if os.environ.get("C2_TIEMOD") else "r59_chem_groups_tie")

    if os.environ.get("C2_LUCKY"):   # 17:58Z: R66-PLAN244's lucky rank kills on group LUCKY_G (tools/chem_ld2_e22lucky.py)
        from chem_ld2_e22lucky import LuckyMixin

        class DcTieX(XMixin, LuckyMixin, dc.DcMixin, tie.TieKill):
            pass
    else:
        class DcTieX(XMixin, dc.DcMixin, tie.TieKill):
            pass
    ntry = max(tries, max((len(v) for v in INJ.values()), default=0))
    kd = DcTieX(an, seed, effort, plan, avoid_pos, modes, V, kap, ntry, gtrials)
    kd._dc_init(free)
    kd._x_init()
    if os.environ.get("C2_LUCKY"):
        kd._lucky_init()
    return kd


dc.make_tie = make_tie

if __name__ == "__main__":
    dc.main()
