#!/usr/bin/env python3
"""invsqrt_fw_chiprog.py -- F3 (ninth holder, 2026-11-14): straight-line AND/XOR programs (F1's minand format) emitted
on FC's kit with a temp budget, and the H1' (m 7) carry rule Chi = C >> 1 by such a program.

program  {"gates": [{"left": [i..], "right": [j..]}], "outputs": [[k..], ..]}: index < n = input, >= n = gate
         (index - n), -1 = constant 1 (outputs only). F1's tool writes factors in reduced echelon form, so each factor
         has an element the other lacks: those are the in-place carriers.
schedule(prog, n, budget)  cheapest reversible pebbling within `budget` temporary wires: compute a gate (1 AND, its
         factors' gates live), erase it by measurement (free, its factors' gates live), XOR a live gate into the
         outputs (free); a gate used only by outputs may instead go straight into each output by one CCX.
emit(c, prog, ins, outs, budget)   outs[k] ^= output k; temporaries allocated and freed; returns the AND count.
apply_h1c(H)   FC2's emit_chi (17 ANDs) -> this module's program for Chi (env CHIPROG = path, default the 8-AND
         program `chi_k8` below; budget = cap - live at the call).
  python3 invsqrt_fw_chiprog.py test [prog.json]   all 128 patterns, random output start values, budgets 3..9
"""
import heapq
import json
import os
import sys

# Chi = C >> 1 of H1' at m 7 (spec DEF7): inputs z0..z3 (top 4 folded bits of T, z0 lowest), l0 = L''_7, l1 = L''_8,
# c0; out0 = c1 = Chi & 1, out1 = c2 = Chi >> 1. Exact on the union of the PRE and POST patterns of all 1,044,433 y
# (F1's minand_sb, 8 ANDs, F3's problem files scratchpad f3h10_chidump.py; k = 7 open).
NAMES = ["z0", "z1", "z2", "z3", "l0", "l1", "c0"]
CHI_K8 = {"gates": [
    {"left": [0, 3, 5], "right": [4, 5]},
    {"left": [1, 2, 5, 7], "right": [3]},
    {"left": [1, 2, 3], "right": [5, 6, 7]},
    {"left": [0, 2, 5, 6], "right": [1, 3, 5, 7]},
    {"left": [1, 5, 7, 10], "right": [2, 3, 6, 9]},
    {"left": [1, 5, 7, 10, 11], "right": [3, 10]},
    {"left": [0, 2, 5, 7, 8, 9, 10, 11, 12], "right": [1, 3, 5, 6, 8, 12]},
    {"left": [0, 7, 8, 9, 10, 11, 12, 13], "right": [2, 6]}],
    "outputs": [[3, 8, 12, 14], [8, 12]]}


# Two stages, 7 ANDs, at most 3 temporaries (fits cap 40's 6 clean wires): c2 alone in 4 ANDs (exact minimum: k = 3
# UNSAT), then c1 from the inputs AND the finished c2 wire in 3 ANDs (k = 2 UNSAT). Valid wherever both outputs
# start at zero (FC2's forward call and the explicit undo's rebuild) and in the mirror of such a call (reverse order,
# both wires holding their values).
CHI_C2 = {"gates": [{"left": [1, 3, 6], "right": [2, 3, 6]}, {"left": [0, 3, 5, 7], "right": [4, 5, 7]},
                    {"left": [2, 6, 7], "right": [5, 6, 8]}, {"left": [3], "right": [5, 7, 8, 9]}],
          "outputs": [[3, 10]]}
CHI_C1G = {"gates": [{"left": [0, 3, 5, 6, 7], "right": [2, 3, 5, 6, 7]}, {"left": [0, 5, 8], "right": [1, 2, 3, 6, 7]},
                     {"left": [0, 2, 5, 6, 8, 9], "right": [1, 3, 7]}],
           "outputs": [[3, 7, 10]]}


# R37's spec six_a (r37_inv_40_patch.spec_of with cub "rom"; LP-feasible): the same inputs (z = top 4 folded bits of
# T', l0 = L'_7, l1 = L'_8, c0), exact on the union of its PRE / POST patterns (105): c2 alone 4 ANDs (k = 3 UNSAT),
# c1 from the inputs and c2 in 3 (k = 2 UNSAT).
CHI6_C2 = {"gates": [{"left": [0, 3, 4], "right": [2, 3, 4, 5, 6]}, {"left": [0, 2, 3, 7], "right": [1, 6]},
                     {"left": [0, 3, 6, 7, 8], "right": [1, 3, 6, 8]}, {"left": [0, 1, 2, 7, 9], "right": [3]}],
           "outputs": [[3, 10]]}
CHI6_C1G = {"gates": [{"left": [0, 3, 7], "right": [1, 2, 3, 5, 7]}, {"left": [1, 2], "right": [5, 6, 8]},
                      {"left": [1, 6, 9], "right": [2, 5, 8]}],
            "outputs": [[3, 7, 10]]}
PROGS = {"def7": (CHI_C2, CHI_C1G), "six_a": (CHI6_C2, CHI6_C1G)}
# six_a joint (c1, c2) in 6 ANDs (k = 6 SAT, 534 s; k = 5 running): one program, outputs only XORed (any start
# value), at most 5 temporaries (g5 goes straight into c1).
CHI6_K6 = {"gates": [{"left": [2, 6], "right": [4, 6]}, {"left": [1], "right": [3]},
                     {"left": [0, 5], "right": [1, 4, 5, 6, 7]}, {"left": [0, 1, 4, 5, 7, 9], "right": [2, 3, 6]},
                     {"left": [1, 2], "right": [3, 10]}, {"left": [0, 4, 5, 7, 8, 9, 10], "right": [3, 10, 11]}],
           "outputs": [[3, 8, 11, 12], [8, 11]]}


def emit_staged(c, ins7, out1, out2, budget, spec="def7"):
    """out2 = c2 (4 ANDs), then out1 = c1 from ins7 + out2 (3 ANDs); out1, out2 must be zero here."""
    p2, p1 = PROGS[spec]
    k = emit(c, p2, ins7, [out2], budget)
    k += emit(c, p1, list(ins7) + [out2], [out1], budget)
    return k


def emit_staged_rev(c, ins7, out1, out2, budget, spec="six_a"):
    """the same XOR into out1, out2 when they HOLD c1, c2 (an erase): c1 first (reads out2 = c2), then c2."""
    p2, p1 = PROGS[spec]
    k = emit(c, p1, list(ins7) + [out2], [out1], budget)
    k += emit(c, p2, ins7, [out2], budget)
    return k


# six_a's cubic table CUB5 (5 folded bits t0..t4 = Tf_5..Tf_9 -> 5 bits) staged, 11 ANDs (best of all orders by a DP
# over the exact stage minima; the joint minimum is >= 6, k = 5 UNSAT): (bit, bits it reads, program on t0..t4 + those).
CUB5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 5, 5, 6, 7, 8, 9, 10, 11, 13, 14, 16, 17, 19]
CUB_STAGES = [
    (1, [], {"gates": [{"left": [0, 1, 2, 3], "right": [4]}, {"left": [0, 1, 2, 4], "right": [3, 4, 5]},
                       {"left": [0, 2, 3], "right": [1, 2, 4, 5]}, {"left": [0, 2, 4, 5, 6, 7], "right": [1, 2, 3, 5, 7]}],
             "outputs": [[4, 5, 8]]}),
    (2, [1], {"gates": [{"left": [1, 5], "right": [2, 3, 4, 5]}, {"left": [1, 3, 5, 6], "right": [4]}],
              "outputs": [[4, 7]]}),
    (3, [1, 2], {"gates": [{"left": [2, 6], "right": [3, 4, 6]}, {"left": [2, 4, 6], "right": [3]}],
                 "outputs": [[7, 8]]}),
    (0, [1, 2, 3], {"gates": [{"left": [0, 2, 3, 4, 5], "right": [1, 2, 3, 5]}, {"left": [3, 6], "right": [4, 7, 8]}],
                    "outputs": [[3, 5, 9]]}),
    (4, [0, 1, 2, 3], {"gates": [{"left": [1, 5, 8], "right": [4, 7, 8]}], "outputs": [[4, 7, 8, 9]]}),
]


# Joint 8-AND program for CUB5 (k = 8 SAT; k = 5 UNSAT; k 6 / 7 open): with W zero, g2 lives on W3 (stored), so it
# computes in 8 ANDs at 6 temporaries (11 at 5) and its exact uncompute (emit_uncompute) costs 7.
CUB8 = {"gates": [{"left": [0, 2, 3, 4], "right": [1, 2, 3, 4]}, {"left": [2], "right": [3, 5]},
                  {"left": [3, 6], "right": [4]}, {"left": [1, 4, 6, 7], "right": [2, 4, 5]},
                  {"left": [2, 4, 7], "right": [6]}, {"left": [0, 3, 5, 6], "right": [1, 2, 3, 6]},
                  {"left": [3, 7, 8], "right": [4, 5, 7, 8, 9]}, {"left": [0, 1, 6, 7, 9, 10, 11], "right": [4, 6, 8, 9, 10, 11]}],
        "outputs": [[2, 3, 5, 6, 7, 8, 9, 10], [2, 4, 9, 10, 11, 12], [4, 5, 7, 8, 9, 11], [7], [6, 9]]}


# The same 11 ANDs in 2 temporaries (order 3, 0, 2, 1, 4; each stage's output wire is zero when it starts, so a gate
# that IS the stage's output lives there): the cheapest staged order within 2 of all 120 [M].
CUB_STAGES2 = [
    (3, [], {"gates": [{"left": [0, 2, 3, 4], "right": [1, 2, 3]}, {"left": [0, 4, 5], "right": [2, 4]},
                       {"left": [0, 2, 3, 5, 6], "right": [4]}], "outputs": [[4, 7]]}),
    (0, [3], {"gates": [{"left": [2, 4, 5], "right": [3]}, {"left": [0, 2, 4, 6], "right": [1, 2, 4, 5, 6]},
                        {"left": [0, 2, 3, 4, 6, 7], "right": [1, 2, 3, 6]}], "outputs": [[3, 4, 5, 6, 8]]}),
    (2, [0, 3], {"gates": [{"left": [0, 4, 5], "right": [1, 3, 4, 6]}, {"left": [0, 2, 3, 5, 7], "right": [4]}],
                 "outputs": [[4, 6, 8]]}),
    (1, [0, 2, 3], {"gates": [{"left": [0, 1, 5, 6], "right": [2, 3, 5, 7]},
                              {"left": [0, 1, 3, 4, 6, 8], "right": [2, 3, 4, 5, 7, 8]}], "outputs": [[3, 5, 6, 7, 9]]}),
    (4, [0, 1, 2, 3], {"gates": [{"left": [1, 5, 8], "right": [4, 7, 8]}], "outputs": [[4, 7, 8, 9]]}),
]


# 10 ANDs within 2 temporaries: stages of 1 or 2 bits solved jointly (exact minima), best by a DP over subsets with
# each stage scheduled within 2 temporaries (output storage): bit 0 (4), bits 3, 4 from bit 0 (3), bits 1, 2 from
# bits 0, 3, 4 (3) [M]. (bits, bits read, program)
CUB_STAGES10 = [
    ([0], [], {"gates": [{"left": [0, 1, 2, 3], "right": [4]}, {"left": [0, 1, 2], "right": [3]},
                         {"left": [0, 3, 5], "right": [1, 4]}, {"left": [1, 7], "right": [3, 5, 7]}],
               "outputs": [[4, 6, 8]]}),
    ([3, 4], [0], {"gates": [{"left": [0, 1, 4, 5], "right": [3, 4, 5]}, {"left": [3], "right": [4, 6]},
                             {"left": [0, 1, 2, 3, 5, 7], "right": [6, 7]}], "outputs": [[8], [6, 7, 8]]}),
    ([1, 2], [0, 3, 4], {"gates": [{"left": [0, 1, 5], "right": [2, 6]}, {"left": [1, 8], "right": [4, 7]},
                                   {"left": [2, 4, 5, 6, 7, 9], "right": [3, 4, 6, 7, 8]}],
                         "outputs": [[5, 7, 8, 10], [6, 9]]}),
]


# 9 ANDs within 4 temporaries: bits 3, 4 jointly (4), then bits 0, 1, 2 from bits 3, 4 (5 = the minimum, k = 4 UNSAT;
# this stage needs 4 temporaries) [M].
CUB_STAGES9 = [
    ([3, 4], [], {"gates": [{"left": [0, 2, 3, 4], "right": [1, 2, 3, 4]}, {"left": [2], "right": [3, 5]},
                            {"left": [3, 6], "right": [4]}, {"left": [2, 3, 4, 5, 7], "right": [6]}],
                  "outputs": [[7], [8]]}),
    ([0, 1, 2], [3, 4], {"gates": [{"left": [0], "right": [1, 3, 5, 6]}, {"left": [1, 2, 6, 7], "right": [3, 4, 5, 6]},
                                   {"left": [1, 4, 5, 6, 8], "right": [3, 5, 6, 7, 8]},
                                   {"left": [0, 3, 5, 8, 9], "right": [1, 5, 7, 8, 9]},
                                   {"left": [1, 4, 6, 8, 9, 10], "right": [3, 4, 5, 6, 7, 8, 10]}],
                         "outputs": [[4, 6, 8, 10], [4, 5, 6, 7, 9, 10], [4, 6, 10, 11]]}),
]


def _emit_cub_fwd(c, tw, W, budget):
    """each stage writes fresh zero output wires, so accumulator storage applies (gates kept on the stage's output
    wires while they hold them alone): same forward count, and the mirror (uncompute) erases those gates free
    (staged 10: uncompute 6 instead of 9) [M]; CUBACC=0 gives the plain single-gate storage."""
    acc = os.environ.get("CUBACC", "1") == "1"
    stages = CUB_STAGES10
    if acc and all(_acc_sched(p, 5 + len(g), budget, 12) is not None for b_, g, p in CUB_STAGES9):
        stages = CUB_STAGES9                       # 9 / 4 within 2 temporaries with accumulator storage [M]
    elif budget >= 4:
        stages = CUB_STAGES9
    k = 0
    for bits, given, prog in stages:
        ins = list(tw) + [W[d] for d in given]
        outs = [W[b] for b in bits]
        if acc:
            sc = _acc_sched(prog, len(ins), budget, 12)
            if sc is not None:
                k += emit_acc(c, prog, ins, outs, budget, sc)
                continue
        sc = schedule(prog, len(ins), budget, zero_outs=True)
        assert sc is not None, "staged table needs 2 temporaries"
        k += emit(c, prog, ins, outs, budget, sc)
    return k


_ACC = {}


def _acc_sched(prog, n, budget, cost_cap):
    key = (json.dumps(prog, sort_keys=True), n, budget, cost_cap)
    if key not in _ACC:
        _ACC[key] = schedule_acc(prog, n, budget, cost_cap)
    return _ACC[key]


def emit_cub(c, tw, W, budget, reverse=False):
    """W[k] ^= bit k of CUB5[t] (t0..t4 = tw), staged 11 ANDs (2 temporaries enough). Forward needs W = 0; reverse
    (W = CUB5[t]) is the kit's exact mirror of a forward emission on the same wires (clears W)."""
    if not reverse:
        return _emit_cub_fwd(c, tw, W, budget)
    save = c.ops
    c.ops = []
    k = _emit_cub_fwd(c, tw, W, budget)
    fwd = c.ops
    c.ops = save
    c.ops.extend(type(c).mirror(fwd))
    return k


def _deps(prog, n):
    out = []
    for g in prog["gates"]:
        out.append(sorted({x - n for x in g["left"] + g["right"] if x >= n}))
    return out


def schedule(prog, n, budget, zero_outs=False):
    """Dijkstra over (live temps, contributed gates): returns (cost, actions) or None. Actions: ("C", g) compute into
    a temp, ("E", g) erase, ("O", g) XOR live g into its outputs, ("D", g) CCX g straight into each of its outputs.
    zero_outs: the outputs start at 0, so a gate that IS an output (output = that gate alone) may be computed into
    the output wire ("S", g: no temporary; it stays live there and is never erased)."""
    if zero_outs:
        return _schedule_store(prog, n, budget)
    G = len(prog["gates"])
    deps = _deps(prog, n)
    users = [set() for _ in range(G)]
    for g, d in enumerate(deps):
        for x in d:
            users[x].add(g)
    outs_of = [set() for _ in range(G)]
    for k, lst in enumerate(prog["outputs"]):
        for x in lst:
            if x >= n:
                outs_of[x - n].add(k)
    need = frozenset(g for g in range(G) if outs_of[g])
    start = (frozenset(), frozenset())
    dist = {start: 0}
    prev = {}
    pq = [(0, 0, start)]
    tie = 0
    while pq:
        d, _, st = heapq.heappop(pq)
        if d != dist.get(st):
            continue
        live, done = st
        if not live and done == need:
            acts = []
            while st != start:
                st, a = prev[st]
                acts.append(a)
            return d, acts[::-1]
        moves = []
        for g in range(G):
            ok = all(x in live for x in deps[g])
            if g in live:
                if ok:
                    moves.append(((live - {g}, done), 0, ("E", g)))
                if g in need and g not in done:
                    moves.append(((live, done | {g}), 0, ("O", g)))
            elif ok:
                if len(live) < budget:
                    moves.append(((live | {g}, done), 1, ("C", g)))
                if g in need and g not in done and not users[g]:
                    moves.append(((live, done | {g}), len(outs_of[g]), ("D", g)))
        for nst, w, a in moves:
            nd = d + w
            if nd < dist.get(nst, 10 ** 9):
                dist[nst] = nd
                prev[nst] = (st, a)
                tie += 1
                heapq.heappush(pq, (nd, tie, nst))
    return None


def _schedule_store(prog, n, budget):
    """schedule() with output storage (outputs start at 0)."""
    G = len(prog["gates"])
    deps = _deps(prog, n)
    users = [set() for _ in range(G)]
    for g, d in enumerate(deps):
        for x in d:
            users[x].add(g)
    outs_of = [set() for _ in range(G)]
    store_of = {}
    for k, lst in enumerate(prog["outputs"]):
        for x in lst:
            if x >= n:
                outs_of[x - n].add(k)
        if len(lst) == 1 and lst[0] >= n and lst[0] - n not in store_of:
            store_of[lst[0] - n] = k
    need = frozenset(g for g in range(G) if outs_of[g])
    start = (frozenset(), frozenset(), frozenset())       # live temps, done contributions, stored gates
    dist = {start: 0}
    prev = {}
    pq = [(0, 0, start)]
    tie = 0
    while pq:
        d, _, st = heapq.heappop(pq)
        if d != dist.get(st):
            continue
        live, done, stored = st
        avail = live | stored
        if not live and done == need:
            acts = []
            while st != start:
                st, a = prev[st]
                acts.append(a)
            return d, acts[::-1]
        moves = []
        for g in range(G):
            ok = all(x in avail for x in deps[g])
            if g in live:
                if ok:
                    moves.append(((live - {g}, done, stored), 0, ("E", g)))
                if g in need and g not in done:
                    moves.append(((live, done | {g}, stored), 0, ("O", g)))
            elif g in stored:
                if g not in done:                          # XOR the stored gate into its other outputs
                    moves.append(((live, done | {g}, stored), 0, ("O", g)))
            elif ok:
                if len(live) < budget:
                    moves.append(((live | {g}, done, stored), 1, ("C", g)))
                if g in store_of and g not in done:        # output k is exactly g: compute it into that wire
                    k = store_of[g]
                    fin = done | {g} if not (outs_of[g] - {k}) else done
                    moves.append(((live, fin, stored | {g}), 1, ("S", g)))
                if g in need and g not in done and not users[g]:
                    moves.append(((live, done | {g}, stored), len(outs_of[g]), ("D", g)))
        for nst, w, a in moves:
            nd = d + w
            if nd < dist.get(nst, 10 ** 9):
                dist[nst] = nd
                prev[nst] = (st, a)
                tie += 1
                heapq.heappush(pq, (nd, tie, nst))
    return None


def schedule_acc(prog, n, budget, cost_cap=40):
    """outputs start at 0 and serve as storage: a gate that is a term of output k may be computed into k while k is
    still zero ("S"), or copied there by CNOT ("O"), and is readable from k until another term is XORed into k.
    Actions: ("C", g) temp, ("E", g) erase temp, ("S", g, k), ("O", g, k), ("D", g, k) CCX into k. Returns
    (cost, actions) or None. The exact uncompute is the kit's mirror of the emission (see emit_acc_uncompute)."""
    G = len(prog["gates"])
    deps = _deps(prog, n)
    nout = len(prog["outputs"])
    terms = [frozenset(x - n for x in lst if x >= n) for lst in prog["outputs"]]
    need = frozenset((g, k) for k in range(nout) for g in terms[k])
    start = (frozenset(), tuple([-1] * nout), frozenset())

    def h(st):                                     # admissible: every pending gate not available costs >= 1
        live_, ost_, done_ = st
        av = live_ | {x for x in ost_ if x >= 0}
        return len({g for (g, k) in need - done_} - av)

    dist = {start: 0}
    prev = {}
    pq = [(h(start), 0, start)]
    tie = 0
    while pq:
        f, _, st = heapq.heappop(pq)
        d = dist.get(st)
        if d is None or f != d + h(st) or d > cost_cap:
            continue
        live, ost, done = st
        if not live and done == need:
            acts = []
            while st != start:
                st, a = prev[st]
                acts.append(a)
            return d, acts[::-1]
        held = {g for g in ost if g >= 0}
        avail = live | held
        moves = []
        for g in range(G):
            ok = all(x in avail for x in deps[g])
            if g in live and ok:
                moves.append(((live - {g}, ost, done), 0, ("E", g)))
            if g in avail:
                for k in range(nout):
                    if g in terms[k] and (g, k) not in done:
                        o2 = list(ost)
                        o2[k] = g if ost[k] == -1 else -2
                        moves.append(((live, tuple(o2), done | {(g, k)}), 0, ("O", g, k)))
            elif ok:
                if len(live) < budget:
                    moves.append(((live | {g}, ost, done), 1, ("C", g)))
                for k in range(nout):
                    if g in terms[k] and (g, k) not in done:
                        o2 = list(ost)
                        if ost[k] == -1:
                            o2[k] = g
                            moves.append(((live, tuple(o2), done | {(g, k)}), 1, ("S", g, k)))
                        if ost[k] >= 0 and ost[k] in deps[g] and ost[k] not in live:
                            continue                       # its factor is read from wire k: no CCX into k
                        o3 = list(ost)
                        o3[k] = g if ost[k] == -1 else -2
                        moves.append(((live, tuple(o3), done | {(g, k)}), 1, ("D", g, k)))
        for nst, w, a in moves:
            nd = d + w
            if nd < dist.get(nst, 10 ** 9) and nd <= cost_cap:
                dist[nst] = nd
                prev[nst] = (st, a)
                tie += 1
                heapq.heappush(pq, (nd + h(nst), tie, nst))
    return None


def emit_acc(c, prog, ins, outs, budget, sched=None):
    """emit with schedule_acc (outputs must be zero here)."""
    n = len(ins)
    if sched is None:
        sched = schedule_acc(prog, n, budget)
    assert sched is not None, "no schedule within %d temporaries" % budget
    cost, acts = sched
    locs = {}                                       # gate -> wires holding exactly it (temps first)
    held_by = {}                                    # output index -> gate it holds exactly (while not mixed)
    temp_of = {}

    def w(x):
        if x < n:
            return ins[x]
        g_ = x - n
        return temp_of[g_] if g_ in temp_of else next(iter(locs[g_]))

    def factors(g):
        L, R = prog["gates"][g]["left"], prog["gates"][g]["right"]
        cl = next(x for x in L if x not in R)
        cr = next(x for x in R if x not in L)
        return L, R, cl, cr

    def form(lst, car):
        for x in lst:
            if x != car:
                c.cx(w(x), w(car))

    def touch(k, g_new):                           # output k receives a term: it holds g_new only if it was zero
        old = held_by.pop(k, None)
        if old is not None:
            locs[old].discard(outs[k])
            if not locs[old]:
                del locs[old]
        if g_new is not None:
            held_by[k] = g_new
            locs.setdefault(g_new, set()).add(outs[k])

    zero = set(range(len(outs)))
    for a in acts:
        g = a[1]
        L, R, cl, cr = factors(g)
        if a[0] == "O":
            k = a[2]
            c.cx(w(n + g), outs[k])
            touch(k, g if k in zero else None)
            zero.discard(k)
            continue
        form(L, cl)
        form(R, cr)
        if a[0] == "C":
            t = c.alloc(1)[0]
            c.and_compute(w(cl), w(cr), t)
            temp_of[g] = t
            locs.setdefault(g, set()).add(t)
        elif a[0] == "E":
            t = temp_of.pop(g)
            c.and_erase(w(cl), w(cr), t)
            locs[g].discard(t)
            if not locs[g]:
                del locs[g]
        elif a[0] == "S":
            k = a[2]
            c.and_compute(w(cl), w(cr), outs[k])
            touch(k, g)
            zero.discard(k)
        else:                                       # "D"
            k = a[2]
            c.ccx(w(cl), w(cr), outs[k])
            touch(k, g if k in zero else None)
            zero.discard(k)
        form(R, cr)
        form(L, cl)
    for k, lst in enumerate(prog["outputs"]):       # input and constant parts
        for x in lst:
            if x == -1:
                c.x(outs[k])
            elif x < n:
                c.cx(ins[x], outs[k])
    assert not temp_of, "temporaries left"
    return cost


def emit_acc_uncompute(c, prog, ins, outs, budget, sched=None):
    """exact inverse of emit_acc (outputs hold the program's values): the kit's mirror of a scratch emission."""
    save = c.ops
    c.ops = []
    k = emit_acc(c, prog, ins, outs, budget, sched)
    fwd = c.ops
    c.ops = save
    c.ops.extend(type(c).mirror(fwd))
    return sum(1 for o in fwd if o[0] in ("ANDE", "CCX"))


def emit(c, prog, ins, outs, budget, sched=None):
    n = len(ins)
    if sched is None:
        sched = schedule(prog, n, budget)
    assert sched is not None, "no schedule within %d temporaries" % budget
    cost, acts = sched
    wire = {}

    def w(x):
        return ins[x] if x < n else wire[x - n]

    def factors(g):
        L, R = prog["gates"][g]["left"], prog["gates"][g]["right"]
        cl = next(x for x in L if x not in R)
        cr = next(x for x in R if x not in L)
        return L, R, cl, cr

    def form(lst, car):
        for x in lst:
            if x != car:
                c.cx(w(x), w(car))

    stored = set()
    for a, g in acts:
        L, R, cl, cr = factors(g)
        if a in ("C", "E", "D", "S"):
            form(L, cl)
            form(R, cr)
            if a == "C":
                t = c.alloc(1)[0]
                c.and_compute(w(cl), w(cr), t)
                wire[g] = t
            elif a == "S":                          # the output that IS this gate (zero now) holds it
                k = next(k for k, lst in enumerate(prog["outputs"]) if lst == [n + g])
                c.and_compute(w(cl), w(cr), outs[k])
                wire[g] = outs[k]
                stored.add(g)
            elif a == "E":
                c.and_erase(w(cl), w(cr), wire[g])
                del wire[g]
            else:
                for k, lst in enumerate(prog["outputs"]):
                    if n + g in lst:
                        c.ccx(w(cl), w(cr), outs[k])
            form(R, cr)
            form(L, cl)
        else:                                       # "O": XOR the live gate into its outputs
            for k, lst in enumerate(prog["outputs"]):
                if n + g in lst and outs[k] != wire[g]:     # a stored gate is already in its own output
                    c.cx(wire[g], outs[k])
    for k, lst in enumerate(prog["outputs"]):       # input and constant parts
        for x in lst:
            if x == -1:
                c.x(outs[k])
            elif x < n:
                c.cx(ins[x], outs[k])
    assert set(wire) == stored, "temporaries left: %s" % sorted(set(wire) - stored)
    return cost


def emit_uncompute(c, prog, ins, outs, budget, zero_outs=True):
    """the exact inverse of emit(..., zero_outs schedule): clears outputs that hold the program's values (the kit's
    mirror of a scratch emission on the same wires: ANDC <-> ANDE)."""
    sc = schedule(prog, len(ins), budget, zero_outs=zero_outs)
    assert sc is not None, "no schedule within %d temporaries" % budget
    save = c.ops
    c.ops = []
    live0, next0 = c.live, c.next
    emit(c, prog, ins, outs, budget, sc)
    fwd = c.ops
    c.ops = save
    c.ops.extend(type(c).mirror(fwd))
    return sc[0]


# six_a joint (c1, c2) in 5 ANDs = the exact minimum (k = 5 SAT, 1,632 s; c1 alone needs 5: k = 4 UNSAT): any start
# value of the outputs; 5 ANDs at 5 temporaries, 6 at 4, 7 at 3.
CHI6_K5 = {"gates": [{"left": [0, 2, 4, 5], "right": [1, 2, 4]}, {"left": [0, 5], "right": [1, 2, 5, 7]},
                     {"left": [1, 6], "right": [2, 7, 8]}, {"left": [2, 6], "right": [9]}, {"left": [2, 9, 10], "right": [3]}],
           "outputs": [[3, 10, 11], [11]]}


def emit_chi6(c, z, l0, l1, c0, out1, out2, outs_zero=True):
    """six_a's carry rule (C >> 1 = c1 + 2 c2 from z = Tf'_6..Tf'_9, l0 = L'_7, l1 = L'_8, c0): out1 ^= c1,
    out2 ^= c2. The cheapest of: the 5-AND minimum (any start value; 5 / 6 / 7 ANDs at 5 / 4 / 3 temporaries), the
    6-AND program, the staged 7-AND form (outs_zero: out1, out2 zero now (forward order) or holding c1, c2)."""
    ins = [z[0], z[1], z[2], z[3], l0, l1, c0]
    budget = c.cap - c.live
    best = None
    for prog in (CHI6_K5, CHI6_K6):
        sc = schedule(prog, 7, budget)
        if sc is not None and (best is None or sc[0] < best[1][0]):
            best = (prog, sc)
    if outs_zero and os.environ.get("CHIACC", "1") == "1":   # accumulator storage (outputs zero): K5 at 5 within
        acc = _acc_sched(CHI6_K5, 7, budget, 7)               # 3 temporaries (plain 7) [M 23:30Z Nov 14]
        if acc is not None and acc[0] <= 7 and (best is None or acc[0] < best[1][0]):
            return emit_acc(c, CHI6_K5, ins, [out1, out2], budget, acc)
    if best is not None and best[1][0] <= 7:
        return emit(c, best[0], ins, [out1, out2], budget, best[1])
    if outs_zero:
        return emit_staged(c, ins, out1, out2, budget, spec="six_a")
    return emit_staged_rev(c, ins, out1, out2, budget, spec="six_a")


def emit_cub_auto(c, tw, W, uncompute=False):
    """six_a's table into W (W zero; uncompute=False) or out of W (W = the table; uncompute=True): CUB8 at >= 5
    temporaries (8 / 7 ANDs at 6, 11 / 10 at 5), else the staged 11-AND form (3 temporaries)."""
    budget = c.cap - c.live
    if budget >= 3 and os.environ.get("CUBACC", "1") == "1":   # CUB8 with accumulator storage: 8 / 4 at 3 [M]
        sc = _acc_sched(CUB8, 5, budget, 9)
        if sc is not None:
            if uncompute:
                return emit_acc_uncompute(c, CUB8, tw, W, budget, sc)
            return emit_acc(c, CUB8, tw, W, budget, sc)
    sc = schedule(CUB8, 5, budget, zero_outs=True)
    if sc is not None and sc[0] < 9:               # CUB8 (8 / 7) when 6 temporaries are free; else the staged form
        if uncompute:
            return emit_uncompute(c, CUB8, tw, W, budget)
        return emit(c, CUB8, tw, W, budget, sc)
    return emit_cub(c, tw, W, budget, reverse=uncompute)


def apply_h1c(H):
    """FC2's emit_chi(c, z, l0, l1, c0, out1, out2) -> the program (out1 ^= c1, out2 ^= c2)."""
    p = os.environ.get("CHIPROG", "")
    progs = [json.load(open(q)) for q in p.split(":")] if p else [CHI_K8]
    orig = H.emit_chi
    log = []

    def emit_chi(c, z, l0, l1, c0, out1, out2):
        budget = c.cap - c.live
        if os.environ.get("CHIMODE", "staged") == "staged" and not p and budget >= 3:
            cost = emit_staged(c, [z[0], z[1], z[2], z[3], l0, l1, c0], out1, out2, budget)
            log.append((budget, cost, "staged"))
            return cost
        best = None
        for prog in progs:                          # the cheapest listed program that fits the budget
            sc = schedule(prog, 7, budget)
            if sc is not None and (best is None or sc[0] < best[1][0]):
                best = (prog, sc)
        if best is None:                            # nothing fits: FC2's rule (17 ANDs)
            log.append((budget, "fc2"))
            return orig(c, z, l0, l1, c0, out1, out2)
        cost = emit(c, best[0], [z[0], z[1], z[2], z[3], l0, l1, c0], [out1, out2], budget, best[1])
        log.append((budget, cost))
        return cost
    H.emit_chi = emit_chi
    H._f3_chi_log = log
    return log


def _ref_chi(z, l, c0):
    """FC2's threshold rule (exact on all y of m 7, DEF7)."""
    l1, l0 = (l >> 1) & 1, l & 1
    t0 = 6 + 3 * c0 + (l1 & l0)
    t1 = 11 + 2 * c0 + ((l1 & l0) if c0 else l1)
    t2 = 16 if c0 else 14 + l1
    return int(z >= t0) + int(z >= t1) + int(z >= t2)


def _test():
    import numpy as np
    sys.path.insert(0, os.path.expanduser("~/qcc/ws/invsqrt_fc/best/w41e/src"))
    import invsqrt_fc_kit as K
    prog = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else CHI_K8
    # occurring patterns: from the problem file if given next to the program, else the threshold rule everywhere
    care = None
    pf = os.path.expanduser("~/qcc/ws/invsqrt_fw/chi/k420/c1.json")
    if os.path.exists(pf):
        d = json.load(open(pf))
        rows = int(d["rows"])
        ins = [int(x["tt"], 16) for x in d["inputs"]]
        care = set()
        for r in range(rows):
            care.add(sum(((ins[i] >> r) & 1) << i for i in range(7)))
    for budget in range(3, 10):
        sc = schedule(prog, 7, budget)
        if sc is None:
            print("budget %d: no schedule" % budget)
            continue
        c = K.Circ(n_fixed=0)
        c.cap = 10 ** 6
        inw = c.alloc(7)
        outw = c.alloc(2)
        live0 = c.live
        st = len(c.ops)
        cost = emit(c, prog, inw, outw, budget, sc)
        ops = c.ops[st:]
        nsh = 128 * 8
        sim = K.Sim(nsh, c.next + 2, seed=budget)
        pat = np.arange(nsh) % 128
        for i in range(7):
            sim.q[inw[i]] = ((pat >> i) & 1).astype(bool)
        rng = np.random.default_rng(budget)
        o0 = rng.integers(0, 2, nsh).astype(bool)
        o1 = rng.integers(0, 2, nsh).astype(bool)
        sim.q[outw[0]] = o0
        sim.q[outw[1]] = o1
        sim.run(ops)
        bad = 0
        for s in range(nsh):
            p = int(pat[s])
            if care is not None and p not in care:
                continue
            z, l0, l1, cc = p & 15, (p >> 4) & 1, (p >> 5) & 1, (p >> 6) & 1
            chi = _ref_chi(z, 2 * l1 + l0, cc)
            if bool(sim.q[outw[0]][s]) != (o0[s] ^ bool(chi & 1)) or bool(sim.q[outw[1]][s]) != (o1[s] ^ bool(chi >> 1)):
                bad += 1
        anc = not sim.q[list(range(live0, c.next))].any() if c.next > live0 else True
        print("budget %d: %d ANDs (expected %.1f), peak temps %d, bad %d, ancillas clean %s, phase %s" % (
            budget, c.toffolis(ops), cost, c.peak - live0, bad, anc, not sim.sign.any()))
    # the staged form: forward from zero outputs, then its mirror must return both outputs to zero
    for budget in (3, 6):
        c = K.Circ(n_fixed=0)
        c.cap = 10 ** 6
        inw = c.alloc(7)
        outw = c.alloc(2)
        live0 = c.live
        st = len(c.ops)
        cost = emit_staged(c, inw, outw[0], outw[1], budget)
        ops = c.ops[st:]
        und = K.Circ.mirror(ops)
        nsh = 128 * 4
        sim = K.Sim(nsh, c.next + 2, seed=budget + 50)
        pat = np.arange(nsh) % 128
        for i in range(7):
            sim.q[inw[i]] = ((pat >> i) & 1).astype(bool)
        sim.run(ops)
        bad = 0
        for s in range(nsh):
            p = int(pat[s])
            if care is not None and p not in care:
                continue
            z, l0, l1, cc = p & 15, (p >> 4) & 1, (p >> 5) & 1, (p >> 6) & 1
            chi = _ref_chi(z, 2 * l1 + l0, cc)
            if bool(sim.q[outw[0]][s]) != bool(chi & 1) or bool(sim.q[outw[1]][s]) != bool(chi >> 1):
                bad += 1
        tmp_clean = not sim.q[list(range(live0, c.next))].any() if c.next > live0 else True
        sim.run(und)
        back = not sim.q[outw].any() and not sim.sign.any()
        print("staged budget %d: %d ANDs (%.1f), peak temps %d, bad %d, temps clean %s, mirror back to zero %s" % (
            budget, c.toffolis(ops), cost, c.peak - live0, bad, tmp_clean, back))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        _test()
    else:
        raise SystemExit(__doc__)
