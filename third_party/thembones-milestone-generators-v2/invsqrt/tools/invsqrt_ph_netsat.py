"""invsqrt_ph: exact search for swap networks of the normaliser with a CP / SAT solver (own code, workstream
invsqrt_ph; OR-tools CP-SAT from a private venv in ~/qcc/ws/invsqrt_ph/venv).

Same model as invsqrt_ph_netsearch.py (reflection idea credited to the best public circuit, manwithacat "HW-25"
= gengsr + maoxin123, read in research/ws_invsqrt_decode.md 3.1 / 3.2), without its restriction that the bits of
an exponent move as one interval: every bit of every exponent is tracked through the layers.

Wires 0..21, exponent e, s = 21 - e. A network is an ordered list of layers; a layer is a set of disjoint swaps.
  bwd (UNNORM on the half thermometer): bits on wires s..21 (leading one on 21) must end on 0..e. Each layer has
      ONE control g(e) = alpha(e >> 1) xor beta (e & 1) (a CNOT sum of half flags and the parity); when it is
      on, all swaps of the layer are done.
  bws: as bwd, but every swap has its own control of that family (CNOTs between swaps are free).
  fwd (NORM): data bits on wires 0..e-1 must end on s..20 (the leading one is cleared by CNOTs). Every swap has
      its own free control over e.
Layer specs (command line, in circuit order): rC = reflection, any of the swaps (q, C - q); pK = patch layer, any
swaps, at most K of them; sK = reflection with a sum chosen by the solver, at most K swaps.
Objective: the number of swaps.

usage: invsqrt_ph_netsat.py bwd|fwd "p2 r25 r21 r20 p2 r18 r12 p2" [seconds=600] [workers=2] [max=60] [out=path]
                            [hint=path of an earlier out file with the same specs] [preset=revfwd]
                            [chains=json of invsqrt_ph_bwdsets.py: reflection layers outside the chain are off]
                            [strict=1: with chains, every chosen swap of a layer in the chain is on]
                            [fix=1: with hint, the swaps of the hint are kept; only others are searched]
                            [leadin=1: the leading one of e = 2k, 2k + 1 enters by CNOTs from the half flags on a
                             free wire at a layer boundary chosen by the solver ("lead_in": [layer, wire] per k)]
                            [leadand=1: the leading one of e = 2k, 2k + 1 may be set by one AND instead ("swaps"
                             in the output then counts swaps + these ANDs, the list is "lead_and")]
"""
import json
import sys
import time

from ortools.sat.python import cp_model

NW = 22


def build(mode, specs, max_swaps, emin=0, lead=True, leadand=False, leadin=False):
    m = cp_model.CpModel()
    # leadand (bwd / bws): the leading one of the exponents 2k, 2k + 1 may be left to one AND of the parity
    # and [e >> 1 = k] after the swaps (sm[k] = 1, price 1) instead of being moved by the swaps.
    sm = [m.NewBoolVar(f"sm{k}") for k in range(11)] if leadand else []
    L = len(specs)
    pairs = []                                   # per layer: list of (a, b)
    for sp in specs:
        if sp[0] == "r":
            c = int(sp[1:])
            pairs.append([(q, c - q) for q in range(NW) if q < c - q <= NW - 1])
        else:
            pairs.append([(a, b) for a in range(NW) for b in range(a + 1, NW)])
    y = [[m.NewBoolVar(f"y{l}_{a}_{b}") for a, b in pairs[l]] for l in range(L)]
    at = []                                      # per layer: wire -> indices of the swaps that touch it
    for l, sp in enumerate(specs):
        d = {q: [] for q in range(NW)}
        for i, (a, b) in enumerate(pairs[l]):
            d[a].append(i)
            d[b].append(i)
        at.append(d)
        for q in range(NW):
            if len(d[q]) > 1:
                m.AddAtMostOne([y[l][i] for i in d[q]])
        if sp[0] in "ps":
            m.Add(sum(y[l]) <= int(sp[1:]))
        if sp[0] == "s":                         # all swaps of the layer have one common sum (a reflection)
            z = {}
            for i, (a, b) in enumerate(pairs[l]):
                if a + b not in z:
                    z[a + b] = m.NewBoolVar(f"z{l}_{a + b}")
                m.AddImplication(y[l][i], z[a + b])
            m.AddAtMostOne(list(z.values()))
    total = sum(sum(v) for v in y) + sum(sm)
    m.Add(total <= max_swaps)
    alpha = beta = None
    act = [[[m.NewBoolVar(f"a{l}_{i}_{e}") for e in range(22)] for i in range(len(pairs[l]))] for l in range(L)]
    if mode == "bwd":
        alpha = [[m.NewBoolVar(f"al{l}_{k}") for k in range(11)] for l in range(L)]
        beta = [m.NewBoolVar(f"be{l}") for l in range(L)]
        for l in range(L):
            for e in range(22):
                if e % 2 == 0:
                    on = alpha[l][e >> 1]
                else:
                    on = m.NewBoolVar(f"on{l}_{e}")
                    m.AddBoolXOr([on, alpha[l][e >> 1], beta[l].Not()])      # on = alpha xor beta
                for i in range(len(pairs[l])):                               # act = y and on
                    m.AddImplication(act[l][i][e], y[l][i])
                    m.AddImplication(act[l][i][e], on)
                    m.AddBoolOr([y[l][i].Not(), on.Not(), act[l][i][e]])
    elif mode == "bws":
        # every swap has its own control of the family alpha(e >> 1) xor beta (e & 1): CNOTs between the swaps
        # are free, so one control per layer is not needed.
        beta = [[m.NewBoolVar(f"bs{l}_{i}") for i in range(len(pairs[l]))] for l in range(L)]
        for l in range(L):
            for i in range(len(pairs[l])):
                for e in range(22):
                    m.AddImplication(act[l][i][e], y[l][i])
                for k in range(11):
                    m.AddBoolXOr([act[l][i][2 * k], act[l][i][2 * k + 1], beta[l][i].Not()])
    else:
        for l in range(L):
            for i in range(len(pairs[l])):
                for e in range(22):
                    m.AddImplication(act[l][i][e], y[l][i])
    done =[[[None] * NW for _ in range(22)] for _ in range(L)]          # some swap at wire q is done for e
    for l in range(L):
        for e in range(22):
            for q in range(NW):
                if at[l][q]:
                    v = m.NewBoolVar(f"d{l}_{e}_{q}")
                    lits = [act[l][i][e] for i in at[l][q]]
                    m.AddBoolOr(lits).OnlyEnforceIf(v)
                    for a_ in lits:
                        m.AddImplication(a_, v)
                    done[l][e][q] = v
    store = {}                                   # e -> bit -> per layer boundary: position variables
    for e in range(emin, 22):
        s = 21 - e
        if mode in ("bwd", "bws"):
            bits = [(i + s, i) for i in range(e)] + ([(21, e)] if lead and not leadin else [])
        else:
            bits = [(i, i + s) for i in range(e)]
        for bi, (q0, q1) in enumerate(bits):
            prev = None
            for l in range(L + 1):
                cur = [m.NewBoolVar(f"x{e}_{bi}_{l}_{q}") for q in range(NW)]
                store.setdefault(e, {}).setdefault(bi, []).append(cur)
                m.AddExactlyOne(cur)
                if l == 0:
                    m.Add(cur[q0] == 1)
                else:
                    ll = l - 1
                    for i, (a, b) in enumerate(pairs[ll]):
                        m.AddBoolOr([prev[a].Not(), act[ll][i][e].Not(), cur[b]])
                        m.AddBoolOr([prev[b].Not(), act[ll][i][e].Not(), cur[a]])
                    for q in range(NW):
                        if done[ll][e][q] is None:
                            m.AddImplication(prev[q], cur[q])
                        else:
                            m.AddBoolOr([prev[q].Not(), done[ll][e][q], cur[q]])
                if l == L:
                    if leadand and lead and bi == e and mode in ("bwd", "bws"):
                        m.Add(cur[q1] == 1).OnlyEnforceIf(sm[e >> 1].Not())
                    else:
                        m.Add(cur[q1] == 1)
                prev = cur
    # leadin (bwd / bws): the leading one need not start on wire 21. It is XORed in by CNOTs from the half
    # flags (free), for both exponents of e >> 1 = k on the same wire inj[k] at the same layer boundary, on a
    # wire that holds no data bit of these exponents there; the swaps after that must bring it to wire e.
    inj = None
    if leadin:
        inj = [[[m.NewBoolVar(f"in{k}_{l}_{q}") for q in range(NW)] for l in range(L + 1)] for k in range(11)]
        for k in range(11):
            m.AddExactlyOne([v for row in inj[k] for v in row])
        for e in range(emin, 22):
            k = e >> 1
            prev = None
            for l in range(L + 1):
                cur = [m.NewBoolVar(f"ld{e}_{l}_{q}") for q in range(NW)]
                m.AddAtMostOne(cur)
                m.Add(sum(cur) == sum(v for l2 in range(l + 1) for v in inj[k][l2]))
                for q in range(NW):
                    m.AddImplication(inj[k][l][q], cur[q])
                    for bi in store.get(e, {}):
                        m.AddImplication(inj[k][l][q], store[e][bi][l][q].Not())
                if l > 0:
                    ll = l - 1
                    for i, (a, b) in enumerate(pairs[ll]):
                        m.AddBoolOr([prev[a].Not(), act[ll][i][e].Not(), cur[b]])
                        m.AddBoolOr([prev[b].Not(), act[ll][i][e].Not(), cur[a]])
                    for q in range(NW):
                        if done[ll][e][q] is None:
                            m.AddImplication(prev[q], cur[q])
                        else:
                            m.AddBoolOr([prev[q].Not(), done[ll][e][q], cur[q]])
                if l == L:
                    if leadand:
                        m.Add(cur[e] == 1).OnlyEnforceIf(sm[k].Not())
                    else:
                        m.Add(cur[e] == 1)
                prev = cur
    m.Minimize(total)
    return m, pairs, y, alpha, beta, act, sm, inj


def main():
    mode = sys.argv[1]
    specs = sys.argv[2].split()
    opt = dict(a.split("=", 1) for a in sys.argv[3:])
    secs, workers, mx = float(opt.get("seconds", 600)), int(opt.get("workers", 2)), int(opt.get("max", 60))
    t0 = time.time()
    m, pairs, y, alpha, beta, act, sm, inj = build(mode, specs, mx, int(opt.get("emin", 0)),
                                                   opt.get("lead", "1") == "1", opt.get("leadand") == "1",
                                                   opt.get("leadin") == "1")
    if opt.get("preset") == "revfwd":
        # Chains of the exponents 8..21 taken from our forward network (tools/invsqrt_ph_net.py) in reverse
        # order: reflections 25, 21, 20, 18, 12 with layer controls that are already sums of half flags and
        # parity (beta = 1 for 21 and 20). Only the exponents 0..7, the leading one and the patches are open.
        use = {25: [8, 9, 10, 11, 14, 15, 16, 17], 21: [9, 11, 12, 15, 17, 18, 20], 20: [9, 10, 13, 15, 16, 19, 20],
               18: [10, 11, 14, 15, 18, 19], 12: [8, 9, 10, 11, 12, 13]}
        for l, sp in enumerate(specs):
            if sp[0] == "r" and int(sp[1:]) in use:
                c = int(sp[1:])
                m.Add(beta[l] == int(c in (21, 20)))
                for k in range(4, 11):
                    m.Add(alpha[l][k] == int(2 * k in use[c]))
    if "chains" in opt:
        # chains of a sum set (invsqrt_ph_bwdsets.py sums=... dump=...): for the exponents listed there, a
        # reflection layer that is not in the chain of the exponent has every swap off; the swaps of the layers
        # in the chain, the patch layers and the other exponents stay free.
        ch = json.load(open(opt["chains"]))["chains"]
        for l, sp in enumerate(specs):
            if sp[0] != "r":
                continue
            for e_s, on in ch.items():
                if int(sp[1:]) not in on:
                    for i in range(len(pairs[l])):
                        m.Add(act[l][i][int(e_s)] == 0)
                elif opt.get("strict") == "1":           # layer in the chain: every chosen swap of it is on
                    for i in range(len(pairs[l])):
                        m.Add(act[l][i][int(e_s)] == y[l][i])
    if "hint" in opt:
        # hint: an earlier network; its layers are matched in order to layers of this run that allow their swaps
        h = json.load(open(opt["hint"]))
        sel, l = set(), 0
        for lay in h["layers"]:
            need = {tuple(p) for p in lay["pairs"]}
            def fits(l):
                if specs[l][0] != lay["spec"][0] or not need <= set(pairs[l]):
                    return False
                return specs[l][0] == "r" or len(need) <= int(specs[l][1:])
            while l < len(specs) and not fits(l):
                l += 1
            assert l < len(specs), "hint does not fit the layer specs"
            sel |= {(l, p) for p in need}
            if mode == "bwd":
                for k in range(11):
                    m.AddHint(alpha[l][k], lay["alpha"][k])
                m.AddHint(beta[l], lay["beta"])
            l += 1
            if mode == "bws" and "alpha" in lay and not isinstance(lay["beta"], list):
                for i, p in enumerate(pairs[l - 1]):
                    if p in need:
                        for e in range(22):
                            m.AddHint(act[l - 1][i][e], lay["alpha"][e >> 1] ^ (lay["beta"] & e & 1))
            elif mode == "bws" and "alpha" in lay:           # hint written by mode bws: controls per swap
                ctl = {tuple(pr): (al, be) for pr, al, be in zip(lay["pairs"], lay["alpha"], lay["beta"])}
                for i, p in enumerate(pairs[l - 1]):
                    if p in need:
                        al, be = ctl[p]
                        for e in range(22):
                            m.AddHint(act[l - 1][i][e], al[e >> 1] ^ (be & e & 1))
        for l in range(len(specs)):
            for i, p in enumerate(pairs[l]):
                m.AddHint(y[l][i], int((l, p) in sel))
                if opt.get("fix") == "1" and (l, p) in sel:      # the swaps of the hint are in the network
                    m.Add(y[l][i] == 1)
    if inj is not None and "hint" in opt:            # hint networks without "lead_in" start it on wire 21
        hl = json.load(open(opt["hint"]))
        where = hl.get("lead_in") or [[0, 21]] * 11
        for k in range(11):
            for l in range(len(specs) + 1):
                for q in range(NW):
                    m.AddHint(inj[k][l][q], int([l, q] == list(where[k])))
            if sm:
                m.AddHint(sm[k], int(k in hl.get("lead_and", [])))
    print("model built in %.0fs" % (time.time() - t0), flush=True)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = secs
    solver.parameters.num_search_workers = workers
    solver.parameters.log_search_progress = False

    class CB(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            print("solution %d swaps%s at %.0fs" % (int(self.ObjectiveValue()), " and ANDs" if sm else "",
                                                   time.time() - t0), flush=True)

    st = solver.Solve(m, CB())
    name = solver.StatusName(st)
    ok = st in (cp_model.OPTIMAL, cp_model.FEASIBLE)
    print("status", name, "objective", solver.ObjectiveValue() if ok else None,
          "bound", solver.BestObjectiveBound(), "time %.0fs" % (time.time() - t0), flush=True)
    if not ok:
        return
    all_layers = []
    for l, sp in enumerate(specs):
        sel = [i for i in range(len(pairs[l])) if solver.Value(y[l][i])]
        lay = {"spec": sp, "pairs": [list(pairs[l][i]) for i in sel]}
        if mode == "bwd":
            lay["alpha"] = [int(solver.Value(v)) for v in alpha[l]]
            lay["beta"] = int(solver.Value(beta[l]))
        elif mode == "bws":                      # per swap: alpha[k] = control for e = 2k, beta = its parity term
            lay["alpha"] = [[int(solver.Value(act[l][i][2 * k])) for k in range(11)] for i in sel]
            lay["beta"] = [int(solver.Value(beta[l][i])) for i in sel]
        else:
            lay["on"] = [[e for e in range(22) if solver.Value(act[l][i][e])] for i in sel]
        all_layers.append(lay)
    net = {"dir": mode, "specs": specs, "status": name, "swaps": int(solver.ObjectiveValue()),
           "lead_and": [k for k in range(11) if sm and solver.Value(sm[k])],
           "lead_in": None if inj is None else [[[l, q] for l in range(len(specs) + 1) for q in range(NW)
                                                 if solver.Value(inj[k][l][q])][0] for k in range(11)],
           "bound": solver.BestObjectiveBound(), "all_layers": all_layers,
           "layers": [lay for lay in all_layers if lay["pairs"]]}
    print("NETSAT " + json.dumps(net))
    if "out" in opt:
        with open(opt["out"], "w") as fh:
            json.dump(net, fh)


if __name__ == "__main__":
    main()
