#!/usr/bin/env python3
"""R91-INV-STEP part 2: executed-gate audit of fc18_M by classically conditioned block.

Two independent counts of the executed Toffoli-like gates (CCX, CCZ, CSWAP) of every conditioned block:
  [P] EXACT: every classical bit of the file is an HMR outcome (a fair coin, independent of the data: kickmix's
      'HMR q b: b := random') or a BIT_INVERT of such bits, so each op's condition is a Boolean function of fair
      coins; it is tracked here as a GF(2) polynomial (ANF) in the coins and its probability counted exactly.
  [M] THE VERIFIER'S OWN COUNTING: kickmix.Simulator(counting_operations=True), as in
      kaggle_data/verifiers/invsqrt/verify_submission_invsqrt.py (imported as a library: spec, reference values,
      bit packing, the phase check), run segment by segment (one segment a block) so the stats delta of a segment
      is that block's executed count. Every shot is checked as the verifier checks it (x restored, scratch 0,
      phase within tol), which proves the segmentation changes nothing.
Usage: python r91_inv_step_exec.py KMB OUTDIR [n_rand=4096] [shots=64] [batch=4096]
"""
import hashlib
import io
import itertools
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, "<HOME>/qcc/data/kickmix")
sys.path.insert(0, "<HOME>/Quantum/kaggle_data/verifiers")
import kickmix as km  # noqa: E402
from invsqrt import verify_submission_invsqrt as V  # noqa: E402

TOF = ("CCX", "CCZ", "CSWAP")
GAP = 8  # unconditioned Clifford ops allowed inside one block (between two ops on the same outer coin)


# ---------- GF(2) polynomials in coin variables: a frozenset of monomials, a monomial = frozenset of coin ids ----------
ONE = frozenset([frozenset()])
ZERO = frozenset()


def padd(a, b):
    return frozenset(a ^ b)


def pmul(a, b):
    out = set()
    for m in a:
        for n in b:
            mn = m | n
            if mn in out:
                out.remove(mn)
            else:
                out.add(mn)
    return frozenset(out)


_pcache = {}


def prob(p):
    """Exact P(p = 1) over independent fair coins (enumeration of the support)."""
    if p in _pcache:
        return _pcache[p]
    vs = sorted(set().union(*p)) if p else []
    if len(vs) > 22:
        raise RuntimeError("support too large: %d" % len(vs))
    idx = {v: i for i, v in enumerate(vs)}
    k = len(vs)
    n = 1 << k
    a = np.arange(n, dtype=np.int64)
    val = np.zeros(n, dtype=bool)
    for m in p:
        t = np.ones(n, dtype=bool)
        for v in m:
            t &= ((a >> idx[v]) & 1).astype(bool)
        val ^= t
    r = float(val.mean())
    _pcache[p] = r
    return r


def parse(text):
    """Ops as dicts: name, args (list of str), cond (bit id or None), line index."""
    ops = []
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        cond = None
        if " if " in s:
            s, c = s.split(" if ")
            cond = int(c.strip()[1:])
        parts = s.split()
        ops.append({"name": parts[0], "args": parts[1:], "cond": cond})
    return ops


def coin_algebra(ops):
    """Per op: the exact probability its condition holds, the outer coin bit, the depth, the condition poly."""
    bits = {}          # bit id -> poly
    stack = []         # list of (bit id, poly)
    out = []
    for i, op in enumerate(ops):
        nm = op["name"]
        sp = ONE
        for _, p in stack:
            sp = pmul(sp, p)
        own = bits.get(op["cond"], ZERO) if op["cond"] is not None else ONE
        cp = pmul(sp, own)
        outer = stack[0][0] if stack else op["cond"]
        rec = {"depth": len(stack), "outer": outer, "cond": op["cond"],
               "stack": [b for b, _ in stack]}
        if nm == "PUSH_CONDITION":
            stack.append((op["cond"], bits.get(op["cond"], ZERO)))
            rec["p"] = None
        elif nm == "POP_CONDITION":
            stack.pop()
            rec["p"] = None
        else:
            rec["p"] = prob(cp) if cp not in (ONE, ZERO) else (1.0 if cp == ONE else 0.0)
            if nm == "HMR":
                b = int(op["args"][1][1:])
                coin = frozenset([frozenset([b])])
                old = bits.get(b, ZERO)
                # executed under cp: b := coin; else b keeps old
                bits[b] = padd(pmul(cp, coin), pmul(padd(ONE, cp), old))
            elif nm == "BIT_INVERT":
                b = int(op["args"][0][1:])
                bits[b] = padd(bits.get(b, ZERO), cp)
        rec["poly"] = cp
        out.append(rec)
    return out, bits


def blocks_of(ops, alg):
    """Conditioned blocks: maximal runs of conditioned ops sharing the outer coin, joined over gaps of at most GAP
    unconditioned Clifford ops. Returns list of dicts with start, end (inclusive op index)."""
    blocks = []
    cur = None
    last_cond_i = None
    for i, (op, a) in enumerate(zip(ops, alg)):
        conditioned = a["depth"] > 0 or op["cond"] is not None or op["name"] in ("PUSH_CONDITION", "POP_CONDITION")
        if op["name"] == "PUSH_CONDITION":
            outer = a["outer"] if a["depth"] > 0 else op["cond"]
        else:
            outer = a["outer"]
        if conditioned:
            if (cur is not None and cur["outer"] == outer and i - last_cond_i - 1 <= GAP
                    and all(ops[j]["name"] not in TOF and ops[j]["name"] != "HMR" for j in range(last_cond_i + 1, i))):
                cur["end"] = i
            else:
                if cur is not None:
                    blocks.append(cur)
                cur = {"start": i, "end": i, "outer": outer}
            last_cond_i = i
    if cur is not None:
        blocks.append(cur)
    return blocks


def main():
    kmb, outdir = sys.argv[1], sys.argv[2]
    n_rand = int(sys.argv[3]) if len(sys.argv) > 3 else 4096
    shots_per_x = int(sys.argv[4]) if len(sys.argv) > 4 else 64
    batch = int(sys.argv[5]) if len(sys.argv) > 5 else 4096
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    circ = km.Circuit.from_file(kmb, format="kmb")
    sha = V._circuit_sha256(circ)
    print("sha256", sha, "ops", len(circ), "qubits", circ.num_qubits, flush=True)
    kmx_path = os.path.join(outdir, "file.kmx")
    circ.to_file(kmx_path, format="kmx")
    text = open(kmx_path).read()
    # round trip: the text must rebuild the same circuit (same kmb sha)
    c2 = km.Circuit(text)
    sha2 = V._circuit_sha256(c2)
    print("text round trip sha", sha2, "SAME" if sha2 == sha else "DIFFERENT", flush=True)
    ops = parse(text)
    lines = [l for l in text.splitlines() if l.strip()]
    assert len(lines) == len(ops)
    alg, bits = coin_algebra(ops)
    exp_tof = sum(a["p"] for op, a in zip(ops, alg) if op["name"] in TOF)
    static_tof = sum(1 for op in ops if op["name"] in TOF)
    print("static Toffoli-like", static_tof, "EXACT expected executed [P]", round(exp_tof, 6), flush=True)
    blocks = blocks_of(ops, alg)
    tb = [b for b in blocks if any(ops[j]["name"] in TOF for j in range(b["start"], b["end"] + 1))]
    print("conditioned blocks", len(blocks), "with Toffoli-like", len(tb), flush=True)

    # ---- segments: cut at every Toffoli block's start and end+1 ----
    cuts = sorted(set([0, len(ops)] + [b["start"] for b in tb] + [b["end"] + 1 for b in tb]))
    segs = [(cuts[k], cuts[k + 1]) for k in range(len(cuts) - 1)]
    # the register header lines (APPEND_TO_REGISTER / REGISTER) sit at the top: keep them in segment 0 only
    seg_circs = [km.Circuit("\n".join(lines[a:b])) for a, b in segs]
    seg_of_block = {}
    for bi, b in enumerate(tb):
        seg_of_block[segs.index((b["start"], b["end"] + 1))] = bi

    spec = V.ProblemSpecInvSqrt.for_n_p(V.N_P)
    rng = np.random.default_rng(int(sha[:16], 16))
    xs_in = np.concatenate([np.arange(1, 4097, dtype=np.uint64),
                            rng.integers(1, spec.max_x + 1, size=n_rand, dtype=np.uint64)])
    shots = np.repeat(xs_in, shots_per_x)
    target, tol = V.reference_values(spec, np.arange(1, spec.max_x + 1, dtype=np.uint64))
    x_reg, scratch = V.validate_registers(circ, spec)
    seeds = np.random.SeedSequence(int(sha, 16), pool_size=8).generate_state(
        4 * ((len(shots) + batch - 1) // batch), dtype=np.uint64).reshape(-1, 4)
    seg_tof = np.zeros(len(segs), dtype=np.float64)
    nsh = 0
    worst = 0.0
    for bi_, start in enumerate(range(0, len(shots), batch)):
        xs = shots[start:start + batch]
        if len(xs) < batch:
            xs = np.concatenate([xs, xs[: batch - len(xs)]])  # pad (counted, harmless)
        sim = km.Simulator(batch_size=batch, debug_print=False, counting_operations=True,
                           seed=int.from_bytes(seeds[bi_], "little"))
        sim.clear_for_shot()
        sim.write_across_shots(x_reg, V.to_bits(xs, spec.width))
        prev = 0
        for k, sc in enumerate(seg_circs):
            sim.do(sc)
            now = int(sim.operation_statistics().num_toffoli_like)
            seg_tof[k] += now - prev
            prev = now
        if not np.array_equal(V.from_bits(sim.read_across_shots(x_reg)), xs):
            raise AssertionError("x was not restored")
        if sim.read_across_shots(scratch).any():
            raise AssertionError("scratch dirty")
        theta = np.array([V.half_turns_to_units(sim.read_shot_phase(i)) for i in range(batch)], dtype=np.uint64)
        allowed = tol[xs - np.uint64(1)]
        error = (theta - target[xs - np.uint64(1)]).view(np.int64)
        used = np.abs(error).view(np.uint64)
        if (used > allowed).any():
            raise AssertionError("phase out of tolerance")
        worst = max(worst, float(np.divide(used, allowed, out=np.zeros(batch), where=allowed > 0).max()))
        nsh += batch
        if bi_ % 16 == 0:
            print("batch", bi_, "shots", nsh, "tof/shot so far", round(seg_tof.sum() / nsh, 4),
                  "t", round(time.time() - t0, 1), flush=True)
    print("ALL SHOTS PASS the verifier's checks: shots", nsh, "worst err/tol", round(worst, 5), flush=True)
    meas_total = seg_tof.sum() / nsh
    print("MEASURED executed Toffoli-like a shot [M]", round(meas_total, 5), flush=True)

    # ---- per block table ----
    rows = []
    for k, (a, b) in enumerate(segs):
        if k not in seg_of_block:
            continue
        blk = tb[seg_of_block[k]]
        tofs = [j for j in range(a, b) if ops[j]["name"] in TOF]
        exact = sum(alg[j]["p"] for j in tofs)
        fr = sorted(set(round(alg[j]["p"], 6) for j in tofs))
        polys = sorted(set(str(sorted(sorted(m) for m in alg[j]["poly"])) for j in tofs))
        names = {}
        for j in tofs:
            names[ops[j]["name"]] = names.get(ops[j]["name"], 0) + 1
        hmr_in = [j for j in range(a, b) if ops[j]["name"] == "HMR"]
        rows.append({"start": a, "end": b - 1, "outer": blk["outer"], "static_tof": len(tofs),
                     "exec_exact": exact, "exec_meas": seg_tof[k] / nsh, "fractions": fr,
                     "cond_polys": polys[:4], "n_polys": len(polys), "gates": names, "n_ops": b - a,
                     "hmr_inside": len(hmr_in), "first_line": lines[a], "lines": lines[a:b]})
    rows.sort(key=lambda r: -r["exec_exact"])
    uncond = sum(1 for op, al in zip(ops, alg) if op["name"] in TOF and al["p"] == 1.0)
    per_op = [{"i": i, "name": op["name"], "p": alg[i]["p"], "depth": alg[i]["depth"], "outer": alg[i]["outer"]}
              for i, op in enumerate(ops)]
    res = {"sha": sha, "static_tof": static_tof, "exact_exec": exp_tof, "meas_exec": meas_total, "shots": nsh,
           "worst": worst, "n_blocks": len(blocks), "n_tof_blocks": len(tb),
           "uncond_tof_ops": uncond, "blocks": rows}
    json.dump(res, open(os.path.join(outdir, "exec_blocks.json"), "w"), indent=1)
    json.dump(per_op, open(os.path.join(outdir, "per_op_p.json"), "w"))
    print("uncond Toffoli-like ops (p = 1)", uncond, "; conditioned static", static_tof - uncond,
          "executed", round(exp_tof - uncond, 4))
    print("TOP 20 blocks by exact executed Toffoli-like a shot:")
    for r in rows[:20]:
        print("ops %d-%d outer b%s static %d exact %.4f meas %.4f fr %s gates %s hmr_in %d polys %d %s" % (
            r["start"], r["end"], r["outer"], r["static_tof"], r["exec_exact"], r["exec_meas"], r["fractions"],
            r["gates"], r["hmr_inside"], r["n_polys"], r["cond_polys"][:2]))
    print("time", round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
