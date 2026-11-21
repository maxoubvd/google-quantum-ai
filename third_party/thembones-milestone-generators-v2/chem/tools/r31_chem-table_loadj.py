"""r31_chem-table_loadj: exact loader prices for other J / word widths / block counts, with OUR compiler.

Design R31-CHEM-TABLE (research/panel/R31-CHEM-TABLE.md). The compiler library tools/chem_loaders_core.py reads the
word width and the words per batch from two module constants (BROT, J = 8, 5); everything else is derived from
the table. This script sets them, re-batches a real table of ours and prices every circuit with the same
deterministic decomposition that the emitter writes (optionally emitting, which re-checks every wire on every
address, and simulating with the verifier's own kickmix check).

Tables are OURS (the entry's words, ~/qcc/ws/chem_ld2/ref_L1p/submission.npz, or a rounding run's rounded.npz).
What is synthetic, and why it is a fair price:
  * "cut": whole 16-address blocks of two-body rows are removed (every remaining block keeps its 16 rows and
    their low addresses, so every per-block code of the table survives). It stands for a table with fewer
    vectors per group; which rows are missing does not matter to the loaders.
  * "carriers": the top bit of every word is rewritten so that top xor second = the slot's quadratic low
    monomial (slot = column mod J). In a real table the same relation is reached by the exact top-bit gauge at
    rounding time; here only its loader price is read (the decoded vectors are NOT preserved by this rewrite).
  * codes tied to the batches of five (the XOR-of-five relations on bits 1 and 2) do not survive a re-batching:
    the J = 6 prices therefore carry NO such relations and are upper bounds for a table rounded with them.

Usage (WSL):
  ~/qcc/.venv/bin/python tools/r31_chem-table_loadj.py <table.npz> --J 6 --brot 8 --cut 5 [--carriers]
        [--effort 1] [--engine auto] [--emit DIR] [--shots 2] [--json out.json] [--ks 0,1,5]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402

QUAD = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]          # the six quadratic monomials of four low bits


def load(path):
    with np.load(path, allow_pickle=False) as z:
        key = "angle_words" if "angle_words" in z.files else "words"
        w = np.array(z[key]).astype(np.int64)
        bits = int(z["bits"]) if "bits" in z.files else 8
    return w, bits


def add_onebody(words, brot, seed=1):
    """Rounding runs store the two-body rows only: put 76 generic rows in front (random words, as dense as the
    real one-body rows; their codes are whatever the table's are NOT, i.e. uncoded: an upper bound)."""
    rng = np.random.default_rng(seed)
    ob = rng.integers(0, 1 << brot, size=(76, words.shape[1]), dtype=np.int64)
    return np.concatenate([ob, words], axis=0)


def cut_blocks(words, n, first_block=6):
    """Remove n whole 16-address blocks of two-body rows, spread over the table."""
    A = words.shape[0]
    full = A // 16
    cand = np.linspace(first_block, full - 1, n).round().astype(int) if n else np.zeros(0, int)
    assert len(set(cand.tolist())) == n and (n == 0 or cand.min() >= 5), cand
    keep = np.ones(A, bool)
    for b in cand:
        keep[16 * b:16 * b + 16] = False
    return words[keep], cand.tolist()


def impose_carriers(words, J, brot):
    """top bit := second bit xor q_slot(l), for every word (slot = column mod J)."""
    w = words.copy()
    A, ncol = w.shape
    l = np.arange(A) & 15
    top, sec = brot - 1, brot - 2
    for t in range(ncol):
        a, b = QUAD[t % J % 6]
        q = ((l >> a) & (l >> b) & 1).astype(np.int64)
        s = (w[:, t] >> sec) & 1
        w[:, t] = (w[:, t] & ~(1 << top)) | ((s ^ q) << top)
    return w


def carrier_report(words, J, brot):
    """How many columns already carry a pure low function in top xor second (as the table stands)."""
    A, ncol = words.shape
    l = np.arange(A) & 15
    n = 0
    for t in range(ncol):
        f = ((words[:, t] >> (brot - 1)) ^ (words[:, t] >> (brot - 2))) & 1
        tt = np.zeros(16, np.int64) - 1
        ok = True
        for x in range(16):
            v = f[l == x]
            if v.size and (v != v[0]).any():
                ok = False
                break
            tt[x] = v[0] if v.size else 0
        n += ok
    return n


def price(words, J, brot, effort=1, engine="auto", lbs=(4, 5, 3), ks=None, emit=None, seed=0):
    core.J, core.BROT = J, brot
    K = -(-words.shape[1] // J)
    res = []
    for k in (range(K) if ks is None else ks):
        last = None
        for lb in lbs:
            try:
                r = core.compile_circuit((words, k, lb, seed, effort, emit is not None, engine))
                break
            except ValueError as exc:
                last = str(exc)
                r = None
        if r is None:
            res.append(dict(k=k, error=last))
            continue
        if emit is not None:
            os.makedirs(emit, exist_ok=True)
            with open(os.path.join(emit, f"angle_k_{k}.kmx"), "w") as fh:
                fh.write(r.pop("kmx"))
        r.pop("log", None)
        res.append(r)
    return res


def simulate(emit_dir, words, shots):
    """The verifier's own simulation and register checks (no chemistry gates)."""
    data = os.environ.get("QCC_DATA", os.path.expanduser("~/qcc/data"))
    sys.path[:0] = [os.path.join(data, "kickmix"), os.path.join(data, "verifiers")]
    from dfthcsos.scoreinfra import circuits as vc
    circ = vc.load_circuits(emit_dir, format="kmx")
    got, layout, tof = vc.recover_words_from_circuits(circ, num_vectors=words.shape[0], num_angles=words.shape[1],
                                                      shots_per_address=shots, seed=12345)
    return dict(words_match=bool(np.array_equal(got, words)), toffoli_both_passes=float(tof),
                selection_bits=len(layout.selection), max_qubits=max(c.num_qubits for c in circ),
                circuits=len(circ))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--J", type=int, default=5)
    ap.add_argument("--brot", type=int, default=0, help="0 = the table's own width")
    ap.add_argument("--cut", type=int, default=0, help="16-address blocks of two-body rows to remove")
    ap.add_argument("--carriers", action="store_true", help="rewrite the top bits to the J-slot carriers")
    ap.add_argument("--effort", type=int, default=1)
    ap.add_argument("--engine", default="auto")
    ap.add_argument("--emit", default=None)
    ap.add_argument("--shots", type=int, default=0)
    ap.add_argument("--json", default=None)
    ap.add_argument("--ks", default=None)
    ap.add_argument("--lb", default="4,5,3")
    args = ap.parse_args()
    t0 = time.time()
    words, bits = load(args.table)
    brot = args.brot or bits
    assert words.max() < (1 << brot), (words.max(), brot)
    if words.shape[0] % 143 == 0:                       # a rounding run: two-body rows only
        words = add_onebody(words, brot)
    removed = []
    if args.cut:
        words, removed = cut_blocks(words, args.cut)
    have = carrier_report(words, args.J, brot)
    if args.carriers:
        words = impose_carriers(words, args.J, brot)
    ks = [int(x) for x in args.ks.split(",")] if args.ks else None
    res = price(words, args.J, brot, args.effort, args.engine, tuple(int(x) for x in args.lb.split(",")), ks,
                args.emit)
    bad = [r for r in res if "error" in r]
    tot = sum(r["toffoli"] for r in res if "error" not in r)
    blocks = -(-words.shape[0] // 16)
    print(f"table {args.table}")
    print(f"rows {words.shape[0]} ({blocks} blocks of 16; removed blocks {removed}), brot {brot}, J {args.J}, "
          f"circuits {-(-75 // args.J)}, carrier columns before rewrite {have}/75, rewrite {args.carriers}, "
          f"effort {args.effort}, engine {args.engine}")
    print(" k lb infra prod total phases carriers clean  outputs(affine,quadratic,generic)")
    for r in res:
        if "error" in r:
            print(f"{r['k']:2d}  DOES NOT FIT: {r['error']}")
        else:
            print(f"{r['k']:2d} {r['lb']:2d} {r['infra']:5d} {r['products']:4d} {r['toffoli']:5d} {r['phases']:4d} "
                  f"{r['carriers']:6d} {r['clean_wires']:5d}   {r['levels']}"
                  + (f"  qubits {r['qubits']}" if "qubits" in r else ""))
    print(f"per pass {tot}{'  (INCOMPLETE: ' + str(len(bad)) + ' circuits do not fit)' if bad else ''}; "
          f"infra {sum(r['infra'] for r in res if 'error' not in r)}, products "
          f"{sum(r['products'] for r in res if 'error' not in r)}; {time.time() - t0:.0f} s")
    out = dict(table=args.table, rows=int(words.shape[0]), blocks=blocks, brot=brot, J=args.J, per_pass=tot,
               incomplete=len(bad), circuits=res, removed=removed, carriers_rewritten=args.carriers,
               effort=args.effort)
    if args.emit and args.shots and not bad and ks is None:
        out["check"] = simulate(args.emit, words, args.shots)
        print("simulation (verifier's kickmix check):", out["check"],
              "| expected Toffolis both passes", 2 * tot)
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
