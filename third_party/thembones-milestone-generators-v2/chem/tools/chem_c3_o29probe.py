"""chem_c3_o29probe.py: builder C3 (fifth holder, Nov 18 07:50Z). LEDGER ROW O29's KILL TEST on today's kept stage: does
ANOTHER VALID BLOCK SET V move the S5T staircase stage A below 149? C2's vset (tools/chem_ld2_c4pub.py, read-only) returns
the FIRST null vector with block 76 of its GF(2) system; every v = that vector XOR any combination of the other null
vectors is equally valid (2^(dim - 1) of them). For each sampled V: the padded table by C2's pad_table (block 76's delta
bits forced by V; check_padded must read 0), then the same two calls my S5T route makes (tools/sh/chem_c3_s5tp.sh):
tools/chem_c3_bsearch.py (searched basis, D 13,14,15) and R66's tools/r66_chem_plan244_stairA.py on the four staircase
splits; the lowest STAGE A with 'outside D: True' is recorded. Box side; one core a V (taskset), pool of NP workers.
    CHEM_BROT=7 CHEM_HB=77 python -I chem_c3_o29probe.py WORDS.npz --k 5 --n 24 --cores 0,1,...,11 --out DIR
"""
import argparse
import json
import os
import random
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

T = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, T)
import chem_loaders_core as core  # noqa: E402
import chem_ld2_wwidth  # noqa: E402,F401
import chem_ld2_c4 as c4  # noqa: E402
import chem_ld2_c4pub as cp  # noqa: E402
from chem_ld2_hb import HB, LAST  # noqa: E402

PY = os.path.expanduser("~/qcc/.venv_work/bin/python")
SPLITS = ("0,2,3|1", "0,1,2|3", "0,1,3|2", "1,2,3|0")


def null_vs(words, k, rel=False):
    F, valid = cp.ct.delta_bits(words, k, core.BROT - 1)
    rows = [np.ones(HB, np.uint8)]
    pidx = np.arange(HB)
    for j in range(7):
        rows.append(((pidx >> j) & 1).astype(np.uint8))
    for b in range(F.shape[0]):
        for l in range(1, 16):
            if valid[LAST, l]:
                rows.append(F[b, :, l] ^ F[b, :, 0])
    if rel:   # 08:20Z: C2's vset_rel rows (chem_ld2_e17gc4.vset_rel): one per wire combination whose mixed parts cancel
        an = core.analyze(c4.pad_table(words, {k: cp.vset(words, k)}), k, 4)
        nw, nb = len(an.out_wires), an.X0.shape[1]
        rows_t = [sum(1 << q for q in range(nw) if (an.G[q] >> o) & 1) for o in range(nb)]
        A = words.shape[0]
        J = core.J
        prev = words[:, J * (k - 1):J * k] if k else np.zeros((A, J), np.int64)
        d = prev ^ words[:, J * k:J * k + J]
        for u in core.nullspace(rows_t, nw):
            col = np.zeros(A, np.uint8)
            for q in core.bits(u):
                j, b = an.out_wires[q]
                col ^= ((d[:, j] >> b) & 1).astype(np.uint8)
            rows.append(np.array([col[16 * p] if 16 * p < A else 0 for p in range(HB)], np.uint8))
    N = [np.asarray(v, np.uint8) for v in cp.ct.nullspace_gf2(np.array(rows, np.uint8))]
    with76 = [v for v in N if v[LAST]]
    rest = [v for v in N if not v[LAST]]
    # make the rest free of block 76 (already) and pick the first with block 76 as the base
    return with76[0], rest + with76[1:]


def tomask(v):
    return sum(1 << p for p in range(HB) if v[p])


def one(args):
    i, V, words, k, core_id, out = args
    d = os.path.join(out, f"v{i}")
    os.makedirs(d, exist_ok=True)
    wp = c4.pad_table(words, {k: V})
    bad = c4.check_padded(wp, {k: V}).get(k, -1)
    np.savez(os.path.join(d, "pad.npz"), angle_words=wp)
    env = dict(os.environ, CHEM_BROT="7", CHEM_HB="77", OMP_NUM_THREADS="1", NUMBA_NUM_THREADS="1")
    pf = os.path.join(d, "pf.json")
    subprocess.run(["taskset", "-c", str(core_id), "nice", "-n", "10", PY, "-I", os.path.join(T, "chem_c3_bsearch.py"),
                    os.path.join(d, "pad.npz"), "--k", str(k), "--D", "13,14,15", "--tries", "300000", "--seed", "1",
                    "--out", pf], env=env, stdout=open(os.path.join(d, "bsearch.log"), "w"), stderr=subprocess.STDOUT,
                   timeout=900)
    best = None
    if os.path.exists(pf):
        for st in SPLITS:
            lg = os.path.join(d, "sa_" + st.replace(",", "").replace("|", "") + ".log")
            try:
                subprocess.run(["taskset", "-c", str(core_id), "nice", "-n", "10", "timeout", "300", PY, "-I",
                                os.path.join(T, "r66_chem_plan244_stairA.py"), os.path.join(d, "pad.npz"), "--k", str(k),
                                "--pfile", pf, "--rep", "4", "--stages", st, "--kills", "2:2", "--seed", "0",
                                "--out", lg[:-4] + ".json"], env=env, stdout=open(lg, "w"), stderr=subprocess.STDOUT,
                               timeout=400)
            except Exception:
                continue
            t = open(lg).read()
            m = re.findall(r"STAGE A (\d+)", t)
            if m and "outside D: True" in t:
                s = min(int(x) for x in m)
                best = s if best is None else min(best, s)
    rec = {"i": i, "V": hex(V), "size": bin(V).count("1"), "check_padded": bad, "stageA": best}
    print(json.dumps(rec), flush=True)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("words")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--n", type=int, default=24)
    ap.add_argument("--cores", default="0,1,2,3,4,5,6,7,8,9,10,11")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--rel", action="store_true", help="add C2's vset_rel wire rows (S-wire tables)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    core.J = 5
    z = np.load(a.words)
    words = np.asarray(z["angle_words"], np.int64)[:1219]
    base, rest = null_vs(words, a.k, a.rel)
    v0 = tomask(base)
    print(f"k{a.k}: null space {1 + len(rest)} vectors (rel {a.rel}); base V equal to plain vset: {v0 == cp.vset(words, a.k)}", flush=True)
    rng = random.Random(a.seed)
    Vs = [v0]
    seen = {v0}
    while len(Vs) < min(a.n, 1 << len(rest)):
        v = base.copy()
        for r in rest:
            if rng.random() < 0.5:
                v ^= r
        m = tomask(v)
        if m not in seen:
            seen.add(m)
            Vs.append(m)
    cores = [int(x) for x in a.cores.split(",")]
    os.makedirs(a.out, exist_ok=True)
    jobs = [(i, V, words, a.k, cores[i % len(cores)], a.out) for i, V in enumerate(Vs)]
    with ThreadPoolExecutor(max_workers=len(cores)) as ex:
        recs = list(ex.map(one, jobs))
    json.dump(recs, open(os.path.join(a.out, "summary.json"), "w"), indent=1)
    got = [r["stageA"] for r in recs if r["stageA"] is not None]
    print(f"k{a.k}: {len(recs)} V's, stage A min {min(got) if got else None}, default V {recs[0]['stageA']}, "
          f"histogram {sorted(set(got))}", flush=True)


if __name__ == "__main__":
    main()
