"""r51_chem_bit1_rp_obn (R51-CHEM-BIT1 re-pricer, 2026-11-16): the 76 one-body rows of a finished table under R51
layouts, read two ways at FIXED h1' (the two-body part of the table held):

  designer  the design file's reading (tools/r51_chem_bit1_obround.py: the rounding tool's decoder on the layout's
            chains, heavy5 placement, ONE placement, RD.HALF = 64 = the full circle, rows 76-79 free); then the
            chain's own signed-range step (tools/r48_chem_q7c_obsigned.py: an exact gauge move of half turns and
            negations) applied to those words, and the chain links counted that the move breaks, per bit.
  chain     the builder's own chain step for a layout (tools/chem_cd7_table7.build_rows7_layout = chem_cd7_oblay:
            rows 76-79 FROZEN in the decoder, every word in -32..31, polish moves out of range rejected, light4 /
            heavy5 placements alternating, best of N), with rows 76-79 given plane-consistent bits on every chain
            first (the layout rounding makes them so; an uncoded table's are not).

    python r51_chem_bit1_rp_obn.py TABLE.npz --spec SPEC [--spec ...] [--specfile F ...] [--mode both|designer|chain]
                                   [--tries 8] [--json OUT.json] [--out-prefix P]

Builders' files are imported read-only; nothing of theirs is edited.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r51_chem_bit1_round as RD  # noqa: E402
import r51_chem_bit1_obround as DOB  # noqa: E402  (sets RD.HALF = 64 at import; set explicitly below)
import chem_cd7_oblay as OL  # noqa: E402
import chem_cd7_table7 as CT  # noqa: E402
import r48_chem_q7c_obsigned as OS  # noqa: E402
from r48_chem_q7c_onebody import (angles_and_mass, check_fcidump, decode_pi, error1, load_fixed_shifts,  # noqa: E402
                                  load_reference, reconstruct_u_vectors, symmetric_from_lower)
from r51_chem_bit1_codemodel import parse_spec  # noqa: E402

LO, HI = -32, 31


def signed(w):
    w = np.mod(np.asarray(w, np.int64), 128)
    return np.where(w >= 64, w - 128, w)


def table_H(path):
    z = np.load(path)
    words = np.asarray(z["angle_words"], np.int64)
    edges = np.asarray(z["tree_edges"], np.int64)
    w = np.asarray(z["W"], np.float64)
    Wo = np.asarray(z["H2_I"], np.float64) - w.sum(axis=1)
    R, B, C = w.shape
    N = edges.shape[0] + 1
    ref = load_reference(check_fcidump())
    sh = load_fixed_shifts()
    S = np.asarray(symmetric_from_lower(sh["Hsym_lower"], N), np.float64)
    beta = float(sh["H1_I"]) - 0.5 * float(np.trace(S))
    eta = int(ref["electrons"])
    tei = np.asarray(ref["tei"], np.float64)
    h1_thc = np.asarray(ref["oei"], np.float64) - 0.5 * np.einsum("prrq->pq", tei)
    del tei
    allv = reconstruct_u_vectors(words / 64.0, edges)
    allv = allv / np.linalg.norm(allv, axis=1, keepdims=True)
    u = allv[N:].reshape(R, B, N)
    k = np.einsum("rc,rbc->rb", Wo, w)
    H = h1_thc - 0.5 * eta * S + beta * np.eye(N) - np.einsum("rb,rbp,rbq->pq", k, u, u, optimize=True)
    e_own = error1(allv[:N], H)[0]
    return words, edges, H, e_own


def link_violations(wq, codes, w76=None, nrows=76):
    """Chain links not affine per 16-row block, counted per bit; block 4 = rows 64-75 (+ rows 76-79 if given)."""
    w = np.mod(np.asarray(wq, np.int64), 128)
    if w76 is not None:
        w = np.vstack([w, np.mod(np.asarray(w76, np.int64), 128)])
    out = {}
    for bit, _, chains in codes:
        for ch in chains:
            for c0, c1 in zip(ch[:-1], ch[1:]):
                for blk in range(5):
                    rows = np.arange(16 * blk, min(16 * blk + 16, w.shape[0]))
                    d = ((w[rows, c0] ^ w[rows, c1]) >> bit) & 1
                    bad = int(not (DOB.CW[:, :len(rows)] == d[None, :]).all(axis=1).any())
                    out[bit] = out.get(bit, 0) + bad
    return out


def obsigned_move(wq, edges):
    """tools/r48_chem_q7c_obsigned.py's exact gauge move on the 76 rows (same steps; returns the new words)."""
    half, full, q = 64, 128, 32
    w0 = np.mod(np.asarray(wq, np.int64), full)
    v0 = OS.decode_pi(w0 / half, edges)
    v0c = np.where(np.abs(v0) < 1e-12, 0.0, v0)
    th, mass = OS.signed_angles(v0c, edges)
    w1 = np.rint(th * half).astype(np.int64)
    dead = mass < 1e-20
    fold = w0 % half
    fold = np.where(fold > q, half - fold, fold)
    w1 = np.where(dead, np.where(fold == q, -q, fold), w1)
    off = float(np.abs(th * half - w1)[~dead].max())
    v1 = OS.decode_pi(w1 / half, edges)
    sgn = np.sign((v0 * v1).sum(1))
    err = float(np.abs(v1 * sgn[:, None] - v0).max())
    return w1, off, err


def repair76(w76, codes):
    """Rows 76-79 (block 4, low addresses 12-15, an affine 2-plane): a link is satisfiable there only if the XOR of
    the four rows' bit is the same on both columns; flip that bit of row 79 on the columns that disagree (XOR on bits
    0..4 keeps the signed range). The layout rounding makes these bits consistent itself."""
    w = signed(w76).copy()
    flips = 0
    for bit, _, chains in sorted(codes, key=lambda x: x[0]):
        for ch in chains:
            par = [int(np.bitwise_xor.reduce((np.mod(w[:, c], 128) >> bit) & 1)) for c in ch]
            for c, p in zip(ch[1:], par[1:]):
                if p != par[0]:
                    w[3, c] = signed(np.mod(w[3, c], 128) ^ (1 << bit))
                    flips += 1
    return w, flips


def designer_rows(H, edges, spec_text):
    _, codes = parse_spec(spec_text, np.arange(75))
    RD.HALF = 64
    tool = None
    if codes:
        deg, chain = RD.layout_arrays([[b_, 1, ch] for b_, d_, chains in codes for ch in chains], (75, 8))
        RD.CHAIN = chain
        tool = (deg, np.zeros((15, 8), np.int64))
    else:
        RD.CHAIN = None
    N = edges.shape[0] + 1
    lam_e, vec = np.linalg.eigh(H)
    th0, _ = angles_and_mass(vec.T, edges)
    row_at = CT._placement(lam_e, N, 0, "heavy5")
    th0, lam_at = th0[row_at], lam_e[row_at]
    wq, nf, nm = DOB.progressive_layout(th0, edges, lam_at, H, 64.0, 300, None, tool)
    step = np.ones(75, np.int64)
    for bit, _, chains in codes:
        for ch in chains:
            for c in ch:
                step[c] = max(step[c], 1 << (bit + 1))
    wq2, e_f = DOB.polish(wq, edges, H, 7, 40, step)
    ws = signed(wq2)
    outside = int(np.sum((ws < LO) | (ws > HI)))
    rows_out = int(np.sum(((ws < LO) | (ws > HI)).any(axis=1)))
    viol0 = link_violations(wq2, codes)
    w1, off, err = obsigned_move(wq2, edges)
    e_after = error1(decode_pi(np.mod(w1, 128) / 64.0, edges), H)[0]
    viol1 = link_violations(w1, codes)
    neg = int(np.sum((signed(w1) != ws) & (np.mod(signed(w1) + ws, 128) == 0) & (ws != 0)))
    changed = int(np.sum(np.mod(w1, 128) != np.mod(wq2, 128)))
    return dict(e1_full_circle_mHa=1e3 * e_f, words_outside=outside, rows_with_words_outside=rows_out,
                violations_before=viol0, words_changed_by_move=changed, words_negated=neg,
                move_exact=(off < 1e-6 and err < 1e-9), e1_after_move_mHa=1e3 * e_after,
                violations_after_move=viol1, links=sum(len(ch) - 1 for _, _, chs in codes for ch in chs))


def chain_rows(H, edges, words, spec_text, tries):
    codes, deg, step, coded = OL.setup(spec_text if spec_text != "none" else "none")
    w76, flips = repair76(words[76:80], codes)
    wq, e, _, info = CT.build_rows7_layout(H, edges, 7, 300, 40, deg, step, codes, tries, w76)
    viol = link_violations(wq, codes, w76)
    ws = signed(wq)
    info.update(dict(e1_best_mHa=1e3 * e, flips76_repair=flips, violations_full=viol,
                     outside_signed_best=int(np.sum((ws < LO) | (ws > HI)))))
    return info, wq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--spec", action="append", default=[])
    ap.add_argument("--specfile", action="append", default=[])
    ap.add_argument("--item", action="append", default=[], help="in order: a spec, or a file ending in .spec")
    ap.add_argument("--mode", default="both", choices=["both", "designer", "chain"])
    ap.add_argument("--tries", type=int, default=8)
    ap.add_argument("--json", default="")
    ap.add_argument("--out-prefix", default="")
    a = ap.parse_args()
    import multiprocessing as mp
    mp.set_start_method("fork", force=True)      # CD7's pool must inherit RD.CHAIN / RD.HALF set by OL.setup
    t0 = time.time()
    words, edges, H, e_own = table_H(a.table)
    print(f"table {a.table}: its own one-body rows error_1 {1e3 * e_own:.4f} mHa ({time.time() - t0:.0f} s)", flush=True)
    specs = [(os.path.basename(x), open(x).read().strip()) if x.endswith(".spec") else (x, x) for x in a.item]
    specs += [(s, s) for s in a.spec] + [(os.path.basename(f), open(f).read().strip()) for f in a.specfile]
    res = {}
    for name, text in specs:
        r = {}
        if a.mode in ("both", "designer"):
            r["designer"] = designer_rows(H, edges, text)
            print(f"[{name}] DESIGNER (full circle, heavy5, one placement): {json.dumps(r['designer'])} "
                  f"({time.time() - t0:.0f} s)", flush=True)
        if a.mode in ("both", "chain"):
            info, wq = chain_rows(H, edges, words, text, a.tries)
            r["chain"] = info
            print(f"[{name}] CHAIN (CD7 build_rows7_layout, best of {a.tries}): {json.dumps(info)} "
                  f"({time.time() - t0:.0f} s)", flush=True)
            if a.out_prefix:
                np.save(f"{a.out_prefix}_{name.replace(':', '_').replace(',', '_')[:40]}.npy", wq)
        res[name] = r
        if a.json:
            with open(a.json, "w") as f:
                json.dump(res, f, indent=1)
    print(f"done ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
