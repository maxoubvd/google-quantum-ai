"""Builder C3 (sixth holder, Nov 18): THE ONE CHAIN's assembly. An entry folder from a table's final words and weights
(SUB.npz: the squeezed or refit submission), loaders 1..14 from the chain, k0 from a verified entry of the SAME table
(its batch 0 words must equal SUB's: owners never touch batch 0). No verification here (tools/sh/chem_c3_chain.sh runs
the local and the official verifier).

  python -I chem_c3_chainasm.py OUT --npz SUB.npz --k0 ENTRY_DIR [--note TEXT] K=file.kmx ...

A loader 1..14 without a K= pair is kept from ENTRY_DIR only if its batches K-1, K are word-identical in ENTRY_DIR's
submission and SUB; otherwise the assembly is refused (exit 4). Circuits are written as given (no count comparison).
SOURCE.txt: per-loader lines, the note, and the chemistry credits block of the lead's s065 word for word.
"""
import os
import shutil
import sys

import numpy as np

DATA = os.environ.get("QCC_DATA", os.path.expanduser("~/qcc/data"))
S65 = "<HOME>/Quantum/submissions/s065_famC_r59b_q28t54chz_c3g8f6rp/SOURCE.txt"
J = 5


def credits():
    first = "  CREDIT (added 2026-11-16 12:40Z, research/midprize/CREDITS_and_lineage.md section 3)"
    last = "loader compiler are our own."
    lines = open(S65, encoding="utf-8").read().splitlines()
    c0 = next(n for n, s in enumerate(lines) if s.startswith("dfthcsos/"))
    i = next(n for n in range(c0, len(lines)) if lines[n].startswith(first))
    j = next(n for n in range(i, len(lines)) if lines[n].rstrip().endswith(last))
    return lines[i:j + 1]


def main():
    args = sys.argv[1:]
    opt = {}
    for key in ("--npz", "--k0", "--note"):
        if key in args:
            i = args.index(key)
            opt[key] = args[i + 1]
            del args[i:i + 2]
    out, pairs = args[0], args[1:]
    sub, k0e, note = opt["--npz"], opt["--k0"], opt.get("--note", "")
    sys.path[:0] = [os.path.join(DATA, "kickmix")]
    import kickmix as km
    ws = np.asarray(np.load(sub)["angle_words"], np.int64)
    we = np.asarray(np.load(os.path.join(k0e, "dfthcsos", "submission.npz"))["angle_words"], np.int64)
    rep = {}
    for p in pairs:
        k, f = p.split("=", 1)
        rep[int(k)] = f
    # 10:30Z: a pair 0=k0.kmx (the chain's own k0, compiled on the table's batch 0) makes --k0 a mere template folder
    if 0 not in rep and (ws.shape != we.shape or not np.array_equal(ws[:, 0:J], we[:, 0:J])):
        print("REFUSED: batch 0 of the k0 entry differs from the submission's")
        sys.exit(4)
    for k in range(1, 15):
        if k not in rep and (ws.shape != we.shape or not np.array_equal(ws[:, J * (k - 1):J * k + J], we[:, J * (k - 1):J * k + J])):
            print(f"REFUSED: loader {k} has no circuit and its batches differ from the k0 entry's")
            sys.exit(4)
    dst = os.path.join(out, "dfthcsos")
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    shutil.copytree(os.path.join(k0e, "dfthcsos"), dst)
    shutil.copyfile(sub, os.path.join(dst, "submission.npz"))
    k0txt = (f"k0 compiled on this table's batch 0 (the project's plain compiler tools/r31_chem-table_loadj.py); folder "
             f"template {k0e}" if 0 in rep else f"k0 from the verified entry {k0e} (same batch 0)")
    lines = [f"C3 ONE-CHAIN entry (tools/sh/chem_c3_chain.sh, tools/chem_c3_chainasm.py): words and weights {sub}; {k0txt}:"]
    tot = 0
    for k in range(15):
        if k in rep:
            c = km.Circuit(open(rep[k]).read())
            if c.num_qubits > 128:
                print(f"REFUSED: loader {k} uses {c.num_qubits} qubits")
                sys.exit(4)
            c.to_file(os.path.join(dst, f"angle_k_{k}.kmb"), format="kmb")
            src = rep[k]
        else:
            c = km.Circuit.from_file(os.path.join(dst, f"angle_k_{k}.kmb"))
            src = f"kept from {os.path.basename(k0e)}"
        n = int(c.operation_statistics().num_toffoli_like)
        tot += n
        lines.append(f"  k{k}: {n} ({src}, qubits {c.num_qubits})")
    lines.append(f"per pass {tot}")
    if note:
        lines.append(note)
    lines += [""] + credits()
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "SOURCE.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines[:17]))
    print(f"PASS_COUNT {tot}")


if __name__ == "__main__":
    main()
