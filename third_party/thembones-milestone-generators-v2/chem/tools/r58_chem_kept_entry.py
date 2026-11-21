"""R58-CHEM-KEPT: a candidate entry = a verified C2 entry with named loaders replaced by kmx circuits (copy of
tools/r55_chem_engine2_rp_entry.py; --equal also replaces a loader of the SAME Toffoli count, for the test entry
that proves a kept-product circuit in the official verifier). Nothing is uploaded.
  python -I r58_chem_kept_entry.py BASE OUT [--equal] K=file.kmx [K=file.kmx ...]
then  ~/qcc/.venv/bin/python tools/qcc_verify.py OUT --tasks chem [--salt alt-salt-1]
"""
import os
import shutil
import sys

DATA = os.environ.get("QCC_DATA", os.path.expanduser("~/qcc/data"))


def main():
    args = sys.argv[1:]
    equal = "--equal" in args
    args = [x for x in args if x != "--equal"]
    base, out, pairs = args[0], args[1], args[2:]
    sys.path[:0] = [os.path.join(DATA, "kickmix")]
    import kickmix as km
    dst = os.path.join(out, "dfthcsos")
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    shutil.copytree(os.path.join(base, "dfthcsos"), dst)
    rep = {int(p.split("=", 1)[0]): p.split("=", 1)[1] for p in pairs}
    lines = [f"R58-CHEM-KEPT candidate: base {base}; loaders replaced by kept-product circuits "
             f"(tools/r58_chem_kept_emit.py):"]
    tb = tn = 0
    for k in range(15):
        c0 = km.Circuit.from_file(os.path.join(dst, f"angle_k_{k}.kmb"))
        n0 = int(c0.operation_statistics().num_toffoli_like)
        tb += n0
        if k not in rep:
            tn += n0
            lines.append(f"  k{k}: kept {n0}")
            continue
        c = km.Circuit(open(rep[k]).read())
        n = int(c.operation_statistics().num_toffoli_like)
        if n > n0 or (n == n0 and not equal) or c.num_qubits > 128:
            tn += n0
            lines.append(f"  k{k}: kept {n0} (candidate {n}, qubits {c.num_qubits})")
            continue
        c.to_file(os.path.join(dst, f"angle_k_{k}.kmb"), format="kmb")
        tn += n
        lines.append(f"  k{k}: {n0} -> {n} ({os.path.basename(rep[k])}, qubits {c.num_qubits})")
    lines.append(f"per pass {tb} -> {tn}")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "SOURCE.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
