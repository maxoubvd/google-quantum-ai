"""R55-CHEM-ENGINE2 re-pricer: a candidate entry = a verified C2 entry with named loaders replaced by kmx circuits.

Copies BASE/dfthcsos into OUT/dfthcsos; for every K=FILE.kmx compiles the text with the verifier's kickmix and
replaces angle_k_K.kmb only if it has FEWER Toffoli-like gates and at most 128 qubits. Writes OUT/SOURCE.txt.
Nothing is uploaded; verify with tools/qcc_verify.py OUT --tasks chem [--salt alt-salt-1].
  python -I r55_chem_engine2_rp_entry.py BASE OUT K=file.kmx [K=file.kmx ...]
"""
import os
import shutil
import sys

DATA = os.environ.get("QCC_DATA", os.path.expanduser("~/qcc/data"))


def main():
    base, out, pairs = sys.argv[1], sys.argv[2], sys.argv[3:]
    sys.path[:0] = [os.path.join(DATA, "kickmix")]
    import kickmix as km
    dst = os.path.join(out, "dfthcsos")
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    shutil.copytree(os.path.join(base, "dfthcsos"), dst)
    rep = {int(p.split("=", 1)[0]): p.split("=", 1)[1] for p in pairs}
    lines = [f"R55-CHEM-ENGINE2 re-pricer candidate: base {base}; replaced loaders (staircase, closure where named):"]
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
        if n >= n0 or c.num_qubits > 128:
            tn += n0
            lines.append(f"  k{k}: kept {n0} (candidate {n}, qubits {c.num_qubits})")
            continue
        c.to_file(os.path.join(dst, f"angle_k_{k}.kmb"), format="kmb")
        tn += n
        lines.append(f"  k{k}: {n0} -> {n} ({os.path.basename(rep[k])}, qubits {c.num_qubits})")
    lines.append(f"per pass {tb} -> {tn}")
    with open(os.path.join(out, "SOURCE.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
