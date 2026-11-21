"""R63-CHEM-WIRED: a candidate entry = a verified C2 entry with named loaders replaced by kmx circuits AND the
angle words replaced by a table rewritten at dead words (copy of tools/r58_chem_kept_entry.py; --equal also replaces
a loader of the SAME Toffoli count; WORDS=file.npz gives the new angle_words: every other array of submission.npz
is kept byte for byte, and the tool refuses a table that differs from the base at a word the decoder reads).
Nothing is uploaded.
  python -I r63_chem_wired_entry.py BASE OUT [--equal] WORDS=table.npz K=file.kmx [K=file.kmx ...]
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
    wfile = [p.split("=", 1)[1] for p in pairs if p.startswith("WORDS=")]
    pairs = [p for p in pairs if not p.startswith("WORDS=")]
    sys.path[:0] = [os.path.join(DATA, "kickmix")]
    import kickmix as km
    dst = os.path.join(out, "dfthcsos")
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    shutil.copytree(os.path.join(base, "dfthcsos"), dst)
    rep = {int(p.split("=", 1)[0]): p.split("=", 1)[1] for p in pairs}
    lines = [f"R63-CHEM-WIRED candidate: base {base}; wired loaders with dead-word kill unknowns "
             f"(tools/r63_chem_wired_dc.py); table rewritten at dead words only: {wfile}"]
    if wfile:
        import numpy as np
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from r62_chem_cache_dead import dead_mask
        sub = os.path.join(dst, "submission.npz")
        with np.load(sub, allow_pickle=False) as z:
            arrs = {k: z[k] for k in z.files}
        new = np.load(wfile[0])["angle_words"]
        old = arrs["angle_words"]
        assert new.shape == old.shape, (new.shape, old.shape)
        dead = dead_mask(np.asarray(old, np.int64), np.asarray(arrs["tree_edges"], np.int64), 7)
        diff = np.asarray(new, np.int64) != np.asarray(old, np.int64)
        assert not (diff & ~dead).any(), "the new table differs at a word the decoder reads"
        arrs["angle_words"] = np.asarray(new).astype(old.dtype)
        np.savez(sub, **arrs)
        lines.append(f"  angle_words: {int(diff.sum())} dead words rewritten (of {int(dead.sum())} dead), "
                     f"per batch {[int(diff[:, 5 * k:5 * k + 5].sum()) for k in range(15)]}")
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
