"""C2 (22nd holder): pick the cheapest emitted circuit per loader in a pulled run folder of tools/sh/chem_ld2_e22dcbox.sh
(owners o<K>/angle_k_<K>.kmx, wired neighbours n<K>/angle_k_<K>.kmx, 30-output neighbours b<K>/*.kmx), counted with the
verifier's kickmix (Toffoli-like gates, at most 128 qubits). Prints one line 'K=path ...' (the entry tool's arguments)
and, on stderr, the counts.
  python -I tools/chem_ld2_e22pick.py RUNDIR
"""
import glob
import os
import re
import sys

DATA = os.environ.get("QCC_DATA", os.path.expanduser("~/qcc/data"))
sys.path[:0] = [os.path.join(DATA, "kickmix")]


def main():
    import kickmix as km
    d = sys.argv[1]
    best = {}
    for f in sorted(glob.glob(os.path.join(d, "[onbra]*", "*.kmx")) + glob.glob(os.path.join(d, "n*", "s*", "*.kmx"))):
        sub = os.path.basename(os.path.dirname(f))
        if re.match(r"^s\d+$", sub):          # R59's g42 driver writes n<K>/s1/angle_k_<K>.kmx
            sub = os.path.basename(os.path.dirname(os.path.dirname(f)))
        m = re.match(r"^[onbra](\d+)$", sub)   # a<K> = an owner rebuilt with its annealed last group
        if not m:
            continue
        k = int(m.group(1))
        try:
            c = km.Circuit(open(f).read())
        except Exception as exc:  # a half-written file
            print(f"skip {f}: {exc}", file=sys.stderr)
            continue
        n = int(c.operation_statistics().num_toffoli_like)
        if c.num_qubits > 128:
            continue
        if k not in best or n < best[k][0]:
            best[k] = (n, f)
    for k in sorted(best):
        print(f"k{k}: {best[k][0]} {best[k][1]}", file=sys.stderr)
    print(" ".join(f"{k}={best[k][1]}" for k in sorted(best)))


if __name__ == "__main__":
    main()
