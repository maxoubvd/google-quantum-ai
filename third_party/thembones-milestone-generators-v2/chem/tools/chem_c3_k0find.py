"""Builder C3 (sixth holder): the cheapest k0 of a table among VERIFIED entries (both box_verify JSONs PASS) whose batch 0
words equal the table's. Prints 'toffoli entry_dir' or 'none'.
  python -I chem_c3_k0find.py TABLE.npz ROOT [ROOT ...]      (ROOT = a folder of entry folders, e.g. ~/qcc/ws/chem_ld2/best)
"""
import glob
import json
import os
import sys

import numpy as np

DATA = os.environ.get("QCC_DATA", os.path.expanduser("~/qcc/data"))


def passed(d):
    n = 0
    for f in glob.glob(os.path.join(d, "box_verify_*.json")):
        try:
            if json.load(open(f)).get("chem", {}).get("status") == "PASS":
                n += 1
        except Exception:
            pass
    return n >= 2


def main():
    sys.path[:0] = [os.path.join(DATA, "kickmix")]
    import kickmix as km
    w = np.asarray(np.load(sys.argv[1])["angle_words"], np.int64)[:, 0:5]
    best = None
    for root in sys.argv[2:]:
        for d in glob.glob(os.path.join(root, "*")):
            s = os.path.join(d, "dfthcsos", "submission.npz")
            if not os.path.exists(s) or not passed(d):
                continue
            try:
                if not np.array_equal(np.asarray(np.load(s)["angle_words"], np.int64)[:, 0:5], w):
                    continue
                n = int(km.Circuit.from_file(os.path.join(d, "dfthcsos", "angle_k_0.kmb")).operation_statistics().num_toffoli_like)
            except Exception:
                continue
            if best is None or n < best[0]:
                best = (n, d)
    print(f"{best[0]} {best[1]}" if best else "none")


if __name__ == "__main__":
    main()
