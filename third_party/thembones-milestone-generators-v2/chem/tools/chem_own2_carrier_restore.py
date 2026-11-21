"""chem_own2_carrier_restore: put the carriers on all 75 columns of an 8-bit table by exact gauge moves (chem_own2;
rank 7 of research/public_frontier_1113.md). CPU, data only.

    python chem_own2_carrier_restore.py IN.npz OUT.npz

For every column whose bit7 XOR bit6 does not follow its slot's low monomial (the monomial every exact column of
that slot follows), apply the loader kit's exact move `chem_loaders_gauge.flip` on exactly the rows where it is
wrong (half turn on that edge, negation of the edge into its parent and of its earlier siblings). The decoded
vectors do not change (checked with `chem_loaders_gauge.check`, up to the sign of each row), so error_2 and error_1
do not change; the bit-0 codes survive (negation never changes bit 0), the XOR-of-five codes on bits 1 and 2 break
where odd words were negated: re-impose them afterwards with `tools/chem_lens_project.py` (it keeps bits 6-7).
OUT.npz = IN.npz with `angle_words` replaced.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_gauge as gauge  # noqa: E402  (loader kit, read-only import)

PAIRS = (3, 5, 6, 9, 10, 12)


def slot_monomials(words):
    rows = words.shape[0]
    low = np.arange(rows) % 16
    top = ((words >> 7) ^ (words >> 6)) & 1
    exact, slot = {}, {}
    for c in range(words.shape[1]):
        for m in PAIRS:
            if np.all(top[:, c] == ((low & m) == m)):
                exact[c] = m
                break
    for s in range(5):
        ms = sorted({exact[c] for c in range(s, words.shape[1], 5) if c in exact})
        assert len(ms) == 1, f"slot {s}: monomials {ms}"
        slot[s] = ms[0]
    return exact, slot


def main():
    src, dst = sys.argv[1], sys.argv[2]
    with np.load(src, allow_pickle=False) as z:
        arrays = {k: np.array(z[k]) for k in z.files}
    w0 = arrays["angle_words"].astype(np.int64)
    edges = arrays["tree_edges"].astype(np.int64)
    exact, slot = slot_monomials(w0)
    missing = [c for c in range(w0.shape[1]) if c not in exact]
    print(f"slot monomials {slot}; columns without an exact carrier: {missing}")
    low = np.arange(w0.shape[0]) % 16
    w = w0.copy()
    for g in missing:
        m = slot[g % 5]
        top = ((w[:, g] >> 7) ^ (w[:, g] >> 6)) & 1
        wrong = top != ((low & m) == m)
        before = w.copy()
        w = np.asarray(gauge.flip(w, edges, g, wrong)).astype(np.int64)
        changed = sorted(set(np.flatnonzero(np.any(w != before, axis=0)).tolist()) - {g})
        print(f"  col {g:2d}: {int(wrong.sum())} rows flipped; other columns changed {changed}")
    # a negated word of value 64 or 192 (lowest set bit = bit 6) toggles its own carrier: repeat on every column
    # until every carrier is exact (each pass flips only the rows still wrong)
    for it in range(20):
        top = ((w >> 7) ^ (w >> 6)) & 1
        bad = [(c, top[:, c] != ((low & slot[c % 5]) == slot[c % 5])) for c in range(w.shape[1])]
        bad = [(c, r) for c, r in bad if r.any()]
        if not bad:
            break
        print(f"  pass {it + 1}: {len(bad)} columns with {sum(int(r.sum()) for _, r in bad)} wrong rows: "
              + ", ".join(f"{c}:{int(r.sum())}" for c, r in bad[:12]))
        for c, r in bad:
            w = np.asarray(gauge.flip(w, edges, c, r)).astype(np.int64)
    dev = gauge.check(w0, w, edges)
    top = ((w >> 7) ^ (w >> 6)) & 1
    ok = all(np.all(top[:, c] == ((low & slot[c % 5]) == slot[c % 5])) for c in range(w.shape[1]))
    bit0 = int(np.sum((w ^ w0) & 1))
    print(f"decoded vectors: largest deviation {dev:.3e} (up to each row's sign); carriers exact on all 75 columns: {ok}; "
          f"bit-0 changes {bit0}; words changed {int(np.sum(w != w0))} of {w.size}")
    assert dev < 1e-9 and ok and bit0 == 0
    arrays["angle_words"] = w.astype(arrays["angle_words"].dtype)
    np.savez(dst, **arrays)
    print(f"wrote {dst}")


if __name__ == "__main__":
    main()
