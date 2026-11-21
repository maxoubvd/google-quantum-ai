"""C2 (22nd holder): writes tools/chem_ld2_e22dcr.py = a pinned COPY of R63-CHEM-WIRED's tools/r63_chem_wired_dc.py
(read-only there) with one change: ginv_masks skips a wire bound by a relation (WR7's P-family bit-1 relations make
some wires sums of others: no dual mask exists) instead of asserting; such wires simply get no free bits. The tool's
own check after the rewrite (wires whose tensor differs from the model) still guards every emitted circuit.
  python3 tools/chem_ld2_e22mkdcr.py
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OLD = """        # check: G[q'] . m = [q' == q]
        for b2, q2 in gen:
            assert ((an.G[q2] & m).bit_count() & 1) == (1 if q2 == q else 0)
        out[q] = (b, m)"""
NEW = """        # check: G[q'] . m = [q' == q]; C2's copy (22nd holder): a wire bound by a relation (bit-1 relations of
        # WR7's P family: the wire is a sum of others) has no such mask: it gets no free bits instead of an assert
        if all(((an.G[q2] & m).bit_count() & 1) == (1 if q2 == q else 0) for b2, q2 in gen):
            out[q] = (b, m)"""
HDR = ("# C2 (22nd holder): PINNED COPY of tools/r63_chem_wired_dc.py (R63-CHEM-WIRED; read-only there), made by\n"
       "# tools/chem_ld2_e22mkdcr.py; one change: ginv_masks skips wires bound by a relation instead of asserting.")


def main():
    src = open(os.path.join(HERE, "r63_chem_wired_dc.py")).read()
    assert src.count(OLD) == 1, "the function changed upstream: re-read it"
    lines = src.replace(OLD, NEW).split("\n")
    i = 1 if lines[0].startswith("#!") else 0
    lines.insert(i, HDR)
    out = os.path.join(HERE, "chem_ld2_e22dcr.py")
    open(out, "w").write("\n".join(lines))
    print(f"written {out}: {len(lines)} lines")


if __name__ == "__main__":
    main()
