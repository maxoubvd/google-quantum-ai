#!/usr/bin/env python3
"""r63_inv_ghost_early.py -- R63-INV-GHOST (designer): END's entry under the cap, at the op-list level.
In the emitted .kmx the retires of END are a run of HMR ops right before the routing.  An X-measurement commutes
with every op that does not touch its wire (or its outcome bit), so each HMR of the run is moved BACK to just
after the last op that touches its wire.  Nothing else changes: same gates, same counts.
  r63_inv_ghost_early.py IN.kmx OUT.kmx q17,q18,q19      (wire numbers as in the op list: the root bit q_j is wire j)
Prints how far each one moved.  Check the result with r63_inv_ghost_chk.py (check and live).
Safe to run: reads IN.kmx, writes OUT.kmx."""
import re
import sys


def main():
    inp, outp, wires = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
    lines = open(inp).read().split("\n")
    isop = [bool(s.strip()) and not s.startswith(("APPEND", "REGISTER", "#")) and not s.startswith("INFO_TAG")
            for s in lines]
    # the run: the first place where the HMRs of all the named wires stand within 12 lines of each other
    hm = {w: [i for i, s in enumerate(lines) if isop[i] and s.split()[0] == "HMR" and s.split()[1] == w]
          for w in wires}
    best = None
    for i0 in hm[wires[0]]:
        idx = [i0]
        for w in wires[1:]:
            near = [i for i in hm[w] if abs(i - i0) <= 12]
            if not near:
                break
            idx.append(near[0])
        else:
            nxt = "\n".join(lines[max(idx): max(idx) + 12])
            if "CSWAP" in nxt:
                best = idx
                break
    assert best is not None, "the run of retires was not found"
    for w, i in sorted(zip(wires, best), key=lambda t: t[1]):
        tok = lines[i].split()
        bit = [t for t in tok if t.startswith("b")][0]
        pw, pb = re.compile(r"\b%s\b" % w), re.compile(r"\b%s\b" % bit)
        j = i
        moved = 0
        while j > 0:
            prev = lines[j - 1]
            k = j - 1
            if isop[k] and (pw.search(prev) or pb.search(prev) or "CONDITION" in prev):
                break                                    # never across a conditioned block's edge
            lines[j - 1], lines[j] = lines[j], lines[j - 1]
            isop[j - 1], isop[j] = isop[j], isop[j - 1]
            j -= 1
            moved += 1 if isop[j + 1] else 0
        print("HMR %s: moved back over %d ops; now right after: %s" % (w, moved, lines[j - 1].strip()[:60]))
    open(outp, "w").write("\n".join(lines))


if __name__ == "__main__":
    main()
