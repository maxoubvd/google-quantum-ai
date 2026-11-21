#!/usr/bin/env python3
"""inv_i39_xpair.py -- I39 (builder, second holder): a PEEPHOLE on a .kmx: remove every pair X w ... X w (unconditioned)
on the same wire w when no op between them touches w.  X X = I on that wire and the ops between do not act on it, so
the circuit is the same (state, phase, Toffoli count); a wire that was set to 1 only for such a span is no longer busy
there (the 39-wire hold file: END's exit flips the 40th id to 1 and back around a Z on another wire).
  python3 inv_i39_xpair.py IN.kmx OUT.kmx [OP_LO OP_HI]      (optional: only pairs whose first X lies in [LO, HI])
Prints the removed pairs (op indices).  Safe to run: reads IN, writes OUT."""
import sys


def main():
    inp, outp = sys.argv[1], sys.argv[2]
    lo = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    hi = int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 9
    lines = open(inp).read().split("\n")
    isop = [bool(s.strip()) and not s.startswith(("APPEND", "REGISTER", "#")) and not s.startswith("INFO_TAG")
            for s in lines]
    opl = [i for i, f in enumerate(isop) if f]                 # op index -> line index
    toks = [lines[i].split() for i in opl]
    last_x = {}                                                 # wire -> op index of an open X (no touch since)
    drop = set()
    pairs = []
    for k, t in enumerate(toks):
        name, args = t[0], t[1:]
        wires = [a for a in args if a.startswith("q")]
        if name == "X" and len(args) == 1 and args[0].startswith("q"):
            w = args[0]
            if w in last_x and lo <= last_x[w] <= hi:
                drop.update((last_x[w], k))
                pairs.append((last_x[w], k, w))
                del last_x[w]
                continue
            last_x[w] = k
            continue
        for w in wires:
            last_x.pop(w, None)
    drop_lines = {opl[k] for k in drop}
    keep = [ln for i, ln in enumerate(lines) if i not in drop_lines]
    open(outp, "w").write("\n".join(keep))
    print("removed %d X pairs: %s" % (len(pairs), pairs[:12]))


if __name__ == "__main__":
    main()
