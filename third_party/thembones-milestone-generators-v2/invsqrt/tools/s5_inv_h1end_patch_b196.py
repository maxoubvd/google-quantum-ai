#!/usr/bin/env python3
"""s5_inv_h1end_patch_b196.py -- S5-INV-H1END: install the b196-stretch change (H1's word un-add with q11 LENT: q11
X-measured with a Clifford repair, the un-add re-emitted at one more clean wire as FC's game, q11 re-made by one AND)
on ANY file whose H1 undo is fc19_b312's (e.g. the base after an integrator installs a find elsewhere).

Method (file level, no generator): the change is the window where the built file NEW differs from BASE (fc19_b312),
located with measurement-bit names blanked (longest common prefix / suffix). In TARGET the bit-blanked op sequence of
BASE's window must occur exactly once; it is replaced by NEW's window with bits renamed:
  - a bit NEW's window reads that was set before the window: NEW's prefix is aligned with BASE's op by op (same ops,
    bit names only differ), then BASE's window with TARGET's (aligned op by op);
  - a bit first set inside NEW's window: a fresh number above TARGET's largest bit;
  - a bit set inside a window and read after it: refused (no such bit in this change; checked).
Writes OUT/invsqrt.kmx, OUT/invsqrt.kmb, OUT/patch.json. Then: tools/r63_inv_ghost_chk.py check TARGET.kmx OUT/invsqrt.kmx 0
(all inputs, phase equal to the target's) and the official verifier on both salts.
  python -I s5_inv_h1end_patch_b196.py TARGET.(kmb|kmx) OUT_DIR [BASE.kmx] [NEW.kmx]
Defaults: BASE = ~/qcc/ws/s5/inv_h1end/base/file.kmx (fc19_b312), NEW = ~/qcc/ws/s5/inv_h1end/best/b196lend/invsqrt.kmx."""
import json
import os
import re
import sys

BIT = re.compile(r"\bb(\d+)\b")
KM = "<HOME>/qcc/data/kickmix"


def read_text(path):
    path = os.path.expanduser(path)
    if path.endswith(".kmb"):
        sys.path.insert(0, KM)
        import kickmix as km
        return km.Circuit.from_file(path, format="kmb").to_text() if hasattr(km.Circuit, "to_text") else _kmb_text(km, path)
    return open(path).read()


def _kmb_text(km, path):
    import tempfile
    c = km.Circuit.from_file(path, format="kmb")
    fd, tmp = tempfile.mkstemp(suffix=".kmx")
    os.close(fd)
    c.to_file(tmp, format="kmx")
    t = open(tmp).read()
    os.remove(tmp)
    return t


def lines_of(text):
    return [l.strip() for l in text.splitlines() if l.strip()]


def blank(ls):
    return [BIT.sub("b", l) for l in ls]


def bits(l):
    return [int(m.group(1)) for m in BIT.finditer(l)]


def align(la, lb):
    """bit map from aligned line lists (same ops, bit names may differ): {bit in a: bit in b}; consistent or error."""
    m, r = {}, {}
    for x, y in zip(la, lb):
        bx, by = bits(x), bits(y)
        assert len(bx) == len(by), (x, y)
        for u, v in zip(bx, by):
            if m.setdefault(u, v) != v or r.setdefault(v, u) != u:
                raise SystemExit("inconsistent bit alignment at %r / %r" % (x, y))
    return m


def window(a, b):
    p = 0
    while p < min(len(a), len(b)) and a[p] == b[p]:
        p += 1
    s = 0
    while s < min(len(a), len(b)) - p and a[len(a) - 1 - s] == b[len(b) - 1 - s]:
        s += 1
    return p, s


def main():
    tgt_path, out = sys.argv[1], os.path.expanduser(sys.argv[2])
    base_path = sys.argv[3] if len(sys.argv) > 3 else "~/qcc/ws/s5/inv_h1end/base/file.kmx"
    new_path = sys.argv[4] if len(sys.argv) > 4 else "~/qcc/ws/s5/inv_h1end/best/b196lend/invsqrt.kmx"
    B, N, T = (lines_of(read_text(x)) for x in (base_path, new_path, tgt_path))
    bB, bN, bT = blank(B), blank(N), blank(T)
    p, s = window(bB, bN)
    wB, wN = B[p:len(B) - s], N[p:len(N) - s]
    print("change window: BASE [%d, %d) %d lines -> NEW [%d, %d) %d lines" % (p, len(B) - s, len(wB), p, len(N) - s, len(wN)))
    # the window of BASE must occur once in TARGET (bit-blanked)
    key = bB[p:len(B) - s]
    hits = [i for i in range(len(bT) - len(key) + 1) if bT[i:i + len(key)] == key]
    if len(hits) != 1:
        raise SystemExit("BASE's window occurs %d times in TARGET (need exactly 1)" % len(hits))
    t0 = hits[0]
    wT = T[t0:t0 + len(key)]
    # bits set before the window (prefix), inside, read after
    def first_set(ls):
        seen = set()
        for l in ls:
            seen.update(bits(l))
        return seen
    pre_N = first_set(N[:p])
    in_N = first_set(wN)
    suf_N = first_set(N[len(N) - s:])
    in_B = first_set(wB)
    pre_B = first_set(B[:p])
    suf_B = first_set(B[len(B) - s:])
    leak_N = (in_N - pre_N) & suf_N                # set inside NEW's window, read after it
    mNB = align(N[:p], B[:p])                      # NEW prefix bits -> BASE prefix bits
    mNBs = align(N[len(N) - s:], B[len(B) - s:])   # NEW suffix bits -> BASE suffix bits
    mBT = align(wB, wT)                            # BASE window bits -> TARGET window bits
    tmax = max((max(bits(l)) for l in T if bits(l)), default=0)
    fresh = {}
    ren = {}
    for u in sorted(in_N):
        if u in pre_N:                             # external: set before the window
            v = mNB.get(u)
            if v is None or v not in mBT:
                raise SystemExit("external bit b%d of NEW's window has no counterpart in BASE's window" % u)
            ren[u] = mBT[v]
        elif u in leak_N:                          # set inside, read after: the suffix's name for it, through BASE
            v = mNBs.get(u)
            if v is None or v not in mBT or v in pre_B:
                raise SystemExit("bit b%d set in NEW's window and read after it has no BASE-window counterpart" % u)
            ren[u] = mBT[v]
        else:
            fresh[u] = tmax + 1 + len(fresh)
            ren[u] = fresh[u]
    wNew = [BIT.sub(lambda m_: "b%d" % ren[int(m_.group(1))], l) for l in wN]
    outl = T[:t0] + wNew + T[t0 + len(key):]
    os.makedirs(out, exist_ok=True)
    open(os.path.join(out, "invsqrt.kmx"), "w").write("\n".join(outl) + "\n")
    sys.path.insert(0, KM)
    import kickmix as km
    c = km.Circuit("\n".join(outl))
    c.to_file(os.path.join(out, "invsqrt.kmb"), format="kmb")
    json.dump(dict(target=tgt_path, base=base_path, new=new_path, window_base=[p, len(B) - s], window_new=[p, len(N) - s],
                   window_target=[t0, t0 + len(key)], external={str(k): v for k, v in ren.items() if k in pre_N},
                   fresh={str(k): v for k, v in fresh.items()}), open(os.path.join(out, "patch.json"), "w"), indent=1)
    print("TARGET window [%d, %d) replaced (%d -> %d lines); %d external bits mapped, %d fresh bits from b%d; wrote %s" % (
        t0, t0 + len(key), len(key), len(wNew), len(ren) - len(fresh), len(fresh), tmax + 1, out))


if __name__ == "__main__":
    main()
