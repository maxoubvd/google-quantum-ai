#!/usr/bin/env python3
"""r41_inv_main_xc.py -- R41-INV-MAIN (designer): EXACT carries in tail rows 3 (and 4) of a divider (in-memory
wrapper on A.fadd, own code; an ARITHMETIC change: the word moves, the LP must be re-fitted).

FC's line-form tail runs its row k (k >= 3) on d's top bits with carry-in v[kk] (kk = k - 2): an approximation of the
carry out of the positions below the window.  Those positions are not on wires, but their bits are known functions
of v and of the tail's own digits, so the exact carry can be made from live wires:
  row 3:  carry = v1 if c3 == c2 else v0 = v1 xor p,  p = (c2 xor c3)(v0 xor v1)          [P, every (v, r), m 7]
  row 4:  carry = MAJ(c2 xor v1 xor c4 xor 1, v2, c'),  c' = v1 if c4 == c3 else v0       [P, with row 3 exact]
(c_k = the control of tail row k, 1 = add).  One temporary AND for row 3, two more for row 4.
The row's lowest target wire L (the tail's leftover, the same wire in every rotation row) receives the carry in its
sum bit, so the leftover stops being affine in (digits, v).  The ANDs are therefore HELD to the tail's last row,
XORed out of L there (L is affine again and FC's relation clears it as before) and then erased by measurement.
Cost: 1 AND forward and 1 in the mirror per row-3 carry (2 more each way for row 4); one clean wire per AND from its
row to the tail's end (DIV1 / DIV2 have 8; MAIN has 2: there this form is for measuring the margin at a larger cap).
  apply(A, cfg)    env R41_XC  = comma list of stages with the exact row-3 carry (MAIN, DIV1, DIV2; default DIV1)
                   env R41_XC4 = comma list of stages with the exact row-4 carry as well (subset of R41_XC)
                   env R41_XT  = comma list of stages whose tail row at taper m - 2 runs at taper m - 1 (what
                                  "sched2": [.., m - 1, m - 1] does in DIV2; for MAIN, whose taper is one knob)
                   env R41_XTK = comma list STAGE.kk: the tail row of STAGE at taper kk runs at taper kk + 1
                   env R41_XCV0 = comma list of stages whose row 3 takes v0 as its carry-in instead of v1 (no AND:
                                  the other free guess; right in the other half of the cases)
                   env R41_XCA4 = comma list of stages whose row 4 takes the digit guess A = c2 xor v1 xor c4 xor 1
                                  as its carry-in instead of v2 (no AND; one clean wire during the row)
"""
import os

ONE = "1"


def apply(A, cfg):
    stages3 = set(t for t in os.environ.get("R41_XC", "DIV1").split(",") if t)
    stages4 = set(t for t in os.environ.get("R41_XC4", "").split(",") if t)
    stages4m = set(t for t in os.environ.get("R41_XC4M", "").split(",") if t)   # row 4: MAJ(A, v2, v1), one AND
    stagest = set(t for t in os.environ.get("R41_XT", "").split(",") if t)       # row kk = m - 2 -> m - 1 (as sched)
    stepk = set((t.split(".")[0], int(t.split(".")[1])) for t in os.environ.get("R41_XTK", "").split(",") if t)
    stagesv0 = set(t for t in os.environ.get("R41_XCV0", "").split(",") if t)
    stagesa4 = set(t for t in os.environ.get("R41_XCA4", "").split(",") if t)
    m = int(cfg["m"])
    v = list(range(22 - m, 22))                    # v = X[22 - m .. 21] (FC's build)
    v0, v1, v2, v3, vtop = v[0], v[1], v[2], v[3], v[m - 1]
    orig_fadd = A.fadd
    st = {"depth": 0, "S": {}, "log": []}

    def fadd(c, ctl, src, tgt, cin=None, force=None):
        if st["depth"]:
            return orig_fadd(c, ctl, src, tgt, cin, force)
        tag = c.tags[-1][1] if getattr(c, "tags", None) else "?"
        key = (id(c), tag)
        st["depth"] += 1
        try:
            if cin == v0 and len(src) and src[0] == v1:          # tail row 2 (taper 0): remember its digit wire
                st["S"][key] = {"w2": ctl}
                return orig_fadd(c, ctl, src, tgt, cin, force)
            if (tag in stagest and m >= 3 and cin == v[m - 2] and len(src) >= 2 and src[0] == vtop and src[1] == ONE
                    and len(tgt) == 3):                          # the row at taper m - 2 run one step faster
                st["log"].append((tag, "fadd", 0, "row at taper m-2 run at taper m-1"))
                return orig_fadd(c, ctl, [ONE], list(tgt[1:]), vtop, force)
            if cin in v and (tag, v.index(cin)) in stepk and len(tgt) >= 2 and len(src) >= 1:
                st["log"].append((tag, "fadd", 0, "row at taper %d run one step faster" % v.index(cin)))
                return orig_fadd(c, ctl, list(src[1:]), list(tgt[1:]), src[0], force)   # the next bit of d is the carry-in
            S = st["S"].get(key)
            if S is None:
                return orig_fadd(c, ctl, src, tgt, cin, force)
            if cin == v1 and len(src) and src[0] == v2 and tag in stagesv0:                     # row 3, free guess v0
                S["w3"] = ctl
                st["log"].append((tag, "fadd", 3, "row 3 carry-in v0"))
                return orig_fadd(c, ctl, src, tgt, v0, force)
            if cin == v2 and len(src) and src[0] == v3 and tag in stagesa4 and "t2" not in S:   # row 4, free guess A
                w2 = S["w2"]                                     # holds c2 xor v0; ctl holds c4
                a = c.alloc(1)[0]
                c.cx(w2, a)
                c.cx(v0, a)
                c.cx(v1, a)
                c.cx(ctl, a)
                c.x(a)                                           # a = c2 xor v1 xor c4 xor 1
                r = orig_fadd(c, ctl, src, tgt, a, force)
                c.x(a)                                           # (ctl still holds c4: the caller's CX comes later)
                c.cx(ctl, a)
                c.cx(v1, a)
                c.cx(v0, a)
                c.cx(w2, a)
                c.free_([a])
                st["log"].append((tag, "fadd", 4, "row 4 carry-in A"))
                return r
            if cin == v1 and len(src) and src[0] == v2 and tag in stages3 and "p" not in S:     # tail row 3
                w2 = S["w2"]                                     # holds c2 xor v0; ctl holds c3
                p = c.alloc(1)[0]
                c.cx(v0, w2)
                c.cx(ctl, w2)                                    # w2 = c2 xor c3
                c.cx(v0, v1)                                     # v1 wire = v0 xor v1
                c.and_compute(w2, v1, p)
                c.cx(v0, v1)
                c.cx(ctl, w2)
                c.cx(v0, w2)
                c.cx(v1, p)                                      # p = the exact carry v1 xor (c2 xor c3)(v0 xor v1)
                r = orig_fadd(c, ctl, src, tgt, p, force)
                c.cx(v1, p)                                      # back to the product
                S.update(p=p, w3=ctl, L=tgt[0])
                st["log"].append((tag, "fadd", 3, "exact carry row 3"))
                return r
            if (cin == v2 and len(src) and src[0] == v3 and tag in stages4m and "p" in S and "q4" not in S):   # row 4
                w2 = S["w2"]                                     # holds c2 xor v0; ctl holds c4
                q4 = c.alloc(1)[0]
                c.cx(v0, w2)
                c.cx(ctl, w2)
                c.x(w2)                                          # w2 = c2 xor c4 xor 1 = A xor v1
                c.cx(v1, v2)                                     # v2 wire = v2 xor v1
                c.and_compute(w2, v2, q4)
                c.cx(v1, v2)
                c.x(w2)
                c.cx(ctl, w2)
                c.cx(v0, w2)
                c.cx(v1, q4)                                     # q4 = MAJ(A, v2, v1)
                r = orig_fadd(c, ctl, src, tgt, q4, force)
                S.update(q4=q4, w4=ctl)                          # (q4 stays the carry while it is held)
                st["log"].append((tag, "fadd", 4, "row 4 carry MAJ(A, v2, v1)"))
                return r
            if (cin == v2 and len(src) and src[0] == v3 and tag in stages4 and "p" in S and "t2" not in S):   # row 4
                w2, w3 = S["w2"], S["w3"]                        # hold c2 xor v0, c3 xor v0; ctl holds c4
                t1 = c.alloc(1)[0]
                t2 = c.alloc(1)[0]
                c.cx(ctl, w3)
                c.cx(v0, w3)                                     # w3 = c3 xor c4
                c.cx(v0, v1)
                c.and_compute(w3, v1, t1)
                c.cx(v0, v1)
                c.cx(v0, w3)
                c.cx(ctl, w3)
                c.cx(v1, t1)                                     # t1 = c' = v1 xor (c3 xor c4)(v0 xor v1)
                c.cx(v0, w2)
                c.cx(v1, w2)
                c.cx(ctl, w2)
                c.x(w2)                                          # w2 = c2 xor v1 xor c4 xor 1 = A
                c.cx(t1, w2)                                     # A xor c'
                c.cx(t1, v2)                                     # v2 wire = v2 xor c'
                c.and_compute(w2, v2, t2)
                c.cx(t1, v2)
                c.cx(t1, w2)
                c.x(w2)
                c.cx(ctl, w2)
                c.cx(v1, w2)
                c.cx(v0, w2)
                c.cx(t1, t2)                                     # t2 = MAJ(A, v2, c') = the exact carry
                r = orig_fadd(c, ctl, src, tgt, t2, force)
                c.cx(t1, t2)                                     # back to the product
                S.update(t1=t1, t2=t2, w4=ctl)
                st["log"].append((tag, "fadd", 4, "exact carry row 4"))
                return r
            if "p" in S and ((cin == vtop and len(src) and src[0] == ONE) or (cin == ONE and not len(src))
                             or (cin == ONE and len(src) and src[0] is None)):   # the tail's last row: settle
                r = orig_fadd(c, ctl, src, tgt, cin, force)
                w2, w3, p, L = S["w2"], S["w3"], S["p"], S["L"]  # w2, w3 hold c2 xor v0, c3 xor v0
                if "t2" in S:
                    t1, t2, w4 = S["t1"], S["t2"], S["w4"]       # w4 holds c4 xor v0 (after the caller's CX);
                    if w4 == ctl:                                # (m = 3: row 4 is the last row, ctl still c4)
                        c.cx(v0, w4)
                    c.cx(t1, t2)
                    c.cx(t2, L)                                  # L loses row 4's exact carry (affine part stays)
                    c.cx(t1, t2)
                    c.cx(w4, w2)
                    c.cx(v1, w2)
                    c.x(w2)
                    c.cx(t1, w2)                                 # w2 = A xor c'
                    c.cx(t1, v2)
                    c.and_erase(w2, v2, t2)
                    c.cx(t1, v2)
                    c.cx(t1, w2)
                    c.x(w2)
                    c.cx(v1, w2)
                    c.cx(w4, w2)
                    c.cx(v1, t1)                                 # t1 back to its product (c3 xor c4)(v0 xor v1)
                    c.cx(w4, w3)                                 # w3 = c3 xor c4
                    c.cx(v0, v1)
                    c.and_erase(w3, v1, t1)
                    c.cx(v0, v1)
                    c.cx(w4, w3)
                    if w4 == ctl:
                        c.cx(v0, w4)
                if "q4" in S:
                    q4, w4 = S["q4"], S["w4"]                    # w4 holds c4 xor v0
                    c.cx(q4, L)                                  # L loses row 4's carry
                    c.cx(v1, q4)                                 # back to the product (A xor v1)(v2 xor v1)
                    c.cx(w4, w2)
                    c.x(w2)                                      # w2 = c2 xor c4 xor 1
                    c.cx(v1, v2)
                    c.and_erase(w2, v2, q4)
                    c.cx(v1, v2)
                    c.x(w2)
                    c.cx(w4, w2)
                c.cx(p, L)                                       # L loses row 3's product: affine in v again
                c.cx(w3, w2)                                     # w2 = c2 xor c3
                c.cx(v0, v1)
                c.and_erase(w2, v1, p)
                c.cx(v0, v1)
                c.cx(w3, w2)
                del st["S"][key]
                return r
            return orig_fadd(c, ctl, src, tgt, cin, force)
        finally:
            st["depth"] -= 1

    A.fadd = fadd
    return st
