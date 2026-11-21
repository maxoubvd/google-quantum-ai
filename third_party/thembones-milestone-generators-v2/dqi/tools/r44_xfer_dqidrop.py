#!/usr/bin/env python3
"""r44_xfer_dqidrop.py -- R44-XFER: the DROP list of a DQI decoder (ops that never fire on valid inputs) and the
file without them.

  sel  OUT.npy ONES1.npz [ONES2.npz ...] [--max-ones 0] [--min-samples 0]
       combine the counts of tools/r44_xfer_dqiones.py; the class = static Toffoli-class ops that execute (act > 0)
       and whose product was 1 on at most --max-ones samples in all; writes their static indices k; prints the class
       by kind and by op-index bin, and the executed mass of the ops by firing rate.
  build IN.kmb CAND.npy OUT.kmb
       delete the listed ops (B5's editor tools/dqi_int_kmbw.py); nothing else changes.
  rate CAND.npy ONES_FRESH.npz ...     the class on FRESH samples: ops of the class that fired, their executed mass.
  sel2 OUT.npz IMP1.npz [IMP2.npz ...] [--ones ONES1.npz,ONES2.npz] [--minpos 30]
       combine the counts of tools/r44_xfer_dqiimp.py (and, for the never-1 class, of r44_xfer_dqiones.py); classes:
       mode 1 product never 1 (delete), 2 product == a (CCX: CX a t; CSWAP: CX a b, CX a t), 3 product == x
       (CCX: CX b t; CSWAP: SWAP b t), 4 product == 1 + a + x; modes 2-4 need 0 violations and >= minpos samples
       with product 1.  Writes k, mode.
  build2 IN.kmb CAND.npz OUT.kmb     the rewrite (ops with their own condition keep it on the inserted gates).
Price of a deleted class: (probability that any of its ops fires on a shot) x 102,400 shots x 0.465% of the cost.
"""
import sys

import numpy as np


def load(files):
    act = ones = None
    n = 0
    for f in files:
        d = np.load(f)
        if act is None:
            act = d["act"].astype(np.int64).copy()
            ones = d["ones"].astype(np.int64).copy()
            t_op, base = d["t_op"], d["base"]
        else:
            assert (d["t_op"] == t_op).all()
            act += d["act"]
            ones += d["ones"]
        n += int(d["samples"])
    return t_op, base, act, ones, n


def main():
    mode = sys.argv[1]
    if mode == "sel":
        out = sys.argv[2]
        args = sys.argv[3:]
        mx = 0
        if "--max-ones" in args:
            i = args.index("--max-ones")
            mx = int(args[i + 1])
            del args[i:i + 2]
        t_op, base, act, ones, n = load(args)
        ex = act / n
        tot = ex.sum()
        cls = (act > 0) & (ones <= mx)
        print("samples %d; static Toffoli-class %d; executed per shot %.1f" % (n, len(t_op), tot))
        print("CLASS (executes, product 1 on <= %d samples): %d ops, %.1f executed per shot = %.3f%% of the count"
              % (mx, int(cls.sum()), ex[cls].sum(), 100 * ex[cls].sum() / tot))
        for b, nm in ((0x31, "CCX"), (0x32, "CSWAP"), (0x30, "CCZ")):
            s = cls & (base == b)
            print("  %-5s %7d ops %9.1f per shot" % (nm, int(s.sum()), ex[s].sum()))
        edges = np.arange(0, int(t_op.max()) + 5_000_000, 5_000_000)
        h = np.histogram(t_op[cls], bins=edges, weights=ex[cls])[0]
        ha = np.histogram(t_op, bins=edges, weights=ex)[0]
        print("  by op index (5M bins), class / all executed per shot:")
        print("   " + " ".join("%d-%dM:%.0f/%.0f" % (edges[i] // 1_000_000, edges[i + 1] // 1_000_000, h[i], ha[i])
                               for i in range(len(h)) if ha[i] > 0))
        rate = np.where(act > 0, ones / np.maximum(act, 1), 0.0)
        print("  executed per shot by firing rate of the op (all ops):")
        lo = 0.0
        for hi in (0.0, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 1e-2, 1e-1, 1.0):
            s = (act > 0) & (rate <= hi) & ((rate > lo) | (hi == 0.0))
            print("   rate %s%g: %7d ops %10.1f per shot" % ("= " if hi == 0 else "<= ", hi, int(s.sum()), ex[s].sum()))
            lo = hi
        np.save(out, np.flatnonzero(cls))
        print("written %s (%d static indices)" % (out, int(cls.sum())))
    elif mode == "rate":
        cand = np.load(sys.argv[2])
        t_op, base, act, ones, n = load(sys.argv[3:])
        ex = act / n
        f = ones[cand] > 0
        print("fresh samples %d: class %d ops, %.1f executed per shot; ops of the class that fired %d (%.1f per shot); "
              "product-1 samples summed over the class %d"
              % (n, len(cand), ex[cand].sum(), int(f.sum()), ex[cand][f].sum(), int(ones[cand].sum())))
        fired = sum(int(np.load(x)["fired"]) for x in sys.argv[3:])
        print("samples with ANY op of the candidate list firing: %d of %d = %.3g per shot -> expected failed shots "
              "per 102,400: %.2f -> cost x %.4f" % (fired, n, fired / n, 102400.0 * fired / n,
                                                     np.exp(0.00465 * 102400.0 * fired / n)))
    elif mode == "sel2":
        out = sys.argv[2]
        args = sys.argv[3:]
        minpos = 30
        onesf = []
        if "--minpos" in args:
            i = args.index("--minpos")
            minpos = int(args[i + 1])
            del args[i:i + 2]
        if "--ones" in args:
            i = args.index("--ones")
            onesf = args[i + 1].split(",")
            del args[i:i + 2]
        cnt = None
        n = 0
        for f in args:
            d = np.load(f)
            if cnt is None:
                cnt = d["cnt"].astype(np.int64).copy()
                t_op, base, raw = d["t_op"], d["base"], d["raw"]
            else:
                assert (d["t_op"] == t_op).all()
                cnt += d["cnt"]
            n += int(d["samples"])
        act, ones = cnt[:, 0], cnt[:, 1].copy()
        ex = act / n
        tot = ex.sum()
        n1 = n
        if onesf:
            _, _, act1, ones1, m1 = load(onesf)
            ones_all = ones + ones1
            n1 = n + m1
        else:
            ones_all = ones
        md = np.zeros(len(t_op), np.int64)
        tw = (base == 0x31) | (base == 0x32)
        md[(act > 0) & (ones_all == 0)] = 1
        A = tw & (md == 0) & (cnt[:, 2] == 0) & (ones >= minpos)
        md[A] = 2
        B = tw & (md == 0) & (cnt[:, 3] == 0) & (ones >= minpos)
        md[B] = 3
        N = tw & (md == 0) & (cnt[:, 4] == 0) & (ones >= minpos)
        md[N] = 4
        print("samples %d (never-1 class on %d); static %d; executed per shot %.1f" % (n, n1, len(t_op), tot))
        for mm, nm in ((1, "product never 1 (delete)"), (2, "product == a"), (3, "product == x"),
                       (4, "product == 1 + a + x")):
            s_ = md == mm
            print("  mode %d %-26s %6d ops %9.1f executed per shot (CCX %.1f, CSWAP %.1f, CCZ %.1f)"
                  % (mm, nm, int(s_.sum()), ex[s_].sum(), ex[s_ & (base == 0x31)].sum(), ex[s_ & (base == 0x32)].sum(),
                     ex[s_ & (base == 0x30)].sum()))
        sel = md > 0
        # a CSWAP that becomes a plain SWAP or CX pair, a CCX that becomes CX: all Toffoli-free
        print("  ALL: %d ops, %.1f executed per shot = %.3f%% of the count; with their own condition bit: %d ops"
              % (int(sel.sum()), ex[sel].sum(), 100 * ex[sel].sum() / tot, int((sel & ((raw & 0x80) != 0)).sum())))
        edges = np.arange(0, int(t_op.max()) + 5_000_000, 5_000_000)
        h = np.histogram(t_op[sel], bins=edges, weights=ex[sel])[0]
        print("  by op index (5M bins): " + " ".join("%d-%dM:%.0f" % (edges[i] // 1_000_000, edges[i + 1] // 1_000_000, h[i])
                                                      for i in range(len(h)) if h[i] > 0))
        # weaker classes, for the record: implications with fewer positives
        for lo, hi in ((1, 5), (5, 30)):
            wk = tw & (md == 0) & ((cnt[:, 2] == 0) | (cnt[:, 3] == 0) | (cnt[:, 4] == 0)) & (ones >= lo) & (ones < hi)
            print("  not taken: 0 violations with %d..%d positives: %d ops, %.1f per shot" % (lo, hi - 1, int(wk.sum()),
                                                                                           ex[wk].sum()))
        np.savez(out, k=np.flatnonzero(sel), mode=md[sel])
        print("written %s" % out)
    elif mode == "build2":
        import os
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import kmb
        import dqi_int_kmbw as kw
        c = kmb.read(sys.argv[2])
        d = np.load(sys.argv[3])
        ops = np.asarray(c.ops)
        toff = ((ops & 0x7C) == 0x30) & (ops != 0x81)
        t_op = np.flatnonzero(toff)
        i3 = np.cumsum(((ops & 0x30) >> 4) == 3) - 1
        ic = np.cumsum((ops & 0x80) != 0) - 1
        q3 = [np.asarray(x) for x in c.q3]
        cond = np.asarray(c.cond)
        delete, inserts = [], {}
        nd = 0
        for k, mm in zip(d["k"].tolist(), d["mode"].tolist()):
            p = int(t_op[k])
            op = int(ops[p])
            base = op & 0x7F
            a, b, t = (int(q3[j][i3[p]]) for j in range(3))
            cf = op & 0x80
            cd = int(cond[ic[p]]) if cf else -1
            g = []
            if mm == 1:
                pass
            elif base == 0x31:
                if mm == 2:
                    g = [(0x21, (a, t))]
                elif mm == 3:
                    g = [(0x21, (b, t))]
                elif mm == 4:
                    g = [(0x11, (t,)), (0x21, (a, t)), (0x21, (b, t))]
            elif base == 0x32:
                if mm == 2:
                    g = [(0x21, (a, b)), (0x21, (a, t))]
                elif mm == 3:
                    g = [(0x22, (b, t))]
                elif mm == 4:
                    g = [(0x22, (b, t)), (0x11, (b,)), (0x11, (t,)), (0x21, (a, b)), (0x21, (a, t))]
            else:
                continue
            delete.append(p)
            nd += 1
            if g:
                inserts[p] = [(code | cf, qs, -1, cd) for code, qs in g]
        k1 = kw.edit(c, sorted(delete), inserts)
        kw.write(sys.argv[4], k1)
        print("written %s: %d Toffoli-class ops rewritten (%d with inserted Clifford gates) of %d static"
              % (sys.argv[4], nd, len(inserts), len(t_op)))
    elif mode == "build":
        import os
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import kmb
        import dqi_int_kmbw as kw
        c = kmb.read(sys.argv[2])
        cand = np.load(sys.argv[3])
        ops = np.asarray(c.ops)
        toff = ((ops & 0x7C) == 0x30) & (ops != 0x81)
        t_op = np.flatnonzero(toff)
        delete = sorted(int(x) for x in t_op[cand])
        k1 = kw.edit(c, delete, {})
        kw.write(sys.argv[4], k1)
        print("written %s: %d Toffoli-class ops deleted of %d static" % (sys.argv[4], len(delete), len(t_op)))


if __name__ == "__main__":
    main()
