#!/usr/bin/env python3
"""invsqrt_s2_core.py -- S2 (sweep of OTHER invsqrt families): shared pieces.

- the input chain (S1's compact chain: k, v, d, h, rc, target, tol for all 3,133,452 inputs) + derived delta, N
- the ideal word W* = 2^40 / sqrt(N), N = 64 n = d^2 + delta (units 2^-37 turn at k = 0)
- the endpoint shared by every family: the public circuit's phase gates rebuilt from a 37-bit word and the angle
  LP (S1's Phase class, imported read-only from tools/invsqrt_sw_model.py)
- the generic cost: Q x (Toffolis + (n_T + n_rot x T_rot(m)) / 4)

Every family module produces the 37-bit word p (before the k rotation) in the public endpoint's format:
p = (W + OFF) mod 2^37 with W the family's integer approximation of W*; OFF is measured here on the public word.
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_sw_model as SM  # noqa: E402  (S1, read-only)

MASK37 = (1 << 37) - 1
S1_WD = os.path.expanduser("~/qcc/ws/invsqrt_sw/wd")


def load(wd=S1_WD):
    ch = SM.load_chain(os.path.join(wd, "chain.npz"))
    ch["delta"] = 64 * ch["h"] + 63 - 8 * ch["rc"]
    ch["N"] = ch["d"] * ch["d"] + ch["delta"]
    ch["Wf"] = np.ldexp(1.0, 40) / np.sqrt(ch["N"].astype(np.float64))
    return ch


def phase(ch, wd=S1_WD):
    return SM.Phase(os.path.join(wd, "desc_pub.npz"), ch)


def rotate(p, k):
    return SM.rotate(p & MASK37, k)


def tper(m, n_rot=33):
    return 0.56 * math.log2(n_rot * 2 ** 64 / (math.pi * m)) + 5.3


def cost(q, tof, m, n_rot=33, n_t=1):
    if m is None or m <= 0:
        return float("inf"), float("nan")
    tr = tper(m, n_rot)
    return q * (tof + (n_t + n_rot * tr) / 4), tr


def evaluate_word(ph, ch, p, lp=True):
    """p: 37-bit words before rotation (int64). Returns dict(m_pa, nfail, m_lp, m, tper)."""
    w = rotate(p.astype(np.int64), ch["k"])
    m_pa, nfail, m_lp, _ = ph.evaluate(w, lp=lp)
    m = max(m_pa, m_lp if m_lp is not None else -1e30)
    return dict(m_pa=m_pa, nfail=nfail, m_lp=m_lp, m=m)


try:
    import numba as _nb

    @_nb.njit(cache=True)
    def _divk(U, d, kk, p_out, ovf, uw, tail):
        """the shipped divider recurrence (rows 9..35, truncation kk for rows 25..35) on a given 15-bit U: word p
        (before rotation), overflow of the 14-bit remainder, remainder wires u, signs of rows 25..36 (bit j-25)."""
        n = U.shape[0]
        for i in range(n):
            Ui = U[i]
            bad = Ui < -16384 or Ui >= 16384
            if Ui < -16384:
                Ui = -16384
            if Ui >= 16384:
                Ui = 16383
            pb = Ui & 32767
            di = d[i]
            V = di // 2
            rr = 4096 - di
            phase = -((1 << 37) - (1 << 29))
            tb = 0
            for j in range(9, 36):
                s = 1 if rr < 0 else 0
                if j >= 25:
                    tb |= s << (j - 25)
                phase -= s << (37 - j)
                if j <= 23:
                    bit = (pb >> (23 - j)) & 1
                    if j == 9:
                        bit = -bit
                    raw = 2 * rr + (2 * s - 1) * di + bit
                elif j == 24:
                    raw = 2 * rr + (2 * s - 1) * di
                else:
                    t = kk[j - 25]
                    if t:
                        A = (((V >> (t + 1)) + ((V >> t) & 1)) << (t + 1)) - (V & 1)
                    else:
                        A = V
                    u = (rr - 1) // 2
                    raw = 2 * (2 * u + s + (2 * s - 1) * A) + 1
                nextu = raw // 2
                if nextu < -8192 or nextu >= 8192:
                    bad = True
                nextu = ((nextu + 8192) & 16383) - 8192
                rr = 2 * nextu + (raw & 1)
            s36 = 1 if rr < 0 else 0
            tb |= s36 << 11
            phase -= s36 << 1
            p_out[i] = (phase - 1) & ((1 << 37) - 1)
            ovf[i] = bad
            uw[i] = ((rr - 1) // 2) & 16383
            tail[i] = tb

    HAVE_NUMBA = True
except Exception:
    HAVE_NUMBA = False


def divide_U(ch, U, kk):
    """shipped divider on a given U: (p, overflow, S rows 25..36, u). numba when available."""
    if HAVE_NUMBA:
        n = len(U)
        p = np.empty(n, dtype=np.int64)
        ovf = np.empty(n, dtype=np.bool_)
        uw = np.empty(n, dtype=np.int64)
        tail = np.empty(n, dtype=np.int64)
        _divk(U.astype(np.int64), ch["d"].astype(np.int64), np.array(kk, dtype=np.int64), p, ovf, uw, tail)
        S = {j: (tail >> (j - 25)) & 1 for j in range(25, 37)}
        return p, ovf, S, uw
    bad = (U < -16384) | (U >= 16384)
    Uc = np.clip(U, -16384, 16383)
    p, ovf = SM.divide(ch["d"], Uc, tuple(kk))
    S, u = SM.divide_signs_u(ch["d"], Uc, tuple(kk))
    return p, ovf | bad, S, u


def pub_word(ch, D, pmask=0, kk=SM.PUB_KK, seed_off=-1):
    """the shipped family's word before rotation (S1's model) and its overflow mask."""
    v = ch["v"]
    C = 2 * D[v]
    U0 = SM.seed_from_C(C, ch["rc"], seed_off)
    U = SM.product(C, U0, ch["h"] + 512, pmask)
    bad = (U < -16384) | (U >= 16384)
    p, ovf = SM.divide(ch["d"], np.clip(U, -16384, 16383), kk)
    return p, ovf | bad


if __name__ == "__main__":
    # word format check: public p against W*
    ch = load()
    D = np.load(os.path.join(S1_WD, "tables", "pub.npy")).astype(np.int64)
    p, ovf = pub_word(ch, D)
    Wf = ch["Wf"]
    diff = p.astype(np.float64) - Wf
    off = np.round(np.median(diff))
    e = diff - off
    print(f"overflow {int(ovf.sum())}; p - W*: median {np.median(diff):.3f}, offset {off:.0f}; "
          f"residual min {e.min():.3f} max {e.max():.3f}")
    tolw = ch["tol"].astype(np.float64) / np.ldexp(1.0, 27) / np.exp2(ch["k"].astype(np.float64))
    print(f"tolerance in word units: min {tolw.min():.3f} max {tolw.max():.3f}; worst |e|/tol "
          f"{np.max(np.abs(e) / tolw):.4f}")
    print(f"p range {p.min()}..{p.max()}; W* range {Wf.min():.1f}..{Wf.max():.1f}")
    # self-check of the numba divider against S1's numpy recurrence (public U, two schedules)
    C = 2 * D[ch["v"]]
    U = SM.product(C, SM.seed_from_C(C, ch["rc"], -1), ch["h"] + 512, 0)
    for kk in (SM.PUB_KK, (0, 0, 1, 2, 4, 5, 6, 7, 8, 9, 10)):
        p1, o1, S1, u1 = divide_U(ch, U, kk)
        p2, o2 = SM.divide(ch["d"], np.clip(U, -16384, 16383), kk)
        S2, u2 = SM.divide_signs_u(ch["d"], np.clip(U, -16384, 16383), kk)
        ok = np.array_equal(p1, p2) and np.array_equal(u1, u2) and all(np.array_equal(S1[j], S2[j])
                                                                        for j in range(25, 37))
        print(f"numba divider == S1 numpy divider for kk {kk}: {ok}")
