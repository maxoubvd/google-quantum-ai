#!/usr/bin/env python3
"""invsqrt_fc2_patch_norm.py -- FC2 patch "norm" for FC's family-C generator (own code).

Replaces FC's barrel normaliser (26 ANDs + 58 CSWAPs forward, 80 in the mirror) by:
  1. thermometer chain P_j = [k >= j] = [x < 2^(22-2j)], j = 1..10: P_1 = AND(~x21, ~x20), then
     P_j = AND(P_(j-1), u_j), u_j = AND(~x_(23-2j), ~x_(22-2j)) erased at once: 19 ANDs (mirror: 9);
  2. both 11-position pair lanes (even bits, odd bits) shifted up by k with the pruned low-bit-first
     butterfly of tools/invsqrt_fc2_lanes.py (17 pairs per lane, one selector wire per pair = a parity of
     the P_j by CNOTs, shared by the two lanes): 34 CSWAPs;
  3. flags F = binary k as parities of the P_j (CNOTs): F0 = P1+..+P10, F1 = P2+P4+P6+P8+P10, F2 = P4+P8,
     F3 = P8 (same wires and meaning as FC's st_norm: F[i] = bit i of k);
  4. the P_j erased from F (k <= 10, so F3 = 1 forces F2 = 0): t = F1 F0 (1 AND), o = F1 or F0 (CNOTs),
     every P_j = affine part xor AND(two live wires): CNOTs then a measurement erase (0 ANDs).
Forward 54 Toffolis, mirror 52 (rebuild of P from F 9, CSWAPs 34, chain 9): 106 both ways against 164.
Mechanism credit: the zero-supported binary butterfly with carry selectors is the public generator's
(stage94_bit_routing_math.py, Syed Asad Ali, Apache-2.0); this file re-derives it for pair lanes.
Use: import invsqrt_fc2_patch_norm as P; P.apply(F8, cfg)  (F8 = FC's invsqrt_fc_g8 module)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_fc2_lanes as LN  # noqa: E402

NET = [(o["stage"], o["a"], o["b"], o["terms"]) for o in LN.network(11, 10)]
NET10 = [(o["stage"], o["a"], o["b"], o["terms"]) for o in LN.network(10, 10)]   # odd lane without its top slot
_OPT = {"cu": False}


def _gf2_solve(cols, target):
    """coefficients a (list of 0/1 over cols) with XOR_i a_i cols[i] = target (int bit vectors), or None."""
    n = len(cols)
    rows = [(cols[i], 1 << i) for i in range(n)]
    basis = []                                         # (vector, combination)
    for v, cmb in rows:
        for bv, bc in basis:
            if v ^ bv < v:
                v, cmb = v ^ bv, cmb ^ bc
        if v:
            basis.append((v, cmb))
            basis.sort(key=lambda t: -t[0])
    v, cmb = target, 0
    for bv, bc in basis:
        if v ^ bv < v:
            v, cmb = v ^ bv, cmb ^ bc
    if v:
        return None
    return [(cmb >> i) & 1 for i in range(n)]


def _code8_tables():
    """truth tables over k = 0..10 (bit k of an int) of the 8-wire code F0..F3, A = F0 F1, B = F0 F2, C = F1 F2,
    D = F0 F1 F2 (plus the constant), and of the thresholds T_j = [k >= j]."""
    full = (1 << 11) - 1
    Fb = [sum(((k >> i) & 1) << k for k in range(11)) for i in range(4)]
    A, B, C = Fb[0] & Fb[1], Fb[0] & Fb[2], Fb[1] & Fb[2]
    D = A & Fb[2]
    code = Fb + [A, B, C, D, full]                     # index 8 = the constant
    T = [full] + [sum(int(k >= j) << k for k in range(11)) for j in range(1, 11)]
    return code, T


def _code8_undo(c, X, F, debts, ids):
    """FC2 norm_code8 (18:25Z): NORM's explicit undo from (x shifted, F = binary k) back to (x, 0), the public
    line's idea re-derived here (P1's log: 8-wire threshold code, 4 CCX): A = F0F1, B = F0F2, C = F1F2, D = A F2
    (4 CCX) span with F every shift selector and T_3..T_8, T_1 + T_9, T_2 + T_10; un-shift (34 CSWAP); then erase:
    c19 -> T_9 by the promise T_1 = 1 + x21 + x20 (x < 3 2^20), debts paid by Z on T_3, T_5, T_7, T_9, c210 erased
    with 2 conditioned ANDs (T_10 = T_9 u_10, T_2 = T_1 u_2), T_9 .. T_4 with 1 each, T_3 with 3 (T_1 u_2 u_3):
    expected 4 + 34 + 5.5 against 6 + 34 + 4.5. Returns the op list."""
    save = c.ops
    c.ops = []
    code_t, T = _code8_tables()
    W = list(F) + list(ids[:4])                        # wires of F0..F3, A, B, C, D
    s, t1, t2, t3 = ids[4:8]
    c.ccx(F[0], F[1], W[4])
    c.ccx(F[0], F[2], W[5])
    c.ccx(F[1], F[2], W[6])
    c.ccx(W[4], F[2], W[7])

    def form(fn_table, host):
        """XOR into wire `host` (zero) the code combination equal to the table fn_table."""
        a = _gf2_solve(code_t, fn_table)
        assert a is not None, "selector not in the code span"
        ops = []
        for i in range(8):
            if a[i]:
                ops.append(("CX", W[i], host, None))
        if a[8]:
            ops.append(("X", host, None))
        return ops
    for op in reversed(NET):                           # un-shift: the forward's network backwards
        _, a_, b_, terms = op
        tab = 0
        for j in terms:
            tab ^= T[j]
        fo = form(tab, s)
        c.ops.extend(fo)
        c.cswap(s, X[2 * a_ + 1], X[2 * b_ + 1])
        c.cswap(s, X[2 * a_], X[2 * b_])
        c.ops.extend(fo[::-1])
    # code basis -> threshold basis in place: wires hold T3..T8, c19 = T1 + T9, c210 = T2 + T10
    targets = [T[3], T[4], T[5], T[6], T[7], T[8], T[1] ^ T[9], T[2] ^ T[10]]
    M = [_gf2_solve(code_t[:8], tg) for tg in targets]
    assert all(m is not None for m in M)
    rows = [sum(m[i] << i for i in range(8)) for m in M]    # wire r must end as XOR of code wires in rows[r]
    # Gauss-Jordan on rows -> identity, recording row ops; apply the inverse sequence as CNOTs
    R_ = list(rows)
    perm = list(range(8))
    rec = []
    for col in range(8):
        piv = next(r for r in range(col, 8) if (R_[r] >> col) & 1)
        if piv != col:
            R_[piv], R_[col] = R_[col], R_[piv]
            perm[piv], perm[col] = perm[col], perm[piv]
            rec.append(("swap", piv, col))
        for r in range(8):
            if r != col and (R_[r] >> col) & 1:
                R_[r] ^= R_[col]
                rec.append(("add", col, r))            # row r += row col
    # E_k..E_1 rows = I  ->  rows = E_1..E_k ; as wire ops (state g = E b), apply E_k first ... E_1 last
    Wn = list(W)                                       # wire currently holding "row slot" i
    for kind, a_, b_ in reversed(rec):
        if kind == "swap":
            Wn[a_], Wn[b_] = Wn[b_], Wn[a_]
        else:
            c.cx(Wn[a_], Wn[b_])                       # slot b_ += slot a_
    th = {3: Wn[0], 4: Wn[1], 5: Wn[2], 6: Wn[3], 7: Wn[4], 8: Wn[5]}
    c19, c210 = Wn[6], Wn[7]
    # c19 -> T9: T1 = NOT x21 AND NOT x20 = 1 + x21 + x20 on every input (x21 x20 = 0)
    c.x(c19); c.cx(X[21], c19); c.cx(X[20], c19)
    th[9] = c19
    for j in (3, 5, 7, 9):                             # norm_debt's debts (-1)^(m_j T_j)
        if j in debts:
            c.z(th[j], debts[j])

    def u_pair(j):
        return X[23 - 2 * j], X[22 - 2 * j]

    def t1_on():                                       # x20 := T1 = 1 + x21 + x20
        c.cx(X[21], X[20]); c.x(X[20])

    def t1_off():
        c.x(X[20]); c.cx(X[21], X[20])
    # c210 = T2 + T10: (-1)^(m T10) = CZ(T9, u_10), (-1)^(m T2) = CZ(T1, u_2), u's conditioned ANDs
    m = c.bit()
    c.ops.append(("HMR", c210, m))
    mc = m                                             # I39 n25: c210's outcome bit, T2's debt carried to T3
    for (j, ctl) in (((10, th[9]),) if _OPT.get("n25") else ((10, th[9]), (2, None))):
        a_, b_ = u_pair(j)
        c.x(a_); c.x(b_)
        c.cand_compute(a_, b_, t1, m)
        if ctl is None:
            t1_on()
            c.cz(X[20], t1)
            t1_off()
        else:
            c.cz(ctl, t1)
        c.cand_erase(a_, b_, t1, m)
        c.x(b_); c.x(a_)
    for j in (9, 8, 7, 6, 5, 4):                       # T_j = T_(j-1) u_j
        m = c.bit()
        c.ops.append(("HMR", th[j], m))
        a_, b_ = u_pair(j)
        c.x(a_); c.x(b_)
        c.cand_compute(a_, b_, t1, m)
        c.cz(th[j - 1], t1)
        c.cand_erase(a_, b_, t1, m)
        c.x(b_); c.x(a_)
    # T3 = T1 u_2 u_3: t2 = T1 NOT x19, t3 = t2 NOT x18 (= T2), t1 = u_3: CZ(t3, t1)
    m = c.bit()
    c.ops.append(("HMR", th[3], m))
    a2, b2 = u_pair(2)
    a3, b3 = u_pair(3)
    for w in (a2, b2, a3, b3):
        c.x(w)
    t1_on()
    c.cand_compute(X[20], a2, t2, m)
    c.cand_compute(t2, b2, t3, m)
    c.cand_compute(a3, b3, t1, m)
    c.cz(t3, t1)
    if _OPT.get("n25"):                                # I39 n25: (-1)^(mc m3 T2): t3 = m3 T2 is on a wire
        c.z(t3, cond=mc)
    c.cand_erase(a3, b3, t1, m)
    c.cand_erase(t2, b2, t3, m)
    c.cand_erase(X[20], a2, t2, m)
    if _OPT.get("n25"):                                # I39 n25: (-1)^(mc (1 - m3) T2) under d = mc AND NOT m3
        assert hasattr(c, "band"), "n25 needs FC's kit extension (invsqrt_fw_fcpol) installed"
        d = c.band(((mc, 1), (m, 0)))
        c.cand_compute(X[20], a2, t2, d)               # t2 = d T1 NOT x19 (t2 is zero again here)
        c.cz(t2, b2)                                   # (-1)^(d T1 NOT x19 NOT x18) = (-1)^(d T2)
        c.cand_erase(X[20], a2, t2, d)
    t1_off()
    for w in (b3, a3, b2, a2):
        c.x(w)
    out = c.ops
    c.ops = save
    return out


def st_norm2(c, X):
    code8 = _OPT.get("code8") and hasattr(c, "blocks")
    if code8:                                          # FC2 norm_code8: NORM is one block with the explicit undo above
        bid8 = len(c.blocks)
        c.blocks[bid8] = None
        c.ops.append(("BBEG", bid8))
    F = c.alloc(4)
    P = [None] + c.alloc(10)
    top = X[2:22]
    for w in top:
        c.x(w)
    cu = _OPT["cu"] and hasattr(c, "blocks")
    if cu:
        # norm_cu (FC2, 16:30Z): the chain sits in a block whose explicit undo (the mirror's) erases P_j by an
        # X-measurement; its phase (-1)^(m P_j) = (-1)^(m P_(j-1) u_j) is repaid by a CZ with u_j computed only
        # when m = 1 (conditioned AND: 1/2 expected): 4.5 instead of the mirror's 9 recomputed u_j.
        bid = len(c.blocks)
        c.blocks[bid] = None
        c.ops.append(("BBEG", bid))
    if _OPT.get("p1free"):
        # norm_p1free (FC2, 21:10Z): x21 and x20 are never both 1 (MAX_X = 3,133,452 < 3 * 2^20), so with the top
        # inverted P_1 = AND(~x21, ~x20) = NOT(~x21 xor ~x20): two CNOTs and an X, no AND (the public line's P_2)
        c.cx(X[21], P[1])
        c.cx(X[20], P[1])
        c.x(P[1])
    else:
        c.and_compute(X[21], X[20], P[1])
    uw = {}
    odd10 = bool(_OPT.get("odd10")) and code8
    if odd10:
        # norm_odd10 (FC2, 22:15Z): per-bit chain Q_j = P_(j-1) AND NOT x_a (a = 23 - 2j), P_j = Q_j AND NOT x_b
        # (the same 2 ANDs per pair as u_j); the leading one, when it sits on an odd bit, is cleared by CNOTs
        # (e_j = P_(j-1) + Q_j) and its parity kept on Lodd, so the odd lane's top slot is always empty and the odd
        # lane takes a 10-position network (15 CSWAPs, not 17); bit 21 = Lodd after the shift. Q_j erased at once:
        # Q_j + P_j = (P_(j-1) + P_j) AND NOT Lodd (the leading pair is the only one where the first factor is 1,
        # and Lodd already holds its e). The explicit undo (norm_code8) un-shifts with the 11-position network: the
        # same map x -> x 4^k, so it is unchanged.
        Lodd = c.alloc(1)[0]
        c.x(Lodd)                                      # top inverted: X[21] = NOT x21, e_1 = x21
        c.cx(X[21], Lodd)                              # Lodd = x21
        c.cx(Lodd, X[21])                              # NOT x21 xor x21 = 1: x21 cleared (inverted frame)
        for j in range(2, 11):
            a, b = X[23 - 2 * j], X[22 - 2 * j]
            q = c.and_compute(P[j - 1], a)             # Q_j = P_(j-1) NOT x_a
            c.and_compute(q, b, P[j])                  # P_j = Q_j NOT x_b
            c.cx(P[j - 1], Lodd)
            c.cx(q, Lodd)                              # Lodd ^= e_j
            c.cx(P[j - 1], a)
            c.cx(q, a)                                 # x_a ^= e_j (clear)
            c.cx(P[j], q)                              # q = (P_(j-1) + P_j) AND NOT Lodd
            c.cx(P[j], P[j - 1])
            c.x(Lodd)
            c.and_erase(P[j - 1], Lodd, q)
            c.x(Lodd)
            c.cx(P[j], P[j - 1])
        # pair 11 (bits 1, 0; k = 10, x in {1, 2, 3}): Q_11 = P_10 NOT x1 (x1 is outside the inverted top), e_11 =
        # P_10 + Q_11 clears x1 when it leads; P_11 = [x = 0] = 0, so Q_11 = P_10 AND NOT Lodd (erased free)
        c.x(X[1])
        q = c.and_compute(P[10], X[1])
        c.x(X[1])
        c.cx(P[10], Lodd)
        c.cx(q, Lodd)
        c.cx(P[10], X[1])
        c.cx(q, X[1])
        c.x(Lodd)
        c.and_erase(P[10], Lodd, q)
        c.x(Lodd)
    else:
        for j in range(2, 11):
            a, b = X[23 - 2 * j], X[22 - 2 * j]
            u = c.and_compute(a, b)
            uw[j] = u
            c.and_compute(P[j - 1], u, P[j])
            c.and_erase(a, b, u)
    if cu and odd10:                                   # (inside norm_code8's block: this inner undo is never used)
        c.ops.append(("BEND", bid))
        c.blocks[bid] = []
    if cu and not odd10:
        c.ops.append(("BEND", bid))
        undo = []
        for j in range(10, 1, -1):
            a, b = X[23 - 2 * j], X[22 - 2 * j]
            m = c.bit()
            undo += [("HMR", P[j], m), ("CANDC", a, b, uw[j], m), ("CZ", P[j - 1], uw[j], m),
                     ("CANDE", a, b, uw[j], m)]
        if _OPT.get("p1free"):
            undo += [("X", P[1], None), ("CX", X[20], P[1], None), ("CX", X[21], P[1], None)]
        else:
            m = c.bit()
            undo += [("HMR", P[1], m), ("CZ", X[21], X[20], m)]
        c.blocks[bid] = undo
    for w in top:
        c.x(w)
    s = c.alloc(1)[0]

    def sel(terms):
        for j in terms:
            if j == 0:
                c.x(s)
            else:
                c.cx(P[j], s)
    if odd10:
        for (_, a, b, terms) in NET:                   # even lane: 11 positions
            sel(terms)
            c.cswap(s, X[2 * a], X[2 * b])
            sel(terms)
        for (_, a, b, terms) in NET10:                 # odd lane: positions 0..9 (bit 21 is empty here)
            sel(terms)
            c.cswap(s, X[2 * a + 1], X[2 * b + 1])
            sel(terms)
        c.cx(Lodd, X[21])                              # bit 21 = the parity of the leading one's position
        c.cx(X[21], Lodd)                              # Lodd = 0
        c.free_([Lodd])
    else:
        for (_, a, b, terms) in NET:
            sel(terms)
            c.cswap(s, X[2 * a], X[2 * b])
            c.cswap(s, X[2 * a + 1], X[2 * b + 1])
            sel(terms)
    c.free_([s])
    for j in range(1, 11):
        c.cx(P[j], F[0])
    for j in (2, 4, 6, 8, 10):
        c.cx(P[j], F[1])
    for j in (4, 8):
        c.cx(P[j], F[2])
    c.cx(P[8], F[3])
    # erase the thermometer from the flags
    if _OPT.get("debt") and hasattr(c, "blocks"):
        # norm_debt (FC2, 17:05Z): P3, P5, P7, P9 X-measured with phase debts (no AND t forward); the block's
        # explicit undo is the mirror of the erase below (it rebuilds P from F with t, as before) followed by Z
        # repairs of the debts on the rebuilt wires: -1 Toffoli forward
        bid = len(c.blocks)
        c.blocks[bid] = None
        c.ops.append(("BBEG", bid))
        t, o = c.alloc(2)
        c.free_([o, t])                               # ids for the undo's temporaries (free there as here)
        mb = {j: c.hmr(P[j]) for j in (3, 5, 7, 9)}
        c.free_([P[9], P[7], P[5], P[3]])
        c.cx(P[2], P[1]); c.cx(F[0], P[1]); c.and_erase(P[2], F[0], P[1])
        for w in (P[4], F[1], P[6], F[3], P[10]):
            c.cx(w, P[2])
        c.cx(F[2], P[4]); c.cx(F[3], P[4])
        c.cx(F[3], P[6]); c.and_erase(F[2], F[1], P[6])
        c.cx(F[3], P[8])
        c.and_erase(F[3], F[1], P[10])
        c.free_([P[8], P[4], P[2]])
        c.ops.append(("BEND", bid))
        X_ = lambda a, b: ("CX", a, b, None)
        orig = ([("ANDC", F[1], F[0], t), X_(F[1], o), X_(F[0], o), X_(t, o),
                 X_(P[2], P[1]), X_(F[0], P[1]), ("ANDE", P[2], F[0], P[1])]
                + [X_(w, P[2]) for w in (P[4], F[1], P[6], F[3], P[10])]
                + [X_(w, P[3]) for w in (P[4], t, P[7], F[3])]
                + [X_(F[2], P[4]), X_(F[3], P[4]), X_(F[3], P[5]), ("ANDE", F[2], o, P[5]),
                   X_(F[3], P[6]), ("ANDE", F[2], F[1], P[6]), X_(F[3], P[7]), ("ANDE", F[2], t, P[7]),
                   X_(F[3], P[8]), ("ANDE", F[3], o, P[9]), ("ANDE", F[3], F[1], P[10]),
                   X_(t, o), X_(F[0], o), X_(F[1], o), ("ANDE", F[1], F[0], t)])
        c.blocks[bid] = c.mirror(orig) + [("Z", P[j], mb[j]) for j in (3, 5, 7, 9)]
        if code8:
            ids = c.alloc(8)                           # wire ids free at NORM's end (the undo's start)
            c.free_(ids)
            c.ops.append(("BEND", bid8))
            c.blocks[bid8] = _code8_undo(c, X, F, mb, ids)
        return F
    t = c.and_compute(F[1], F[0])
    o = c.alloc(1)[0]
    for w in (F[1], F[0], t):
        c.cx(w, o)                                   # o = F1 or F0
    c.cx(P[2], P[1]); c.cx(F[0], P[1]); c.and_erase(P[2], F[0], P[1])     # P1 = P2 or F0
    # k <= 10: F3 = 1 forces F2 = 0 and F1 F0 = 0, so P4 F1 = F2 F1 + F3 F1 and P4 t = F2 t: CNOTs only
    for w in (P[4], F[1], P[6], F[3], P[10]):                              # P2 = P4 + F1 + (P6 + F3) + P10
        c.cx(w, P[2])
    for w in (P[4], t, P[7], F[3]):                                        # P3 = P4 + t + (P7 + F3)
        c.cx(w, P[3])
    c.cx(F[2], P[4]); c.cx(F[3], P[4])                                    # P4 = F2 + F3
    c.cx(F[3], P[5]); c.and_erase(F[2], o, P[5])                          # P5 = F3 + F2 (F1 or F0)
    c.cx(F[3], P[6]); c.and_erase(F[2], F[1], P[6])                       # P6 = F3 + F2 F1
    c.cx(F[3], P[7]); c.and_erase(F[2], t, P[7])                          # P7 = F3 + F2 F1 F0
    c.cx(F[3], P[8])                                                      # P8 = F3
    c.and_erase(F[3], o, P[9])                                            # P9 = F3 (F1 or F0)
    c.and_erase(F[3], F[1], P[10])                                        # P10 = F3 F1
    for w in (t, F[0], F[1]):
        c.cx(w, o)
    c.free_([P[8], P[4], P[3], P[2], o])
    c.and_erase(F[1], F[0], t)
    return F


def check_identities():
    for k in range(11):
        f = [(k >> i) & 1 for i in range(4)]
        p = [None] + [int(k >= j) for j in range(1, 11)]
        t = f[1] & f[0]
        o = f[1] | f[0]
        assert p[1] == p[2] | f[0] and p[2] == p[4] | f[1] and p[3] == p[4] | t
        assert p[4] == f[2] ^ f[3] and p[5] == f[3] ^ (f[2] & o) and p[6] == f[3] ^ (f[2] & f[1])
        assert p[7] == f[3] ^ (f[2] & t) and p[8] == f[3] and p[9] == f[3] & o and p[10] == f[3] & f[1]
        assert f[0] == sum(p[1:]) % 2 and f[1] == (p[2] + p[4] + p[6] + p[8] + p[10]) % 2
        assert p[2] == p[4] ^ f[1] ^ p[6] ^ f[3] ^ p[10] and p[3] == p[4] ^ t ^ p[7] ^ f[3]
        assert f[2] == (p[4] + p[8]) % 2 and f[3] == p[8]
    return True


def apply(F8, cfg):
    done = F8.__dict__.setdefault("_fc2_done", set())
    _OPT["cu"] = bool(cfg.get("norm_cu"))
    _OPT["debt"] = bool(cfg.get("norm_debt"))
    _OPT["code8"] = bool(cfg.get("norm_code8"))
    _OPT["p1free"] = bool(cfg.get("norm_p1free"))
    _OPT["odd10"] = bool(cfg.get("norm_odd10"))
    _OPT["n25"] = bool(cfg.get("_i39n25"))            # I39 (tools/inv_i39_patch_norm25.py)
    assert not _OPT["code8"] or _OPT["debt"], "norm_code8 is written for the norm_debt forward"
    if "norm" in done:
        return
    check_identities()
    F8.G.st_norm = st_norm2
    done.add("norm")
