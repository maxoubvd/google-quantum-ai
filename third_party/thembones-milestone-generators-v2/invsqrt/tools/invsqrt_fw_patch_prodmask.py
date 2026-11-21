# Copyright 2026 Syed Asad Ali. SPDX-License-Identifier: Apache-2.0. Credits: Opus and Astra.
# _product_carry / _clear_product_carry (candidate01/lifecycle.py) and product_inverse (patches/patch_seedinv.py) below
# are copies of the public invsqrt source (Apache-2.0) with the changes marked F3, by builder F3 of our fork
# (invsqrt_fw, 2026-11-13).
"""F3 patch prodmask: ANY set of product rows 0..7 rounded plainly (MASK bit i = row i plain).

Public row i (s = 8 - i) adds C[s+1:] with the multiplexed carry c_i = (q_i == p) ? C[s-1] : C[s], p = C[8], and
row 9's aliased carry (q_9 == p): the ten rows then sum exactly to floor(F / 2). Plain row: carry-in C[s] (data wire,
no CCX, no owned wire). Parity rule for any mask (F3 derivation, 2026-11-13): write a_m = C[8 - m] (a_8 = C[0] = 0);
the multiplexed rows form runs c = [L_c, U_c] between plain rows, X_c = a_{L_c} ^ a_{U_c + 1}. Run c uses
    p_c = NOT a_{L_c} ^ X_1 ^ ... ^ X_{c-1}        (applied to q by CNOTs and one X: no wire, no Toffoli)
and row 9 uses parity 0 (its carry is NOT q). Each run then telescopes to a residual r_c = s_c X_c with s_c chosen so
the running sum stays in {-1, 0}, and the rows sum exactly to floor(F' / 2), F' = the public F with each plain row's
term replaced by 2 A_j, A_j = (C >> (s+1)) + C[s]: U = U0 - F' // 2 (the sweep's coarse-row model) for every mask.
For MASK = rows 0..N-1 this is the same arithmetic as patch prodplain (other parity polarity, same sums).
Place right AFTER patch_seedinv.py.
"""

MASK = 0xFF   # bit i set: product row i (0..7) rounded plainly


def _runs():
    """F3: runs of multiplexed rows among 0..7 and each run's parity as (C indices, constant)."""
    rows = [i for i in range(8) if not (MASK >> i) & 1]
    runs = []
    for i in rows:
        if runs and runs[-1][1] == i - 1:
            runs[-1][1] = i
        else:
            runs.append([i, i])
    par = {}
    acc = set()                                   # C indices of X_1 ^ ... ^ X_{c-1}
    for L, U in runs:
        p = set(acc) ^ {8 - L}                    # a_L = C[8 - L]
        for i in range(L, U + 1):
            par[i] = (sorted(x for x in p if x != 0), 1)      # NOT a_L ^ ...: constant 1
        x = {8 - L} ^ {8 - (U + 1)}               # X_c = a_L ^ a_{U+1}
        acc ^= {v for v in x if v != 0}
    return par


def apply(ctx):
    bl = ctx.baseline
    lc = ctx.components['lifecycle']
    km = bl.km
    PAR = _runs()

    def parity_on(builder, C, q, i):
        idx, const = PAR[i]
        for k in idx:
            builder.cx(C[k], q)
        if const:
            builder.x(q)

    def parity_off(builder, C, q, i):
        idx, const = PAR[i]
        if const:
            builder.x(q)
        for k in reversed(idx):
            builder.cx(C[k], q)

    def _product_carry(builder, C, q, i):
        shift = 8 - i
        if 0 <= i < 8 and (MASK >> i) & 1:                    # F3: plain rounding, the data wire C[s]
            assert C[shift] is not False
            return C[shift]
        if shift < 0:                                         # F3: row 9 parity 0: carry = NOT q
            builder.x(q); return q
        if shift == 0:
            assert C[0] is False; return False
        carry = builder.alloc_qubits(1)[0]; upper = C[shift]; lower = C[shift - 1]
        if lower is False:
            parity_on(builder, C, q, i); builder.ccx(upper, q, carry); parity_off(builder, C, q, i)
        else:
            parity_on(builder, C, q, i)                       # F3: parity applied before upper is mixed into lower
            builder.cx(upper, lower); builder.x(q)
            builder.ccx(lower, q, carry)
            builder.x(q); builder.cx(upper, lower)
            parity_off(builder, C, q, i)
            builder.cx(upper, carry)
        return carry

    def _clear_product_carry(builder, C, q, i, carry):
        shift = 8 - i
        if 0 <= i < 8 and (MASK >> i) & 1:                    # F3: nothing owned
            assert carry is C[shift]
            return
        if shift < 0:
            builder.x(q); return
        if shift == 0:
            assert carry is False; return
        upper = C[shift]; lower = C[shift - 1]
        if lower is False:
            parity_on(builder, C, q, i); builder.del_and(upper, q, carry); parity_off(builder, C, q, i)
        else:
            builder.cx(upper, carry)
            parity_on(builder, C, q, i)
            builder.cx(upper, lower); builder.x(q)
            builder.del_and(lower, q, carry)
            builder.x(q); builder.cx(upper, lower)
            parity_off(builder, C, q, i)
        builder.free(carry)

    installed = lc.product_inverse                            # patch_seedinv's closure (this patch follows it)
    g = installed.__globals__
    assert 'STATE' in g and 'product_inverse' in g
    STATE = g['STATE']

    def product_inverse(builder, v, D, word, token, dirty, *, baseline):
        D, word, dirty = list(D), list(word), list(dirty)
        assert len(D) == 13 and len(word) == 16 and len(token.rows) == 10
        entry = builder.num_allocated_qubits
        outcomes = {}
        for i in (0, 1, 3):
            old = v[i]
            assert isinstance(old, km.q) and dirty.count(old) == 1
            m = builder.alloc_bits(1)[0]
            builder.hmr(old, m)
            builder.free(old)
            outcomes[i] = m
            v[i] = None
            dirty.remove(old)
        C = [False] + D
        builder.x(word[15])
        for i in reversed(range(10)):
            debt = token.rows[i]
            q = word[6 + i]
            target = word[:6 + i]
            source = C[9 - i:] if i < 9 else C
            for t in target:
                builder.cx(q, t)
            carry = baseline._product_carry(builder, C, q, i)

            def erase(value, i=i, q=q):
                baseline._clear_product_carry(builder, C, q, i, value)

            def regenerate(i=i, q=q):
                return baseline._product_carry(builder, C, q, i)

            plain = i < 8 and (MASK >> i) & 1
            released = baseline.product_ar.inverse_add(
                builder, source, target, debt, carry_in=carry,
                dirty_wires=dirty + word[6:],
                carry_lifecycle=(erase, regenerate) if (i < 8 and not plain) else None)   # F3: none for plain rows
            if not released:
                baseline._clear_product_carry(builder, C, q, i, carry)
            for t in target:
                builder.cx(q, t)
            assert debt.consumed and not debt.measurements
        builder.x(word[15])
        assert builder.num_allocated_qubits == entry - 3
        STATE['outcomes'] = outcomes
        return word[:6], word[6:]

    bl._product_carry = _product_carry
    bl._clear_product_carry = _clear_product_carry
    g['product_inverse'] = product_inverse
