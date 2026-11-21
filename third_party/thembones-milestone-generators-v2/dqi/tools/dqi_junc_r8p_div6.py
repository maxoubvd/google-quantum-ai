"""dqi_junc_r8p_div6.py -- the un-synthesis division's 6-bit inverse without retained coordinates.

The division (exp/round7/root/norm_deferred_flag_return.divide_from_norms, as patched by r6n / latem) inverts the
divisor's 6-bit subnorm M in place with exp/frontier/a14/decode_retained19.inverse6_retained19 (a searched 19-product
plan, 14 live qubits with M, two source coordinates retained for norm_cache.rebuild), the peak of un-synthesis at 905
(ip8: 864 + virtual leader 12 + sweep leader 12 + N 12 + 14 + flags 3). Here, through the two names the fold divider
reads from its own globals (decode_inverse_retained's environment):
  produce_retained(b, M): the in-place GF(64) inverse of workstream dqi_junc's 12-bit line (tools/dqi_junc_div.py
      inv6_inplace: GF(8) norm 6, its inverse 3, conj, the two 3-bit halves scaled through a 3-bit target), built on
      the kit's own functions: 21 Toffolis, 6 private (12 live with M); returns an all-zero transcript mask (the kit's
      cz(M, mask) becomes nothing) and its own masks in place of the retained coordinates;
  rebuild_retained(b, N, kept, k6): M = Norm6(N) recomputed into a fresh register (one 6-bit norm), then the
      inverse's masks (on conj(M_old)'s halves, linear in M_old) paid by CZ on it.
Width at the inverse: 12 live with M instead of 14. Toffolis: 21 + the norm against 19 + the cached rebuild.
"""
import sys

STATS = dict(produce=0, rebuild=0)


def _make(G):
    ar = G['ar']
    norm = G['norm_xor_noscratch']
    mul = G['mul_xor_noscratch']
    inv3 = sys.modules['exp.noscratch.kernels'].inverse3_noscratch
    TO6, FRM6, K6 = ar.tower(6)

    def measured(b, reg):
        m = b.alloc_bits(len(reg))
        b.hmr(reg, m)
        return m

    def produce(b, M):
        STATS['produce'] += 1
        ar.linear_inplace(b, M, TO6)
        N3 = b.alloc_qubits(3)
        norm(b, M, N3, K6)                         # Norm of M over GF(8)
        P = inv3(b, N3)                            # 1/Norm
        ar.erase_inverse(b, P, N3)                 # N3 = 1/P: free
        b.free(N3)
        b.cx(M[3:], M[:3])                         # conj(M) = (a + c, c)
        masks = []
        for half in (M[:3], M[3:]):
            t = b.alloc_qubits(3)
            mul(b, half, P, t)
            masks.append(measured(b, half))
            b.swap(half, t)
            b.free(t)
        ar._erase_norm(b, M, P, K6)                # P = Norm(new M): free
        b.free(P)
        ar.linear_inplace(b, M, FRM6)
        zero = b.alloc_bits(6)
        b.bit_store(zero, False)
        return zero, masks

    def rebuild(b, N, kept, k6):
        STATS['rebuild'] += 1
        M = b.alloc_qubits(6)
        norm(b, N, M, k6)                          # M_old = Norm6(N) (N is in the tower basis here)
        ar.linear_inplace(b, M, TO6)
        b.cx(M[3:], M[:3])
        b.cz(M[:3], kept[0])
        b.cz(M[3:], kept[1])
        b.cx(M[3:], M[:3])
        ar.linear_inplace(b, M, FRM6)
        for m in kept:
            b.free(m)
        return M
    return produce, rebuild


def install(km, module):
    walk = sys.modules['walk_bound']
    done, n = set(), 0
    for fn in walk.functions(module.decode):
        G = fn.__globals__
        if id(G) in done or 'produce_retained' not in G or 'rebuild_retained' not in G:
            continue
        if not fn.__code__.co_filename.endswith('decode_inverse_retained.py'):
            continue
        done.add(id(G))
        produce, rebuild = _make(G)
        G['produce_retained'] = produce
        G['rebuild_retained'] = rebuild
        n += 1
    assert n, 'the fold divider with produce_retained / rebuild_retained was not found'
    import atexit
    atexit.register(lambda: print('dqb_div6', STATS, file=sys.stderr))
    return dict(div6=n)
