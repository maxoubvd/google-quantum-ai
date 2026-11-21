"""Patch for the q26 generator: hoisted-L prefix, J = 3 (workstream dqi_build).

Identity (tools/dqi_build_prefix_check.py): for a later node g the L-products of chain steps 1 and 2 equal
    L0 += U_g (x + H_0)^2 + b_2,   U_g = [a_1 + b_1 (g^2 + k_0) + b_2 (a_1^2 + b_1^2 g^2)] / q_0(g)^2,
applied to the jet before the chain (forward) or after it (return). So:
  * chain steps with head 1 or head 2 become N-only products (99 instead of 182..192);
  * one block before the first homogeneous group (and the same block after its return) computes the two words
    Pa = b_2 a_1^2, Pb = b_2 b_1^2 and does one 99-Toffoli product per homogeneous group.
Usable with dqi_extract_gen.py / dqi_build_stage_means.py (--patch) or imported by an entry file.
"""
import sys

NONLY = frozenset({1, 2})

CHAIN_SRC = '''
def chain_step(b, head, N, D, L, q, ct, c, reverse=False, first=False):
    if first:
        return core.external_chain_step(b, head, N, D, L, q, ct, c, reverse, first)
    if reverse:
        core.scalepair(b, N, core.inv(q))
        b.swap(N, D)
        core.scalepair(b, L, core.inv(q))
    if head[0][0].id // 24 in _dqb_nonly:
        _dqb_mul_head(b, head, D, N, ct)
    else:
        a, h = head
        core.cx(b, h, a, ct)
        accumulator(b, core.combine(head), D, list(N) + list(L), joint_tensor(c))
        core.cx(b, h, a, ct)
    if not reverse:
        core.scalepair(b, L, q)
        b.swap(N, D)
        core.scalepair(b, N, q)
'''

EXTRACTION_SRC = '''
def external_extraction(b,words,reverse=False,token=None):
    history=[]if not reverse else list(token)
    for gi in reversed(range(17))if reverse else range(17):
        g=groups[gi];u=combine(words[g['word']:g['word']+2]);v=combine(words[g['word']+2:g['word']+4])
        a,d=words[g['word']:g['word']+2]
        if not reverse:
            if gi==min(EXTERNAL_HOMOGENEOUS):_dqb_block(b,words,EXTERNAL_HOMOGENEOUS)
            nt=external_node_chain(b,words,gi,normalizer=None)if gi in EXTERNAL_HOMOGENEOUS else external_ordinary_node(b,words,gi)
            cx(b,d,a,g['c']);b.cx(d,v[:12]);local=inverse24_forward190(b,v,K)
            cx(b,words[g['word']+3],words[g['word']+2],g['c'])
            history.append((nt,local))
        else:
            nt,local=history.pop();cx(b,words[g['word']+3],words[g['word']+2],g['c'])
            inverse24_reverse209(b,v,local,K);b.cx(d,v[:12]);cx(b,d,a,g['c'])
            if gi in EXTERNAL_HOMOGENEOUS:external_node_chain(b,words,gi,True,nt,None)
            else:external_ordinary_node(b,words,gi,True,nt)
            if gi==min(EXTERNAL_HOMOGENEOUS):_dqb_block(b,words,EXTERNAL_HOMOGENEOUS)
    assert not reverse or not history
    return history
'''


def make_block(G):
    """The hoisted block, built from the primitives of one final_decoder-style namespace G."""
    ar, cx, scalepair, combine = G['ar'], G['cx'], G['scalepair'], G['combine']
    groups, ks, K, mul, inv = G['groups'], G['ks'], G['K'], G['mul'], G['inv']
    square12, square24, mul24 = G['square12'], G['external_square24'], G['mul24_xor99']

    def block(b, words, homogeneous):
        with b.mark('dqb_prefix3'):
            a0, h0, a1, b1, b2 = words[0], words[1], words[2], words[3], words[5]
            pa, pb = b.alloc_qubits(12), b.alloc_qubits(12)
            for src, dst in ((a1, pa), (b1, pb)):
                square12(b, src)
                ar.mul_xor0(b, b2, src, dst)
                square12(b, src, True)
            m = combine([a1, b1])
            for gi in sorted(homogeneous):
                g = groups[gi]
                ct = g['c']
                N = combine(words[g['word']:g['word'] + 2])
                L = combine(words[g['word'] + 2:g['word'] + 4])
                e = ct ^ g['kap']
                q0 = ks[0] ^ g['kap']
                q0sq = mul(q0, q0)
                cx(b, h0, a0, ct); b.cx(combine([a0, h0]), N); cx(b, h0, a0, ct)
                square24(b, N)
                b.cx(pa, a1); cx(b, b1, a1, ks[0]); b.cx(pb, b1); cx(b, b1, a1, e)
                scalepair(b, L, q0sq)
                mul24(b, m, N, L, K)
                scalepair(b, L, inv(q0sq))
                cx(b, b1, a1, e); b.cx(pb, b1); cx(b, b1, a1, ks[0]); b.cx(pa, a1)
                square24(b, N, True)
                cx(b, h0, a0, ct); b.cx(combine([a0, h0]), N); cx(b, h0, a0, ct)
                b.cx(b2, L[:12])
            for src, dst in ((b1, pb), (a1, pa)):
                square12(b, src)
                ar.erase_product(b, b2, src, dst)
                square12(b, src, True)
            b.free(pb); b.free(pa)
    return block


def install(km, module):
    walk = sys.modules['walk_bound']
    fns = walk.functions(module.decode)
    scope = {}
    exec(compile(CHAIN_SRC, __file__ + ':chain', 'exec'), scope)
    chain_code = scope['chain_step'].__code__
    scope = {}
    exec(compile(EXTRACTION_SRC, __file__ + ':extraction', 'exec'), scope)
    extraction_code = scope['external_extraction'].__code__
    need = ('ar', 'cx', 'scalepair', 'combine', 'groups', 'ks', 'K', 'mul', 'inv', 'square12',
            'external_square24', 'mul24_xor99')
    mul_head = next(f for f in fns if f.__name__ == 'external_mul_head' and all(k in f.__globals__ for k in need))
    block = make_block(mul_head.__globals__)
    log = dict(chain=0, extraction=0)
    for f in fns:
        if f.__name__ == 'chain_step' and 'joint_tensor' in f.__code__.co_names:
            assert f.__closure__ is None
            f.__globals__.update(_dqb_nonly=NONLY, _dqb_mul_head=mul_head)
            f.__code__ = chain_code
            log['chain'] += 1
        elif f.__name__ == 'external_extraction' and 'inverse24_forward190' in f.__code__.co_names:
            assert f.__closure__ is None
            f.__globals__['_dqb_block'] = block
            f.__code__ = extraction_code
            log['extraction'] += 1
    assert log['chain'] and log['extraction'], log
    return log
