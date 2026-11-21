"""dqi_var_route.py -- B7 (dqi_var), fourth holder: a layout-free final routing for the line's decode (patch token
`var_route`, install LAST). Zero Toffolis: it only changes which free wires a few data bits pass through.

The trap it removes (S3's recipe 4; B5's q8u, my q5 + KEEP 4..6): the kit's decode calls
`route_words(b, heads + real, S)` between un-synthesis and unextract to put the data back on the input words S, and
`common.reclaim` asserts that every S wire not holding data is free at that moment. A register that is still live
there (a parked flag, G, a pole flag) and that a changed allocation order placed on an S wire breaks the build.
Here route_words sends each data bit whose S wire is still held to a spare free wire instead, writes that spare
into the word list S in place (unextract, uninitial and restore_masks run on the list, so they follow it), and
records the move; right after the decode's last step (baseline.restore) every recorded bit is swapped onto its S
wire, which is free by then, and the spare is freed. Every other bit is routed exactly as the kit routes it.
Env: DQI_VAR_ROUTE_LOG=1 prints the moves.
"""
import os
import sys

STATS = dict(calls=0, held=0, final_swaps=0)
# pending moves per builder (scratch decodes of other patches run on their own builders and may stop early)
_PENDING = {}


def _soft_route_words(b, source_words, target_words, km):
    """The kit's route_words (exp/round7/a05/raw_return.py) with held target wires replaced by spares."""
    STATS['calls'] += 1
    source = [q for w in source_words for q in w]
    target = [q for w in target_words for q in w]
    source_set = set(source)
    target_set = set(target)
    missing = [q for q in target if q not in source_set]
    extras = [q for q in source if q not in target_set]
    assert len(missing) == len(extras)
    if missing:
        wanted = {str(q) for q in missing}
        available = b.alloc_qubits(b.finish_circuit().num_qubits - b.num_allocated_qubits)
        have = {str(q): q for q in available}
        held = [q for q in missing if str(q) not in have]
        spares = [q for q in available if str(q) not in wanted]
        assert len(spares) >= len(held), ('not enough spare wires', len(spares), len(held))
        subst = {}
        for q, s in zip(held, spares):
            subst[str(q)] = s
            _PENDING.setdefault(id(b), []).append((s, str(q)))
        keep = {str(q) for q in missing if str(q) in have} | {str(s) for s in subst.values()}
        b.free(km.array([q for q in available if str(q) not in keep]))
        if held:
            STATS['held'] += len(held)
            # the word list S follows the spares
            for wi, w in enumerate(target_words):
                qs = [subst.get(str(q), q) for q in w]
                if any(str(q) in subst for q in w):
                    target_words[wi] = km.array(qs)
            target = [subst.get(str(q), q) for q in target]
            missing = [subst.get(str(q), q) for q in missing]
            if os.environ.get('DQI_VAR_ROUTE_LOG') == '1':
                print('dqi_var_route: held', [str(q) for q in held], '-> spares', [str(s) for s in subst.values()],
                      flush=True)
    labels = {q: q for q in source + missing}
    places = {q: q for q in source + missing}
    for dest, wanted_q in zip(target, source):
        where = places[wanted_q]
        if where != dest:
            old = labels[dest]
            b.swap(dest, where)
            labels[dest], labels[where] = wanted_q, old
            places[wanted_q], places[old] = dest, where
    if extras:
        b.free(extras)


def _final_swaps(b, km):
    """After the decode's last step: every recorded data bit goes from its spare to its own S wire."""
    pending = _PENDING.pop(id(b), [])
    if not pending:
        return
    available = b.alloc_qubits(b.finish_circuit().num_qubits - b.num_allocated_qubits)
    have = {str(q): q for q in available}
    need = {wid for _, wid in pending}
    assert need <= set(have), ('S wires still held at the end', sorted(need - set(have)))
    for spare, wid in pending:
        b.swap(spare, have[wid])
        STATS['final_swaps'] += 1
    b.free(km.array([q for q in available if str(q) not in need]))
    b.free(km.array([s for s, _ in pending]))


class _BaselineProxy:
    """Delegates to the decode's baseline object; its restore() runs the final swaps after the original."""

    def __init__(self, inner, km):
        object.__setattr__(self, '_inner', inner)
        object.__setattr__(self, '_km', km)

    def __getattr__(self, name):
        attr = getattr(object.__getattribute__(self, '_inner'), name)
        if name != 'restore':
            return attr
        km = object.__getattribute__(self, '_km')

        def restore(b, *a, **k):
            out = attr(b, *a, **k)
            _final_swaps(b, km)
            return out
        return restore


def install(km_unused, ns):
    import kickmix as km
    dec = ns.decode
    g = dec.__globals__
    assert 'route_words' in g and 'g' in g and 'baseline' in g['g'], 'var_route: not the kit decode it expects'
    g['route_words'] = lambda b, src, tgt: _soft_route_words(b, src, tgt, km)
    g['g']['baseline'] = _BaselineProxy(g['g']['baseline'], km)
    return dict(route='soft', stats=STATS, decode=getattr(dec, '__module__', '?'))
