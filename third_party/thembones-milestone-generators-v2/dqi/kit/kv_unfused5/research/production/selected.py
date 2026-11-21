"""Selected back half, assembled without displaced intermediate installers.

All arithmetic helpers remain in their frozen research source files. The
private-graph copy between the forward assembly and the reverse additions is
intentional: nested closure dictionaries and their shared caches must retain
the same ownership as the measured composition.
"""
from pathlib import Path
import inspect
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'research/agents' / name)
                for name in ('a12', 'a11', 'r06_joint_top', 'a09',
                             'a08', 'a07', 'a01')]
import hankel_compact
import hankel_cache_even
import joint_cached
import forward_norm_cache
import known_zero_inverse
import middle_fusion
import blocked_three_noscratch
import root_fused_three
import scalar_compact_v3
from rebind_graph import private_rebind

p = hankel_compact.p
cv = p.install_dot_cache(p.install_base_zero(p.install_undo_zeros(
    p.install_child_inverse(p.install_retained_gaps(p.make_converter())))))
# Z2 kv_unfused3: the forward descent is NOT fused (no a07 middle fusion, no level-14 three-frame block); the
# root inverse stays ephemeral (a07's retire / rebuild, with the full root pre action). The arithmetic helpers of
# the three-frame block stay in the scope (a14's reverse fusion and a18 use payload_middle).
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location('ephemeral_root_only', str(ROOT / 'research/agents/a07z/ephemeral_root_block.py'))
ephemeral_root_only = _ilu.module_from_spec(_spec); _spec.loader.exec_module(ephemeral_root_only)
cells = dict(zip(cv.__code__.co_freevars, cv.__closure__))
short = cells['core'].cell_contents.__globals__
scope = short['old'].direct_post.__globals__
scope.update(ordinary_solve=scope['look_solve'],
             convolution=lambda b, x, y, out: scope['dk'].conv(b, x, y, out, len(x) + len(y) - 1),
             phase_classical_parity=blocked_three_noscratch.phase_classical_parity)
for fun in (blocked_three_noscratch.mask_word, blocked_three_noscratch.erase_mask, blocked_three_noscratch.copy_flag,
            blocked_three_noscratch.phase_flag, blocked_three_noscratch.payload_middle,
            blocked_three_noscratch.pre_partial):
    scope[fun.__name__] = p.clone(fun, scope)
ephemeral_root_only.install(cv, blocked_three_noscratch.p)
cv = hankel_cache_even.install(cv)
cv = joint_cached.install(cv, root_sources=18)
cv = private_rebind(cv, p.base.replacements, p.ROOT)
decode = scalar_compact_v3.make_decoder(converter=cv, reference_decode=p.decode)
forward_norm_cache.install(scalar_compact_v3.converter_of(decode), levels=(11, 5))
known_zero_inverse.install(scalar_compact_v3.converter_of(decode))

sys.path[:0] = [str(ROOT / 'research/agents' / name)
                for name in ('a14', 'r08_weight_handoff', 'a13')]
import reverse_fusion
import reverse_caches
import weight_handoff_v3

decode = private_rebind(decode, {}, p.ROOT)
cv = scalar_compact_v3.converter_of(decode)
scope = scalar_compact_v3.scope_of(cv)
scope['payload_middle'] = p.clone(weight_handoff_v3.payload_middle, scope)
reverse_caches.install(cv)
import reverse_root_lifetime  # Z1 kv_unfused5: a14 reverse fusion NOT installed, a12 kept
reverse_root_lifetime.install(cv)

from inverse_character_width import install_decoder
install_decoder(decode)
