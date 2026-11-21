"""dqi_junc_r8_entry.py -- the submitted decoder's generator (best_r8: tools/dqi_build_entry.py, read-only, same logic)
plus workstream dqi_junc's patches for the 24-bit in-place line.

decode(b, S, E) for tools/dqi_convert_build.py. Environment as tools/dqi_build_entry.py (DQI_BUILD_SOURCE,
DQI_BUILD_FIRST, DQI_BUILD_PATCHES, DQI_BUILD_R4JT) and DQI_JUNC_R8P: comma list of my patch modules
`dqi_junc_r8p_<name>` installed after the build line's patches, in order.
"""
import importlib
import os
import sys
import types
from pathlib import Path

ROOT = Path(os.environ.get('DQI_BUILD_SOURCE', os.path.expanduser('~/qcc/ws/dqi_build/src/source'))).resolve()
sys.path[:0] = [str(ROOT), str(ROOT / 'exp/round7/root'), str(ROOT / 'exp/next/carry')]
import decode_next_best as _base                                      # noqa: E402
import flag_carry                                                      # noqa: E402

flag_carry.FIRST = int(os.environ.get('DQI_BUILD_FIRST', '3'))
CARRY_LOG = flag_carry.install(_base.decode)
from decode_trace_merged_completed_return import decode               # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
PATCH_LOG = {}
for _name in [p for p in os.environ.get('DQI_BUILD_PATCHES', 'prefix5').split(',') if p]:
    _mod = importlib.import_module('dqi_build_patch_' + _name)
    PATCH_LOG[_name] = _mod.install(None, types.SimpleNamespace(decode=decode))
for _name in [p for p in os.environ.get('DQI_JUNC_R8P', '').split(',') if p]:
    _mod = importlib.import_module('dqi_junc_r8p_' + _name)
    PATCH_LOG['junc_' + _name] = _mod.install(None, types.SimpleNamespace(decode=decode))

__all__ = ['decode']
