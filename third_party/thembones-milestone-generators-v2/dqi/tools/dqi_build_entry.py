"""Our DQI generator entry: the public kit's q26 carry decoder plus the dqi_build patches.

decode(b, S, E) for tools/dqi_convert_build.py. The kit source tree is read from DQI_BUILD_SOURCE
(default ~/qcc/ws/dqi_build/src/source, a verbatim copy of the public bundle's `source` directory).
Parameters (environment, with the shipped defaults):
  DQI_BUILD_FIRST   first group on the homogeneous route (kit: 5; shipped: 3, the minimum for prefix5)
  DQI_BUILD_PATCHES comma list of patch modules to install, in order (default: prefix5; prefix3 is the
                    earlier J = 3 block and must not be combined with it)
"""
import importlib
import os
import sys
import types
from pathlib import Path

ROOT = Path(os.environ.get('DQI_BUILD_SOURCE', os.path.expanduser('~/qcc/ws/dqi_build/src/source'))).resolve()
sys.path[:0] = [str(ROOT), str(ROOT / 'exp/round7/root'), str(ROOT / 'exp/next/carry')]
import decode_next_best as _base
import flag_carry

flag_carry.FIRST = int(os.environ.get('DQI_BUILD_FIRST', '3'))
CARRY_LOG = flag_carry.install(_base.decode)
from decode_trace_merged_completed_return import decode

sys.path.insert(0, str(Path(__file__).resolve().parent))
PATCH_LOG = {}
for _name in [p for p in os.environ.get('DQI_BUILD_PATCHES', 'prefix5').split(',') if p]:
    _mod = importlib.import_module('dqi_build_patch_' + _name)
    PATCH_LOG[_name] = _mod.install(None, types.SimpleNamespace(decode=decode))

__all__ = ['decode']
