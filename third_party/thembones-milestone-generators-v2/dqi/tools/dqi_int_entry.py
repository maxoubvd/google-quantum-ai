"""dqi_int_entry.py -- workstream dqi_int (builder B5): the 24-bit decoder of builder J (tools/dqi_junc_r8_entry.py)
with x12d's converter-side patches (tools/dqi_x12d_r8_entry.py, DQI_X12D_R8P), both run unchanged (read-only), then
my patches, in order:
  DQI_INT_P   comma list of my patch modules `dqi_int_<name>` (each has install(km, ns) like J's dqi_junc_r8p_*)
decode(b, S, E) for tools/dqi_convert_build.py.
"""
import importlib
import importlib.util
import os
import sys
import types

_TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _TOOLS)
_spec = importlib.util.spec_from_file_location('dqi_x12d_r8_entry', os.path.join(_TOOLS, 'dqi_x12d_r8_entry.py'))
_x12d = importlib.util.module_from_spec(_spec)
sys.modules['dqi_x12d_r8_entry'] = _x12d
_spec.loader.exec_module(_x12d)
decode = _x12d.decode

PATCH_LOG = dict(_x12d.PATCH_LOG)
for _name in [p for p in os.environ.get('DQI_INT_P', '').split(',') if p]:
    _mod = importlib.import_module('dqi_int_' + _name)
    PATCH_LOG['int_' + _name] = _mod.install(None, types.SimpleNamespace(decode=decode))
print('dqi_int_entry patches:', {k: v for k, v in PATCH_LOG.items() if k.startswith('int_')}, flush=True)

__all__ = ['decode']
