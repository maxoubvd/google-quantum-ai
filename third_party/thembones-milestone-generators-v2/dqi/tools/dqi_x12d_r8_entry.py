"""dqi_x12d_r8_entry.py -- builder J's 24-bit entry (tools/dqi_junc_r8_entry.py, read-only, run unchanged with its own
environment) plus my converter-side patches (workstream dqi_x12d), installed after J's, in order:
  DQI_X12D_R8P   comma list of my patch modules `dqi_x12d_<name>` (each has install(km, ns), as J's dqi_junc_r8p_*)
decode(b, S, E) for tools/dqi_convert_build.py.
"""
import importlib
import importlib.util
import os
import sys
import types

_TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _TOOLS)
_spec = importlib.util.spec_from_file_location('dqi_junc_r8_entry', os.path.join(_TOOLS, 'dqi_junc_r8_entry.py'))
_junc = importlib.util.module_from_spec(_spec)
sys.modules['dqi_junc_r8_entry'] = _junc
_spec.loader.exec_module(_junc)
decode = _junc.decode

PATCH_LOG = dict(_junc.PATCH_LOG)
for _name in [p for p in os.environ.get('DQI_X12D_R8P', '').split(',') if p]:
    _mod = importlib.import_module('dqi_x12d_' + _name)
    PATCH_LOG['x12d_' + _name] = _mod.install(None, types.SimpleNamespace(decode=decode))
print('dqi_x12d_r8_entry patches:', {k: v for k, v in PATCH_LOG.items() if k.startswith('x12d_')}, flush=True)

__all__ = ['decode']
