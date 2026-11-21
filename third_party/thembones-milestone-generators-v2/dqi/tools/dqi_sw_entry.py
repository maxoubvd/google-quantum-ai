"""dqi_sw_entry.py -- S3 sweep: B5's decode entry (tools/dqi_int_entry.py, read-only logic copied) with one addition:
a token `sw_<name>` in DQI_INT_P installs my patch module tools/dqi_sw_<name>.py at that position (so that a patch can
sit before a3, whose SYNC capture must see it). Every other token loads B5's `dqi_int_<name>` as his entry does.
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
    if _name.startswith('sw_'):
        _mod = importlib.import_module('dqi_sw_' + _name[3:])
    elif _name.startswith('var_'):          # B7's patch modules (read-only), e.g. var_keep with DQI_VAR_KEEP_ADD
        _mod = importlib.import_module('dqi_var_' + _name[4:])
    elif _name.startswith('pub_'):          # D1's patch modules (read-only), e.g. pub_pfd with DQI_PUB_PFD_RINF
        _mod = importlib.import_module('dqi_pub_' + _name[4:])
    else:
        _mod = importlib.import_module('dqi_int_' + _name)
    PATCH_LOG['int_' + _name] = _mod.install(None, types.SimpleNamespace(decode=decode))
print('dqi_sw_entry patches:', {k: v for k, v in PATCH_LOG.items() if k.startswith('int_')}, flush=True)

__all__ = ['decode']
