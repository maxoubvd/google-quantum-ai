"""dqi_pub_entry.py -- builder D1: B5's decode (tools/dqi_int_entry.py with DQI_INT_P / DQI_INT_FP, read-only) plus
D1's patches, in order: DQI_PUB_P = comma list of tools/dqi_pub_<name>.py modules with install(km, ns)
(default 'a3', the projective finish). decode(b, S, E) for tools/dqi_convert_build.py.
"""
import importlib
import importlib.util
import os
import sys
import types

_TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _TOOLS)
_spec = importlib.util.spec_from_file_location('dqi_int_entry', os.path.join(_TOOLS, 'dqi_int_entry.py'))
_int = importlib.util.module_from_spec(_spec)
sys.modules['dqi_int_entry'] = _int
_spec.loader.exec_module(_int)
decode = _int.decode

PATCH_LOG = dict(_int.PATCH_LOG)
for _name in [p for p in os.environ.get('DQI_PUB_P', 'a3').split(',') if p]:
    _mod = importlib.import_module('dqi_pub_' + _name)
    PATCH_LOG['pub_' + _name] = _mod.install(None, types.SimpleNamespace(decode=decode))
print('dqi_pub_entry patches:', {k: v for k, v in PATCH_LOG.items() if k.startswith('pub_')}, flush=True)

__all__ = ['decode']
