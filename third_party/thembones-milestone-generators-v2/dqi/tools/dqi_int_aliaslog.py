"""dqi_int_aliaslog.py -- workstream dqi_int (B5, Nov 14): prints Z1's dqi_zn_alias LOG when the build process ends.

Z1's tools/dqi_zn_alias.py (R36-DQI-LOC's cross-word alias) counts at run time how many descents took the short conv
(conv_short), the root alias (root_aliased) and patch B (b_conv_short); its installer returns only the install-time
counts. This module registers an atexit print of the live LOG so a full decode's build.log shows them.
Install right after dqi_zn_alias in DQI_INT_EXT. It emits no gate and changes no circuit.
"""
import atexit
import sys


def install(km, ns):
    mod = sys.modules.get('dqi_zn_alias')
    if mod is None:
        return {'alias': 'absent'}
    atexit.register(lambda: print('dqi_int_aliaslog', dict(mod.LOG), flush=True))
    return {'alias_log_at_install': dict(mod.LOG)}
