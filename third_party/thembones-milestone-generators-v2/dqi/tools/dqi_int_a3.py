"""dqi_int_a3.py -- workstream dqi_int (B5): builder D1's transplant a3 (the projective finish in the converter and the
scan's node polynomial P -> P * M; tools/dqi_pub_a3.py + tools/dqi_pub_finish.py, read-only) as a DQI_INT_P patch of
tools/dqi_int_entry.py. Mechanism credit: the projective finish of the public DQI decoder (Syed Asad Ali,
Apache-2.0), re-implemented by D1 on our kit's kernels. Install it after flagpark and retl:
DQI_INT_P=retj,polepark2,flagpark,retl,div6f,a3.
DQI_INT_A3_SNAP=<dir>: load a frozen copy of D1's two files from <dir> (provenance of a build; dqi_pub_a3 puts its
own directory first on sys.path, so dqi_pub_finish comes from the same copy); default: D1's live files in tools/.
"""
import importlib
import importlib.util
import os
import sys


def install(km, ns):
    # (Nov 14 20:35Z) DQI_INT_A3_SNAP_OVR=<dir> takes precedence (the pinned recipes set DQI_INT_A3_SNAP themselves);
    # unset: unchanged.
    snap = os.environ.get('DQI_INT_A3_SNAP_OVR') or os.environ.get('DQI_INT_A3_SNAP')
    if snap:
        spec = importlib.util.spec_from_file_location('dqi_pub_a3', os.path.join(snap, 'dqi_pub_a3.py'))
        mod = importlib.util.module_from_spec(spec)
        sys.modules['dqi_pub_a3'] = mod
        spec.loader.exec_module(mod)
        assert os.path.dirname(os.path.abspath(sys.modules['dqi_pub_finish'].__file__)) == os.path.abspath(snap)
    else:
        mod = importlib.import_module('dqi_pub_a3')
    # DQI_INT_A3_FDIV=dt42|dt (with DQI_PUB_A3_DIV=dr): B6's wider direct divisions (309.89 at 42 / 294.93 at 48
    # ancillas, tools/dqi_kern_div.py) in the finish's reverse in place of dr (315.00 at 36). Only the finish calls
    # dqi_kern_div.divide_direct, so rebinding it there changes nothing else.
    # DQI_INT_A3_FDIV=auto: per call, the cheapest form whose width fits under DQI_INT_WCAP (default 902) at that
    # moment (live qubits at entry + its ancillas): dt (48), dt42 (42), else dr (36). Choices printed at exit.
    fdiv = os.environ.get('DQI_INT_A3_FDIV')
    if fdiv:
        import dqi_kern_div as KD
        if fdiv == 'auto':
            import atexit
            cap = int(os.environ.get('DQI_INT_WCAP', '902'))
            forms = [(48, KD.divide_direct_tower, 'dt'), (42, KD.divide_direct_tower42, 'dt42'),
                     (0, KD.divide_direct, 'dr')]
            seen = []

            def chosen(b, x, y, K, g):
                n = b.num_allocated_qubits
                w, fn, name = next(f for f in forms if n + f[0] <= cap or f[0] == 0)
                seen.append((n, name))
                return fn(b, x, y, K, g)
            KD.divide_direct = chosen
            atexit.register(lambda: print('dqi_int_a3 fdiv auto (live at entry, form):', seen, flush=True))
        else:
            KD.divide_direct = {'dt42': KD.divide_direct_tower42, 'dt': KD.divide_direct_tower}[fdiv]
    return {'pub_a3': mod.install(km, ns), 'snap': snap or 'live', 'fdiv': fdiv or 'dr'}
