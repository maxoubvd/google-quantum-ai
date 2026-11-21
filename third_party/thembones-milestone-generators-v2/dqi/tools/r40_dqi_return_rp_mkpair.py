"""r40_dqi_return_rp_mkpair.py -- R40-DQI-RETURN re-pricer (Nov 15): a MEASUREMENT DEVICE, not a line change.

D1's tools/dqi_pub_pair.py (7228f649, the line's: pin q8u) rewrites the a3 v5 core's text in core_pair_source() and
asserts that the delta loop `for j in _D1_J: ... delta)` + `b.cx(delta,top)` occurs exactly once. Its replacement
(_DELTA_NEW) only differs from the old text when DQI_PUB_PAIR_DELTA_PROG is set (make_delta); the line sets no delta
program, so in the line the delta-loop rewrite is a no-op. K1's J / J2 snapshots (tools/dqi_k1_a3snap_k1j*/) moved
the delta loop, so the assert stops a line build with them. This copy keeps every rewrite and assert when a delta
program is set and, without one, skips ONLY the delta-loop rewrite (FWD / REV pair rewrites and their asserts kept).
Control: the base build (a3snap_h3) with this copy must reproduce best_q22ch's raw file byte for byte (b725cc46).

    python3 r40_dqi_return_rp_mkpair.py <source dqi_pub_pair.py (7228f649)> <out dqi_pub_pair.py>
"""
import hashlib
import sys

SRC_SHA = '7228f64944a2fb26'
OLD = '''    assert src.count(_DELTA_OLD) == 1, 'a3 v5 delta loop changed'
    return src.replace(_FWD_OLD, _FWD_NEW).replace(_REV_OLD, _REV_NEW).replace(_DELTA_OLD, _DELTA_NEW)
'''
NEW = '''    if not os.environ.get('DQI_PUB_PAIR_DELTA_PROG'):
        # R40-DQI-RETURN re-pricer's copy: without a delta program _DELTA_NEW runs exactly _DELTA_OLD's ops
        LOG['r40rp_delta_rewrite'] = 'skipped (no delta program): old text count %d' % src.count(_DELTA_OLD)
        return src.replace(_FWD_OLD, _FWD_NEW).replace(_REV_OLD, _REV_NEW)
    assert src.count(_DELTA_OLD) == 1, 'a3 v5 delta loop changed'
    return src.replace(_FWD_OLD, _FWD_NEW).replace(_REV_OLD, _REV_NEW).replace(_DELTA_OLD, _DELTA_NEW)
'''


def main(src, out):
    raw = open(src, 'rb').read()
    sha = hashlib.sha256(raw).hexdigest()
    assert sha.startswith(SRC_SHA), 'source differs: %s' % sha[:16]
    text = raw.decode('utf-8').replace('\r\n', '\n')
    assert text.count(OLD) == 1, 'core_pair_source text not found once'
    assert '\nimport os' in text or '\nimport os,' in text or 'import os\n' in text, 'module does not import os'
    new = text.replace(OLD, NEW)
    with open(out, 'w', encoding='utf-8', newline='\n') as f:
        f.write(new)
    print('wrote', out, hashlib.sha256(new.encode('utf-8')).hexdigest()[:16], 'from', sha[:16])


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
