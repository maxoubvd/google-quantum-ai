#!/usr/bin/env python3
"""dqi_int_spec901.py -- B5 (17th holder, Nov 16): port a park SPEC of tools/dqi_int_park901b.py from the build it was
written for to another build of the same converter (wires renamed, ops shifted), through the two builds' 902 censuses
(tools/r48_dqi_ret2_rp_idle.py 'park wire' lines) and the new build's flag checks (tools/r48_dqi_ret2_rp_flagchk.py).

Matching: a census park is keyed by (birth stack, closing stack, gap length); the k-th old park of a key (by gap start)
maps to the k-th new park of that key.  Spec lines:
  park f lt nt W        -> the matched new park's f, lt, nt and the word the new flag check names for it;
  open / trans / remake -> each op is an endpoint (or endpoint - 1) of one of the old flag's census gaps: mapped to the
                           same endpoint of the matched new gap; the flag and word likewise;
  hmr op w              -> dropped unless given again by --hmr <op> <wire> (the new build's linpark answer).
Prints, per new park line, the flag check's verdict for the close ('EXACT ON SAMPLES' or not).
Usage: python dqi_int_spec901.py <old census> <old spec> <new census> <new flagchk> <out spec> [--hmr <op> <wire>]...
"""
import re
import sys

PARK = re.compile(r'park wire\s+(\d+)\s+gap (\d+)\.\.(\d+) \((\d+) ops\)\s+covers (\d+) runs\s+born (\S+) in \[([^\]]*)\]'
                  r'\s+closes at \[([^\]]*)\]')
FWORD = re.compile(r'park f\s+(\d+) gap (\d+)\.\.(\d+) born \S+ \[[^\]]*\]: word ([\d,]+) ')
RUN = re.compile(r'^(\w+)\s+(\d+)\.\.(\d+)\s+\d+ ops [C-] \[')
FCHK = re.compile(r'^f\s+(\d+) gap\s+(\d+)\.\.(\d+)\s+writes\s+(-?\d+) \| (.*)$')


def census(path):
    out = []
    for line in open(path, errors='replace'):
        m = PARK.search(line)
        if m:
            f, lt, nt, ln, cov, born, bst, cst = m.groups()
            out.append(dict(f=int(f), lt=int(lt), nt=int(nt), len=int(ln), born=born, key=(bst, cst, int(ln))))
    return out


def flagchk(path):
    word, verdict = {}, {}
    for line in open(path, errors='replace'):
        m = FWORD.search(line)
        if m:
            f, lt, nt, w = m.groups()
            word[(int(f), int(lt), int(nt))] = [int(x) for x in w.split(',')]
            continue
        m = FCHK.match(line.strip())
        if m:
            f, lt, nt, _wr, rest = m.groups()
            verdict[(int(f), int(lt), int(nt))] = rest
    return word, verdict


def main():
    args = sys.argv[1:]
    hmrs = []
    while '--hmr' in args:
        i = args.index('--hmr')
        hmrs.append((int(args[i + 1]), int(args[i + 2])))
        del args[i:i + 3]
    oc, ospec, nc, nfc, out = args
    old, new = census(oc), census(nc)
    word, verdict = flagchk(nfc)
    bykey = {}
    for p in sorted(new, key=lambda p: p['lt']):
        bykey.setdefault(p['key'], []).append(p)
    okey = {}
    for p in sorted(old, key=lambda p: p['lt']):
        okey.setdefault(p['key'], []).append(p)
    match = {}
    used = set()
    for key, lst in okey.items():
        cand = bykey.get(key, [])
        for k, p in enumerate(lst):
            q = cand[k] if k < len(cand) else None
            match[(p['f'], p['lt'], p['nt'])] = q
            if q is not None:
                used.add((q['f'], q['lt'], q['nt']))
    # fallback for a gap whose length or closing kernel changed: the same birth stack, the same rank among the parks
    # of that birth stack (by gap start), not used by a strict match
    def birth_rank(lst):
        d = {}
        for p in sorted(lst, key=lambda p: p['lt']):
            d.setdefault(p['key'][0], []).append(p)
        return d
    ob, nb = birth_rank(old), birth_rank(new)
    for k0, q0 in list(match.items()):
        if q0 is not None:
            continue
        p = [x for x in old if (x['f'], x['lt'], x['nt']) == k0][0]
        r = ob[p['key'][0]].index(p)
        cand = nb.get(p['key'][0], [])
        q = cand[r] if r < len(cand) else None
        if q is not None and 'EXACT ON SAMPLES' not in verdict.get((q['f'], q['lt'], q['nt']), ''):
            print('fallback REFUSED (re-test not exact): %d %d..%d -> %d %d..%d' % (
                p['f'], p['lt'], p['nt'], q['f'], q['lt'], q['nt']))
            q = None
        if q is not None and (q['f'], q['lt'], q['nt']) not in used:
            match[k0] = q
            used.add((q['f'], q['lt'], q['nt']))
            print('fallback match by birth rank: %d %d..%d (%d ops) -> %d %d..%d (%d ops)' % (
                p['f'], p['lt'], p['nt'], p['len'], q['f'], q['lt'], q['nt'], q['len']))
    # endpoint map of every old gap -> new gap
    ends = {}
    for (f, lt, nt), q in match.items():
        if q is None:
            continue
        ends.setdefault(f, []).append(((lt, nt), (q['lt'], q['nt']), q['f']))

    def map_op(f, op):
        for (lt, nt), (qlt, qnt), qf in ends.get(f, []):
            for o, n in ((lt, qlt), (nt, qnt), (lt - 1, qlt - 1), (nt - 1, qnt - 1), (lt + 1, qlt + 1)):
                if op == o:
                    return n, qf, (qf, qlt, qnt)
        raise SystemExit('op %d of flag %d is not at an endpoint of its census gaps' % (op, f))

    def word_of(gk):
        w = word.get(gk)
        assert w and len(w) == 12, ('no 12-wire word for new gap', gk)
        return w

    lines, report = [], []
    for raw in open(ospec):
        s = raw.split('#')[0].split()
        if not s:
            continue
        kind = s[0]
        if kind == 'park':
            f, lt, nt = int(s[1]), int(s[2]), int(s[3])
            q = match.get((f, lt, nt))
            if q is None:
                report.append('SKIPPED old park %s: no match in the new census' % ' '.join(s[:4]))
                continue
            key = (q['f'], q['lt'], q['nt'])
            w = word.get(key)
            assert w and len(w) == 12, ('no 12-wire word', key)
            lines.append('park %d %d %d %s' % (q['f'], q['lt'], q['nt'], ','.join(map(str, w))))
            report.append('park %d %d..%d -> %d %d..%d | %s' % (f, lt, nt, key[0], key[1], key[2],
                                                               verdict.get(key, 'NO VERDICT')[:150]))
        elif kind == 'open':
            op, f, name = int(s[1]), int(s[2]), s[3]
            n, qf, _gk = map_op(f, op)
            lines.append('open %d %d %s' % (n, qf, name))
        elif kind == 'trans':
            o1, o2, f = int(s[1]), int(s[2]), int(s[3])
            n1, qf, gk = map_op(f, o1)
            n2, _, _ = map_op(f, o2)
            lines.append('trans %d %d %d %s' % (n1, n2, qf, ','.join(map(str, word_of(gk)))))
        elif kind == 'remake':
            op, f, name = int(s[1]), int(s[2]), s[4]
            n, qf, gk = map_op(f, op)
            lines.append('remake %d %d %s %s' % (n, qf, ','.join(map(str, word_of(gk))), name))
        elif kind == 'hmr':
            report.append('dropped: %s (give --hmr for the new build)' % raw.strip())
        else:
            raise SystemExit('unknown spec line ' + raw)
    for op, w in hmrs:
        lines.append('hmr %d %d' % (op, w))
    # coverage: every run at the census level must lie inside a park gap, the K interval (open .. remake) of its flag,
    # or follow an hmr of a wire; runs left uncovered get the new census's cover parks that contain them and whose
    # re-test at the close is exact on samples (flag check)
    runs = []
    for ln in open(nc, errors='replace'):
        m = RUN.match(ln)
        if m:
            runs.append((m.group(1), int(m.group(2)), int(m.group(3))))
    ivs = []
    opens = {}
    for ln in lines:
        t = ln.split()
        if t[0] == 'park':
            ivs.append((int(t[2]), int(t[3])))
        elif t[0] == 'open':
            opens[(int(t[2]), t[3])] = int(t[1])
        elif t[0] == 'remake':
            o = opens.get((int(t[2]), t[4]))
            if o is not None:
                ivs.append((o, int(t[1])))
        elif t[0] == 'hmr':
            ivs.append((int(t[1]), int(t[1]) + 150000))     # X's run follows within 92k ops
    def covered(ra, rb, extra=()):
        return any(lo < ra and rb < hi for lo, hi in list(ivs) + list(extra))
    added = []
    for st, ra, rb in runs:
        if covered(ra, rb, added):
            continue
        cand = [p for p in new if p['lt'] < ra and rb < p['nt']]
        cand = [p for p in cand if 'RE-TEST AT THE CLOSE EXACT ON SAMPLES' in verdict.get((p['f'], p['lt'], p['nt']), '')]
        if hmrs and st == 'unconvert' and not cand:
            pass
        if not cand:
            report.append('UNCOVERED run %s %d..%d: no exact park in the new census' % (st, ra, rb))
            continue
        p = max(cand, key=lambda p: sum(1 for s2, a2, b2 in runs if p['lt'] < a2 and b2 < p['nt']))
        w = word.get((p['f'], p['lt'], p['nt']))
        if not w or len(w) != 12:
            report.append('UNCOVERED run %s %d..%d: park %d without a 12-wire word' % (st, ra, rb, p['f']))
            continue
        lines.append('park %d %d %d %s' % (p['f'], p['lt'], p['nt'], ','.join(map(str, w))))
        added.append((p['lt'], p['nt']))
        report.append('ADDED park %d %d..%d (born %s) for run %s %d..%d' % (p['f'], p['lt'], p['nt'], p['born'], st, ra, rb))
    report.append('coverage: %d runs, %d added parks, %d uncovered' % (
        len(runs), len(added), sum(1 for r in report if r.startswith('UNCOVERED'))))
    with open(out, 'w') as fh:
        fh.write('# ported by dqi_int_spec901.py from %s (%s -> %s)\n' % (ospec, oc, nc))
        fh.write('\n'.join(lines) + '\n')
    print('\n'.join(report))
    print('wrote %s: %d lines' % (out, len(lines)))


if __name__ == '__main__':
    main()
