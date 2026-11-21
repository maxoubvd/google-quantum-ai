"""F3 (invsqrt_fw): circuit-level wire moves on a kickmix text circuit (.kmx), checked exactly, then relabelled.

A move (erase / rebuild) frees wire w over an idle interval:
  - after op t1 (its last value touch before the interval): 'HMR qw bM' (M a fresh classical bit);
  - every phase gate on w (Z, CZ, CCZ, Z_POW) between t1 and the rebuild point moves to right after the rebuild
    (w holds the same value again there; the move is checked by the exact simulation);
  - before op tr: the rebuild program (lines that compute w's old value into qw from live wires), then
    'Z qw if bM', which repays the sign (-1)^(m w) the HMR left.
Ops are numbered as in the original file (0-based, the order of invsqrt_fw_sim.Circ).

    python invsqrt_fw_rw.py <in.kmx> <moves.json> <out_prefix>      writes <out_prefix>.kmx (+ .kmb via kickmix)
    python invsqrt_fw_rw.py relabel <in.kmb> <out_prefix> [Q]        smallest-width relabel (+ final SWAPs for x)

moves.json: [{"w": 9, "t1": 6983, "tr": 9840, "program": ["X q34", "CCX q20 q34 q9", "X q34", "CX q43 q9"]}, ...]
"""
import json
import re
import subprocess
import sys

sys.path.insert(0, "<HOME>/Quantum/tools")

KMPY = "PYTHONPATH=$HOME/qcc/data/kickmix $HOME/qcc/.venv/bin/python"
PHASE_NAMES = {"Z", "CZ", "CCZ", "Z_POW"}
QRE = re.compile(r"\bq(\d+)\b")


def read_kmx(path):
    lines = open(path).read().splitlines()
    head = [l for l in lines if l.startswith(("APPEND_TO_REGISTER", "REGISTER"))]
    body = [l for l in lines if not l.startswith(("APPEND_TO_REGISTER", "REGISTER")) and l.strip()]
    return head, body


def qubits(line):
    return [int(x) for x in QRE.findall(line)]


def max_bit(body):
    m = -1
    for l in body:
        for x in re.findall(r"\bb(\d+)\b", l):
            m = max(m, int(x))
    return m


def apply_moves(body, moves):
    nb = max_bit(body) + 1
    ins_after = {}
    ins_before = {}
    drop = set()
    for mv in moves:
        w, t1, tr = mv["w"], mv["t1"], mv["tr"]
        bit = nb
        nb += 1
        moved = []
        for i in range(t1 + 1, tr):
            l = body[i]
            name = l.split()[0]
            if w in qubits(l):
                if name in PHASE_NAMES:
                    moved.append(l)
                    drop.add(i)
                elif name not in ("PUSH_CONDITION", "POP_CONDITION"):
                    raise ValueError(f"move q{w}: value touch inside the interval at op {i}: {l}")
        # a moved gate keeps its condition: inside PUSH blocks the gate itself carries no 'if'; re-attach it
        fixed = []
        for i in sorted(j for j in drop if t1 < j < tr and w in qubits(body[j])):
            l = body[i]
            conds = []
            depth = []
            for j in range(0, i):
                pass
            fixed.append((i, l))
        ins_after.setdefault(t1, []).append(f"HMR q{w} b{bit}")
        prog = list(mv["program"]) + [f"Z q{w} if b{bit}"]
        ins_before.setdefault(tr, []).extend(prog + [("MOVED", i) for i, l in fixed])
    return ins_after, ins_before, drop


def condition_stack_at(body):
    """For every op index, the list of bits of the PUSH_CONDITION blocks open at that op."""
    st = []
    out = []
    for l in body:
        out.append(list(st))
        t = l.split()
        if t[0] == "PUSH_CONDITION":
            st.append(t[-1])
        elif t[0] == "POP_CONDITION":
            st.pop()
    return out


def emit(body, ins_after, ins_before, drop):
    stk = condition_stack_at(body)
    out = []
    for i, l in enumerate(body):
        for x in ins_before.get(i, []):
            if isinstance(x, tuple):
                j = x[1]
                src = body[j]
                conds = stk[j]
                if " if " in src:
                    conds = conds + [src.split(" if ")[1].strip()]
                    src = src.split(" if ")[0]
                if len(conds) == 0:
                    out.append(src)
                elif len(conds) == 1:
                    out.append(f"{src} if {conds[0]}")
                else:
                    for cnd in conds[:-1]:
                        out.append(f"PUSH_CONDITION if {cnd}")
                    out.append(f"{src} if {conds[-1]}")
                    for cnd in conds[:-1]:
                        out.append("POP_CONDITION")
            else:
                out.append(x)
        if i not in drop:
            out.append(l)
        for x in ins_after.get(i, []):
            out.append(x)
    return out


def to_kmb(kmx, kmb):
    cmd = f"{KMPY} -c \"import kickmix as km; km.Circuit.from_file('{kmx}').to_file('{kmb}', format='kmb')\""
    subprocess.run(["bash", "-lc", cmd], check=True)


def main(inp, moves_json, prefix):
    head, body = read_kmx(inp)
    moves = json.load(open(moves_json))
    a, b, d = apply_moves(body, moves)
    out = emit(body, a, b, d)
    open(prefix + ".kmx", "w").write("\n".join(head + out) + "\n")
    to_kmb(prefix + ".kmx", prefix + ".kmb")
    print(f"wrote {prefix}.kmx/.kmb: {len(out)} ops ({len(out) - len(body):+d}), moves {len(moves)}")


def relabel(kmb, prefix, Q=None):
    """Smallest-width relabel: every segment (wire non-zero on some input) gets a physical wire, greedy by start."""
    import numpy as np
    import invsqrt_fw_sim as F
    c = F.Circ(kmb)
    s = F.Sim(c, hmr="zero", phase=False, track=True)
    s.run_to(len(c.ops))
    xr = c.regs["x"]
    nz0 = set(xr)
    ev = {}
    for i, w, z in s.events:
        ev.setdefault(i, []).append((w, z))
    # text of the input circuit, op aligned
    subprocess.run(["bash", "-lc", f"{KMPY} -c \"import kickmix as km; open('{prefix}.in.kmx','w').write(str(km.Circuit.from_file('{kmb}', format='kmb')))\""], check=True)
    head, body = read_kmx(prefix + ".in.kmx")
    assert len(body) == len(c.ops), (len(body), len(c.ops))
    phys = {w: w for w in xr}     # logical -> physical for live segments
    free = sorted(set(range(200)) - set(xr))
    live = set(xr)
    out = []
    peak = len(xr)
    for i, l in enumerate(body):
        qs = qubits(l)
        before = set(live)
        for w, z in ev.get(i, []):
            if z:
                live.add(w)
            else:
                live.discard(w)
        used_now = {phys[w] for w in qs if w in phys}
        tmp = {}
        for w in qs:
            if w not in phys:
                # becomes live at this op, or a zero operand: a free physical wire not used by this op
                p = next(f for f in free if f not in used_now and f not in tmp.values())
                tmp[w] = p
        for w, p in tmp.items():
            if w in live:
                phys[w] = p
                free.remove(p)
        m = dict(phys)
        m.update(tmp)
        out.append(QRE.sub(lambda mo: f"q{m[int(mo.group(1))]}", l))
        for w in list(phys):
            if w not in live:
                free.append(phys.pop(w))
                free.sort()
        peak = max(peak, len(live | before))
    # x back to its register wires with free SWAPs
    for k, w in enumerate(xr):
        p = phys[w]
        if p != w:
            # swap physical p and w; whoever logical sits on physical w moves to p
            other = [lw for lw, pp in phys.items() if pp == w]
            out.append(f"SWAP q{p} q{w}")
            if other:
                phys[other[0]] = p
            phys[w] = w
    nq = 1 + max(int(x) for l in out for x in QRE.findall(l))
    open(prefix + ".kmx", "w").write("\n".join(head + out) + "\n")
    to_kmb(prefix + ".kmx", prefix + ".kmb")
    print(f"relabelled: peak live {peak}, physical wires {nq}")


if __name__ == "__main__":
    if sys.argv[1] == "relabel":
        relabel(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else None)
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
