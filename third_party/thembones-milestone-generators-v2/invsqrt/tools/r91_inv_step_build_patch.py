#!/usr/bin/env python3
"""R91-INV-C2-SPLITDEBT: build a split-debt patch of fc18_M (kmx text -> kmb).

A patch for debt b (spec JSON, a list of entries {"b", "h_op", "t_op", "g_terms", "h_terms"}):
  - DELETE every op conditioned on b (" if b<b>"), and every PUSH_CONDITION if b<b> ... POP_CONDITION span with its
    contents (b's own repair: its Toffoli-like gates and the Cliffords of its measured temporaries);
  - INSERT right after op h_op (b's HMR) one Clifford under b per g term: [] -> NEG, [w] -> Z qw, [w, v] -> CZ qw qv;
  - INSERT right after op t_op the same for the h terms.
Ops are numbered as kickmix's instruction index (lines after the 23 register lines, as the census tools count).
Checks: b's HMR is at h_op; the op at t_op is a top-level unconditioned Toffoli-like gate; b is used nowhere else.
Usage: python r91_inv_step_build_patch.py SRC.kmx SPEC.json OUTDIR  (writes OUTDIR/invsqrt.kmx, invsqrt.kmb, patch.json)
"""
import json
import os
import sys

sys.path.insert(0, "<HOME>/qcc/data/kickmix")
sys.path.insert(0, os.path.expanduser("~/qcc/data/kickmix"))
import kickmix as km  # noqa: E402


def term_op(term, b):
    if len(term) == 0:
        return "NEG if b%d" % b
    if len(term) == 1:
        return "Z q%d if b%d" % (term[0], b)
    return "CZ q%d q%d if b%d" % (term[0], term[1], b)


def main():
    src, specp, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    spec = json.load(open(specp))
    lines = open(src).read().split("\n")
    while lines and lines[-1] == "":
        lines.pop()
    head = [i for i, s in enumerate(lines) if s.startswith("APPEND_TO_REGISTER") or s.startswith("REGISTER")]
    nhead = len(head)
    assert head == list(range(nhead)), "register lines must come first"
    ops = lines[nhead:]
    delete = set()
    after = {}
    report = []
    for e in spec:
        b, h_op, t_op = e["b"], e["h_op"], e["t_op"]
        tag = " if b%d" % b
        assert ops[h_op].startswith("HMR ") and ops[h_op].endswith(" b%d" % b), ops[h_op]
        # depth and condition of t_op
        depth = 0
        for i in range(t_op + 1):
            if ops[i].startswith("PUSH_CONDITION"):
                depth += 1
            elif ops[i].startswith("POP_CONDITION"):
                depth -= 1
        assert depth == 0 and " if " not in ops[t_op] and ops[t_op].split()[0] in ("CCX", "CCZ", "CSWAP"), ops[t_op]
        dels = []
        i = 0
        while i < len(ops):
            s = ops[i]
            if s == "PUSH_CONDITION" + tag:
                j, d = i, 0
                while True:
                    if ops[j].startswith("PUSH_CONDITION"):
                        d += 1
                    elif ops[j].startswith("POP_CONDITION"):
                        d -= 1
                    dels.append(j)
                    if d == 0:
                        break
                    j += 1
                i = j + 1
                continue
            if s.endswith(tag):
                assert i > h_op, "b%d used before its HMR at op %d" % (b, i)
                dels.append(i)
            i += 1
        # b must not be read by anything else (BIT_INVERT etc.)
        for i, s in enumerate(ops):
            if i in dels or i == h_op:
                continue
            assert (" b%d " % b) not in (s + " ") or s.startswith("HMR"), "b%d also used at op %d: %s" % (b, i, s)
        delete.update(dels)
        after.setdefault(h_op, []).extend(term_op(t, b) for t in e["g_terms"])
        after.setdefault(t_op, []).extend(term_op(t, b) for t in e["h_terms"])
        report.append({"b": b, "deleted_ops": [[i, ops[i]] for i in dels],
                       "inserted_at_H": after[h_op], "inserted_at_t": after[t_op], "h_op": h_op, "t_op": t_op,
                       "t_gate": ops[t_op]})
        print("b%d: delete %d ops %s; insert %d after op %d, %d after op %d (%s)" % (
            b, len(dels), [ops[i] for i in dels], len(e["g_terms"]), h_op, len(e["h_terms"]), t_op, ops[t_op]), flush=True)
    new_ops = []
    for i, s in enumerate(ops):
        if i not in delete:
            new_ops.append(s)
        new_ops.extend(after.get(i, []))
    text = "\n".join(lines[:nhead] + new_ops) + "\n"
    os.makedirs(outdir, exist_ok=True)
    open(os.path.join(outdir, "invsqrt.kmx"), "w").write(text)
    c = km.Circuit(text)
    c.to_file(os.path.join(outdir, "invsqrt.kmb"), format="kmb")
    c2 = km.Circuit.from_file(os.path.join(outdir, "invsqrt.kmb"), format="kmb")
    assert c2 == c
    print("ops", len(ops), "->", len(new_ops), "qubits", c.num_qubits, "bits", c.num_bits, flush=True)
    json.dump({"src": src, "spec": spec, "report": report, "ops_before": len(ops), "ops_after": len(new_ops)},
              open(os.path.join(outdir, "patch.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
