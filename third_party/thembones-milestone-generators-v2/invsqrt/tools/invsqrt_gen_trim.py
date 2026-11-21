"""invsqrt_gen: adapt a lookup program made for a block spec WITH its tied outputs to the spec of the same block
without them (builder option tin: the builder writes those outputs itself). Own tool.

usage: invsqrt_gen_trim.py <program.json> <spec.json> <out.json>
The steps whose target is an output that the new spec does not have are dropped; if a later step reads such a
wire the program cannot be adapted and nothing is written. The result carries the new spec's hash; it is
checked like every program when the builder replays it (every output against its truth table).
"""
import json
import sys

prog, spec, out = (json.load(open(sys.argv[1])), json.load(open(sys.argv[2])), sys.argv[3])
keep = {o["name"] for o in spec["outputs"]} | {i["name"] for i in spec["inputs"]}
steps, gone = [], set()
for st in prog["steps"]:
    tgt = st[1]
    used = {n for f in st[2:] for n in f if n != "1"}
    if used & gone:
        sys.exit(f"cannot adapt: step {st[:2]} reads {sorted(used & gone)}, which the new spec does not have")
    if tgt in keep or (tgt.startswith("s") and tgt[1:].isdigit()):
        steps.append(st)
    else:
        gone.add(tgt)
ns = len({s[1] for s in steps if s[1].startswith("s") and s[1][1:].isdigit()})
assert ns <= spec["max_scratch"], f"needs {ns} scratch wires, the spec allows {spec['max_scratch']}"
new = dict(block=spec["block"], spec_sha=spec["spec_sha"], n_scratch=ns, steps=steps,
           note=f"adapted from a program for spec {prog.get('spec_sha')}: outputs {sorted(gone)} left to the builder")
json.dump(new, open(out, "w"))
print(f"{spec['block']}: {sum(1 for s in steps if s[0] == 'and')} ANDs, {ns} scratch, dropped outputs {sorted(gone)}")
