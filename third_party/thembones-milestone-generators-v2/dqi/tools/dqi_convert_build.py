"""Build a DQI decoder from a generator module exposing decode(b, S, E).

Usage: python dqi_convert_build.py <decode_module.py> <out.kmb> [--lists]

Same shape as the public carry-q26-892 build.py: kickmix stands in for the kit's builder module
`circ`, registers S0..S69 (qubits) and E0..E4094 (bits) are named in the file. Writes
<out>.build.json with digest and counts.
"""
import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.expanduser("~/qcc/data/kickmix"))
import kickmix as km

sys.modules["circ"] = km

entry = Path(sys.argv[1]).resolve()
target = Path(sys.argv[2]).resolve()
if "--lists" in sys.argv:
    aq, ab = km.CircuitBuilder.alloc_qubits, km.CircuitBuilder.alloc_bits
    km.CircuitBuilder.alloc_qubits = lambda self, n: list(aq(self, n))
    km.CircuitBuilder.alloc_bits = lambda self, n: list(ab(self, n))
start = time.time()
sys.path.insert(0, str(entry.parent))
spec = importlib.util.spec_from_file_location("decoder_under_test", entry)
module = importlib.util.module_from_spec(spec)
sys.modules["decoder_under_test"] = module
spec.loader.exec_module(module)
print("imported in", round(time.time() - start, 1), "s", flush=True)
b = km.CircuitBuilder()
S = [b.alloc_qubits(12) for _ in range(70)]
E = [b.alloc_bits(12) for _ in range(4095)]
module.decode(b, S, E)
print("decoded in", round(time.time() - start, 1), "s", flush=True)
metadata = []
for i in range(70):
    metadata.append(f'REGISTER r{i} "S{i}"')
    metadata.extend(f"APPEND_TO_REGISTER q{12 * i + k} r{i}" for k in range(12))
for j in range(4095):
    metadata.append(f'REGISTER r{70 + j} "E{j}"')
    metadata.extend(f"APPEND_TO_REGISTER b{12 * j + k} r{70 + j}" for k in range(12))
c = km.Circuit("\n".join(metadata)) + b.finish_circuit()
c.to_file(str(target), format="kmb")
with open(target, "rb") as f:
    digest = hashlib.file_digest(f, "sha256").hexdigest()
record = dict(sha256=digest, Q=c.num_qubits, T_static=c.operation_statistics().num_toffoli_like, num_ops=len(c),
              num_bits=c.num_bits, seconds=round(time.time() - start, 1))
target.with_suffix(".build.json").write_text(json.dumps(record, indent=1) + "\n")
print(json.dumps(record), flush=True)
