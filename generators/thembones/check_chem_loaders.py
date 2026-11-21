#!/usr/bin/env python3
# New file added by the team (Maxime Bendavid), 2026-11-21, plan step S2.26.
"""Rebuild ThemBones' 15 chemistry loaders from their .kmx sources and compare with the set's SHA-256 manifest.

    ~/venvs/tb313/bin/python generators/thembones/check_chem_loaders.py   (inside WSL, after setup_sandbox.sh)

The loaders' .kmx and the table (submission.npz e05fb18f...) are carried as data by ThemBones (randomized searches and
a GPU fit, not replayable); this checks the part that is rebuildable: kmx -> kmb through kickmix's writer.
Time: 0.6 s for the 15 loaders.
"""
import hashlib
import os
import sys
import time

R = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DS = os.path.join(R, "third_party", "thembones-milestone-generators-v2")
sys.path.insert(0, os.path.join(R, "google-verifiable-quantum-circuit-compilation", "kickmix"))
import kickmix as km  # noqa: E402

want = {}
for line in open(os.path.join(DS, "SET_MANIFEST.sha256")):
    h, p = line.split()
    want[p.lstrip("./")] = h
t0, ok = time.perf_counter(), 0
for k in range(15):
    c = km.Circuit.from_file(os.path.join(DS, "chem", "loaders_kmx", f"angle_k_{k}.kmx"), format="kmx")
    out = f"/tmp/tb_angle_k_{k}.kmb"
    c.to_file(out, format="kmb")
    got = hashlib.sha256(open(out, "rb").read()).hexdigest()
    ok += got == want[f"dfthcsos/angle_k_{k}.kmb"]
    print(f"angle_k_{k}: {got[:16]} {'MATCH' if got == want[f'dfthcsos/angle_k_{k}.kmb'] else 'DIFF'} ({c.num_qubits} qubits)")
print(f"{ok}/15 byte-identical in {time.perf_counter() - t0:.1f} s")
sys.exit(0 if ok == 15 else 1)
