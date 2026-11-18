"""Build submission.zip from a circuit folder and run the checks the Kaggle metric runs.

Usage (WSL, venv with kickmix):
    python tools/package.py SRC_DIR OUT_ZIP          # check, then write the zip
    python tools/package.py --check-only SRC_DIR     # checks only
    python tools/package.py --self-test SRC_DIR      # checks that planted defects are rejected

SRC_DIR must hold invsqrt.kmb|kmx, dqi.kmb|kmx and dfthcsos/{angle_k_0..n-1.kmb|kmx, submission.npz}.

Checks (forum topic 751951 and the competition overview): exactly one invsqrt, one dqi and one
dfthcsos/ folder and no stray circuit file anywhere; angle_k_0..n-1 without gaps; every circuit has
max(num_qubits, num_bits) - 1 < 699,050; every register is named; ANGLE circuits use at most 128
qubits; submission.npz holds U, W, H2_I, H1_I, Hsym_lower as little-endian floats and tree_edges as
an integer array. The zip is deterministic (sorted members, fixed timestamps) and re-read after
writing to confirm every member's SHA-256.

Cost: O(total bytes) for hashing and compression plus one kickmix load per circuit (the 455 MB dqi
circuit loads in about 5 s). Memory: one circuit at a time.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import sys
import tempfile
import time
import zipfile

import numpy as np

MAX_QUBIT_OR_BIT_ID = 2**30 // 1536  # 699,050: the metric requires max(num_qubits, num_bits) - 1 < this
ANGLE_MAX_QUBITS = 128
CIRCUIT_EXT = (".kmb", ".kmx")
NPZ_FLOAT_KEYS = ("U", "W", "H2_I", "H1_I", "Hsym_lower")
ZIP_TIME = (1980, 1, 1, 0, 0, 0)


class CheckError(Exception):
    pass


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def collect_members(src: str) -> dict[str, str]:
    """Map zip member name -> source path, enforcing the file layout. O(number of files)."""
    members: dict[str, str] = {}
    stray = []
    for root, _, files in os.walk(src):
        for f in files:
            rel = os.path.relpath(os.path.join(root, f), src).replace(os.sep, "/")
            if rel.endswith(CIRCUIT_EXT) or rel.endswith(".npz"):
                members[rel] = os.path.join(root, f)
            else:
                stray.append(rel)
    if stray:
        raise CheckError(f"unexpected non-circuit files (they would be zipped): {sorted(stray)[:5]}")
    for task in ("invsqrt", "dqi"):
        hits = [m for m in members if os.path.basename(m).split(".")[0] == task]
        if hits != [f"{task}.kmb"] and hits != [f"{task}.kmx"]:
            raise CheckError(f"need exactly one top-level {task}.kmb or {task}.kmx, found {hits}")
    angle = {}
    for m in members:
        if m in ("invsqrt.kmb", "invsqrt.kmx", "dqi.kmb", "dqi.kmx", "dfthcsos/submission.npz"):
            continue
        g = re.fullmatch(r"dfthcsos/angle_k_(\d+)\.(kmb|kmx)", m)
        if not g:
            raise CheckError(f"stray file in the submission: {m}")
        k = int(g.group(1))
        if k in angle:
            raise CheckError(f"angle_k_{k} present twice")
        angle[k] = m
    if "dfthcsos/submission.npz" not in members:
        raise CheckError("missing dfthcsos/submission.npz")
    if not angle or sorted(angle) != list(range(len(angle))):
        raise CheckError(f"angle_k indices must be 0..n-1 without gaps, got {sorted(angle)}")
    return members


def check_npz(path: str) -> None:
    with np.load(path, allow_pickle=False) as z:
        for key in NPZ_FLOAT_KEYS + ("tree_edges",):
            if key not in z.files:
                raise CheckError(f"submission.npz lacks {key}")
        for key in NPZ_FLOAT_KEYS:
            a = z[key]
            if a.dtype.kind != "f":
                raise CheckError(f"submission.npz {key} must be a float array, got {a.dtype}")
            if a.dtype.byteorder == ">" or (a.dtype.byteorder == "=" and sys.byteorder != "little"):
                raise CheckError(f"submission.npz {key} must be little-endian, got {a.dtype.str}")
        t = z["tree_edges"]
        if t.dtype.kind not in "iu":
            raise CheckError(f"submission.npz tree_edges must be integer, got {t.dtype}")


def check_circuit(member: str, path: str) -> str:
    import kickmix as km

    c = km.Circuit.from_file(path)
    top = max(c.num_qubits, c.num_bits) - 1
    if top >= MAX_QUBIT_OR_BIT_ID:
        raise CheckError(f"{member}: max(num_qubits, num_bits) - 1 = {top} >= {MAX_QUBIT_OR_BIT_ID}")
    names = list(c.register_data.keys())
    if len(names) != c.num_registers or any(not n for n in names):
        raise CheckError(f"{member}: {c.num_registers} registers but {len(names)} named ones")
    if member.startswith("dfthcsos/") and c.num_qubits > ANGLE_MAX_QUBITS:
        raise CheckError(f"{member}: {c.num_qubits} qubits > {ANGLE_MAX_QUBITS}")
    return f"{c.num_qubits} qubits, {c.num_bits} bits, {c.num_registers} registers"


def check(src: str, verbose: bool = True) -> dict[str, str]:
    members = collect_members(src)
    check_npz(members["dfthcsos/submission.npz"])
    for m in sorted(members):
        if m.endswith(CIRCUIT_EXT):
            info = check_circuit(m, members[m])
            if verbose:
                print(f"  ok  {m:28s} {info}")
    return members


def write_zip(members: dict[str, str], out_zip: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(out_zip)), exist_ok=True)
    tmp = out_zip + ".tmp"
    with zipfile.ZipFile(tmp, "w") as z:
        for m in sorted(members):
            info = zipfile.ZipInfo(m, date_time=ZIP_TIME)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            with open(members[m], "rb") as f, z.open(info, "w", force_zip64=True) as g:
                shutil.copyfileobj(f, g, 1 << 22)
    with zipfile.ZipFile(tmp) as z:  # re-read: every member must match its source byte for byte
        for m in sorted(members):
            h = hashlib.sha256()
            with z.open(m) as f:
                for block in iter(lambda: f.read(1 << 22), b""):
                    h.update(block)
            want = sha256_file(members[m])
            if h.hexdigest() != want:
                raise CheckError(f"zip member {m} does not match its source")
            print(f"  {want}  {m}")
    os.replace(tmp, out_zip)


def self_test(src: str) -> None:
    """Plant one defect at a time in a light copy of SRC (symlinks for big files) and expect rejection."""
    members = collect_members(src)
    n_angle = sum(1 for m in members if "angle_k_" in m)

    def make_copy(d: str) -> None:
        os.makedirs(os.path.join(d, "dfthcsos"))
        for m, p in members.items():
            os.symlink(os.path.abspath(p), os.path.join(d, m))

    cases = {
        "gap in angle_k": lambda d: os.rename(
            os.path.join(d, f"dfthcsos/angle_k_{n_angle - 1}.kmb"), os.path.join(d, f"dfthcsos/angle_k_{n_angle}.kmb")),
        "leftover backup copy": lambda d: os.symlink(os.path.abspath(members["dqi.kmb"]), os.path.join(d, "dfthcsos/dqi.kmb")),
        "big-endian npz": lambda d: _write_big_endian_npz(members["dfthcsos/submission.npz"], os.path.join(d, "dfthcsos/submission.npz")),
        "missing invsqrt": lambda d: os.remove(os.path.join(d, "invsqrt.kmb")),
    }
    failures = 0
    for name, plant in cases.items():
        with tempfile.TemporaryDirectory() as d:
            make_copy(d)
            plant(d)
            try:
                collect_members(d)
                check_npz(os.path.join(d, "dfthcsos/submission.npz"))
                print(f"  FAIL  planted '{name}' was accepted")
                failures += 1
            except CheckError as e:
                print(f"  ok    planted '{name}' rejected: {e}")
    if failures:
        raise SystemExit(f"self-test: {failures} planted defects accepted")
    print("self-test PASS")


def _write_big_endian_npz(src: str, dst: str) -> None:
    os.remove(dst)
    with np.load(src) as z:
        arrays = {k: (z[k].astype(z[k].dtype.newbyteorder(">")) if k == "W" else z[k]) for k in z.files}
    np.savez(dst, **arrays)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src")
    ap.add_argument("out_zip", nargs="?")
    ap.add_argument("--check-only", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    t0 = time.perf_counter()
    try:
        if a.self_test:
            self_test(a.src)
        else:
            print("checks:")
            members = check(a.src)
            if not a.check_only:
                if not a.out_zip:
                    ap.error("OUT_ZIP is required unless --check-only")
                print(f"writing {a.out_zip}:")
                write_zip(members, a.out_zip)
                print(f"zip sha256 {sha256_file(a.out_zip)}  size {os.path.getsize(a.out_zip):,} bytes")
    except CheckError as e:
        raise SystemExit(f"REJECTED: {e}")
    print(f"done in {time.perf_counter() - t0:.1f} s")


if __name__ == "__main__":
    main()
