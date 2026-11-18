"""Static reader and fast cost statistics for kickmix `.kmb` circuits (no quantum simulation).

Usage:
    python tools/kmb.py FILE.kmb [FILE2.kmb ...] [--shots 1048576] [--no-mean] [--profile] [--json]

What it gives, per circuit: qubits, bits, instruction counts, static Toffoli-like count (CCX + CCZ + CSWAP),
Z_POW / HMR counts, the **expected number of Toffoli-like gates executed per shot**, and optionally an
allocated-width profile and per-stage Toffolis (stages from INFO_TAG instructions, when the circuit has them).

Why the expected Toffoli count needs no quantum simulation: in kickmix, classical bits are written only by HMR
(a fair coin, independent of the data), BIT_STORE0/1 and BIT_INVERT, and conditions read only classical bits.
So whether a conditioned gate fires is a function of the coins and of the circuit's classical input bits alone.
`mean_toffolis` replays that classical trajectory on many shots at once, 64 shots per uint64 word.
Classical input bits (bits listed in registers) are fair coins by default; this is exact for invsqrt and
chemistry (no classical inputs) but NOT the dqi verifier's distribution of the E registers, so for dqi use
the verifier or kickmix-engine for the cost (`--no-mean` skips it).

The result is a Monte Carlo estimate with its standard error (se); the official verifiers also report a
sampled mean (invsqrt: 3,133,452 shots; dqi: 102,400), so agreement is checked within the combined error.

File layout (kmb version 2): magic, header (version, num_qubits, num_bits, num_registers, num_ops, 8 packet
counts), register packets, then the instruction stream as struct-of-arrays: one opcode byte per instruction and
one operand array per operand kind (3-qubit, 2-qubit, 1-qubit, target bit, condition bit, angle), each consumed
in stream order; INFO_TAG texts last. Layout as documented by GrantPrater's `vqcc_ir.py` (Apache-2.0, dataset
grantprater/vqcc-dqi-61c8348d-merged-closings); this file is our own implementation.

Complexity: reading is O(file size) with memory-mapped arrays (no copy); counts are O(num_ops) vectorized;
`mean_toffolis` is O(E × W) where E = classical events (writers, PUSH/POP, Toffoli-likes) and W = shots / 64
(invsqrt: E ≈ 5k, W = 16,384 → well under a second); the width profile is O(touches log touches).
"""
from __future__ import annotations

import argparse
import json
import struct
import sys
import time
from dataclasses import dataclass, field

import numpy as np

MAGIC = bytes.fromhex("d750c7d5c329d326e3cc9f6834f2b8bf")
HAS_COND = 0x80
BASE_NAMES = {0x00: "NEG", 0x01: "POP_CONDITION", 0x10: "Z", 0x11: "X", 0x12: "R", 0x13: "Z_POW", 0x50: "HMR",
              0x20: "CZ", 0x21: "CX", 0x22: "SWAP", 0x30: "CCZ", 0x31: "CCX", 0x32: "CSWAP",
              0x40: "BIT_STORE0", 0x41: "BIT_STORE1", 0x42: "BIT_INVERT", 0x06: "INFO_TAG",
              0x07: "DEBUG_PRINT_EMPTY", 0x47: "DEBUG_PRINT_C", 0x17: "DEBUG_PRINT_Q"}
PUSH = 0x81
POP = 0x01
TOFFOLI_BASE = (0x30, 0x31, 0x32)
WRITER_BASE = (0x50, 0x40, 0x41, 0x42)
RESET_BASE = (0x50, 0x12)  # HMR and R leave the qubit in |0>
PACKETS = ("qqq0", "qqq1", "qqq2", "qq0", "qq1", "q0", "b0", "bc", "angles")  # payload types 2..10


def op_name(t: int) -> str:
    if t == PUSH:
        return "PUSH_CONDITION"
    return BASE_NAMES.get(t & 0x7F, f"?{t:#x}") + ("_IF" if t & HAS_COND else "")


def _lut(pred) -> np.ndarray:
    return np.array([bool(pred(t)) for t in range(256)])


L3Q = _lut(lambda t: t & 0x30 == 0x30)
L2Q = _lut(lambda t: t & 0x30 == 0x20)
L1Q = _lut(lambda t: t & 0x30 == 0x10)
L1B = _lut(lambda t: t & 0x40)
LCOND = _lut(lambda t: t & HAS_COND)
LANGLE = _lut(lambda t: t & 0x7F == 0x13)
LTAG = _lut(lambda t: t == 0x06)
LTOF = _lut(lambda t: (t & 0x7F) in TOFFOLI_BASE)
LWRITER = _lut(lambda t: (t & 0x7F) in WRITER_BASE and t != PUSH)
LRESET = _lut(lambda t: (t & 0x7F) in RESET_BASE)
LEVENT = LWRITER | LTOF | _lut(lambda t: t in (PUSH, POP))
_PACKET_LUT = {"qqq0": L3Q, "qqq1": L3Q, "qqq2": L3Q, "qq0": L2Q, "qq1": L2Q, "q0": L1Q,
               "b0": L1B, "bc": LCOND, "angles": LANGLE}


@dataclass
class Kmb:
    path: str
    num_qubits: int
    num_bits: int
    types: np.ndarray
    arrays: dict
    registers: list = field(default_factory=list)  # (name, qubit ids, bit ids)
    tags: list = field(default_factory=list)       # (instruction index, text)

    @classmethod
    def read(cls, path: str) -> "Kmb":
        buf = np.memmap(path, dtype=np.uint8, mode="r")
        if bytes(buf[:16]) != MAGIC:
            raise ValueError(f"{path}: not a .kmb file")
        version, nq, nb, nreg, nops = struct.unpack("<IIIIQ", bytes(buf[16:40]))
        if version != 2:
            raise ValueError(f"{path}: kmb version {version} not supported")
        header_counts = struct.unpack("<8Q", bytes(buf[40:104]))
        pos = 104
        regs = []
        for _ in range(nreg):
            ptype, plen = struct.unpack("<QQ", bytes(buf[pos:pos + 16]))
            if ptype != 0:
                raise ValueError(f"{path}: register packet expected, got {ptype}")
            name_len, n_ids = struct.unpack("<II", bytes(buf[pos + 16:pos + 24]))
            name = bytes(buf[pos + 24:pos + 24 + name_len]).decode()
            ids = np.frombuffer(buf, dtype="<u4", count=n_ids, offset=pos + 24 + name_len)
            qubit = ids >> 31 == 1
            regs.append((name, (ids[qubit] & 0x7FFFFFFF).astype(np.int64), ids[~qubit].astype(np.int64)))
            pos += 16 + plen

        def packet(expect: int):
            nonlocal pos
            ptype, plen = struct.unpack("<QQ", bytes(buf[pos:pos + 16]))
            if ptype != expect:
                raise ValueError(f"{path}: packet {ptype}, expected {expect}")
            off = pos + 16
            pos += 16 + plen
            return off, plen

        off, plen = packet(1)
        types = np.frombuffer(buf, dtype=np.uint8, count=nops, offset=off)
        arrays = {}
        for k, name in enumerate(PACKETS):
            off, plen = packet(2 + k)
            if name == "angles":
                arrays[name] = np.frombuffer(buf, dtype="<u8", count=plen // 8, offset=off).reshape(-1, 2)
            else:
                arrays[name] = np.frombuffer(buf, dtype="<u4", count=plen // 4, offset=off)
        off, plen = packet(11)
        n_tags = struct.unpack("<Q", bytes(buf[off:off + 8]))[0]
        ends = np.frombuffer(buf, dtype="<u8", count=n_tags, offset=off + 8)
        chars = bytes(buf[off + 8 + 8 * n_tags:off + plen])
        tag_idx = np.flatnonzero(LTAG[types])
        starts = np.concatenate(([0], ends[:-1])) if n_tags else []
        tags = [(int(i), chars[int(a):int(b)].decode(errors="replace")) for i, a, b in zip(tag_idx, starts, ends)]
        c = cls(path, nq, nb, types, arrays, regs, tags)
        c._validate(header_counts)
        return c

    def _validate(self, header_counts) -> None:
        bc = np.bincount(self.types, minlength=256)
        for name, lut in _PACKET_LUT.items():
            if len(self.arrays[name]) != int(bc[lut].sum()):
                raise ValueError(f"{self.path}: {name} array length does not match the opcodes")
        computed = (int(bc[L3Q].sum()), int(bc[L2Q].sum()), int(bc[L1Q].sum()), int(bc[L1B].sum()),
                    int(bc[LCOND].sum()), int(bc[LANGLE].sum()), int(bc[LTAG].sum()))
        if tuple(header_counts[:7]) != computed:
            raise ValueError(f"{self.path}: header counts {header_counts[:7]} != computed {computed}")

    @property
    def num_ops(self) -> int:
        return len(self.types)

    def static_counts(self) -> dict:
        bc = np.bincount(self.types, minlength=256)
        return {op_name(t): int(n) for t, n in enumerate(bc) if n}

    def toffoli_like(self) -> int:
        return int(LTOF[self.types].sum())

    def classical_inputs(self) -> set:
        return {int(b) for _, _, bits in self.registers for b in bits}

    def _events(self):
        """Classical events in stream order: (opcode, condition bit or -1, target bit or -1). O(num_ops)."""
        t = self.types
        ev = np.flatnonzero(LEVENT[t])
        kind = t[ev]
        cpos = np.cumsum(LCOND[t], dtype=np.int64)[ev] - 1
        bpos = np.cumsum(L1B[t], dtype=np.int64)[ev] - 1
        bc, b0 = self.arrays["bc"].astype(np.int64), self.arrays["b0"].astype(np.int64)
        cbit = np.where(LCOND[kind], bc[np.clip(cpos, 0, max(len(bc) - 1, 0))] if len(bc) else 0, -1)
        tbit = np.where(L1B[kind], b0[np.clip(bpos, 0, max(len(b0) - 1, 0))] if len(b0) else 0, -1)
        return ev, kind.tolist(), cbit.astype(np.int64).tolist(), tbit.astype(np.int64).tolist()

    def mean_toffolis(self, shots: int = 1 << 20, seed: int = 1, per_event: bool = False):
        """Expected Toffoli-like gates executed per shot, by replaying the classical trajectory bit-parallel.

        Returns (mean, standard error, fire_counts or None). fire_counts[i] = number of the `shots` shots on
        which the i-th Toffoli-like of the stream fired (per_event=True). Cost O(events × shots/64)."""
        W = max(shots // 64, 2)
        rng = np.random.default_rng(seed)
        free = self.classical_inputs()
        ones = np.full(W, np.uint64(0xFFFFFFFFFFFFFFFF))
        zero = np.zeros(W, dtype=np.uint64)
        bits: dict[int, np.ndarray] = {}
        stack = [ones]
        tally = np.zeros(W, dtype=np.int64)
        fires = [] if per_event else None

        def coin():
            return rng.integers(0, 1 << 64, size=W, dtype=np.uint64, endpoint=False)

        def get(b):
            v = bits.get(b)
            if v is None:
                v = bits[b] = coin() if b in free else zero
            return v

        _, kinds, cbits, tbits = self._events()
        for op, cb, tb in zip(kinds, cbits, tbits):
            if op == PUSH:
                stack.append(stack[-1] & get(cb))
                continue
            if op == POP:
                stack.pop()
                continue
            g = stack[-1] & get(cb) if cb >= 0 else stack[-1]
            base = op & 0x7F
            if base in TOFFOLI_BASE:
                f = np.bitwise_count(g)
                tally += f
                if per_event:
                    fires.append(int(f.sum()))
            elif base == 0x50:   # HMR: fresh coin where the gate fires
                old = get(tb)
                bits[tb] = (coin() & g) | (old & ~g)
            elif base == 0x40:   # BIT_STORE0
                bits[tb] = get(tb) & ~g
            elif base == 0x41:   # BIT_STORE1
                bits[tb] = get(tb) | g
            else:                # BIT_INVERT
                bits[tb] = get(tb) ^ g
        per_word = tally / 64.0
        return float(per_word.mean()), float(per_word.std(ddof=1) / np.sqrt(W)), fires

    def touches(self):
        """(instruction index, qubit, is_reset) for every qubit operand, in stream order. O(num_ops)."""
        t = self.types
        parts = []
        for lut, names in ((L3Q, ("qqq0", "qqq1", "qqq2")), (L2Q, ("qq0", "qq1")), (L1Q, ("q0",))):
            idx = np.flatnonzero(lut[t])
            for name in names:
                parts.append((idx, self.arrays[name].astype(np.int64), LRESET[t[idx]] & (name == "q0")))
        op = np.concatenate([p[0] for p in parts])
        q = np.concatenate([p[1] for p in parts])
        reset = np.concatenate([p[2] for p in parts])
        order = np.argsort(op, kind="stable")
        return op[order], q[order], reset[order]

    def width_profile(self):
        """Allocated width along the stream: a qubit counts from its first non-reset touch after a reset (or
        from the start, for register qubits) until the HMR/R that resets it. An upper bound on the qubits that
        can be non-zero. Returns (width per instruction as int32, peak). O(touches log touches)."""
        op, q, reset = self.touches()
        n = self.num_ops
        delta = np.zeros(n + 1, dtype=np.int64)
        order = np.lexsort((op, q))  # by qubit, then stream position
        oq, oo, orr = q[order], op[order], reset[order]
        new_q = np.r_[True, oq[1:] != oq[:-1]]
        prev_reset = np.r_[True, orr[:-1]] | new_q  # a touch after a reset (or the first one) opens an interval
        reg_q = {int(x) for _, qs, _ in self.registers for x in qs}
        is_reg = np.isin(oq, np.fromiter(reg_q, dtype=np.int64, count=len(reg_q))) if reg_q else np.zeros(len(oq), bool)
        opens = prev_reset & ~orr & ~(new_q & is_reg)
        np.add.at(delta, oo[opens], 1)
        np.add.at(delta, oo[orr & ~new_q & ~prev_reset] + 1, -1)  # a reset closes the open interval
        width = np.cumsum(delta[:n]) + len(reg_q)
        return width.astype(np.int32), int(width.max()) if n else 0

    def stages(self):
        """Instruction ranges between consecutive INFO_TAGs: [(start, end, text of the opening tag)]."""
        cuts = [(i, txt) for i, txt in self.tags]
        if not cuts or cuts[0][0] != 0:
            cuts = [(0, "(start)")] + cuts
        return [(a, b, txt) for (a, txt), (b, _) in zip(cuts, cuts[1:] + [(self.num_ops, "")])]


def summarize(path: str, shots: int, mean: bool, profile: bool) -> dict:
    t0 = time.perf_counter()
    c = Kmb.read(path)
    counts = c.static_counts()
    out = {"file": path, "num_qubits": c.num_qubits, "num_bits": c.num_bits, "num_ops": c.num_ops,
           "toffoli_like_static": c.toffoli_like(), "z_pow": counts.get("Z_POW", 0) + counts.get("Z_POW_IF", 0),
           "hmr": counts.get("HMR", 0) + counts.get("HMR_IF", 0), "info_tags": len(c.tags),
           "registers": len(c.registers), "classical_inputs": len(c.classical_inputs()), "counts": counts}
    if mean:
        m, se, fires = c.mean_toffolis(shots, per_event=profile)
        out.update(mean_toffolis=m, mean_toffolis_se=se, shots=shots,
                   mean_note="classical inputs drawn as fair coins" if out["classical_inputs"] else "exact model")
    if profile:
        width, peak = c.width_profile()
        out["allocated_width_peak"] = peak
        if c.tags:
            tof_pos = np.flatnonzero(LTOF[c.types])
            per = np.asarray(fires, dtype=np.float64) / max(shots // 64, 2) / 64 if mean else None
            rows = []
            for a, b, txt in c.stages():
                sel = (tof_pos >= a) & (tof_pos < b)
                rows.append({"stage": txt[:60], "ops": b - a, "toffoli_static": int(sel.sum()),
                             "toffoli_mean": float(per[sel].sum()) if per is not None else None,
                             "width_max": int(width[a:b].max()) if b > a else 0})
            out["stages"] = rows
    out["seconds"] = round(time.perf_counter() - t0, 3)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--shots", type=int, default=1 << 20)
    ap.add_argument("--no-mean", action="store_true", help="skip the expected Toffoli count (use for dqi)")
    ap.add_argument("--profile", action="store_true", help="allocated-width peak and per-stage Toffolis")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    for f in a.files:
        s = summarize(f, a.shots, not a.no_mean, a.profile)
        if a.json:
            print(json.dumps(s))
            continue
        line = f"{f}: Q={s['num_qubits']} bits={s['num_bits']} ops={s['num_ops']:,} toffoli_static={s['toffoli_like_static']:,}"
        if "mean_toffolis" in s:
            line += f" mean={s['mean_toffolis']:.4f}±{s['mean_toffolis_se']:.4f}"
        if "allocated_width_peak" in s:
            line += f" width_peak={s['allocated_width_peak']}"
        print(line + f" z_pow={s['z_pow']} hmr={s['hmr']} tags={s['info_tags']} [{s['seconds']} s]")
        for r in s.get("stages", []):
            print(f"    {r}")


if __name__ == "__main__":
    sys.exit(main())
