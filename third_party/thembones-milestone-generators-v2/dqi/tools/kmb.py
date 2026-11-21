"""Independent numpy reader for kickmix .kmb circuits (struct-of-arrays binary format).

Gives direct array access to a circuit for analysis and rewriting, which the kickmix
Python API does not expose. Layout per research/05_kickmix_and_techniques.md section A3.
"""
from __future__ import annotations

import dataclasses

import numpy as np

MAGIC = bytes.fromhex("d750c7d5c329d326e3cc9f6834f2b8bf")
ID_MASK = 0x1FFFFFFF

# Base opcodes (condition flag 0x80 stripped).
NEG, POP, INFO_TAG, DEBUG = 0x00, 0x01, 0x06, 0x07
Z, X, R, Z_POW = 0x10, 0x11, 0x12, 0x13
CZ, CX, SWAP = 0x20, 0x21, 0x22
CCZ, CCX, CSWAP = 0x30, 0x31, 0x32
BIT0, BIT1, BITINV, HMR = 0x40, 0x41, 0x42, 0x50
PUSH = 0x81  # full opcode (always carries a condition bit)
NAMES = {0x00: "NEG", 0x01: "POP_CONDITION", 0x06: "INFO_TAG", 0x07: "DEBUG_PRINT", 0x17: "DEBUG_PRINT_Q",
         0x47: "DEBUG_PRINT_B", 0x10: "Z", 0x11: "X", 0x12: "R", 0x13: "Z_POW", 0x20: "CZ", 0x21: "CX",
         0x22: "SWAP", 0x30: "CCZ", 0x31: "CCX", 0x32: "CSWAP", 0x40: "BIT_STORE0", 0x41: "BIT_STORE1",
         0x42: "BIT_INVERT", 0x50: "HMR"}


@dataclasses.dataclass
class Kmb:
    num_qubits: int
    num_bits: int
    registers: list  # (name, ids array, is_qubit array)
    ops: np.ndarray  # u8 per op (with 0x80 condition flag)
    q3: tuple  # three u32 arrays
    q2: tuple  # two u32 arrays
    q1: np.ndarray
    bit_t: np.ndarray
    cond: np.ndarray
    angles: np.ndarray  # (n, 2) u64: w0, w1
    tags: list  # strings, one per INFO_TAG op, in order

    @property
    def base(self) -> np.ndarray:
        """Opcode with the condition flag removed, except PUSH_CONDITION which keeps 0x81."""
        b = self.ops & 0x7F
        b = b.astype(np.uint8)
        b[self.ops == PUSH] = PUSH
        return b

    @property
    def nq(self) -> np.ndarray:
        n = (self.ops & 0x30) >> 4
        n = n.astype(np.uint8)
        return n

    @property
    def is_toffoli(self) -> np.ndarray:
        return (self.ops & 0x7C) == 0x30  # 0x30, 0x31, 0x32 with or without the condition flag

    @property
    def has_cond(self) -> np.ndarray:
        return (self.ops & 0x80) != 0

    def class_index(self, k: int) -> np.ndarray:
        """For every op, its index into the k-qubit operand arrays (valid where nq == k)."""
        return np.cumsum(self.nq == k) - 1

    def tag_positions(self) -> np.ndarray:
        return np.flatnonzero(self.ops == INFO_TAG)

    def histogram(self) -> dict:
        vals, counts = np.unique(self.ops, return_counts=True)
        out = {}
        for v, c in zip(vals.tolist(), counts.tolist()):
            name = "PUSH_CONDITION" if v == PUSH else NAMES.get(v & 0x7F, hex(v)) + ("_IF" if v & 0x80 else "")
            out[name] = c
        return out


def read(path: str, mmap: bool = True) -> Kmb:
    data = np.memmap(path, dtype=np.uint8, mode="r") if mmap else np.fromfile(path, dtype=np.uint8)
    if bytes(data[:16]) != MAGIC:
        raise ValueError("not a .kmb file")
    version, num_qubits, num_bits, num_registers = np.frombuffer(data[16:32], dtype="<u4").tolist()
    if version != 2:
        raise ValueError(f"unsupported kmb version {version}")
    counts = np.frombuffer(data[32:104], dtype="<u8").tolist()
    pos = 104
    packets: dict[int, list] = {}
    n = len(data)
    while pos < n:
        ptype, plen = np.frombuffer(data[pos:pos + 16], dtype="<u8").tolist()
        pos += 16
        packets.setdefault(ptype, []).append((pos, plen))
        pos += plen

    def arr(ptype: int, dtype: str) -> np.ndarray:
        if ptype not in packets:
            return np.zeros(0, dtype=dtype)
        off, plen = packets[ptype][0]
        return np.frombuffer(data, dtype=dtype, count=plen // np.dtype(dtype).itemsize, offset=off)

    registers = []
    for off, plen in packets.get(0, []):
        name_len, reg_len = np.frombuffer(data[off:off + 8], dtype="<u4").tolist()
        name = bytes(data[off + 8:off + 8 + name_len]).decode("utf-8", "replace")
        ids = np.frombuffer(data, dtype="<u4", count=reg_len, offset=off + 8 + name_len)
        registers.append((name, ids & ID_MASK, (ids >> 31) == 1))
    ops = arr(1, "u1")
    tags = []
    if 11 in packets:
        off, plen = packets[11][0]
        ntag = int(np.frombuffer(data[off:off + 8], dtype="<u8")[0])
        ends = np.frombuffer(data, dtype="<u8", count=ntag, offset=off + 8)
        base = off + 8 + 8 * ntag
        start = 0
        for e in ends.tolist():
            tags.append(bytes(data[base + start:base + e]).decode("utf-8", "replace"))
            start = e
    ang = arr(10, "<u8")
    return Kmb(num_qubits, num_bits, registers, ops,
               (arr(2, "<u4"), arr(3, "<u4"), arr(4, "<u4")), (arr(5, "<u4"), arr(6, "<u4")),
               arr(7, "<u4"), arr(8, "<u4"), arr(9, "<u4"), ang.reshape(-1, 2), tags)
