"""midprize_nb_reroll.py -- builder NB (Nov 18): DQI re-roll tickets as a pure function of (file A, label).

A re-roll is file A plus ONE trailing INFO_TAG op (opcode 0x06) whose tag text is 'reroll:<label>'; no gate changes.
This module rebuilds a ticket from file A, or strips a ticket back to file A, by patching the .kmb bytes directly
(standard library only: no kickmix, no numpy, no team module), so the public notebook can carry it verbatim.

.kmb layout (kickmix format version 2): 16-byte magic | <IIII version, qubits, bits, registers | <9Q counts
(ops, 3-qubit, 2-qubit, 1-qubit, bit targets, conditions, angles, tags, tag bytes) | packets '<QQ type, length'
+ payload. Packet 1 holds one byte per op; packet 11 holds '<Q n' + n u64 end offsets + the tag bytes.
An INFO_TAG op has no qubit, bit or condition operand, so only the op count, the tag count, the tag-byte count and
packets 1 and 11 change.

CLI (python -I):
  strip  <ticket.kmb> <fileA.kmb>               -> writes file A, prints the removed label and both sha256
  append <fileA.kmb> <label> <out.kmb>          -> writes the ticket, prints its sha256
  check  <fileA.kmb> <label>=<sha256> [...]     -> rebuilds each ticket in memory, compares sha256, prints timings
"""
import hashlib
import struct
import sys
import time

# ---- BEGIN NOTEBOOK CODE (build_notebook.py copies the lines between these markers into the notebook) ----
KMB_MAGIC = bytes.fromhex("d750c7d5c329d326e3cc9f6834f2b8bf")
INFO_TAG_OP = 0x06
REROLL_PREFIX = "reroll:"


def kmb_layout(buf):
    """Header fields and the packet list [(type, payload offset, payload length)] of a .kmb held in `buf`."""
    mv = memoryview(buf)
    if bytes(mv[:16]) != KMB_MAGIC:
        raise ValueError("not a .kmb file (magic)")
    head = list(struct.unpack_from("<IIII", mv, 16))
    counts = list(struct.unpack_from("<9Q", mv, 32))
    if head[0] != 2:
        raise ValueError(f"unsupported .kmb version {head[0]}")
    packets, pos, n = [], 104, len(mv)
    while pos < n:
        ptype, plen = struct.unpack_from("<QQ", mv, pos)
        packets.append((ptype, pos + 16, plen))
        pos += 16 + plen
    if pos != n:
        raise ValueError("truncated .kmb (packet lengths overrun the file)")
    kinds = [p[0] for p in packets]
    if kinds.count(1) != 1 or kinds.count(11) != 1:
        raise ValueError(f"expected exactly one ops packet and one tags packet, found types {kinds}")
    return head, counts, packets


def _tags_packet(mv, off, plen):
    ntag = struct.unpack_from("<Q", mv, off)[0]
    ends = list(struct.unpack_from(f"<{ntag}Q", mv, off + 8))
    data = bytes(mv[off + 8 + 8 * ntag: off + plen])
    if (ends[-1] if ends else 0) != len(data):
        raise ValueError("tags packet: end offsets do not match its data")
    return ends, data


def _emit(head, counts, packets, mv, ops_extra, tags_payload, sink):
    """Write header + packets to sink(bytes-like); packet 1 gets `ops_extra` appended (b'' or one op byte, or None to
    drop its last byte); packet 11 is replaced by `tags_payload`."""
    sink(KMB_MAGIC)
    sink(struct.pack("<IIII", *head))
    sink(struct.pack("<9Q", *counts))
    for ptype, off, plen in packets:
        if ptype == 1:
            body = mv[off: off + plen - 1] if ops_extra is None else mv[off: off + plen]
            extra = b"" if ops_extra is None else ops_extra
            sink(struct.pack("<QQ", 1, len(body) + len(extra)))
            sink(body)
            if extra:
                sink(extra)
        elif ptype == 11:
            sink(struct.pack("<QQ", 11, len(tags_payload)))
            sink(tags_payload)
        else:
            sink(mv[off - 16: off + plen])


def reroll_from_file_a(buf_a, label, out_path=None):
    """Ticket = file A + INFO_TAG 'reroll:<label>'. Returns its sha256 (hex); writes it to out_path if given."""
    mv = memoryview(buf_a)
    head, counts, packets = kmb_layout(buf_a)
    tag = (REROLL_PREFIX + label).encode("utf-8")
    _, toff, tlen = [p for p in packets if p[0] == 11][0]
    ends, data = _tags_packet(mv, toff, tlen)
    ends.append((ends[-1] if ends else 0) + len(tag))
    tags_payload = struct.pack("<Q", len(ends)) + struct.pack(f"<{len(ends)}Q", *ends) + data + tag
    counts = list(counts)
    counts[0] += 1          # ops
    counts[7] += 1          # tags
    counts[8] += len(tag)   # tag bytes
    h = hashlib.sha256()
    fh = open(out_path, "wb") if out_path else None
    try:
        def sink(b):
            h.update(b)
            if fh:
                fh.write(b)
        _emit(head, counts, packets, mv, bytes([INFO_TAG_OP]), tags_payload, sink)
    finally:
        if fh:
            fh.close()
    return h.hexdigest()


def file_a_from_reroll(buf_t, out_path=None, collect=None):
    """Inverse: drop the trailing INFO_TAG 'reroll:<label>'. Returns (label, sha256 of file A); writes file A to out_path
    if given; appends its byte chunks to the list `collect` if given."""
    mv = memoryview(buf_t)
    head, counts, packets = kmb_layout(buf_t)
    _, ooff, olen = [p for p in packets if p[0] == 1][0]
    if mv[ooff + olen - 1] != INFO_TAG_OP:
        raise ValueError("the last op is not an INFO_TAG")
    _, toff, tlen = [p for p in packets if p[0] == 11][0]
    ends, data = _tags_packet(mv, toff, tlen)
    start = ends[-2] if len(ends) > 1 else 0
    last = data[start:].decode("utf-8")
    if not last.startswith(REROLL_PREFIX):
        raise ValueError(f"the last tag is not a re-roll tag: {last!r}")
    tag = last.encode("utf-8")
    ends = ends[:-1]
    tags_payload = struct.pack("<Q", len(ends)) + (struct.pack(f"<{len(ends)}Q", *ends) if ends else b"") + data[:start]
    counts = list(counts)
    counts[0] -= 1
    counts[7] -= 1
    counts[8] -= len(tag)
    h = hashlib.sha256()
    fh = open(out_path, "wb") if out_path else None
    try:
        def sink(b):
            h.update(b)
            if fh:
                fh.write(b)
            if collect is not None:
                collect.append(bytes(b))
        _emit(head, counts, packets, mv, None, tags_payload, sink)
    finally:
        if fh:
            fh.close()
    return last[len(REROLL_PREFIX):], h.hexdigest()
# ---- END NOTEBOOK CODE ----


def _read(path):
    with open(path, "rb") as fh:
        return fh.read()


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 23), b""):
            h.update(b)
    return h.hexdigest()


def main(argv):
    cmd = argv[1]
    if cmd == "strip":
        t0 = time.time()
        label, sha_a = file_a_from_reroll(_read(argv[2]), argv[3])
        print(f"label {label!r}  file A sha256 {sha_a}  ticket sha256 {_sha(argv[2])}  ({time.time() - t0:.1f} s)")
    elif cmd == "append":
        t0 = time.time()
        sha = reroll_from_file_a(_read(argv[2]), argv[3], argv[4])
        print(f"ticket {argv[3]!r} sha256 {sha}  ({time.time() - t0:.1f} s)")
    elif cmd == "check":
        t0 = time.time()
        buf = _read(argv[2])
        print(f"file A read {len(buf) / 2**20:.1f} MiB ({time.time() - t0:.1f} s), sha256 {hashlib.sha256(buf).hexdigest()}")
        bad = 0
        for item in argv[3:]:
            label, want = item.split("=")
            t1 = time.time()
            got = reroll_from_file_a(buf, label)
            ok = got == want
            bad += not ok
            print(f"  {label}: {got[:16]} want {want[:16]} {'EQUAL' if ok else 'DIFFERENT'} ({time.time() - t1:.1f} s)")
        print(f"check: {len(argv) - 3 - bad}/{len(argv) - 3} equal, total {time.time() - t0:.1f} s")
        return 1 if bad else 0
    else:
        raise SystemExit(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
