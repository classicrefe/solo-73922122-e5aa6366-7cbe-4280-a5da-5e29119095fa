#!/usr/bin/env python3
"""Generate a demo lane whose sample reads overlap and can be merged.

The output is shaped exactly like ``make_demo_fastq.py`` (paired gzip FASTQ
with the two demo barcodes) so it can be submitted through POST /jobs, and
after demultiplexing the sampleA bucket contains pairs engineered for the
overlap-merge API: R2' (the reverse complement of R2) shares a unique 30-base
ungapped overlap with the tail of every R1. A few non-overlapping pairs are
included so merged and unmerged counts are both visible.

Usage:
    python scripts/make_merge_demo.py OUT_DIR [pairs]
"""
import gzip
import random
import sys
from pathlib import Path

SAMPLES = [
    ("sampleA", "AAAAAAAA", "CCCCCCCC"),
    ("sampleB", "GGGGGGGG", "TTTTTTTT"),
]

_OVERLAP = 30
_INSERT = 50          # post-barcode insert length (inside the 40..150 window)


def _overlap_pair(rng: random.Random, name: bytes):
    """R2' shares a unique 30-base overlap with R1's tail, with a T tail."""
    shared = bytes(rng.choice(b"ACGT") for _ in range(_OVERLAP))
    r1_insert = bytes(rng.choice(b"ACGT")
                      for _ in range(_INSERT - _OVERLAP)) + shared
    r2_rc = shared + b"T" * (_INSERT - _OVERLAP)
    r2_insert = r2_rc.translate(
        bytes.maketrans(b"ACGT", b"TGCA"))[::-1]
    return ((name, r1_insert, b"I" * _INSERT),
            (name, r2_insert, b"J" * _INSERT))


def _nonoverlap_pair(rng: random.Random, name: bytes):
    r1 = bytes(rng.choice(b"AC") for _ in range(_INSERT))
    r2 = bytes(rng.choice(b"GT") for _ in range(_INSERT))
    return (name, r1, b"I" * _INSERT), (name, r2, b"J" * _INSERT)


def main(out_dir: str, pairs: int = 60) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(2026)
    r1_records: list[tuple[bytes, bytes, bytes]] = []
    r2_records: list[tuple[bytes, bytes, bytes]] = []
    idx = 0

    def add(b1: str, b2: str, insert1: bytes, insert2: bytes,
            q1: bytes, q2: bytes) -> None:
        nonlocal idx
        b1b, b2b = b1.encode(), b2.encode()
        r1_records.append(
            (f"read{idx}".encode(), b1b + insert1, b"I" * len(b1b) + q1))
        r2_records.append(
            (f"read{idx}".encode(), b2b + insert2, b"J" * len(b2b) + q2))
        idx += 1

    # Majority sampleA: mostly mergeable, a couple deliberately non-overlapping.
    n_merge = max(1, int(pairs * 0.8))
    for i in range(n_merge):
        (_, i1, q1), (_, i2, q2) = _overlap_pair(rng, f"m{i}".encode())
        add(SAMPLES[0][1], SAMPLES[0][2], i1, i2, q1, q2)
    for i in range(max(1, pairs - n_merge)):
        (_, i1, q1), (_, i2, q2) = _nonoverlap_pair(rng, f"u{i}".encode())
        add(SAMPLES[0][1], SAMPLES[0][2], i1, i2, q1, q2)
    # A handful of sampleB reads so the lane also has a second sample bucket.
    for i in range(max(1, pairs // 10)):
        (_, i1, q1), (_, i2, q2) = _overlap_pair(rng, f"b{i}".encode())
        add(SAMPLES[1][1], SAMPLES[1][2], i1, i2, q1, q2)

    for mate, records, fname in (
            (1, r1_records, "R1.fastq.gz"), (2, r2_records, "R2.fastq.gz")):
        with gzip.open(out / fname, "wb") as g:
            for rname, seq, qual in records:
                g.write(b"@" + rname + b" mate" + str(mate).encode() + b"\n")
                g.write(seq + b"\n+\n" + qual + b"\n")
    print(f"wrote {len(r1_records)} paired records to {out} "
          f"({n_merge} mergeable sampleA pairs)")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 60)
