#!/usr/bin/env python3
"""Generate a small demo lane (paired gzip FASTQ) for the /docs workflow.

Usage:
    python scripts/make_demo_fastq.py OUT_DIR [pairs]

Writes OUT_DIR/R1.fastq.gz and OUT_DIR/R2.fastq.gz containing a mix of:
  * exact barcode hits for two samples,
  * one-substitution hits (still assignable, unique),
  * a two-substitution read and a random read (quarantine).
"""
import gzip
import sys
from pathlib import Path

SAMPLES = [
    ("sampleA", "AAAAAAAA", "CCCCCCCC"),
    ("sampleB", "GGGGGGGG", "TTTTTTTT"),
]


def main(out_dir: str, pairs: int = 60) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    r1_records: list[tuple[bytes, bytes]] = []
    r2_records: list[tuple[bytes, bytes]] = []
    idx = 0

    def add(b1: bytes, b2: bytes) -> None:
        nonlocal idx
        insert1 = b"ACGTACGT"
        insert2 = b"TGCATGCA"
        s1, s2 = b1 + insert1, b2 + insert2
        r1_records.append((f"read{idx}".encode(), s1, b"I" * len(s1)))
        r2_records.append((f"read{idx}".encode(), s2, b"J" * len(s2)))
        idx += 1

    per = max(1, pairs // 4)
    for _ in range(per):  # exact hits
        add(SAMPLES[0][1].encode(), SAMPLES[0][2].encode())
        add(SAMPLES[1][1].encode(), SAMPLES[1][2].encode())
    # one substitution on R2 -> still sampleA, unique
    add(b"AAAAAAAA", b"CCCGCCCC")
    # one substitution on R1 -> still sampleB, unique
    add(b"GGGGGGGT", b"TTTTTTTT")
    # two substitutions -> quarantine
    add(b"AAAAAATT", b"CCCCCCCG")
    # random-ish -> quarantine
    add(b"NNNNNNNN", b"ACGTACGT")

    for mate, records, name in (
            (1, r1_records, "R1.fastq.gz"), (2, r2_records, "R2.fastq.gz")):
        with gzip.open(out / name, "wb") as g:
            for rname, seq, qual in records:
                g.write(b"@" + rname + b" mate" + str(mate).encode() + b"\n")
                g.write(seq + b"\n+\n" + qual + b"\n")
    print(f"wrote {len(r1_records)} paired records to {out}")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 60)
