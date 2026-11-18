"""Unit tests for barcode assignment rules and the decompressed budget."""
import gzip
import os
import resource
import time
from pathlib import Path

import pytest

from app import processing


KEYS = [
    processing.SampleKey("s1", b"AAAAAAAACCCCCCCC"),
    processing.SampleKey("s2", b"GGGGGGGGTTTTTTTT"),
]


def test_exact_and_single_substitution_assignment():
    assert processing._best_sample(b"AAAAAAAAT", b"CCCCCCCCG", KEYS).sample_id == "s1"
    # one substitution across both mates combined
    assert processing._best_sample(b"AAAAAAAT", b"CCCCCCCC", KEYS).sample_id == "s1"
    assert processing._best_sample(b"AAAAAAAA", b"CCCCCCCG", KEYS).sample_id == "s1"
    assert processing._best_sample(b"GGGGGGGG", b"TTTTTTTT", KEYS).sample_id == "s2"


def test_two_substitutions_quarantined():
    assert processing._best_sample(b"AAAAAATT", b"CCCCCCCC", KEYS) is None
    assert processing._best_sample(b"AAAAAAAT", b"CCCCCCCG", KEYS) is None


def test_tie_at_distance_one_quarantined():
    # Two barcodes exactly one base apart: a read on the middle base is
    # distance 1 from both -> not unique -> quarantine.
    keys = [
        processing.SampleKey("a", b"AAAAAAAACCCCCCCC"),
        processing.SampleKey("b", b"AAAAAAACCCCCCCCC"),
    ]
    # probe = AAAAAAA C CCCCCCCC? Build explicitly:
    probe1 = b"AAAAAAAC"  # 1 diff vs a's r1, 0 vs b's r1
    probe2 = b"CCCCCCCC"
    # distance to a: 1; distance to b r1 0 + r2 same -> 0: unique b
    assert processing._best_sample(probe1, probe2, keys).sample_id == "b"
    # A genuinely tied probe (distance equal to both).
    tied = [
        processing.SampleKey("a", b"AAAAAAAACCCCCCCC"),
        processing.SampleKey("b", b"CAAAAAAACCCCCCCC"),
    ]
    # probe differing equally:
    probe = b"GAAAAAAA", b"CCCCCCCC"
    assert processing._best_sample(*probe, tied) is None


def test_short_read_quarantined():
    assert processing._best_sample(b"AAAAAAA", b"CCCCCCCC", KEYS) is None


def _write_lane(tmp: Path, pairs: int, seq_len: int = 40):
    r1 = tmp / "R1.fastq.gz"
    r2 = tmp / "R2.fastq.gz"
    seq1 = b"A" * seq_len
    seq2 = b"C" * seq_len
    qual = b"I" * seq_len
    with gzip.open(r1, "wb") as g:
        for i in range(pairs):
            g.write(b"@r%d\n%s\n+\n%s\n" % (i, seq1, qual))
    with gzip.open(r2, "wb") as g:
        for i in range(pairs):
            g.write(b"@r%d\n%s\n+\n%s\n" % (i, seq2, qual))
    return r1, r2


def test_size_budget_enforced(tmp_path):
    # Each record pair is ~2*(1+10+2+40) = ~106 bytes; make it exceed budget.
    r1, r2 = _write_lane(tmp_path, pairs=5, seq_len=40)
    with pytest.raises(processing.SizeExceededError):
        list(processing._paired_records(r1, r2, size_budget=128))


def test_streaming_memory_on_large_lane(tmp_path, data_dir):
    """A lane with ~100+ MiB decompressed stays in small RSS (constant memory)."""
    pairs = 500_000  # ~104 MiB uncompressed, streamed record by record
    r1, r2 = _write_lane(tmp_path, pairs=pairs, seq_len=40)
    assert r1.stat().st_size < 10 * 1024 * 1024  # highly compressible gzip

    rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    count = 0
    t0 = time.time()
    for _ in processing._paired_records(
            r1, r2, size_budget=512 * 1024 * 1024):
        count += 1
        if count == pairs:
            break
    elapsed = time.time() - t0
    rss_after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    assert count == pairs
    # ru_maxrss is KiB on Linux; growth must be a small fraction of input.
    growth_mib = (rss_after - rss_before) / 1024
    print(f"\nstreamed {pairs} pairs ({elapsed:.1f}s), "
          f"RSS growth {growth_mib:.1f} MiB")
    assert growth_mib < 64, f"RSS grew {growth_mib:.1f} MiB during streaming"
