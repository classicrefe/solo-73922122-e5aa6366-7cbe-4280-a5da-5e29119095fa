"""Overlap merging of paired-end reads into single-end reads.

For every validated R1/R2 pair the R2 mate is reverse-complemented and the
algorithm searches for the *unique* ungapped overlap between the tail of R1
and the head of R2' of length at least 12 with a mismatch rate of at most 5%.
No merge happens when there is no qualifying position or more than one.

Quality rules inside an overlap:
  * same base: Phred scores are added and capped at 40;
  * different base: the higher-quality base wins; on a tie the base becomes
    ``N`` carrying the lower quality (equal on a tie).

Everything streams record by record, so a lane whose decompressed size reaches
the 512 MiB budget still fits inside a 128 MiB container.
"""
from __future__ import annotations

import errno
import fcntl
import gzip
import json
import os
import shutil
import time
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from . import config, db, processing

# bytes.translate table: complement, preserving N.
_COMPLEMENT = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")

# 2-bit codes used by the packed-integer overlap prefilter. N is coded as A:
# the prefilter therefore only ever under-counts mismatches, so every
# candidate it accepts is re-verified exactly on the raw bytes (which handle
# N correctly).
_CODE_TABLE = bytes(
    (b"ACGT".index(bytes([b])) if b in b"ACGT" else 0) for b in range(256))


class MergeError(Exception):
    """Any condition that fails the whole merge job (bad lane, I/O, ...)."""


@dataclass
class MergeResult:
    merged: bool
    seq: bytes = b""
    qual: bytes = b""
    overlap_len: int = 0


def reverse_complement(seq: bytes) -> bytes:
    return seq.translate(_COMPLEMENT)[::-1]


def _max_mismatches(overlap_len: int) -> int:
    """Largest mismatch count still inside the inclusive 5% budget."""
    return (overlap_len * config.MERGE_MAX_MISMATCH_PCT) // 100


def _pack(seq: bytes) -> int:
    """Pack bases into 2 bits each (base 0 at the least-significant end).

    A single small integer holds a whole read (<=150 bases => <=300 bits).
    """
    code = _CODE_TABLE
    value = 0
    for i in range(len(seq)):
        value |= code[seq[i]] << (2 * i)
    return value


def find_overlap(r1: bytes, r2_rc: bytes) -> int | None:
    """Return the unique qualifying overlap length, else None.

    A candidate of length ``L`` aligns ``r1[-L:]`` with ``r2_rc[:L]``.
    Candidates are ungapped, L in 12..min(len(r1), len(r2_rc)), with a
    mismatch rate <= 5%. A packed 2-bit integer Hamming prefilter (big-int
    operations implemented in C) rejects most candidates; because N is coded
    as A the prefilter can only under-count mismatches, so every survivor is
    re-verified exactly on the raw bytes before acceptance. Multiple
    qualifying positions or none at all mean "do not merge".
    """
    n1 = len(r1)
    max_len = min(n1, len(r2_rc))
    if max_len < config.MERGE_MIN_OVERLAP:
        return None
    p1 = _pack(r1)
    p2 = _pack(r2_rc)
    qualifying: list[int] = []
    for length in range(config.MERGE_MIN_OVERLAP, max_len + 1):
        budget = _max_mismatches(length)
        field_mask = (1 << (2 * length)) - 1
        # Mask after XOR: tail1's bits above 2L are shift padding, not data.
        xor = (((p1 >> (2 * (n1 - length))) ^ p2) & field_mask)
        low = field_mask // 3  # repeating 01 over the packed width
        if (((xor & low) << 1) | (xor & (low << 1))).bit_count() > budget:
            continue
        # Exact verification on the raw bytes (handles N correctly).
        diffs = 0
        for x, y in zip(r1[n1 - length:], r2_rc[:length]):
            if x != y:
                diffs += 1
                if diffs > budget:
                    break
        if diffs <= budget:
            qualifying.append(length)
            if len(qualifying) > 1:
                return None
    return qualifying[0] if len(qualifying) == 1 else None


def _overlap_qualities(q1: int, q2: int) -> int:
    """Combine two Phred scores inside an overlap, capped at MERGE_QUAL_CAP."""
    return min(q1 + q2, config.MERGE_QUAL_CAP)


def merge_pair(r1: processing.Record, r2: processing.Record) -> MergeResult:
    """Merge one validated pair. Returns merged=False when no unique overlap."""
    r2_seq_rc = reverse_complement(r2.seq)
    length = find_overlap(r1.seq, r2_seq_rc)
    if length is None:
        return MergeResult(merged=False)

    # R2 was reverse-complemented; its quality string is simply reversed.
    q2_rev = r2.qual[::-1]
    merged_seq = bytearray()
    merged_qual = bytearray()

    # 1. R1-only prefix.
    merged_seq += r1.seq[:-length]
    merged_qual += r1.qual[:-length]

    # 2. Overlap region, position by position.
    r1_tail = r1.seq[-length:]
    q1_tail = r1.qual[-length:]
    r2_head = r2_seq_rc[:length]
    q2_head = q2_rev[:length]
    for b1, b2, c1, c2 in zip(r1_tail, r2_head, q1_tail, q2_head):
        if b1 == b2:
            merged_seq.append(b1)
            merged_qual.append(33 + _overlap_qualities(c1 - 33, c2 - 33))
        else:
            if c1 > c2:
                merged_seq.append(b1)
                merged_qual.append(c1)
            elif c2 > c1:
                merged_seq.append(b2)
                merged_qual.append(c2)
            else:
                # Tie: ambiguity base N with the lower (equal) quality.
                merged_seq.append(ord("N"))
                merged_qual.append(min(c1, c2))

    # 3. R2' tail not covered by R1.
    merged_seq += r2_seq_rc[length:]
    merged_qual += q2_rev[length:]
    return MergeResult(
        merged=True,
        seq=bytes(merged_seq),
        qual=bytes(merged_qual),
        overlap_len=length,
    )


# ---------------- streaming lane iteration ----------------

def _check_mate_lengths(r1: processing.Record, r2: processing.Record,
                        lane: str, index: int) -> None:
    for mate, rec in (("R1", r1), ("R2", r2)):
        n = len(rec.seq)
        if n < config.MERGE_MIN_READ_LEN or n > config.MERGE_MAX_READ_LEN:
            raise MergeError(
                f"lane {lane}: record {index} {mate} length {n} is outside "
                f"the allowed {config.MERGE_MIN_READ_LEN}.."
                f"{config.MERGE_MAX_READ_LEN} bases")


def iter_lane_pairs(r1_path: Path, r2_path: Path, lane: str
                    ) -> Iterator[tuple[processing.Record, processing.Record]]:
    """Yield length-checked, name-verified R1/R2 pairs of one published lane.

    Reuses the demux FASTQ validator (gzip integrity, FASTQ framing, Phred
    range, name/order equality, equal record counts) and additionally enforces
    the 40..150 base window. The decompressed total is capped at 512 MiB.
    """
    index = 0
    try:
        for r1, r2 in processing._paired_records(
                r1_path, r2_path, config.MAX_TOTAL_UNCOMPRESSED):
            _check_mate_lengths(r1, r2, lane, index)
            yield r1, r2
            index += 1
    except processing.FastqFormatError as exc:
        # Any corrupt/ill-formed lane fails the whole merge job, never a
        # partial file.
        raise MergeError(f"lane {lane}: {exc}") from exc
    except processing.SizeExceededError as exc:
        raise MergeError(f"lane {lane}: {exc}") from exc


# ---------------- job lifecycle ----------------

def _merge_work_dir(job_id: str) -> Path:
    return config.MERGE_WORK_DIR / job_id


def _merge_lock_path(job_id: str) -> Path:
    return _merge_work_dir(job_id) / "lock"


def _acquire_merge_lock(job_id: str):
    wd = _merge_work_dir(job_id)
    wd.mkdir(parents=True, exist_ok=True)
    fh = open(_merge_lock_path(job_id), "w")
    try:
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as exc:
        fh.close()
        if exc.errno in (errno.EACCES, errno.EAGAIN):
            raise RuntimeError(
                f"merge job {job_id} is locked by another instance") from exc
        raise
    fh.write(f"{os.getpid()}\n")
    fh.flush()
    return fh


def _lane_sample_files(conn, lane: str, sample_id: str) -> tuple[Path, Path]:
    row = conn.execute(
        "SELECT r1_file, r2_file FROM lane_counts "
        "WHERE lane=? AND sample_id=?", (lane, sample_id)).fetchone()
    if row is None:
        raise MergeError(
            f"sample {sample_id!r} is absent from published lane {lane!r}")
    r1_path, r2_path = Path(row["r1_file"]), Path(row["r2_file"])
    root = config.PUBLISHED_DIR.resolve()
    for p in (r1_path, r2_path):
        rp = p.resolve()
        if root not in rp.parents or not rp.is_file():
            raise MergeError(f"lane {lane}: published file missing ({p.name})")
    return r1_path, r2_path


def _expected_pairs(conn, lanes: list[str], sample_id: str) -> int:
    """Published total of pairs that will be presented to the merger."""
    total = 0
    found = False
    for lane in lanes:
        row = conn.execute(
            "SELECT read_pairs FROM lane_counts WHERE lane=? AND sample_id=?",
            (lane, sample_id)).fetchone()
        if row is not None:
            found = True
            total += row["read_pairs"]
    if not found:
        raise MergeError(
            f"sample {sample_id!r} exists in none of the snapshot lanes")
    return total


class _MergeWriter:
    """One gzip output, written record by record."""

    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self._raw = open(path, "wb")
        self._gz = gzip.GzipFile(
            filename="", mode="wb", fileobj=self._raw, compresslevel=6)

    def write(self, name: bytes, seq: bytes, qual: bytes) -> None:
        self._gz.write(b"@" + name[1:] + b"\n")
        self._gz.write(seq + b"\n+\n" + qual + b"\n")

    def close(self) -> None:
        self._gz.close()
        self._raw.close()


def process_merge_job(job_row, instance_id: str) -> None:
    """Stream every snapshot lane for one sample into a single merged FASTQ."""
    job_id = job_row["merge_job_id"]
    snapshot_id = job_row["snapshot_id"]
    sample_id = job_row["sample_id"]
    lanes = list(json.loads(job_row["lanes_json"]))  # frozen snapshot order

    lock = _acquire_merge_lock(job_id)
    work = _merge_work_dir(job_id)
    staged = work / "output"
    tmp_publish = config.MERGES_DIR / f".tmp.{job_id}"
    final_publish = config.MERGES_DIR / job_id
    conn = db.get_conn()
    try:
        if staged.exists():
            shutil.rmtree(staged)
        if tmp_publish.exists():
            shutil.rmtree(tmp_publish)
        if final_publish.exists():
            shutil.rmtree(final_publish)

        expected = _expected_pairs(conn, lanes, sample_id)
        out_path = staged / "merged.fastq.gz"
        writer = _MergeWriter(out_path)
        merged_count = 0
        unmerged_count = 0
        seen_pairs = 0
        try:
            for lane in lanes:
                r1_path, r2_path = _lane_sample_files(conn, lane, sample_id)
                for r1, r2 in iter_lane_pairs(r1_path, r2_path, lane):
                    result = merge_pair(r1, r2)
                    if result.merged:
                        writer.write(r1.name, result.seq, result.qual)
                        merged_count += 1
                    else:
                        unmerged_count += 1
                    seen_pairs += 1
                    if seen_pairs % config.HEARTBEAT_FLUSH_EVERY == 0:
                        db.merge_heartbeat(job_id, instance_id)
        finally:
            writer.close()

        # Cross-check against the published counts before anything is visible.
        if seen_pairs != expected:
            raise MergeError(
                f"pair count mismatch: scanned {seen_pairs} but published "
                f"counts say {expected}")
        if merged_count + unmerged_count != expected:
            raise MergeError("merged/unmerged totals do not add up")

        # Verify the staged gzip is a complete, readable stream.
        try:
            with gzip.open(out_path, "rb") as verify:
                while verify.read(65536):
                    pass
        except (OSError, EOFError, zlib.error) as exc:
            raise MergeError(f"merged output is not a valid gzip: {exc}") from exc

        tmp_publish.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(staged), str(tmp_publish))
        os.rename(tmp_publish, final_publish)

        db.complete_merge_job(
            job_id, instance_id,
            str(final_publish / "merged.fastq.gz"),
            merged_count, unmerged_count, expected)
    except Exception:
        # A failed merge must never leave a half file or staging debris.
        # The work directory itself (and its lock) stay until the lock is
        # released; startup recovery removes any directory left by a crash.
        shutil.rmtree(staged, ignore_errors=True)
        shutil.rmtree(tmp_publish, ignore_errors=True)
        shutil.rmtree(final_publish, ignore_errors=True)
        raise
    finally:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        finally:
            lock.close()
        # Success path: all scratch is gone now that the lock is released.
        if db.get_merge_job(job_id) is not None \
                and db.get_merge_job(job_id)["status"] == "completed":
            shutil.rmtree(work, ignore_errors=True)


# ---------------- startup recovery ----------------

def _lock_is_free(path: Path) -> bool:
    if not path.exists():
        return True
    fh = open(path, "r")
    try:
        try:
            fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            if exc.errno in (errno.EACCES, errno.EAGAIN):
                return False
            raise
        fcntl.flock(fh.fileno(), fcntl.LOCK_UN)
        return True
    finally:
        fh.close()


def _remove_merge_scratch(job_id: str) -> None:
    wd = _merge_work_dir(job_id)
    if wd.exists():
        shutil.rmtree(wd, ignore_errors=True)
    tmp = config.MERGES_DIR / f".tmp.{job_id}"
    if tmp.exists():
        shutil.rmtree(tmp, ignore_errors=True)


def recover_merges_on_startup() -> list[str]:
    """Mirror demux recovery: flock protects live work, DB rows protect output."""
    notes: list[str] = []
    conn = db.get_conn()

    for job in conn.execute(
            "SELECT * FROM merge_jobs WHERE status='processing'").fetchall():
        job_id = job["merge_job_id"]
        if not _lock_is_free(_merge_lock_path(job_id)):
            continue  # a live instance owns this merge
        notes.append(f"recovering interrupted merge job {job_id}")
        conn.execute(
            "UPDATE merge_jobs SET status='failed', instance_id=NULL,"
            "error='interrupted (owner process disappeared)', updated_at=? "
            "WHERE merge_job_id=?", (db.utcnow(), job_id))
        _remove_merge_scratch(job_id)

    live = {
        r["merge_job_id"] for r in conn.execute(
            "SELECT merge_job_id FROM merge_jobs "
            "WHERE status='processing'").fetchall()
    }
    if config.MERGE_WORK_DIR.exists():
        for wd in config.MERGE_WORK_DIR.iterdir():
            if not wd.is_dir() or wd.name in live:
                continue
            lock = wd / "lock"
            if lock.exists() and not _lock_is_free(lock):
                continue
            if time.time() - wd.stat().st_mtime < 5:
                continue
            notes.append(f"removing orphan merge work directory {wd.name}")
            shutil.rmtree(wd, ignore_errors=True)

    if config.MERGES_DIR.exists():
        finished = {
            r["merge_job_id"] for r in conn.execute(
                "SELECT merge_job_id FROM merge_jobs WHERE status='completed'"
            ).fetchall()
        }
        for child in config.MERGES_DIR.iterdir():
            if child.name.startswith(".tmp."):
                job_id = child.name[len(".tmp."):]
                job = conn.execute(
                    "SELECT status FROM merge_jobs WHERE merge_job_id=?",
                    (job_id,)).fetchone()
                if job is not None and job["status"] == "processing" \
                        and not _lock_is_free(_merge_lock_path(job_id)):
                    continue
                notes.append(
                    f"removing unfinished merge staging {child.name}")
                shutil.rmtree(child, ignore_errors=True)
            elif child.name not in finished:
                notes.append(
                    f"removing unpublished merge output {child.name}")
                shutil.rmtree(child, ignore_errors=True)
    return notes
