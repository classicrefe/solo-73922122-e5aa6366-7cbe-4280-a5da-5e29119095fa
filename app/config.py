"""Runtime configuration driven by environment variables.

All state lives under one shared data directory so two container instances
mounting the same volume can cooperate.
"""
from __future__ import annotations

import os
from pathlib import Path


def _default_data_dir() -> Path:
    candidate = Path(os.environ.get("DEMUX_DATA_DIR", "/data"))
    try:
        candidate.mkdir(parents=True, exist_ok=True)
        probe = candidate / ".write_test"
        probe.touch()
        probe.unlink()
        return candidate.resolve()
    except OSError:
        # Fallback for non-container local runs/tests.
        fallback = Path(os.environ.get("TMPDIR", "/tmp")) / "demux-data"
        fallback.mkdir(parents=True, exist_ok=True)
        return fallback.resolve()


DATA_DIR = _default_data_dir()
DB_PATH = DATA_DIR / "demux.db"
INCOMING_DIR = DATA_DIR / "incoming"
WORK_DIR = DATA_DIR / "work"
PUBLISHED_DIR = DATA_DIR / "published"
# Overlap-merge jobs get their own scratch/output trees beside the demux ones.
MERGE_WORK_DIR = WORK_DIR / "merges"
MERGES_DIR = PUBLISHED_DIR / "merges"

# Hard limits from the task specification.
MAX_SAMPLES = 8
BARCODE_LEN = 8
MAX_TOTAL_UNCOMPRESSED = 512 * 1024 * 1024  # 512 MiB across both mates
HEARTBEAT_FLUSH_EVERY = 256  # records between heartbeat writes
POLL_INTERVAL = float(os.environ.get("DEMUX_POLL_INTERVAL", "1.0"))

# Overlap-merge rules.
MERGE_MIN_READ_LEN = 40          # mates outside 40..150 bases fail the lane
MERGE_MAX_READ_LEN = 150
MERGE_MIN_OVERLAP = 12           # shortest accepted ungapped overlap
MERGE_MAX_MISMATCH_PCT = 5       # mismatch rate ceiling (percent, inclusive)
MERGE_QUAL_CAP = 40              # Phred cap when summing matching-base qualities

for _d in (DATA_DIR, INCOMING_DIR, WORK_DIR, PUBLISHED_DIR,
           MERGE_WORK_DIR, MERGES_DIR):
    _d.mkdir(parents=True, exist_ok=True)
