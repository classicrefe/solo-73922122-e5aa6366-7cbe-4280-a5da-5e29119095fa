"""Shared pytest fixtures: isolated data dir and in-process API client."""
import gzip
import os
import sys
import tempfile
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


@pytest.fixture()
def data_dir(tmp_path, monkeypatch):
    d = tmp_path / "data"
    d.mkdir()
    monkeypatch.setenv("DEMUX_DATA_DIR", str(d))
    # Reload configuration before any app import of the paths.
    import importlib
    from app import config as cfg
    importlib.reload(cfg)
    for sub in (cfg.INCOMING_DIR, cfg.WORK_DIR, cfg.PUBLISHED_DIR):
        sub.mkdir(parents=True, exist_ok=True)
    from app import db as dbmod
    importlib.reload(dbmod)
    dbmod.init_db()
    yield d


@pytest.fixture()
def client(data_dir):
    # Import after env/reload so paths point at the tmp data dir.
    from fastapi.testclient import TestClient
    from app import main as mainmod
    with TestClient(mainmod.app) as c:
        yield c


def write_fastq_gz(path: Path, records):
    """records: list of (name_no_at, seq, qual)."""
    with gzip.open(path, "wb") as g:
        for name, seq, qual in records:
            g.write(b"@" + name + b"\n")
            g.write(seq + b"\n")
            g.write(b"+\n")
            g.write(qual + b"\n")


def make_pair(r1_records, r2_records, tmp: Path,
              name1="R1.fastq.gz", name2="R2.fastq.gz"):
    p1 = tmp / name1
    p2 = tmp / name2
    write_fastq_gz(p1, r1_records)
    write_fastq_gz(p2, r2_records)
    return p1, p2


def q_for(seq: bytes, ch=b"I") -> bytes:
    return ch * len(seq)


@pytest.fixture()
def helpers():
    return type("H", (), {
        "write_fastq_gz": staticmethod(write_fastq_gz),
        "make_pair": staticmethod(make_pair),
        "q_for": staticmethod(q_for),
    })
