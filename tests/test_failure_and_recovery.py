"""Corruption, retry, snapshot isolation and crash recovery tests."""
import gzip
import io
import json
import time
import zlib
from pathlib import Path

from fastapi.testclient import TestClient

BARCODES = [
    {"sample_id": "s1", "barcode_r1": "AAAAAAAA", "barcode_r2": "CCCCCCCC"},
]


def _post(client, lane, barcodes, r1_path, r2_path):
    with open(r1_path, "rb") as f1, open(r2_path, "rb") as f2:
        return client.post(
            "/jobs",
            data={"lane": lane, "barcodes": json.dumps(barcodes)},
            files={"r1": ("R1.fastq.gz", f1, "application/gzip"),
                   "r2": ("R2.fastq.gz", f2, "application/gzip")},
        )


def _wait(client, job_id, timeout=15):
    deadline = time.time() + timeout
    while time.time() < deadline:
        row = client.get(f"/jobs/{job_id}").json()
        if row["status"] in ("published", "failed"):
            return row
        time.sleep(0.1)
    raise AssertionError("timeout waiting for job")


def _good_pair(tmp: Path):
    r1 = tmp / "R1.fastq.gz"
    r2 = tmp / "R2.fastq.gz"
    with gzip.open(r1, "wb") as g:
        g.write(b"@x\nAAAAAAAATT\n+\nIIIIIIIIII\n")
    with gzip.open(r2, "wb") as g:
        g.write(b"@x\nCCCCCCCCGG\n+\nJJJJJJJJJJ\n")
    return r1, r2


def test_corrupt_gzip_fails_whole_lane(client, tmp_path):
    r1, r2 = _good_pair(tmp_path)
    r1.write_bytes(r1.read_bytes()[:-6])  # truncated gzip footer/stream
    resp = _post(client, "badGz", BARCODES, r1, r2)
    job = _wait(client, resp.json()["job"]["job_id"])
    assert job["status"] == "failed"
    assert "format" in job["error"].lower() or "gzip" in job["error"].lower()
    assert client.get("/lanes/badGz").status_code == 404


def test_mismatched_names_fail(client, tmp_path):
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    with gzip.open(r1, "wb") as g:
        g.write(b"@x\nAAAAAAAATT\n+\nIIIIIIIIII\n")
    with gzip.open(r2, "wb") as g:
        g.write(b"@y\nCCCCCCCCGG\n+\nJJJJJJJJJJ\n")  # different read name
    resp = _post(client, "badNames", BARCODES, r1, r2)
    job = _wait(client, resp.json()["job"]["job_id"])
    assert job["status"] == "failed"
    assert "name" in job["error"]


def test_unequal_record_counts_fail(client, tmp_path):
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    with gzip.open(r1, "wb") as g:
        g.write(b"@x\nAAAAAAAATT\n+\nIIIIIIIIII\n"
                b"@y\nAAAAAAAATT\n+\nIIIIIIIIII\n")
    with gzip.open(r2, "wb") as g:
        g.write(b"@x\nCCCCCCCCGG\n+\nJJJJJJJJJJ\n")
    resp = _post(client, "badCount", BARCODES, r1, r2)
    job = _wait(client, resp.json()["job"]["job_id"])
    assert job["status"] == "failed"
    assert "different numbers" in job["error"]


def test_retry_failed_job_succeeds(client, tmp_path):
    r1, r2 = _good_pair(tmp_path)
    bad_r1 = tmp_path / "badR1.fastq.gz"
    bad_r1.write_bytes(r1.read_bytes()[:-6])
    resp = _post(client, "retryLane", BARCODES, bad_r1, r2)
    job_id = resp.json()["job"]["job_id"]
    failed = _wait(client, job_id)
    assert failed["status"] == "failed"

    # Retry endpoint refuses non-failed states and requeues failed ones.
    assert client.post(f"/jobs/{job_id}/retry").status_code == 200
    # Requeue once more would now be a conflict (already queued).
    assert client.post(f"/jobs/{job_id}/retry").status_code == 409
    bad_again = _wait(client, job_id)
    assert bad_again["status"] == "failed"  # input still corrupt

    # Repair the retained upload, then retry: same lane/input identity is
    # preserved (retry, not a new conflicting submission).
    from app import db as dbmod
    raw = dbmod.get_job(job_id)
    Path(raw["r1_incoming"]).write_bytes(r1.read_bytes())
    assert client.post(f"/jobs/{job_id}/retry").status_code == 200
    job = _wait(client, job_id)
    assert job["status"] == "published"
    assert client.get("/lanes/retryLane").status_code == 200


def test_snapshot_excludes_later_lanes(client, helpers, tmp_path):
    # Lane A published, snapshot taken, lane B afterwards must not appear.
    r1 = [(b"a", b"AAAAAAAATT", b"I" * 10)]
    r2 = [(b"a", b"CCCCCCCCGG", b"J" * 10)]
    d1 = tmp_path / "a"; d1.mkdir()
    p1, p2 = helpers.make_pair(r1, r2, d1)
    resp = _post(client, "laneA", BARCODES, p1, p2)
    _wait(client, resp.json()["job"]["job_id"])

    snap = client.post("/snapshots").json()
    assert snap["lanes"] == ["laneA"]

    d2 = tmp_path / "b"; d2.mkdir()
    q1, q2 = helpers.make_pair(r1, r2, d2)
    resp = _post(client, "laneB", BARCODES, q1, q2)
    _wait(client, resp.json()["job"]["job_id"])

    same = client.get(f"/snapshots/{snap['snapshot_id']}").json()
    assert same["lanes"] == ["laneA"]
    headers_resp = client.get(
        f"/snapshots/{snap['snapshot_id']}/samples/s1/r1")
    assert headers_resp.headers["x-snapshot-lanes"] == "laneA"
    # A sample that exists in no snapshot lane yields 404.
    assert client.get(
        f"/snapshots/{snap['snapshot_id']}/samples/ghost/r1").status_code == 404
    assert client.get("/snapshots/nope").status_code == 404


def test_recovery_clears_interrupted_remnants(data_dir, helpers, tmp_path):
    import importlib
    from app import db as dbmod, processing, main as mainmod

    # Simulate a crashed processing job: DB row + scratch work dir + lock,
    # no live holder.
    r1, r2 = _good_pair(tmp_path)
    upload = data_dir / "incoming" / "deadjob"
    upload.mkdir(parents=True)
    up1, up2 = upload / "R1.fastq.gz", upload / "R2.fastq.gz"
    up1.write_bytes(r1.read_bytes())
    up2.write_bytes(r2.read_bytes())

    dbmod.create_job(
        {"job_id": "deadjob", "lane": "deadLane", "input_hash": "x",
         "barcode_table": "[]", "r1_incoming": str(up1),
         "r2_incoming": str(up2), "created_at": dbmod.utcnow()},
        [{"sample_id": "s1", "barcode_r1": "AAAAAAAA",
          "barcode_r2": "CCCCCCCC"}],
    )
    conn = dbmod.get_conn()
    conn.execute(
        "UPDATE jobs SET status='processing', instance_id='ghost', "
        "heartbeat_at=? WHERE job_id='deadjob'",
        ("2000-01-01T00:00:00+00:00",))
    work = data_dir / "work" / "deadjob"
    work.mkdir(parents=True)
    (work / "lock").write_text("ghost\n")
    (work / "output").mkdir()

    # An already-published lane must be left untouched by recovery.
    pub = data_dir / "published" / "lanes" / "liveLane"
    pub.mkdir(parents=True)
    (pub / "keep.txt").write_text("durable")
    conn.execute(
        "INSERT INTO lane_publishes (lane, job_id, published_at, seq) "
        "VALUES ('liveLane','j1',?,1)", (dbmod.utcnow(),))

    notes = processing.recover_on_startup()
    assert any("deadjob" in n for n in notes)
    assert not work.exists()
    assert (pub / "keep.txt").exists()
    row = conn.execute(
        "SELECT status FROM jobs WHERE job_id='deadjob'").fetchone()
    assert row["status"] == "failed"

    # A second recovery run is a no-op (idempotent).
    assert processing.recover_on_startup() == []


def test_recovery_never_removes_live_instance_work(data_dir):
    import fcntl
    from app import db as dbmod, processing

    work = data_dir / "work" / "livejob"
    work.mkdir(parents=True)
    lock_fh = open(work / "lock", "w")
    fcntl.flock(lock_fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    try:
        dbmod.create_job(
            {"job_id": "livejob", "lane": "liveLane2", "input_hash": "x",
             "barcode_table": "[]",
             "r1_incoming": str(data_dir / "incoming/x/R1.fastq.gz"),
             "r2_incoming": str(data_dir / "incoming/x/R2.fastq.gz"),
             "created_at": dbmod.utcnow()},
            [{"sample_id": "s1", "barcode_r1": "AAAAAAAA",
              "barcode_r2": "CCCCCCCC"}],
        )
        dbmod.get_conn().execute(
            "UPDATE jobs SET status='processing', instance_id='live', "
            "heartbeat_at=? WHERE job_id='livejob'",
            ("2000-01-01T00:00:00+00:00",))
        notes = processing.recover_on_startup()
        assert work.exists()  # flock owner protects it
        assert not any("livejob" in n for n in notes)
    finally:
        fcntl.flock(lock_fh.fileno(), fcntl.LOCK_UN)
        lock_fh.close()
