"""Crash-recovery behaviour specific to overlap-merge jobs."""
import gzip
import random
import time

from tests.test_merge_flow import (
    BARCODES, _build_merging_lane, _post_job, _publish_lane, _wait_merge,
)


def test_interrupted_merge_is_failed_and_scratch_cleared(client, tmp_path):
    from app import db as dbmod, merging, config as cfg

    pairs = _build_merging_lane(tmp_path, "Lrec", n_merge=2, n_unmerged=1,
                                seed=21, sub="rec")
    _publish_lane(client, tmp_path, "Lrec", pairs, "rec")
    snap = client.post("/snapshots").json()
    mid = client.post(
        "/merge/jobs",
        json={"snapshot_id": snap["snapshot_id"], "sample_id": "s1"}
    ).json()["job"]["merge_job_id"]
    _wait_merge(client, mid)  # let it finish first to obtain a known job row

    # Recreate the state a hard kill would leave behind.
    conn = dbmod.get_conn()
    conn.execute(
        "UPDATE merge_jobs SET status='processing', instance_id='ghost', "
        "heartbeat_at=? WHERE merge_job_id=?",
        ("2000-01-01T00:00:00+00:00", mid))
    work = cfg.MERGE_WORK_DIR / mid
    work.mkdir(parents=True, exist_ok=True)
    (work / "lock").write_text("ghost\n")
    (work / "output").mkdir(exist_ok=True)
    staged_gz = work / "output" / "merged.fastq.gz"
    with gzip.open(staged_gz, "wb") as g:
        g.write(b"@half\nACGT\n+\nIIII\n")
    tmp_pub = cfg.MERGES_DIR / f".tmp.{mid}"
    tmp_pub.mkdir(parents=True, exist_ok=True)

    notes = merging.recover_merges_on_startup()
    assert any(mid in n for n in notes)
    assert not work.exists()
    assert not tmp_pub.exists()
    row = conn.execute(
        "SELECT status, error FROM merge_jobs WHERE merge_job_id=?",
        (mid,)).fetchone()
    assert row["status"] == "failed"
    assert "interrupted" in row["error"]

    # A second recovery run is a no-op.
    assert merging.recover_merges_on_startup() == []


def test_recovery_never_touches_live_merge(client, tmp_path):
    import fcntl
    from app import db as dbmod, merging, config as cfg

    pairs = _build_merging_lane(tmp_path, "Llive", n_merge=1, n_unmerged=0,
                                seed=22, sub="live")
    _publish_lane(client, tmp_path, "Llive", pairs, "live")
    snap = client.post("/snapshots").json()
    mid = client.post(
        "/merge/jobs",
        json={"snapshot_id": snap["snapshot_id"], "sample_id": "s1"}
    ).json()["job"]["merge_job_id"]

    work = cfg.MERGE_WORK_DIR / mid
    work.mkdir(parents=True, exist_ok=True)
    lock_fh = open(work / "lock", "w")
    fcntl.flock(lock_fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    try:
        dbmod.get_conn().execute(
            "UPDATE merge_jobs SET status='processing', instance_id='live', "
            "heartbeat_at=? WHERE merge_job_id=?",
            ("2000-01-01T00:00:00+00:00", mid))
        notes = merging.recover_merges_on_startup()
        assert work.exists()           # flock owner protects it
        assert not any(mid in n for n in notes)
        row = dbmod.get_conn().execute(
            "SELECT status FROM merge_jobs WHERE merge_job_id=?",
            (mid,)).fetchone()
        assert row["status"] == "processing"
    finally:
        fcntl.flock(lock_fh.fileno(), fcntl.LOCK_UN)
        lock_fh.close()


def test_completed_merge_output_survives_recovery(client, tmp_path):
    from app import merging, config as cfg

    pairs = _build_merging_lane(tmp_path, "Lkeep", n_merge=2, n_unmerged=0,
                                seed=23, sub="keep")
    _publish_lane(client, tmp_path, "Lkeep", pairs, "keep")
    snap = client.post("/snapshots").json()
    mid = client.post(
        "/merge/jobs",
        json={"snapshot_id": snap["snapshot_id"], "sample_id": "s1"}
    ).json()["job"]["merge_job_id"]
    done = _wait_merge(client, mid)
    assert done["status"] == "completed"
    out = cfg.MERGES_DIR / mid / "merged.fastq.gz"
    assert out.is_file()

    # Recovery must leave completed output untouched.
    assert merging.recover_merges_on_startup() == []
    assert out.is_file()
    assert client.get(f"/merge/jobs/{mid}/download").status_code == 200
