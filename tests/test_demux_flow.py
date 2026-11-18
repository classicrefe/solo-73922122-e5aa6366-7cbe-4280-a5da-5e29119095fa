"""End-to-end happy path: assignment, trimming, quarantine, counts, files."""
import gzip
import time

from fastapi.testclient import TestClient


BARCODES = [
    {"sample_id": "s1", "barcode_r1": "AAAAAAAA", "barcode_r2": "CCCCCCCC"},
    {"sample_id": "s2", "barcode_r1": "GGGGGGGG", "barcode_r2": "TTTTTTTT"},
]


def _post(client, lane, barcodes, r1_path, r2_path):
    with open(r1_path, "rb") as f1, open(r2_path, "rb") as f2:
        return client.post(
            "/jobs",
            data={"lane": lane, "barcodes": __import__("json").dumps(barcodes)},
            files={"r1": ("R1.fastq.gz", f1, "application/gzip"),
                   "r2": ("R2.fastq.gz", f2, "application/gzip")},
        )


def _wait_published(client, job_id, timeout=15):
    deadline = time.time() + timeout
    while time.time() < deadline:
        r = client.get(f"/jobs/{job_id}")
        status = r.json()["status"]
        if status in ("published", "failed"):
            return r.json()
        time.sleep(0.1)
    raise AssertionError("job did not finish in time")


def test_full_demux_flow(client, helpers, tmp_path):
    # s1 exact; s1 with one substitution on R2; s2 exact; unknown -> quarantine;
    # ambiguous (equidistant to s1 and s2) -> quarantine.
    r1 = [
        (b"r1", b"AAAAAAAATTTT", b"I" * 12),                 # s1 exact
        (b"r2", b"AAAAAAAATTTT", b"I" * 12),                 # s1 1 sub on R2
        (b"r3", b"GGGGGGGGAAAA", b"I" * 12),                 # s2 exact
        (b"r4", b"NNNNNNNNAAAA", b"I" * 12),                 # unknown
        (b"r5", b"AGAGAGAGAAAA", b"I" * 12),                 # ambiguous
    ]
    r2 = [
        (b"r1", b"CCCCCCCCGGGG", b"J" * 12),
        (b"r2", b"CCCGCCCCGGGG", b"J" * 12),  # 1 substitution
        (b"r3", b"TTTTTTTTGGGG", b"J" * 12),
        (b"r4", b"ACGTACGTGGGG", b"J" * 12),
        (b"r5", b"CTCTCTCTGGGG", b"J" * 12),  # 4 diffs vs each -> quarantine
    ]
    p1, p2 = helpers.make_pair(r1, r2, tmp_path)
    resp = _post(client, "lane1", BARCODES, p1, p2)
    assert resp.status_code == 201, resp.text
    job = resp.json()["job"]
    assert job["status"] == "queued"

    done = _wait_published(client, job["job_id"])
    assert done["status"] == "published", done

    lane = client.get("/lanes/lane1").json()
    assert lane["total_read_pairs"] == 5
    by_sample = {c["sample_id"]: c["read_pairs"] for c in lane["counts"]}
    assert by_sample["s1"] == 2
    assert by_sample["s2"] == 1
    assert by_sample["__quarantine__"] == 2

    # All expected files are downloadable and contain trimmed records.
    r = client.get("/lanes/lane1/samples/s1/r1")
    assert r.status_code == 200
    with gzip.open(__import__("io").BytesIO(r.content)) as g:
        lines = g.read().splitlines()
    assert lines[1] == b"TTTT"  # 8-base barcode trimmed from both seq...
    assert lines[3] == b"IIII"  # ...and the matching quality chars

    r = client.get("/lanes/lane1/samples/s1/r2")
    with gzip.open(__import__("io").BytesIO(r.content)) as g:
        lines = g.read().splitlines()
    assert lines[1] == b"GGGG"
    assert lines[3] == b"JJJJ"

    # Quarantine keeps the full, untrimmed records.
    r = client.get("/lanes/lane1/quarantine/r1")
    with gzip.open(__import__("io").BytesIO(r.content)) as g:
        lines = g.read().splitlines()
    assert b"NNNNNNNNAAAA" in lines
    assert len(lines) == 8  # 2 quarantined pairs x 4 lines

    # Snapshot fixes the published lane set.
    snap = client.post("/snapshots").json()
    assert snap["lanes"] == ["lane1"]
    r = client.get(f"/snapshots/{snap['snapshot_id']}/samples/s1/r1")
    assert r.status_code == 200 and len(r.content) > 0
    r = client.get(f"/snapshots/{snap['snapshot_id']}/quarantine/r2")
    assert r.status_code == 200


def test_duplicate_submission_reused_and_conflicting_rejected(
        client, helpers, tmp_path):
    r1 = [(b"a", b"AAAAAAAATTTT", b"I" * 12)]
    r2 = [(b"a", b"CCCCCCCCGGGG", b"J" * 12)]
    p1, p2 = helpers.make_pair(r1, r2, tmp_path)
    first = _post(client, "dupLane", BARCODES, p1, p2)
    assert first.status_code == 201
    job_id = first.json()["job"]["job_id"]
    _wait_published(client, job_id)

    # Identical resubmission reuses the same job (201 body carries the flag).
    again = _post(client, "dupLane", BARCODES, p1, p2)
    assert again.status_code == 201
    assert again.json()["reused"] is True
    assert again.json()["job"]["job_id"] == job_id

    # Different FASTQ content for the same lane is rejected.
    r1b = [(b"a", b"AAAAAAAATTTG", b"I" * 12)]
    p1b, p2b = helpers.make_pair(r1b, r2, tmp_path / "x" if False else tmp_path,
                                 "R1b.fastq.gz", "R2b.fastq.gz")
    conflict = _post(client, "dupLane", BARCODES, p1b, p2b)
    assert conflict.status_code == 409


def test_validation_errors(client, helpers, tmp_path):
    r1 = [(b"a", b"AAAAAAAATTTT", b"I" * 12)]
    r2 = [(b"a", b"CCCCCCCCGGGG", b"J" * 12)]
    p1, p2 = helpers.make_pair(r1, r2, tmp_path)

    # Too many samples.
    many = [
        {"sample_id": f"s{i}", "barcode_r1": "AAAAAAAA", "barcode_r2": "CCCCCCCC"}
        for i in range(9)
    ]
    resp = _post(client, "bad", many, p1, p2)
    assert resp.status_code == 400

    # Bad barcode alphabet/length.
    bad = [{"sample_id": "s1", "barcode_r1": "AAAAAAN_",
            "barcode_r2": "CCCCCCCC"}]
    resp = _post(client, "bad", bad, p1, p2)
    assert resp.status_code == 400

    # Duplicate barcode pairs.
    dup = [
        {"sample_id": "s1", "barcode_r1": "AAAAAAAA", "barcode_r2": "CCCCCCCC"},
        {"sample_id": "s2", "barcode_r1": "AAAAAAAA", "barcode_r2": "CCCCCCCC"},
    ]
    resp = _post(client, "bad", dup, p1, p2)
    assert resp.status_code == 400
