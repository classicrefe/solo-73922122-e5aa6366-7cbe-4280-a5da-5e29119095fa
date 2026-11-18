"""End-to-end tests for the snapshot overlap-merge API."""
import gzip
import io
import json
import random
import time

BARCODES = [
    {"sample_id": "s1", "barcode_r1": "AAAAAAAA", "barcode_r2": "CCCCCCCC"},
    {"sample_id": "s2", "barcode_r1": "GGGGGGGG", "barcode_r2": "TTTTTTTT"},
]


def _post_job(client, lane, r1_path, r2_path, barcodes=BARCODES):
    with open(r1_path, "rb") as f1, open(r2_path, "rb") as f2:
        return client.post(
            "/jobs",
            data={"lane": lane, "barcodes": json.dumps(barcodes)},
            files={"r1": ("R1.fastq.gz", f1, "application/gzip"),
                   "r2": ("R2.fastq.gz", f2, "application/gzip")})


def _wait(client, job_id, path="/jobs/{}", timeout=20):
    deadline = time.time() + timeout
    while time.time() < deadline:
        row = client.get(path.format(job_id)).json()
        if row["status"] in ("published", "failed", "completed"):
            return row
        time.sleep(0.05)
    raise AssertionError("timeout waiting for job")


def _wait_merge(client, merge_id, timeout=20):
    deadline = time.time() + timeout
    while time.time() < deadline:
        row = client.get(f"/merge/jobs/{merge_id}").json()
        if row["status"] in ("completed", "failed"):
            return row
        time.sleep(0.05)
    raise AssertionError("timeout waiting for merge job")


def _overlap_pair(rng, name, overlap=30, total=50, mismatch=0,
                  qual1=b"I", qual2=b"J"):
    """Build one overlapping R1/R2 pair (R2 given in raw orientation).

    R2' (reverse-complemented R2) shares exactly ``overlap`` bases with R1's
    tail; its remaining bases are T so no shorter/longer alignment qualifies,
    making the overlap unique.
    """
    shared = bytes(rng.choice(b"ACGT") for _ in range(overlap))
    if mismatch:
        shared = bytearray(shared)
        for pos in rng.sample(range(overlap), mismatch):
            shared[pos] = rng.choice(
                [b for b in b"ACGT" if b != shared[pos]])
        shared = bytes(shared)
    r1_extra = bytes(rng.choice(b"ACGT")
                     for _ in range(total - overlap))
    r1_seq = r1_extra + shared
    r2_rc = shared + b"T" * (total - overlap)
    r2_seq = r2_rc.translate(bytes.maketrans(b"ACGT", b"TGCA"))[::-1]
    return ((name, r1_seq, qual1 * total),
            (name, r2_seq, qual2 * total))


def _nonoverlap_pair(rng, name, total=50):
    r1 = bytes(rng.choice(b"AC") for _ in range(total))       # A/C only
    r2 = bytes(rng.choice(b"GT") for _ in range(total))       # G/T only
    return (name, r1, b"I" * total), (name, r2, b"J" * total)


def _write_pairs(tmp, pairs1, pairs2, sub="lane"):
    d = tmp / sub
    d.mkdir(parents=True, exist_ok=True)
    p1, p2 = d / "R1.fastq.gz", d / "R2.fastq.gz"
    with gzip.open(p1, "wb") as g:
        for name, seq, qual in pairs1:
            g.write(b"@" + name + b"\n" + seq + b"\n+\n" + qual + b"\n")
    with gzip.open(p2, "wb") as g:
        for name, seq, qual in pairs2:
            g.write(b"@" + name + b"\n" + seq + b"\n+\n" + qual + b"\n")
    return p1, p2


def _publish_lane(client, tmp, lane, pairs, sub):
    # Demux needs the dual-index barcodes at the read heads; after trimming,
    # the published sample reads keep exactly the engineered insert.
    b1, b2 = b"AAAAAAAA", b"CCCCCCCC"
    r1_recs, r2_recs = [], []
    for (name, s1, q1), (_, s2, q2) in pairs:
        r1_recs.append((name, b1 + s1, b"I" * len(b1) + q1))
        r2_recs.append((name, b2 + s2, b"J" * len(b2) + q2))
    p1, p2 = _write_pairs(tmp, r1_recs, r2_recs, sub)
    resp = _post_job(client, lane, p1, p2)
    assert resp.status_code == 201, resp.text
    job = _wait(client, resp.json()["job"]["job_id"])
    assert job["status"] == "published", job
    return lane


def _build_merging_lane(tmp, lane, n_merge=4, n_unmerged=2, seed=1, sub="l1"):
    rng = random.Random(seed)
    pairs = []
    for i in range(n_merge):
        pairs.append(_overlap_pair(rng, f"m{i}".encode(), overlap=30))
    for i in range(n_unmerged):
        pairs.append(_nonoverlap_pair(rng, f"u{i}".encode()))
    return pairs


def test_merge_happy_path_counts_and_output(client, tmp_path):
    pairs = _build_merging_lane(tmp_path, "L1")
    _publish_lane(client, tmp_path, "L1", pairs, "l1")
    snap = client.post("/snapshots").json()

    resp = client.post("/merge/jobs",
                       json={"snapshot_id": snap["snapshot_id"],
                             "sample_id": "s1"})
    assert resp.status_code == 201, resp.text
    assert resp.json()["reused"] is False
    merge_id = resp.json()["job"]["merge_job_id"]
    done = _wait_merge(client, merge_id)
    assert done["status"] == "completed", done
    assert done["merged_count"] == 4
    assert done["unmerged_count"] == 2
    assert done["total_pairs"] == 6
    assert done["lanes"] == ["L1"]
    assert done["download_available"] is True

    dl = client.get(f"/merge/jobs/{merge_id}/download")
    assert dl.status_code == 200
    assert dl.headers["content-type"] == "application/gzip"
    with gzip.open(io.BytesIO(dl.content)) as g:
        lines = g.read().splitlines()
    # Four merged single-end records => 16 FASTQ lines.
    assert len(lines) == 4 * 4
    # Every merged sequence/quality pair has equal length and every name came
    # from an R1 name of the mergeable subset.
    names = {lines[i][1:] for i in range(0, len(lines), 4)}
    assert names == {b"m0", b"m1", b"m2", b"m3"}
    for i in range(0, len(lines), 4):
        assert len(lines[i + 1]) == len(lines[i + 3])
        assert lines[i + 2] == b"+"


def test_merge_read_names_match_demux_output(client, tmp_path):
    pairs = _build_merging_lane(tmp_path, "L2", n_merge=3, n_unmerged=0,
                                seed=11, sub="l2")
    _publish_lane(client, tmp_path, "L2", pairs, "l2")
    # Reference: the demux R1 output names (barcode trimming keeps names).
    ref = client.get("/lanes/L2/samples/s1/r1")
    with gzip.open(io.BytesIO(ref.content)) as g:
        ref_names = {line[1:] for line in g.read().splitlines()
                     if line.startswith(b"@")}
    snap = client.post("/snapshots").json()
    mid = client.post("/merge/jobs",
                      json={"snapshot_id": snap["snapshot_id"],
                            "sample_id": "s1"}).json()["job"]["merge_job_id"]
    _wait_merge(client, mid)
    with gzip.open(io.BytesIO(
            client.get(f"/merge/jobs/{mid}/download").content)) as g:
        merged_names = {line[1:] for line in g.read().splitlines()
                        if line.startswith(b"@")}
    assert merged_names == ref_names


def test_merge_concatenates_all_snapshot_lanes(client, tmp_path):
    p_a = _build_merging_lane(tmp_path, "A", n_merge=3, n_unmerged=1,
                              seed=2, sub="a")
    p_b = _build_merging_lane(tmp_path, "B", n_merge=2, n_unmerged=2,
                              seed=3, sub="b")
    _publish_lane(client, tmp_path, "A", p_a, "a")
    _publish_lane(client, tmp_path, "B", p_b, "b")
    snap = client.post("/snapshots").json()
    assert snap["lanes"] == ["A", "B"]
    mid = client.post("/merge/jobs",
                      json={"snapshot_id": snap["snapshot_id"],
                            "sample_id": "s1"}).json()["job"]["merge_job_id"]
    done = _wait_merge(client, mid)
    assert done["status"] == "completed", done
    assert done["total_pairs"] == 8       # 4 per lane
    assert done["merged_count"] == 5
    assert done["unmerged_count"] == 3
    assert done["lanes"] == ["A", "B"]


def test_snapshot_isolation_lanes_published_later_excluded(client, tmp_path):
    pairs = _build_merging_lane(tmp_path, "early", n_merge=2, n_unmerged=0,
                                seed=4, sub="e")
    _publish_lane(client, tmp_path, "early", pairs, "e")
    snap = client.post("/snapshots").json()

    later = _build_merging_lane(tmp_path, "later", n_merge=5, n_unmerged=0,
                                seed=5, sub="l")
    _publish_lane(client, tmp_path, "later", later, "l")

    mid = client.post("/merge/jobs",
                      json={"snapshot_id": snap["snapshot_id"],
                            "sample_id": "s1"}).json()["job"]["merge_job_id"]
    done = _wait_merge(client, mid)
    assert done["total_pairs"] == 2       # only the frozen lane
    assert done["lanes"] == ["early"]


def test_merge_unknown_snapshot_and_missing_sample_are_404(client, tmp_path):
    pairs = _build_merging_lane(tmp_path, "L3", sub="l3")
    _publish_lane(client, tmp_path, "L3", pairs, "l3")
    snap = client.post("/snapshots").json()

    r = client.post("/merge/jobs",
                    json={"snapshot_id": "nope", "sample_id": "s1"})
    assert r.status_code == 404
    r = client.post("/merge/jobs",
                    json={"snapshot_id": snap["snapshot_id"],
                          "sample_id": "ghost"})
    assert r.status_code == 404
    r = client.post("/merge/jobs",
                    json={"snapshot_id": snap["snapshot_id"], "sample_id": ""})
    assert r.status_code in (400, 422)  # explicit regex 400 or pydantic 422


def test_duplicate_merge_request_is_idempotent(client, tmp_path):
    pairs = _build_merging_lane(tmp_path, "L4", sub="l4")
    _publish_lane(client, tmp_path, "L4", pairs, "l4")
    snap = client.post("/snapshots").json()
    body = {"snapshot_id": snap["snapshot_id"], "sample_id": "s1"}
    first = client.post("/merge/jobs", json=body).json()
    second = client.post("/merge/jobs", json=body).json()
    assert second["reused"] is True
    assert second["job"]["merge_job_id"] == first["job"]["merge_job_id"]
    listed = client.get("/merge/jobs").json()
    assert len(listed) == 1


def test_corrupt_lane_fails_whole_merge_without_output(client, tmp_path):
    pairs = _build_merging_lane(tmp_path, "L5", sub="l5")
    _publish_lane(client, tmp_path, "L5", pairs, "l5")
    snap = client.post("/snapshots").json()
    mid = client.post(
        "/merge/jobs",
        json={"snapshot_id": snap["snapshot_id"], "sample_id": "s1"}
    ).json()["job"]["merge_job_id"]

    # Corrupt the published R1 gzip of the snapshot lane in place.
    from app import db as dbmod
    row = dbmod.get_conn().execute(
        "SELECT r1_file FROM lane_counts WHERE lane='L5' AND sample_id='s1'"
    ).fetchone()
    p = row["r1_file"]
    with open(p, "rb") as fh:
        blob = fh.read()
    with open(p, "wb") as fh:
        fh.write(blob[:-8])  # truncated gzip stream

    done = _wait_merge(client, mid)
    assert done["status"] == "failed"
    assert "gzip" in done["error"].lower() or "format" in done["error"].lower()
    assert client.get(f"/merge/jobs/{mid}/download").status_code == 409
    # No half file, staging or temporary publish artifacts remain.
    from app import config as cfg
    assert not (cfg.MERGES_DIR / mid).exists()
    assert not (cfg.MERGES_DIR / f".tmp.{mid}").exists()
    assert not (cfg.MERGE_WORK_DIR / mid / "output").exists()


def test_overlong_read_fails_merge(client, tmp_path):
    rng = random.Random(6)
    b1, b2 = b"AAAAAAAA", b"CCCCCCCC"

    def _post_with_barcodes(lane, pairs, sub):
        r1 = [(n, b1 + s, b"I" * len(b1) + q) for (n, s, q), _ in pairs]
        r2 = [(n, b2 + s, b"J" * len(b2) + q) for _, (n, s, q) in pairs]
        p1, p2 = _write_pairs(tmp_path, r1, r2, sub)
        _post_job(client, lane, p1, p2)
        deadline = time.time() + 20
        while time.time() < deadline:
            if client.get(f"/lanes/{lane}").status_code == 200:
                return
            time.sleep(0.05)
        raise AssertionError(f"lane {lane} never published")

    good = _overlap_pair(rng, b"ok", overlap=30, total=50)
    long_r1, long_r2 = _overlap_pair(rng, b"long", overlap=30, total=151)
    _post_with_barcodes("Llong", [good, (long_r1, long_r2)], "long")
    snap = client.post("/snapshots").json()
    mid = client.post(
        "/merge/jobs",
        json={"snapshot_id": snap["snapshot_id"], "sample_id": "s1"}
    ).json()["job"]["merge_job_id"]
    done = _wait_merge(client, mid)
    assert done["status"] == "failed"
    assert "151" in done["error"] and "outside" in done["error"]


def test_short_read_fails_merge(client, tmp_path):
    rng = random.Random(7)
    b1, b2 = b"AAAAAAAA", b"CCCCCCCC"
    short_r1, short_r2 = _overlap_pair(rng, b"short", overlap=20, total=39)
    r1 = [(short_r1[0], b1 + short_r1[1], b"I" * 8 + short_r1[2])]
    r2 = [(short_r2[0], b2 + short_r2[1], b"J" * 8 + short_r2[2])]
    p1, p2 = _write_pairs(tmp_path, r1, r2, "short")
    _post_job(client, "Lshort", p1, p2)
    deadline = time.time() + 20
    while time.time() < deadline:
        if client.get("/lanes/Lshort").status_code == 200:
            break
        time.sleep(0.05)
    snap = client.post("/snapshots").json()
    mid = client.post(
        "/merge/jobs",
        json={"snapshot_id": snap["snapshot_id"], "sample_id": "s1"}
    ).json()["job"]["merge_job_id"]
    done = _wait_merge(client, mid)
    assert done["status"] == "failed"
    assert "39" in done["error"] and "outside" in done["error"]


def test_failed_merge_retries_and_completes(client, tmp_path):
    pairs = _build_merging_lane(tmp_path, "L6", sub="l6")
    _publish_lane(client, tmp_path, "L6", pairs, "l6")
    snap = client.post("/snapshots").json()
    mid = client.post(
        "/merge/jobs",
        json={"snapshot_id": snap["snapshot_id"], "sample_id": "s1"}
    ).json()["job"]["merge_job_id"]

    from app import db as dbmod
    row = dbmod.get_conn().execute(
        "SELECT r1_file FROM lane_counts WHERE lane='L6' AND sample_id='s1'"
    ).fetchone()
    target = row["r1_file"]
    original = open(target, "rb").read()
    with open(target, "wb") as fh:
        fh.write(original[:-8])
    done = _wait_merge(client, mid)
    assert done["status"] == "failed"
    # Retry is refused while inputs stay broken...
    assert client.post(f"/merge/jobs/{mid}/retry").status_code == 200
    assert _wait_merge(client, mid)["status"] == "failed"
    # ...then succeeds once the published file is repaired.
    with open(target, "wb") as fh:
        fh.write(original)
    assert client.post(f"/merge/jobs/{mid}/retry").status_code == 200
    fixed = _wait_merge(client, mid)
    assert fixed["status"] == "completed"
    assert fixed["merged_count"] == 4


def test_published_count_tampering_fails_merge(client, tmp_path):
    pairs = _build_merging_lane(tmp_path, "L7", n_merge=2, n_unmerged=0,
                                seed=7, sub="l7")
    _publish_lane(client, tmp_path, "L7", pairs, "l7")
    snap = client.post("/snapshots").json()
    from app import db as dbmod
    # Lie about how many pairs the published files contain.
    dbmod.get_conn().execute(
        "UPDATE lane_counts SET read_pairs=read_pairs+1 "
        "WHERE lane='L7' AND sample_id='s1'")
    mid = client.post(
        "/merge/jobs",
        json={"snapshot_id": snap["snapshot_id"], "sample_id": "s1"}
    ).json()["job"]["merge_job_id"]
    done = _wait_merge(client, mid)
    assert done["status"] == "failed"
    assert "count" in done["error"].lower()
    assert client.get(f"/merge/jobs/{mid}/download").status_code == 409


def test_lanes_without_sample_are_skipped(client, tmp_path):
    # Lane A has s1; lane B is published with a different sample set so s1 is
    # absent there. The merge must skip B and process only A's pairs.
    p_a = _build_merging_lane(tmp_path, "mA", n_merge=3, n_unmerged=1,
                              seed=31, sub="ma")
    _publish_lane(client, tmp_path, "mA", p_a, "ma")
    only_s2 = [
        {"sample_id": "s2", "barcode_r1": "GGGGGGGG",
         "barcode_r2": "TTTTTTTT"},
    ]
    rng = random.Random(32)
    pairs_b = []
    for i in range(2):
        r1, r2 = _overlap_pair(rng, f"x{i}".encode(), overlap=30)
        pairs_b.append((r1, r2))
    r1_recs = [(n, b"GGGGGGGG" + s, b"I" * 8 + q) for (n, s, q), _ in pairs_b]
    r2_recs = [(n, b"TTTTTTTT" + s, b"J" * 8 + q) for _, (n, s, q) in pairs_b]
    p1, p2 = _write_pairs(tmp_path, r1_recs, r2_recs, "mb")
    resp = _post_job(client, "mB", p1, p2, barcodes=only_s2)
    assert resp.status_code == 201
    _wait(client, resp.json()["job"]["job_id"])

    snap = client.post("/snapshots").json()
    assert snap["lanes"] == ["mA", "mB"]
    mid = client.post(
        "/merge/jobs",
        json={"snapshot_id": snap["snapshot_id"], "sample_id": "s1"}
    ).json()["job"]["merge_job_id"]
    done = _wait_merge(client, mid)
    assert done["status"] == "completed", done
    assert done["total_pairs"] == 4       # only lane mA's s1 bucket
    assert done["lanes"] == ["mA"]        # mB recorded as skipped
    assert done["merged_count"] == 3


def test_other_requests_stay_available_during_merge(client, tmp_path):
    # A large lane keeps the worker busy; health and queries must stay live.
    rng = random.Random(8)
    pairs = [_overlap_pair(rng, f"b{i}".encode(), overlap=30)
             for i in range(20000)]
    _publish_lane(client, tmp_path, "big", pairs, "big")
    snap = client.post("/snapshots").json()
    mid = client.post(
        "/merge/jobs",
        json={"snapshot_id": snap["snapshot_id"], "sample_id": "s1"}
    ).json()["job"]["merge_job_id"]
    # The event loop answers immediately while the background thread merges.
    assert client.get("/health").status_code == 200
    assert client.get("/lanes").status_code == 200
    assert client.get("/lanes/big").status_code == 200
    dl = client.get("/lanes/big/samples/s1/r1")
    assert dl.status_code == 200 and len(dl.content) > 0
    done = _wait_merge(client, mid, timeout=60)
    assert done["status"] == "completed"
    assert done["merged_count"] == 20000
