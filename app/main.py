"""FastAPI application: demultiplexing submission, query, download API."""
from __future__ import annotations

import hashlib
import json
import logging
import re
import shutil
import uuid
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, RedirectResponse, StreamingResponse

from . import config, db, processing
from .worker import JobWorker, new_instance_id

_IO_CHUNK = 1024 * 1024

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s %(message)s",
)
log = logging.getLogger("demux.api")

_LANE_RE = re.compile(r"^[A-Za-z0-9_.-]{1,64}$")
_SAMPLE_RE = re.compile(r"^[A-Za-z0-9_.-]{1,64}$")
_BASE_RE = re.compile(r"^[ACGT]{%d}$" % config.BARCODE_LEN)
_RESERVED_PREFIX = "__"


def _validate_barcodes(raw: str) -> list[dict[str, str]]:
    try:
        data = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise HTTPException(400, f"barcodes is not valid JSON: {exc}")
    if not isinstance(data, list) or not (1 <= len(data) <= config.MAX_SAMPLES):
        raise HTTPException(
            400, f"barcodes must be a list of 1..{config.MAX_SAMPLES} samples")
    seen_samples: set[str] = set()
    seen_combos: set[tuple[str, str]] = set()
    normalized: list[dict[str, str]] = []
    for item in data:
        if not isinstance(item, dict):
            raise HTTPException(400, "each barcode entry must be an object")
        sid = item.get("sample_id")
        b1 = item.get("barcode_r1")
        b2 = item.get("barcode_r2")
        if not isinstance(sid, str) or not _SAMPLE_RE.match(sid):
            raise HTTPException(400, f"invalid sample_id: {sid!r}")
        if sid.startswith(_RESERVED_PREFIX):
            raise HTTPException(400, f"sample_id {sid!r} uses a reserved prefix")
        if not isinstance(b1, str) or not _BASE_RE.match(b1):
            raise HTTPException(
                400, f"sample {sid}: barcode_r1 must be {config.BARCODE_LEN} A/C/G/T")
        if not isinstance(b2, str) or not _BASE_RE.match(b2):
            raise HTTPException(
                400, f"sample {sid}: barcode_r2 must be {config.BARCODE_LEN} A/C/G/T")
        if sid in seen_samples:
            raise HTTPException(400, f"duplicate sample_id: {sid}")
        combo = (b1, b2)
        if combo in seen_combos:
            raise HTTPException(400, f"duplicate barcode pair on sample {sid}")
        seen_samples.add(sid)
        seen_combos.add(combo)
        normalized.append(
            {"sample_id": sid, "barcode_r1": b1, "barcode_r2": b2})
    return normalized


async def _stream_to_disk(upload: UploadFile, dest: Path,
                          hasher: hashlib._Hash) -> int:
    """Spool an upload in fixed chunks; constant memory regardless of size."""
    total = 0
    with open(dest, "wb") as out:
        while True:
            chunk = await upload.read(_IO_CHUNK)
            if not chunk:
                break
            hasher.update(chunk)
            out.write(chunk)
            total += len(chunk)
    return total


def _job_public(row) -> dict[str, Any]:
    out = {k: row[k] for k in row.keys()}
    # Do not leak full server paths; expose only whether uploads are retained.
    out["uploads_retained"] = (
        Path(row["r1_incoming"]).is_file() and Path(row["r2_incoming"]).is_file())
    out.pop("r1_incoming", None)
    out.pop("r2_incoming", None)
    if row["status"] == "published":
        out["samples"] = [
            dict(c) for c in db.lane_counts(row["lane"])]
    else:
        out["samples"] = [dict(s) for s in db.get_samples(row["job_id"])]
    return out


def _safe_published_path(path_str: str) -> Path:
    root = config.PUBLISHED_DIR.resolve()
    p = Path(path_str).resolve()
    if root not in p.parents or not p.is_file():
        raise HTTPException(404, "file not found")
    return p


@asynccontextmanager
async def lifespan(app: FastAPI):
    db.init_db()
    notes = processing.recover_on_startup()
    for note in notes:
        log.info("startup recovery: %s", note)
    instance_id = new_instance_id()
    worker = JobWorker(instance_id)
    app.state.instance_id = instance_id
    app.state.worker = worker
    worker.start()
    log.info("instance %s ready, data dir %s", instance_id, config.DATA_DIR)
    try:
        yield
    finally:
        worker.stop()


app = FastAPI(
    title="Sequencing Lane Demultiplexing API",
    version="1.0.0",
    description=(
        "Submit a lane with up to 8 dual-index samples and paired gzip "
        "FASTQ files. Reads are assigned by combined R1+R2 8-base barcodes "
        "with at most one substitution and a unique candidate; everything "
        "else is quarantined. Barcodes and qualities are trimmed on "
        "assignment. Published lanes, per-sample files, quarantine files and "
        "counts become visible atomically."
    ),
    lifespan=lifespan,
)


@app.get("/", include_in_schema=False)
def root():
    return RedirectResponse("/docs")


@app.get("/health")
def health():
    worker: JobWorker = app.state.worker
    return {
        "status": "ok",
        "instance_id": app.state.instance_id,
        "current_job": worker.current_job,
        "queued_jobs": db.get_conn().execute(
            "SELECT COUNT(*) c FROM jobs WHERE status='queued'").fetchone()["c"],
    }


# ---------------- jobs ----------------

@app.post("/jobs", status_code=201)
async def submit_job(
    lane: str = Form(...),
    barcodes: str = Form(...),
    r1: UploadFile = File(..., description="R1 gzip FASTQ"),
    r2: UploadFile = File(..., description="R2 gzip FASTQ"),
):
    if not _LANE_RE.match(lane or ""):
        raise HTTPException(400, "lane must match [A-Za-z0-9_.-]{1,64}")
    samples = _validate_barcodes(barcodes)

    existing = db.get_job_by_lane(lane)
    upload_id = uuid.uuid4().hex
    upload_dir = config.INCOMING_DIR / upload_id
    upload_dir.mkdir(parents=True)
    hasher = hashlib.sha256()
    # Barcode table is part of the identity hash, canonicalised.
    hasher.update(json.dumps(samples, sort_keys=True).encode())
    try:
        r1_path = upload_dir / "R1.fastq.gz"
        r2_path = upload_dir / "R2.fastq.gz"
        n1 = await _stream_to_disk(r1, r1_path, hasher)
        n2 = await _stream_to_disk(r2, r2_path, hasher)
    finally:
        await r1.close()
        await r2.close()

    if n1 == 0 or n2 == 0:
        shutil.rmtree(upload_dir, ignore_errors=True)
        raise HTTPException(400, "R1/R2 uploads must be non-empty gzip files")

    # Marker lets startup recovery distinguish finished uploads from an
    # upload still streaming on another instance.
    (upload_dir / "COMPLETE").touch()

    input_hash = hasher.hexdigest()

    if existing is not None:
        shutil.rmtree(upload_dir, ignore_errors=True)
        if existing["input_hash"] != input_hash:
            raise HTTPException(
                409,
                f"lane {lane} already exists with different input; "
                "same lane requires identical barcodes and FASTQ data")
        return {
            "reused": True,
            "job": _job_public(existing),
        }

    job_id = uuid.uuid4().hex
    try:
        db.create_job(
            {
                "job_id": job_id,
                "lane": lane,
                "input_hash": input_hash,
                "barcode_table": json.dumps(samples),
                "r1_incoming": str(r1_path),
                "r2_incoming": str(r2_path),
                "created_at": db.utcnow(),
            },
            samples,
        )
    except Exception:
        # Concurrent submission for the same lane lost the UNIQUE race.
        shutil.rmtree(upload_dir, ignore_errors=True)
        rival = db.get_job_by_lane(lane)
        if rival is not None and rival["input_hash"] == input_hash:
            return {"reused": True, "job": _job_public(rival)}
        raise HTTPException(
            409, f"lane {lane} already exists with different input")
    return {"reused": False, "job": _job_public(db.get_job(job_id))}


@app.get("/jobs")
def get_jobs():
    return [_job_public(r) for r in db.list_jobs()]


@app.get("/jobs/{job_id}")
def get_job(job_id: str):
    row = db.get_job(job_id)
    if row is None:
        raise HTTPException(404, "job not found")
    return _job_public(row)


@app.post("/jobs/{job_id}/retry")
def retry_job(job_id: str):
    row = db.get_job(job_id)
    if row is None:
        raise HTTPException(404, "job not found")
    if row["status"] != "failed":
        raise HTTPException(409, f"job is {row['status']}, only failed jobs retry")
    for key in ("r1_incoming", "r2_incoming"):
        if not Path(row[key]).is_file():
            raise HTTPException(
                409, "original uploads are gone; resubmit the lane to retry")
    if not db.requeue_job(job_id):
        raise HTTPException(409, "job could not be requeued")
    return _job_public(db.get_job(job_id))


# ---------------- published lanes ----------------

@app.get("/lanes")
def list_lanes():
    return [
        {
            "lane": r["lane"],
            "job_id": r["job_id"],
            "published_at": r["published_at"],
            "seq": r["seq"],
        }
        for r in db.published_lanes()
    ]


@app.get("/lanes/{lane}")
def get_lane(lane: str):
    pub = db.get_conn().execute(
        "SELECT * FROM lane_publishes WHERE lane=?", (lane,)).fetchone()
    if pub is None:
        raise HTTPException(404, "lane not published (or unknown)")
    counts = [dict(r) for r in db.lane_counts(lane)]
    total = sum(c["read_pairs"] for c in counts)
    return {
        "lane": lane,
        "job_id": pub["job_id"],
        "published_at": pub["published_at"],
        "total_read_pairs": total,
        "counts": counts,
    }


def _bucket_file(lane: str, bucket: str, mate: str) -> Path:
    if mate not in ("r1", "r2"):
        raise HTTPException(400, "mate must be r1 or r2")
    row = db.get_conn().execute(
        "SELECT * FROM lane_counts WHERE lane=? AND sample_id=?",
        (lane, bucket)).fetchone()
    if row is None:
        raise HTTPException(404, f"bucket {bucket!r} not found in lane {lane}")
    return _safe_published_path(row[f"{mate}_file"])


@app.get("/lanes/{lane}/samples/{sample_id}/{mate}")
def download_sample(lane: str, sample_id: str, mate: str):
    path = _bucket_file(lane, sample_id, mate)
    return FileResponse(
        path, media_type="application/gzip",
        filename=f"{lane}.{sample_id}.{mate}.fastq.gz")


@app.get("/lanes/{lane}/quarantine/{mate}")
def download_quarantine(lane: str, mate: str):
    path = _bucket_file(lane, processing.QUARANTINE, mate)
    return FileResponse(
        path, media_type="application/gzip",
        filename=f"{lane}.quarantine.{mate}.fastq.gz")


# ---------------- snapshots ----------------

@app.post("/snapshots", status_code=201)
def create_snapshot():
    lanes = [r["lane"] for r in db.published_lanes()]
    return db.create_snapshot(uuid.uuid4().hex, lanes)


@app.get("/snapshots")
def get_snapshots():
    return db.list_snapshots()


@app.get("/snapshots/{snapshot_id}")
def get_snapshot(snapshot_id: str):
    snap = db.get_snapshot(snapshot_id)
    if snap is None:
        raise HTTPException(404, "snapshot not found")
    return snap


def _snapshot_bucket_files(snapshot_id: str, sample_id: str, mate: str):
    if mate not in ("r1", "r2"):
        raise HTTPException(400, "mate must be r1 or r2")
    snap = db.get_snapshot(snapshot_id)
    if snap is None:
        raise HTTPException(404, "snapshot not found")
    files: list[Path] = []
    conn = db.get_conn()
    for lane in snap["lanes"]:
        row = conn.execute(
            "SELECT * FROM lane_counts WHERE lane=? AND sample_id=?",
            (lane, sample_id)).fetchone()
        if row is None:
            continue  # sample absent from this lane's table
        files.append(_safe_published_path(row[f"{mate}_file"]))
    return snap, files


def _iter_concat(files: list[Path]):
    # gzip member concatenation is a valid gzip stream per RFC 1952.
    for path in files:
        with open(path, "rb") as fh:
            while True:
                chunk = fh.read(_IO_CHUNK)
                if not chunk:
                    break
                yield chunk


@app.get("/snapshots/{snapshot_id}/samples/{sample_id}/{mate}")
def snapshot_sample(snapshot_id: str, sample_id: str, mate: str):
    snap, files = _snapshot_bucket_files(snapshot_id, sample_id, mate)
    if not files:
        raise HTTPException(404, f"sample {sample_id!r} absent from snapshot")
    headers = {"X-Snapshot-Lanes": ",".join(snap["lanes"])}
    return StreamingResponse(
        _iter_concat(files), media_type="application/gzip",
        headers={**headers,
                 "content-disposition":
                     f'attachment; filename="{sample_id}.{mate}.fastq.gz"'})


@app.get("/snapshots/{snapshot_id}/quarantine/{mate}")
def snapshot_quarantine(snapshot_id: str, mate: str):
    snap, files = _snapshot_bucket_files(
        snapshot_id, processing.QUARANTINE, mate)
    if not files:
        raise HTTPException(404, "quarantine bucket absent from snapshot")
    headers = {"X-Snapshot-Lanes": ",".join(snap["lanes"])}
    return StreamingResponse(
        _iter_concat(files), media_type="application/gzip",
        headers={**headers,
                 "content-disposition":
                     f'attachment; filename="quarantine.{mate}.fastq.gz"'})
