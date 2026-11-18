#!/usr/bin/env python3
"""Two-instance integration test using real uvicorn processes.

Simulates the docker-compose deployment:
  * two uvicorn workers on different ports, one shared data directory;
  * lanes submitted while one instance is busy still get accepted (queued)
    and query/download stay responsive;
  * jobs are claimed by exactly one instance each and work is shared;
  * duplicate submissions reuse the task, conflicting input is rejected;
  * snapshots fix the published lane set.

Hard-crash recovery (SIGKILL mid-lane) is covered by the pytest suite
(tests/test_failure_and_recovery.py) because some sandboxes forbid kill().
When signals are allowed here, the optional crash section runs as well.

Run:  python scripts/integration_two_instances.py
"""
from __future__ import annotations

import json
import os
import shutil
import socket
import subprocess
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scripts.make_demo_fastq import main as make_demo  # noqa: E402


def free_port() -> int:
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    return port


PORT_A, PORT_B = free_port(), free_port()
BASE_A = f"http://127.0.0.1:{PORT_A}"
BASE_B = f"http://127.0.0.1:{PORT_B}"


def http(method: str, url: str, data: bytes | None = None,
         ctype: str = "application/json", raw_path: tuple | None = None):
    if raw_path is not None:
        boundary, body = raw_path  # multipart built manually, no extra deps
        req = urllib.request.Request(
            url, data=body, method=method,
            headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    else:
        req = urllib.request.Request(
            url, data=data, method=method,
            headers={"Content-Type": ctype} if data else {})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            return resp.status, resp.read()
    except urllib.error.HTTPError as exc:
        return exc.code, exc.read()


def wait_health(port: int, timeout: float = 20) -> None:
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            code, _ = http("GET", f"http://127.0.0.1:{port}/health")
            if code == 200:
                return
        except OSError:
            pass
        time.sleep(0.3)
    raise RuntimeError(f"port {port} never became healthy")


def start_instance(port: int, data_dir: Path) -> subprocess.Popen:
    env = {
        **os.environ,
        "DEMUX_DATA_DIR": str(data_dir),
        "DEMUX_POLL_INTERVAL": "0.2",
    }
    log = open(data_dir / f"uvicorn-{port}.log", "ab")
    proc = subprocess.Popen(
        [sys.executable, "-m", "uvicorn", "app.main:app",
         "--host", "127.0.0.1", "--port", str(port)],
        cwd=str(ROOT), env=env, stdout=log, stderr=subprocess.STDOUT)
    try:
        wait_health(port)
    except Exception:
        log.close()
        raise
    return proc


def stop(proc: subprocess.Popen) -> bool:
    """Return True if the process was signalled; False if the sandbox denies."""
    try:
        proc.terminate()
        proc.wait(timeout=15)
        return True
    except PermissionError:
        return False


def multipart(lane: str, barcodes: list, r1: Path, r2: Path) -> tuple:
    boundary = "----demuxboundary"
    parts: list[bytes] = []

    def field(name: str, value: str) -> None:
        parts.append(f"--{boundary}\r\nContent-Disposition: form-data; "
                     f'name="{name}"\r\n\r\n{value}\r\n'.encode())

    def file_field(name: str, path: Path) -> None:
        parts.append(
            f"--{boundary}\r\nContent-Disposition: form-data; name=\"{name}\"; "
            f'filename="{path.name}"\r\nContent-Type: application/gzip\r\n\r\n'
            .encode())
        parts.append(path.read_bytes())
        parts.append(b"\r\n")

    field("lane", lane)
    field("barcodes", json.dumps(barcodes))
    file_field("r1", r1)
    file_field("r2", r2)
    parts.append(f"--{boundary}--\r\n".encode())
    return boundary, b"".join(parts)


def wait_job(base: str, job_id: str, timeout: float = 90) -> dict:
    deadline = time.time() + timeout
    while time.time() < deadline:
        code, body = http("GET", f"{base}/jobs/{job_id}")
        row = json.loads(body)
        if row["status"] in ("published", "failed"):
            return row
        time.sleep(0.2)
    raise RuntimeError("job timeout")


def main() -> int:
    data_dir = Path(os.environ.get("TMPDIR", "/tmp")) / f"demux-it-{os.getpid()}"
    shutil.rmtree(data_dir, ignore_errors=True)
    data_dir.mkdir(parents=True)
    fixture_dir = data_dir / "fixture"
    make_demo(str(fixture_dir), pairs=20_000)
    r1, r2 = fixture_dir / "R1.fastq.gz", fixture_dir / "R2.fastq.gz"
    barcodes = [
        {"sample_id": "sampleA", "barcode_r1": "AAAAAAAA",
         "barcode_r2": "CCCCCCCC"},
        {"sample_id": "sampleB", "barcode_r1": "GGGGGGGG",
         "barcode_r2": "TTTTTTTT"},
    ]
    procs: list[subprocess.Popen] = []

    print("[1] start two instances sharing one volume")
    a = start_instance(PORT_A, data_dir)
    b = start_instance(PORT_B, data_dir)
    procs = [a, b]

    try:
        print("[2] submit 8 lanes; work must be split across both instances")
        job_ids = []
        for i in range(8):
            boundary, body = multipart(f"laneIT{i}", barcodes, r1, r2)
            code, resp = http("POST", f"{BASE_A}/jobs",
                              raw_path=(boundary, body))
            assert code == 201, (code, resp)
            job_ids.append(json.loads(resp)["job"]["job_id"])
        owners = {}
        for jid in job_ids:
            row = wait_job(BASE_A, jid)
            assert row["status"] == "published", row
            owners[jid] = row["instance_id"]
        distinct = set(owners.values())
        print(f"    {len(distinct)} distinct owners across 8 lanes: {distinct}")
        assert len(distinct) == 2, "both instances must claim jobs"

        print("[3] new submission/query/download stay responsive during load")
        code, _ = http("GET", f"{BASE_A}/health")
        assert code == 200
        boundary, body = multipart("laneITExtra", barcodes, r1, r2)
        code, resp = http("POST", f"{BASE_B}/jobs", raw_path=(boundary, body))
        assert code == 201, resp
        jid = json.loads(resp)["job"]["job_id"]
        assert wait_job(BASE_B, jid)["status"] == "published"

        for base in (BASE_A, BASE_B):
            for path in ("/lanes/laneIT0/samples/sampleA/r1",
                         "/lanes/laneIT0/samples/sampleB/r2",
                         "/lanes/laneIT0/quarantine/r1"):
                code, body = http("GET", base + path)
                assert code == 200 and len(body) > 0, (base, path, code)
        code, lane = http("GET", f"{BASE_A}/lanes/laneIT0")
        counts = {c["sample_id"]: c["read_pairs"]
                  for c in json.loads(lane)["counts"]}
        print("    laneIT0 counts:", counts)
        assert counts["__quarantine__"] == 2

        print("[4] snapshot fixes the lane list (later lane excluded)")
        code, snap = http("POST", f"{BASE_A}/snapshots")
        snap_id = json.loads(snap)["snapshot_id"]
        before = set(json.loads(snap)["lanes"])
        boundary, body = multipart("laneITLater", barcodes, r1, r2)
        code, resp = http("POST", f"{BASE_A}/jobs", raw_path=(boundary, body))
        assert wait_job(BASE_A, json.loads(resp)["job"]["job_id"])[
            "status"] == "published"
        code, snap2 = http("GET", f"{BASE_B}/snapshots/{snap_id}")
        assert set(json.loads(snap2)["lanes"]) == before
        print(f"    snapshot keeps {len(before)} lanes; later lane excluded")

        print("[5] duplicate lane reuses task; different input rejected (409)")
        boundary, body = multipart("laneIT0", barcodes, r1, r2)
        code, resp = http("POST", f"{BASE_B}/jobs", raw_path=(boundary, body))
        assert code == 201 and json.loads(resp)["reused"] is True
        bad_barcodes = [{**barcodes[0], "barcode_r1": "AAAAAAAC"}]
        boundary, body = multipart("laneIT0", bad_barcodes, r1, r2)
        code, resp = http("POST", f"{BASE_B}/jobs", raw_path=(boundary, body))
        assert code == 409, (code, resp)
    finally:
        signalled = [stop(p) for p in procs]

    if not all(signalled):
        print("\n[6] SKIPPED hard-crash section: sandbox denies process signals")
        print("    Crash recovery is covered by: pytest tests/test_failure_and_recovery.py")
        print("\nINTEGRATION CHECKS PASSED (crash section skipped here)")
        return 0

    print("[6] hard-crash recovery (SIGKILL mid-lane, restart cleanup)")
    c = start_instance(PORT_A, data_dir)
    try:
        boundary, body = multipart("laneCrash", barcodes, r1, r2)
        code, resp = http("POST", f"{BASE_A}/jobs", raw_path=(boundary, body))
        crash_jid = json.loads(resp)["job"]["job_id"]
        c.kill()
        c.wait(timeout=15)
        os.environ["DEMUX_DATA_DIR"] = str(data_dir)
        import importlib
        from app import config as cfg, db as dbmod
        importlib.reload(cfg)
        importlib.reload(dbmod)
        dbmod.init_db()
        dbmod.get_conn().execute(
            "UPDATE jobs SET status='processing', instance_id='ghost', "
            "heartbeat_at=? WHERE job_id=?",
            ("2000-01-01T00:00:00+00:00", crash_jid))
        wd = cfg.WORK_DIR / crash_jid
        wd.mkdir(parents=True, exist_ok=True)
        (wd / "lock").write_text("ghost\n")
        (wd / "output").mkdir(exist_ok=True)

        d = start_instance(PORT_B, data_dir)
        try:
            deadline = time.time() + 15
            body = b""
            while time.time() < deadline:
                code, body = http("GET", f"{BASE_B}/jobs/{crash_jid}")
                if json.loads(body)["status"] == "failed":
                    break
                time.sleep(0.3)
            row = json.loads(body)
            assert row["status"] == "failed" and "interrupted" in row["error"]
            code, _ = http("GET", f"{BASE_B}/lanes/laneCrash")
            assert code == 404  # no half-published lane
            code, _ = http("POST", f"{BASE_B}/jobs/{crash_jid}/retry")
            assert code == 200
            assert wait_job(BASE_B, crash_jid)["status"] == "published"
        finally:
            stop(d)
    finally:
        if c.poll() is None:
            try:
                c.kill()
            except PermissionError:
                pass

    print("\nALL INTEGRATION CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
