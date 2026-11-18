"""Background worker: one lane processed at a time per instance.

Runs in a daemon thread inside the uvicorn process. While a slow lane is
being processed the event loop keeps serving uploads, queries and downloads;
the poller simply does not claim another job until the current one finishes.
"""
from __future__ import annotations

import logging
import threading
import time
import uuid

from . import config, db, processing

log = logging.getLogger("demux.worker")


def new_instance_id() -> str:
    return f"{uuid.getnode():x}-{uuid.uuid4().hex[:8]}"


class JobWorker(threading.Thread):
    def __init__(self, instance_id: str):
        super().__init__(name=f"demux-worker-{instance_id[:8]}", daemon=True)
        self.instance_id = instance_id
        self._stop = threading.Event()
        self.current_job: str | None = None

    def stop(self) -> None:
        self._stop.set()

    def run(self) -> None:
        log.info("worker %s started", self.instance_id)
        while not self._stop.is_set():
            try:
                job = db.claim_next_job(self.instance_id)
            except Exception:
                log.exception("failed to claim a job")
                job = None
            if job is None:
                self._stop.wait(config.POLL_INTERVAL)
                continue
            self.current_job = job["job_id"]
            log.info("instance %s claimed job %s lane %s",
                     self.instance_id, job["job_id"], job["lane"])
            try:
                processing.process_job(job, self.instance_id)
                log.info("job %s published", job["job_id"])
            except processing.FastqFormatError as exc:
                log.warning("job %s failed on format: %s", job["job_id"], exc)
                self._fail(job, f"FASTQ format error: {exc}")
            except processing.SizeExceededError as exc:
                log.warning("job %s failed on size: %s", job["job_id"], exc)
                self._fail(job, str(exc))
            except Exception:
                log.exception("job %s crashed", job["job_id"])
                self._fail(job, "internal processing error")
            finally:
                self.current_job = None

    def _fail(self, job, message: str) -> None:
        try:
            db.mark_failed(job["job_id"], self.instance_id, message)
        except Exception:
            log.exception("could not mark job %s failed", job["job_id"])
