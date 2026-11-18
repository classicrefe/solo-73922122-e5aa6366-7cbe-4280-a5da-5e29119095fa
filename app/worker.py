"""Background worker: one job (demux or merge) at a time per instance.

Runs in a daemon thread inside the uvicorn process. While a slow job is being
processed the event loop keeps serving uploads, queries and downloads; the
poller simply does not claim another job until the current one finishes.
Demultiplexing jobs are claimed before merge jobs so the ingestion path is
never starved, but both queues drain on every loop.
"""
from __future__ import annotations

import json
import logging
import threading
import time
import uuid

from . import config, db, merging, processing

log = logging.getLogger("demux.worker")


def new_instance_id() -> str:
    return f"{uuid.getnode():x}-{uuid.uuid4().hex[:8]}"


class JobWorker(threading.Thread):
    def __init__(self, instance_id: str):
        super().__init__(name=f"demux-worker-{instance_id[:8]}", daemon=True)
        self.instance_id = instance_id
        self._stop = threading.Event()
        self.current_job: str | None = None
        self.current_merge_job: str | None = None

    def stop(self) -> None:
        self._stop.set()

    def run(self) -> None:
        log.info("worker %s started", self.instance_id)
        while not self._stop.is_set():
            try:
                job = db.claim_next_job(self.instance_id)
                kind = "demux"
                if job is None:
                    job = db.claim_next_merge_job(self.instance_id)
                    kind = "merge"
            except Exception:
                log.exception("failed to claim a job")
                job = None
            if job is None:
                self._stop.wait(config.POLL_INTERVAL)
                continue
            if kind == "demux":
                self._run_demux(job)
            else:
                self._run_merge(job)

    def _run_demux(self, job) -> None:
        self.current_job = job["job_id"]
        log.info("instance %s claimed job %s lane %s",
                 self.instance_id, job["job_id"], job["lane"])
        try:
            processing.process_job(job, self.instance_id)
            log.info("job %s published", job["job_id"])
        except processing.FastqFormatError as exc:
            log.warning("job %s failed on format: %s", job["job_id"], exc)
            self._fail_demux(job, f"FASTQ format error: {exc}")
        except processing.SizeExceededError as exc:
            log.warning("job %s failed on size: %s", job["job_id"], exc)
            self._fail_demux(job, str(exc))
        except Exception:
            log.exception("job %s crashed", job["job_id"])
            self._fail_demux(job, "internal processing error")
        finally:
            self.current_job = None

    def _run_merge(self, job) -> None:
        self.current_merge_job = job["merge_job_id"]
        log.info(
            "instance %s claimed merge job %s (snapshot %s, sample %s, "
            "%d lanes)", self.instance_id, job["merge_job_id"],
            job["snapshot_id"], job["sample_id"],
            len(json.loads(job["lanes_json"])))
        try:
            merging.process_merge_job(job, self.instance_id)
            log.info("merge job %s completed", job["merge_job_id"])
        except merging.MergeError as exc:
            log.warning("merge job %s failed: %s",
                        job["merge_job_id"], exc)
            self._fail_merge(job, str(exc))
        except (processing.FastqFormatError, processing.SizeExceededError) as exc:
            log.warning("merge job %s failed on lane data: %s",
                        job["merge_job_id"], exc)
            self._fail_merge(job, str(exc))
        except Exception:
            log.exception("merge job %s crashed", job["merge_job_id"])
            self._fail_merge(job, "internal merge error")
        finally:
            self.current_merge_job = None

    def _fail_demux(self, job, message: str) -> None:
        try:
            db.mark_failed(job["job_id"], self.instance_id, message)
        except Exception:
            log.exception("could not mark job %s failed", job["job_id"])

    def _fail_merge(self, job, message: str) -> None:
        try:
            db.mark_merge_failed(
                job["merge_job_id"], self.instance_id, message)
        except Exception:
            log.exception(
                "could not mark merge job %s failed", job["merge_job_id"])
