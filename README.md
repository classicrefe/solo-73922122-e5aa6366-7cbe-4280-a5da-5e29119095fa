# Sequencing Lane Demultiplexing & Read-Merging API

FastAPI + SQLite service that splits mixed paired-end gzip FASTQ uploads into
per-sample files and merges a snapshot sample's overlapping paired-end reads
into single-end reads. Designed to run as **two container instances on one
shared volume** with a **128 MiB memory limit**.

## What it does

* Accepts a lane number, a dual-index barcode table (1–8 samples, two
  8-base barcodes each) and paired `R1`/`R2` gzip FASTQ files.
  Decompressed input across both mates is capped at **512 MiB**.
* Validates both mates: valid gzip, well-formed FASTQ, equal record counts,
  identical read names in the same order, matching sequence/quality lengths.
  Any corruption fails the **whole lane**.
* Assigns a pair to a sample when the combined R1+R2 barcode matches with
  **at most one substitution** and the candidate sample is **unique**.
  Everything else (no hit, >1 substitution, or a tie) goes to quarantine.
* Assigned reads have both 8-base barcodes and their corresponding quality
  characters trimmed; quarantined reads are kept intact.
* A lane becomes visible only when **all** per-sample paired files, the
  quarantine paired files and the count table exist (atomic publish:
  files staged, renamed, then a single DB transaction).
* Same lane + identical input (hash of barcode table and both files) reuses
  the existing task; same lane with different input is rejected with `409`.
* Failed tasks can be retried (`POST /jobs/{id}/retry`); uploads are kept
  until a lane is published.
* Snapshots freeze the published lane list at creation time. Sample/
  quarantine downloads through a snapshot only ever concatenate files from
  that fixed list — lanes published later never mix in.

## Overlap read merging (`POST /merge/jobs`)

For one sample within a snapshot, the service reads that sample's R1/R2 from
**exactly the snapshot's frozen lane list** (lanes published during or after
the run never mix in) and merges each pair into a single-end read:

* Only mates of **40..150 bases** are accepted; any mate outside that window,
  any corrupt gzip/FASTQ, mismatched read names or unequal record counts fail
  the **whole merge job**.
* R2 is reverse-complemented, then the merger looks for the **unique ungapped
  overlap** between the R1 tail and the R2' head of **at least 12 bases** with
  a **mismatch rate of at most 5%** (inclusive). No qualifying position, or
  more than one, means the pair is left **unmerged** (counted, not written).
* Overlap quality rules — same base: Phred scores are **added and capped at
  40**; different bases: the **higher-quality base wins**, and on an equal
  quality the base becomes **`N`** carrying the lower (equal) quality.
* Output is one gzip FASTQ containing only the merged single-end reads, plus
  the **`merged_count` and `unmerged_count`** totals. Read names, FASTQ
  framing and the published per-lane pair counts are cross-checked per pair
  and before publish; a count mismatch fails the job.
* A failed merge never exposes a half file: staging and temporary output are
  removed, no download row is committed, and the job can be retried
  (`POST /merge/jobs/{id}/retry`).
* Merging streams record by record (verified flat RSS on hundreds of MiB),
  runs in the same single background worker as demultiplexing, and never
  blocks uploads, queries or downloads on the event loop.

### Merge demo

```bash
python scripts/make_merge_demo.py /tmp/merge-demo     # engineered overlaps
# Submit the two generated gzip files through POST /jobs like any lane, then:
curl -X POST http://localhost:8001/snapshots          # freeze the lane list
curl -X POST http://localhost:8001/merge/jobs \
  -H 'content-type: application/json' \
  -d '{"snapshot_id":"<id>","sample_id":"sampleA"}'
curl http://localhost:8001/merge/jobs/<merge_job_id>            # counts
curl -OJ http://localhost:8001/merge/jobs/<merge_job_id>/download
```

## Running (Docker Compose, two instances on a shared volume)

```bash
docker compose up --build
# instance A: http://localhost:8001   instance B: http://localhost:8002
# OpenAPI demo UI:  http://localhost:8001/docs  (also on 8002)
```

Both containers mount the named volume `demux-shared` at `/data`, are capped
at 128 MiB (`mem_limit`), and automatically recover on restart.

### Quick demo

```bash
python scripts/make_demo_fastq.py /tmp/demo          # makes R1/R2.fastq.gz
# In http://localhost:8001/docs -> POST /jobs:
#   lane: demoLane
#   barcodes: [{"sample_id":"sampleA","barcode_r1":"AAAAAAAA","barcode_r2":"CCCCCCCC"},
#              {"sample_id":"sampleB","barcode_r1":"GGGGGGGG","barcode_r2":"TTTTTTTT"}]
#   r1/r2: the generated gzip files
curl http://localhost:8001/lanes/demoLane             # counts
curl -OJ http://localhost:8001/lanes/demoLane/samples/sampleA/r1
curl -OJ http://localhost:8001/lanes/demoLane/quarantine/r1
```

## API surface

| Method | Path | Purpose |
|---|---|---|
| POST | `/jobs` | Submit a lane (multipart: `lane`, `barcodes` JSON, `r1`, `r2`) |
| GET | `/jobs`, `/jobs/{id}` | Task status (`queued`/`processing`/`published`/`failed`) |
| POST | `/jobs/{id}/retry` | Requeue a failed task |
| GET | `/lanes`, `/lanes/{lane}` | Published lanes and per-bucket counts |
| GET | `/lanes/{lane}/samples/{sample_id}/{r1\|r2}` | Download a sample mate |
| GET | `/lanes/{lane}/quarantine/{r1\|r2}` | Download a quarantine mate |
| POST | `/snapshots` | Freeze the current published lane list |
| GET | `/snapshots`, `/snapshots/{id}` | List/inspect snapshots |
| GET | `/snapshots/{id}/samples/{sample_id}/{mate}` | Snapshot-scoped concatenated download |
| GET | `/snapshots/{id}/quarantine/{mate}` | Snapshot-scoped quarantine download |
| POST | `/merge/jobs` | Queue an overlap merge `{snapshot_id, sample_id}` (idempotent) |
| GET | `/merge/jobs`, `/merge/jobs/{id}` | Merge job status and merged/unmerged counts |
| POST | `/merge/jobs/{id}/retry` | Requeue a failed merge job |
| GET | `/merge/jobs/{id}/download` | Download the merged gzip FASTQ (completed jobs) |
| GET | `/health` | Liveness, instance id, current job, queue depth |

## Concurrency and crash safety

* **SQLite WAL** with `BEGIN IMMEDIATE` job claiming: exactly one instance
  ever picks up a given job; the other keeps serving requests.
* Each instance processes one lane at a time in a background thread, so a
  slow lane never blocks uploads, queries or downloads on the event loop.
* A job runs under an exclusive `flock` on its work directory. Startup
  recovery removes scratch/staging debris only when no lock is held, and
  never deletes directories that have a published-lane DB row — so another
  instance's in-flight work and all published data are preserved.
* Publish order: write/verify outputs in `work/`, move to
  `published/lanes/.tmp.<job>`, atomic `rename` to the final lane directory,
  then insert metadata in one DB transaction. Inputs are deleted only after
  the commit, so a crash at any point leaves no half-visible lane.
* Overlap merges follow the same pattern under their own
  `work/merges/<job>` flock and `published/merges/<job>` output: a crashed
  merge is marked failed and its scratch reaped on startup, while a live
  instance's flock and every completed merge row are always preserved.

## Memory

All I/O is streamed record-by-record (1 MiB network chunks, line-based gzip
reads, one open gzip writer per bucket). No file is ever buffered whole.
The automated streaming test processes hundreds of MiB of decompressed
data with RSS growth below 64 MiB.

## Layout

```
app/            FastAPI app, SQLite layer, streaming splitter, worker,
                overlap merger (merging.py)
scripts/        demo FASTQ generators, two-instance integration check
tests/          pytest suite (happy path, corruption, retry, recovery,
                barcode rules, size budget, streaming memory, merge
                algorithm/flow/recovery)
Dockerfile      python:3.11-slim, non-root user, /data volume, healthcheck
docker-compose.yml  two 128 MiB instances sharing one named volume
```

## Tests

```bash
python -m pytest tests/ -q
python scripts/integration_two_instances.py   # real two-process check
```
