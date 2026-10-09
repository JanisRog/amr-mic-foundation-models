#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
import random
import re
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS = ",".join([
    "run_accession",
    "sample_accession",
    "experiment_accession",
    "study_accession",
    "library_layout",
    "instrument_platform",
    "instrument_model",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "submitted_ftp",
    "submitted_md5",
    "submitted_bytes",
])


def safe_name(text: str) -> str:
    text = str(text).strip()
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    text = text.strip("._")
    return text or "sample"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        first = handle.readline().rstrip("\n")
        if not first:
            return []
        handle.seek(0)
        first_fields = first.split("\t")
        if first_fields[0] in {"sample_id", "Sample_IDs"}:
            rows = []
            for row in csv.DictReader(handle, delimiter="\t"):
                rows.append({
                    "sample_id": (row.get("sample_id") or row.get("Sample_IDs") or "").strip(),
                    "ena_sample": (row.get("ena_sample") or row.get("ENA_SAMPLE") or "").strip(),
                    "status": (row.get("status") or "").strip(),
                    "reason": (row.get("reason") or "").strip(),
                })
            return rows

        rows = []
        reader = csv.reader(handle, delimiter="\t")
        for parts in reader:
            if not parts:
                continue
            rows.append({
                "sample_id": parts[0].strip() if len(parts) > 0 else "",
                "ena_sample": parts[1].strip() if len(parts) > 1 else "",
                "status": parts[2].strip() if len(parts) > 2 else "",
                "reason": parts[3].strip() if len(parts) > 3 else "",
            })
        return rows


def load_row(metadata_tsv: Path, task_id: int) -> dict[str, str]:
    rows = read_rows(metadata_tsv)
    if task_id < 1 or task_id > len(rows):
        raise SystemExit(f"No metadata row for SLURM_ARRAY_TASK_ID={task_id}; n_rows={len(rows)}")
    return rows[task_id - 1]


def fetch_run_report(ena_sample: str) -> list[dict[str, str]]:
    params = {
        "accession": ena_sample,
        "result": "read_run",
        "fields": FIELDS,
        "format": "tsv",
    }
    url = ENA_FILEREPORT + "?" + urllib.parse.urlencode(params)
    last_err: Exception | None = None
    for attempt in range(1, 8):
        try:
            with urllib.request.urlopen(url, timeout=90) as response:
                text = response.read().decode("utf-8", errors="replace")
            lines = [x for x in text.splitlines() if x.strip()]
            if not lines:
                return []
            reader = csv.DictReader(lines, delimiter="\t")
            return [{k: (v or "").strip() for k, v in row.items()} for row in reader]
        except Exception as exc:  # noqa: BLE001
            last_err = exc
            time.sleep(min(120, 10 * attempt) + random.randint(0, 10))
    raise RuntimeError(f"ENA read_run API failed for {ena_sample}: {last_err}")


def ftp_to_https(raw: str) -> str:
    u = raw.strip()
    if not u:
        return ""
    if u.startswith("ftp://ftp.sra.ebi.ac.uk/"):
        return "https://ftp.sra.ebi.ac.uk/" + u.split("ftp://ftp.sra.ebi.ac.uk/", 1)[1]
    if u.startswith("ftp://ftp.ebi.ac.uk/"):
        return "https://ftp.ebi.ac.uk/" + u.split("ftp://ftp.ebi.ac.uk/", 1)[1]
    if u.startswith("ftp.sra.ebi.ac.uk/") or u.startswith("ftp.ebi.ac.uk/"):
        return "https://" + u
    if u.startswith("http://"):
        return "https://" + u.split("http://", 1)[1]
    return u


def split_values(value: str) -> list[str]:
    return [x.strip() for x in re.split(r"[;,]", value or "") if x.strip()]


def is_fastq_url(url: str) -> bool:
    lower = url.lower()
    return lower.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq")) or "fastq" in Path(lower).name


def read_index_from_name(name: str, pos: int) -> str:
    lower = name.lower()
    if re.search(r"(^|[_\-.])r?1([_\-.]|\.f(ast)?q)", lower):
        return "1"
    if re.search(r"(^|[_\-.])r?2([_\-.]|\.f(ast)?q)", lower):
        return "2"
    return str(pos)


def collect_fastq_records(report: list[dict[str, str]]) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    for row in report:
        run = row.get("run_accession", "run") or "run"
        urls = split_values(row.get("fastq_ftp", ""))
        md5s = split_values(row.get("fastq_md5", ""))
        bytes_ = split_values(row.get("fastq_bytes", ""))
        source = "fastq_ftp"
        if not urls:
            urls = [u for u in split_values(row.get("submitted_ftp", "")) if is_fastq_url(u)]
            md5s = split_values(row.get("submitted_md5", ""))
            bytes_ = split_values(row.get("submitted_bytes", ""))
            source = "submitted_ftp"
        for i, raw in enumerate(urls, start=1):
            url = ftp_to_https(raw)
            if not is_fastq_url(url):
                continue
            name = Path(urllib.parse.urlparse(url).path).name or f"{run}_{i}.fastq.gz"
            records.append({
                "run_accession": run,
                "url": url,
                "filename": name,
                "read_index": read_index_from_name(name, i),
                "md5": md5s[i - 1] if i - 1 < len(md5s) else "",
                "bytes": bytes_[i - 1] if i - 1 < len(bytes_) else "",
                "source_field": source,
                "library_layout": row.get("library_layout", ""),
                "instrument_platform": row.get("instrument_platform", ""),
                "instrument_model": row.get("instrument_model", ""),
            })
    return records


def run(cmd: list[str]) -> None:
    print("[CMD]", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def clean_tmp(path: Path) -> None:
    try:
        if path.exists():
            path.unlink()
    except OSError:
        pass


def python_download(url: str, tmp: Path) -> None:
    with urllib.request.urlopen(url, timeout=1200) as response, tmp.open("wb") as out:
        shutil.copyfileobj(response, out)


def download(url: str, dest: Path, force: bool = False) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0 and not force:
        return
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    last_err: Exception | None = None
    for attempt in range(1, 11):
        clean_tmp(tmp)
        try:
            if shutil.which("curl"):
                run([
                    "curl", "-L", "--fail",
                    "--retry", "12", "--retry-delay", "20", "--retry-max-time", "1800",
                    "--connect-timeout", "60", "--max-time", "1800",
                    "-o", str(tmp), url,
                ])
            elif shutil.which("wget"):
                run(["wget", "--tries=12", "--waitretry=20", "--timeout=1800", "-O", str(tmp), url])
            else:
                python_download(url, tmp)
            if not tmp.exists() or tmp.stat().st_size == 0:
                raise RuntimeError("download produced empty file")
            tmp.replace(dest)
            return
        except Exception as exc:  # noqa: BLE001
            last_err = exc
            clean_tmp(tmp)
            wait = min(300, 20 * attempt) + random.randint(0, 20)
            print(f"[WARN] download_attempt={attempt} failed wait_seconds={wait} error={exc}", flush=True)
            time.sleep(wait)
    clean_tmp(tmp)
    try:
        python_download(url, tmp)
        if not tmp.exists() or tmp.stat().st_size == 0:
            raise RuntimeError("python fallback produced empty file")
        tmp.replace(dest)
        return
    except Exception as exc:  # noqa: BLE001
        clean_tmp(tmp)
        raise RuntimeError(f"download failed after retries: {last_err}; python_fallback={exc}") from exc


def write_status(path: Path, data: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(data.keys()))
        writer.writeheader()
        writer.writerow(data)


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    metadata_tsv = Path(os.environ["METADATA_TSV"])
    outroot = Path(os.environ.get("OUTROOT") or os.path.join(os.environ["MTB"], "reads"))
    task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", "1"))
    force = os.environ.get("FORCE", "0") == "1"

    row = load_row(metadata_tsv, task_id)
    sample_id = row.get("sample_id", "")
    ena_sample = row.get("ena_sample", "")
    sample_safe = safe_name(sample_id)

    sample_dir = outroot / "fastq_by_sample" / sample_safe
    status_path = outroot / "status" / f"{sample_safe}.tsv"
    manifest_path = outroot / "manifests" / f"{sample_safe}.fastq_manifest.tsv"

    base_status = {
        "sample_id": sample_id,
        "sample_safe": sample_safe,
        "ena_sample": ena_sample,
        "status": "",
        "reason": "",
        "n_runs": "0",
        "n_fastq_files": "0",
        "sample_dir": str(sample_dir),
        "manifest": str(manifest_path),
    }

    if not ena_sample:
        base_status.update({"status": "failed", "reason": "missing_ena_sample"})
        write_status(status_path, base_status)
        print(json.dumps(base_status, indent=2))
        raise SystemExit(2)

    report = fetch_run_report(ena_sample)
    records = collect_fastq_records(report)
    if not records:
        base_status.update({"status": "missing", "reason": "no_read_fastq_url", "n_runs": str(len(report))})
        write_status(status_path, base_status)
        print(json.dumps(base_status, indent=2))
        raise SystemExit(2)

    manifest_rows: list[dict[str, str]] = []
    try:
        for rec in records:
            run_acc = safe_name(rec["run_accession"])
            dest = sample_dir / run_acc / rec["filename"]
            download(rec["url"], dest, force=force)
            manifest_rows.append({
                "sample_id": sample_id,
                "ena_sample": ena_sample,
                "run_accession": rec["run_accession"],
                "read_index": rec["read_index"],
                "library_layout": rec["library_layout"],
                "instrument_platform": rec["instrument_platform"],
                "instrument_model": rec["instrument_model"],
                "source_field": rec["source_field"],
                "url": rec["url"],
                "md5": rec["md5"],
                "bytes": rec["bytes"],
                "fastq": str(dest),
            })
        write_manifest(manifest_path, manifest_rows)
    except Exception as exc:  # noqa: BLE001
        if manifest_rows:
            write_manifest(manifest_path, manifest_rows)
        base_status.update({
            "status": "failed",
            "reason": f"download_failed: {exc}",
            "n_runs": str(len({r.get("run_accession", "") for r in records})),
            "n_fastq_files": str(len(manifest_rows)),
        })
        write_status(status_path, base_status)
        print(json.dumps(base_status, indent=2))
        raise SystemExit(2) from exc

    base_status.update({
        "status": "complete",
        "reason": "",
        "n_runs": str(len({r.get("run_accession", "") for r in records})),
        "n_fastq_files": str(len(manifest_rows)),
    })
    write_status(status_path, base_status)
    print(json.dumps(base_status, indent=2))


if __name__ == "__main__":
    main()
