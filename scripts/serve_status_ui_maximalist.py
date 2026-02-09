#!/usr/bin/env python3
from __future__ import annotations

import datetime as dt
import json
import os
import re
import subprocess
import sys
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LOGS = ROOT / "logs"
ORCH_LOG = LOGS / "mason_orchestrator.log"
FULL_LOG = LOGS / "mason_full_run.log"
ANNOTATED_TSV = ROOT / "results" / "annotated_variants.tsv"
CANDIDATE_TSV = ROOT / "results" / "carrier_candidates_couple.tsv"
SAMPLES_TSV = ROOT / "config" / "samples.tsv"
STATE_DIR = ROOT / "data" / "interim" / "mason" / "state"
ALIGN_DIR = ROOT / "data" / "interim" / "mason" / "alignment"
REF_FA = ROOT / "data" / "refs" / "grch38" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
START_SCRIPT = ROOT / "scripts" / "start_mason_local.sh"
STOP_SCRIPT = ROOT / "scripts" / "stop_mason_local.sh"
PORT = int(os.environ.get("STATUS_UI_MAX_PORT", "8788"))


def run(cmd: str) -> str:
    try:
        out = subprocess.check_output(cmd, shell=True, text=True, stderr=subprocess.DEVNULL)
        return out.strip()
    except Exception:
        return ""


def run_script(path: Path) -> tuple[bool, str]:
    if not path.exists():
        return False, f"script not found: {path}"
    try:
        proc = subprocess.run(
            [str(path)],
            cwd=str(ROOT),
            text=True,
            capture_output=True,
            check=False,
        )
    except Exception as exc:
        return False, f"failed to execute {path.name}: {exc}"

    output = (proc.stdout or "").strip()
    err = (proc.stderr or "").strip()
    msg = output if output else (err if err else f"{path.name} exited {proc.returncode}")
    ok = proc.returncode == 0
    return ok, msg


def pid_of(pattern: str) -> str:
    return run(f"pgrep -f \"{pattern}\" | head -n1")


def pid_of_any(patterns: list[str]) -> str:
    for pat in patterns:
        pid = pid_of(pat)
        if pid:
            return pid
    return ""


def read_tail(path: Path, n: int) -> str:
    if not path.exists():
        return ""
    try:
        lines = path.read_text(errors="replace").splitlines()[-n:]
        return "\n".join(lines)
    except Exception:
        return ""


def parse_elapsed_to_seconds(text: str) -> int:
    text = text.strip()
    if not text:
        return 0
    days = 0
    if "-" in text:
        d, text = text.split("-", 1)
        try:
            days = int(d)
        except ValueError:
            days = 0
    parts = []
    for p in text.split(":"):
        try:
            parts.append(int(p))
        except ValueError:
            parts.append(0)
    if len(parts) == 3:
        h, m, s = parts
    elif len(parts) == 2:
        h, m, s = 0, parts[0], parts[1]
    else:
        h, m, s = 0, 0, parts[0] if parts else 0
    return days * 86400 + h * 3600 + m * 60 + s


def fmt_dur(sec: int) -> str:
    if sec < 0:
        sec = 0
    d, rem = divmod(sec, 86400)
    h, rem = divmod(rem, 3600)
    m, s = divmod(rem, 60)
    if d:
        return f"{d}d {h:02d}h {m:02d}m {s:02d}s"
    return f"{h:02d}h {m:02d}m {s:02d}s"


def clamp01(v: float) -> float:
    if v < 0:
        return 0.0
    if v > 1:
        return 1.0
    return v


def file_size(path: Path) -> int:
    try:
        if path.exists():
            return path.stat().st_size
    except Exception:
        pass
    return 0


def state_done(name: str) -> bool:
    return (STATE_DIR / f"{name}.done").exists()


def parse_samples() -> dict[str, str]:
    values: dict[str, str] = {}
    if not SAMPLES_TSV.exists():
        return values
    try:
        lines = SAMPLES_TSV.read_text(errors="replace").splitlines()
        for row in lines[1:]:
            cols = row.split("\t")
            if len(cols) < 5:
                continue
            sid = cols[0]
            if sid == "mason":
                values["mason_r1"] = cols[2]
                values["mason_r2"] = cols[3]
            if sid == "hannah":
                values["hannah_status"] = cols[4]
    except Exception:
        return values
    return values


def tier_counts(path: Path) -> dict[str, int]:
    counts = {"tier_a": 0, "tier_b": 0, "tier_c": 0, "not_reportable": 0}
    if not path.exists():
        return counts
    try:
        with path.open("r", errors="replace") as fh:
            next(fh, None)
            for raw in fh:
                cols = raw.rstrip("\n").split("\t")
                if len(cols) < 17:
                    continue
                tier = cols[16]
                if tier == "Tier A - High confidence candidate":
                    counts["tier_a"] += 1
                elif tier == "Tier B - Needs orthogonal confirmation":
                    counts["tier_b"] += 1
                elif tier == "Tier C - Low confidence / likely artifact":
                    counts["tier_c"] += 1
                elif tier == "Not reportable":
                    counts["not_reportable"] += 1
    except Exception:
        return counts
    return counts


def candidate_count(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        return max(0, len(path.read_text(errors="replace").splitlines()) - 1)
    except Exception:
        return 0


def last_phase_name(full_log_text: str) -> str:
    matches = re.findall(r"\[(?:\d{4}-\d{2}-\d{2} [^\]]+)\]\s+(Phase\s+\d+:\s+[^\n]+)", full_log_text)
    return matches[-1] if matches else ""


def phase3_alignment_fraction_from_sizes(mason_r1: Path, mason_r2: Path, raw_bam: Path) -> tuple[float, str]:
    in_total = file_size(mason_r1) + file_size(mason_r2)
    raw_size = file_size(raw_bam)
    if in_total <= 0:
        return 0.0, "input size unavailable"

    # Empirical envelope for BAM output growth from compressed FASTQ inputs.
    expected_final = int(in_total * 1.28)
    frac = clamp01(raw_size / max(expected_final, 1))
    msg = f"raw BAM {raw_size/1e9:.1f}G / est final {expected_final/1e9:.1f}G"
    return frac, msg


def phase3_sort_fraction_from_sizes(src: Path, out: Path) -> tuple[float, str]:
    src_size = file_size(src)
    out_size = file_size(out)
    if src_size <= 0:
        return 0.0, "source BAM size unavailable"
    frac = clamp01(out_size / src_size)
    return frac, f"{out.name} {out_size/1e9:.1f}G / src {src_size/1e9:.1f}G"


def progress_payload() -> dict[str, object]:
    now = dt.datetime.now().astimezone()
    now_txt = now.strftime("%Y-%m-%d %H:%M:%S %Z")
    sample_info = parse_samples()

    mason_r1 = Path(sample_info["mason_r1"]) if sample_info.get("mason_r1") else Path("/__missing_mason_r1__")
    mason_r2 = Path(sample_info["mason_r2"]) if sample_info.get("mason_r2") else Path("/__missing_mason_r2__")

    orch_pid = pid_of_any([
        "run_mason_orchestrator.sh",
        "bash scripts/run_mason_orchestrator.sh",
    ])
    bwa_index_pid = pid_of_any([
        f"bwa index {REF_FA}",
        "bwa index data/refs/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        "bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    ])
    full_pid = pid_of_any([
        "run_mason_full_analysis.sh",
        "bash scripts/run_mason_full_analysis.sh",
    ])

    # Phase/subphase process probes for finer progress.
    bwa_mem_pid = pid_of_any(["bwa mem -t", "bwa mem "])
    sam_view_pid = pid_of_any(["samtools view -@", "samtools view "])
    sort_name_pid = pid_of_any(["samtools sort -n", "mason.name.bam"])
    fixmate_pid = pid_of("samtools fixmate")
    sort_pos_pid = pid_of_any(["samtools sort -@", "mason.pos.bam"])
    markdup_pid = pid_of("samtools markdup")
    sam_index_pid = pid_of_any(["samtools index", "mason.markdup.bam.bai"])
    bcftools_pid = pid_of("bcftools")
    snpeff_pid = pid_of_any(["snpEff", "java .*snpEff"])

    orch_text = ORCH_LOG.read_text(errors="replace") if ORCH_LOG.exists() else ""
    full_text = FULL_LOG.read_text(errors="replace") if FULL_LOG.exists() else ""
    full_tail = read_tail(FULL_LOG, 40)

    failed = ("Operation timed out" in full_tail) or ("ERROR:" in full_tail and not full_pid)

    idx_files = {
        ".amb": REF_FA.with_suffix(REF_FA.suffix + ".amb"),
        ".ann": REF_FA.with_suffix(REF_FA.suffix + ".ann"),
        ".bwt": REF_FA.with_suffix(REF_FA.suffix + ".bwt"),
        ".pac": REF_FA.with_suffix(REF_FA.suffix + ".pac"),
        ".sa": REF_FA.with_suffix(REF_FA.suffix + ".sa"),
    }
    idx_state = {k: (p.exists() and file_size(p) > 0) for k, p in idx_files.items()}
    idx_score = int(100 * sum(1 for v in idx_state.values() if v) / 5)

    raw_bam = ALIGN_DIR / "mason.raw.bam"
    name_bam = ALIGN_DIR / "mason.name.bam"
    fixmate_bam = ALIGN_DIR / "mason.fixmate.bam"
    pos_bam = ALIGN_DIR / "mason.pos.bam"
    markdup_bam = ALIGN_DIR / "mason.markdup.bam"

    phase = "Idle"
    phase_detail = "No active pipeline process"
    phase_progress = 0.0
    phase_eta = "Unknown"
    overall_progress = 0.0
    stage = "IDLE VOID"

    # Overall weighting by phase
    phase_weights = {
        "Phase 1": (2, 10),
        "Phase 2": (10, 24),
        "Phase 3": (24, 64),
        "Phase 4": (64, 76),
        "Phase 5": (76, 89),
        "Phase 6": (89, 96),
        "Phase 7": (96, 99),
    }

    if bwa_index_pid:
        stage = "INDEX RAGE"
        phase = "Phase 2"
        phase_detail = "Building GRCh38 BWA index"
        iter_match = re.findall(r"BWTIncConstructFromPacked\]\s+(\d+)\s+iterations done", orch_text[-400000:])
        iter_done = int(iter_match[-1]) if iter_match else 0
        if iter_done > 0:
            phase_progress = clamp01(iter_done / 688.0)
            etime = run(f"ps -o etime= -p {bwa_index_pid}")
            elapsed = parse_elapsed_to_seconds(etime)
            eta_idx = int(elapsed * max(1, 688 - iter_done) / iter_done)
            phase_eta = f"~{fmt_dur(eta_idx)}"
        else:
            phase_progress = 0.02
            phase_eta = "estimating"
    elif full_pid:
        stage = "PIPELINE OVERDRIVE"
        last_phase = last_phase_name(full_text)
        if last_phase.startswith("Phase 1"):
            phase = "Phase 1"
            phase_progress = 0.95 if state_done("phase1_qc") else 0.55
            phase_detail = "FASTQ integrity + FastQC/MultiQC"
        elif last_phase.startswith("Phase 2"):
            phase = "Phase 2"
            phase_progress = 0.95 if state_done("phase2_ref") else 0.45
            phase_detail = "Reference setup"
        elif last_phase.startswith("Phase 3") or bwa_mem_pid or sam_view_pid or sort_name_pid or fixmate_pid or sort_pos_pid or markdup_pid:
            phase = "Phase 3"
            if bwa_mem_pid or sam_view_pid:
                frac, msg = phase3_alignment_fraction_from_sizes(mason_r1, mason_r2, raw_bam)
                phase_progress = 0.76 * frac
                phase_detail = f"BWA MEM stream -> raw BAM ({msg})"
            elif sort_name_pid:
                frac, msg = phase3_sort_fraction_from_sizes(raw_bam, name_bam)
                phase_progress = 0.76 + (0.10 * frac)
                phase_detail = f"Name-sort BAM ({msg})"
            elif fixmate_pid:
                frac, msg = phase3_sort_fraction_from_sizes(name_bam, fixmate_bam)
                phase_progress = 0.86 + (0.06 * frac)
                phase_detail = f"Fixmate ({msg})"
            elif sort_pos_pid:
                frac, msg = phase3_sort_fraction_from_sizes(fixmate_bam if file_size(fixmate_bam) else name_bam, pos_bam)
                phase_progress = 0.92 + (0.04 * frac)
                phase_detail = f"Position sort ({msg})"
            elif markdup_pid:
                frac, msg = phase3_sort_fraction_from_sizes(pos_bam if file_size(pos_bam) else raw_bam, markdup_bam)
                phase_progress = 0.96 + (0.03 * frac)
                phase_detail = f"Mark duplicates ({msg})"
            elif sam_index_pid:
                phase_progress = 0.995
                phase_detail = "Indexing final BAM"
            else:
                phase_progress = 0.80
                phase_detail = "Alignment/post-processing running"
            etime = run(f"ps -o etime= -p {full_pid}")
            elapsed = parse_elapsed_to_seconds(etime)
            phase_eta = f"~{fmt_dur(max(1800, 12 * 3600 - elapsed))}"
        elif last_phase.startswith("Phase 4") or bcftools_pid:
            phase = "Phase 4"
            phase_progress = 0.55
            phase_detail = "Variant calling (bcftools)"
        elif last_phase.startswith("Phase 5") or snpeff_pid:
            phase = "Phase 5"
            phase_progress = 0.55
            phase_detail = "Annotation (snpEff + ClinVar)"
        elif last_phase.startswith("Phase 6"):
            phase = "Phase 6"
            phase_progress = 0.7
            phase_detail = "Building result tables"
        elif last_phase.startswith("Phase 7"):
            phase = "Phase 7"
            phase_progress = 0.6
            phase_detail = "Report generation"
        else:
            phase = "Initializing"
            phase_progress = 0.1
            phase_detail = "Pipeline process alive; awaiting phase banner"

        if phase in phase_weights:
            lo, hi = phase_weights[phase]
            overall_progress = lo + (hi - lo) * clamp01(phase_progress)
        else:
            overall_progress = max(2.0, 20.0 * clamp01(phase_progress))

        phase_eta = phase_eta if phase_eta != "Unknown" else "estimating"
    elif "Full pipeline exit status: 0" in orch_text:
        stage = "RUN LANDED"
        phase = "Completed"
        phase_detail = "Pipeline finished successfully"
        phase_progress = 1.0
        phase_eta = "Done"
        overall_progress = 100.0
    elif failed:
        stage = "RUN FAILED"
        phase = "Failed"
        phase_detail = "Write timeout in latest run; inspect full log tail"
        phase_progress = 0.0
        phase_eta = "Stopped"
        overall_progress = 0.0
    elif orch_pid:
        stage = "BOOTSTRAPPING"
        phase = "Init"
        phase_detail = "Orchestrator alive; waiting on worker process"
        phase_progress = 0.1
        phase_eta = "Preparing"
        overall_progress = 1.0

    tiers = tier_counts(ANNOTATED_TSV)
    couple_candidates = candidate_count(CANDIDATE_TSV)

    return {
        "now": now_txt,
        "stage": stage,
        "detail": f"{phase}: {phase_detail}",
        "eta": phase_eta,
        "progress": round(clamp01(overall_progress / 100.0) * 100.0, 2),
        "phase": phase,
        "phase_progress": round(100.0 * clamp01(phase_progress), 2),
        "phase_detail": phase_detail,
        "orch_pid": orch_pid or "none",
        "bwa_pid": bwa_index_pid or "none",
        "full_pid": full_pid or "none",
        "orchestrator_log": read_tail(ORCH_LOG, 26),
        "full_log": read_tail(FULL_LOG, 26),
        "idx": idx_state,
        "idx_score": idx_score,
        "tiers": tiers,
        "couple_candidates": couple_candidates,
        "has_active_run": bool(bwa_index_pid or full_pid),
        "hannah_status": sample_info.get("hannah_status", "unknown"),
    }


def control_action(action: str) -> dict[str, object]:
    action = (action or "").strip().lower()
    if action == "start":
        ok, msg = run_script(START_SCRIPT)
    elif action == "stop":
        ok, msg = run_script(STOP_SCRIPT)
    elif action == "restart":
        ok_stop, msg_stop = run_script(STOP_SCRIPT)
        ok_start, msg_start = run_script(START_SCRIPT)
        ok = ok_start
        msg = f"stop: {msg_stop}\nstart: {msg_start}" if ok_stop else f"stop (non-fatal): {msg_stop}\nstart: {msg_start}"
    else:
        return {
            "ok": False,
            "message": f"unsupported action: {action}",
            "status": progress_payload(),
        }

    safe_note = "No data deletion: controls do not remove interim files or checkpoints."
    return {
        "ok": ok,
        "message": f"{msg}\n{safe_note}",
        "status": progress_payload(),
    }


def html_page() -> str:
    return """<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
  <title>Mason Monitor MAX</title>
  <style>
    :root {
      --bg0: #07040f;
      --bg1: #130824;
      --ink: #f8fbff;
      --muted: #b8bfd8;
      --acid: #30f2a2;
      --hot: #ff4fd8;
      --sun: #ffbf3c;
      --sky: #4ad1ff;
      --line: rgba(255,255,255,0.16);
    }
    * { box-sizing: border-box; }
    html {
      min-height: 100%;
      background-color: #07040f;
      overscroll-behavior-y: none;
      overflow-x: hidden;
      background:
        radial-gradient(1100px 700px at -10% -20%, #ff4fd830 0, transparent 55%),
        radial-gradient(900px 640px at 110% -10%, #30f2a22c 0, transparent 58%),
        radial-gradient(1200px 900px at 50% 120%, #4ad1ff26 0, transparent 62%),
        repeating-linear-gradient(135deg, #ffffff07 0 2px, transparent 2px 16px),
        linear-gradient(160deg, var(--bg0), var(--bg1));
    }
    body {
      margin: 0;
      min-height: 100vh;
      background-color: #07040f;
      position: relative;
      overscroll-behavior-y: none;
      -webkit-overflow-scrolling: touch;
      overflow-x: hidden;
      overflow-y: auto;
      font-family: \"Arial Black\", \"Impact\", \"Avenir Next\", sans-serif;
      color: var(--ink);
      background:
        radial-gradient(1100px 700px at -10% -20%, #ff4fd830 0, transparent 55%),
        radial-gradient(900px 640px at 110% -10%, #30f2a22c 0, transparent 58%),
        radial-gradient(1200px 900px at 50% 120%, #4ad1ff26 0, transparent 62%),
        repeating-linear-gradient(135deg, #ffffff07 0 2px, transparent 2px 16px),
        linear-gradient(160deg, var(--bg0), var(--bg1));
    }
    body::before {
      content: "";
      position: fixed;
      inset: 0;
      z-index: -1;
      background:
        radial-gradient(1100px 700px at -10% -20%, #ff4fd830 0, transparent 55%),
        radial-gradient(900px 640px at 110% -10%, #30f2a22c 0, transparent 58%),
        radial-gradient(1200px 900px at 50% 120%, #4ad1ff26 0, transparent 62%),
        repeating-linear-gradient(135deg, #ffffff07 0 2px, transparent 2px 16px),
        linear-gradient(160deg, var(--bg0), var(--bg1));
      pointer-events: none;
    }
    .wrap {
      max-width: 1320px;
      margin: 0 auto;
      display: grid;
      gap: 12px;
      min-height: 100vh;
      padding: 14px;
    }
    .hero {
      border: 2px solid var(--line);
      border-radius: 20px;
      padding: 14px;
      background: linear-gradient(95deg, #ffffff08, #ffffff02);
      backdrop-filter: blur(1px);
      box-shadow: 0 0 0 2px #ffffff08 inset, 0 10px 50px #00000066;
    }
    .h {
      margin: 10px 0 4px;
      line-height: 1;
      font-size: clamp(30px, 7vw, 88px);
      letter-spacing: 0.03em;
      text-transform: uppercase;
      text-shadow: 0 0 18px #30f2a255, 0 0 46px #ff4fd866;
    }
    .sub { color: var(--muted); font-size: 13px; letter-spacing: 0.03em; }
    .grid { display: grid; gap: 12px; grid-template-columns: repeat(3, minmax(0, 1fr)); align-items: start; }
    .card {
      border: 2px solid var(--line);
      border-radius: 16px;
      background: linear-gradient(150deg, #ffffff10, #ffffff03);
      padding: 12px;
      min-width: 0;
    }
    .k { font-size: 11px; text-transform: uppercase; letter-spacing: 0.09em; color: #d4d9ea; }
    .v { margin-top: 8px; font-size: 30px; line-height: 1.06; }
    .tiny { margin-top: 8px; font-size: 13px; color: var(--muted); }
    .meter {
      margin-top: 10px;
      height: 16px;
      border-radius: 999px;
      border: 1px solid #ffffff4a;
      overflow: hidden;
      background: #0000005c;
    }
    .fill {
      height: 100%;
      width: 0%;
      background: linear-gradient(90deg, var(--acid), var(--sky), var(--hot));
      box-shadow: 0 0 18px #4ad1ff88;
      transition: width 0.45s ease;
    }
    .phase-fill {
      background: linear-gradient(90deg, var(--sun), var(--hot));
      box-shadow: 0 0 18px #ffbf3c88;
    }
    .strip { display: grid; gap: 10px; grid-template-columns: repeat(4, minmax(0, 1fr)); }
    .chip {
      border: 2px solid var(--line);
      border-radius: 14px;
      padding: 10px;
      background: #ffffff08;
      min-height: 90px;
    }
    .num { font-size: 32px; margin-top: 4px; }
    .table {
      width: 100%; border-collapse: collapse; font-size: 13px; margin-top: 7px;
    }
    .table th, .table td {
      border-bottom: 1px solid #ffffff26; padding: 7px 6px; text-align: left;
    }
    .logs { display: grid; gap: 12px; grid-template-columns: repeat(2, minmax(0, 1fr)); }
    pre {
      margin: 8px 0 0;
      padding: 12px;
      border: 2px solid #ffffff1f;
      border-radius: 12px;
      background: #04040add;
      min-height: 250px;
      max-height: 350px;
      overflow: auto;
      white-space: pre-wrap;
      word-break: break-word;
      color: #d9e1ff;
      font-family: Menlo, Consolas, monospace;
      font-size: 12px;
      line-height: 1.35;
    }
    .ok { color: #7bffb8; }
    .no { color: #ffb2ca; }
    .controls { display: flex; gap: 8px; margin-top: 10px; flex-wrap: wrap; }
    .btn {
      border: 1px solid #ffffff4a;
      background: #ffffff12;
      color: #f8fbff;
      border-radius: 10px;
      font-size: 12px;
      font-weight: 700;
      letter-spacing: 0.03em;
      text-transform: uppercase;
      padding: 8px 10px;
      cursor: pointer;
    }
    .btn:hover { background: #ffffff22; }
    .btn:disabled { opacity: 0.5; cursor: not-allowed; }
    .ctrl-msg {
      margin-top: 10px;
      white-space: pre-wrap;
      font-size: 12px;
      color: #d9e1ff;
      border: 1px dashed #ffffff33;
      border-radius: 10px;
      padding: 8px;
      max-height: 120px;
      overflow: auto;
    }
    .ctrl-help {
      margin-top: 8px;
      font-size: 11px;
      line-height: 1.35;
      color: #d4d9ea;
      border: 1px dashed #ffffff2a;
      border-radius: 10px;
      padding: 8px;
      white-space: pre-wrap;
    }
    @media (max-width: 980px) {
      .grid { grid-template-columns: 1fr; }
      .strip { grid-template-columns: 1fr 1fr; }
      .logs { grid-template-columns: 1fr; }
    }
  </style>
</head>
<body>
  <div class=\"wrap\">
    <section class=\"hero\">
      <div class=\"h\" id=\"stage\">Loading</div>
      <div class=\"sub\" id=\"detail\">Waiting for status payload</div>
    </section>

    <section class=\"grid\">
      <div class=\"card\">
        <div class=\"k\">Run Momentum</div>
        <div class=\"v\" id=\"eta\">ETA --</div>
        <div class=\"tiny\">Overall progress</div>
        <div class=\"meter\"><div class=\"fill\" id=\"fill\"></div></div>
        <div class=\"tiny\"><b id=\"progress\">0%</b></div>
        <div class=\"tiny\" style=\"margin-top:10px\">Current phase progress</div>
        <div class=\"meter\"><div class=\"fill phase-fill\" id=\"phaseFill\"></div></div>
        <div class=\"tiny\"><b id=\"phase\">Unknown</b> · <b id=\"phaseProgress\">0%</b></div>
        <div class=\"tiny\" id=\"phaseDetail\">--</div>
      </div>
      <div class=\"card\">
        <div class=\"k\">Process IDs</div>
        <div class=\"tiny\">orchestrator: <b id=\"orch\">none</b></div>
        <div class=\"tiny\">bwa index: <b id=\"bwa\">none</b></div>
        <div class=\"tiny\">full run: <b id=\"full\">none</b></div>
        <div class=\"tiny\">hannah status: <b id=\"hannah\">unknown</b></div>
        <div class=\"tiny\">refresh: <b id=\"now\">--</b></div>
        <div class=\"controls\">
          <button class=\"btn\" id=\"btnStart\" onclick=\"control('start')\">Start</button>
          <button class=\"btn\" id=\"btnStop\" onclick=\"control('stop')\">Stop</button>
          <button class=\"btn\" id=\"btnRestart\" onclick=\"control('restart')\">Restart</button>
        </div>
        <div class=\"ctrl-help\">Start: launches orchestrator if not running.\nStop: terminates active run processes.\nRestart: Stop then Start.\nNo data deletion: these controls do NOT remove interim BAMs or checkpoint .done files.</div>
        <div class=\"ctrl-msg\" id=\"ctrlMsg\">No control action yet.</div>
      </div>
      <div class=\"card\">
        <div class=\"k\">Index Completeness</div>
        <div class=\"v\"><span id=\"idxscore\">0</span>%</div>
        <table class=\"table\">
          <thead><tr><th>File</th><th>State</th></tr></thead>
          <tbody id=\"idxrows\"></tbody>
        </table>
      </div>
    </section>

    <section class=\"strip\">
      <div class=\"chip\"><div class=\"k\">Tier A</div><div class=\"num\" id=\"ta\">0</div></div>
      <div class=\"chip\"><div class=\"k\">Tier B</div><div class=\"num\" id=\"tb\">0</div></div>
      <div class=\"chip\"><div class=\"k\">Tier C</div><div class=\"num\" id=\"tc\">0</div></div>
      <div class=\"chip\"><div class=\"k\">Couple Rows</div><div class=\"num\" id=\"cc\">0</div></div>
    </section>

    <section class=\"logs\">
      <div class=\"card\">
        <div class=\"k\">Orchestrator Tail</div>
        <pre id=\"orchlog\"></pre>
      </div>
      <div class=\"card\">
        <div class=\"k\">Full Run Tail</div>
        <pre id=\"fulllog\"></pre>
      </div>
    </section>
  </div>

  <script>
    const byId = (id) => document.getElementById(id);
    const setText = (id, value) => byId(id).textContent = value;
    const controlButtons = ["btnStart", "btnStop", "btnRestart"];
    const asNum = (value) => {
      const n = Number(value);
      return Number.isFinite(n) ? n : 0;
    };
    const fmtPct = (value) => `${asNum(value).toFixed(2)}%`;

    function render(data) {
      setText("stage", data.stage || "Unknown");
      setText("detail", data.detail || "No detail");
      setText("eta", `ETA ${data.eta || "Unknown"}`);
      setText("progress", fmtPct(data.progress));
      setText("phase", data.phase || "Unknown");
      setText("phaseProgress", fmtPct(data.phase_progress));
      setText("phaseDetail", data.phase_detail || "--");
      setText("orch", data.orch_pid || "none");
      setText("bwa", data.bwa_pid || "none");
      setText("full", data.full_pid || "none");
      setText("hannah", data.hannah_status || "unknown");
      setText("now", data.now || "--");
      setText("idxscore", String(data.idx_score || 0));

      byId("fill").style.width = `${asNum(data.progress)}%`;
      byId("phaseFill").style.width = `${asNum(data.phase_progress)}%`;

      const idx = data.idx || {};
      byId("idxrows").innerHTML = [".amb", ".ann", ".bwt", ".pac", ".sa"]
        .map((k) => {
          const ok = !!idx[k];
          return `<tr><td>${k}</td><td class=\"${ok ? "ok" : "no"}\">${ok ? "yes" : "no"}</td></tr>`;
        })
        .join("");

      const tiers = data.tiers || {};
      setText("ta", String(tiers.tier_a || 0));
      setText("tb", String(tiers.tier_b || 0));
      setText("tc", String(tiers.tier_c || 0));
      setText("cc", String(data.couple_candidates || 0));

      byId("orchlog").textContent = data.orchestrator_log || "";
      byId("fulllog").textContent = data.full_log || "";
    }

    function setControlsEnabled(enabled) {
      controlButtons.forEach((id) => {
        const el = byId(id);
        if (el) el.disabled = !enabled;
      });
    }

    async function control(action) {
      if (action === "restart") {
        const ok = confirm("Restart = Stop then Start. This does NOT delete interim data or checkpoints. Continue?");
        if (!ok) return;
      }
      if (action === "stop") {
        const ok = confirm("Stop will terminate active run processes only. Interim data and checkpoints are kept. Continue?");
        if (!ok) return;
      }
      setControlsEnabled(false);
      setText("ctrlMsg", `Running action: ${action} ...`);
      try {
        const res = await fetch("/api/control", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ action })
        });
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        const payload = await res.json();
        const ok = payload.ok ? "OK" : "FAILED";
        setText("ctrlMsg", `${ok}: ${payload.message || "no message"}`);
        if (payload.status) render(payload.status);
      } catch (err) {
        setText("ctrlMsg", `FAILED: ${err.message}`);
      } finally {
        setControlsEnabled(true);
      }
    }

    async function tick() {
      try {
        const res = await fetch('/api/status', { cache: 'no-store' });
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        render(await res.json());
      } catch (err) {
        render({
          stage: "IDLE VOID",
          detail: `status fetch failed: ${err.message}`,
          eta: "Unknown",
          progress: 0,
          phase: "Unknown",
          phase_progress: 0,
          phase_detail: "--",
          orch_pid: "none",
          bwa_pid: "none",
          full_pid: "none",
          hannah_status: "unknown",
          now: "--",
          idx_score: 0,
          idx: { ".amb": false, ".ann": false, ".bwt": false, ".pac": false, ".sa": false },
          tiers: { tier_a: 0, tier_b: 0, tier_c: 0 },
          couple_candidates: 0,
          orchestrator_log: "",
          full_log: ""
        });
      }
    }

    tick();
    setInterval(tick, 5000);
  </script>
</body>
</html>
"""


class Handler(BaseHTTPRequestHandler):
    def _send(self, code: int, content: bytes, ctype: str) -> None:
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(content)))
        self.end_headers()
        self.wfile.write(content)

    def do_GET(self) -> None:
        if self.path in ("/", "/index.html"):
            self._send(200, html_page().encode("utf-8"), "text/html; charset=utf-8")
            return
        if self.path == "/api/status":
            self._send(200, json.dumps(progress_payload()).encode("utf-8"), "application/json; charset=utf-8")
            return
        self._send(404, b"not found", "text/plain; charset=utf-8")

    def do_POST(self) -> None:
        if self.path != "/api/control":
            self._send(404, b"not found", "text/plain; charset=utf-8")
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError:
            length = 0
        raw = self.rfile.read(length) if length > 0 else b"{}"
        try:
            body = json.loads(raw.decode("utf-8"))
        except Exception:
            body = {}
        action = str(body.get("action", "")).strip().lower()
        payload = control_action(action)
        status = 200 if payload.get("ok") else 400
        self._send(status, json.dumps(payload).encode("utf-8"), "application/json; charset=utf-8")

    def do_HEAD(self) -> None:
        if self.path in ("/", "/index.html"):
            body = html_page().encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            return
        if self.path == "/api/status":
            body = json.dumps(progress_payload()).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            return
        self.send_response(404)
        self.end_headers()

    def log_message(self, fmt: str, *args) -> None:
        return


def main() -> None:
    global PORT
    if len(sys.argv) > 1:
        try:
            PORT = int(sys.argv[1])
        except ValueError:
            pass
    server = HTTPServer(("127.0.0.1", PORT), Handler)
    print(f"Maximalist monitor on http://127.0.0.1:{PORT}/")
    server.serve_forever()


if __name__ == "__main__":
    main()
