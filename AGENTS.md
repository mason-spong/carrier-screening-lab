# AGENTS.md

This repository is a local, research-grade genomics workspace for building a **carrier-screening candidate list** from raw FASTQ files for two samples: `mason` and `hannah`.

## Start Here

When a coding agent opens this repo, read these in order:

1. `AGENTS.md` (this file)
2. `agents/PROJECT_CONTEXT.md`
3. `agents/QUALITY_STANDARDS.md`
4. `agents/WORKFLOW_PLAN.md`
5. `config/samples.tsv`

## Project Goal

Produce a high-confidence, reproducible, and auditable candidate carrier-screening output that includes:

- Candidate variants with evidence and confidence tiering
- Couple-level interpretation signals (shared gene carrier status risk)
- Explicit limitations and blind spots
- Clinical confirmation guidance for actionable candidates

Scope is **technical analysis and prioritization**, not diagnosis.

## Non-Negotiables

- Protect privacy: treat all genomic files as highly sensitive.
- Reproducibility first: pin tool/database versions and capture provenance.
- No silent assumptions: all thresholds and filters must be explicit.
- No clinical claims from raw FASTQ or unvalidated calls.
- Any actionable finding must be marked `research-only pending clinical confirmation`.

## Repository Conventions

- Track in git: code, configs, docs, lightweight manifests.
- Do not track in git: raw/intermediate/final genomic payloads (`.gitignore` already enforces this).
- Keep paths stable and sample-centric:
  - `data/raw/mason/`
  - `data/raw/hannah/`

## Execution Style

- Build pipeline in phases: smoke test first, then full run.
- Fail fast on QC gate violations.
- Prefer clear tabular outputs (`.tsv`) plus a concise report (`.md`).
- Preserve intermediate logs needed for auditability.

## Long-Run Utilities

Use these scripts for Mason-only execution and monitoring:

- `scripts/run_mason_full_analysis.sh`
  - End-to-end Mason technical run (QC, reference prep, alignment, calling, annotation, reports).
  - Uses checkpoint markers in `data/interim/mason/state/` for resumable phases.
- `scripts/run_mason_orchestrator.sh`
  - Supervises long execution.
  - Ensures full BWA reference index exists (`.amb/.ann/.bwt/.pac/.sa`) before launching full run.
  - Intended to be run with `caffeinate` for overnight stability.
- `scripts/start_mason_local.sh`
  - One-command local launcher for orchestrator with `nohup` + `caffeinate`.
  - Writes PID to `logs/mason_local_runner.pid` and validates process startup.
- `scripts/stop_mason_local.sh`
  - One-command local stop helper for orchestrator/pipeline worker processes.
  - Removes stale `logs/mason_local_runner.pid`.
- `scripts/serve_status_ui_maximalist.py`
  - Local read-only status dashboard with stage, phase/subphase progress, ETA, index-file status, and log tails.
  - Includes local run controls (`Start`, `Stop`, `Restart`) on the web page.
- `scripts/run_status_ui_maximalist_local.sh`
  - Convenience launcher for the local status UI server.

Primary logs for long runs:

- `logs/mason_orchestrator.log`
- `logs/mason_full_run.log`
- `logs/mason_local_runner.log` (if launched via local `nohup`)

## Execution-Context Rules (Critical)

To avoid localhost/process-lifecycle failures between agent and user sessions:

1. User-facing web UI must be started from the user's own terminal session.
   - Use `scripts/run_status_ui_maximalist_local.sh`.
2. Long-running analysis process of record must be started from the user's own terminal.
   - Use `scripts/start_mason_local.sh`.
3. Never assume `127.0.0.1` is shared between execution contexts.
   - Validate from the same context as the consumer (the user's browser).
4. Never treat BWA indexing as complete unless all files exist and are non-empty:
   - `.amb`, `.ann`, `.bwt`, `.pac`, `.sa`
5. Treat PID files as hints only.
   - Always confirm with `ps` before declaring a process live.

See `agents/OPERATIONS_RETROSPECTIVE.md` for detailed root-cause analysis and standard operating sequence.

## Required Deliverables

At minimum, create and maintain:

- `reports/run_summary.md`
- `results/qc_metrics.tsv`
- `results/alignment_metrics.tsv`
- `results/coverage_metrics.tsv`
- `results/annotated_variants.tsv` (or split by sample + joint)
- `results/carrier_candidates_couple.tsv`
- `reports/limitations_and_followup.md`

## Quality Bar

Use `agents/QUALITY_STANDARDS.md` as the source of truth for:

- QC pass/fail gates
- Variant confidence tiers
- Reportability rules
- Known blind spots (CNV/SV/repeats/pseudogene complexity)

## If Blocked

If data, tools, or references are missing:

1. Record the blocker in `reports/run_summary.md`.
2. Propose the smallest next action to unblock.
3. Do not fabricate outputs.
