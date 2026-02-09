# Carrier Screening Lab

Research-grade local pipeline workspace for joint carrier-screening candidate analysis.

Agent-facing operating docs are in:
- `AGENTS.md`
- `agents/README.md`

## Samples
- `mason`: FASTQ files present.
- `hannah`: waiting for FASTQ upload.

## Layout
- `config/` sample sheet and run settings
- `workflow/` pipeline definitions
- `scripts/` helper scripts
- `data/raw/<sample>/` raw FASTQ inputs (gitignored)
- `data/interim/` temporary/intermediate artifacts (gitignored)
- `data/processed/` stable processed artifacts (gitignored)
- `results/` callsets and tables (gitignored)
- `reports/` human-readable reports (gitignored)

## Next steps
1. Add Hannah FASTQs under `data/raw/hannah/`.
2. Install tooling (`fastqc`, `multiqc`, `fastp`, `bwa-mem2`, `samtools`, `bcftools`).
3. Run a subset smoke test before full run.

## Long-run Execution and Monitoring
- Full Mason-only run: `scripts/run_mason_full_analysis.sh`
- Supervising orchestrator: `scripts/run_mason_orchestrator.sh`
- One-command local runner: `scripts/start_mason_local.sh`
- One-command local stop: `scripts/stop_mason_local.sh`
- Local dashboard server: `scripts/serve_status_ui_maximalist.py`
- Local dashboard launcher: `scripts/run_status_ui_maximalist_local.sh`

Typical local usage:
1. Start run supervisor with power guard:
   - `./scripts/start_mason_local.sh`
2. Start local status UI in a separate terminal:
   - `./scripts/run_status_ui_maximalist_local.sh 8788`
3. Open:
   - `http://127.0.0.1:8788/`
4. Optional: use UI controls (`Start`, `Stop`, `Restart`) from the web page.

Important:
- Launch both commands from your own terminal session (not an external agent session) so process lifetime and localhost visibility match your browser context.
