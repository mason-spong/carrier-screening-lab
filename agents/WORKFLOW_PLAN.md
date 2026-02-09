# Workflow Plan

## Phase 0: Setup

1. Confirm sample sheet in `config/samples.tsv`.
2. Install/verify required tools.
3. Pin reference bundle and annotation database versions.

## Phase 1: Raw QC

1. Run FASTQ integrity checks.
2. Run `FastQC` and aggregate `MultiQC`.
3. Summarize read count, read length, and quality distributions.

Output:
- `results/qc_metrics.tsv`
- `reports/run_summary.md` (QC section)

## Phase 2: Alignment and Post-Processing

1. Align reads to GRCh38.
2. Sort/index alignments.
3. Mark duplicates.
4. Recalibrate base qualities if workflow uses BQSR.
5. Compute alignment and depth metrics.

Output:
- `results/alignment_metrics.tsv`
- `results/coverage_metrics.tsv`

## Phase 3: Variant Calling

1. Per-sample variant calling in gVCF mode.
2. Joint genotyping for both samples.
3. Apply quality filters with explicit thresholds.

Output:
- Joint filtered VCF (gitignored artifact path in `results/`)
- Filter summary in `reports/run_summary.md`

## Phase 4: Annotation

1. Annotate effects/transcripts.
2. Add clinical evidence fields (ClinVar review/significance).
3. Add population frequencies.
4. Normalize and generate analysis-friendly TSV.

Output:
- `results/annotated_variants.tsv`

## Phase 5: Carrier Candidate Extraction

1. Restrict to carrier-screening relevant genes/logic.
2. Apply couple-level inheritance filters.
3. Assign confidence tier (`Tier A/B/C/Not reportable`).
4. Produce confirmation guidance.

Output:
- `results/carrier_candidates_couple.tsv`
- `reports/limitations_and_followup.md`

## Phase 6: Final Reporting

1. Produce concise run summary with:
  - inputs
  - versions
  - pass/fail gates
  - candidate counts by tier
  - explicit limitations
2. Mark all actionable candidates as research-only pending clinical confirmation.

Output:
- `reports/run_summary.md`

