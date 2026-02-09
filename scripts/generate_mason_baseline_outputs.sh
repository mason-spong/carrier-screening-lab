#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

RUN_TS_UTC="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
RUN_TS_LOCAL="$(date +"%Y-%m-%d %H:%M:%S %Z")"

mkdir -p results reports logs

MASON_R1="$(awk -F '\t' '$1=="mason" {print $3}' config/samples.tsv)"
MASON_R2="$(awk -F '\t' '$1=="mason" {print $4}' config/samples.tsv)"
HANNAH_STATUS="$(awk -F '\t' '$1=="hannah" {print $5}' config/samples.tsv)"

if [[ -z "$MASON_R1" || -z "$MASON_R2" ]]; then
  echo "Failed to resolve mason FASTQ paths from config/samples.tsv" >&2
  exit 1
fi

if [[ ! -f "$MASON_R1" || ! -f "$MASON_R2" ]]; then
  echo "Mason FASTQ files not found: $MASON_R1 $MASON_R2" >&2
  exit 1
fi

have_tool() {
  if command -v "$1" >/dev/null 2>&1; then
    echo "yes"
  else
    echo "no"
  fi
}

tool_version() {
  local tool="$1"
  if command -v "$tool" >/dev/null 2>&1; then
    case "$tool" in
      bwa)
        LC_ALL=C LANG=C bwa 2>&1 | awk 'NF {print; exit}'
        ;;
      fastqc)
        LC_ALL=C LANG=C fastqc --version 2>&1 | head -n 1
        ;;
      snpEff)
        LC_ALL=C LANG=C snpEff -version 2>&1 | head -n 1
        ;;
      vep)
        if LC_ALL=C LANG=C vep --help >/dev/null 2>&1; then
          echo "vep runnable"
        else
          echo "vep present but not runnable"
        fi
        ;;
      *)
        LC_ALL=C LANG=C "$tool" --version 2>&1 | head -n 1
        ;;
    esac
  else
    echo "missing"
  fi
}

# Lightweight FASTQ sanity check: validate first 10k records structure and read lengths.
fastq_sample_metrics() {
  local fq="$1"
  set +o pipefail
  gunzip -c "$fq" | awk '
    NR > 40000 { exit }
    {
      mod = NR % 4
      if (mod == 1) {
        if (substr($0,1,1) != "@") bad_header++
      } else if (mod == 3) {
        if (substr($0,1,1) != "+") bad_plus++
      } else if (mod == 2) {
        seq_len = length($0)
      } else if (mod == 0) {
        qual_len = length($0)
        if (seq_len != qual_len) bad_len++
        if (min_len == 0 || seq_len < min_len) min_len = seq_len
        if (seq_len > max_len) max_len = seq_len
        sum_len += seq_len
        reads++
      }
    }
    END {
      if (reads == 0) {
        printf "reads=0\tmean_len=NA\tmin_len=NA\tmax_len=NA\tformat_ok=no\n"
      } else {
        fmt = (bad_header == 0 && bad_plus == 0 && bad_len == 0) ? "yes" : "no"
        printf "reads=%d\tmean_len=%.2f\tmin_len=%d\tmax_len=%d\tformat_ok=%s\n", reads, sum_len/reads, min_len, max_len, fmt
      }
    }
  '
  set -o pipefail
}

R1_BYTES="$(stat -f '%z' "$MASON_R1")"
R2_BYTES="$(stat -f '%z' "$MASON_R2")"
R1_SIZE_H="$(ls -lh "$MASON_R1" | awk '{print $5}')"
R2_SIZE_H="$(ls -lh "$MASON_R2" | awk '{print $5}')"

R1_METRICS="$(fastq_sample_metrics "$MASON_R1")"
R2_METRICS="$(fastq_sample_metrics "$MASON_R2")"

R1_READS="$(echo "$R1_METRICS" | awk -F '\t' '{print $1}' | cut -d '=' -f2)"
R1_MEAN_LEN="$(echo "$R1_METRICS" | awk -F '\t' '{print $2}' | cut -d '=' -f2)"
R1_MIN_LEN="$(echo "$R1_METRICS" | awk -F '\t' '{print $3}' | cut -d '=' -f2)"
R1_MAX_LEN="$(echo "$R1_METRICS" | awk -F '\t' '{print $4}' | cut -d '=' -f2)"
R1_FMT_OK="$(echo "$R1_METRICS" | awk -F '\t' '{print $5}' | cut -d '=' -f2)"

R2_READS="$(echo "$R2_METRICS" | awk -F '\t' '{print $1}' | cut -d '=' -f2)"
R2_MEAN_LEN="$(echo "$R2_METRICS" | awk -F '\t' '{print $2}' | cut -d '=' -f2)"
R2_MIN_LEN="$(echo "$R2_METRICS" | awk -F '\t' '{print $3}' | cut -d '=' -f2)"
R2_MAX_LEN="$(echo "$R2_METRICS" | awk -F '\t' '{print $4}' | cut -d '=' -f2)"
R2_FMT_OK="$(echo "$R2_METRICS" | awk -F '\t' '{print $5}' | cut -d '=' -f2)"

FASTQC_PRESENT="$(have_tool fastqc)"
MULTIQC_PRESENT="$(have_tool multiqc)"
FASTP_PRESENT="$(have_tool fastp)"
BWA_PRESENT="$(have_tool bwa-mem2)"
if [[ "$BWA_PRESENT" == "no" ]]; then
  BWA_PRESENT="$(have_tool bwa)"
fi
SAMTOOLS_PRESENT="$(have_tool samtools)"
BCFTOOLS_PRESENT="$(have_tool bcftools)"
SNPEFF_PRESENT="$(have_tool snpEff)"
VEP_PRESENT="$(have_tool vep)"

REF_PRESENT="no"
if find data/refs -type f ! -name '.gitkeep' | grep -q .; then
  REF_PRESENT="yes"
fi

# Gate outcomes
GATE_FASTQ_SANITY="fail"
if [[ "$R1_FMT_OK" == "yes" && "$R2_FMT_OK" == "yes" ]]; then
  GATE_FASTQ_SANITY="pass"
fi

# Full integrity (stream complete) not run in this baseline pass.
GATE_FASTQ_INTEGRITY="incomplete"

GATE_ALIGNMENT_READY="fail"
if [[ "$BWA_PRESENT" == "yes" && "$SAMTOOLS_PRESENT" == "yes" && "$REF_PRESENT" == "yes" ]]; then
  GATE_ALIGNMENT_READY="pass"
fi

GATE_VARIANT_READY="fail"
if [[ "$BCFTOOLS_PRESENT" == "yes" && "$GATE_ALIGNMENT_READY" == "pass" ]]; then
  GATE_VARIANT_READY="pass"
fi

GATE_COUPLE_READY="fail"
if [[ "$HANNAH_STATUS" == "ready" ]]; then
  GATE_COUPLE_READY="pass"
fi

cat > results/qc_metrics.tsv <<TSV
run_timestamp_utc\tsample_id\tread_pair\tfile_path\tmetric\tvalue\tstatus\tnotes
$RUN_TS_UTC\tmason\tR1\t$MASON_R1\tfile_size_bytes\t$R1_BYTES\tmeasured\tcompressed_size=$R1_SIZE_H
$RUN_TS_UTC\tmason\tR2\t$MASON_R2\tfile_size_bytes\t$R2_BYTES\tmeasured\tcompressed_size=$R2_SIZE_H
$RUN_TS_UTC\tmason\tR1\t$MASON_R1\tformat_sanity_sampled_reads\t$R1_READS\tmeasured\tfirst_10k_reads
$RUN_TS_UTC\tmason\tR2\t$MASON_R2\tformat_sanity_sampled_reads\t$R2_READS\tmeasured\tfirst_10k_reads
$RUN_TS_UTC\tmason\tR1\t$MASON_R1\tmean_read_length\t$R1_MEAN_LEN\tmeasured\tfirst_10k_reads
$RUN_TS_UTC\tmason\tR2\t$MASON_R2\tmean_read_length\t$R2_MEAN_LEN\tmeasured\tfirst_10k_reads
$RUN_TS_UTC\tmason\tR1\t$MASON_R1\tmin_read_length\t$R1_MIN_LEN\tmeasured\tfirst_10k_reads
$RUN_TS_UTC\tmason\tR2\t$MASON_R2\tmin_read_length\t$R2_MIN_LEN\tmeasured\tfirst_10k_reads
$RUN_TS_UTC\tmason\tR1\t$MASON_R1\tmax_read_length\t$R1_MAX_LEN\tmeasured\tfirst_10k_reads
$RUN_TS_UTC\tmason\tR2\t$MASON_R2\tmax_read_length\t$R2_MAX_LEN\tmeasured\tfirst_10k_reads
$RUN_TS_UTC\tmason\tR1\t$MASON_R1\tformat_sanity\t$R1_FMT_OK\t$GATE_FASTQ_SANITY\tFASTQ_structure_check_on_first_10k_reads
$RUN_TS_UTC\tmason\tR2\t$MASON_R2\tformat_sanity\t$R2_FMT_OK\t$GATE_FASTQ_SANITY\tFASTQ_structure_check_on_first_10k_reads
$RUN_TS_UTC\tmason\tboth\t$MASON_R1,$MASON_R2\tfastq_integrity_full_stream\tnot_run\t$GATE_FASTQ_INTEGRITY\tfull_gzip_stream_validation_not_executed_in_baseline
TSV

cat > results/alignment_metrics.tsv <<TSV
run_timestamp_utc\tsample_id\tmetric\tvalue\tstatus\tnotes
$RUN_TS_UTC\tmason\talignment_executed\tno\tblocked\trequires_bwa_or_bwa-mem2_plus_samtools_and_GRCh38_reference
$RUN_TS_UTC\tmason\tbwa_available\t$BWA_PRESENT\t$([[ "$BWA_PRESENT" == "yes" ]] && echo pass || echo fail)\ttool_check
$RUN_TS_UTC\tmason\tsamtools_available\t$SAMTOOLS_PRESENT\t$([[ "$SAMTOOLS_PRESENT" == "yes" ]] && echo pass || echo fail)\ttool_check
$RUN_TS_UTC\tmason\treference_bundle_present\t$REF_PRESENT\t$([[ "$REF_PRESENT" == "yes" ]] && echo pass || echo fail)\tdata/refs currently missing GRCh38 fasta+indexes
TSV

cat > results/coverage_metrics.tsv <<TSV
run_timestamp_utc\tsample_id\tmetric\tvalue\tstatus\tnotes
$RUN_TS_UTC\tmason\tcoverage_assessment_executed\tno\tblocked\trequires_aligned_bam_and_samtools_depth_or_equivalent
$RUN_TS_UTC\tmason\tcoverage_gate\tnot_evaluable\tfail\talignment_not_available
TSV

cat > results/annotated_variants.tsv <<'TSV'
sample_id	chrom	pos	ref	alt	gene	transcript	variant_hgvs	zygosity	genotype	depth	genotype_quality	allele_balance	clinvar_significance	clinvar_review_status	gnomad_af	confidence_tier	inclusion_reason	confirmation_recommendation	analysis_status	notes
TSV

cat > results/carrier_candidates_couple.tsv <<'TSV'
pair_id	sample_id	partner_id	gene	transcript	variant_hgvs	chrom	pos	ref	alt	zygosity	genotype	depth	genotype_quality	allele_balance	clinvar_significance	clinvar_review_status	gnomad_af	inheritance_model	couple_risk_flag	confidence_tier	inclusion_reason	confirmation_recommendation
TSV

cat > reports/limitations_and_followup.md <<MD
# Limitations and Follow-up

## Technical blind spots (must disclose)
- Repeat expansions are not assessed by this baseline workflow.
- Complex CNVs/gene conversions are not assessed.
- Difficult paralog/pseudogene loci are at elevated mapping risk.
- Structural variants beyond short-read small-variant callers are not assessed.
- Low/uneven coverage loci cannot be evaluated until alignment+depth are available.

## Current run-specific blockers
- Reference bundle in data/refs/ is missing (GRCh38 fasta/indexes and annotation resources).
- Annotation databases are not configured yet (snpEff genomes and VEP cache/FASTA are still required).
- VEP command-line tool availability: $VEP_PRESENT ($(tool_version vep)).
- snpEff command-line tool availability: $SNPEFF_PRESENT ($(tool_version snpEff)).
- Partner sample hannah remains pending in config/samples.tsv and data/raw/hannah/.

## Follow-up actions (minimum to unblock)
1. Add pinned GRCh38 reference assets to data/refs/ and document versions in run summary.
2. Configure annotation datasets (snpEff database and VEP cache/FASTA if VEP is used).
3. Import Hannah FASTQs and update config/samples.tsv to ready with real filenames.
4. Re-run full workflow including full FASTQ integrity scan, alignment, variant calling, annotation, and couple-level filtering.

## Confirmation guidance policy
- Any future actionable candidate must be tagged: research-only pending clinical confirmation.
- Confirmation should use an orthogonal clinical assay (e.g., validated targeted sequencing, CNV assay when relevant, and genetics professional review).
MD

cat > reports/run_summary.md <<MD
# Run Summary

## Run metadata
- Run timestamp (UTC): $RUN_TS_UTC
- Run timestamp (local): $RUN_TS_LOCAL
- Mode: Mason-only baseline execution (single-sample technical pass)

## Input manifest
- Mason R1 path: $MASON_R1
- Mason R1 compressed size: $R1_SIZE_H ($R1_BYTES bytes)
- Mason R2 path: $MASON_R2
- Mason R2 compressed size: $R2_SIZE_H ($R2_BYTES bytes)
- Hannah sample status from config/samples.tsv: $HANNAH_STATUS

## Tool and reference availability
| Component | Availability | Version/Note |
|---|---|---|
| fastqc | $FASTQC_PRESENT | $(tool_version fastqc) |
| multiqc | $MULTIQC_PRESENT | $(tool_version multiqc) |
| fastp | $FASTP_PRESENT | $(tool_version fastp) |
| bwa/bwa-mem2 | $BWA_PRESENT | $(tool_version bwa-mem2) / $(tool_version bwa) |
| samtools | $SAMTOOLS_PRESENT | $(tool_version samtools) |
| bcftools | $BCFTOOLS_PRESENT | $(tool_version bcftools) |
| snpEff | $SNPEFF_PRESENT | $(tool_version snpEff) |
| vep | $VEP_PRESENT | $(tool_version vep) |
| GRCh38 reference bundle (data/refs/) | $REF_PRESENT | requires fasta + indexes + annotation DBs |

## QC and analysis gate outcomes
| Gate | Outcome | Evidence |
|---|---|---|
| FASTQ format sanity (sampled first 10k reads/file) | $GATE_FASTQ_SANITY | R1 format_ok=$R1_FMT_OK, R2 format_ok=$R2_FMT_OK |
| FASTQ full-stream gzip integrity | $GATE_FASTQ_INTEGRITY | full stream check not executed in this baseline |
| Alignment gate (tools + reference) | $GATE_ALIGNMENT_READY | bwa_available=$BWA_PRESENT, samtools_available=$SAMTOOLS_PRESENT, ref_present=$REF_PRESENT |
| Variant calling gate | $GATE_VARIANT_READY | bcftools_available=$BCFTOOLS_PRESENT and alignment gate not satisfied |
| Couple-level interpretation gate | $GATE_COUPLE_READY | Hannah sample status is $HANNAH_STATUS |

Overall status: **analysis incomplete**

## Candidate counts by tier
- Tier A: 0 (not evaluated; variant calling not run)
- Tier B: 0 (not evaluated; variant calling not run)
- Tier C: 0 (not evaluated; variant calling not run)
- Not reportable: 0 (not evaluated; variant calling not run)

## Top candidate table
No candidates available because variant calling/annotation was blocked.

## Deliverables generated in this run
- results/qc_metrics.tsv
- results/alignment_metrics.tsv
- results/coverage_metrics.tsv
- results/annotated_variants.tsv (header only)
- results/carrier_candidates_couple.tsv (header only)
- reports/limitations_and_followup.md
- reports/run_summary.md

## Limitations and unresolved blockers
- Reference bundle and annotation datasets are not configured, so full alignment/annotation remains blocked.
- Couple-level risk logic cannot be completed without Hannah FASTQs.
- Full FASTQ stream integrity validation remains pending.

## Clinical confirmation disclaimer
Any actionable future finding must be considered **research-only pending clinical confirmation**.
No diagnostic claim is made from this baseline run.
MD

cat > logs/tool_availability.tsv <<TSV
run_timestamp_utc\ttool\tavailable\tversion
$RUN_TS_UTC\tfastqc\t$FASTQC_PRESENT\t$(tool_version fastqc)
$RUN_TS_UTC\tmultiqc\t$MULTIQC_PRESENT\t$(tool_version multiqc)
$RUN_TS_UTC\tfastp\t$FASTP_PRESENT\t$(tool_version fastp)
$RUN_TS_UTC\tbwa-mem2\t$(have_tool bwa-mem2)\t$(tool_version bwa-mem2)
$RUN_TS_UTC\tbwa\t$(have_tool bwa)\t$(tool_version bwa)
$RUN_TS_UTC\tsamtools\t$SAMTOOLS_PRESENT\t$(tool_version samtools)
$RUN_TS_UTC\tbcftools\t$BCFTOOLS_PRESENT\t$(tool_version bcftools)
$RUN_TS_UTC\tsnpEff\t$SNPEFF_PRESENT\t$(tool_version snpEff)
$RUN_TS_UTC\tvep\t$VEP_PRESENT\t$(tool_version vep)
TSV

echo "Baseline deliverables generated at $RUN_TS_UTC"
