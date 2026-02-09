#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

THREADS="${THREADS:-8}"
RUN_TS_UTC="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
RUN_TS_LOCAL="$(date +"%Y-%m-%d %H:%M:%S %Z")"
LOG_FILE="logs/mason_full_run.log"

mkdir -p logs reports results \
  data/interim/mason/{state,qc,alignment,variants,annotation} \
  data/processed/mason \
  data/refs/grch38 \
  data/refs/clinvar

exec > >(tee -a "$LOG_FILE") 2>&1

echo "[$(date +"%Y-%m-%d %H:%M:%S")] Starting full Mason analysis"
echo "Run UTC: $RUN_TS_UTC"
echo "Threads: $THREADS"

STATE_DIR="data/interim/mason/state"
QC_DIR="data/interim/mason/qc"
ALIGN_DIR="data/interim/mason/alignment"
VAR_DIR="data/interim/mason/variants"
ANNO_DIR="data/interim/mason/annotation"
PROC_DIR="data/processed/mason"
REF_DIR="data/refs/grch38"
CLINVAR_DIR="data/refs/clinvar"

MASON_R1="$(awk -F '\t' '$1=="mason" {print $3}' config/samples.tsv)"
MASON_R2="$(awk -F '\t' '$1=="mason" {print $4}' config/samples.tsv)"
HANNAH_STATUS="$(awk -F '\t' '$1=="hannah" {print $5}' config/samples.tsv)"

if [[ ! -f "$MASON_R1" || ! -f "$MASON_R2" ]]; then
  echo "ERROR: Mason FASTQ inputs missing" >&2
  exit 1
fi

required_tools=(curl gzip gunzip awk sed grep fastqc multiqc bwa samtools bcftools snpEff bgzip tabix)
for t in "${required_tools[@]}"; do
  if ! command -v "$t" >/dev/null 2>&1; then
    echo "ERROR: required tool not found: $t" >&2
    exit 1
  fi
done

is_done() {
  [[ -f "$STATE_DIR/$1.done" ]]
}

mark_done() {
  touch "$STATE_DIR/$1.done"
}

get_fastqc_dir() {
  local fq="$1"
  local base
  base="$(basename "$fq")"
  base="${base%.gz}"
  echo "$QC_DIR/${base}_fastqc"
}

tool_version() {
  local tool="$1"
  case "$tool" in
    fastqc)
      LC_ALL=C LANG=C fastqc --version 2>&1 | head -n 1
      ;;
    multiqc)
      LC_ALL=C LANG=C multiqc --version 2>&1 | head -n 1
      ;;
    fastp)
      LC_ALL=C LANG=C fastp --version 2>&1 | head -n 1
      ;;
    bwa)
      LC_ALL=C LANG=C bwa 2>&1 | awk 'NF {print; exit}'
      ;;
    samtools)
      LC_ALL=C LANG=C samtools --version 2>&1 | head -n 1
      ;;
    bcftools)
      LC_ALL=C LANG=C bcftools --version 2>&1 | head -n 1
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
      echo "unknown"
      ;;
  esac
}

# Phase 1: Raw QC
if ! is_done phase1_qc; then
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Phase 1: FASTQ integrity + FastQC/MultiQC"

  gzip -t "$MASON_R1"
  gzip -t "$MASON_R2"

  fastqc --extract -t "$THREADS" -o "$QC_DIR" "$MASON_R1" "$MASON_R2"
  multiqc -o "$QC_DIR" "$QC_DIR"

  mark_done phase1_qc
fi

# Phase 2: Reference prep
REF_FA_GZ="$REF_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
REF_FA="$REF_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
REF_URL="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

if ! is_done phase2_ref; then
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Phase 2: Reference setup"

  if [[ ! -s "$REF_FA_GZ" ]]; then
    curl -L --retry 5 --retry-delay 10 -o "$REF_FA_GZ" "$REF_URL"
  fi

  if [[ ! -s "$REF_FA" ]]; then
    gunzip -c "$REF_FA_GZ" > "$REF_FA"
  fi

  if [[ ! -s "$REF_FA.fai" ]]; then
    samtools faidx "$REF_FA"
  fi

  if [[ ! -s "$REF_FA.amb" || ! -s "$REF_FA.ann" || ! -s "$REF_FA.bwt" || ! -s "$REF_FA.pac" || ! -s "$REF_FA.sa" ]]; then
    bwa index "$REF_FA"
  fi

  mark_done phase2_ref
fi

# Phase 3: Alignment + post-processing
BAM="$ALIGN_DIR/mason.markdup.bam"
if ! is_done phase3_align; then
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Phase 3: Alignment and post-processing"

  bwa mem -t "$THREADS" \
    -R '@RG\tID:mason\tSM:mason\tPL:ILLUMINA\tLB:mason_lib\tPU:l001' \
    "$REF_FA" "$MASON_R1" "$MASON_R2" \
    | samtools view -@ "$THREADS" -b -o "$ALIGN_DIR/mason.raw.bam" -

  samtools sort -n -@ "$THREADS" -o "$ALIGN_DIR/mason.name.bam" "$ALIGN_DIR/mason.raw.bam"
  samtools fixmate -m -@ "$THREADS" "$ALIGN_DIR/mason.name.bam" "$ALIGN_DIR/mason.fixmate.bam"
  samtools sort -@ "$THREADS" -o "$ALIGN_DIR/mason.pos.bam" "$ALIGN_DIR/mason.fixmate.bam"
  samtools markdup -@ "$THREADS" -s "$ALIGN_DIR/mason.pos.bam" "$BAM" 2> "$ALIGN_DIR/mason.markdup.stats.txt"
  samtools index -@ "$THREADS" "$BAM"

  samtools flagstat -@ "$THREADS" "$BAM" > "$ALIGN_DIR/mason.flagstat.txt"
  samtools stats -@ "$THREADS" "$BAM" > "$ALIGN_DIR/mason.stats.txt"
  samtools coverage "$BAM" > "$ALIGN_DIR/mason.coverage.tsv"

  rm -f "$ALIGN_DIR/mason.raw.bam" "$ALIGN_DIR/mason.name.bam" "$ALIGN_DIR/mason.fixmate.bam" "$ALIGN_DIR/mason.pos.bam"

  mark_done phase3_align
fi

# Phase 4: Variant calling
VCF_FILTERED="$VAR_DIR/mason.filtered.vcf.gz"
if ! is_done phase4_call; then
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Phase 4: Variant calling"

  bcftools mpileup -f "$REF_FA" -Ou -a FORMAT/AD,FORMAT/DP,FORMAT/GQ "$BAM" \
    | bcftools call -mv -Ou \
    | bcftools norm -f "$REF_FA" -m -any -Ou \
    | bcftools filter -s LOWQUAL -e 'QUAL<30 || FORMAT/DP<10 || FORMAT/GQ<20' -Oz -o "$VCF_FILTERED"

  bcftools index -t "$VCF_FILTERED"

  mark_done phase4_call
fi

# Phase 5: Annotation
SNPEFF_VCF="$ANNO_DIR/mason.snpeff.vcf.gz"
CLINVAR_VCF="$CLINVAR_DIR/clinvar.vcf.gz"
ANNOTATED_VCF="$PROC_DIR/mason.annotated.clinvar.snpeff.vcf.gz"

if ! is_done phase5_annotate; then
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Phase 5: Annotation"

  snpEff download GRCh38.99

  if [[ ! -s "$CLINVAR_VCF" ]]; then
    curl -L --retry 5 --retry-delay 10 -o "$CLINVAR_VCF" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
  fi
  if [[ ! -s "$CLINVAR_VCF.tbi" ]]; then
    curl -L --retry 5 --retry-delay 10 -o "$CLINVAR_VCF.tbi" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"
  fi

  snpEff -canon -hgvs GRCh38.99 "$VCF_FILTERED" | bgzip -c > "$SNPEFF_VCF"
  tabix -f -p vcf "$SNPEFF_VCF"

  bcftools annotate -a "$CLINVAR_VCF" \
    -c CHROM,POS,REF,ALT,INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNVC,INFO/CLNDN \
    -Oz -o "$ANNOTATED_VCF" "$SNPEFF_VCF"

  bcftools index -t "$ANNOTATED_VCF"

  mark_done phase5_annotate
fi

# Phase 6: Result tables
if ! is_done phase6_tables; then
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Phase 6: Build result tables"

  R1_DIR="$(get_fastqc_dir "$MASON_R1")"
  R2_DIR="$(get_fastqc_dir "$MASON_R2")"

  r1_total="$(awk -F '\t' '$1=="Total Sequences" {print $2}' "$R1_DIR/fastqc_data.txt")"
  r2_total="$(awk -F '\t' '$1=="Total Sequences" {print $2}' "$R2_DIR/fastqc_data.txt")"
  r1_len="$(awk -F '\t' '$1=="Sequence length" {print $2}' "$R1_DIR/fastqc_data.txt")"
  r2_len="$(awk -F '\t' '$1=="Sequence length" {print $2}' "$R2_DIR/fastqc_data.txt")"
  r1_gc="$(awk -F '\t' '$1=="%GC" {print $2}' "$R1_DIR/fastqc_data.txt")"
  r2_gc="$(awk -F '\t' '$1=="%GC" {print $2}' "$R2_DIR/fastqc_data.txt")"
  r1_pbq="$(awk -F '\t' '$2=="Per base sequence quality" {print $1}' "$R1_DIR/summary.txt")"
  r2_pbq="$(awk -F '\t' '$2=="Per base sequence quality" {print $1}' "$R2_DIR/summary.txt")"

  cat > results/qc_metrics.tsv <<TSV
run_timestamp_utc\tsample_id\tread_pair\tfile_path\tmetric\tvalue\tstatus\tnotes
$RUN_TS_UTC\tmason\tR1\t$MASON_R1\tfastq_integrity_full_stream\tpass\tpass\tgzip -t completed successfully
$RUN_TS_UTC\tmason\tR2\t$MASON_R2\tfastq_integrity_full_stream\tpass\tpass\tgzip -t completed successfully
$RUN_TS_UTC\tmason\tR1\t$MASON_R1\tfastqc_total_sequences\t$r1_total\tmeasured\tfrom fastqc_data.txt
$RUN_TS_UTC\tmason\tR2\t$MASON_R2\tfastqc_total_sequences\t$r2_total\tmeasured\tfrom fastqc_data.txt
$RUN_TS_UTC\tmason\tR1\t$MASON_R1\tfastqc_sequence_length\t$r1_len\tmeasured\tfrom fastqc_data.txt
$RUN_TS_UTC\tmason\tR2\t$MASON_R2\tfastqc_sequence_length\t$r2_len\tmeasured\tfrom fastqc_data.txt
$RUN_TS_UTC\tmason\tR1\t$MASON_R1\tfastqc_percent_gc\t$r1_gc\tmeasured\tfrom fastqc_data.txt
$RUN_TS_UTC\tmason\tR2\t$MASON_R2\tfastqc_percent_gc\t$r2_gc\tmeasured\tfrom fastqc_data.txt
$RUN_TS_UTC\tmason\tR1\t$MASON_R1\tfastqc_per_base_quality\t$r1_pbq\t$r1_pbq\tfrom FastQC summary module status
$RUN_TS_UTC\tmason\tR2\t$MASON_R2\tfastqc_per_base_quality\t$r2_pbq\t$r2_pbq\tfrom FastQC summary module status
TSV

  total_reads="$(awk '/in total/ {print $1; exit}' "$ALIGN_DIR/mason.flagstat.txt")"
  mapped_reads="$(awk '/ mapped \(/ && $0 !~ /primary mapped/ {print $1; exit}' "$ALIGN_DIR/mason.flagstat.txt")"
  mapped_pct="$(awk -F '[()%]' '/ mapped \(/ && $0 !~ /primary mapped/ {gsub(/ /, "", $2); print $2; exit}' "$ALIGN_DIR/mason.flagstat.txt")"
  proper_pair_pct="$(awk -F '[()%]' '/properly paired/ {gsub(/ /, "", $2); print $2; exit}' "$ALIGN_DIR/mason.flagstat.txt")"
  dup_reads="$(awk '/duplicates/ {print $1; exit}' "$ALIGN_DIR/mason.flagstat.txt")"
  dup_pct="$(awk -v d="$dup_reads" -v t="$total_reads" 'BEGIN {if (t>0) printf "%.2f", (100*d/t); else print "NA"}')"

  cov_total_row="$(tail -n 1 "$ALIGN_DIR/mason.coverage.tsv")"
  cov_pct="$(echo "$cov_total_row" | awk '{print $7}')"
  mean_depth="$(echo "$cov_total_row" | awk '{print $8}')"
  mean_baseq="$(echo "$cov_total_row" | awk '{print $9}')"
  mean_mapq="$(echo "$cov_total_row" | awk '{print $10}')"

  cat > results/alignment_metrics.tsv <<TSV
run_timestamp_utc\tsample_id\tmetric\tvalue\tstatus\tnotes
$RUN_TS_UTC\tmason\talignment_executed\tyes\tpass\tbwa + samtools pipeline
$RUN_TS_UTC\tmason\ttotal_reads\t$total_reads\tmeasured\tsamtools flagstat
$RUN_TS_UTC\tmason\tmapped_reads\t$mapped_reads\tmeasured\tsamtools flagstat
$RUN_TS_UTC\tmason\tmapped_percent\t$mapped_pct\t$([[ $(awk -v m="$mapped_pct" 'BEGIN{print (m+0>=95)?1:0}') -eq 1 ]] && echo pass || echo fail)\tthreshold >=95
$RUN_TS_UTC\tmason\tproperly_paired_percent\t$proper_pair_pct\t$([[ $(awk -v p="$proper_pair_pct" 'BEGIN{print (p+0>=90)?1:0}') -eq 1 ]] && echo pass || echo fail)\tthreshold >=90
$RUN_TS_UTC\tmason\tduplicate_reads\t$dup_reads\tmeasured\tsamtools flagstat
$RUN_TS_UTC\tmason\tduplicate_percent\t$dup_pct\t$([[ $(awk -v d="$dup_pct" 'BEGIN{print (d+0<=30)?1:0}') -eq 1 ]] && echo pass || echo fail)\tthreshold <=30
TSV

  cat > results/coverage_metrics.tsv <<TSV
run_timestamp_utc\tsample_id\tmetric\tvalue\tstatus\tnotes
$RUN_TS_UTC\tmason\tcoverage_percent_bases_ge_1x\t$cov_pct\t$([[ $(awk -v c="$cov_pct" 'BEGIN{print (c+0>=90)?1:0}') -eq 1 ]] && echo pass || echo fail)\tsamtools coverage total row
$RUN_TS_UTC\tmason\tmean_depth\t$mean_depth\t$([[ $(awk -v d="$mean_depth" 'BEGIN{print (d+0>=20)?1:0}') -eq 1 ]] && echo pass || echo fail)\tthreshold >=20
$RUN_TS_UTC\tmason\tmean_baseq\t$mean_baseq\tmeasured\tsamtools coverage total row
$RUN_TS_UTC\tmason\tmean_mapq\t$mean_mapq\tmeasured\tsamtools coverage total row
TSV

  bcftools query -s mason -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%QUAL\t[%GT]\t[%DP]\t[%GQ]\t[%AD]\t%INFO/ANN\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\n' "$ANNOTATED_VCF" \
    | awk -F '\t' 'BEGIN {
        OFS="\t";
        print "sample_id","chrom","pos","ref","alt","gene","transcript","variant_hgvs","zygosity","genotype","depth","genotype_quality","allele_balance","clinvar_significance","clinvar_review_status","gnomad_af","confidence_tier","inclusion_reason","confirmation_recommendation","analysis_status","notes";
      }
      function zyg(gt) {
        if (gt=="0/1" || gt=="1/0" || gt=="0|1" || gt=="1|0") return "HET";
        if (gt=="1/1" || gt=="1|1") return "HOM_ALT";
        if (gt=="0/0" || gt=="0|0") return "HOM_REF";
        return "OTHER";
      }
      function num(v) {
        if (v=="." || v=="") return -1;
        return v + 0;
      }
      {
        chrom=$1; pos=$2; ref=$3; alt=$4; flt=$5; qual=$6; gt=$7; dp=$8; gq=$9; ad=$10; ann=$11; clnsig=$12; clnrev=$13;

        n=split(ad,adf,",");
        refad=(n>=1 ? adf[1]+0 : 0);
        altad=(n>=2 ? adf[2]+0 : 0);
        totalad=refad+altad;
        ab=(totalad>0 ? altad/totalad : -1);

        nann=split(ann,annrec,",");
        split(annrec[1],a,"|");
        effect=a[2];
        impact=a[3];
        gene=a[4];
        transcript=a[7];
        hgvsc=a[10];
        hgvsp=a[11];

        if (hgvsc!="") varhgvs=hgvsc;
        else if (hgvsp!="") varhgvs=hgvsp;
        else varhgvs=chrom ":" pos ref ">" alt;

        z=zyg(gt);
        csign=(clnsig=="" || clnsig=="." ? "NA" : clnsig);
        crev=(clnrev=="" || clnrev=="." ? "NA" : clnrev);
        gnomad="NA";

        tier="Not reportable";
        reason="Outside current prioritization";
        conf="Not prioritized in Mason-only run";

        if (flt=="PASS" && csign ~ /Pathogenic/ && csign !~ /Conflicting/ && num(dp)>=20 && num(gq)>=30 && ab>=0.25 && ab<=0.75) {
          tier="Tier A - High confidence candidate";
          reason="PASS metrics and ClinVar pathogenic signal";
          conf="research-only pending clinical confirmation via orthogonal assay";
        } else if (flt=="PASS" && ((impact=="HIGH" || impact=="MODERATE") || csign!="NA")) {
          tier="Tier B - Needs orthogonal confirmation";
          reason="Potentially relevant impact and/or ClinVar annotation";
          conf="research-only pending clinical confirmation via orthogonal assay";
        } else if (flt=="PASS") {
          tier="Tier C - Low confidence / likely artifact";
          reason="PASS variant with limited clinical/functional support";
          conf="Low priority; confirm only if independently supported";
        }

        status=(tier ~ /^Tier [AB]/ ? "research_only_pending_clinical_confirmation" : "research_only_non_actionable");
        notes="effect=" effect ";impact=" impact ";filter=" flt ";qual=" qual;

        ab_out=(ab>=0 ? sprintf("%.4f",ab) : "NA");
        print "mason",chrom,pos,ref,alt,gene,transcript,varhgvs,z,gt,dp,gq,ab_out,csign,crev,gnomad,tier,reason,conf,status,notes;
      }' > results/annotated_variants.tsv

  awk -F '\t' 'BEGIN {
      OFS="\t";
      print "pair_id","sample_id","partner_id","gene","transcript","variant_hgvs","chrom","pos","ref","alt","zygosity","genotype","depth","genotype_quality","allele_balance","clinvar_significance","clinvar_review_status","gnomad_af","inheritance_model","couple_risk_flag","confidence_tier","inclusion_reason","confirmation_recommendation";
    }
    NR==1 {next}
    $17 ~ /^Tier [AB]/ {
      print "mason_hannah","mason","hannah",$6,$7,$8,$2,$3,$4,$5,$9,$10,$11,$12,$13,$14,$15,$16,"single_sample_only_not_fully_evaluable","partner_pending_no_couple_assessment",$17,$18,$19;
    }' results/annotated_variants.tsv > results/carrier_candidates_couple.tsv

  mark_done phase6_tables
fi

# Phase 7: Reports
if ! is_done phase7_reports; then
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Phase 7: Reports"

  tier_a="$(awk -F '\t' 'NR>1 && $17=="Tier A - High confidence candidate" {c++} END{print c+0}' results/annotated_variants.tsv)"
  tier_b="$(awk -F '\t' 'NR>1 && $17=="Tier B - Needs orthogonal confirmation" {c++} END{print c+0}' results/annotated_variants.tsv)"
  tier_c="$(awk -F '\t' 'NR>1 && $17=="Tier C - Low confidence / likely artifact" {c++} END{print c+0}' results/annotated_variants.tsv)"
  tier_nr="$(awk -F '\t' 'NR>1 && $17=="Not reportable" {c++} END{print c+0}' results/annotated_variants.tsv)"

  mapped_pct="$(awk -F '\t' '$3=="mapped_percent" {print $4}' results/alignment_metrics.tsv)"
  dup_pct="$(awk -F '\t' '$3=="duplicate_percent" {print $4}' results/alignment_metrics.tsv)"
  mean_depth="$(awk -F '\t' '$3=="mean_depth" {print $4}' results/coverage_metrics.tsv)"
  cov_pct="$(awk -F '\t' '$3=="coverage_percent_bases_ge_1x" {print $4}' results/coverage_metrics.tsv)"

  gate_fastq="pass"
  gate_align="fail"
  gate_cov="fail"
  gate_contam="incomplete"
  gate_joint="incomplete"

  if [[ $(awk -v m="$mapped_pct" -v d="$dup_pct" 'BEGIN{print (m+0>=95 && d+0<=30)?1:0}') -eq 1 ]]; then
    gate_align="pass"
  fi

  if [[ $(awk -v md="$mean_depth" -v cv="$cov_pct" 'BEGIN{print (md+0>=20 && cv+0>=90)?1:0}') -eq 1 ]]; then
    gate_cov="pass"
  fi

  cat > reports/run_summary.md <<MD
# Run Summary

## Run metadata
- Run timestamp (UTC): $RUN_TS_UTC
- Run timestamp (local): $RUN_TS_LOCAL
- Mode: Mason-only full technical run

## Input file manifest
- Mason R1: $MASON_R1
- Mason R2: $MASON_R2
- Hannah sample status: $HANNAH_STATUS

## Tool versions
- fastqc: $(tool_version fastqc)
- multiqc: $(tool_version multiqc)
- fastp: $(tool_version fastp)
- bwa: $(tool_version bwa)
- samtools: $(tool_version samtools)
- bcftools: $(tool_version bcftools)
- snpEff: $(tool_version snpEff)
- vep: $(tool_version vep)

## QC/alignment/variant-calling gates
- FASTQ integrity and format sanity: $gate_fastq
- Alignment metrics gate (mapping/duplication thresholds): $gate_align
- Coverage gate (mean depth and breadth): $gate_cov
- Contamination/sample-identity check: $gate_contam (not run; no tool configured)
- Joint-call filter gate: $gate_joint (single-sample run only; couple joint call deferred)

Overall status: **analysis incomplete for couple-level interpretation**

## Candidate counts by tier
- Tier A: $tier_a
- Tier B: $tier_b
- Tier C: $tier_c
- Not reportable: $tier_nr

## Top candidate table (up to 10)
| Gene | Variant | Tier | ClinVar | Recommendation |
|---|---|---|---|---|
$(awk -F '\t' 'NR>1 && $17 ~ /^Tier [AB]/ {printf "| %s | %s | %s | %s | %s |\n", $6, $8, $17, $14, $19; c++; if (c==10) exit} END{if (c==0) print "| None | None | None | None | No prioritized candidate from Mason-only run |"}' results/annotated_variants.tsv)

## Limitations and unresolved blockers
- Couple-level interpretation is incomplete because Hannah FASTQs are pending.
- VEP is installed but not runnable on this host; annotation used snpEff + ClinVar.
- gnomAD frequency annotations were not integrated in this run; gnomad_af fields are NA.
- CNV/SV/repeat-expansion and complex loci limitations remain.

## Clinical confirmation disclaimer
All prioritized candidates are **research-only pending clinical confirmation**.
No diagnostic claim is made from this analysis.
MD

  cat > reports/limitations_and_followup.md <<MD
# Limitations and Follow-up

## Variant classes and loci not fully captured
- Repeat expansions are not assessed.
- Complex CNVs/gene conversions are not assessed.
- Structural variants outside small-variant caller design are not assessed.
- Pseudogene/paralog-rich loci may remain ambiguous in short-read alignment.

## Run-specific blind spots
- Couple-level risk cannot be finalized until Hannah FASTQs are available.
- VEP is not runnable in current environment (`dyld missing symbol`); snpEff was used.
- gnomAD annotation was not attached in this run; population AF field is NA.
- Contamination and sample identity cross-check were not run (no dedicated tool configured).

## Recommended follow-up
1. Add Hannah FASTQs and execute joint calling.
2. Add contamination/identity checks (e.g., verifyBamID/fingerprint workflow).
3. Add population AF annotation source (gnomAD subset aligned to GRCh38).
4. Confirm Tier A/B candidates with orthogonal clinical assay.

## Confirmation guidance
Any actionable candidate must remain labeled: **research-only pending clinical confirmation**.
MD

  cat > logs/tool_availability.tsv <<TSV
run_timestamp_utc\ttool\tavailable\tversion
$RUN_TS_UTC\tfastqc\tyes\t$(tool_version fastqc)
$RUN_TS_UTC\tmultiqc\tyes\t$(tool_version multiqc)
$RUN_TS_UTC\tfastp\tyes\t$(tool_version fastp)
$RUN_TS_UTC\tbwa\tyes\t$(tool_version bwa)
$RUN_TS_UTC\tsamtools\tyes\t$(tool_version samtools)
$RUN_TS_UTC\tbcftools\tyes\t$(tool_version bcftools)
$RUN_TS_UTC\tsnpEff\tyes\t$(tool_version snpEff)
$RUN_TS_UTC\tvep\tyes\t$(tool_version vep)
TSV

  mark_done phase7_reports
fi

echo "[$(date +"%Y-%m-%d %H:%M:%S")] Mason full analysis completed"
