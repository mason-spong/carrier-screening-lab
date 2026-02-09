# Output Specification

## `results/carrier_candidates_couple.tsv`

Required columns:

- `pair_id`
- `sample_id`
- `partner_id`
- `gene`
- `transcript`
- `variant_hgvs` (or canonical variant string used consistently)
- `chrom`
- `pos`
- `ref`
- `alt`
- `zygosity`
- `genotype`
- `depth`
- `genotype_quality`
- `allele_balance`
- `clinvar_significance`
- `clinvar_review_status`
- `gnomad_af`
- `inheritance_model`
- `couple_risk_flag`
- `confidence_tier`
- `inclusion_reason`
- `confirmation_recommendation`

## `reports/run_summary.md`

Must include:

1. Run timestamp and tool versions
2. Input file manifest
3. QC/alignment/variant-calling gate outcomes
4. Candidate counts by tier
5. Top candidate table (brief)
6. Limitations and unresolved blockers
7. Clinical confirmation disclaimer

## `reports/limitations_and_followup.md`

Must include:

- Variant classes not well captured by this pipeline
- Genes/loci with known technical complexity
- Recommended follow-up test types for high-impact unresolved risk

