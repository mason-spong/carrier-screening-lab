# Quality Standards

## Standard of Evidence

Prioritize variants with converging support:

- Strong technical quality (depth, genotype quality, allele balance)
- Credible clinical annotation (ClinVar review strength, not just assertion count)
- Population plausibility (gnomAD frequency consistent with severe Mendelian disease model)
- Inheritance consistency at couple level

## Mandatory QC Gates (Before Interpretation)

1. FASTQ integrity and format sanity pass.
2. Alignment metrics within expected range (mapping rate, duplicate burden).
3. Coverage acceptable for target interpretation scope.
4. Contamination and sample-identity checks are acceptable.
5. Joint-call VCF passes defined filters.

If any gate fails, report as `analysis incomplete` rather than over-interpreting.

## Candidate Tiering

- `Tier A - High confidence candidate`
  - Good technical metrics
  - Plausible consequence/inheritance
  - Strong or moderately strong clinical evidence
- `Tier B - Needs orthogonal confirmation`
  - Potentially relevant but evidence conflicts, marginal quality, or uncertain interpretation
- `Tier C - Low confidence / likely artifact`
  - Weak technical support, problematic locus, or low clinical relevance
- `Not reportable`
  - Benign/likely benign, common polymorphism, or outside scope

## Coupled Carrier Logic

- Autosomal recessive focus: flag when both partners carry qualifying variants in same gene.
- X-linked focus: flag relevant maternal carrier candidates.
- Always include caveat when phasing or variant class coverage is incomplete.

## Known Blind Spots (Must Be Stated)

- Repeat expansions
- Complex CNVs and gene conversions
- Difficult paralog/pseudogene regions
- Structural variants beyond caller design
- Regions with low or uneven coverage

## Reporting Rules

- Every reported candidate must include:
  - Gene
  - Variant representation
  - Zygosity
  - Technical quality fields used
  - Annotation summary
  - Confidence tier
  - Why it was included
  - Confirmation recommendation

