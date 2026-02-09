# Project Context

## Objective

Generate a **couple-focused carrier-screening candidate list** from short-read FASTQ data for:

- `mason` (available)
- `hannah` (pending FASTQ import)

This is a research workflow intended to prioritize candidates for clinical follow-up.

## Current Data State

- Mason FASTQ files are present in `data/raw/mason/`.
- Hannah placeholders exist in `config/samples.tsv` and must be replaced with real file names once downloaded.

## Analysis Philosophy

- Use conservative filtering and transparent evidence.
- Minimize false confidence.
- Make technical limitations explicit.

## Out of Scope

- Diagnostic conclusions
- Treatment recommendations
- Returning “negative means no risk” statements

## Success Criteria

- End-to-end reproducible run from raw FASTQ to candidate list.
- Clear confidence tiers with rationale.
- Ready-to-hand-off package for genetics professional review.

