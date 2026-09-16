# Public Benchmark v0.1 Freeze Record

This record freezes the identity of the public benchmark cohort before
public-cohort ranking performance is evaluated.

## Source

- Source: Monarch Initiative Phenopacket Store
- Release: 0.1.27
- Strict eligible source cases: 8,610

## Frozen cohort

- Total cases: 500
- Development cases: 300
- Locked evaluation cases: 200
- Unique causal genes: 500
- Unique publications: 500
- Development/locked publication overlap: 0
- Development/locked exact causal variant overlap: 0
- Development/locked causal gene overlap: 0

The v0.1 design therefore evaluates cross-publication and cross-gene
generalization. All 200 locked-evaluation causal genes are unseen in
the development cohort.

## Canonical manifest fingerprint

The truth-bearing cohort manifest is not committed to Git.

Canonical serialization uses sorted JSON keys, compact separators,
UTF-8 encoding, and a final newline.

SHA-256:

a8c51d43b0300e0ea12de95e92388cbb8a4cc493b2255e0e7ac833ee4a284d29

Any cohort manifest with a different digest is not the frozen public
benchmark v0.1 cohort.

## Semantic boundaries

The cohort manifest contains benchmark truth and MUST NOT be used as
ranking input.

Benchmark truth MUST NOT enter CandidateCase construction, evidence
integration, interpretation, candidate scoring, rank aggregation, or
any other pre-ranking computation.

Truth may be joined only after ranking for benchmark evaluation.

This freeze record does not constitute clinical validation.
