# Public Case Benchmark Protocol v0.1

## Purpose

This benchmark evaluates relative candidate prioritization performance
on public, solved, phenotype-annotated rare-disease cases.

It is an analysis benchmark. It is not clinical validation and does not
measure diagnostic accuracy, pathogenicity probability, or ACMG/AMP
classification performance.

## Source

Primary source:

- Monarch Initiative Phenopacket Store
- Frozen source release: 0.1.27

The original source phenopackets may contain solved diagnoses and causal
variant interpretations. These truth-bearing fields MUST NOT enter the
ranking-time case representation.

## Scope

Analysis Core v1 evaluates germline:

- SNVs
- small INDELs
- small MNVs

Structural variants, copy-number variants, repeat expansions, and other
variant classes are outside the primary v0.1 benchmark scope and should
be evaluated separately.

## Eligibility

A public case is eligible only when all of the following are true:

1. The phenopacket contains a SOLVED interpretation.
2. At least three present HPO terms are available.
3. At least one genomic interpretation is explicitly CAUSATIVE.
4. At least one causative variant has an explicit VCF record.
5. The VCF record explicitly declares hg38.
6. Chromosome, position, reference, and alternate allele are available.
7. The causative variant is an SNV, small INDEL, or small MNV.
8. A causal gene symbol is explicitly available from geneContext.
9. Source provenance can be retained.

Absent/excluded HPO terms are retained when available.

Cases with two usable causative variants remain eligible. Both variants
must remain represented in benchmark truth.

## Canonical Variant Identity

The VCF representation is the canonical genomic identity:

    chrom:pos:ref:alt

HGVS expressions are retained as source/provenance information but are
not used as a replacement for canonical VCF identity.

## Truth Separation

Benchmark truth MUST be physically and logically separated from ranking
input.

Ranking-time input:

    case.json

Benchmark-only truth:

    truth.json

The ranking-time representation MUST NOT contain:

- causal candidate identifiers
- is_causal flags
- solved diagnosis labels used as truth
- expected ranks
- benchmark labels
- any equivalent truth-bearing field

Ranking functions MUST NOT accept benchmark truth as an argument.

Truth may be loaded only after ranking results have been frozen for
benchmark evaluation.

## Frozen Evidence

Candidate evidence used for the reproducible benchmark is stored as a
frozen production snapshot.

Every candidate snapshot must retain:

- snapshot date
- pipeline version
- evidence source versions
- canonical production CandidateCase state

A frozen-evidence benchmark and a future current-evidence reanalysis are
different experiments and must not be conflated.

## Cohort Size

The target public benchmark contains 500 cases:

- 300 development cases
- 200 locked evaluation cases

The development cohort may be used for engineering, debugging, failure
analysis, and explicitly documented algorithm development.

The locked evaluation cohort must not be used to tune ranking behavior.

## Cohort Construction

The cohort is selected deterministically from the strict eligible pool.

Selection should maximize biological and provenance diversity rather
than randomly reproduce the frequency distribution of the source store.

Primary objectives:

1. maximize causal-gene diversity
2. maximize publication diversity
3. avoid domination by large source cohorts
4. preserve meaningful inheritance and variant diversity
5. retain both single-variant and multi-variant solved cases

Exact duplicate cases and exact duplicate causal variant sets must not
appear across development and locked evaluation.

## Split Leakage

Development and locked evaluation should have no publication overlap
whenever the requested cohort sizes can be achieved under that
constraint.

Exact causal variant overlap across the two splits is prohibited.

Gene overlap is permitted and measured rather than universally
prohibited.

Locked evaluation results must therefore be stratified into:

- genes seen in development
- genes unseen in development

This measures generalization without redefining the primary cohort
after results are observed.

## Ranking Strategies

The benchmark compares frozen production strategies, including:

- weighted variant-first ranking
- weighted exact-phenotype ranking
- weighted semantic-phenotype ranking
- Borda aggregation
- experimental geometric-mean rank aggregation

Benchmark truth is not available to any ranking strategy.

## Primary Metrics

Primary case-level retrieval metrics:

- mean reciprocal rank (MRR)
- Hit@1
- Hit@3
- Hit@5
- causal rank distribution

## Secondary Metrics

Secondary analyses include:

- unique Hit@k
- top-tie rate
- evidence coverage
- paired strategy rank changes
- seen-gene versus unseen-gene performance
- variant-type stratification
- inheritance/allelic-state stratification
- negative-HPO availability stratification
- failure taxonomy

For cases with multiple causal variants, additional recovery metrics
must distinguish finding any causal variant from recovering all required
causal variants.

## Interpretation Boundaries

Benchmark scores represent candidate retrieval/prioritization
performance only.

They are not:

- pathogenicity probabilities
- ACMG/AMP classifications
- diagnoses
- clinical recommendations
- proof of clinical validity
- proof of clinical utility

## Reproducibility

The following must be versioned or recorded for every benchmark run:

- benchmark protocol version
- source dataset release
- cohort manifest
- development/locked split assignment
- pipeline version
- ranking engine versions
- evidence snapshot metadata
- benchmark artifact schema version

Changes to cohort eligibility or split construction after observing
locked evaluation performance require a new benchmark protocol version.

## Cohort v0.1 Diversity Freeze

Before any ranking performance was evaluated, the deterministic v0.1
cohort construction was dry-run against the strict eligible source pool.

The strict eligible pool contained:

- 8,610 cases
- 675 causal genes
- 1,656 publications

The deterministic 500-case cohort selected for the v0.1 design contains:

- 500 cases
- 500 unique causal genes
- 500 unique publications
- 84 multi-variant solved cases
- 405 cases with at least one absent/excluded HPO term

The development/locked design is therefore intentionally a
cross-gene, cross-publication evaluation:

- development: 300 cases / 300 unique causal genes / 300 publications
- locked evaluation: 200 cases / 200 unique causal genes / 200 publications
- publication overlap: 0
- exact causal variant overlap: 0
- causal gene overlap: 0

Accordingly, all 200 locked-evaluation causal genes are unseen in the
development cohort.

This property was established before ranking performance on the public
cohort was observed. It is part of the v0.1 benchmark design rather than
a post hoc performance filter.

The v0.1 benchmark therefore emphasizes breadth and cross-gene
generalization rather than reproducing the frequency distribution of
the Phenopacket Store.

Future benchmark versions may add a separate repeated-gene cohort to
measure within-gene generalization. Such a cohort must be reported
separately and must not replace or retroactively modify the locked v0.1
evaluation cohort.
