"""Tests for post-ranking IGV review manifests."""

from __future__ import annotations

import pytest

from genomics_platform.evidence.variant_evidence_contract import (
    VariantIdentity,
)
from genomics_platform.review.igv_manifest import (
    build_igv_manifest,
    build_igv_manifest_from_variant_key,
    parse_variant_key,
)
from genomics_platform.review.review_contract import (
    ReviewStatus,
    TechnicalReviewRecord,
)


def test_manifest_from_structured_identity():
    identity = VariantIdentity(
        assembly="GRCh38",
        chrom="chr17",
        pos=43071077,
        ref="G",
        alt="A",
    )

    manifest = build_igv_manifest(
        candidate_id="candidate-001",
        identity=identity,
        alignment_path="sample.bam",
        alignment_index_path="sample.bam.bai",
        variant_path="sample.vcf.gz",
        variant_index_path="sample.vcf.gz.tbi",
        flank_bp=100,
    )

    assert (
        manifest.variant_key
        == "chr17:43071077:G:A"
    )

    assert (
        manifest.locus
        == "chr17:43070977-43071177"
    )

    assert (
        manifest.genome_build
        == "GRCh38"
    )


def test_indel_window_includes_reference_span():
    identity = VariantIdentity(
        assembly="GRCh38",
        chrom="chr1",
        pos=1000,
        ref="ATGC",
        alt="A",
    )

    manifest = build_igv_manifest(
        candidate_id="deletion",
        identity=identity,
        alignment_path="sample.bam",
        alignment_index_path="sample.bam.bai",
        variant_path="sample.vcf.gz",
        variant_index_path="sample.vcf.gz.tbi",
        flank_bp=10,
    )

    assert manifest.locus_start == 990
    assert manifest.locus_end == 1013


def test_window_never_starts_below_one():
    identity = VariantIdentity(
        assembly="GRCh38",
        chrom="chr1",
        pos=5,
        ref="A",
        alt="G",
    )

    manifest = build_igv_manifest(
        candidate_id="near-start",
        identity=identity,
        alignment_path="sample.bam",
        alignment_index_path="sample.bam.bai",
        variant_path="sample.vcf.gz",
        variant_index_path="sample.vcf.gz.tbi",
        flank_bp=100,
    )

    assert manifest.locus_start == 1


def test_ranked_candidate_key_fallback():
    manifest = (
        build_igv_manifest_from_variant_key(
            candidate_id="candidate-002",
            variant_key="chr20:1000:C:T",
            assembly="GRCh38",
            alignment_path="sample.bam",
            alignment_index_path=(
                "sample.bam.bai"
            ),
            variant_path="sample.vcf.gz",
            variant_index_path=(
                "sample.vcf.gz.tbi"
            ),
        )
    )

    assert manifest.chromosome == "chr20"
    assert manifest.position == 1000
    assert manifest.reference_allele == "C"
    assert manifest.alternate_allele == "T"


def test_invalid_variant_key_is_rejected():
    with pytest.raises(
        ValueError,
        match="chrom:pos:ref:alt",
    ):
        parse_variant_key(
            "chr20:1000:C",
            assembly="GRCh38",
        )


def test_invalid_variant_position_is_rejected():
    with pytest.raises(
        ValueError,
        match="integer",
    ):
        parse_variant_key(
            "chr20:not-a-position:C:T",
            assembly="GRCh38",
        )


def test_review_defaults_to_not_reviewed():
    record = TechnicalReviewRecord(
        candidate_id="candidate-001",
        variant_key="chr17:43071077:G:A",
    )

    assert (
        record.status
        is ReviewStatus.NOT_REVIEWED
    )

    assert record.reviewed_at is None


def test_completed_review_has_timestamp():
    record = TechnicalReviewRecord.reviewed(
        candidate_id="candidate-001",
        variant_key="chr17:43071077:G:A",
        status=ReviewStatus.SUPPORTS_CALL,
        reviewer="reviewer-01",
        notes="Read support visually consistent.",
    )

    assert (
        record.status
        is ReviewStatus.SUPPORTS_CALL
    )

    assert record.reviewed_at is not None


def test_review_cannot_affect_ranking():
    record = TechnicalReviewRecord.reviewed(
        candidate_id="candidate-001",
        variant_key="chr17:43071077:G:A",
        status=(
            ReviewStatus.SUSPECTED_ARTIFACT
        ),
    )

    data = record.to_dict()

    assert (
        data["semantic_boundaries"][
            "affects_ranking"
        ]
        is False
    )

    assert (
        data["semantic_boundaries"][
            "pathogenicity_probability"
        ]
        is False
    )

    assert (
        data["semantic_boundaries"][
            "diagnosis"
        ]
        is False
    )


def test_manifest_declares_review_boundaries():
    identity = VariantIdentity(
        assembly="GRCh38",
        chrom="chr2",
        pos=2000,
        ref="A",
        alt="T",
    )

    manifest = build_igv_manifest(
        candidate_id="candidate-003",
        identity=identity,
        alignment_path="sample.bam",
        alignment_index_path="sample.bam.bai",
        variant_path="sample.vcf.gz",
        variant_index_path="sample.vcf.gz.tbi",
    )

    boundaries = (
        manifest.to_dict()[
            "review_semantics"
        ]
    )

    assert (
        boundaries["post_ranking"]
        is True
    )

    assert (
        boundaries["affects_ranking"]
        is False
    )

    assert (
        boundaries[
            "pathogenicity_classification"
        ]
        is False
    )


def test_cram_track_is_supported():
    identity = VariantIdentity(
        assembly="GRCh38",
        chrom="chr3",
        pos=3000,
        ref="G",
        alt="C",
    )

    manifest = build_igv_manifest(
        candidate_id="cram-case",
        identity=identity,
        alignment_path="sample.cram",
        alignment_index_path="sample.cram.crai",
        variant_path="sample.vcf.gz",
        variant_index_path="sample.vcf.gz.tbi",
        alignment_format="cram",
    )

    assert (
        manifest.alignment.format
        == "cram"
    )


def test_manifest_does_not_require_real_files():
    identity = VariantIdentity(
        assembly="GRCh38",
        chrom="chr4",
        pos=4000,
        ref="T",
        alt="C",
    )

    manifest = build_igv_manifest(
        candidate_id="manifest-only",
        identity=identity,
        alignment_path="/not/mounted/sample.bam",
        alignment_index_path=(
            "/not/mounted/sample.bam.bai"
        ),
        variant_path="/not/mounted/sample.vcf.gz",
        variant_index_path=(
            "/not/mounted/sample.vcf.gz.tbi"
        ),
    )

    assert (
        manifest.alignment.path
        == "/not/mounted/sample.bam"
    )
