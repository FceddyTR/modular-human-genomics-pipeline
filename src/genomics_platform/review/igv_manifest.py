"""IGV review manifest generation for ranked variants.

The manifest describes files and genomic loci required by an IGV
client. It does not inspect reads, classify variants, or alter ranking.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genomics_platform.evidence.variant_evidence_contract import (
    VariantIdentity,
)


SCHEMA_VERSION = "1.0"
DEFAULT_FLANK_BP = 100


@dataclass(frozen=True)
class IGVAlignmentTrack:
    path: str
    index_path: str
    format: str = "bam"

    def __post_init__(self) -> None:
        fmt = self.format.lower()

        if fmt not in {
            "bam",
            "cram",
        }:
            raise ValueError(
                "Alignment format must be "
                "'bam' or 'cram'."
            )

        if not self.path.strip():
            raise ValueError(
                "Alignment path is required."
            )

        if not self.index_path.strip():
            raise ValueError(
                "Alignment index path is required."
            )

        object.__setattr__(
            self,
            "format",
            fmt,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "index_path": self.index_path,
            "format": self.format,
        }


@dataclass(frozen=True)
class IGVVariantTrack:
    path: str
    index_path: str
    format: str = "vcf"

    def __post_init__(self) -> None:
        if not self.path.strip():
            raise ValueError(
                "Variant track path is required."
            )

        if not self.index_path.strip():
            raise ValueError(
                "Variant track index path is required."
            )

    def to_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "index_path": self.index_path,
            "format": self.format,
        }


@dataclass(frozen=True)
class IGVReviewManifest:
    candidate_id: str
    variant_key: str
    genome_build: str
    chromosome: str
    position: int
    reference_allele: str
    alternate_allele: str
    locus_start: int
    locus_end: int
    alignment: IGVAlignmentTrack
    variants: IGVVariantTrack
    schema_version: str = SCHEMA_VERSION

    @property
    def locus(self) -> str:
        return (
            f"{self.chromosome}:"
            f"{self.locus_start}-"
            f"{self.locus_end}"
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "candidate_id": self.candidate_id,
            "variant_key": self.variant_key,
            "genome_build": self.genome_build,
            "chromosome": self.chromosome,
            "position": self.position,
            "reference_allele": (
                self.reference_allele
            ),
            "alternate_allele": (
                self.alternate_allele
            ),
            "locus_start": self.locus_start,
            "locus_end": self.locus_end,
            "locus": self.locus,
            "alignment": (
                self.alignment.to_dict()
            ),
            "variants": (
                self.variants.to_dict()
            ),
            "review_semantics": {
                "post_ranking": True,
                "technical_visual_review": True,
                "affects_ranking": False,
                "pathogenicity_classification": False,
                "diagnosis": False,
            },
        }


def parse_variant_key(
    variant_key: str,
    *,
    assembly: str,
) -> VariantIdentity:
    """Parse the canonical chrom:pos:ref:alt key.

    This fallback exists for contracts such as RankedCandidate that
    intentionally retain only the canonical variant key.
    """

    value = variant_key.strip()

    parts = value.split(
        ":",
        maxsplit=3,
    )

    if len(parts) != 4:
        raise ValueError(
            "variant_key must have canonical "
            "chrom:pos:ref:alt format."
        )

    chrom, pos_text, ref, alt = parts

    if not chrom:
        raise ValueError(
            "Variant chromosome is required."
        )

    try:
        pos = int(pos_text)
    except ValueError as exc:
        raise ValueError(
            "Variant position must be an integer."
        ) from exc

    if pos < 1:
        raise ValueError(
            "Variant position must be >= 1."
        )

    if not ref or not alt:
        raise ValueError(
            "Reference and alternate alleles "
            "are required."
        )

    return VariantIdentity(
        assembly=assembly,
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
    )


def build_igv_manifest(
    *,
    candidate_id: str,
    identity: VariantIdentity,
    alignment_path: str | Path,
    alignment_index_path: str | Path,
    variant_path: str | Path,
    variant_index_path: str | Path,
    flank_bp: int = DEFAULT_FLANK_BP,
    alignment_format: str = "bam",
) -> IGVReviewManifest:
    """Build an IGV manifest from structured variant identity."""

    candidate_id = candidate_id.strip()

    if not candidate_id:
        raise ValueError(
            "candidate_id is required."
        )

    if flank_bp < 0:
        raise ValueError(
            "flank_bp cannot be negative."
        )

    # Include the full reference allele span for deletions/MNVs.
    variant_end = (
        identity.pos
        + max(
            len(identity.ref),
            1,
        )
        - 1
    )

    locus_start = max(
        1,
        identity.pos - flank_bp,
    )

    locus_end = (
        variant_end + flank_bp
    )

    alignment = IGVAlignmentTrack(
        path=str(alignment_path),
        index_path=str(
            alignment_index_path
        ),
        format=alignment_format,
    )

    variants = IGVVariantTrack(
        path=str(variant_path),
        index_path=str(
            variant_index_path
        ),
    )

    return IGVReviewManifest(
        candidate_id=candidate_id,
        variant_key=identity.variant_key,
        genome_build=identity.assembly,
        chromosome=identity.chrom,
        position=identity.pos,
        reference_allele=identity.ref,
        alternate_allele=identity.alt,
        locus_start=locus_start,
        locus_end=locus_end,
        alignment=alignment,
        variants=variants,
    )


def build_igv_manifest_from_variant_key(
    *,
    candidate_id: str,
    variant_key: str,
    assembly: str,
    alignment_path: str | Path,
    alignment_index_path: str | Path,
    variant_path: str | Path,
    variant_index_path: str | Path,
    flank_bp: int = DEFAULT_FLANK_BP,
    alignment_format: str = "bam",
) -> IGVReviewManifest:
    """Build a manifest from a ranked candidate's canonical key."""

    identity = parse_variant_key(
        variant_key,
        assembly=assembly,
    )

    return build_igv_manifest(
        candidate_id=candidate_id,
        identity=identity,
        alignment_path=alignment_path,
        alignment_index_path=(
            alignment_index_path
        ),
        variant_path=variant_path,
        variant_index_path=(
            variant_index_path
        ),
        flank_bp=flank_bp,
        alignment_format=alignment_format,
    )
