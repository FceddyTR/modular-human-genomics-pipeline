"""Truth-aware Phenopacket discovery for public benchmark construction.

This module is used only to determine benchmark eligibility and extract
benchmark truth. Its output is not a ranking-time case representation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


PUBLIC_CASE_PROTOCOL_VERSION = "0.1"
PHENOPACKET_STORE_RELEASE = "0.1.27"

SUPPORTED_ASSEMBLY = "hg38"
SUPPORTED_VARIANT_TYPES = frozenset(
    {
        "SNV",
        "INDEL",
        "MNV",
    }
)


@dataclass(frozen=True)
class PublicCausalVariant:
    gene_symbol: str
    hgnc_id: str | None
    chrom: str
    pos: int
    ref: str
    alt: str
    variant_type: str
    allelic_state: str | None = None

    def __post_init__(self) -> None:
        gene_symbol = self.gene_symbol.strip()

        if not gene_symbol:
            raise ValueError(
                "gene_symbol is required."
            )

        if self.pos <= 0:
            raise ValueError(
                "pos must be positive."
            )

        if not self.chrom.strip():
            raise ValueError(
                "chrom is required."
            )

        if not self.ref:
            raise ValueError(
                "ref is required."
            )

        if not self.alt:
            raise ValueError(
                "alt is required."
            )

        if (
            self.variant_type
            not in SUPPORTED_VARIANT_TYPES
        ):
            raise ValueError(
                "Unsupported variant_type."
            )

        object.__setattr__(
            self,
            "gene_symbol",
            gene_symbol.upper(),
        )

    @property
    def variant_key(self) -> str:
        return (
            f"{self.chrom}:"
            f"{self.pos}:"
            f"{self.ref}:"
            f"{self.alt}"
        )


@dataclass(frozen=True)
class PublicCaseEligibility:
    case_id: str
    eligible: bool
    present_hpo_terms: tuple[str, ...]
    absent_hpo_terms: tuple[str, ...]
    causal_variants: tuple[
        PublicCausalVariant, ...
    ]
    exclusion_reasons: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.case_id.strip():
            raise ValueError(
                "case_id is required."
            )

        if self.eligible:
            if self.exclusion_reasons:
                raise ValueError(
                    "Eligible cases cannot have "
                    "exclusion reasons."
                )

            if len(self.present_hpo_terms) < 3:
                raise ValueError(
                    "Eligible cases require at "
                    "least three present HPO terms."
                )

            if not self.causal_variants:
                raise ValueError(
                    "Eligible cases require at "
                    "least one causal variant."
                )

    @property
    def causal_genes(self) -> tuple[str, ...]:
        return tuple(
            sorted(
                {
                    variant.gene_symbol
                    for variant
                    in self.causal_variants
                }
            )
        )


def _deduplicate(
    values: list[str],
) -> tuple[str, ...]:
    return tuple(
        dict.fromkeys(values)
    )


def _phenotypes(
    packet: dict[str, Any],
) -> tuple[
    tuple[str, ...],
    tuple[str, ...],
]:
    present: list[str] = []
    absent: list[str] = []

    for feature in packet.get(
        "phenotypicFeatures",
        [],
    ):
        if not isinstance(feature, dict):
            continue

        term = feature.get("type")

        if not isinstance(term, dict):
            continue

        term_id = term.get("id")

        if not (
            isinstance(term_id, str)
            and term_id.startswith("HP:")
            and len(term_id) == 10
            and term_id[3:].isdigit()
        ):
            continue

        if feature.get("excluded") is True:
            absent.append(term_id)
        else:
            present.append(term_id)

    return (
        _deduplicate(present),
        _deduplicate(absent),
    )


def _solved_genomic_interpretations(
    packet: dict[str, Any],
):
    for interpretation in packet.get(
        "interpretations",
        [],
    ):
        if not isinstance(
            interpretation,
            dict,
        ):
            continue

        if (
            interpretation.get(
                "progressStatus"
            )
            != "SOLVED"
        ):
            continue

        diagnosis = interpretation.get(
            "diagnosis"
        )

        if not isinstance(
            diagnosis,
            dict,
        ):
            continue

        genomic = diagnosis.get(
            "genomicInterpretations",
            [],
        )

        if not isinstance(genomic, list):
            continue

        yield from genomic


def _variant_type(
    ref: str,
    alt: str,
) -> str | None:
    if (
        ref.startswith("<")
        or alt.startswith("<")
        or "[" in alt
        or "]" in alt
        or alt == "*"
    ):
        return None

    if len(ref) == 1 and len(alt) == 1:
        return "SNV"

    if max(len(ref), len(alt)) > 50:
        return None

    if len(ref) != len(alt):
        return "INDEL"

    return "MNV"


def _causal_variant(
    genomic: Any,
) -> PublicCausalVariant | None:
    if not isinstance(genomic, dict):
        return None

    if (
        str(
            genomic.get(
                "interpretationStatus",
                "",
            )
        ).upper()
        != "CAUSATIVE"
    ):
        return None

    variant_interpretation = (
        genomic.get(
            "variantInterpretation"
        )
    )

    if not isinstance(
        variant_interpretation,
        dict,
    ):
        return None

    descriptor = (
        variant_interpretation.get(
            "variationDescriptor"
        )
    )

    if not isinstance(
        descriptor,
        dict,
    ):
        return None

    vcf = descriptor.get("vcfRecord")

    if not isinstance(vcf, dict):
        return None

    if (
        str(
            vcf.get(
                "genomeAssembly",
                "",
            )
        ).lower()
        != SUPPORTED_ASSEMBLY
    ):
        return None

    chrom = vcf.get("chrom")
    pos = vcf.get("pos")
    ref = vcf.get("ref")
    alt = vcf.get("alt")

    if not (
        isinstance(chrom, str)
        and chrom
        and isinstance(pos, int)
        and not isinstance(pos, bool)
        and pos > 0
        and isinstance(ref, str)
        and ref
        and isinstance(alt, str)
        and alt
    ):
        return None

    variant_type = _variant_type(
        ref,
        alt,
    )

    if variant_type is None:
        return None

    gene_context = descriptor.get(
        "geneContext"
    )

    if not isinstance(
        gene_context,
        dict,
    ):
        return None

    gene_symbol = gene_context.get(
        "symbol"
    )

    if not (
        isinstance(gene_symbol, str)
        and gene_symbol.strip()
    ):
        return None

    hgnc_id = gene_context.get(
        "valueId"
    )

    if not isinstance(hgnc_id, str):
        hgnc_id = None

    allelic_state = descriptor.get(
        "allelicState"
    )

    allelic_state_value = None

    if isinstance(allelic_state, dict):
        value = (
            allelic_state.get("label")
            or allelic_state.get("id")
        )

        if isinstance(value, str):
            allelic_state_value = value

    return PublicCausalVariant(
        gene_symbol=gene_symbol,
        hgnc_id=hgnc_id,
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        variant_type=variant_type,
        allelic_state=allelic_state_value,
    )


def evaluate_phenopacket(
    packet: dict[str, Any],
) -> PublicCaseEligibility:
    """Evaluate one source packet without creating ranking input."""

    if not isinstance(packet, dict):
        raise TypeError(
            "packet must be a dictionary."
        )

    case_id = packet.get("id")

    if not (
        isinstance(case_id, str)
        and case_id.strip()
    ):
        raise ValueError(
            "Phenopacket id is required."
        )

    present, absent = _phenotypes(
        packet
    )

    reasons: list[str] = []

    if len(present) < 3:
        reasons.append(
            "present_hpo_lt_3"
        )

    solved = list(
        _solved_genomic_interpretations(
            packet
        )
    )

    if not solved:
        reasons.append(
            "no_solved_genomic_interpretation"
        )

    causal_variants = tuple(
        variant
        for genomic in solved
        if (
            variant := _causal_variant(
                genomic
            )
        )
        is not None
    )

    if not causal_variants:
        reasons.append(
            "no_usable_causative_small_variant"
        )

    return PublicCaseEligibility(
        case_id=case_id,
        eligible=not reasons,
        present_hpo_terms=present,
        absent_hpo_terms=absent,
        causal_variants=causal_variants,
        exclusion_reasons=tuple(
            dict.fromkeys(reasons)
        ),
    )
