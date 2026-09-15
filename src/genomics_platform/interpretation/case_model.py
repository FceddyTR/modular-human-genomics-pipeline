"""Stable candidate case model for downstream ranking.

Case Model v0.1 binds:
- a VariantEvidenceContract,
- its InterpretationProfile,
- optional patient/case context,
- explicit semantic boundaries.

It does NOT:
- rank candidates,
- assign tiers,
- classify variants,
- perform ACMG/AMP classification,
- calculate pathogenicity probability,
- assign a diagnosis.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable, Optional

from genomics_platform.evidence.variant_evidence_contract import (
    VariantEvidenceContract,
)
from genomics_platform.interpretation.interpretation_profile import (
    InterpretationProfile,
)


SCHEMA_VERSION = "1.0"
MODEL_VERSION = "0.1.0"


def _normalize_hpo_terms(
    values: Optional[Iterable[str]],
) -> tuple[str, ...]:
    if values is None:
        return ()

    if isinstance(values, str):
        values = [values]

    output: list[str] = []
    seen: set[str] = set()

    for value in values:
        if not isinstance(value, str):
            raise TypeError(
                "HPO terms must be strings."
            )

        term = value.strip().upper()

        if (
            not term.startswith("HP:")
            or len(term) != 10
            or not term[3:].isdigit()
        ):
            raise ValueError(
                f"Invalid HPO identifier: {value!r}"
            )

        if term not in seen:
            seen.add(term)
            output.append(term)

    return tuple(output)


def _normalize_optional_text(
    value: Optional[str],
) -> Optional[str]:
    if value is None:
        return None

    if not isinstance(value, str):
        raise TypeError(
            "Expected string or None."
        )

    value = value.strip()

    return value or None


@dataclass(frozen=True)
class CaseContext:
    """Patient/case information relevant to interpretation."""

    present_hpo_terms: tuple[str, ...] = ()
    absent_hpo_terms: tuple[str, ...] = ()
    proband_sex: Optional[str] = None
    genotype_state: Optional[str] = None
    de_novo_status: Optional[bool] = None
    family_context: Optional[str] = None

    def __post_init__(self) -> None:
        present = _normalize_hpo_terms(
            self.present_hpo_terms
        )
        absent = _normalize_hpo_terms(
            self.absent_hpo_terms
        )

        overlap = set(present) & set(absent)

        if overlap:
            raise ValueError(
                "HPO terms cannot be both present "
                "and absent: "
                + ", ".join(sorted(overlap))
            )

        object.__setattr__(
            self,
            "present_hpo_terms",
            present,
        )

        object.__setattr__(
            self,
            "absent_hpo_terms",
            absent,
        )

        sex = _normalize_optional_text(
            self.proband_sex
        )

        if sex is not None:
            sex = sex.lower()

            if sex not in {
                "male",
                "female",
                "unknown",
            }:
                raise ValueError(
                    "proband_sex must be male, "
                    "female, unknown, or None."
                )

        object.__setattr__(
            self,
            "proband_sex",
            sex,
        )

        genotype = _normalize_optional_text(
            self.genotype_state
        )

        if genotype is not None:
            genotype = genotype.lower()

            if genotype not in {
                "heterozygous",
                "homozygous_alt",
                "hemizygous",
                "unknown",
            }:
                raise ValueError(
                    "Unsupported genotype_state: "
                    f"{genotype!r}"
                )

        object.__setattr__(
            self,
            "genotype_state",
            genotype,
        )

        if (
            self.de_novo_status is not None
            and not isinstance(
                self.de_novo_status,
                bool,
            )
        ):
            raise TypeError(
                "de_novo_status must be bool "
                "or None."
            )

        object.__setattr__(
            self,
            "family_context",
            _normalize_optional_text(
                self.family_context
            ),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "present_hpo_terms": list(
                self.present_hpo_terms
            ),
            "absent_hpo_terms": list(
                self.absent_hpo_terms
            ),
            "proband_sex": self.proband_sex,
            "genotype_state": (
                self.genotype_state
            ),
            "de_novo_status": (
                self.de_novo_status
            ),
            "family_context": (
                self.family_context
            ),
        }


@dataclass(frozen=True)
class CandidateCase:
    """Stable unit consumed by future ranking engines."""

    candidate_id: str
    evidence_contract: VariantEvidenceContract
    interpretation_profile: InterpretationProfile
    case_context: CaseContext = field(
        default_factory=CaseContext
    )
    gene_symbol: Optional[str] = None
    schema_version: str = SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not isinstance(
            self.evidence_contract,
            VariantEvidenceContract,
        ):
            raise TypeError(
                "evidence_contract must be a "
                "VariantEvidenceContract."
            )

        if not isinstance(
            self.interpretation_profile,
            InterpretationProfile,
        ):
            raise TypeError(
                "interpretation_profile must be an "
                "InterpretationProfile."
            )

        if not isinstance(
            self.case_context,
            CaseContext,
        ):
            raise TypeError(
                "case_context must be a CaseContext."
            )

        candidate_id = (
            _normalize_optional_text(
                self.candidate_id
            )
        )

        if candidate_id is None:
            raise ValueError(
                "candidate_id cannot be empty."
            )

        object.__setattr__(
            self,
            "candidate_id",
            candidate_id,
        )

        gene = _normalize_optional_text(
            self.gene_symbol
        )

        if gene is not None:
            gene = gene.upper()

        object.__setattr__(
            self,
            "gene_symbol",
            gene,
        )

        contract_key = (
            self.evidence_contract
            .identity
            .variant_key
        )

        profile_key = (
            self.interpretation_profile
            .variant_key
        )

        if contract_key != profile_key:
            raise ValueError(
                "Variant identity mismatch between "
                "evidence contract and "
                "interpretation profile: "
                f"{contract_key!r} != "
                f"{profile_key!r}"
            )

    @property
    def variant_key(self) -> str:
        return (
            self.evidence_contract
            .identity
            .variant_key
        )

    @property
    def available_dimensions(
        self,
    ) -> tuple[str, ...]:
        return tuple(
            self.interpretation_profile
            .available_dimensions
        )

    @property
    def error_dimensions(
        self,
    ) -> tuple[str, ...]:
        return tuple(
            self.interpretation_profile
            .error_dimensions
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": (
                self.schema_version
            ),
            "model": {
                "name": "candidate_case",
                "version": MODEL_VERSION,
            },
            "candidate_id": (
                self.candidate_id
            ),
            "variant_key": self.variant_key,
            "gene_symbol": self.gene_symbol,
            "case_context": (
                self.case_context.to_dict()
            ),
            "evidence_contract": (
                self.evidence_contract.to_dict()
            ),
            "interpretation_profile": (
                self.interpretation_profile
                .to_dict()
            ),
            "summary": {
                "available_dimensions": list(
                    self.available_dimensions
                ),
                "error_dimensions": list(
                    self.error_dimensions
                ),
            },
            "semantic_boundaries": {
                "candidate_ranking": False,
                "candidate_tier": False,
                "variant_classification": False,
                "acmg_amp_classification": False,
                "pathogenicity_probability": False,
                "diagnosis": False,
            },
        }


def build_candidate_case(
    *,
    candidate_id: str,
    evidence_contract: VariantEvidenceContract,
    interpretation_profile: InterpretationProfile,
    gene_symbol: Optional[str] = None,
    present_hpo_terms: Optional[
        Iterable[str]
    ] = None,
    absent_hpo_terms: Optional[
        Iterable[str]
    ] = None,
    proband_sex: Optional[str] = None,
    genotype_state: Optional[str] = None,
    de_novo_status: Optional[bool] = None,
    family_context: Optional[str] = None,
) -> CandidateCase:
    """Build a validated CandidateCase."""

    context = CaseContext(
        present_hpo_terms=(
            _normalize_hpo_terms(
                present_hpo_terms
            )
        ),
        absent_hpo_terms=(
            _normalize_hpo_terms(
                absent_hpo_terms
            )
        ),
        proband_sex=proband_sex,
        genotype_state=genotype_state,
        de_novo_status=de_novo_status,
        family_context=family_context,
    )

    return CandidateCase(
        candidate_id=candidate_id,
        evidence_contract=evidence_contract,
        interpretation_profile=(
            interpretation_profile
        ),
        case_context=context,
        gene_symbol=gene_symbol,
    )
