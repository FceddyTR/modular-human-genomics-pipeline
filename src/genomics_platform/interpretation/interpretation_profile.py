from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Optional


SCHEMA_VERSION = "1.0"

VALID_DIMENSION_STATUSES = {
    "AVAILABLE",
    "PARTIAL",
    "UNAVAILABLE",
    "ERROR",
}


@dataclass(frozen=True)
class InterpretationObservation:
    """
    One descriptive observation derived from normalized evidence.

    Observations are not ACMG/AMP criteria and do not independently
    establish pathogenicity.
    """

    code: str
    label: str
    value: Any = None
    source: Optional[str] = None
    details: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class EvidenceDimension:
    """
    Descriptive interpretation of one evidence dimension.

    The status describes evidence availability, not biological or
    clinical significance.
    """

    name: str
    status: str
    observations: tuple[InterpretationObservation, ...] = field(
        default_factory=tuple
    )
    limitations: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        normalized_status = self.status.upper()

        if normalized_status not in VALID_DIMENSION_STATUSES:
            raise ValueError(
                f"Invalid interpretation dimension status "
                f"for {self.name}: {self.status}"
            )

        object.__setattr__(
            self,
            "status",
            normalized_status,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "status": self.status,
            "observation_count": len(self.observations),
            "observations": [
                observation.to_dict()
                for observation in self.observations
            ],
            "limitations": list(self.limitations),
        }


@dataclass(frozen=True)
class InterpretationProfile:
    """
    Source-normalized descriptive interpretation profile.

    This object intentionally contains no:
      - pathogenicity classification,
      - ACMG/AMP classification,
      - candidate tier,
      - ranking score,
      - pathogenicity probability.

    Ranking and classification belong to later policy layers.
    """

    variant_key: str

    population: EvidenceDimension
    clinical: EvidenceDimension
    functional: EvidenceDimension

    phenotype: Optional[EvidenceDimension] = None
    inheritance: Optional[EvidenceDimension] = None
    gene_disease: Optional[EvidenceDimension] = None

    schema_version: str = SCHEMA_VERSION

    @property
    def error_dimensions(self) -> list[str]:
        dimensions = (
            self.population,
            self.clinical,
            self.functional,
            self.phenotype,
            self.inheritance,
            self.gene_disease,
        )

        return [
            dimension.name
            for dimension in dimensions
            if dimension is not None
            and dimension.status == "ERROR"
        ]

    @property
    def available_dimensions(self) -> list[str]:
        dimensions = (
            self.population,
            self.clinical,
            self.functional,
            self.phenotype,
            self.inheritance,
            self.gene_disease,
        )

        return [
            dimension.name
            for dimension in dimensions
            if dimension is not None
            and dimension.status in {
                "AVAILABLE",
                "PARTIAL",
            }
        ]

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "variant_key": self.variant_key,
            "dimensions": {
                "population": self.population.to_dict(),
                "clinical": self.clinical.to_dict(),
                "functional": self.functional.to_dict(),
                "phenotype": (
                    self.phenotype.to_dict()
                    if self.phenotype
                    else None
                ),
                "inheritance": (
                    self.inheritance.to_dict()
                    if self.inheritance
                    else None
                ),
                "gene_disease": (
                    self.gene_disease.to_dict()
                    if self.gene_disease
                    else None
                ),
            },
            "available_dimensions": self.available_dimensions,
            "error_dimensions": self.error_dimensions,
            "semantics": {
                "profile": (
                    "Descriptive interpretation profile derived "
                    "from normalized evidence."
                ),
                "status": (
                    "Dimension status describes evidence availability "
                    "and completeness, not pathogenicity."
                ),
                "observations": (
                    "Observations describe evidence signals and must "
                    "not be interpreted as ACMG/AMP criteria unless "
                    "a separate validated policy layer explicitly "
                    "implements such rules."
                ),
                "ranking": (
                    "This profile does not rank variants or assign "
                    "candidate tiers."
                ),
                "classification": (
                    "This profile does not classify variants as "
                    "pathogenic, likely pathogenic, benign, or "
                    "likely benign."
                ),
            },
        }
