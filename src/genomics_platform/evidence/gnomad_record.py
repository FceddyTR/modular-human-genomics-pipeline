from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Optional


SCHEMA_VERSION = "1.0"


def safe_af(ac: Optional[int], an: Optional[int]) -> Optional[float]:
    """
    Calculate allele frequency only when AC/AN are usable.

    AN == 0 or missing data yields None rather than an inferred frequency.
    """
    if ac is None or an is None or an <= 0:
        return None
    return ac / an


@dataclass(frozen=True)
class PopulationFrequency:
    id: str
    ac: Optional[int] = None
    an: Optional[int] = None
    homozygote_count: Optional[int] = None

    @property
    def af(self) -> Optional[float]:
        return safe_af(self.ac, self.an)

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["af"] = self.af
        return data


@dataclass(frozen=True)
class DatasetFrequency:
    ac: Optional[int] = None
    an: Optional[int] = None
    homozygote_count: Optional[int] = None
    populations: tuple[PopulationFrequency, ...] = field(default_factory=tuple)

    @property
    def af(self) -> Optional[float]:
        return safe_af(self.ac, self.an)

    def to_dict(self) -> dict[str, Any]:
        return {
            "ac": self.ac,
            "an": self.an,
            "af": self.af,
            "homozygote_count": self.homozygote_count,
            "populations": [p.to_dict() for p in self.populations],
        }


@dataclass(frozen=True)
class GnomADRecord:
    assembly: str
    chrom: str
    pos: int
    ref: str
    alt: str

    lookup_status: str

    exome: Optional[DatasetFrequency] = None
    genome: Optional[DatasetFrequency] = None

    source: str = "gnomAD"
    source_dataset: str = "gnomad_r4"
    source_release: Optional[str] = None
    access_method: str = "gnomAD Browser GraphQL API"
    schema_version: str = SCHEMA_VERSION

    error: Optional[str] = None

    @property
    def variant_key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @property
    def gnomad_variant_id(self) -> str:
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "source": self.source,
            "source_dataset": self.source_dataset,
            "source_release": self.source_release,
            "access_method": self.access_method,
            "assembly": self.assembly,
            "variant_key": self.variant_key,
            "gnomad_variant_id": self.gnomad_variant_id,
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "lookup_status": self.lookup_status,
            "exome": self.exome.to_dict() if self.exome else None,
            "genome": self.genome.to_dict() if self.genome else None,
            "error": self.error,
            "semantics": {
                "not_found": "No gnomAD record returned; this does not imply AF=0 or rarity.",
                "frequency": "Population frequency is evidence, not pathogenicity classification.",
                "datasets": "Exome and genome observations are preserved separately.",
                "populations": (
                    "Population IDs are preserved as returned by gnomAD; "
                    "no grpmax or ancestry interpretation is inferred from arbitrary subgroups."
                ),
            },
        }
