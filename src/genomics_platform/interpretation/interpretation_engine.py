from __future__ import annotations

from genomics_platform.evidence.variant_evidence_contract import (
    VariantEvidenceContract,
)
from genomics_platform.interpretation.clinical_interpreter import (
    interpret_clinical,
)
from genomics_platform.interpretation.functional_interpreter import (
    interpret_functional,
)
from genomics_platform.interpretation.interpretation_profile import (
    InterpretationProfile,
)
from genomics_platform.interpretation.population_interpreter import (
    interpret_population,
)


ENGINE_VERSION = "0.1.0"


def interpret_variant(
    contract: VariantEvidenceContract,
) -> InterpretationProfile:
    """
    Convert a normalized VariantEvidenceContract into a descriptive
    InterpretationProfile.

    Current v0.1 dimensions:
      - population
      - clinical
      - functional

    Future dimensions:
      - phenotype
      - inheritance
      - gene-disease validity

    This engine intentionally does NOT:
      - implement ACMG/AMP classification,
      - assign pathogenicity probabilities,
      - assign candidate tiers,
      - rank variants,
      - convert source annotations into diagnoses.

    Those functions belong to explicit downstream policy and ranking
    layers.
    """

    if not isinstance(contract, VariantEvidenceContract):
        raise TypeError(
            "interpret_variant expects a VariantEvidenceContract"
        )

    return InterpretationProfile(
        variant_key=contract.identity.variant_key,
        population=interpret_population(contract),
        clinical=interpret_clinical(contract),
        functional=interpret_functional(contract),
    )
