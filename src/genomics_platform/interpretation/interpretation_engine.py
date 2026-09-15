from __future__ import annotations

from typing import Any, Dict, Optional

from genomics_platform.evidence.variant_evidence_contract import (
    VariantEvidenceContract,
)
from genomics_platform.interpretation.clinical_interpreter import (
    interpret_clinical,
)
from genomics_platform.interpretation.functional_interpreter import (
    interpret_functional,
)
from genomics_platform.interpretation.gene_disease_interpreter import (
    interpret_gene_disease,
)
from genomics_platform.interpretation.interpretation_profile import (
    InterpretationProfile,
)
from genomics_platform.interpretation.population_interpreter import (
    interpret_population,
)


ENGINE_VERSION = "0.2.0"


def interpret_variant(
    contract: VariantEvidenceContract,
    gene_evidence_context: Optional[Dict[str, Any]] = None,
) -> InterpretationProfile:
    """
    Convert a normalized VariantEvidenceContract into a descriptive
    InterpretationProfile.

    Current dimensions:
      - population
      - clinical
      - functional
      - optional gene-disease context

    Future dimensions:
      - phenotype
      - inheritance

    Gene-disease evidence is supplied as already-built context.
    This keeps the interpretation layer independent from databases,
    files, network services, and evidence-source storage details.

    This engine intentionally does NOT:
      - implement ACMG/AMP classification,
      - assign pathogenicity probabilities,
      - assign candidate tiers,
      - rank variants,
      - convert gene-disease evidence into variant pathogenicity,
      - convert source annotations into diagnoses.
    """

    if not isinstance(contract, VariantEvidenceContract):
        raise TypeError(
            "interpret_variant expects a VariantEvidenceContract"
        )

    gene_disease = None

    if gene_evidence_context is not None:
        gene_disease = interpret_gene_disease(
            gene_evidence_context
        )

    return InterpretationProfile(
        variant_key=contract.identity.variant_key,
        population=interpret_population(contract),
        clinical=interpret_clinical(contract),
        functional=interpret_functional(contract),
        gene_disease=gene_disease,
    )
