"""Top-level descriptive interpretation engine.

The engine converts a stable Variant Evidence Contract into an
Interpretation Profile.

Optional case-level evidence contexts can be supplied for additional
descriptive interpretation dimensions.

This engine does NOT perform:
- ACMG/AMP classification,
- pathogenicity probability estimation,
- candidate ranking,
- diagnosis,
- treatment recommendation.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, Optional

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
from genomics_platform.interpretation.phenotype_interpreter import (
    interpret_phenotype,
)
from genomics_platform.interpretation.population_interpreter import (
    interpret_population,
)


ENGINE_VERSION = "0.3.0"


def interpret_variant(
    contract: VariantEvidenceContract,
    gene_evidence_context: Optional[
        Dict[str, Any]
    ] = None,
    present_hpo_terms: Optional[
        Iterable[str]
    ] = None,
    absent_hpo_terms: Optional[
        Iterable[str]
    ] = None,
) -> InterpretationProfile:
    """Build a descriptive Interpretation Profile.

    Parameters
    ----------
    contract:
        Stable Variant Evidence Contract.

    gene_evidence_context:
        Optional pre-built Gene Evidence Bundle.

    present_hpo_terms:
        Optional HPO terms observed in the patient.

    absent_hpo_terms:
        Optional HPO terms explicitly assessed and absent.

    Notes
    -----
    Phenotype interpretation is only performed when at least one
    phenotype input is explicitly supplied.

    Therefore:

    phenotype=None
        means phenotype interpretation was not requested.

    phenotype.status == "UNAVAILABLE"
        means phenotype interpretation was requested but the required
        context was unavailable.

    This distinction is intentional and important for future
    case-level reasoning and ranking.
    """

    if not isinstance(
        contract,
        VariantEvidenceContract,
    ):
        raise TypeError(
            "contract must be a VariantEvidenceContract."
        )

    population = interpret_population(
        contract
    )

    clinical = interpret_clinical(
        contract
    )

    functional = interpret_functional(
        contract
    )

    gene_disease = None

    if gene_evidence_context is not None:
        gene_disease = interpret_gene_disease(
            gene_evidence_context
        )

    phenotype = None

    phenotype_requested = (
        present_hpo_terms is not None
        or absent_hpo_terms is not None
    )

    if phenotype_requested:
        phenotype = interpret_phenotype(
            gene_evidence_bundle=(
                gene_evidence_context
            ),
            present_hpo_terms=present_hpo_terms,
            absent_hpo_terms=absent_hpo_terms,
        )

    return InterpretationProfile(
        variant_key=contract.identity.variant_key,
        population=population,
        clinical=clinical,
        functional=functional,
        phenotype=phenotype,
        inheritance=None,
        gene_disease=gene_disease,
    )
