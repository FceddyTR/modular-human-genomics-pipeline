"""Descriptive inheritance interpretation.

Inheritance Interpreter v0.1 combines:
- gene/disease inheritance evidence from a Gene Evidence Bundle,
- optional case-level genotype context.

The interpreter is deliberately conservative.

It does NOT:
- classify variants,
- establish pathogenicity,
- rank candidates,
- diagnose disease,
- infer phase,
- infer compound heterozygosity,
- infer de novo status from missing parental data.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Set

from genomics_platform.interpretation.interpretation_profile import (
    EvidenceDimension,
    InterpretationObservation,
)


ENGINE_VERSION = "0.1.0"


AD = "AUTOSOMAL_DOMINANT"
AR = "AUTOSOMAL_RECESSIVE"
XL = "X_LINKED"
MT = "MITOCHONDRIAL"
MONOALLELIC_AUTOSOMAL = "AUTOSOMAL_DOMINANT_OR_MONOALLELIC"
UNKNOWN = "UNKNOWN"


VALID_GENOTYPES = {
    "heterozygous",
    "homozygous_alt",
    "hemizygous",
    "unknown",
}


def _normalize_genotype(
    genotype_state: Optional[str],
) -> Optional[str]:
    if genotype_state is None:
        return None

    if not isinstance(genotype_state, str):
        raise TypeError(
            "genotype_state must be a string or None."
        )

    value = genotype_state.strip().lower()

    if value not in VALID_GENOTYPES:
        raise ValueError(
            "Unsupported genotype_state: "
            f"{genotype_state!r}"
        )

    return value


def _normalize_inheritance_text(
    value: str,
) -> str:
    text = value.strip().lower()

    if not text:
        return UNKNOWN

    if (
        "autosomal dominant" in text
        or text == "ad"
    ):
        return AD

    if (
        "autosomal recessive" in text
        or text == "ar"
    ):
        return AR

    if "x-linked" in text or "x linked" in text:
        return XL

    if (
        "mitochondrial" in text
        or "maternal inheritance" in text
    ):
        return MT

    if (
        "monoallelic" in text
        and "autosomal" in text
    ):
        return MONOALLELIC_AUTOSOMAL

    return UNKNOWN


def _extract_inheritance_evidence(
    gene_evidence_bundle: Dict[str, Any],
) -> List[Dict[str, Any]]:
    output: List[Dict[str, Any]] = []

    records = gene_evidence_bundle.get(
        "records",
        [],
    )

    if not isinstance(records, list):
        return output

    for record in records:
        if not isinstance(record, dict):
            continue

        evidence = record.get("evidence")

        if not isinstance(evidence, dict):
            continue

        raw_values: List[str] = []

        inheritance = evidence.get(
            "inheritance",
            [],
        )

        if isinstance(inheritance, list):
            raw_values.extend(
                value
                for value in inheritance
                if isinstance(value, str)
            )

        phenotype_ids = evidence.get(
            "phenotype_ids",
            [],
        )

        description = evidence.get(
            "evidence_description"
        )

        # HP:0000006 is Autosomal dominant inheritance.
        if (
            isinstance(phenotype_ids, list)
            and "HP:0000006" in phenotype_ids
        ):
            raw_values.append(
                "Autosomal dominant inheritance"
            )

        for raw_value in raw_values:
            output.append(
                {
                    "source": evidence.get(
                        "source"
                    ),
                    "gene_symbol": evidence.get(
                        "gene_symbol"
                    ),
                    "disease_id": evidence.get(
                        "disease_id"
                    ),
                    "raw_inheritance": raw_value,
                    "normalized_inheritance": (
                        _normalize_inheritance_text(
                            raw_value
                        )
                    ),
                    "source_record_id": (
                        evidence.get(
                            "source_record_id"
                        )
                    ),
                    "evidence_description": (
                        description
                    ),
                    "disease_normalization": (
                        record.get(
                            "disease_normalization"
                        )
                    ),
                }
            )

    return output


def _compatibility_for_model(
    model: str,
    genotype_state: Optional[str],
) -> str:
    """Return descriptive genotype/model compatibility."""

    if genotype_state is None:
        return "UNKNOWN"

    if genotype_state == "unknown":
        return "UNKNOWN"

    if model in {
        AD,
        MONOALLELIC_AUTOSOMAL,
    }:
        if genotype_state == "heterozygous":
            return "COMPATIBLE"

        if genotype_state == "homozygous_alt":
            return "PARTIAL"

        return "UNKNOWN"

    if model == AR:
        if genotype_state == "homozygous_alt":
            return "COMPATIBLE"

        if genotype_state == "heterozygous":
            # A single heterozygous small variant is not enough
            # to establish a recessive genotype. A second allele,
            # phase, CNV/SV, or other event may still exist.
            return "PARTIAL"

        return "UNKNOWN"

    if model == XL:
        if genotype_state in {
            "heterozygous",
            "hemizygous",
        }:
            return "PARTIAL"

        return "UNKNOWN"

    if model == MT:
        return "UNKNOWN"

    return "UNKNOWN"


def interpret_inheritance(
    gene_evidence_bundle: Optional[
        Dict[str, Any]
    ],
    genotype_state: Optional[str] = None,
    proband_sex: Optional[str] = None,
    de_novo_status: Optional[bool] = None,
) -> EvidenceDimension:
    """Interpret inheritance evidence descriptively."""

    genotype = _normalize_genotype(
        genotype_state
    )

    if proband_sex is not None:
        if not isinstance(proband_sex, str):
            raise TypeError(
                "proband_sex must be a string or None."
            )

        proband_sex = (
            proband_sex.strip().lower()
        )

        if proband_sex not in {
            "male",
            "female",
            "unknown",
        }:
            raise ValueError(
                "proband_sex must be male, female, "
                "unknown, or None."
            )

    if (
        de_novo_status is not None
        and not isinstance(
            de_novo_status,
            bool,
        )
    ):
        raise TypeError(
            "de_novo_status must be bool or None."
        )

    if gene_evidence_bundle is None:
        return EvidenceDimension(
            name="inheritance",
            status="UNAVAILABLE",
            observations=[],
            limitations=[
                "No gene evidence context was supplied.",
                "Inheritance compatibility was not evaluated.",
                "No candidate ranking or variant classification was performed.",
            ],
        )

    if not isinstance(
        gene_evidence_bundle,
        dict,
    ):
        raise TypeError(
            "gene_evidence_bundle must be a dictionary or None."
        )

    evidence_records = (
        _extract_inheritance_evidence(
            gene_evidence_bundle
        )
    )

    observations: List[
        InterpretationObservation
    ] = []

    observations.append(
        InterpretationObservation(
            code="INHERITANCE_CASE_CONTEXT",
            label="Case inheritance context",
            value={
                "genotype_state": genotype,
                "proband_sex": proband_sex,
                "de_novo_status": de_novo_status,
            },
            source="patient_context",
            details={
                "phase_inferred": False,
                "compound_heterozygosity_inferred": False,
                "candidate_ranking": False,
                "variant_classification": False,
            },
        )
    )

    normalized_models: Set[str] = set()

    compatibility_values: List[str] = []

    for item in evidence_records:
        model = item[
            "normalized_inheritance"
        ]

        normalized_models.add(model)

        compatibility = (
            _compatibility_for_model(
                model,
                genotype,
            )
        )

        compatibility_values.append(
            compatibility
        )

        observations.append(
            InterpretationObservation(
                code="INHERITANCE_SOURCE_EVIDENCE",
                label="Inheritance source evidence",
                value={
                    "source": item["source"],
                    "gene_symbol": (
                        item["gene_symbol"]
                    ),
                    "disease_id": (
                        item["disease_id"]
                    ),
                    "raw_inheritance": (
                        item[
                            "raw_inheritance"
                        ]
                    ),
                    "normalized_inheritance": (
                        model
                    ),
                    "genotype_compatibility": (
                        compatibility
                    ),
                },
                source=item["source"],
                details={
                    "source_record_id": (
                        item[
                            "source_record_id"
                        ]
                    ),
                    "evidence_description": (
                        item[
                            "evidence_description"
                        ]
                    ),
                    "disease_normalization": (
                        item[
                            "disease_normalization"
                        ]
                    ),
                    "variant_pathogenicity": False,
                    "candidate_ranking": False,
                    "diagnosis": False,
                },
            )
        )

    observations.append(
        InterpretationObservation(
            code="INHERITANCE_SUMMARY",
            label="Inheritance interpretation summary",
            value={
                "normalized_models": sorted(
                    normalized_models
                ),
                "source_evidence_count": len(
                    evidence_records
                ),
                "compatibility_states": sorted(
                    set(compatibility_values)
                ),
            },
            source="inheritance_interpreter",
            details={
                "candidate_ranking": False,
                "variant_classification": False,
                "pathogenicity_probability": False,
            },
        )
    )

    if not evidence_records:
        status = "UNAVAILABLE"
    elif genotype is None or genotype == "unknown":
        status = "PARTIAL"
    elif "COMPATIBLE" in compatibility_values:
        status = "AVAILABLE"
    else:
        status = "PARTIAL"

    limitations = [
        "Inheritance interpretation is descriptive and does not establish variant pathogenicity.",
        "Genotype compatibility does not establish causality or diagnosis.",
        "An incompatible or partial observation must not automatically exclude a candidate.",
        "Phase is not inferred.",
        "Compound heterozygosity is not inferred from a single variant.",
        "Copy-number, structural, repeat, mosaic, and complex alleles are not evaluated by this v0.1 interpreter.",
        "X-linked interpretation requires richer sex-aware and locus-aware logic in future versions.",
        "Mitochondrial interpretation requires heteroplasmy and maternal-lineage context in future versions.",
        "De novo status is preserved when supplied but is not independently verified here.",
        "PanelApp monoallelic terminology is preserved conservatively and is not automatically treated as definitive autosomal-dominant disease inheritance.",
        "No ACMG/AMP classification, pathogenicity probability, candidate tier, or ranking is produced.",
    ]

    return EvidenceDimension(
        name="inheritance",
        status=status,
        observations=observations,
        limitations=limitations,
    )
