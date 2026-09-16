from genomics_platform.evidence.evidence_bundle import EvidenceBundle
from genomics_platform.evidence.evidence_record import EvidenceRecord


bundle = EvidenceBundle(gene_symbol="TESTGENE")

bundle.add(
    EvidenceRecord(
        source="GenCC",
        evidence_type="gene_disease_validity",
        gene_symbol="TESTGENE",
        gene_id="HGNC:999999",
        disease_id="MONDO:9999999",
        disease_name="Synthetic inflammatory disorder",
        inheritance=["autosomal dominant"],
        classification="Strong",
        source_record_id="SGC-TEST",
        source_version="1",
        license="CC0-1.0",
    )
)

bundle.add(
    EvidenceRecord(
        source="PanelApp",
        evidence_type="expert_gene_panel",
        gene_symbol="TESTGENE",
        disease_id="MONDO:9999999",
        disease_name="Synthetic inflammatory disorder",
        inheritance=["monoallelic"],
        confidence="green",
        source_record_id="TEST-PANEL",
    )
)

bundle.add(
    EvidenceRecord(
        source="HPO",
        evidence_type="gene_phenotype",
        gene_symbol="TESTGENE",
        disease_id="MONDO:9999999",
        phenotype_ids=[
            "HP:0000001",
            "HP:0000002",
        ],
    )
)

data = bundle.to_dict()

summary = data["summary"]

print("Evidence assertions:")
print(f"  records       = {summary['record_count']}")
print(f"  diseases      = {len(summary['disease_ids'])}")
print(f"  phenotypes    = {len(summary['phenotype_ids'])}")
print(f"  sources       = {summary['sources']}")

assert summary["record_count"] == 3
assert summary["sources"]["GenCC"] == 1
assert summary["sources"]["PanelApp"] == 1
assert summary["sources"]["HPO"] == 1
assert len(summary["disease_ids"]) == 1
assert len(summary["phenotype_ids"]) == 2

assert data["scope"]["clinical_classification"] is False
assert data["scope"]["candidate_ranking"] is False

print()
print("==============================================")
print(" PASS: EVIDENCE BUNDLE")
print("==============================================")
