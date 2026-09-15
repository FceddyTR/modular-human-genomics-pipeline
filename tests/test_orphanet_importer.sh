#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

cat > "$TMP/orphanet.xml" <<'XML'
<?xml version="1.0" encoding="UTF-8"?>
<JDBOR date="2026-06-23 07:57:31">
  <Availability>
    <Licence>
      <ShortIdentifier>CC-BY-4.0</ShortIdentifier>
    </Licence>
  </Availability>

  <DisorderList count="1">
    <Disorder id="17601">
      <OrphaCode>166024</OrphaCode>
      <Name lang="en">Synthetic rare disease</Name>

      <DisorderGeneAssociationList count="1">
        <DisorderGeneAssociation>
          <SourceOfValidation>22587682[PMID]</SourceOfValidation>

          <Gene id="20160">
            <Name lang="en">kinesin family member 7</Name>
            <Symbol>KIF7</Symbol>

            <GeneType id="25993">
              <Name lang="en">gene with protein product</Name>
            </GeneType>

            <ExternalReferenceList count="3">
              <ExternalReference>
                <Source>Ensembl</Source>
                <Reference>ENSG00000166813</Reference>
              </ExternalReference>

              <ExternalReference>
                <Source>HGNC</Source>
                <Reference>30497</Reference>
              </ExternalReference>

              <ExternalReference>
                <Source>OMIM</Source>
                <Reference>611254</Reference>
              </ExternalReference>
            </ExternalReferenceList>
          </Gene>

          <DisorderGeneAssociationType>
            <Name lang="en">Disease-causing germline mutation(s) in</Name>
          </DisorderGeneAssociationType>

          <DisorderGeneAssociationStatus>
            <Name lang="en">Assessed</Name>
          </DisorderGeneAssociationStatus>
        </DisorderGeneAssociation>
      </DisorderGeneAssociationList>
    </Disorder>
  </DisorderList>
</JDBOR>
XML

python -m src.genomics_platform.evidence.adapters.orphanet_importer \
  --input-xml "$TMP/orphanet.xml" \
  --output-json "$TMP/orphanet.json"

python -m src.genomics_platform.evidence.orphanet_store build \
  --input-xml "$TMP/orphanet.xml" \
  --database "$TMP/orphanet.sqlite"

python - "$TMP/orphanet.json" "$TMP/orphanet.sqlite" <<'PY'
import json
import sys

from src.genomics_platform.evidence.orphanet_store import (
    query_disease,
    query_gene,
)

json_path = sys.argv[1]
database = sys.argv[2]

with open(json_path, encoding="utf-8") as handle:
    data = json.load(handle)

assert data["summary"]["record_count"] == 1
assert data["summary"]["unique_genes"] == 1
assert data["summary"]["unique_diseases"] == 1

record = data["records"][0]

assert record["source"] == "Orphanet"
assert record["evidence_type"] == "gene_disease"
assert record["gene_symbol"] == "KIF7"
assert record["gene_id"] == "HGNC:30497"
assert record["disease_id"] == "ORPHA:166024"
assert record["classification"] is None
assert record["confidence"] == "Assessed"
assert record["citations"] == ["PMID:22587682"]
assert record["license"] == "CC-BY-4.0"
assert record["source_release"] == "2026-06-23 07:57:31"

gene_records = query_gene(
    database,
    "kif7",
)

assert len(gene_records) == 1
assert gene_records[0]["disease_id"] == "ORPHA:166024"

disease_records = query_disease(
    database,
    "166024",
)

assert len(disease_records) == 1
assert disease_records[0]["gene_symbol"] == "KIF7"

print("Orphanet synthetic regression: PASS")
PY
