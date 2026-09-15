#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

INPUT="${ROOT}/tests/fixtures/mondo/mondo.synthetic.json"
DB="${ROOT}/tests/fixtures/mondo/mondo.synthetic.sqlite"

rm -f "${DB}"

PYTHONPATH="${ROOT}" python \
  "${ROOT}/src/genomics_platform/evidence/adapters/mondo_importer.py" \
  build \
  --input-json "${INPUT}" \
  --database "${DB}"

python - "${DB}" <<'PY'
import sqlite3
import sys

db = sys.argv[1]

con = sqlite3.connect(db)

try:
    diseases = con.execute(
        "SELECT COUNT(*) FROM disease"
    ).fetchone()[0]

    xrefs = con.execute(
        "SELECT COUNT(*) FROM xref"
    ).fetchone()[0]

    row = con.execute(
        """
        SELECT mondo_id
        FROM xref
        WHERE source_id = 'OMIM:999001'
        """
    ).fetchone()

finally:
    con.close()

print()
print("MONDO assertions:")
print("  diseases =", diseases)
print("  xrefs    =", xrefs)
print("  mapping  =", row[0])

assert diseases == 2
assert xrefs == 3
assert row[0] == "MONDO:0001234"

print()
print("==============================================")
print(" PASS: MONDO IMPORTER")
print("==============================================")
PY

PYTHONPATH="${ROOT}" python \
  "${ROOT}/src/genomics_platform/evidence/adapters/mondo_importer.py" \
  query \
  --database "${DB}" \
  --disease-id OMIM:999001
