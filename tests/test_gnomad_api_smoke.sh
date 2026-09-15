#!/usr/bin/env bash
set -euo pipefail

export PYTHONPATH="${PYTHONPATH:-}:$(pwd)/src"

python - <<'PY'
from genomics_platform.evidence.adapters.gnomad_api_adapter import query_variant

r = query_variant("1", 55051215, "G", "GA")

assert r.lookup_status == "FOUND", r.to_dict()
assert r.exome is not None
assert r.genome is not None

assert r.exome.ac is not None
assert r.exome.an is not None
assert r.exome.af is not None

assert r.genome.ac is not None
assert r.genome.an is not None
assert r.genome.af is not None

assert len(r.exome.populations) > 0
assert len(r.genome.populations) > 0

print("PASS: GNOMAD API SMOKE v0.1")
print("variant:", r.variant_key)
print("exome AF:", r.exome.af)
print("genome AF:", r.genome.af)
print("exome populations:", len(r.exome.populations))
print("genome populations:", len(r.genome.populations))
PY
