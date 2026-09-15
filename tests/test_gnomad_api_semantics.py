from genomics_platform.evidence.adapters.gnomad_api_adapter import query_variant


def test_api_error_is_not_not_found():
    record = query_variant(
        "1",
        55051215,
        "G",
        "GA",
        api_url="http://127.0.0.1:1/api",
        timeout=1,
    )

    assert record.lookup_status == "ERROR"
    assert record.error is not None
    assert record.exome is None
    assert record.genome is None


def test_not_found_record_has_no_frequency_evidence():
    # Synthetic semantic assertion:
    # NOT_FOUND must never be represented as AF=0.
    from genomics_platform.evidence.gnomad_record import GnomADRecord

    record = GnomADRecord(
        assembly="GRCh38",
        chrom="1",
        pos=1,
        ref="A",
        alt="T",
        lookup_status="NOT_FOUND",
    )

    data = record.to_dict()

    assert data["lookup_status"] == "NOT_FOUND"
    assert data["exome"] is None
    assert data["genome"] is None
    assert "does not imply AF=0" in data["semantics"]["not_found"]


def test_graphql_variant_not_found_is_not_found(monkeypatch):
    import io
    import json

    from genomics_platform.evidence.adapters import gnomad_api_adapter

    payload = {
        "errors": [
            {
                "message": "Variant not found"
            }
        ]
    }

    class FakeResponse:
        def __enter__(self):
            return io.BytesIO(json.dumps(payload).encode())

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr(
        gnomad_api_adapter.urllib.request,
        "urlopen",
        lambda *args, **kwargs: FakeResponse(),
    )

    record = gnomad_api_adapter.query_variant(
        "1",
        66926,
        "AG",
        "A",
    )

    assert record.lookup_status == "NOT_FOUND"
    assert record.error is None
    assert record.exome is None
    assert record.genome is None


def test_graphql_service_overloaded_remains_error(monkeypatch):
    import io
    import json

    from genomics_platform.evidence.adapters import gnomad_api_adapter

    payload = {
        "errors": [
            {
                "message": "Service overloaded"
            }
        ]
    }

    class FakeResponse:
        def __enter__(self):
            return io.BytesIO(json.dumps(payload).encode())

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr(
        gnomad_api_adapter.urllib.request,
        "urlopen",
        lambda *args, **kwargs: FakeResponse(),
    )

    record = gnomad_api_adapter.query_variant(
        "1",
        66926,
        "AG",
        "A",
    )

    assert record.lookup_status == "ERROR"
    assert record.error == "GraphQL: Service overloaded"
    assert record.exome is None
    assert record.genome is None
