"""Build the local HPO ontology database."""

from __future__ import annotations

import argparse
import json

from genomics_platform.phenotype.hpo_ontology_store import (
    build_ontology_database,
)


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--obo",
        required=True,
    )

    parser.add_argument(
        "--database",
        required=True,
    )

    args = parser.parse_args()

    result = build_ontology_database(
        obo_path=args.obo,
        database_path=args.database,
    )

    print(
        json.dumps(
            result,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
