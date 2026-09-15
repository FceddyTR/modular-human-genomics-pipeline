"""Build disease-corpus HPO information content."""

from __future__ import annotations

import argparse
import json

from genomics_platform.phenotype.hpo_information_content import (
    build_information_content_database,
)


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--associations",
        required=True,
    )

    parser.add_argument(
        "--ontology",
        required=True,
    )

    parser.add_argument(
        "--output",
        required=True,
    )

    args = parser.parse_args()

    result = (
        build_information_content_database(
            association_database=(
                args.associations
            ),
            ontology_database=(
                args.ontology
            ),
            output_database=(
                args.output
            ),
        )
    )

    print(
        json.dumps(
            result,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
