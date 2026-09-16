"""Controlled Borda tie and coverage stress cases."""

from genomics_platform.ranking.benchmark.controlled_suite import (
    ControlledCase,
    ControlledCaseSpec,
)
from genomics_platform.ranking.benchmark.synthetic_factory import (
    make_synthetic_candidate,
)


def _case(case_id, description, expected_behavior, candidates):
    return ControlledCase(
        spec=ControlledCaseSpec(
            case_id=case_id,
            family="borda_ties",
            description=description,
            causal_candidate_ids=("causal",),
            expected_behavior=expected_behavior,
            tags=("synthetic", "controlled", "stress", "borda", "ties"),
        ),
        candidates=tuple(candidates),
    )


def exact_two_way_tie():
    return _case(
        "borda-exact-two-way-tie-001",
        "Two candidates have intentionally identical evidence.",
        "Borda should preserve a genuine evidence tie.",
        (
            make_synthetic_candidate(
                candidate_id="causal",
                gene="BT1C",
                position=79001,
            ),
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="BT1D",
                position=79002,
            ),
        ),
    )


def exact_three_way_tie():
    return _case(
        "borda-exact-three-way-tie-001",
        "Three candidates have intentionally identical evidence.",
        "Competition ranking should preserve the top tie.",
        (
            make_synthetic_candidate(
                candidate_id="decoy-a",
                gene="BT2A",
                position=79011,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="BT2C",
                position=79012,
            ),
            make_synthetic_candidate(
                candidate_id="decoy-b",
                gene="BT2B",
                position=79013,
            ),
        ),
    )


def sparse_vs_complete():
    return _case(
        "borda-sparse-vs-complete-001",
        "Sparse and complete candidates test evidence coverage correction.",
        "Sparse ballot participation must not gain missingness inflation.",
        (
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="BT3D",
                position=79021,
                af=None,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="BT3C",
                position=79022,
                af=1e-5,
            ),
        ),
    )


def one_signal_sparse():
    return _case(
        "borda-one-signal-sparse-001",
        "Sparse candidate has one strong signal but missing population evidence.",
        "Coverage-aware Borda should expose the completeness penalty.",
        (
            make_synthetic_candidate(
                candidate_id="decoy",
                gene="BT4D",
                position=79031,
                consequence="frameshift_variant",
                impact="HIGH",
                significance="Uncertain_significance",
                af=None,
            ),
            make_synthetic_candidate(
                candidate_id="causal",
                gene="BT4C",
                position=79032,
                consequence="missense_variant",
                impact="MODERATE",
                significance="Uncertain_significance",
                af=1e-5,
            ),
        ),
    )


BUILDERS = (
    exact_two_way_tie,
    exact_three_way_tie,
    sparse_vs_complete,
    one_signal_sparse,
)
