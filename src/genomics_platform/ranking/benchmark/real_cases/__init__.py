"""Public and de-identified case benchmark contracts.

Real-case ranking inputs are intentionally separated from benchmark
truth so causal labels cannot enter candidate prioritization.
"""

from .contract import (
    CaseProvenance,
    RealCaseInput,
)
from .manifest import (
    load_case_manifest,
    load_truth_manifest,
)
from .loader import load_real_case_input
from .ranking import (
    RealCaseRankingBundle,
    RealCaseRankingRun,
    rank_real_case,
)
from .validation import (
    validate_case_truth_pair,
    validate_real_case_input,
)

__all__ = (
    "CaseProvenance",
    "RealCaseInput",
    "load_case_manifest",
    "load_truth_manifest",
    "load_real_case_input",
    "RealCaseRankingBundle",
    "RealCaseRankingRun",
    "rank_real_case",
    "validate_case_truth_pair",
    "validate_real_case_input",
)
