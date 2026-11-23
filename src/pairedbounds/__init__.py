"""pairedbounds: hidden-bias sensitivity and internal-consistency checks for
1:1 matched binary outcomes reported without the paired 2x2 table."""

from .sensitivity import rho, worst_case_p, sensitivity_value
from .mcnemar import mcnemar_p, mcnemar_chi2, CONVENTIONS
from .recovery import (
    Reported, Contrast, base_tables, matches_reported, feasible_tables,
    attainable_range,
)

__all__ = [
    "rho", "worst_case_p", "sensitivity_value",
    "mcnemar_p", "mcnemar_chi2", "CONVENTIONS",
    "Reported", "Contrast", "base_tables", "matches_reported",
    "feasible_tables", "attainable_range",
]
__version__ = "1.0.0"
