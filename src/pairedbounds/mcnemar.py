"""Two-sided McNemar p-values under the conventions met in published reports.

b: discordant pairs with the event in the treated member only.
c: discordant pairs with the event in the control member only.

uncorrected : chi2 = (c - b)^2 / (b + c), 1 df
corrected   : chi2 = (|c - b| - 1)^2 / (b + c), 1 df (Edwards; R mcnemar.test default)
exact       : two-sided exact conditional binomial test, min(1, 2 * tail)
midp        : two-sided mid-p (Fagerland et al. 2013)
"""
from __future__ import annotations

from scipy.stats import binom, chi2

CONVENTIONS = ("uncorrected", "corrected", "exact", "midp")


def mcnemar_chi2(b: int, c: int, convention: str) -> float | None:
    D = b + c
    if D == 0:
        return None
    d = abs(c - b)
    if convention == "uncorrected":
        return d * d / D
    if convention == "corrected":
        return max(d - 1, 0) ** 2 / D
    return None


def mcnemar_p(b: int, c: int, convention: str) -> float:
    D = b + c
    if D == 0:
        return 1.0
    if convention in ("uncorrected", "corrected"):
        return float(chi2.sf(mcnemar_chi2(b, c, convention), 1))
    k = max(b, c)
    if convention == "exact":
        return float(min(1.0, 2.0 * binom.sf(k - 1, D, 0.5)))
    if convention == "midp":
        return float(min(1.0, 2.0 * (binom.sf(k, D, 0.5) + 0.5 * binom.pmf(k, D, 0.5))))
    raise ValueError(f"unknown convention {convention!r}")
