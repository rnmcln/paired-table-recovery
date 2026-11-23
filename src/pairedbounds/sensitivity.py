"""Conventional (Rosenbaum) sensitivity analysis for McNemar's statistic.

For 1:1 matched pairs and a binary outcome, only discordant pairs carry
information. Let D be the number of discordant pairs and c the number in which
only the control member had the adverse event (pairs favouring treatment).
Under Fisher's sharp null and a bias of at most Gamma in the within-pair odds of
treatment, the upper bound on the one-sided p-value is obtained by treating
each discordant pair as favouring treatment independently with probability
rho = Gamma / (1 + Gamma) (Rosenbaum 1987; Rosenbaum 2020, ch. 3). The
sensitivity value Gamma-hat is the smallest Gamma at which this bound exceeds
alpha (Zhao 2019).

Three tail conventions are provided: the exact binomial tail, the mid-p tail
and the large-sample normal approximation. The exact tail is used for the
primary analysis; the other two are sensitivity analyses.
"""
from __future__ import annotations

import math

from scipy.stats import binom, norm

GAMMA_MAX = 25.0


def rho(gamma: float) -> float:
    """Worst-case probability that a discordant pair favours treatment."""
    return gamma / (1.0 + gamma)


def worst_case_p(D: int, c: int, gamma: float, method: str = "exact") -> float:
    """Upper bound on the one-sided p-value at bias gamma.

    D: discordant pairs; c: discordant pairs favouring treatment.
    method: 'exact' | 'midp' | 'normal'.
    """
    D, c = int(D), int(c)
    if D <= 0:
        return 1.0
    p = rho(gamma)
    if method == "exact":
        return float(binom.sf(c - 1, D, p))
    if method == "midp":
        return float(binom.sf(c, D, p) + 0.5 * binom.pmf(c, D, p))
    if method == "normal":
        mean, var = D * p, D * p * (1.0 - p)
        return float(norm.sf((c - mean) / math.sqrt(var)))
    raise ValueError(f"unknown method {method!r}")


def sensitivity_value(D: int, c: int, alpha: float = 0.05, method: str = "exact",
                      gamma_max: float = GAMMA_MAX, tol: float = 1e-6) -> float:
    """Sensitivity value Gamma-hat for a discordant split (D, c).

    Returns 1.0 if the one-sided test does not reject at Gamma = 1, and
    gamma_max if it still rejects at gamma_max. The worst-case p-value is
    non-decreasing in Gamma, so bisection is valid.
    """
    f = lambda g: worst_case_p(D, c, g, method)  # noqa: E731
    if f(1.0) > alpha:
        return 1.0
    if f(gamma_max) <= alpha:
        return gamma_max
    lo, hi = 1.0, gamma_max
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if f(mid) <= alpha:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
