import math
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from pairedbounds import (Contrast, Reported, base_tables, feasible_tables,  # noqa: E402
                          mcnemar_p, rho, sensitivity_value, worst_case_p)


def scratch(D, c, alpha=0.05, gmax=25.0):
    def p(g):
        r = g / (1 + g)
        return sum(math.comb(D, k) * r ** k * (1 - r) ** (D - k) for k in range(c, D + 1))
    if p(1.0) > alpha:
        return 1.0
    lo, hi = 1.0, gmax
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if p(mid) <= alpha else (lo, mid)
    return (lo + hi) / 2


@pytest.mark.parametrize("D,c", [(26, 26), (24, 20), (38, 32), (11, 10), (122, 110)])
def test_sensitivity_matches_scratch(D, c):
    assert abs(sensitivity_value(D, c) - scratch(D, c)) < 1e-3


def test_closed_form_all_favouring():
    # all D pairs favour treatment: rho^D = alpha  =>  Gamma = r / (1 - r)
    r = 0.05 ** (1 / 26)
    assert abs(sensitivity_value(26, 26) - r / (1 - r)) < 1e-4


def test_hammond_benchmark():
    assert abs(sensitivity_value(122, 110, method="normal") - 5.59) < 0.01


@pytest.mark.parametrize("g", [1.0, 2.0, 7.5])
def test_rho(g):
    assert abs(rho(g) - g / (1 + g)) < 1e-12


def test_worst_case_p_monotone():
    ps = [worst_case_p(30, 22, g) for g in (1, 1.5, 2, 3, 5)]
    assert all(b >= a for a, b in zip(ps, ps[1:]))


def test_null_returns_one():
    assert sensitivity_value(40, 20) == 1.0
    assert sensitivity_value(0, 0) == 1.0


def test_exact_more_conservative_than_midp():
    for D, c in [(8, 7), (15, 12), (24, 20)]:
        assert sensitivity_value(D, c, method="exact") <= sensitivity_value(D, c, method="midp") + 1e-9


def test_mcnemar_conventions():
    assert abs(mcnemar_p(0, 4, "uncorrected") - 0.0455003) < 1e-6
    assert abs(mcnemar_p(0, 4, "corrected") - 0.1336144) < 1e-6
    assert abs(mcnemar_p(0, 4, "exact") - 0.125) < 1e-12
    assert mcnemar_p(5, 5, "exact") == 1.0


def test_base_tables_complete_identity():
    ct = Contrast("x", "s", "o", 50, "t", "c", 3, 9)
    tabs = base_tables(ct)
    assert tabs == [(b, b + 6) for b in range(4)]


def test_base_tables_respects_concordant_cap():
    # n - C limits b: with n=10, T=8, C=9 we need f >= 0 => b <= 1
    ct = Contrast("x", "s", "o", 10, "t", "c", 8, 9)
    assert max(b for b, _ in base_tables(ct)) == 1


def test_incomplete_with_and_without_location_assumption():
    ct = Contrast("x", "s", "o", 20, "t", "c", 1, 4, treated_known=15, control_known=12)
    assume = base_tables(ct, True)
    free = base_tables(ct, False)
    assert assume == [(0, 3), (1, 4)]
    assert set(assume) <= set(free)
    assert (0, 0) in free


def test_reported_rounding_and_thresholds():
    r = Reported("p", "p", "=", 0.04, 2)
    assert r.holds(0.0449) and r.holds(0.035) and not r.holds(0.0451)
    assert Reported("p", "p", "<", 0.001).holds(0.0009)


def test_hamou_singleton():
    ct = Contrast("h", "Hamou", "rec", 126, "statin", "none", 33, 41)
    rep = Reported("chi2", "chi2", "=", 1.02, 2, ("uncorrected", "corrected"))
    tabs = feasible_tables(ct, (rep,))
    assert [(b, c) for b, c, _ in tabs] == [(20, 28)]
    assert tabs[0][2]["chi2"] == ["corrected"]


def test_unattainable_p_detected():
    ct = Contrast("m", "Mohanty", "reop", 27, "aligned", "malaligned", 2, 6)
    rep = Reported("p", "p", "=", 0.0412, 4)
    assert feasible_tables(ct, (rep,)) == []
