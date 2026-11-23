#!/usr/bin/env python3
"""Validation of the sensitivity and McNemar routines (results/validation_log.txt).

1. External benchmark: Hammond's smoking data (122 discordant pairs, 110
   favouring the hypothesis); the sensitivity value under the normal
   approximation is 5.59 (Lin, Fogarty and Xu 2026, Section 6.1).
2. From-scratch binomial recomputation agrees with the library (exact tail).
3. McNemar p-values agree with statsmodels for all four conventions.
4. Negative control: a balanced discordant split returns 1.0.
5. Monotonicity of the sensitivity value in the number of favouring pairs.
6. Hamou et al.: the continuity-corrected statistic of the single compatible
   table reproduces the printed chi-square of 1.02 and p of 0.312.
"""
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from pairedbounds import mcnemar_chi2, mcnemar_p, sensitivity_value  # noqa: E402
from statsmodels.stats.contingency_tables import mcnemar as sm_mcnemar  # noqa: E402

OUT = os.path.join(ROOT, "results", "validation_log.txt")


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


def main():
    log, ok = [], True

    def check(name, cond, detail):
        nonlocal ok
        ok = ok and cond
        log.append(f"[{'pass' if cond else 'FAIL'}] {name}: {detail}")

    g_norm = sensitivity_value(122, 110, method="normal")
    g_exact = sensitivity_value(122, 110, method="exact")
    check("Hammond benchmark (normal)", abs(g_norm - 5.59) < 0.01,
          f"{g_norm:.3f} vs published 5.59; exact tail {g_exact:.3f}")
    for D, c in [(24, 20), (16, 16), (11, 10), (100, 70), (8, 6)]:
        a, b = sensitivity_value(D, c), scratch(D, c)
        check(f"scratch D={D} c={c}", abs(a - b) < 1e-3, f"library {a:.4f}, scratch {b:.4f}")
    for b_, c_ in [(20, 28), (1, 10), (3, 19), (31, 49), (0, 4)]:
        tab = [[0, b_], [c_, 0]]
        refs = {"uncorrected": sm_mcnemar(tab, exact=False, correction=False).pvalue,
                "corrected": sm_mcnemar(tab, exact=False, correction=True).pvalue,
                "exact": sm_mcnemar(tab, exact=True).pvalue}
        for conv, ref in refs.items():
            mine = mcnemar_p(b_, c_, conv)
            check(f"McNemar {conv} b={b_} c={c_}", abs(mine - ref) < 1e-10,
                  f"{mine:.6g} vs statsmodels {ref:.6g}")
    check("negative control", sensitivity_value(100, 50) == 1.0, "balanced split gives 1.0")
    vals = [sensitivity_value(40, c) for c in range(20, 41, 4)]
    check("monotonicity", all(y >= x for x, y in zip(vals, vals[1:])),
          ", ".join(f"{v:.2f}" for v in vals))
    x = mcnemar_chi2(20, 28, "corrected")
    check("Hamou reproduction", round(x, 2) == 1.02 and round(mcnemar_p(20, 28, "corrected"), 3) == 0.312,
          f"chi2 {x:.4f}, p {mcnemar_p(20, 28, 'corrected'):.4f}")
    log.append(f"ALL CHECKS PASSED: {ok}")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        f.write("\n".join(log) + "\n")
    print("\n".join(log))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
