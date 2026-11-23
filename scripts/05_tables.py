#!/usr/bin/env python3
"""Compact result tables (results/table_*.csv) used for reporting."""
import csv
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")


def read(name):
    with open(os.path.join(RES, name)) as f:
        return list(csv.DictReader(f))


def rng(lo, hi, k=2):
    lo, hi = float(lo), float(hi)
    f = lambda v: ">25" if v >= 24.999 else f"{v:.{k}f}"  # noqa: E731
    return f(lo) if abs(lo - hi) < 5e-3 else f"{f(lo)}-{f(hi)}"


def main():
    summ = read("sensitivity_summary.csv")
    rows = []
    for s in summ:
        if s["n_tables"] == "0":
            alt = {k: "none" for k in ("exact_0.05", "midp_0.05", "normal_0.05", "exact_0.025")}
        else:
            alt = {k: rng(s[f"gamma_{k}_min"], s[f"gamma_{k}_max"])
                   for k in ("exact_0.05", "midp_0.05", "normal_0.05", "exact_0.025")}
        rows.append({"cid": s["cid"], "role": s["role"], "study": s["study"], "outcome": s["outcome"],
                     "scenario": s["scenario"], "n_tables": s["n_tables"],
                     "b": s.get("b_range", ""), "D": s.get("D_range", ""),
                     "gamma_exact_alpha0.05": alt["exact_0.05"],
                     "gamma_midp_alpha0.05": alt["midp_0.05"],
                     "gamma_normal_alpha0.05": alt["normal_0.05"],
                     "gamma_exact_alpha0.025": alt["exact_0.025"]})
    with open(os.path.join(RES, "table_sensitivity_all.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote table_sensitivity_all.csv ({len(rows)} rows)")


if __name__ == "__main__":
    main()
