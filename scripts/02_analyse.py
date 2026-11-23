#!/usr/bin/env python3
"""Compatible paired tables, sensitivity values and consistency checks.

Outputs (results/):
  compatible_tables.csv  one row per contrast x scenario x compatible table
  sensitivity_summary.csv one row per contrast x scenario (ranges over tables)
  consistency.csv        attainability of each reported statistic
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
from contrasts import CONTRASTS  # noqa: E402
from pairedbounds import (CONVENTIONS, attainable_range, feasible_tables,  # noqa: E402
                          mcnemar_p, sensitivity_value)
import numpy as np  # noqa: E402
from scipy.stats import chi2_contingency, fisher_exact  # noqa: E402


def unpaired_p(ct):
    """Two-sided p-values of unpaired tests on the marginal counts (diagnostic only)."""
    n1 = ct.treated_known or ct.n_pairs
    n2 = ct.control_known or ct.n_pairs
    t = np.array([[ct.treated_events, n1 - ct.treated_events],
                  [ct.control_events, n2 - ct.control_events]])
    return {"pearson": float(chi2_contingency(t, correction=False)[1]),
            "yates": float(chi2_contingency(t, correction=True)[1]),
            "fisher": float(fisher_exact(t)[1])}

RES = os.path.join(ROOT, "results")
os.makedirs(RES, exist_ok=True)

SV_SPECS = [("exact", 0.05), ("midp", 0.05), ("normal", 0.05), ("exact", 0.025)]


def scenarios(ct):
    """(scenario label, constraints, location assumption) for one contrast."""
    out = []
    assumptions = [False] if ct.complete else [False, True]
    for a in assumptions:
        tag = "" if ct.complete else (" [events in complete pairs]" if a else " [no location assumption]")
        out.append(("counts only" + tag, (), a))
        for r in ct.reported:
            out.append((f"counts + {r.label}" + tag, (r,), a))
        if len(ct.reported) > 1:
            out.append(("counts + all reported" + tag, tuple(ct.reported), a))
    return out


def fmt(v):
    return ">25" if v >= 25.0 - 1e-9 else f"{v:.2f}"


def main():
    rows_t, rows_s, rows_c = [], [], []
    for ct in CONTRASTS:
        for label, cons, assume in scenarios(ct):
            tabs = feasible_tables(ct, cons, assume)
            svs = {}
            for b, c, hits in tabs:
                rec = {"cid": ct.cid, "study": ct.study, "outcome": ct.outcome, "role": ct.role,
                       "scenario": label, "b": b, "c": c, "D": b + c,
                       "matched_conventions": "; ".join(f"{k}: {','.join(v)}" for k, v in hits.items())}
                for m, a in SV_SPECS:
                    v = sensitivity_value(b + c, c, alpha=a, method=m)
                    rec[f"gamma_{m}_{a}"] = round(v, 4)
                    svs.setdefault((m, a), []).append(v)
                for conv in CONVENTIONS:
                    rec[f"p_{conv}"] = round(mcnemar_p(b, c, conv), 6)
                rows_t.append(rec)
            s = {"cid": ct.cid, "study": ct.study, "outcome": ct.outcome, "role": ct.role,
                 "n_pairs": ct.n_pairs, "treated_events": ct.treated_events,
                 "control_events": ct.control_events, "scenario": label,
                 "n_tables": len(tabs)}
            if tabs:
                bs = [t[0] for t in tabs]; cs = [t[1] for t in tabs]; ds = [t[0] + t[1] for t in tabs]
                s.update(b_range=f"{min(bs)}-{max(bs)}", c_range=f"{min(cs)}-{max(cs)}",
                         D_range=f"{min(ds)}-{max(ds)}")
                for (m, a), vals in svs.items():
                    s[f"gamma_{m}_{a}_min"] = round(min(vals), 4)
                    s[f"gamma_{m}_{a}_max"] = round(max(vals), 4)
                lo, hi = min(svs[("exact", 0.05)]), max(svs[("exact", 0.05)])
                s["gamma_exact_0.05"] = fmt(lo) if abs(lo - hi) < 5e-3 else f"{fmt(lo)}-{fmt(hi)}"
            else:
                s["gamma_exact_0.05"] = "no compatible table"
            rows_s.append(s)
        # consistency: is each reported statistic attainable, and under which convention?
        assumes = [False] if ct.complete else [False, True]
        for r in ct.reported:
            for a in assumes:
                tabs = feasible_tables(ct, (r,), a)
                convs = sorted({cv for _, _, h in tabs for cv in h[r.label]})
                rng = {cv: attainable_range(ct, cv, a) for cv in CONVENTIONS}
                up = unpaired_p(ct)
                up_hits = [k for k, v in up.items() if r.kind == "p" and r.holds(v)]
                rows_c.append({
                    "cid": ct.cid, "study": ct.study, "outcome": ct.outcome,
                    "reported": r.label,
                    "location_assumption": "" if ct.complete else ("events in complete pairs" if a else "none"),
                    "attainable": bool(tabs), "n_tables": len(tabs),
                    "conventions_reproducing": ",".join(convs),
                    **{f"p_range_{cv}": f"{rng[cv][0]:.4g}-{rng[cv][1]:.4g}" for cv in CONVENTIONS},
                    **{f"unpaired_p_{k}": round(v, 4) for k, v in up.items()},
                    "unpaired_tests_reproducing": ",".join(up_hits),
                })

    def write(name, rows):
        cols = []
        for r in rows:
            for k in r:
                if k not in cols:
                    cols.append(k)
        with open(os.path.join(RES, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)

    write("compatible_tables.csv", rows_t)
    write("sensitivity_summary.csv", rows_s)
    write("consistency.csv", rows_c)
    for s in rows_s:
        print(f"{s['cid']:20s} {s['scenario'][:62]:62s} n={s['n_tables']:4d} "
              f"D={s.get('D_range', ''):8s} G={s['gamma_exact_0.05']}")
    print()
    for c in rows_c:
        print(f"{c['cid']:20s} {c['reported']:28s} {c['location_assumption']:24s} "
              f"attainable={c['attainable']} conv={c['conventions_reproducing']} unpaired={c['unpaired_tests_reproducing']}")


if __name__ == "__main__":
    main()
