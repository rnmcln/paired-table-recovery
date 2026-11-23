"""Enumerate paired 2x2 tables compatible with what a matched study reports.

Notation (1:1 matched pairs, adverse binary event; the treated arm is the one
hypothesised to have fewer events):

                         control: event   control: no event
    treated: event             e                 b
    treated: no event          c                 f

The sensitivity analysis needs only b and c. Published reports usually give
the marginal event counts T = e + b (treated) and C = e + c (control), so that
c - b = C - T, together with a test statistic or p-value.

Base set (marginal counts only)
-------------------------------
Complete outcome data: b = 0..min(T, n - C), c = b + (C - T); the upper limit
keeps e >= 0 and f >= 0.

Incomplete outcome data (marginals counted over patients with known outcome,
the McNemar test over pairs in which both outcomes are known): the events of
the complete pairs are not reported. Two base sets are available:
  * no location assumption: any b <= T, c <= C with b + c not exceeding the
    largest possible number of complete pairs;
  * events-in-complete-pairs assumption: as for complete data, c - b = C - T.

Reported-statistic constraint
-----------------------------
A table is compatible with a reported statistic if, under at least one
declared test convention, its statistic agrees with the printed value at the
printed precision (or satisfies a printed threshold such as p < 0.001).
If no table in the base set is compatible, the reported numbers are mutually
inconsistent under the declared conventions.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .mcnemar import CONVENTIONS, mcnemar_chi2, mcnemar_p


@dataclass(frozen=True)
class Reported:
    label: str                  # e.g. "Table 2 p" or "text p"
    kind: str                   # 'p' | 'chi2'
    relation: str               # '=' | '<' | '>'
    value: float
    decimals: int | None = None
    conventions: tuple = CONVENTIONS

    def holds(self, x: float) -> bool:
        if self.relation == "=":
            half = 0.5 * 10.0 ** (-self.decimals)
            return abs(x - self.value) <= half + 1e-12
        if self.relation == "<":
            return x < self.value
        if self.relation == ">":
            return x > self.value
        raise ValueError(self.relation)


@dataclass
class Contrast:
    cid: str
    study: str
    outcome: str
    n_pairs: int
    treated: str
    control: str
    treated_events: int
    control_events: int
    treated_known: int | None = None
    control_known: int | None = None
    reported: list = field(default_factory=list)
    role: str = "main"

    @property
    def complete(self) -> bool:
        return self.treated_known is None and self.control_known is None


def base_tables(ct: Contrast, assume_events_in_complete_pairs: bool = False):
    """Return the list of (b, c) admitted by the marginal counts."""
    T, C, n = ct.treated_events, ct.control_events, ct.n_pairs
    diff = C - T
    if ct.complete:
        out = []
        for b in range(0, min(T, n - C) + 1):
            c = b + diff
            if 0 <= c <= C:
                out.append((b, c))
        return out
    max_complete = min(ct.treated_known, ct.control_known)
    if assume_events_in_complete_pairs:
        return [(b, b + diff) for b in range(0, T + 1)
                if 0 <= b + diff <= C and 2 * b + diff <= max_complete]
    return [(b, c) for b in range(0, T + 1) for c in range(0, C + 1)
            if b + c <= max_complete]


def matches_reported(b: int, c: int, rep: Reported) -> list:
    """Conventions under which table (b, c) reproduces the reported value."""
    hits = []
    for conv in rep.conventions:
        if rep.kind == "chi2":
            x = mcnemar_chi2(b, c, conv)
            if x is None:
                continue
        else:
            x = mcnemar_p(b, c, conv)
        if rep.holds(x):
            hits.append(conv)
    return hits


def feasible_tables(ct: Contrast, constraints=(), assume_events_in_complete_pairs=False):
    """Tables in the base set compatible with every constraint in `constraints`.

    Returns a list of (b, c, {label: [conventions]}).
    """
    out = []
    for b, c in base_tables(ct, assume_events_in_complete_pairs):
        hits = {}
        ok = True
        for rep in constraints:
            h = matches_reported(b, c, rep)
            if not h:
                ok = False
                break
            hits[rep.label] = h
        if ok:
            out.append((b, c, hits))
    return out


def attainable_range(ct: Contrast, convention: str, assume_events_in_complete_pairs=False):
    """Smallest and largest two-sided McNemar p attainable over the base set."""
    ps = [mcnemar_p(b, c, convention)
          for b, c in base_tables(ct, assume_events_in_complete_pairs) if b + c > 0]
    return (min(ps), max(ps)) if ps else (None, None)
