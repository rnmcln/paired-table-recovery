# paired-table-recovery

Hidden-bias sensitivity values and internal-consistency checks for 1:1 matched
binary outcomes that are published without the paired 2×2 table.

## Summary

For a 1:1 matched comparison of a binary outcome, McNemar's test and the
Rosenbaum sensitivity analysis depend only on the discordant pairs: *b* (event
in the treated member only) and *c* (event in the control member only). Matched
studies usually report the marginal event counts and a p-value instead. The
marginal counts fix *c* − *b* and bound *b*; a reported McNemar statistic or
p-value, read at its printed precision under each standard test convention,
narrows the set further. This repository enumerates every paired table
compatible with what a report prints and returns

* the range of the conventional sensitivity value Γ̂ (the smallest bias in the
  within-pair odds of treatment at which the one-sided test no longer rejects
  at α = 0.05) across the compatible tables; and
* whether each reported test result can be reproduced by any compatible table
  (an empty set indicates that the printed numbers are mutually inconsistent
  under the conventions considered).

Test conventions: uncorrected and continuity-corrected McNemar χ², exact
conditional and mid-p binomial tests. Sensitivity values use the exact binomial
tail (primary), with mid-p and normal-approximation tails and one-sided
α = 0.025 as sensitivity analyses. Outcomes with incomplete follow-up are
handled either without any assumption about which pairs contain the events or
under the assumption that all reported events fall in complete pairs.

## Worked contrasts (exact tail, one-sided α = 0.05)

| Contrast | Printed information | Compatible tables | Γ̂ |
|---|---|---|---|
| Schoenfeld 2019 (STRIDES), non-ambulatory at 6 months | 4/154 vs 20/125; p < 0.001 | 28 (D 10–24) | 1.93–6.19 |
| STRIDES, death by 6 months | 41 vs 59 of 201; p = 0.04 (table), p < 0.001 (text) | 8 and 7 respectively; none for both | 1.06–1.12 or 1.80–5.52 |
| STRIDES, death by 1 year | 75 vs 102 of 201; p = 0.007 (table), p = 0.11 (text) | 4; none for p = 0.11 | 1.21–1.25 |
| Mohanty 2023, reoperation by 2 years | 2 vs 6 of 27; p = 0.0412 | 3 from counts; none for the p-value | 1.00 |
| Hamou 2026, recurrence | 33 vs 41 of 126; χ² = 1.02, p = 0.312 | 1 (corrected χ², D = 48) | 1.00 |

Supplementary contrasts (Mohanty 2024, Bak 2024, Shultz 2018, Pang 2026) and
all sensitivity analyses are in `results/`.

## Layout

```
src/pairedbounds/    sensitivity values, McNemar conventions, table enumeration
scripts/
  contrasts.py       contrasts with values as printed in each source
  01_screen.py       prespecified PubMed candidate screen (needs NCBI_EMAIL)
  02_analyse.py      compatible tables, sensitivity values, consistency checks
  03_validate.py     benchmark and cross-checks -> results/validation_log.txt
  04_figures.py      figures (written to figures/, not version-controlled)
  05_tables.py       summary table of all scenarios and sensitivity analyses
data/
  sources/extracted_values.csv   every input number with its location in the source
  screen/                        query, records retrieved, screening decisions
results/             CSV outputs and validation log
tests/               unit tests
```

## Reproduce

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
make all     # validate, analyse, tables, figures
make test
```

The analysis is deterministic (no random number generation). Re-running the
candidate screen (`make screen`) requires `NCBI_EMAIL` and optionally
`NCBI_API_KEY` in the environment; the record list from the search of
23 November 2026 is included.

## Validation

`scripts/03_validate.py` reproduces the published sensitivity value of 5.59 for
Hammond's smoking data (122 discordant pairs, 110 favouring the hypothesis,
normal approximation), cross-checks the exact tail against a from-scratch
binomial computation and the McNemar p-values against `statsmodels`, and
reproduces the printed χ² of 1.02 and p of 0.312 for the Hamou contrast.

## Data

All inputs are summary statistics from published articles. No individual
patient data are used.

## Licence and citation

MIT licence. See `CITATION.cff`.
