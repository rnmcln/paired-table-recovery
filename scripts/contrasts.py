"""Matched binary contrasts analysed, with the values exactly as printed.

Every number below is transcribed from the source article; the location is
given in data/sources/extracted_values.csv. The treated arm is the arm
hypothesised to have fewer adverse events.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from pairedbounds import Contrast, Reported  # noqa: E402

CONTRASTS = [
    # Schoenfeld et al. 2019 (STRIDES), Cancer 125:2631-2637, Table 2.
    # Ambulatory status known for 154 operative and 125 non-operative patients.
    Contrast("strides_ambulation", "Schoenfeld 2019 (STRIDES)", "Non-ambulatory at 6 months",
             201, "operative", "non-operative", 4, 20, treated_known=154, control_known=125,
             reported=[Reported("Table 2 and text: p<0.001", "p", "<", 0.001)]),
    Contrast("strides_mort6", "Schoenfeld 2019 (STRIDES)", "Death by 6 months",
             201, "operative", "non-operative", 41, 59,
             reported=[Reported("Table 2: p=0.04", "p", "=", 0.04, 2),
                       Reported("text: p<0.001", "p", "<", 0.001)]),
    Contrast("strides_mort12", "Schoenfeld 2019 (STRIDES)", "Death by 1 year",
             201, "operative", "non-operative", 75, 102,
             reported=[Reported("Table 2: p=0.007", "p", "=", 0.007, 3),
                       Reported("text: p=0.11", "p", "=", 0.11, 2)]),
    # Mohanty et al. 2023, J Neurosurg Spine 39:301-310, Table 7 (matched columns).
    Contrast("mohanty_crsvah", "Mohanty 2023 (CrSVA-H)", "Reoperation by 2 years",
             27, "aligned (CrSVA-H <2 cm)", "malaligned", 2, 6,
             reported=[Reported("Table 7: p=0.0412", "p", "=", 0.0412, 4)]),
    # Hamou et al. 2026, Brain Spine 6:106048, Results (matched cohort).
    Contrast("hamou_recurrence", "Hamou 2026", "Recurrence requiring reoperation",
             126, "pre-existing statin", "no statin", 33, 41,
             reported=[Reported("text: chi2=1.02", "chi2", "=", 1.02, 2,
                                ("uncorrected", "corrected")),
                       Reported("text: p=0.312", "p", "=", 0.312, 3)]),
    # Diagnostic only: the printed percentages (3%, 16%) applied to all 201 pairs, which ignores
    # the denominators of patients with known status printed in Table 2. Not used for inference.
    Contrast("strides_ambulation_pct201", "Schoenfeld 2019 (STRIDES)",
             "Non-ambulatory at 6 months, percentages applied to 201 pairs (diagnostic)",
             201, "operative", "non-operative", 6, 32, role="diagnostic"),
    # Supplementary contrasts from the candidate screen.
    Contrast("mohanty_pill", "Mohanty 2023 (PI-LL)", "Reoperation by 2 years",
             38, "PI-LL <10", "PI-LL >=10", 2, 6,
             reported=[Reported("Table 6: p=0.2207", "p", "=", 0.2207, 4)], role="supplement"),
    Contrast("bak_ais", "Bak 2024", "No AIS-grade improvement",
             283, "direct admission", "transfer", 124, 152,
             reported=[Reported("Table 2: p=.024", "p", "=", 0.024, 3)], role="supplement"),
    Contrast("shultz_transfusion", "Shultz 2018", "Blood transfusion",
             255, "TDA", "ALIF", 1, 10,
             reported=[Reported("text: p=0.007", "p", "=", 0.007, 3)], role="supplement"),
    Contrast("pang_improvement", "Pang 2026", "No neurological improvement",
             167, "glucocorticoid", "no glucocorticoid", 49, 48,
             reported=[Reported("abstract: exact p=1.000", "p", "=", 1.0, 3, ("exact",))],
             role="supplement"),
]
