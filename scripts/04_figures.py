#!/usr/bin/env python3
"""Figure 1 (main contrasts) and Figure S1 (supplementary contrasts).

Each row is one set of paired tables compatible with what a report prints:
the marginal counts alone, or the counts together with one reported test
result. Ticks mark the sensitivity value of each compatible table (exact
tail, one-sided alpha = 0.05); the line spans their range. Rows with no
compatible table are marked as such. Output: figures/ (not version-controlled).
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FixedLocator, NullLocator  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

NORD = {"ink": "#2E3440", "muted": "#4C566A", "rule": "#D8DEE9",
        "blue": "#5E81AC", "red": "#BF616A", "frost": "#81A1C1"}
plt.rcParams.update({
    "font.family": "monospace",
    "font.monospace": ["DejaVu Sans Mono", "Menlo", "Courier New"],
    "font.size": 8, "axes.edgecolor": NORD["muted"], "axes.labelcolor": NORD["ink"],
    "xtick.color": NORD["muted"], "ytick.color": NORD["ink"], "text.color": NORD["ink"],
    "svg.fonttype": "none", "pdf.fonttype": 42,
})

MAIN = [
    ("strides_ambulation", "STRIDES: non-ambulatory at 6 months (4/154 vs 20/125)", [
        ("counts only [events in complete pairs]", "counts only, events in complete pairs"),
        ("counts + Table 2 and text: p<0.001 [no location assumption]", "counts + reported p<0.001"),
    ]),
    ("strides_mort6", "STRIDES: death by 6 months (41 vs 59 of 201)", [
        ("counts only", "counts only"),
        ("counts + Table 2: p=0.04", "counts + Table 2 p=0.04"),
        ("counts + text: p<0.001", "counts + text p<0.001"),
        ("counts + all reported", "both reported p-values"),
    ]),
    ("strides_mort12", "STRIDES: death by 1 year (75 vs 102 of 201)", [
        ("counts only", "counts only"),
        ("counts + Table 2: p=0.007", "counts + Table 2 p=0.007"),
        ("counts + text: p=0.11", "counts + text p=0.11"),
    ]),
    ("mohanty_crsvah", "Mohanty 2023: reoperation by 2 years (2 vs 6 of 27)", [
        ("counts only", "counts only"),
        ("counts + Table 7: p=0.0412", "counts + reported p=0.0412"),
    ]),
    ("hamou_recurrence", "Hamou 2026: recurrence (33 vs 41 of 126)", [
        ("counts only", "counts only"),
        ("counts + all reported", "counts + reported χ²=1.02, p=0.312"),
    ]),
]
SUPP = [
    ("mohanty_pill", "Mohanty 2024: reoperation by 2 years (2 vs 6 of 38)", [
        ("counts only", "counts only"), ("counts + Table 6: p=0.2207", "counts + reported p=0.2207")]),
    ("bak_ais", "Bak 2024: no AIS-grade improvement (124 vs 152 of 283)", [
        ("counts only", "counts only"), ("counts + Table 2: p=.024", "counts + reported p=.024")]),
    ("shultz_transfusion", "Shultz 2018: transfusion (1 vs 10 of 255)", [
        ("counts only", "counts only"), ("counts + text: p=0.007", "counts + reported p=0.007")]),
    ("pang_improvement", "Pang 2026: no improvement (49 vs 48 of 167)", [
        ("counts only", "counts only"), ("counts + abstract: exact p=1.000", "counts + reported exact p=1.000")]),
]


def load():
    tabs = {}
    with open(os.path.join(RES, "compatible_tables.csv")) as f:
        for r in csv.DictReader(f):
            tabs.setdefault((r["cid"], r["scenario"]), []).append(float(r["gamma_exact_0.05"]))
    return tabs


def draw(spec, name, xmax=10.0):
    tabs = load()
    rows = []
    for cid, title, scen in spec:
        rows.append(("title", title, None))
        for key, lab in scen:
            rows.append(("row", lab, tabs.get((cid, key), [])))
    h = 0.24 * len(rows) + 0.7
    fig, ax = plt.subplots(figsize=(7.2, h))
    y = len(rows)
    yl = []
    for kind, lab, vals in rows:
        y -= 1
        if kind == "title":
            ax.text(-0.9, y, lab, ha="left", va="center", fontsize=8, fontweight="bold",
                    transform=ax.get_yaxis_transform())
            continue
        ax.text(-0.02, y, lab, ha="right", va="center", fontsize=7.5, color=NORD["muted"],
                transform=ax.get_yaxis_transform())
        if not vals:
            ax.text(1.03, y, "no compatible table", va="center", fontsize=7.5, color=NORD["red"])
            continue
        lo, hi = min(vals), max(vals)
        col = NORD["muted"] if lab.startswith("counts only") else NORD["blue"]
        ax.plot([lo, hi], [y, y], color=col, lw=1.2, solid_capstyle="butt")
        ax.scatter(vals, [y] * len(vals), marker="|", s=40, color=col, lw=0.9, zorder=3)
        txt = f"{lo:.2f}" if hi - lo < 5e-3 else f"{lo:.2f}-{hi:.2f}"
        ax.text(min(hi, xmax) * 1.06, y, txt, va="center", fontsize=7.5, color=col)
    ax.set_xscale("log")
    ax.set_xlim(0.97, xmax * 1.6)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.xaxis.set_major_locator(FixedLocator([1, 1.5, 2, 3, 5, 10]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels(["1", "1.5", "2", "3", "5", "10"])
    ax.set_yticks([])
    for v in (1, 2):
        ax.axvline(v, color=NORD["rule"], lw=0.8, zorder=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_bounds(1, xmax)
    ax.set_xlabel(r"sensitivity value $\hat{\Gamma}$ (log scale)", fontsize=8)
    fig.subplots_adjust(left=0.47, right=0.98, top=0.98, bottom=0.6 / h)
    for ext in ("pdf", "png", "tiff"):
        kw = {"dpi": 600} if ext in ("png", "tiff") else {}
        if ext == "tiff":
            kw["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), **kw)
    plt.close(fig)
    print("wrote", name)


if __name__ == "__main__":
    draw(MAIN, "figure1_compatible_sensitivity")
    draw(SUPP, "figureS1_supplementary_contrasts")
