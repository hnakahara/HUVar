"""Shared style for the revised manuscript figures (Genome Medicine revision).

Sequential values (0-1 metrics) use a single blue ramp, light -> dark.
Categorical tool identity uses a fixed-order 6-slot palette (validated for
colour-vision deficiency; three slots are below 3:1 contrast on white, so every
categorical mark also carries a distinct marker shape and a legend).
"""
from __future__ import annotations

import matplotlib
from matplotlib.colors import LinearSegmentedColormap

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update({
    "font.family": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.size": 9,
    "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": "#0b0b0b",
    "xtick.color": "#52514e",
    "ytick.color": "#52514e",
    "axes.titleweight": "bold",
    "axes.titlesize": 10,
})

INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
NA_FILL = "#f0efec"

BLUE_RAMP = ["#e8f1fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
             "#256abf", "#184f95", "#0d366b"]
SEQ = LinearSegmentedColormap.from_list("seq_blue", BLUE_RAMP)

# fixed categorical order (InterVar, BIAS-2015, BIAS-2015 no PP5/BP6,
# HUVar default, HUVar AM+SpliceAI, HUVar REVEL+SpliceAI)
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
MARKERS = ["o", "s", "D", "^", "v", "P"]


def text_color(v: float) -> str:
    """Ink on light cells, white on dark cells."""
    return "white" if v >= 0.62 else INK


def save(fig, stem: str, formats=("png", "pdf")) -> None:
    for f in formats:
        fig.savefig(f"{stem}.{f}", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"wrote {stem}.{{{','.join(formats)}}}")
