"""Minimal, self-contained plotting helpers so the figure scripts run anywhere."""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

GREY = "#4d4d4d"


def apply_figure_style(sizes: tuple[int, int, int] = (8, 7, 6)) -> None:
    """Role-mapped font ladder, open frame, outward ticks, 300-dpi save."""
    base, secondary, tick = sizes
    mpl.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": base,
            "axes.titlesize": base,
            "axes.labelsize": base,
            "figure.titlesize": base,
            "legend.fontsize": secondary,
            "xtick.labelsize": tick,
            "ytick.labelsize": tick,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titlelocation": "left",
            "axes.titleweight": "regular",
            "xtick.direction": "out",
            "ytick.direction": "out",
            "legend.frameon": False,
            "axes.grid": False,
        }
    )


def panel_letter(ax, letter: str, dx: float = -0.18, dy: float = 1.02) -> None:
    """Bold panel letter outside the top-left of the axes."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontweight="bold",
            fontsize=mpl.rcParams["font.size"] + 1, va="bottom", ha="left")


def strip_with_median(ax, groups, values, colors, order=None, jitter: float = 0.12,
                      seed: int = 0):
    """Jittered points per group with a bold horizontal median tick."""
    groups = np.asarray(groups)
    values = np.asarray(values, dtype=float)
    order = list(colors) if order is None else list(order)
    rng = np.random.default_rng(seed)
    for x, g in enumerate(order):
        v = values[groups == g]
        ax.scatter(x + rng.uniform(-jitter, jitter, v.size), v, s=26,
                   c=colors[g], edgecolors="white", linewidths=0.5, zorder=3)
        ax.plot([x - 0.26, x + 0.26], [np.median(v)] * 2, color=GREY,
                lw=1.8, solid_capstyle="butt", zorder=4)
    ax.set_xlim(-0.6, len(order) - 0.4)
    return ax
