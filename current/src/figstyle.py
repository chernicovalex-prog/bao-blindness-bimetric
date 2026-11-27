"""Shared style and helpers for the paper figures (src/make_figures.py).

Computer-Modern look without requiring a LaTeX installation (mathtext 'cm').
Colours: first slots of the validated categorical palette; identity is never
carried by colour alone (direct labels, markers, line styles).
"""
import pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
FIGDIR = ROOT / "figures"
RESDIR = ROOT / "results"

INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#dcdbd5"
BLUE = "#2a78d6"     # BAO
ORANGE = "#eb6834"   # SNe
AQUA = "#169a6b"     # CMB (slightly darker aqua for white paper)
VIOLET = "#6250d6"   # joint
RED = "#c8342f"
GREY = "#8a8983"

PROBE = {"BAO": BLUE, "SNe": ORANGE, "CMB": AQUA, "joint": VIOLET}


def apply_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["cmr10", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.formatter.use_mathtext": True,
        "axes.unicode_minus": False,
        "font.size": 8.5,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "legend.fontsize": 7.5,
        "legend.frameon": False,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.edgecolor": INK2,
        "axes.labelcolor": INK,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "text.color": INK,
        "axes.linewidth": 0.7,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 3.2,
        "ytick.major.size": 3.2,
        "xtick.minor.size": 1.8,
        "ytick.minor.size": 1.8,
        "lines.linewidth": 1.6,
        "lines.markersize": 4.5,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "grid.alpha": 0.8,
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def panel_label(ax, s, x=0.02, y=0.96, **kw):
    ax.text(x, y, s, transform=ax.transAxes, ha="left", va="top",
            fontsize=9, fontweight="bold", **kw)


def save(fig, name):
    FIGDIR.mkdir(exist_ok=True)
    fig.savefig(FIGDIR / (name + ".pdf"))
    fig.savefig(FIGDIR / (name + ".png"), dpi=200)
    plt.close(fig)


class Checks:
    """Numbers annotated on / read off the figures must equal the paper's."""

    def __init__(self):
        self.lines = []
        self.failed = 0

    def close(self, label, got, want, tol, rel=False):
        err = abs(got - want) / (abs(want) if rel else 1.0)
        ok = err <= tol
        self.lines.append("%s  %s: got %.6g, paper %.6g (%s err %.2e <= %.1e)"
                          % ("PASS" if ok else "FAIL", label, got, want,
                             "rel" if rel else "abs", err, tol))
        self.failed += (not ok)

    def true(self, label, cond):
        self.lines.append("%s  %s" % ("PASS" if cond else "FAIL", label))
        self.failed += (not cond)

    def write(self, path):
        txt = "\n".join(self.lines) + "\n" + (
            "ALL FIGURE CHECKS PASSED (%d)\n" % len(self.lines) if not self.failed
            else "%d FIGURE CHECKS FAILED\n" % self.failed)
        pathlib.Path(path).write_text(txt)
        print(txt)
        return self.failed
