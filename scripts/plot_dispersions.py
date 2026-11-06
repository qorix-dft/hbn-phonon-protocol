#!/usr/bin/env python3
"""
Panels (a)-(d): pristine phonon dispersions, DFT against the MLFF, each with
its density of states.

    python plot_dispersions.py

Reads band.yaml and total_dos.dat from the branch tree. Everything
configurable is in the SETTINGS block.

On styling: seaborn is used for the palette and context scaling only. Its
figure styles (despined axes, background grids, sans-serif defaults) are built
for statistical graphics; a band structure wants boxed axes with inward ticks
and no grid, which is what journals expect. Set STYLE = "seaborn" to see the
alternative.
"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

# ----------------------------------------------------------------------
# SETTINGS
# ----------------------------------------------------------------------
PANELS = [("bulk",      "4x4x2", "Bulk"),
          ("monolayer", "8x8x1", "Monolayer"),
          ("bilayer",   "6x6x1", "Bilayer"),
          ("trilayer",  "4x4x1", "Trilayer")]

REF_BRANCH = "ai"        # DFT
CMP_BRANCH = "ml_ml"     # MLFF at its own equilibrium

THZ_TO_MEV = 4.13567
E_MIN, E_MAX = -5, 205   # 205 so the bulk LO at 197 meV is not clipped
Y_MAJOR = 50

STYLE = "publication"    # or "seaborn"
C_REF = "#1a1a1a"
C_CMP = "#b2182b"
LW_REF, LW_CMP = 1.35, 1.15
ALPHA_CMP = 0.9
DPI = 400
FIGSIZE = (9.6, 7.8)
DOS_FRAC = 0.30
OUT = "figure_dispersions"

# ----------------------------------------------------------------------
try:
    import seaborn as sns
    sns.set_context("paper", font_scale=1.15)
    if STYLE == "seaborn":
        sns.set_style("whitegrid")
    # take only the muted red; seaborn's grey is too weak for the reference
    # curve, which should stay near-black so the two are easy to separate
    C_CMP = sns.color_palette("deep")[3]
except ImportError:
    pass

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "font.size": 11,
    "axes.linewidth": 1.0, "axes.grid": False,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 5, "ytick.major.size": 5,
    "xtick.minor.size": 2.8, "ytick.minor.size": 2.8,
    "xtick.major.width": 1.0, "ytick.major.width": 1.0,
    "legend.frameon": False, "savefig.bbox": "tight",
})


def load_band(system, branch, dim):
    import yaml
    p = os.path.join(f"{system}_{branch}", dim, "band.yaml")
    if not os.path.exists(p):
        return None, None, None
    d = yaml.safe_load(open(p))
    f = np.array([[b["frequency"] for b in q["band"]]
                  for q in d["phonon"]]) * THZ_TO_MEV
    x = np.array([q["distance"] for q in d["phonon"]])
    seg = d.get("segment_nqpoint")
    return x, f, (list(np.cumsum(seg) - 1) if seg else None)


def load_dos(system, branch, dim):
    p = os.path.join(f"{system}_{branch}", dim, "total_dos.dat")
    if not os.path.exists(p):
        return None, None
    a = np.loadtxt(p, comments="#")
    return a[:, 0] * THZ_TO_MEV, a[:, 1]


fig = plt.figure(figsize=FIGSIZE)
# a narrow spacer column keeps the left DOS panel clear of the right
# dispersion; without it the two column groups run together
gs = fig.add_gridspec(2, 5,
                      width_ratios=[1, DOS_FRAC, 0.18, 1, DOS_FRAC],
                      hspace=0.30, wspace=0.10,
                      left=0.07, right=0.985, top=0.92, bottom=0.07)
COL = {0: 0, 1: 3}          # panel column -> gridspec column

for k, (system, dim, label) in enumerate(PANELS):
    r, c = divmod(k, 2)
    ax = fig.add_subplot(gs[r, COL[c]])
    axd = fig.add_subplot(gs[r, COL[c] + 1], sharey=ax)

    xr, fr, ends = load_band(system, REF_BRANCH, dim)
    xc, fc, _ = load_band(system, CMP_BRANCH, dim)
    if xr is None:
        ax.text(0.5, 0.5, f"{system}_{REF_BRANCH}/{dim}\nmissing",
                ha="center", va="center", transform=ax.transAxes)
        continue

    for b in range(fr.shape[1]):
        ax.plot(xr, fr[:, b], color=C_REF, lw=LW_REF, zorder=2,
                solid_capstyle="round")
    if xc is not None:
        for b in range(fc.shape[1]):
            ax.plot(xc, fc[:, b], color=C_CMP, lw=LW_CMP, alpha=ALPHA_CMP,
                    zorder=3, solid_capstyle="round")

    ax.set_ylim(E_MIN, E_MAX)
    ax.set_xlim(xr[0], xr[-1])
    ax.yaxis.set_major_locator(MultipleLocator(Y_MAJOR))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.axhline(0, color="0.65", ls=(0, (4, 4)), lw=0.8, zorder=1)
    if ends:
        for e in ends[:-1]:
            ax.axvline(xr[e], color="0.65", ls=(0, (4, 4)), lw=0.8, zorder=1)
        ax.set_xticks([xr[0]] + [xr[e] for e in ends])
        ax.set_xticklabels([r"$\Gamma$", "M", "K", r"$\Gamma$"])
    ax.tick_params(axis="x", which="minor", bottom=False, top=False)
    if c == 0:                       # only the left column carries the label
        ax.set_ylabel(r"$E$ (meV)")
    # panel label above the axes, so it never sits on a branch
    ax.set_title(f"({'abcd'[k]}) {label}", loc="left", fontsize=11.5, pad=6)

    er, dr = load_dos(system, REF_BRANCH, dim)
    ec, dc = load_dos(system, CMP_BRANCH, dim)
    if er is not None:
        axd.fill_betweenx(er, 0, dr, color=C_REF, alpha=0.13, lw=0)
        axd.plot(dr, er, color=C_REF, lw=LW_REF)
    if ec is not None:
        axd.plot(dc, ec, color=C_CMP, lw=LW_CMP, alpha=ALPHA_CMP)
    if er is None and ec is None:
        axd.text(0.5, 0.5, "no\nDOS", ha="center", va="center",
                 fontsize=8, color="0.5", transform=axd.transAxes)
    axd.set_xticks([])
    axd.tick_params(labelleft=False)
    axd.yaxis.set_minor_locator(AutoMinorLocator(2))
    axd.set_xlabel("DOS", fontsize=10, labelpad=3)
    axd.axhline(0, color="0.65", ls=(0, (4, 4)), lw=0.8)
    axd.set_xlim(left=0)
    axd.grid(False)

# one legend for the whole figure, clear of every panel
h = [plt.Line2D([], [], color=C_REF, lw=LW_REF, label="Ab initio"),
     plt.Line2D([], [], color=C_CMP, lw=LW_CMP, label="MLFF")]
fig.legend(handles=h, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.985),
           fontsize=11, handlelength=1.8, columnspacing=2.2)

fig.savefig(f"{OUT}.png", dpi=DPI)
fig.savefig(f"{OUT}.pdf")
print(f"wrote {OUT}.png and {OUT}.pdf")
