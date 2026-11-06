#!/usr/bin/env python3
"""
Panels (e)-(h): phonon density of states of the defective trilayer supercells,
DFT against the MLFF.

    python plot_defect_dos.py

Input files are two columns, frequency in THz and DOS, as written by phonopy.
Expected in the working directory:

    phdos_cbcn_ai.dat   phdos_cbcn_ml.dat      -> panel (e)
    phdos_cn_ai.dat     phdos_cn_ml.dat        -> panel (f)
    phdos_cb_ai.dat     phdos_cb_ml.dat        -> panel (g)
    phdos_tri_ai.dat    phdos_tri_ml.dat       -> panel (h), pristine reference

A missing pair is drawn as an empty panel with a note, so the figure can be
built before every file exists.

IMPORTANT -- what the shading means
-----------------------------------
The coloured bands are FIXED ENERGY WINDOWS, not an eigenvector projection.
They mark the regions where in-plane and out-of-plane character dominate; they
do not measure it. Describe them that way in the caption, or compute the real
projection by weighting each mode by |e_z|^2 against |e_x|^2 + |e_y|^2 from a
mesh run with EIGENVECTORS = .TRUE.
"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

# ----------------------------------------------------------------------
# SETTINGS
# ----------------------------------------------------------------------
PANELS = [("cbcn", r"$\mathrm{C_B C_N}$"),
          ("cn",   r"$\mathrm{C_N}$"),
          ("cb",   r"$\mathrm{C_B}$"),
          ("tri",  r"tri-$h$BN")]

THZ_TO_MEV = 4.13567
E_MIN, E_MAX = 0, 192
X_MAJOR = 25

# window edges in THz, as in the original figure; character alternates
# out-of-plane / in-plane / out-of-plane / in-plane
EDGES_THZ = [0, 10, 16, 24, 48]
CHARACTER = ["out", "in", "out", "in"]

C_REF = "#1a1a1a"
C_CMP = "#c44e52"
C_IN, C_OUT = "#4c72b0", "#dbaa00"
FILL_ALPHA = 0.20
LW_REF, LW_CMP = 1.5, 1.3

NORMALISE = True         # scale each panel to its own DFT maximum
# Ignore everything below this when setting the scale. A Gamma-sampled
# supercell puts its three acoustic modes at exactly zero, which the smearing
# turns into a spike far taller than any real feature; normalising to it
# flattens the rest of the panel. The pristine cell shows this most strongly
# because its acoustic modes are not shifted by a defect.
NORM_EMIN = 5.0          # meV
Y_HEADROOM = 1.12        # top of the y-axis, in units of the scale
DPI = 400
FIGSIZE = (7.6, 8.4)
OUT = "figure_defect_dos"

# ----------------------------------------------------------------------
try:
    import seaborn as sns
    sns.set_context("paper", font_scale=1.15)
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
    "legend.frameon": False, "savefig.bbox": "tight",
})


def load(tag, branch):
    p = f"phdos_{tag}_{branch}.dat"
    if not os.path.exists(p):
        return None, None
    a = np.loadtxt(p, comments="#")
    e = a[:, 0] * THZ_TO_MEV
    m = e >= -1e-9                    # phonopy writes a small negative tail
    return e[m], a[m, 1]


fig, axes = plt.subplots(len(PANELS), 1, figsize=FIGSIZE, sharex=True)
fig.subplots_adjust(hspace=0.12, left=0.11, right=0.98, top=0.965, bottom=0.075)

for k, (tag, label) in enumerate(PANELS):
    ax = axes[k]
    er, dr = load(tag, "ai")
    ec, dc = load(tag, "ml")

    if er is None and ec is None:
        ax.text(0.5, 0.5, f"phdos_{tag}_*.dat not found",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=9, color="0.5")
    else:
        if NORMALISE and er is not None:
            m = er > NORM_EMIN
            scale = dr[m].max() if m.any() else dr.max()
        else:
            scale = 1.0
        if er is not None:
            # shade the energy windows under the DFT curve
            for j in range(len(EDGES_THZ) - 1):
                lo = EDGES_THZ[j] * THZ_TO_MEV
                hi = EDGES_THZ[j + 1] * THZ_TO_MEV
                m = (er >= lo) & (er <= hi)
                if m.any():
                    ax.fill_between(er[m], 0, dr[m] / scale, lw=0,
                                    color=C_IN if CHARACTER[j] == "in" else C_OUT,
                                    alpha=FILL_ALPHA, zorder=1)
            ax.plot(er, dr / scale, color=C_REF, lw=LW_REF, zorder=3)
        if ec is not None:
            ax.plot(ec, dc / scale, color=C_CMP, lw=LW_CMP, zorder=4)

    ax.set_xlim(E_MIN, E_MAX)
    if NORMALISE:
        ax.set_ylim(0, Y_HEADROOM)      # clips the acoustic spike at E ~ 0
    else:
        ax.set_ylim(bottom=0)
    ax.xaxis.set_major_locator(MultipleLocator(X_MAJOR))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.set_yticklabels([])
    ax.set_ylabel("DOS (a.u.)")
    ax.text(0.015, 0.90, f"({'efgh'[k]}) {label}", transform=ax.transAxes,
            va="top", ha="left", fontsize=12)

    if k == 0:
        h = [plt.Line2D([], [], color=C_REF, lw=LW_REF, label="Ab initio"),
             plt.Line2D([], [], color=C_CMP, lw=LW_CMP, label="MLFF"),
             plt.Rectangle((0, 0), 1, 1, fc=C_IN, alpha=FILL_ALPHA,
                           label="In-plane"),
             plt.Rectangle((0, 0), 1, 1, fc=C_OUT, alpha=FILL_ALPHA,
                           label="Out-of-plane")]
        ax.legend(handles=h, loc="upper right", ncol=2, fontsize=9.5,
                  handlelength=1.5, columnspacing=1.4, labelspacing=0.35)

axes[-1].set_xlabel(r"$E$ (meV)")

fig.savefig(f"{OUT}.png", dpi=DPI)
fig.savefig(f"{OUT}.pdf")
print(f"wrote {OUT}.png and {OUT}.pdf")
