#!/usr/bin/env python3
"""
MLFF vs ab initio phonon dispersion benchmark for trilayer h-BN.

Reads phonopy band-style .dat files (col 0 = q-path coordinate,
cols 1..N = branch frequencies in THz) and produces:
  - error statistics resolved by branch, manifold and energy window
  - a two-panel comparison figure (dispersion overlay + signed residual)
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# ----------------------------------------------------------------------
# SETTINGS
# ----------------------------------------------------------------------
F_AI        = "/mnt/user-data/uploads/phdispai.dat"
F_ML        = "/mnt/user-data/uploads/phdispmlff.dat"

THZ_TO_MEV  = 4.13567          # 1 THz -> meV
MEV_TO_CM   = 8.06554          # 1 meV -> cm^-1

# High-symmetry point indices along the stored q-path
HSP         = {"$\\Gamma$": 0, "M": 39, "K": 79, "$\\Gamma$ ": 119}

# Mode-manifold grouping (0-based branch indices) for trilayer h-BN,
# 6 atoms/cell -> 18 branches
MANIFOLDS = [
    ("Acoustic",              list(range(0, 3))),
    ("Interlayer shear",      list(range(3, 7))),
    ("Interlayer breathing",  list(range(7, 9))),
    ("ZO",                    list(range(9, 12))),
    ("LO/TO",                 list(range(12, 18))),
]

OUT_PNG     = "/mnt/user-data/outputs/trilayer_phonon_benchmark.png"
OUT_PDF     = "/mnt/user-data/outputs/trilayer_phonon_benchmark.pdf"

C_AI        = "#1a1a1a"        # ab initio
C_ML        = "#99180f"        # MLFF
LW          = 1.3
DPI         = 300
FIGSIZE     = (7.0, 6.4)
HEIGHT_RATIO = (2.6, 1.0)

# ----------------------------------------------------------------------
# STYLE
# ----------------------------------------------------------------------
mpl.rcParams.update({
    "font.family":       "serif",
    "font.serif":        ["DejaVu Serif"],
    "mathtext.fontset":  "dejavuserif",
    "font.size":         10,
    "axes.linewidth":    0.9,
    "xtick.direction":   "in",
    "ytick.direction":   "in",
    "xtick.top":         True,
    "ytick.right":       True,
    "xtick.major.size":  5,
    "ytick.major.size":  5,
    "xtick.minor.size":  2.8,
    "ytick.minor.size":  2.8,
    "xtick.major.width": 0.9,
    "ytick.major.width": 0.9,
    "legend.frameon":    False,
})

# ----------------------------------------------------------------------
# LOAD
# ----------------------------------------------------------------------
ai = np.loadtxt(F_AI)
ml = np.loadtxt(F_ML)
assert ai.shape == ml.shape, "grids differ"

q = ai[:, 0]
A = ai[:, 1:] * THZ_TO_MEV
M = ml[:, 1:] * THZ_TO_MEV
D = M - A
nq, nb = A.shape

x = np.arange(nq)


def stats(d):
    return (np.sqrt((d ** 2).mean()), np.abs(d).mean(),
            np.abs(d).max(), d.mean())


# ----------------------------------------------------------------------
# REPORT
# ----------------------------------------------------------------------
print(f"{nb} branches, {nq} q-points\n")

print("=== Gamma-point frequencies ===")
print(f"{'branch':>7} {'manifold':>21} {'AI (meV)':>10} {'MLFF':>9} "
      f"{'d':>8} {'rel %':>8} {'AI cm-1':>10} {'ML cm-1':>9}")
name_of = {}
for nm, idx in MANIFOLDS:
    for b in idx:
        name_of[b] = nm
for b in range(nb):
    rel = 100 * D[0, b] / A[0, b] if abs(A[0, b]) > 0.5 else np.nan
    print(f"{b+1:>7} {name_of[b]:>21} {A[0,b]:>10.3f} {M[0,b]:>9.3f} "
          f"{D[0,b]:>+8.3f} {rel:>8.2f} "
          f"{A[0,b]*MEV_TO_CM:>10.2f} {M[0,b]*MEV_TO_CM:>9.2f}")

print("\n=== Error by manifold (full path, meV) ===")
print(f"{'manifold':>21} {'RMSE':>7} {'MAE':>7} {'max|d|':>8} {'mean':>8}")
r, m, x_, s = stats(D)
print(f"{'ALL':>21} {r:>7.3f} {m:>7.3f} {x_:>8.3f} {s:>+8.3f}")
for nm, idx in MANIFOLDS:
    r, m, x_, s = stats(D[:, idx])
    print(f"{nm:>21} {r:>7.3f} {m:>7.3f} {x_:>8.3f} {s:>+8.3f}")

print("\n=== Error by energy window (binned on ab initio energy) ===")
for lo, hi in [(0, 20), (20, 50), (50, 110), (110, 150), (150, 200)]:
    msk = (A >= lo) & (A < hi)
    if msk.sum():
        print(f"{lo:>4}-{hi:<4} meV  n={msk.sum():>5}  "
              f"RMSE={np.sqrt((D[msk]**2).mean()):.3f}  "
              f"MAE={np.abs(D[msk]).mean():.3f}  "
              f"rel.MAE={100*np.abs(D[msk]).mean()/A[msk].mean():.2f}%")

print("\n=== Degeneracy of the interlayer shear doublets at Gamma ===")
for nm, pair in [("shear-1", (3, 4)), ("shear-2", (5, 6))]:
    print(f" {nm}: AI split = {A[0,pair[1]]-A[0,pair[0]]:.4f} meV, "
          f"MLFF split = {M[0,pair[1]]-M[0,pair[0]]:.4f} meV")

print("\n=== LO/TO spread at Gamma (tests for non-analytic correction) ===")
print(f" AI   spread = {A[0,12:].max()-A[0,12:].min():.4f} meV")
print(f" MLFF spread = {M[0,12:].max()-M[0,12:].min():.4f} meV")

print("\n=== Flexural (ZA) branch near Gamma ===")
print(f"{'q':>9} {'AI b1':>9} {'MLFF b1':>9}")
for i in range(6):
    print(f"{q[i]:>9.5f} {A[i,0]:>9.3f} {M[i,0]:>9.3f}")
for nm, X in [("AI", A), ("MLFF", M)]:
    sl = slice(3, 15)
    w = X[sl, 0]
    if (w > 0).all():
        n = np.polyfit(np.log(q[sl]), np.log(w), 1)[0]
        print(f" {nm}: fitted exponent n = {n:.2f}  (harmonic 2D theory: 2)")

# ----------------------------------------------------------------------
# FIGURE
# ----------------------------------------------------------------------
fig, (ax, axr) = plt.subplots(
    2, 1, figsize=FIGSIZE, sharex=True,
    gridspec_kw={"height_ratios": HEIGHT_RATIO, "hspace": 0.06})

for b in range(nb):
    ax.plot(x, A[:, b], color=C_AI, lw=LW, zorder=2)
    ax.plot(x, M[:, b], color=C_ML, lw=LW, zorder=3)
ax.plot([], [], color=C_AI, lw=LW, label="Ab initio")
ax.plot([], [], color=C_ML, lw=LW, label="MLFF")

ax.set_ylabel(r"$E$ (meV)")
ax.set_ylim(-5, 180)
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.axhline(0, color="0.55", ls=(0, (5, 4)), lw=0.8, zorder=1)
ax.legend(loc="upper left", fontsize=9, bbox_to_anchor=(0.015, 0.99),
          handlelength=1.6, labelspacing=0.35)

for b in range(nb):
    axr.plot(x, D[:, b], color=C_ML, lw=0.8, alpha=0.75)
axr.axhline(0, color="0.55", ls=(0, (5, 4)), lw=0.8)
axr.set_ylabel(r"$\Delta E$ (meV)")
axr.set_ylim(-3.2, 2.6)
axr.yaxis.set_minor_locator(AutoMinorLocator(2))

for a in (ax, axr):
    for lab, i in HSP.items():
        if 0 < i < nq - 1:
            a.axvline(i, color="0.55", ls=(0, (5, 4)), lw=0.8, zorder=1)
    a.set_xlim(0, nq - 1)

axr.set_xticks(list(HSP.values()))
axr.set_xticklabels([k.strip() for k in HSP.keys()])
axr.tick_params(axis="x", which="minor", top=False, bottom=False)

rmse = np.sqrt((D ** 2).mean())
axr.text(0.30, 0.06, f"RMSE = {rmse:.2f} meV",
         transform=axr.transAxes, fontsize=9, va="bottom")

fig.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight")
fig.savefig(OUT_PDF, bbox_inches="tight")
print(f"\nwrote {OUT_PNG}\nwrote {OUT_PDF}")
