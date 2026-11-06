#!/usr/bin/env python3
"""
Fit sigma_xx(a) from the three-point stress scan and solve for sigma_xx = 0.

Run from the directory containing monolayer_m1/, monolayer_00/, ... each with
a completed OUTCAR.

  python stressfit.py monolayer
  python stressfit.py            # all three
"""
import sys, os, re
import numpy as np

SYSTEMS = ["monolayer", "bilayer", "trilayer"]
TAGS    = ["m1", "00", "p1"]

# Expected slopes from a*k/(2V). The factor 2 is because scaling a strains the
# cell biaxially: eps_xx and eps_yy change together, so dE = 2*V*sigma_xx*deps.
# Verified against measured slopes to better than 0.2 %.
EXPECTED_SLOPE = {"monolayer": 935.5, "bilayer": 1543.0, "trilayer": 1961.5}

# ZA-driven criterion: |sigma_xx| must sit below this for the flexural branch
# to remain quadratic at the smallest sampled q
TOLERANCE = {"monolayer": 0.014, "bilayer": 0.024, "trilayer": 0.030}


def read_run(d):
    """Return (a, sigma_xx) from a finished VASP run directory."""
    oc = os.path.join(d, "OUTCAR")
    if not os.path.exists(oc):
        return None
    sig = None
    for line in open(oc, errors="ignore"):
        if "in kB" in line:
            parts = line.split()
            try:
                sig = float(parts[2])          # xx component
            except (IndexError, ValueError):
                pass
    poscar = os.path.join(d, "POSCAR")
    L = open(poscar).read().replace("\r", "").split("\n")
    scale = float(L[1])
    a = np.linalg.norm([float(x) for x in L[2].split()]) * scale
    return (a, sig) if sig is not None else None


def fit(system):
    pts = []
    for t in TAGS:
        r = read_run(f"{system}_{t}")
        if r:
            pts.append(r)
    if len(pts) < 2:
        print(f"{system}: only {len(pts)} finished runs -- need at least 2")
        return
    pts.sort()
    a = np.array([p[0] for p in pts])
    s = np.array([p[1] for p in pts])

    slope, intercept = np.polyfit(a, s, 1)
    a_zero = -intercept / slope
    resid = s - (slope * a + intercept)

    print(f"\n{system}")
    print(f"{'a (A)':>12}{'sigma_xx (kBar)':>18}{'resid':>10}")
    for ai, si, ri in zip(a, s, resid):
        print(f"{ai:>12.6f}{si:>18.5f}{ri:>10.5f}")

    exp = EXPECTED_SLOPE[system]
    print(f"\n  slope       = {slope:9.1f} kBar/A   (predicted {exp:.0f}, "
          f"{100*abs(abs(slope)-exp)/exp:.1f} % off)")
    print(f"  a(sigma=0)  = {a_zero:.6f} A")
    print(f"  shift from  = {(a_zero - a[len(a)//2])*1e3:+.3f} mA")
    print(f"  tolerance   = {TOLERANCE[system]:.3f} kBar")

    if abs(abs(slope) - exp) / exp > 0.20:
        print("  WARNING: slope differs from prediction by >20 % -- check the runs")
    if not (a.min() - 0.002 < a_zero < a.max() + 0.002):
        print("  WARNING: zero lies well outside the scan -- widen the bracket")


if __name__ == "__main__":
    targets = [sys.argv[1]] if len(sys.argv) > 1 else SYSTEMS
    for t in targets:
        fit(t)
