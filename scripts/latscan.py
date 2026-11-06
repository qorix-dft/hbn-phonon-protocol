#!/usr/bin/env python3
"""
In-plane lattice constant scan for h-BN slabs.

Replaces ISIF=3 relaxation, which fails on cells containing vacuum because
the sigma_zz stress component is numerical noise referenced to a volume
dominated by empty space (ZBRENT bracketing failure).

Usage
-----
  generate:  python latscan.py gen   <system>
  inspect:   python latscan.py check <system>    # what was found, and where
  fit:       python latscan.py fit   <system>

<system> in {monolayer, bilayer, trilayer}
"""
import sys, os, re
import numpy as np

# ----------------------------------------------------------------------
# SETTINGS
# ----------------------------------------------------------------------
ENCUT_NOTE = 650            # run the scan at the PRODUCTION cutoff, not higher
A_SPAN   = 0.06             # +/- range in Angstrom
N_POINTS = 7
VACUUM   = 15.157           # vacuum gap, Angstrom
ROOT     = "latscan"

# Centre each scan on that system's own AI result, and start from the
# AI-relaxed interlayer spacing. The ionic relaxation will move the spacing
# anyway, but starting close saves steps and keeps it from wandering.
A_CENTRE = {"monolayer": 2.508960, "bilayer": 2.507926, "trilayer": 2.507510}
D_INTER  = {"monolayer": 3.28229155, "bilayer": 3.3597, "trilayer": 3.3384}

B_SITE = np.array([1/3, 2/3])
N_SITE = np.array([2/3, 1/3])

NLAYERS = {"monolayer": 1, "bilayer": 2, "trilayer": 3}


def write_poscar(path, a, nlayers, title, d_inter):
    a1 = np.array([a, 0.0, 0.0])
    a2 = np.array([-a/2, a*np.sqrt(3)/2, 0.0])
    c  = (nlayers - 1) * d_inter + VACUUM
    z0 = (c - (nlayers - 1) * d_inter) / 2
    rows = []
    for L in range(nlayers):
        zc = z0 + L * d_inter
        bxy, nxy = (B_SITE, N_SITE) if L % 2 == 0 else (N_SITE, B_SITE)
        rows.append(('B', bxy[0], bxy[1], zc/c))
        rows.append(('N', nxy[0], nxy[1], zc/c))
    rows.sort(key=lambda r: r[0])
    nb = sum(1 for r in rows if r[0] == 'B')
    with open(path, "w") as f:
        f.write(f"{title}\n1.0\n")
        for v in (a1, a2, np.array([0, 0, c])):
            f.write(f"  {v[0]:22.16f}{v[1]:22.16f}{v[2]:22.16f}\n")
        f.write(f"  B  N\n  {nb:3d} {len(rows)-nb:3d}\nDirect\n")
        for r in rows:
            f.write(f"  {r[1]:22.16f}{r[2]:22.16f}{r[3]:22.16f}\n")


def gen(system):
    nl = NLAYERS[system]
    centre = A_CENTRE[system]
    avals = np.linspace(centre - A_SPAN, centre + A_SPAN, N_POINTS)
    base = os.path.join(ROOT, system)
    os.makedirs(base, exist_ok=True)
    for a in avals:
        d = os.path.join(base, f"a_{a:.5f}")
        os.makedirs(d, exist_ok=True)
        write_poscar(os.path.join(d, "POSCAR"), a, nl,
                     f"hBN {system} | a={a:.5f} | vacuum {VACUUM} A",
                     D_INTER[system])
    print(f"{system}: {N_POINTS} cells in {base}/")
    print(f"  a from {avals[0]:.5f} to {avals[-1]:.5f} A")
    if nl == 1:
        print("  monolayer has 0 free ionic DOF -> run static (IBRION=-1, NSW=0)")
    else:
        print("  run ionic relaxation only (IBRION=1, ISIF=2) -- NOT ISIF=3")


def read_energy(d):
    """Total energy (eV) from a finished run, whichever file has it.

    An MLFF run (ML_MODE = run) does no SCF, so the OUTCAR does not carry the
    'free  energy   TOTEN' line an ab initio run writes. vasprun.xml records
    e_fr_energy the same way in both cases, so read that first.
    """
    # 1. vasprun.xml -- format independent, works for AI and MLFF alike
    v = os.path.join(d, "vasprun.xml")
    if os.path.exists(v):
        pat = re.compile(r'e_fr_energy"?\s*>\s*([-+0-9.eEdD]+)')
        e = None
        for line in open(v, errors="ignore"):
            m = pat.search(line)
            if m:
                try:
                    e = float(m.group(1).replace("D", "E"))
                except ValueError:
                    pass
        if e is not None:
            return e, "vasprun.xml"

    # 2. OUTCAR -- several spellings, including the MLFF variants
    oc = os.path.join(d, "OUTCAR")
    if os.path.exists(oc):
        pats = [("free  energy   TOTEN", -2),
                ("ML energy  without entropy", -1),
                ("ML free energy   TOTEN", -2),
                ("free energy    TOTEN", -2)]
        for pat, idx in pats:
            e = None
            for line in open(oc, errors="ignore"):
                if pat in line:
                    try:
                        e = float(line.split()[idx])
                    except (IndexError, ValueError):
                        pass
            if e is not None:
                return e, f"OUTCAR ({pat.strip()})"

    # 3. OSZICAR -- last resort
    oz = os.path.join(d, "OSZICAR")
    if os.path.exists(oz):
        e = None
        for line in open(oz, errors="ignore"):
            if "E0=" in line:
                try:
                    e = float(line.split("E0=")[1].split()[0])
                except (IndexError, ValueError):
                    pass
        if e is not None:
            return e, "OSZICAR (E0)"

    return None, None


def check(system):
    """Report what is present in each scan directory and where the energy
    was found. Run this when `fit` says it has no finished runs."""
    base = os.path.join(ROOT, system)
    if not os.path.isdir(base):
        print(f"no directory {base}")
        return
    print(f"\n{system}")
    print(f"{'directory':>14}{'vasprun':>9}{'OUTCAR':>8}{'OSZICAR':>9}"
          f"{'energy (eV)':>15}   source")
    for d in sorted(os.listdir(base)):
        wd = os.path.join(base, d)
        if not os.path.isdir(wd):
            continue
        e, src = read_energy(wd)
        print(f"{d:>14}"
              f"{'yes' if os.path.exists(os.path.join(wd,'vasprun.xml')) else '-':>9}"
              f"{'yes' if os.path.exists(os.path.join(wd,'OUTCAR')) else '-':>8}"
              f"{'yes' if os.path.exists(os.path.join(wd,'OSZICAR')) else '-':>9}"
              + (f"{e:>15.6f}   {src}" if e is not None
                 else f"{'--':>15}   not found"))


def fit(system, deg=4):
    """Fit E(a). Default quartic: E(a) has significant cubic anharmonicity over
    a +/-0.04 A window, and a parabola biases the minimum by ~0.6 mA."""
    base = os.path.join(ROOT, system)
    a_list, e_list = [], []
    for d in sorted(os.listdir(base)):
        wd = os.path.join(base, d)
        if not os.path.isdir(wd):
            continue
        e, _ = read_energy(wd)
        if e is not None:
            a_list.append(float(d.split("_")[1]))
            e_list.append(e)
    if len(a_list) < deg + 2:
        print(f"only {len(a_list)} finished runs -- need at least {deg+2}")
        print(f"run:  python {os.path.basename(__file__)} check {system}")
        return
    a = np.array(a_list); E = np.array(e_list)

    # centre the abscissa: a Vandermonde in a ~ 2.5 over a 0.08 A window is
    # badly conditioned at degree 4 and makes np.polyfit's covariance estimate
    # go non-positive-semidefinite.
    a0 = a.mean()
    x = a - a0

    def minimum(p):
        r = [xx.real for xx in np.roots(np.polyder(p))
             if abs(xx.imag) < 1e-9 and x.min() < xx.real < x.max()]
        return min(r, key=lambda xx: np.polyval(p, xx)) + a0 if r else None

    p2, pn = np.polyfit(x, E, 2), np.polyfit(x, E, deg)
    a2, an = minimum(p2), minimum(pn)
    resid = E - np.polyval(pn, x)

    print(f"\n{system}: fitted from {len(a)} points, degree {deg}")
    print(f"{'a (A)':>10}{'E (eV)':>16}{'resid (meV)':>14}")
    for ai, Ei, ri in zip(a, E, resid):
        print(f"{ai:>10.5f}{Ei:>16.6f}{ri*1e3:>14.4f}")

    print(f"\n  a_opt (deg {deg}) = {an:.5f} A        <-- use this")
    print(f"  a_opt (quadratic) = {a2:.5f} A   (biased by {(a2-an)*1e3:+.2f} mA)")
    print(f"  curvature at minimum = {np.polyval(np.polyder(pn, 2), an - a0):.4f} eV/A^2")
    print(f"  max |residual| = {np.abs(resid).max()*1e3:.4f} meV")

    if an is None or not (a.min() < an < a.max()):
        print("  WARNING: minimum outside the scan range -- widen A_SPAN")
    if np.abs(resid).max() * 1e3 > 0.05:
        print("  WARNING: residual > 0.05 meV even at this order -- check convergence")

    # bootstrap uncertainty
    try:
        p_, cov = np.polyfit(x, E, deg, cov=True)
        samp = []
        for _ in range(2000):
            ps = np.random.multivariate_normal(p_, cov)
            m = minimum(ps)
            if m is not None:
                samp.append(m)
        if samp:
            print(f"  bootstrap uncertainty = +/- {np.std(samp)*1e3:.3f} mA")
    except Exception:
        pass


if __name__ == "__main__":
    cmds = {"gen": gen, "fit": fit, "check": check}
    if len(sys.argv) < 3 or sys.argv[1] not in cmds or sys.argv[2] not in NLAYERS:
        print(__doc__); sys.exit(1)
    cmds[sys.argv[1]](sys.argv[2])
