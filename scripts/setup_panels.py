#!/usr/bin/env python3
"""
Write band.conf and mesh.conf into every directory needed for panels (a)-(d),
then print the commands to run. No VASP -- this is all post-processing from
FORCE_SETS.

    python setup_panels.py            # write configs
    python setup_panels.py --check    # report what is present

Bulk gets NAC = .TRUE.; the slabs do not. In two dimensions the LO-TO
splitting vanishes at Gamma, so applying phonopy's 3D Ewald sum to a slab
would produce a splitting that is slab-slab interaction, not physics.
"""
import sys, os

PANELS = [("bulk",      "4x4x2", True),
          ("monolayer", "8x8x1", False),
          ("bilayer",   "6x6x1", False),
          ("trilayer",  "4x4x1", False)]
BRANCHES = ["ai", "ml_ml"]

BAND = """DIM = {a} {b} {c}
PRIMITIVE_AXES = 1 0 0  0 1 0  0 0 1
BAND = 0.0 0.0 0.0  0.5 0.0 0.0  0.3333 0.3333 0.0  0.0 0.0 0.0
BAND_LABELS = $\\Gamma$ M K $\\Gamma$
BAND_POINTS = 101
FC_SYMMETRY = .TRUE.
SYMMETRY = .TRUE.
SYMPREC = 1E-4
EIGENVECTORS = .TRUE.
{nac}"""

MESH = """DIM = {a} {b} {c}
PRIMITIVE_AXES = 1 0 0  0 1 0  0 0 1
MP = {mp}
FC_SYMMETRY = .TRUE.
SYMMETRY = .TRUE.
SYMPREC = 1E-4
DOS = .TRUE.
DOS_RANGE = 0 200 0.2
SIGMA = 0.5
{nac}"""


def main(check=False):
    todo = []
    for system, dim, nac in PANELS:
        d = [int(x) for x in dim.split("x")]
        mp = "30 30 10" if system == "bulk" else "40 40 1"
        nac_line = "NAC = .TRUE.\n" if nac else ""
        for br in BRANCHES:
            wd = os.path.join(f"{system}_{br}", dim)
            if not os.path.isdir(wd):
                print(f"  MISSING  {wd}")
                continue
            if not os.path.exists(os.path.join(wd, "FORCE_SETS")):
                print(f"  no FORCE_SETS in {wd}")
                continue
            if nac and not os.path.exists(os.path.join(wd, "BORN")):
                print(f"  *** {wd}: NAC requested but no BORN file. Copy it in,")
                print("      or the LO branch collapses onto TO in panel (a).")
            if check:
                have = [f for f in ("band.yaml", "total_dos.dat")
                        if os.path.exists(os.path.join(wd, f))]
                print(f"  {wd:<30}{', '.join(have) if have else 'nothing yet'}")
                continue
            with open(os.path.join(wd, "band.conf"), "w") as f:
                f.write(BAND.format(a=d[0], b=d[1], c=d[2], nac=nac_line))
            with open(os.path.join(wd, "mesh.conf"), "w") as f:
                f.write(MESH.format(a=d[0], b=d[1], c=d[2], mp=mp,
                                    nac=nac_line))
            todo.append(wd)

    if check or not todo:
        return
    print(f"\nwrote band.conf and mesh.conf into {len(todo)} directories\n")
    print("run:\n")
    for wd in todo:
        print(f"  (cd {wd} && phonopy -c POSCAR-unitcell band.conf"
              f" && phonopy -c POSCAR-unitcell mesh.conf)")
    print("\nthen:  python plot_dispersions.py")
    print("\nBAND_POINTS is 101 for a smooth figure -- higher than the 41 used")
    print("for the convergence checks, and free since it is only a Fourier")
    print("interpolation of force constants you already have.")


if __name__ == "__main__":
    main("--check" in sys.argv)
