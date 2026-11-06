#!/usr/bin/env python3
"""
Assemble the Gamma-point table from the branch tree.

    python table1.py                       # plain text, all systems
    python table1.py --latex               # LaTeX rows for Table 1
    python table1.py --cmp ml_ai           # compare against the other branch
    python table1.py --dim 4x4x1           # pin the supercell (rarely useful)

Which MLFF branch to report
---------------------------
Default is `ml_ml`: the force field at its OWN relaxed geometry. That is what
it predicts when handed the material, so it is the number that belongs beside
an experimental value.

`ml_ai` puts the force field on the DFT geometry, isolating the force-constant
error alone. Useful as a diagnostic, not as a prediction -- and the two differ
by more than rounding. For bulk the curvature error is 13.0 cm^-1 while the
total is 9.1, because the geometric error partially CANCELS it. A force field
benchmarked at its own equilibrium can therefore look more accurate than its
interatomic forces are. Report `ml_ml` in the table and the decomposition in
the text.

Reads FORCE_SETS via the phonopy API, so no FORCE_CONSTANTS file and no rerun.
Prints DFT and MLFF side by side in cm^-1 at the format the table wants:
DFT [MLFF].

Irrep labels are entered by hand below -- phonopy reports them per run, and
they do not change between branches, so there is nothing to parse. Update
IRREPS if a supercell or system changes.

Acoustic modes come out as one 3-fold degenerate set at 0.000, which phonopy
cannot split ("Not found"). The assignment is fixed by symmetry:
  D_6h  A_2u (z) + E_1u (x,y)      D_3h  A_2'' + E'
  D_3d  A_2u     + E_u             (bulk needs --nac for the LO row)
"""
import sys, os, glob
import numpy as np

THZ_TO_CM = 33.35641
BRANCHES = ["ai", "ml_ai", "ml_ml"]

# branch index (0-based, inclusive start) -> (mode label, irrep)
IRREPS = {
    "monolayer": [(0, 2, "acoustic", "A_2'' + E'"),
                  (3, 3, "ZO",       "A_2''"),
                  (4, 5, "TO/LO",    "E'")],
    "bilayer":   [(0, 2, "acoustic",  "A_2u + E_u"),
                  (3, 4, "shear",     "E_g"),
                  (5, 5, "breathing", "A_1g"),
                  (6, 6, "ZO",        "A_2u"),
                  (7, 7, "ZO2",       "A_1g"),
                  (8, 9, "TO/LO",     "E_g"),
                  (10, 11, "TO/LO",   "E_u")],
    "trilayer":  [(0, 2, "acoustic",    "A_2'' + E'"),
                  (3, 4, "shear-1",     "E''"),
                  (5, 6, "shear-2",     "E'"),
                  (7, 7, "breathing-1", "A_1'"),
                  (8, 8, "breathing-2", "A_2''"),
                  (9, 9, "ZO",          "A_2''"),
                  (10, 10, "ZO",        "A_2''"),
                  (11, 11, "ZO",        "A_1'"),
                  (12, 13, "TO/LO",     "E'"),
                  (14, 15, "TO/LO",     "E''"),
                  (16, 17, "TO/LO",     "E'")],
    # with NAC along [1 0 0] the E_1u doublet splits into TO (perpendicular
    # to q) and LO (along q); E_2g carries no dipole and stays degenerate
    "bulk":      [(0, 2, "acoustic",  "A_2u + E_1u"),
                  (3, 4, "shear",     "E_2g"),
                  (5, 5, "ZO'",       "B_1g"),
                  (6, 6, "ZO",        "A_2u"),
                  (7, 7, "ZO2",       "B_1g"),
                  (8, 8, "TO",        "E_1u"),
                  (9, 10, "E_2g",     "E_2g"),
                  (11, 11, "LO",      "E_1u")],
}


def gamma(system, branch, dim=None):
    import phonopy
    base = f"{system}_{branch}"
    if not os.path.isdir(base):
        return None, None, False
    tags = []
    for t in os.listdir(base):
        if os.path.exists(os.path.join(base, t, "FORCE_SETS")):
            try:
                d = tuple(int(x) for x in t.split("x"))
            except ValueError:
                continue
            if len(d) == 3 and (dim is None or t == dim):
                tags.append((np.prod(d), t))
    if not tags:
        return None, None, False
    tag = sorted(tags)[-1][1]
    wd = os.path.join(base, tag)
    born = os.path.join(wd, "BORN")
    kw = dict(phonopy_yaml=os.path.join(wd, "phonopy_disp.yaml"),
              force_sets_filename=os.path.join(wd, "FORCE_SETS"),
              produce_fc=True, symmetrize_fc=True, log_level=0)
    nac = os.path.exists(born)
    if nac:
        kw.update(born_filename=born, is_nac=True)
    ph = phonopy.load(**kw)
    if nac:
        # LO-TO splitting is direction dependent at q = 0; take the in-plane
        # limit, which is the branch to compare against the ~1610 cm^-1 LO
        ph.run_qpoints([[0, 0, 0]], nac_q_direction=[1, 0, 0])
    else:
        ph.run_qpoints([[0, 0, 0]])
    return ph.qpoints.frequencies[0] * THZ_TO_CM, tag, nac


if __name__ == "__main__":
    latex = "--latex" in sys.argv
    dim = None
    if "--dim" in sys.argv:
        dim = sys.argv[sys.argv.index("--dim") + 1]
    cmp_branch = "ml_ml"
    if "--cmp" in sys.argv:
        cmp_branch = sys.argv[sys.argv.index("--cmp") + 1]
    if cmp_branch not in BRANCHES[1:]:
        print(f"--cmp must be one of {BRANCHES[1:]}")
        sys.exit(1)
    print(f"comparing _ai against _{cmp_branch}"
          + ("   (force field at its own equilibrium)" if cmp_branch == "ml_ml"
             else "   (force field on the DFT geometry -- diagnostic only)"))

    for system in ["bulk", "monolayer", "bilayer", "trilayer"]:
        ai, tag, nac = gamma(system, "ai", dim)
        if ai is None:
            print(f"\n{system}: no _ai data"
                  + (f" at DIM {dim}" if dim else ""))
            continue
        # Pin the comparison branch to the SAME supercell as _ai. Taking the
        # largest available in each branch independently silently pairs
        # different DIMs -- which is how _ai at 4x4x2 ended up compared
        # against 5x5x4 and produced a shear error 10 % too large.
        ml, mltag, mlnac = gamma(system, cmp_branch, tag)
        if ml is None:
            ml, mltag, mlnac = gamma(system, cmp_branch, None)
            if ml is not None:
                print(f"\n*** {system}: no _{cmp_branch} at {tag}; using "
                      f"{mltag} instead. Frequencies are converged to")
                print("    <0.05 cm^-1 across these sizes, but match them if"
                      " you can.")
        if ml is not None and nac != mlnac:
            print(f"*** {system}: NAC on for _ai={nac}, _{cmp_branch}={mlnac}."
                  " Copy BORN into both directories.")

        print(f"\n{'='*70}\n{system}   ({tag})"
              + ("   NAC on" if nac else "   no NAC")
              + ("" if ml is not None else f"   [no _{cmp_branch} -- DFT only]"))
        if not latex:
            print(f"{'branch':>9}{'mode':>13}{'DFT':>11}{'MLFF':>11}"
                  f"{'err %':>9}{'irrep':>13}")
        for lo, hi, mode, irrep in IRREPS.get(system, []):
            if hi >= ai.size:
                continue
            a = ai[lo]
            m = ml[lo] if ml is not None and lo < ml.size else None
            rng = f"{lo+1}" if lo == hi else f"{lo+1}-{hi+1}"
            if latex:
                val = (f"${a:.2f}\\,[{m:.2f}]$" if m is not None
                       else f"${a:.2f}$")
                print(f"  {mode:<12} & {val} & (${irrep}$) \\\\")
            else:
                err = f"{100*(m/a-1):+8.1f}%" if (m is not None and a > 1) else f"{'--':>9}"
                mm = f"{m:>11.2f}" if m is not None else f"{'--':>11}"
                print(f"{rng:>9}{mode:>13}{a:>11.2f}{mm}{err}{irrep:>13}")

    print("\n  Reported MLFF values are at the branch selected with --cmp")
    print(f"  (default ml_ml). Copy BORN into <system>_<branch>/<DIM>/ for")
    print("  every branch you compare, or the bulk LO row will not split.")
    print("\n  NAC is applied automatically wherever a BORN file sits in the")
    print("  branch directory. Without it the bulk E_1u and E_2g manifolds")
    print("  stay degenerate at Gamma and the LO row is wrong.")
    print("  Omit --dim to take the largest supercell available per system;")
    print("  a single --dim value cannot match slabs (N N 1) and bulk (N N M).")
