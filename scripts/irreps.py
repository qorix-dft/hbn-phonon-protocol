#!/usr/bin/env python3
"""
Gamma-point irreducible representations for Table 1.

    python irreps.py                    # every system found, _ai branch
    python irreps.py trilayer           # one system
    python irreps.py trilayer ml_ml     # a different branch

Builds force constants from FORCE_SETS via the phonopy API -- no FORCE_CONSTANTS
file needed, and no new VASP. Symmetrised to match FC_SYMMETRY in band.conf.

Expected label sets. Anything outside these is wrong:

  bulk       D_6h   A_1g A_2u B_1g B_2g E_1g E_1u E_2g E_2u
  monolayer  D_3h   A_1' A_2' E' A_1'' A_2'' E''
  bilayer    D_3d   A_1g A_2u E_g E_u
  trilayer   D_3h   A_1' A_2' E' A_1'' A_2'' E''

E''' does not exist in any point group -- do not carry it over from the draft.
If the monolayer returns D_6h labels, that calculation is not a monolayer.
"""
import sys, os, glob
import numpy as np

THZ_TO_CM = 33.35641

EXPECTED = {
    "bulk":      ("D_6h", "P6_3/mmc"),
    "monolayer": ("D_3h", "P-6m2"),
    "bilayer":   ("D_3d", "P-3m1"),
    "trilayer":  ("D_3h", "P-6m2"),
}


def dims_for(system, branch):
    base = f"{system}_{branch}"
    if not os.path.isdir(base):
        return []
    out = []
    for tag in os.listdir(base):
        try:
            d = tuple(int(x) for x in tag.split("x"))
        except ValueError:
            continue
        if len(d) == 3 and os.path.exists(
                os.path.join(base, tag, "FORCE_SETS")):
            out.append((tag, d))
    return sorted(out, key=lambda t: np.prod(t[1]))


def run(system, branch, tag):
    import phonopy
    wd = os.path.join(f"{system}_{branch}", tag)
    ph = phonopy.load(phonopy_yaml=os.path.join(wd, "phonopy_disp.yaml"),
                      force_sets_filename=os.path.join(wd, "FORCE_SETS"),
                      produce_fc=True, symmetrize_fc=True, log_level=0)
    ph.set_irreps([0, 0, 0], degeneracy_tolerance=1e-4)
    ir = ph.irreps

    pg = ir.get_point_group_symbol() if hasattr(ir, "get_point_group_symbol") \
        else getattr(ir, "_pointgroup_symbol", "?")
    exp_pg, exp_sg = EXPECTED.get(system, ("?", "?"))
    ok = exp_pg.replace("_", "").lower() in str(pg).replace("_", "").lower() \
        or str(pg) in exp_pg

    print(f"\n  {system}_{branch}/{tag}")
    print(f"    point group: {pg}   (expect {exp_pg} / {exp_sg})"
          f"   {'OK' if ok else '*** MISMATCH ***'}")

    freqs = ir.ddm.get_eigenvalues() if hasattr(ir, "ddm") else None
    labels = getattr(ir, "_ir_labels", None)
    band_idx = getattr(ir, "band_indices", None)

    if labels is None or band_idx is None:
        print("    (labels unavailable from this phonopy version -- run")
        print("     `phonopy --irreps=\"0 0 0\" --fc-symmetry` on the command")
        print("     line in that directory instead)")
        return

    ph.run_qpoints([[0, 0, 0]])
    w = ph.qpoints.frequencies[0] * THZ_TO_CM
    print(f"    {'branches':>12}{'freq (cm^-1)':>15}{'irrep':>10}{'deg':>5}")
    for bi, lab in zip(band_idx, labels):
        idx = list(bi)
        rng = f"{idx[0]+1}" if len(idx) == 1 else f"{idx[0]+1}-{idx[-1]+1}"
        print(f"    {rng:>12}{w[idx[0]]:>15.3f}{str(lab):>10}{len(idx):>5}")


if __name__ == "__main__":
    systems = ([sys.argv[1]] if len(sys.argv) > 1
               else sorted({p.rsplit("_", 1)[0].replace("_ml", "")
                            for p in glob.glob("*_ai") if os.path.isdir(p)}))
    branch = sys.argv[2] if len(sys.argv) > 2 else "ai"
    for s in systems:
        s = s.replace("_ai", "")
        dims = dims_for(s, branch)
        if not dims:
            print(f"\n  {s}_{branch}: nothing with FORCE_SETS")
            continue
        tag = dims[-1][0]        # largest supercell available
        try:
            run(s, branch, tag)
        except Exception as e:
            print(f"\n  {s}_{branch}/{tag}: {type(e).__name__}: {e}")
            print("    fall back to the command line:")
            print(f"      cd {s}_{branch}/{tag}")
            print(f"      phonopy -c POSCAR-unitcell --dim=\"{tag.replace('x',' ')}\" \\")
            print(f"              --irreps=\"0 0 0\" --tolerance=1e-4 --fc-symmetry")
