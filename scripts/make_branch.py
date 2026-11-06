#!/usr/bin/env python3
"""
Build one <system>_<branch>/<DIM>/ directory: supercell, displacements, and a
disp-NNN/ per displacement with all VASP inputs.

    python make_branch.py <system> <branch> <DIM> <unitcell POSCAR>

  <branch> in {ai, ml_ai, ml_ml}
  <DIM>    like 4x4x1 or 4x4x2

Examples
--------
    python make_branch.py bulk      ai    4x4x2 cells_zerostress/POSCAR_bulk.vasp
    python make_branch.py monolayer ml_ml 4x4x1 CONTCAR_monolayer
    python make_branch.py bilayer   ml_ml 4x4x1 CONTCAR_bilayer
    python make_branch.py trilayer  ml_ml 4x4x1 CONTCAR_trilayer
    python make_branch.py bulk      ml_ml 4x4x2 CONTCAR_bulk_ml

Needs in the working directory: INCAR_ai and INCAR_ml (the Stage 3 blocks),
POTCAR, kpoints/, and ML_FF for the MLFF branches.

Existing disp-NNN/ with a finished vasprun.xml are left untouched.
"""
import sys, os, glob, shutil, subprocess

KPT_DIR = "kpoints"
BRANCHES = {"ai": "INCAR_ai", "ml_ai": "INCAR_ml", "ml_ml": "INCAR_ml"}


def run_done(d):
    v = os.path.join(d, "vasprun.xml")
    if not os.path.exists(v) or os.path.getsize(v) < 1000:
        return False
    with open(v, "rb") as f:
        f.seek(max(0, os.path.getsize(v) - 400))
        return b"</modeling>" in f.read()


def main(system, branch, dim, cellfile):
    if branch not in BRANCHES:
        print(f"branch must be one of {list(BRANCHES)}")
        return
    incar = BRANCHES[branch]
    d = [int(x) for x in dim.split("x")]
    wd = os.path.join(f"{system}_{branch}", dim)
    os.makedirs(wd, exist_ok=True)

    # 1. unit cell into place
    shutil.copy(cellfile, os.path.join(wd, "POSCAR-unitcell"))

    # 2. supercell + symmetry-reduced displacements
    #    writes SPOSCAR, phonopy_disp.yaml, POSCAR-001 ... POSCAR-NNN
    for exe in ("phonopy", "phonopy-init"):      # v4 moved -d to phonopy-init
        r = subprocess.run([exe, "-d", f"--dim={d[0]} {d[1]} {d[2]}",
                            "-c", "POSCAR-unitcell"],
                           cwd=wd, capture_output=True, text=True)
        if r.returncode == 0:
            break
    else:
        print(r.stderr.strip()[:300])
        return

    poscars = sorted(glob.glob(os.path.join(wd, "POSCAR-[0-9]*")))
    print(f"{system}_{branch}/{dim}: {len(poscars)} displacements")

    # 3. band.conf for the post-processing step
    with open(os.path.join(wd, "band.conf"), "w") as f:
        f.write(f"""DIM = {d[0]} {d[1]} {d[2]}
PRIMITIVE_AXES = 1 0 0  0 1 0  0 0 1
BAND = 0.0 0.0 0.0  0.5 0.0 0.0  0.3333 0.3333 0.0  0.0 0.0 0.0
BAND_LABELS = $\\Gamma$ M K $\\Gamma$
BAND_POINTS = 41
FC_SYMMETRY = .TRUE.
SYMMETRY = .TRUE.
SYMPREC = 1E-4
EIGENVECTORS = .TRUE.
""")

    # 4. one directory per displacement, with inputs
    kpt = os.path.join(KPT_DIR, f"KPOINTS_{system}_{dim}")
    need = [incar, "POTCAR", kpt]
    missing = [f for f in need if not os.path.exists(f)]
    if missing:
        print(f"  missing: {', '.join(missing)} -- POSCARs written, inputs not")
        return
    ml = "ML_FF" if branch.startswith("ml") else None
    if ml and not os.path.exists(ml):
        print("  warning: no ML_FF in the working directory")
        ml = None

    made = kept = 0
    for p in poscars:
        num = os.path.basename(p).split("-")[1]
        dd = os.path.join(wd, f"disp-{num}")
        if run_done(dd):
            kept += 1
            continue
        os.makedirs(dd, exist_ok=True)
        shutil.copy(p, os.path.join(dd, "POSCAR"))
        shutil.copy(incar, os.path.join(dd, "INCAR"))
        shutil.copy("POTCAR", dd)
        shutil.copy(kpt, os.path.join(dd, "KPOINTS"))
        if ml:
            shutil.copy(ml, dd)
        made += 1
    msg = f"  {made} disp-NNN/ prepared"
    if kept:
        msg += f", {kept} already finished and left untouched"
    print(msg + f"  (INCAR {incar}, KPOINTS {os.path.basename(kpt)})")
    print(f"\n  run VASP in {wd}/disp-*/ then:")
    print(f"    cd {wd}")
    print(f"    phonopy -f disp-*/vasprun.xml")
    print(f"    phonopy -c POSCAR-unitcell band.conf")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)
    main(*sys.argv[1:])
