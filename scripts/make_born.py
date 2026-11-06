#!/usr/bin/env python3
"""
Build the phonopy BORN file from a VASP LCALCEPS / LEPSILON run, and check it.

    python make_born.py                 # reads ./vasprun.xml, writes ./BORN
    python make_born.py <dir>           # from another directory
    python make_born.py <dir> --check   # report only, write nothing

Equivalent one-liner if you trust the run:

    phonopy-vasp-born > BORN

This script exists because that one-liner writes whatever it finds without
telling you whether it is sensible. The checks below matter:

  charge neutrality   sum of Z* over all atoms must vanish. VASP does not
                      enforce it, and a finite residual propagates straight
                      into the LO-TO splitting. phonopy can impose it, but you
                      want to know how large it was.
  dielectric tensor   for bulk h-BN, eps_inf is ~4.8-5.1 in plane and
                      ~2.6-2.9 along c. Note this is the ELECTRONIC tensor;
                      the static constant is larger (~3.4 out of plane) because
                      it carries the ionic contribution. A near-isotropic
                      tensor means the run did not resolve the anisotropy --
                      h-BN is strongly hyperbolic.
  Born charges        Z*(B) ~ +2.7 in plane, Z*(N) ~ -2.7. Out-of-plane
                      components are much smaller, ~ +/-0.8.

BORN file format
  line 1        conversion factor (14.400 for the phonopy default units)
  line 2        dielectric tensor, 9 values: xx xy xz yx yy yz zx zy zz
  line 3..      Z* for each symmetry-INEQUIVALENT atom, 9 values each,
                in the order phonopy's primitive cell lists them
"""
import sys, os, re
import numpy as np

FACTOR = 14.400


def parse_vasprun(path):
    """Return (epsilon 3x3, born N x 3 x 3) from vasprun.xml."""
    txt = open(path, errors="ignore").read()

    def grab_matrix(name):
        m = re.search(rf'<varray name="{name}"[^>]*>(.*?)</varray>',
                      txt, re.S)
        if not m:
            return None
        rows = re.findall(r'<v>([^<]*)</v>', m.group(1))
        return np.array([[float(x) for x in r.split()] for r in rows])

    eps = grab_matrix("epsilon")
    if eps is None:
        eps = grab_matrix("epsilon_scf")
    if eps is None:
        # LCALCEPS writes it under a different tag in some versions
        eps = grab_matrix("epsilon_ion")

    m = re.search(r'<array name="born_charges".*?</array>', txt, re.S)
    born = None
    if m:
        sets = re.findall(r'<set[^>]*>(.*?)</set>', m.group(0), re.S)
        mats = []
        for s in sets:
            rows = re.findall(r'<v>([^<]*)</v>', s)
            if len(rows) == 3:
                mats.append([[float(x) for x in r.split()] for r in rows])
        if mats:
            born = np.array(mats)
    return eps, born


def parse_outcar(path):
    """Fallback: BORN EFFECTIVE CHARGES block and the dielectric tensor."""
    lines = open(path, errors="ignore").read().split("\n")
    eps = None
    for i, l in enumerate(lines):
        if "MACROSCOPIC STATIC DIELECTRIC TENSOR" in l:
            eps = np.array([[float(x) for x in lines[i + j].split()]
                            for j in (2, 3, 4)])
    born = []
    for i, l in enumerate(lines):
        if "BORN EFFECTIVE CHARGES" in l:
            j = i + 2
            born = []
            while j < len(lines):
                if not lines[j].strip().startswith("ion"):
                    break
                mat = [[float(x) for x in lines[j + k].split()[1:4]]
                       for k in (1, 2, 3)]
                born.append(mat)
                j += 4
    return eps, (np.array(born) if born else None)


def report(eps, born, species=None):
    print("\n=== dielectric tensor (electronic, eps_inf) ===")
    for r in eps:
        print("   " + "  ".join(f"{x:9.5f}" for x in r))
    print(f"   in-plane  xx, yy : {eps[0,0]:.4f}, {eps[1,1]:.4f}"
          "   (bulk h-BN ~4.8-5.1)")
    print(f"   out-of-plane  zz : {eps[2,2]:.4f}"
          "                (bulk h-BN ~2.6-2.9, electronic)")
    aniso = eps[0, 0] / eps[2, 2] if eps[2, 2] else 0
    print(f"   anisotropy xx/zz : {aniso:.2f}"
          "                (bulk h-BN ~1.7-1.8)")
    if aniso < 1.2:
        print("   *** nearly isotropic -- suspicious for a layered material")

    print("\n=== Born effective charges ===")
    for i, z in enumerate(born):
        lab = species[i] if species and i < len(species) else f"atom {i+1}"
        print(f"   {lab}:  xx {z[0,0]:+8.4f}   yy {z[1,1]:+8.4f}"
              f"   zz {z[2,2]:+8.4f}")

    tot = born.sum(axis=0)
    resid = np.abs(tot).max()
    print("\n=== charge neutrality (sum of Z* must vanish) ===")
    for r in tot:
        print("   " + "  ".join(f"{x:+9.5f}" for x in r))
    print(f"   max |residual| = {resid:.5f}")
    if resid < 0.01:
        print("   good -- well below anything that would affect LO-TO")
    elif resid < 0.05:
        print("   acceptable; phonopy will symmetrise it away")
    else:
        print("   *** large. Check ISYM (PEAD usually needs ISYM=0) and the")
        print("       k-point density before trusting the LO-TO splitting.")


def write_born(path, eps, born, dedup=True):
    out = born
    if dedup:
        # keep the first occurrence of each distinct Z* tensor: atoms related
        # by symmetry have identical tensors and phonopy wants only one of each
        keep, seen = [], []
        for z in born:
            if not any(np.allclose(z, s, atol=1e-6) for s in seen):
                seen.append(z)
                keep.append(z)
        out = np.array(keep)
    with open(path, "w") as f:
        f.write(f"{FACTOR}\n")
        f.write(" ".join(f"{x:.8f}" for x in eps.flatten()) + "\n")
        for z in out:
            f.write(" ".join(f"{x:.8f}" for x in z.flatten()) + "\n")
    print(f"\nwrote {path}  ({len(out)} independent atoms "
          f"of {len(born)} total)")
    print("  NOTE: phonopy expects one line per symmetry-INEQUIVALENT atom of")
    print("  the primitive cell. For bulk h-BN (4 atoms, 2 inequivalent) keep")
    print("  the B and N lines only -- drop the duplicates if all 4 are here.")


if __name__ == "__main__":
    d = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else "."
    eps = born = None
    v = os.path.join(d, "vasprun.xml")
    if os.path.exists(v):
        eps, born = parse_vasprun(v)
    if eps is None or born is None:
        o = os.path.join(d, "OUTCAR")
        if os.path.exists(o):
            e2, b2 = parse_outcar(o)
            eps = eps if eps is not None else e2
            born = born if born is not None else b2
    if eps is None or born is None:
        print("could not find the dielectric tensor and/or Born charges.")
        print("check that the run finished and that LCALCEPS or LEPSILON was set:")
        print("  grep -A5 'MACROSCOPIC STATIC DIELECTRIC' OUTCAR")
        print("  grep -A12 'BORN EFFECTIVE CHARGES' OUTCAR")
        sys.exit(1)

    report(eps, born)
    if "--check" not in sys.argv:
        write_born(os.path.join(d, "BORN"), eps, born)
