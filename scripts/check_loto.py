#!/usr/bin/env python3
"""
Check which bulk modes split under the non-analytic correction, and along
which direction.

    python check_loto.py --write     # write the three conf files
    python check_loto.py             # compare, once the runs have finished

Run from bulk_ai/4x4x2/ (or bulk_ml_ml/4x4x2/) with BORN present.

The idea
--------
The LO-TO splitting comes from the macroscopic electric field a mode's dipole
sets up. A mode is pushed up only when its dipole lies ALONG q:

  no NAC             nothing splits; the reference
  q in-plane [100]   E_1u has an in-plane dipole -> one component parallel to
                     q becomes LO. A_2u's dipole is along c, perpendicular to
                     q, so it stays transverse. E_2g is the out-of-phase
                     Davydov partner: the two layers' dipoles cancel, there is
                     no macroscopic polarisation, and it never splits.
  q along c [001]    now A_2u is the longitudinal one and should rise; E_1u
                     becomes purely transverse and should not move.

If that pattern holds, the assignment in the table is right. If E_2g moves
under either direction, something is wrong with BORN or the mode ordering.
"""
import sys, os
import numpy as np

THZ_TO_CM = 33.35641

CONFS = {
    "nonac":   ("", "0.02 0.0 0.0"),
    "inplane": ("NAC = .TRUE.\n", "0.02 0.0 0.0"),
    "alongc":  ("NAC = .TRUE.\n", "0.0 0.0 0.02"),
}

TEMPLATE = """DIM = 4 4 2
PRIMITIVE_AXES = 1 0 0  0 1 0  0 0 1
BAND = 0.0 0.0 0.0  {end}
BAND_POINTS = 11
FC_SYMMETRY = .TRUE.
SYMPREC = 1E-4
{nac}"""

# branch index (0-based) -> label, for the 4-atom bulk cell
LABELS = {3: "shear E_2g", 5: "ZO' B_1g", 6: "ZO A_2u", 7: "ZO2 B_1g",
          8: "in-plane optical", 9: "in-plane optical",
          10: "in-plane optical", 11: "in-plane optical"}


def write():
    for name, (nac, end) in CONFS.items():
        with open(f"loto_{name}.conf", "w") as f:
            f.write(TEMPLATE.format(nac=nac, end=end))
    print("wrote loto_nonac.conf, loto_inplane.conf, loto_alongc.conf\n")
    print("run each, saving band.yaml under a distinct name:\n")
    for name in CONFS:
        print(f"  phonopy -c POSCAR-unitcell loto_{name}.conf"
              f" && mv band.yaml band_{name}.yaml")
    print("\nthen:  python check_loto.py")


def gamma(name):
    import yaml
    p = f"band_{name}.yaml"
    if not os.path.exists(p):
        return None
    d = yaml.safe_load(open(p))
    return np.array([b["frequency"]
                     for b in d["phonon"][0]["band"]]) * THZ_TO_CM


def compare():
    f = {k: gamma(k) for k in CONFS}
    missing = [k for k, v in f.items() if v is None]
    if missing:
        print(f"missing band_{{{','.join(missing)}}}.yaml"
              " -- run `python check_loto.py --write` first")
        return
    n = min(len(v) for v in f.values())
    print(f"{'br':>4}{'no NAC':>12}{'q in-plane':>13}{'q along c':>12}"
          f"{'moves?':>26}")
    for i in range(3, n):
        a, b, c = f["nonac"][i], f["inplane"][i], f["alongc"][i]
        db, dc = b - a, c - a
        tag = []
        if abs(db) > 1.0:
            tag.append("splits for q in-plane")
        if abs(dc) > 1.0:
            tag.append("splits for q along c")
        note = "; ".join(tag) if tag else "no dipole along either"
        print(f"{i+1:>4}{a:>12.2f}{b:>13.2f}{c:>12.2f}   {note}")

    print("\nexpected, if the assignment is right:")
    print("  in-plane optical : one component rises for q in-plane (LO),")
    print("                     nothing moves for q along c")
    print("  ZO A_2u          : rises for q along c, unmoved for q in-plane")
    print("  E_2g, B_1g       : never move -- no macroscopic dipole")


if __name__ == "__main__":
    (write if "--write" in sys.argv else compare)()
