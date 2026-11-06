#!/usr/bin/env python3
"""
Stage 5 validation. Runs on the band.yaml files you already have -- no new VASP.

    python stage5_checks.py                 # every system and branch found
    python stage5_checks.py trilayer        # one system

Checks, in the order they were originally problems:

  branch count        should equal 3 x (atoms in the primitive cell)
  acoustic at Gamma   |omega| < 0.02 meV once FC_SYMMETRY imposes the ASR
  imaginary near Gamma  none within 0.05 1/A
  doublet splitting   symmetry-required degeneracies, < 0.01 meV
  ZA exponent         omega ~ q^n near Gamma, n = 2 for a slab

The ZA check applies to slabs only. In a 3D crystal the out-of-plane acoustic
branch disperses linearly along an in-plane q, so bulk is skipped.
"""
import sys, os, glob
import numpy as np

THZ_TO_CM = 33.35641
CM_TO_MEV = 1.0 / 8.06554
BRANCHES = ["ai", "ml_ai", "ml_ml"]
NLAYER = {"monolayer": 1, "bilayer": 2, "trilayer": 3, "bulk": 2}
NATOM = {"monolayer": 2, "bilayer": 4, "trilayer": 6, "bulk": 4}

TOL_ACOUSTIC = 0.02      # meV
TOL_SPLIT = 0.01         # meV

# phonopy writes band.yaml distances using reciprocal vectors WITHOUT the 2*pi
# factor, so Gamma-M is ~0.23 in those units, not 1.445 1/A. Convert before
# applying any criterion stated in 1/A.
Q_IMAG = 0.05            # 1/A, physical
A_INPLANE = 2.5085       # for the unit conversion


def systems_present():
    out = set()
    order = sorted(BRANCHES, key=len, reverse=True)
    for p in glob.glob("*_*"):
        if not os.path.isdir(p):
            continue
        for b in order:
            if p.endswith("_" + b):
                out.add(p[: -len(b) - 1])
                break
    return sorted(out)


def load_band(path):
    import yaml
    d = yaml.safe_load(open(path))
    f = np.array([[b["frequency"] for b in p["band"]] for p in d["phonon"]])
    q = np.array([p["distance"] for p in d["phonon"]])
    return f * THZ_TO_CM, q


def check(system, branch, tag, path):
    f, q = load_band(path)
    nb = f.shape[1]
    mev = f * CM_TO_MEV
    out = []

    exp_nb = 3 * NATOM[system]
    out.append(("branch count", f"{nb}", nb == exp_nb, f"expect {exp_nb}"))

    ac = mev[0, :3]
    out.append(("acoustic at Gamma", f"{np.abs(ac).max():.4f} meV",
                np.abs(ac).max() < TOL_ACOUSTIC, f"< {TOL_ACOUSTIC}"))

    # phonopy distance -> physical 1/A
    q_phys = q * 2 * np.pi
    near = q_phys < Q_IMAG
    if near.any():
        worst = mev[near].min()
        out.append(("imaginary near Gamma", f"{worst:+.4f} meV",
                    worst > -1e-6, f"none within {Q_IMAG} 1/A"))
    # a tolerance, not zero: the closing Gamma has acoustic modes at ~1e-15
    neg = mev.min(axis=1) < -1e-3
    if neg.any():
        last = q_phys[neg].max()
        out.append(("imaginary extends to", f"{last:.4f} 1/A", last < Q_IMAG,
                    f"< {Q_IMAG}"))

    # Shear doublets only. An N-layer slab has N-1 of them, bulk has 1; the
    # breathing modes that follow are NOT degenerate, so pairing blindly two
    # at a time compares a breathing mode against the ZO branch above it.
    nl = NLAYER[system]
    n_doublet = 1 if system == "bulk" else nl - 1
    if n_doublet > 0:
        worst = 0.0
        for i in range(n_doublet):
            lo, hi = 3 + 2 * i, 4 + 2 * i
            if hi < nb:
                worst = max(worst, abs(mev[0, hi] - mev[0, lo]))
        out.append((f"doublet splitting ({n_doublet})", f"{worst:.5f} meV",
                    worst < TOL_SPLIT, f"< {TOL_SPLIT}"))

    if system != "bulk":
        # fit only where ZA is positive, and say how many points that left
        w, qq = f[1:20, 0], q_phys[1:20]
        m = (w > 0) & (qq > 0)
        if m.sum() >= 4:
            n = np.polyfit(np.log(qq[m]), np.log(w[m]), 1)[0]
            out.append(("ZA exponent", f"{n:.2f} ({m.sum()} pts)",
                        1.9 < n < 2.1, "1.9 - 2.1"))
        else:
            out.append(("ZA exponent", f"only {m.sum()} positive pts", False,
                        "need 4"))

    print(f"\n  {branch:<7}{tag}")
    for name, val, ok, crit in out:
        print(f"    {name:<24}{val:>18}   {'PASS' if ok else 'FAIL':<6}  {crit}")
    if system != "bulk":
        print(f"    {'ZA near Gamma:':<24}", end="")
        for i in range(1, 7):
            print(f"{f[i,0]*CM_TO_MEV:8.3f}", end="")
        print(f"\n    {'  at q (1/A):':<24}", end="")
        for i in range(1, 7):
            print(f"{q_phys[i]:8.4f}", end="")
        print()
    return

    print(f"\n  {branch:<7}{tag}")
    for name, val, ok, crit in out:
        print(f"    {name:<24}{val:>18}   {'PASS' if ok else 'FAIL':<6}  {crit}")


if __name__ == "__main__":
    # --path <dir> <system>: check a single band.yaml anywhere on disk, e.g.
    # one of the supercell_conv/ directories from Stage 2
    if "--path" in sys.argv:
        i = sys.argv.index("--path")
        d, syst = sys.argv[i + 1], sys.argv[i + 2]
        p = os.path.join(d, "band.yaml")
        if not os.path.exists(p):
            print(f"no band.yaml in {d}")
            sys.exit(1)
        print(f"\n{'='*66}\n{syst}   {d}")
        check(syst, "path", os.path.basename(d.rstrip("/")), p)
        sys.exit(0)

    found = systems_present()
    targets = [sys.argv[1]] if len(sys.argv) > 1 else found
    for s in targets:
        if s not in found:
            print(f"unknown system: {s}")
            continue
        print(f"\n{'='*66}\n{s}")
        for b in BRANCHES:
            base = f"{s}_{b}"
            if not os.path.isdir(base):
                continue
            for tag in sorted(os.listdir(base)):
                p = os.path.join(base, tag, "band.yaml")
                if os.path.exists(p):
                    check(s, b, tag, p)
