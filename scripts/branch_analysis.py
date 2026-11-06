#!/usr/bin/env python3
"""
Analyse the three-branch phonon tree.

Layout expected in the working directory:

    <system>_ai      /<DIM>/     AI-relaxed geometry,   AI forces
    <system>_ml_ai   /<DIM>/     AI-relaxed geometry,   MLFF forces
    <system>_ml_ml   /<DIM>/     MLFF-relaxed geometry, MLFF forces

  e.g.  bulk_ai/3x3x2/   trilayer_ml_ml/4x4x1/

Each <DIM> directory needs phonopy_disp.yaml and FORCE_SETS. Gamma frequencies
are computed directly from FORCE_SETS via the phonopy API, so band.yaml is not
required; if band.yaml is present the full-path comparison is added too.

Commands
--------
    python branch_analysis.py status                # what exists
    python branch_analysis.py gamma   bulk          # Gamma table, all branches
    python branch_analysis.py drift   bulk          # Stage 2d: c- and in-plane
                                                    #   drift WITHIN one branch
    python branch_analysis.py compare trilayer      # Stage 6a: AI vs ML at
                                                    #   matched DIM

Omit the system name to act on all systems found.
"""
import sys, os, glob, re
import numpy as np

THZ_TO_CM = 33.35641
BRANCHES  = ["ai", "ml_ai", "ml_ml"]
BRANCH_LABEL = {"ai": "AI geom / AI force",
                "ml_ai": "AI geom / ML force",
                "ml_ml": "ML geom / ML force"}

# number of layers in the primitive cell -> rigid-layer branch count
NLAYER = {"monolayer": 1, "bilayer": 2, "trilayer": 3, "bulk": 2}

# reference values for the bulk shear mode, for context in `drift`
BULK_SHEAR_REF = {"experiment": 51.8, "DFT (Table 1)": 51.99,
                  "MLFF 3x3x2": 35.304, "MLFF 5x5x2": 36.164,
                  "MLFF 5x5x3": 35.530, "MLFF 5x5x4": 34.920}


def parse_dim(tag):
    try:
        d = tuple(int(x) for x in tag.split("x"))
        return d if len(d) == 3 and all(v > 0 for v in d) else None
    except ValueError:
        return None


def systems_present():
    out = set()
    # longest branch suffix first: "bulk_ml_ai" ends with "_ai" too, and must
    # not be read as a system called "bulk_ml"
    order = sorted(BRANCHES, key=len, reverse=True)
    for p in glob.glob("*_*"):
        if not os.path.isdir(p):
            continue
        for b in order:
            if p.endswith("_" + b):
                out.add(p[: -len(b) - 1])
                break
    return sorted(out)


def dims_present(system, branch):
    base = f"{system}_{branch}"
    if not os.path.isdir(base):
        return []
    out = []
    for tag in os.listdir(base):
        d = parse_dim(tag)
        if d and os.path.isdir(os.path.join(base, tag)):
            out.append((tag, d))
    return sorted(out, key=lambda t: (np.prod(t[1]), t[1]))


def has_forces(wd):
    return (os.path.exists(os.path.join(wd, "phonopy_disp.yaml"))
            and os.path.exists(os.path.join(wd, "FORCE_SETS")))


def gamma(wd, symmetrize=True):
    """Gamma frequencies in cm^-1, built from FORCE_SETS."""
    try:
        import phonopy
    except ImportError:
        return None
    if not has_forces(wd):
        return None
    try:
        ph = phonopy.load(phonopy_yaml=os.path.join(wd, "phonopy_disp.yaml"),
                          force_sets_filename=os.path.join(wd, "FORCE_SETS"),
                          produce_fc=True, symmetrize_fc=symmetrize,
                          log_level=0)
        ph.run_qpoints([[0, 0, 0]])
        return ph.qpoints.frequencies[0] * THZ_TO_CM
    except Exception as e:
        print(f"    ! {wd}: {type(e).__name__}: {e}")
        return None


def band_path(wd):
    """Full-path frequencies in cm^-1 from band.yaml, or None."""
    by = os.path.join(wd, "band.yaml")
    if not os.path.exists(by):
        return None
    import yaml
    d = yaml.safe_load(open(by))
    return np.array([[b["frequency"] for b in p["band"]]
                     for p in d["phonon"]]) * THZ_TO_CM


def rigid_slice(system):
    nl = NLAYER[system]
    return slice(3, 3 + (0 if nl == 1 else 3 * (nl - 1)))


# ----------------------------------------------------------------------
def cmd_status(systems):
    print(f"{'system':<12}{'branch':<8}{'DIM':>9}{'FORCE_SETS':>12}"
          f"{'band.yaml':>11}{'disp done':>11}")
    for s in systems:
        for b in BRANCHES:
            for tag, _ in dims_present(s, b):
                wd = os.path.join(f"{s}_{b}", tag)
                nd = len(glob.glob(os.path.join(wd, "disp-[0-9]*")))
                done = sum(
                    os.path.exists(os.path.join(x, "vasprun.xml"))
                    for x in glob.glob(os.path.join(wd, "disp-[0-9]*")))
                print(f"{s:<12}{b:<8}{tag:>9}"
                      f"{'yes' if os.path.exists(os.path.join(wd,'FORCE_SETS')) else 'no':>12}"
                      f"{'yes' if os.path.exists(os.path.join(wd,'band.yaml')) else 'no':>11}"
                      f"{f'{done}/{nd}':>11}")


# ----------------------------------------------------------------------
def cmd_gamma(systems):
    for s in systems:
        rs = rigid_slice(s)
        nshow = max(rs.stop, 6)
        print(f"\n{s}   (rigid-layer branches {rs.start+1}-{rs.stop})")
        print(f"{'branch':<8}{'DIM':>9}" +
              "".join(f"{'b'+str(i+1):>10}" for i in range(3, nshow)))
        for b in BRANCHES:
            for tag, _ in dims_present(s, b):
                f = gamma(os.path.join(f"{s}_{b}", tag))
                if f is None:
                    continue
                print(f"{b:<8}{tag:>9}" +
                      "".join(f"{f[i]:>10.3f}"
                              for i in range(3, min(nshow, f.size))))


# ----------------------------------------------------------------------
def cmd_drift(systems):
    """Stage 2d: does the drift appear WITHIN a single branch?

    Compares only along one axis at a time. A drift present in `ai` is
    physical; one present only in the MLFF branches is a force-field artefact
    and belongs in Stage 7, not in supercell convergence.
    """
    for s in systems:
        rs = rigid_slice(s)
        if rs.stop == rs.start:
            print(f"\n{s}: no rigid-layer modes")
            continue
        print(f"\n{s}")
        for b in BRANCHES:
            dims = dims_present(s, b)
            data = []
            for tag, d in dims:
                f = gamma(os.path.join(f"{s}_{b}", tag))
                if f is not None:
                    data.append((tag, d, f))
            if len(data) < 2:
                continue
            print(f"\n  branch {b}  ({BRANCH_LABEL[b]})")

            for axis, name, key in ((2, "c series", lambda d: d[0]),
                                    (0, "in-plane series", lambda d: d[2])):
                groups = {}
                for tag, d, f in data:
                    groups.setdefault(key(d), []).append((tag, d, f))
                for gk, grp in sorted(groups.items()):
                    if len(grp) < 2:
                        continue
                    grp.sort(key=lambda r: r[1][axis])
                    fixed = ("in-plane = %dx%d" % (gk, gk) if axis == 2
                             else "c repetition = %d" % gk)
                    print(f"    {name}, {fixed}")
                    for i in range(len(grp) - 1):
                        (t0, _, f0), (t1, _, f1) = grp[i], grp[i + 1]
                        d0 = f0[rs]
                        d1 = f1[rs][:d0.size]
                        print(f"      {t0} -> {t1}:  "
                              + "  ".join(f"{a:8.3f} -> {c:8.3f} "
                                          f"({c-a:+6.3f})"
                                          for a, c in zip(d0, d1)))

        if s == "bulk":
            print("\n  reference values for the bulk E_2g shear (cm^-1):")
            for k, v in BULK_SHEAR_REF.items():
                print(f"    {k:<16}{v:8.3f}")
            print("\n  If the AI branch is flat in c while the MLFF branches")
            print("  drift, the drift is a force-field artefact -> Stage 7.")
            print("  If AI drifts too, both axes need converging jointly.")


# ----------------------------------------------------------------------
def cmd_compare(systems):
    """Stage 6a: error decomposition at matched DIM.

      ml_ai vs ai      pure curvature error of the force field
      ml_ml vs ml_ai   geometric contribution
      ml_ml vs ai      what the manuscript reports
    """
    pairs = [("ml_ai", "ai", "curvature error"),
             ("ml_ml", "ml_ai", "geometric contribution"),
             ("ml_ml", "ai", "total, as reported")]
    for s in systems:
        rs = rigid_slice(s)
        print(f"\n{s}")
        common = set(t for t, _ in dims_present(s, "ai"))
        for b in BRANCHES[1:]:
            common &= set(t for t, _ in dims_present(s, b))
        if not common:
            print("  no DIM present in all three branches")
            continue
        for tag in sorted(common, key=lambda t: np.prod(parse_dim(t))):
            print(f"\n  DIM {tag}")
            g = {b: gamma(os.path.join(f"{s}_{b}", tag)) for b in BRANCHES}
            p = {b: band_path(os.path.join(f"{s}_{b}", tag)) for b in BRANCHES}
            print(f"{'comparison':>22}{'Gamma rigid':>14}{'Gamma optical':>15}"
                  f"{'path RMSE':>12}")
            for hi, lo, label in pairs:
                if g[hi] is None or g[lo] is None:
                    continue
                n = min(g[hi].size, g[lo].size)
                d = g[hi][:n] - g[lo][:n]
                rig = np.abs(d[rs]).max() if rs.stop > rs.start else 0.0
                opt = np.abs(d[rs.stop:]).max() if d.size > rs.stop else 0.0
                if p[hi] is not None and p[lo] is not None:
                    m = min(p[hi].shape[0], p[lo].shape[0])
                    k = min(p[hi].shape[1], p[lo].shape[1])
                    rmse = np.sqrt(((p[hi][:m, :k] - p[lo][:m, :k]) ** 2).mean())
                    rtxt = f"{rmse:>12.3f}"
                else:
                    rtxt = f"{'--':>12}"
                print(f"{label:>22}{rig:>14.3f}{opt:>15.3f}{rtxt}")
        print("\n  If the rigid-mode error is mostly in the geometric row,")
        print("  the fix is the relaxation, not retraining.")


# ----------------------------------------------------------------------
if __name__ == "__main__":
    cmds = {"status": cmd_status, "gamma": cmd_gamma,
            "drift": cmd_drift, "compare": cmd_compare}
    if len(sys.argv) < 2 or sys.argv[1] not in cmds:
        print(__doc__)
        sys.exit(1)
    found = systems_present()
    if not found:
        print("no <system>_{ai,ml_ai,ml_ml} directories in the current folder")
        sys.exit(1)
    targets = [sys.argv[2]] if len(sys.argv) > 2 else found
    for t in targets:
        if t not in found:
            print(f"unknown system: {t}  (found: {', '.join(found)})")
            continue
    cmds[sys.argv[1]]([t for t in targets if t in found])
