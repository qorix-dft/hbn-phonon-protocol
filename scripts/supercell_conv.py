#!/usr/bin/env python3
"""
Stage 2: supercell convergence for the phonon force constants, MLFF only.

Run every command from the WORKING DIRECTORY described below -- both `gen` and
`fit` build paths relative to the current directory.

    workdir/
    |-- supercell_conv.py          <- this script
    |-- cells_zerostress/          <- Stage 1e geometries
    |   |-- POSCAR_monolayer.vasp
    |   |-- POSCAR_bilayer.vasp
    |   |-- POSCAR_trilayer.vasp
    |   `-- POSCAR_bulk.vasp       <- copy your ISIF=3 CONTCAR here
    |-- kpoints/
    `-- supercell_conv/            <- created by `gen`
        `-- trilayer/
            |-- 3x3x1/
            |-- 4x4x1/
            |-- 5x5x1/
            `-- 6x6x1/

Workflow
--------
    python supercell_conv.py status trilayer    # what exists, what has run
    python supercell_conv.py gen    trilayer    # build displacement sets
    python supercell_conv.py setup  trilayer    # make disp-NNN/ dirs + inputs
    # run VASP (MLFF) in every supercell_conv/trilayer/*/disp-*/
    python supercell_conv.py fit    trilayer    # compare commensurate sizes

Omit the system name to act on all four. Add --force to rebuild.

Nothing already on disk is overwritten. `gen` skips any DIM directory that
already contains phonopy_disp.yaml; `setup` skips any disp-NNN/ with a finished
vasprun.xml; `fit` picks up every DIM directory it finds, including sizes you
created by hand that are not in the settings list. So adding a new size to
BULK_DIMS -- or making one yourself -- never disturbs completed runs.

`setup` needs INCAR, POTCAR and the ML_FF file in the working directory; it
copies them into each disp-NNN/ along with the matching KPOINTS.

Slabs use DIM = N N 1 -- never repeat along c, that replicates the vacuum.
Bulk uses N N M and converges the c repetition separately.

Rigid-layer modes converge slowest, because interlayer force constants are
long-ranged relative to their magnitude. Converging on the optical branches
alone is how the low-energy region gets under-converged unnoticed, and that is
exactly where the MLFF benchmark currently fails -- so they are reported in
their own column.
"""
import sys, os, glob, shutil, subprocess, itertools
import numpy as np

# ----------------------------------------------------------------------
# SETTINGS
# ----------------------------------------------------------------------
SLAB_DIMS = [(3, 3, 1), (4, 4, 1), (5, 5, 1), (6, 6, 1)]
BULK_DIMS = [(3, 3, 2), (4, 4, 2), (5, 5, 2), (5, 5, 3), (5, 5, 4), (6, 6, 2)]

CELL_DIR   = "cells_zerostress"
KPT_DIR    = "kpoints"
ROOT       = "supercell_conv"
THZ_TO_CM  = 33.35641

TOL_OPT    = 0.5    # cm^-1, Gamma optical
TOL_ZB     = 1.0    # cm^-1, M and K
TOL_RIGID  = 0.5    # cm^-1, interlayer shear and breathing

NLAYER = {"monolayer": 1, "bilayer": 2, "trilayer": 3, "bulk": 2}

BAND_CONF = """DIM = {a} {b} {c}
PRIMITIVE_AXES = 1 0 0  0 1 0  0 0 1
BAND = 0.0 0.0 0.0  0.5 0.0 0.0  0.3333 0.3333 0.0  0.0 0.0 0.0
BAND_LABELS = $\\Gamma$ M K $\\Gamma$
BAND_POINTS = 41
FC_SYMMETRY = .TRUE.
SYMMETRY = .TRUE.
SYMPREC = 1E-4
EIGENVECTORS = .TRUE.
"""


def dims_for(system):
    return BULK_DIMS if system == "bulk" else SLAB_DIMS


def cell_path(system):
    return os.path.join(CELL_DIR, f"POSCAR_{system}.vasp")


# ----------------------------------------------------------------------
def parse_tag(tag):
    """'5x5x3' -> (5, 5, 3), or None if it is not a DIM directory name."""
    try:
        d = tuple(int(x) for x in tag.split("x"))
        return d if len(d) == 3 and all(v > 0 for v in d) else None
    except ValueError:
        return None


def discover(system):
    """Every DIM directory present on disk, plus any from the settings list
    that have not been created yet. Sorted by supercell size.

    Directories you created by hand are picked up automatically, so nothing
    already run is ignored or overwritten."""
    found = {}
    base = os.path.join(ROOT, system)
    if os.path.isdir(base):
        for tag in os.listdir(base):
            d = parse_tag(tag)
            if d and os.path.isdir(os.path.join(base, tag)):
                found[d] = tag
    for d in dims_for(system):
        found.setdefault(d, "x".join(map(str, d)))
    return [(found[d], d) for d in sorted(found, key=lambda t: (np.prod(t), t))]


def is_built(wd):
    """Has phonopy -d already been run here?"""
    return os.path.exists(os.path.join(wd, "phonopy_disp.yaml"))


def run_done(d):
    """Has VASP finished in this disp directory?"""
    v = os.path.join(d, "vasprun.xml")
    if not os.path.exists(v) or os.path.getsize(v) < 1000:
        return False
    with open(v, "rb") as f:
        f.seek(max(0, os.path.getsize(v) - 400))
        return b"</modeling>" in f.read()


def status(system):
    print(f"\n{system}")
    print(f"{'DIM':>10}{'built':>8}{'disp':>7}{'done':>7}{'FORCE_SETS':>12}"
          f"{'band.yaml':>11}")
    for tag, d in discover(system):
        wd = os.path.join(ROOT, system, tag)
        if not os.path.isdir(wd):
            print(f"{tag:>10}{'-':>8}{'-':>7}{'-':>7}{'-':>12}{'-':>11}")
            continue
        disp = sorted(glob.glob(os.path.join(wd, "disp-[0-9]*")))
        done = sum(run_done(x) for x in disp)
        print(f"{tag:>10}{'yes' if is_built(wd) else 'no':>8}{len(disp):>7}"
              f"{done:>7}"
              f"{'yes' if os.path.exists(os.path.join(wd,'FORCE_SETS')) else 'no':>12}"
              f"{'yes' if os.path.exists(os.path.join(wd,'band.yaml')) else 'no':>11}")


# ----------------------------------------------------------------------
def gen(system, force=False):
    cell = cell_path(system)
    if not os.path.exists(cell):
        print(f"missing {cell}")
        if system == "bulk":
            print("  copy your ISIF=3 CONTCAR to this path")
        return
    for tag, d in discover(system):
        wd = os.path.join(ROOT, system, tag)
        if is_built(wd) and not force:
            n = len(glob.glob(os.path.join(wd, "POSCAR-[0-9]*")))
            done = sum(run_done(x)
                       for x in glob.glob(os.path.join(wd, "disp-[0-9]*")))
            print(f"{system} {tag}: already built ({n} displacements, "
                  f"{done} runs finished) -- skipping")
            continue
        os.makedirs(wd, exist_ok=True)
        shutil.copy(cell, os.path.join(wd, "POSCAR-unitcell"))
        # phonopy v4 moved the -d setup operation to phonopy-init
        for exe in ("phonopy", "phonopy-init"):
            r = subprocess.run([exe, "-d", f"--dim={d[0]} {d[1]} {d[2]}",
                                "-c", "POSCAR-unitcell"],
                               cwd=wd, check=False,
                               capture_output=True, text=True)
            if r.returncode == 0:
                break
            if "phonopy-init" not in (r.stderr or ""):
                print(r.stderr.strip()[:200])
                break
        with open(os.path.join(wd, "band.conf"), "w") as f:
            f.write(BAND_CONF.format(a=d[0], b=d[1], c=d[2]))
        n = len(glob.glob(os.path.join(wd, "POSCAR-[0-9]*")))
        print(f"{system} {tag}: {n} displacements -> {wd}")
    print(f"\nnext:  python {os.path.basename(__file__)} setup {system}")


# ----------------------------------------------------------------------
def setup(system, force=False):
    """Create disp-NNN/ directories with POSCAR, INCAR, KPOINTS, POTCAR.

    Directories with a finished vasprun.xml are left untouched unless
    force=True, so re-running this never disturbs completed work."""
    need = ["INCAR", "POTCAR"]
    missing = [f for f in need if not os.path.exists(f)]
    if missing:
        print(f"missing in working directory: {', '.join(missing)}")
        return
    ml = "ML_FF" if os.path.exists("ML_FF") else None
    if ml is None:
        print("warning: no ML_FF in working directory -- MLFF runs will fail")

    for tag, d in discover(system):
        wd = os.path.join(ROOT, system, tag)
        if not is_built(wd):
            print(f"{system} {tag}: not built yet -- run `gen` first")
            continue
        kpt = os.path.join(KPT_DIR, f"KPOINTS_{system}_{tag}")
        if not os.path.exists(kpt):
            print(f"{system} {tag}: missing {kpt} -- skipping")
            continue
        made = kept = 0
        for p in sorted(glob.glob(os.path.join(wd, "POSCAR-[0-9]*"))):
            num = os.path.basename(p).split("-")[1]
            dd = os.path.join(wd, f"disp-{num}")
            if run_done(dd) and not force:
                kept += 1
                continue
            os.makedirs(dd, exist_ok=True)
            shutil.copy(p, os.path.join(dd, "POSCAR"))
            shutil.copy("INCAR", dd)
            shutil.copy("POTCAR", dd)
            shutil.copy(kpt, os.path.join(dd, "KPOINTS"))
            if ml:
                shutil.copy(ml, dd)
            made += 1
        msg = f"{system} {tag}: {made} disp-NNN/ prepared"
        if kept:
            msg += f", {kept} already finished and left untouched"
        print(msg + f"  (KPOINTS {os.path.basename(kpt)})")
    print("\nrun VASP in each disp-NNN/, then in each supercell directory:")
    print("  phonopy -f disp-*/vasprun.xml")
    print("  phonopy -c POSCAR-unitcell band.conf")
    print(f"\nthen:  python {os.path.basename(__file__)} fit {system}")


# ----------------------------------------------------------------------
def hsp_frequencies(bandyaml):
    """Frequencies (cm^-1) at Gamma, M and K, located from segment_nqpoint."""
    import yaml
    d = yaml.safe_load(open(bandyaml))
    freqs = np.array([[b["frequency"] for b in p["band"]]
                      for p in d["phonon"]]) * THZ_TO_CM
    seg = d.get("segment_nqpoint")
    if not seg:
        n = len(freqs) // 3
        seg = [n, n, n]
    ends = np.cumsum(seg) - 1
    return freqs[0], freqs[ends[0]], freqs[ends[1]]


def gamma_asr_test(wd):
    """Gamma frequencies with and without FC_SYMMETRY.

    FC_SYMMETRY enforces the acoustic sum rule by redistributing the residual
    onto the diagonal. That residual accumulates over pairs, so the correction
    grows with supercell size -- which shifts weak modes (small interlayer
    force constants) while leaving stiff ones alone. If a mode drifts with
    supercell size under symmetrisation but not without it, the drift is the
    redistribution, not convergence, and no larger supercell will fix it.
    """
    try:
        import phonopy
    except ImportError:
        return None
    yml = os.path.join(wd, "phonopy_disp.yaml")
    fsf = os.path.join(wd, "FORCE_SETS")
    if not (os.path.exists(yml) and os.path.exists(fsf)):
        return None
    out = {}
    for sym in (True, False):
        try:
            ph = phonopy.load(phonopy_yaml=yml, force_sets_filename=fsf,
                              produce_fc=True, symmetrize_fc=sym, log_level=0)
            ph.run_qpoints([[0, 0, 0]])
            out[sym] = ph.qpoints.frequencies[0] * THZ_TO_CM
        except Exception:
            return None
    return out


def commensurate(dim, point):
    """Is a high-symmetry point exactly represented by this supercell?

    Finite displacement with DIM = N N M gives exact frequencies only at
    q = (i/N, j/N, k/M). M = (1/2,0,0) needs N even; K = (1/3,1/3,0) needs N
    divisible by 3. Everything else is Fourier interpolation of truncated
    force constants, so comparing an exact value against an interpolated one
    measures the interpolation, not the convergence.
    """
    N = dim[0]
    return {"M": N % 2 == 0, "K": N % 3 == 0}[point]


def fc_decay(wd):
    """Largest |Phi| in the outer shell of the supercell, as a fraction of the
    global maximum. Commensurability-independent, so this is the primary test.

    `phonopy band.conf` builds force constants in memory and only writes
    FORCE_CONSTANTS when asked, so fall back to rebuilding them from FORCE_SETS
    via the phonopy API. Symmetrised to match FC_SYMMETRY in band.conf.
    """
    try:
        import phonopy
    except ImportError:
        return None, "phonopy not importable"

    yml = os.path.join(wd, "phonopy_disp.yaml")
    if not os.path.exists(yml):
        return None, "phonopy_disp.yaml not found"

    kw = dict(phonopy_yaml=yml, produce_fc=True, symmetrize_fc=True, log_level=0)
    fcf, fsf = (os.path.join(wd, "FORCE_CONSTANTS"),
                os.path.join(wd, "FORCE_SETS"))
    if os.path.exists(fcf):
        kw["force_constants_filename"] = fcf
    elif os.path.exists(fsf):
        kw["force_sets_filename"] = fsf
    else:
        return None, "neither FORCE_CONSTANTS nor FORCE_SETS"

    try:
        ph = phonopy.load(**kw)
        fc = ph.force_constants
        sc = ph.supercell
    except Exception as e:
        return None, f"{type(e).__name__}: {e}"
    if fc is None:
        return None, "no force constants produced"

    # Minimum-image distances. Rounding fractional coordinates is INVALID for
    # non-orthogonal cells -- in a hexagonal supercell it overestimates by up
    # to 33 %. Search the neighbouring lattice translations explicitly.
    pos = sc.scaled_positions
    raw = pos[None, :, :] - pos[:, None, :]
    shifts = np.array(list(itertools.product([-1, 0, 1], repeat=3)))
    dist = np.linalg.norm((raw[..., None, :] + shifts) @ sc.cell,
                          axis=-1).min(axis=-1)
    if fc.shape[0] != fc.shape[1]:                 # compact form
        dist = dist[ph.primitive.p2s_map]

    # In-plane separation is the quantity that matters: DIM = N N 1 truncates
    # only in-plane. Using the full 3D distance confounds the test with slab
    # thickness -- r_max comes out as sqrt(rho_max^2 + t^2), so a thick slab's
    # outer shell holds weak interlayer pairs and looks spuriously converged.
    inpl = np.linalg.norm(
        (raw[..., None, :] + shifts) @ sc.cell, axis=-1)
    disp3 = (raw[..., None, :] + shifts) @ sc.cell
    k = np.linalg.norm(disp3, axis=-1).argmin(axis=-1)
    ii, jj = np.indices(k.shape)
    vec = disp3[ii, jj, k]
    rho = np.linalg.norm(vec[..., :2], axis=-1)
    if fc.shape[0] != fc.shape[1]:
        rho = rho[ph.primitive.p2s_map]

    norm = np.linalg.norm(fc.reshape(fc.shape[0], fc.shape[1], 9), axis=-1)
    rmax = rho.max()
    outer = norm[rho > 0.85 * rmax]
    if outer.size == 0:
        return None, "no pairs in the outer shell"

    # power-law exponent of the in-plane tail. hBN is polar, so a long-ranged
    # dipole-dipole tail (|Phi| ~ r^-2.5) is physical and cannot be converged
    # away in real space -- that is what the non-analytic correction handles.
    bins, vals = [], []
    for lo, hi in [(0.30, 0.50), (0.50, 0.70), (0.70, 0.85), (0.85, 1.0)]:
        m = (rho > lo * rmax) & (rho <= hi * rmax)
        if m.any():
            bins.append(0.5 * (lo + hi) * rmax)
            vals.append(norm[m].max() / norm.max())
    expo = (np.polyfit(np.log(bins), np.log(vals), 1)[0]
            if len(bins) >= 3 else float("nan"))
    note = (f"rho_max = {rmax:5.2f} A  tail ~ r^{expo:+.2f}  "
            + " ".join(f"{x:.1e}" for x in vals))

    # Out-of-plane decay. The in-plane metric is structurally blind to the c
    # direction: it takes max |Phi| per in-plane bin, and that maximum is
    # attained by an intra-layer pair whose value does not change when more
    # layers are stacked. Bulk (DIM = N N M) truncates along c as well, so it
    # needs its own tail.
    dz = np.abs(vec[..., 2])
    if fc.shape[0] != fc.shape[1]:
        dz = np.abs(vec[ph.primitive.p2s_map][..., 2])
    zmax = dz.max()
    if zmax > 1.0:
        zb, zv = [], []
        for lo, hi in [(0.30, 0.50), (0.50, 0.70), (0.70, 0.85), (0.85, 1.0)]:
            m = (dz > lo * zmax) & (dz <= hi * zmax)
            if m.any():
                zb.append(0.5 * (lo + hi) * zmax)
                zv.append(norm[m].max() / norm.max())
        zexpo = (np.polyfit(np.log(zb), np.log(zv), 1)[0]
                 if len(zb) >= 3 else None)
        zout = norm[dz > 0.85 * zmax]
        note += (f"\n{'':>16}{zout.max()/norm.max():>11.2e}   "
                 f"dz_max  = {zmax:5.2f} A  tail ~ z^"
                 + (f"{zexpo:+.2f}" if zexpo is not None else " n/a ") + "  "
                 + " ".join(f"{x:.1e}" for x in zv))
    return outer.max() / norm.max(), note


def fit(system):
    nl = NLAYER[system]
    n_rigid = 0 if nl == 1 else 3 * (nl - 1)
    rows = []
    for tag, d in discover(system):
        wd = os.path.join(ROOT, system, tag)
        by = os.path.join(wd, "band.yaml")
        if os.path.exists(by):
            rows.append((tag, d, hsp_frequencies(by), fc_decay(wd)))
    if len(rows) < 2:
        print(f"{system}: only {len(rows)} finished -- need at least 2")
        return

    label = (f"rigid-layer modes = branches 4-{3+n_rigid}"
             if n_rigid else "no rigid-layer modes; watch ZA instead")
    print(f"\n{system}  ({label})")

    print("\n  [1] Gamma -- exact for every supercell, so always comparable")
    print("      compared only within single-axis series: changing the")
    print("      in-plane size and the c repetition at the same time is")
    print("      not a convergence test.")

    def series(rows, key, label):
        groups = {}
        for r in rows:
            groups.setdefault(key(r), []).append(r)
        for gk in sorted(groups):
            grp = sorted(groups[gk], key=lambda r: np.prod(r[1]))
            if len(grp) < 2:
                continue
            print(f"\n      {label(gk)}")
            print(f"{'compare':>18}{'optical':>11}{'rigid':>10}{'verdict':>10}")
            for i in range(len(grp) - 1):
                (t0, _, f0, _), (t1, _, f1, _) = grp[i], grp[i + 1]
                n = min(f0[0].size, f1[0].size)
                dG = np.abs(f1[0][:n] - f0[0][:n])
                opt = dG[3 + n_rigid:].max() if dG.size > 3 + n_rigid else 0.0
                rig = dG[3:3 + n_rigid].max() if n_rigid else 0.0
                ok = opt < TOL_OPT and rig < TOL_RIGID
                print(f"{t0 + ' -> ' + t1:>18}{opt:>11.3f}{rig:>10.3f}"
                      f"{'PASS' if ok else 'FAIL':>10}")

    series(rows, lambda r: r[1][2], lambda k: f"in-plane series, c repetition = {k}")
    if any(r[1][2] > 1 for r in rows) and len({r[1][2] for r in rows}) > 1:
        series(rows, lambda r: r[1][0], lambda k: f"c series, in-plane = {k}x{k}")

    print("\n  [1b] absolute Gamma frequencies (cm^-1) -- read the trend, not")
    print("       just the differences. A monotonic drift is convergence; an")
    print("       oscillation at the 0.1-1 cm^-1 level is symmetrisation noise")
    print("       in weak modes, not something a bigger supercell will fix.")
    nshow = 3 + (n_rigid if n_rigid else 3)
    print(f"{'DIM':>12}" + "".join(f"{'b'+str(i+1):>10}"
                                   for i in range(3, nshow)))
    for tag, _, f, _ in rows:
        print(f"{tag:>12}" + "".join(f"{f[0][i]:>10.3f}"
                                     for i in range(3, min(nshow, f[0].size))))
    if n_rigid:
        arr = np.array([f[0][3:3 + n_rigid] for _, _, f, _ in rows])
        spread = arr.max(axis=0) - arr.min(axis=0)
        print(f"{'spread':>12}" + "".join(f"{x:>10.3f}" for x in spread))

    print("\n  [1c] acoustic sum rule test at Gamma")
    print("       FC_SYMMETRY redistributes the ASR residual onto the diagonal,")
    print("       and that residual grows with supercell size. A mode that")
    print("       drifts with symmetrisation on but is stable with it off is")
    print("       showing the redistribution, not a convergence error.")
    hdr = False
    for tag, _, _, _ in rows:
        g = gamma_asr_test(os.path.join(ROOT, system, tag))
        if g is None:
            continue
        if not hdr:
            print(f"{'DIM':>12}{'sym':>6}" + "".join(
                f"{'b'+str(i+1):>10}" for i in range(3, 3 + max(n_rigid, 3))))
            hdr = True
        for sym in (True, False):
            print(f"{tag if sym else '':>12}{'on' if sym else 'off':>6}"
                  + "".join(f"{g[sym][i]:>10.3f}"
                            for i in range(3, min(3 + max(n_rigid, 3),
                                                  g[sym].size))))

    for pt, idx in (("M", 1), ("K", 2)):
        ok_rows = [r for r in rows if commensurate(r[1], pt)]
        skipped = [r[0] for r in rows if not commensurate(r[1], pt)]
        print(f"\n  [2] {pt} -- commensurate supercells only"
              + (f"  (skipped: {', '.join(skipped)})" if skipped else ""))
        if len(ok_rows) < 2:
            print("      need two commensurate sizes to compare")
            continue
        for i in range(len(ok_rows) - 1):
            (t0, _, f0, _), (t1, _, f1, _) = ok_rows[i], ok_rows[i + 1]
            n = min(f0[idx].size, f1[idx].size)
            dd = np.abs(f1[idx][:n] - f0[idx][:n]).max()
            print(f"{t0 + ' -> ' + t1:>16}{dd:>11.3f}"
                  f"{'PASS' if dd < TOL_ZB else 'FAIL':>20}")

    print("\n  [3] in-plane force-constant decay at the supercell boundary")
    print("      (outer 15 % of IN-PLANE separations, as a fraction of max |Phi|)")
    print("      diagnostic, not pass/fail: hBN is polar, so a long-ranged")
    print("      dipole tail ~r^-2.5 is physical and is handled by NAC.")
    for t, _, _, (dec, note) in rows:
        if dec is None:
            print(f"{t:>16}        --   ({note})")
        else:
            print(f"{t:>16}{dec:>11.2e}   {note}")

    print(f"\n  tolerances: Gamma optical {TOL_OPT}, zone boundary {TOL_ZB}, "
          f"rigid-layer {TOL_RIGID} cm^-1")


# ----------------------------------------------------------------------
if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--force"]
    force = "--force" in sys.argv
    if len(args) < 1 or args[0] not in ("gen", "setup", "fit", "status"):
        print(__doc__)
        sys.exit(1)
    cmd = args[0]
    targets = [args[1]] if len(args) > 1 else list(NLAYER)
    for t in targets:
        if t not in NLAYER:
            print(f"unknown system: {t}")
            continue
        if cmd == "gen":
            gen(t, force)
        elif cmd == "setup":
            setup(t, force)
        elif cmd == "status":
            status(t)
        else:
            fit(t)
