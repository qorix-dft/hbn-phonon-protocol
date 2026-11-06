# Reference: tags, flags and commands

---

## INCAR tags that matter for phonons

One tag per line. VASP's in-line separator is `;`, not whitespace
([Pitfall 17](02-pitfalls.md#17-vasps-in-incar-separator-is--not-whitespace)).

### Shared block

```
PREC    = Accurate
ENCUT   = 650          ! the cutoff your MLFF was trained at; same everywhere
EDIFF   = 1E-8         ! not 1E-6 — force-constant errors are amplified 100x
NELM    = 200
NELMIN  = 5
ALGO    = Normal
ISYM    = 2
SYMPREC = 1E-4         ! never lower
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.      ! reciprocal-space projectors: forces are cleaner
IVDW    = 12           ! DFT-D3(BJ); closes DFPT, see LCALCEPS below
ADDGRID = .TRUE.
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

### Per-task additions

| Task | Tags | File |
|---|---|---|
| Bulk cell relaxation | `IBRION = 2`, `ISIF = 3`, `NSW = 200`, `POTIM = 0.3`, `EDIFFG = -1E-3` | `INCAR.relax_bulk` |
| Static (zero ionic DOF) | `IBRION = -1`, `NSW = 0` | `INCAR.static` |
| Slab, ions only | `IBRION = 1`, `ISIF = 2`, `NSW = 100`, `POTIM = 0.2`, `EDIFFG = -1E-3` | `INCAR.relax_slab` |
| Displacement, ab initio | `IBRION = -1`, `NSW = 0` | `INCAR.disp_ai` |
| Displacement, MLFF | `+ ML_LMLFF = .TRUE.`, `ML_MODE = run` | `INCAR.disp_mlff` |
| Born charges, ε | `LCALCEPS = .TRUE.`, `LPEAD = .TRUE.`, `IBRION = 6`, `NFREE = 2`, `POTIM = 0.015` | `INCAR.born` |

### Tags with sharp edges

| Tag | Note |
|---|---|
| `IBRION = 1` vs `2` | For slabs use `1` (RMM-DIIS). `IBRION = 2` is conjugate gradient, whose line search is ZBRENT — the thing that aborts on a vacuum cell. |
| `ISIF = 3` | Never on a cell containing vacuum. |
| `LEPSILON` | Incompatible with `IVDW`. Use `LCALCEPS` + `LPEAD`. |
| `IBRION = 7` / `8` | DFPT; same incompatibility. |
| `LPHON_READ_FORCE_CONSTANTS` | Makes every other INCAR tag inert. |
| `LPHON_DISPERSION` | VASP's internal dispersion. Works, but gives you no control over symmetrisation or sum rules; prefer phonopy. |
| `LPHON_POLAR` | VASP's NAC equivalent. Same reasoning. |
| `ML_MODE = run` | No SCF, so `OUTCAR` has no `TOTEN` line. Parse `vasprun.xml`. |
| `NGXF`/`NGYF`/`NGZF` | Doubling these is the FFT-grid test for whether a residual is numerical. |
| `ADDGRID` | Keep it on and keep it consistent across branches. |

---

## phonopy

### Version differences

| Task | phonopy 2.x | phonopy 4.x |
|---|---|---|
| Generate displacements | `phonopy -d --dim="4 4 1" -c POSCAR` | `phonopy-init --dim="4 4 1" -c POSCAR` |
| Collect forces | `phonopy -f disp-*/vasprun.xml` | same |

`scripts/make_branch.py` tries both so it works either way.

### Conf-file tags

```
DIM = 4 4 1                     # must be the supercell actually displaced
PRIMITIVE_AXES = 1 0 0  0 1 0  0 0 1
BAND = 0 0 0  0.5 0 0  0.3333333 0.3333333 0  0 0 0
BAND_LABELS = $\Gamma$ M K $\Gamma$
BAND_POINTS = 101
FC_SYMMETRY = .TRUE.            # permutation + translational (ASR). NOT rotational.
SYMMETRY = .TRUE.
SYMPREC = 1E-4
EIGENVECTORS = .TRUE.
NAC = .TRUE.                    # needs a BORN file; 3D systems only
Q_DIRECTION = 1 0 0             # conf tag only — there is no CLI flag
MP = 40 40 1                    # DOS mesh
DOS = .TRUE.
DOS_RANGE = 0 200 0.1
SIGMA = 0.5                     # identical across any panels you compare
FORCE_CONSTANTS = READ          # reads FORCE_CONSTANTS instead of FORCE_SETS
```

### Command cheat-sheet

```bash
# symmetry of a cell
phonopy --symmetry --tolerance=1e-4 -c POSCAR

# displacements, then forces
phonopy -d --dim="6 6 1" -c CONTCAR
phonopy -f disp-*/vasprun.xml

# dispersion
phonopy -p band.conf                      # with plot
phonopy band.conf                          # band.yaml only

# DOS
phonopy -p mesh.conf

# Gamma irreps  (no --readfc unless you actually wrote FORCE_CONSTANTS)
phonopy -c POSCAR-unitcell --irreps="0 0 0" --tolerance=1e-4 --dim="6 6 1"

# Born charges from a finished LCALCEPS run
phonopy-vasp-born > BORN.raw

# export force constants
phonopy --writefc band.conf
```

### `band.yaml` gotchas

| Field | Note |
|---|---|
| `distance` | **No 2π factor.** Multiply by 2π for a physical Å⁻¹ wavevector. |
| `segment_nqpoint` | Use this to find segment boundaries. Do not compute the index of M or K as `2*seg-2` or any other guess. |
| `frequency` | THz by default. 1 THz = 4.13567 meV = 33.3564 cm⁻¹. |

### Unit conversions

| From | To | Multiply by |
|---|---|---|
| THz | meV | 4.135667 |
| THz | cm⁻¹ | 33.35641 |
| meV | cm⁻¹ | 8.065544 |

---

## The `BORN` file

Format: first line is a conversion factor (phonopy's default `14.399652`),
second line is the dielectric tensor as nine numbers, then one nine-number Z*
block per **symmetry-independent** atom, in the order phonopy lists them.

```
14.399652
  4.802  0.000  0.000  0.000  4.802  0.000  0.000  0.000  2.749
  2.713  0.000  0.000  0.000  2.713  0.000  0.000  0.000 -0.730
 -2.713  0.000  0.000  0.000 -2.713  0.000  0.000  0.000  0.730
```

Two lines of Z* for bulk h-BN — one B, one N — not four
([Pitfall 9](02-pitfalls.md#9-born-wants-only-the-symmetry-independent-atoms)).
`scripts/make_born.py` produces this from a finished run.

Sanity checks before using it: Σ Z* = 0 to better than 0.01; ε compared against
the **electronic** constant, not the static one.

---

## Scripts

| Script | Purpose |
|---|---|
| `make_cells.py` | Build slab stacks from a bulk primitive cell and a target vacuum |
| `latscan.py` | E(a) scan replacing `ISIF = 3`; `gen` / `fit`; quartic fit |
| `stressfit.py` | σ_xx(a) → zero-stress lattice constant; validates the a·k/2V slope |
| `supercell_conv.py` | `status` / `gen` / `setup` / `fit`; commensurability, decay, ASR test |
| `make_branch.py` | Build one `<system>_<branch>/<DIM>/` tree end to end |
| `branch_analysis.py` | `status` / `gamma` / `drift` / `compare` across the branch tree |
| `make_born.py` | Parse, validate and deduplicate Born charges into `BORN` |
| `check_loto.py` | Per-direction LO–TO verification |
| `irreps.py` | Γ irreps across the tree, checked against expected label sets |
| `stage5_checks.py` | Branch count, ASR, imaginary modes, doublet splitting, ZA exponent |
| `za_fit.py` | Fit ω² = A q⁴ − c; report κ and the residual offset |
| `table1.py` | Γ-frequency table; `--latex`, `--cmp {ml_ai,ml_ml}` |
| `phonon_benchmark.py` | Full-path RMSE/MAE between two dispersion files |
| `setup_panels.py` | Write `band.conf` and `mesh.conf` into panel directories |
| `plot_dispersions.py` | Dispersion panels with DOS sidebars |
| `plot_defect_dos.py` | Defect DOS panels with energy-window shading |
| `phonon_triage.sh` | Audit an inherited run tree for salvage |

Scripts carrying material-specific constants at the top: `latscan.py`,
`stressfit.py`, `stage5_checks.py`, `table1.py`, `setup_panels.py`,
`plot_dispersions.py`.
