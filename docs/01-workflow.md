# The staged protocol

Each stage has an exit criterion. Do not start a stage until the previous one
meets it — phonons are second derivatives, and an error in the geometry is
amplified before you ever see a band structure.

Commands below assume you are at the repository root and that `scripts/` is on
your path or invoked as `python scripts/<name>.py`.

---

## Stage 0 — Cells

**Goal.** A primitive cell for each system whose symmetry spglib reports
correctly.

```bash
for s in bulk monolayer bilayer trilayer; do
    echo "== $s"
    phonopy --symmetry --tolerance=1e-4 -c examples/hBN/cells/POSCAR_$s.vasp | head -3
done
```

**Exit criterion.** The space group and operation count match what the point
group demands.

| System | Atoms | Modes | Space group | Ops | Point group | NAC |
|---|---|---|---|---|---|---|
| Bulk | 4 | 12 | `P6_3/mmc` (#194) | 24 | D_6h | **yes** |
| Monolayer | 2 | 6 | `P-6m2` (#187) | 12 | D_3h | no |
| Bilayer | 4 | 12 | `P-3m1` (#164) | 12 | D_3d | no |
| Trilayer | 6 | 18 | `P-6m2` (#187) | 12 | D_3h | no |

**Rules.**

- Hexagonal primitive cell, not orthorhombic
  ([Pitfall 1](02-pitfalls.md#1-orthorhombic-cells-lose-the-3-fold)).
- `SYMPREC` / `--tolerance` never below `1e-4`.
- Vacuum for slabs: ≥ 4× the interlayer spacing. We used 15.157 Å for an
  interlayer of 3.34 Å and verified convergence against 12–28 Å.
- Carry no selective-dynamics flags into any phonon POSCAR.

`scripts/make_cells.py` builds the slab stack from a bulk primitive cell and a
target vacuum.

### If you also need a defect reference

A pristine reference for defect DOS panels must match the defect supercells in
**cell, q-sampling and `SIGMA`** — which generally means computing the pristine
system twice: once in the primitive cell for the dispersion, once in the defect
supercell at Γ only for the DOS comparison.

---

## Stage 1 — Geometry

**Goal.** In-plane stress low enough that the flexural branch is in the
bending-dominated regime at your smallest sampled q.

### Why the threshold is what it is

In-plane tension σ_2D changes the flexural branch from ω ∝ q² to ω ∝ q below

$$q_c = \sqrt{\sigma_{2D}/\kappa}$$

With κ ≈ 0.9 eV per sheet and a smallest sampled q of 0.037 Å⁻¹, anything above
~0.03 kBar puts your first few q-points in the wrong regime. Compute your own
threshold from your own κ and q_min; do not inherit ours.

| System | Target \|σ_xx\| (kBar) | Achieved |
|---|---|---|
| Monolayer | 0.014 | 0.00024 |
| Bilayer | 0.024 | 0.00038 |
| Trilayer | 0.030 | 0.00107 |

### Procedure

| System | Free ionic DOF | Method | INCAR |
|---|---|---|---|
| Bulk | 0 | `ISIF = 3` | `INCAR.relax_bulk` |
| Monolayer | 0 | static only | `INCAR.static` |
| Bilayer | 2 | ions only, `ISIF = 2` | `INCAR.relax_slab` |
| Trilayer | 2 | ions only, `ISIF = 2` | `INCAR.relax_slab` |

Bulk relaxes normally. Slabs do not
([Pitfall 2](02-pitfalls.md#2-isif--3-aborts-on-a-slab-with-vacuum)):

```bash
python scripts/latscan.py gen monolayer   # writes a_*/ directories
# run VASP in each a_*/
python scripts/latscan.py fit monolayer   # quartic fit, centred abscissa

python scripts/stressfit.py               # three statics per system, sigma_xx -> 0
```

`latscan.py` carries per-system constants at the top of the file — the scan
centre, the interlayer spacing and the window width. Edit those for a new
material.

`stressfit.py` validates its own fit against a predicted slope of **a·k/2V**
(the factor 2 because scaling `a` strains the cell biaxially,
[Pitfall 15](02-pitfalls.md#15-the-stress-slope-is-ak2v-not-akv)). A measured
slope far from the prediction means something else is wrong.

**Exit criterion.** |σ_xx| below the threshold for every system, and the three
consistency checks below passing.

### Free consistency checks worth doing

These cost nothing and catch a surprising amount:

- **Additivity.** In-plane stiffness should scale 1 : 2 : 3 with layer count.
  Ours did exactly (38.48 / 77.17 / 115.54 eV/Å²).
- **Interlayer spacing monotonic.** bilayer > trilayer > bulk — fewer
  neighbours, weaker binding, larger spacing. Ours: 3.3597 / 3.3384 / 3.2751 Å.
- **Interlayer binding per interface** roughly constant, with a small excess for
  the trilayer from second-neighbour terms. Ours: −93.4 meV per interface,
  trilayer/bilayer ratio 2.088 against an expected 2.0.

A violation of any of these is a geometry bug, found before you spend a single
displacement run on it.

---

## Stage 2 — Supercell convergence

**Goal.** A supercell large enough that the frequencies you will report are
converged at **commensurate** q-points.

```bash
python scripts/supercell_conv.py status trilayer
python scripts/supercell_conv.py gen    trilayer
python scripts/supercell_conv.py setup  trilayer
# run VASP in every supercell_conv/<system>/*/disp-*/
python scripts/supercell_conv.py fit    trilayer
```

The script is non-destructive by design: `gen` skips directories that already
exist, `setup` skips finished runs, and `fit` discovers every `DIM` directory
including ones you made by hand.

It also runs:

- `commensurate(dim, point)` on every comparison, so M and K are only compared
  where they are exactly representable
  ([Pitfall 4](02-pitfalls.md#4-incommensurate-q-points));
- an explicit minimum-image search for the in-plane force-constant decay, rather
  than rounding fractional coordinates
  ([Pitfall 5](02-pitfalls.md#5-minimum-image-by-rounding-fails-in-a-hexagonal-lattice));
- `gamma_asr_test()`, comparing Γ frequencies with `FC_SYMMETRY` on and off, which
  tells you how much of your acoustic behaviour is the sum rule doing work.

Do not expect the force constants to fall below ~10⁻³ in a polar material
([Pitfall 12](02-pitfalls.md#12-force-constants-do-not-decay-to-zero-in-a-polar-material)).

### KPOINTS

Scale inversely with `DIM`, holding roughly constant spacing (~0.24 Å⁻¹ here).
Files in [`examples/hBN/kpoints/`](../examples/hBN/kpoints/).

| System | Unit cell | 3×3 | 4×4 | 5×5 | 6×6 |
|---|---|---|---|---|---|
| Slabs | 12×12×1 | 4×4×1 | 3×3×1 | 3×3×1 | 2×2×1 |
| Bulk | 12×12×4 | 4×4×2 | 3×3×2 | 3×3×2 | 2×2×2 |

Use the **same** KPOINTS for the ab initio and force-field branches
([Pitfall 7](02-pitfalls.md#7-mismatched-encut-between-stages)).

**Exit criterion.** Frequencies at the q-points you report change by less than
your target precision between the two largest supercells.

Our production sizes: bulk 4×4×2, monolayer 8×8×1, bilayer 6×6×1, trilayer
6×6×1.

---

## Stage 3 — Displacements

**Goal.** Force constants for each branch of the benchmark, differing in exactly
one controlled respect. See
[`04-benchmark-design.md`](04-benchmark-design.md) for why there are three.

```
<system>_ai       AI-relaxed geometry,   AI forces          -> FC_ai
<system>_ml_ai    AI-relaxed geometry,   MLFF forces        -> FC_ml@ai
<system>_ml_ml    MLFF-relaxed geometry, MLFF forces        -> FC_ml
```

```bash
python scripts/make_branch.py <system> <branch> <DIM>   # builds one tree
python scripts/branch_analysis.py status
```

`make_branch.py` builds a `<system>_<branch>/<DIM>/` directory end to end:
supercell, displacements, one `disp-NNN/` per displacement with INCAR, KPOINTS,
POTCAR and `ML_FF` where needed, plus a `band.conf`. It tries both `phonopy -d`
and `phonopy-init` so it works across phonopy 2.x and 4.x.

By hand:

```bash
phonopy -d --dim="6 6 1" -c CONTCAR     # phonopy-init in v4+
# ... run VASP in disp-*/ ...
phonopy -f disp-*/vasprun.xml           # writes FORCE_SETS
```

**Generate the displacement pattern once per geometry** and share the
`POSCAR-{number}` files between `_ai` and `_ml_ai`. That is what makes the
comparison conclusive: the only remaining variable is the force engine.

**Rules.**

- `IBRION = -1`, `NSW = 0`. phonopy applies the displacement; VASP must not move
  the ions. Any ionic step in these OUTCARs invalidates the force constants.
- `EDIFF = 1E-8`
  ([Pitfall 16](02-pitfalls.md#16-ediff--1e-6-is-too-loose-for-force-constants)).
- The MLFF branch is the same INCAR plus `ML_LMLFF = .TRUE.` and
  `ML_MODE = run`.

### MLFF geometries for the `_ml_ml` branch

`_ml_ml` needs the force field's **own** equilibrium. Reusing the ab initio
geometry makes `FC_ml` identical to `FC_ml@ai` and the geometric contribution —
the thing this branch exists to measure — zero by construction.

Repeat Stage 1 with `INCAR.relax_bulk_mlff` and `INCAR.relax_slab_mlff`. Slabs
keep the E(a) scan: `ISIF = 3` fails on a vacuum slab for reasons unrelated to
the force engine. Bulk goes straight to `ISIF = 3`.

The force field has no plane-wave basis and therefore no Pulay stress, so the
E(a) minimum and σ_xx = 0 coincide. Run both anyway — agreement to well under a
mÅ confirms it; disagreement means the MLFF stress tensor needs understanding
before `FC_ml` can be trusted.

Cost is negligible: ~21 MLFF statics plus three relaxations.

---

## Stage 4 — Born charges

Bulk (3D) only. Skip entirely for slabs
([Pitfall 10](02-pitfalls.md#10-q_direction-is-a-conf-file-tag-not-a-command-line-flag)).

```bash
# INCAR.born on the bulk primitive cell, then:
phonopy-vasp-born > BORN.raw
python scripts/make_born.py          # deduplicate, validate, write BORN
python scripts/check_loto.py         # per-direction verification
```

`make_born.py` parses `vasprun.xml`/`OUTCAR`, checks charge neutrality and ε
against expectation, and **keeps only the symmetry-independent atoms**
([Pitfall 9](02-pitfalls.md#9-born-wants-only-the-symmetry-independent-atoms)).

`check_loto.py` writes three conf files — no NAC, q in-plane, q along c — and
reports which modes move under which. Each mode should split only where its
dipole lies along **q**.

**Exit criterion.** Charge neutrality < 0.01; ε in a physically sensible range
for the *electronic* constant
([Pitfall 21](02-pitfalls.md#21-ε_zz-electronic-versus-static)); the
direction-resolved splitting pattern matches the selection rules.

Our values: ε_∞ = 4.802 (in-plane) / 2.749 (out-of-plane), Z*(B) = +2.713,
neutrality 0.00000; in-plane LO–TO 235.33 cm⁻¹ against 243 measured (−3.2%).

---

## Stage 5 — Post-processing and validation

```bash
python scripts/stage5_checks.py                        # whole tree
python scripts/stage5_checks.py --path <dir> <system>  # one band.yaml
python scripts/za_fit.py monolayer_ai/*/
python scripts/irreps.py
```

Conf files in [`examples/hBN/conf/`](../examples/hBN/conf/). `FC_SYMMETRY =
.TRUE.` applies the permutation and translational sum rules; VASP applies
neither. It does not apply rotational invariance
([Pitfall 11](02-pitfalls.md#11-the-acoustic-sum-rule-is-not-rotational-invariance)).

| Check | Criterion | Script |
|---|---|---|
| Branch count | 3 × N_atoms | `stage5_checks.py` |
| Acoustic modes at Γ | \|ω\| < 0.02 meV | `stage5_checks.py` |
| Imaginary modes near Γ | none within 0.05 Å⁻¹ | `stage5_checks.py` |
| Degenerate doublet splitting | < 0.01 meV | `stage5_checks.py` |
| ZA exponent (monolayer) | 1.9 – 2.1 | `za_fit.py` |
| Γ irrep labels | match the point group | `irreps.py` |
| Tolerance sweep 1e-4 → 1e-2 | Γ moves < 0.05 meV | manual |

`za_fit.py` fits ω² = A q⁴ − c and reports the flexural rigidity κ. A nonzero
`c` that survives a doubled FFT grid is the rotational-invariance residual, not
a bug. Ours: κ = 0.83 eV (DFT), 0.955 eV (MLFF, +17%).

For a DOS, use `mesh.conf` (MP 40×40×1 for slabs, 30×30×10 for bulk). For an
in-plane/out-of-plane decomposition, weight eigenvectors by Cartesian component
— phonopy's `PDOS` gives atom projection, not direction projection.

---

## Stage 6 — Analysis, tables and figures

```bash
python scripts/branch_analysis.py gamma   trilayer
python scripts/branch_analysis.py drift   bulk
python scripts/branch_analysis.py compare trilayer
python scripts/phonon_benchmark.py <ref.dat> <cmp.dat>
python scripts/table1.py --latex --cmp ml_ml
python scripts/setup_panels.py
python scripts/plot_dispersions.py
python scripts/plot_defect_dos.py
```

| Comparison | What it isolates |
|---|---|
| `ml_ai` vs `ai` | curvature error at fixed geometry |
| `ml_ml` vs `ml_ai` | geometric contribution |
| `ml_ml` vs `ai` | the error a user of the force field actually sees |

`table1.py` pins the comparison branch to the reference branch's `DIM`
([Pitfall 14](02-pitfalls.md#14-comparing-two-branches-at-different-supercell-sizes))
and prints which it used. `phonon_benchmark.py` gives full-path RMSE and MAE
between two dispersion files.

**Be explicit about the branch in every caption.** A table computed from
`_ml_ai` and a figure computed from `_ml_ml` are answering different questions,
and the mismatch is invisible to a reader.

---

## Triage: inheriting someone else's runs

```bash
bash scripts/phonon_triage.sh <directory>
```

Finds salvageable VASP runs, audits `ENCUT` consistency, and greps the symmetry
and stress lines. Useful before deciding whether to rerun from scratch — which,
in our case, was the right answer.
