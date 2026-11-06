# hbn-phonon-protocol

A reproducible workflow for computing phonon dispersions of layered materials
with VASP + phonopy, and for benchmarking a machine-learned force field (MLFF)
against the ab initio reference it was trained on.

Hexagonal boron nitride — bulk, monolayer, bilayer and trilayer — is carried
through end to end as a worked example, with every number in
[`examples/hBN/`](examples/hBN/) reproducible from the inputs in this repository.

**The short version of why this exists.** A phonon dispersion is easy to produce
and hard to produce correctly. Most of the failure modes are silent: the
calculation finishes, the bands look plausible, and the error hides in a
degeneracy that should be exact, an acoustic branch with the wrong exponent, or
a comparison between two supercells that do not share a q-point. This repository
records the checks that catch those, and the fourteen specific traps that cost
us time, in [`docs/02-pitfalls.md`](docs/02-pitfalls.md).

---

## Contents

| Path | What it is |
|---|---|
| [`docs/01-workflow.md`](docs/01-workflow.md) | The staged protocol, start to finish |
| [`docs/02-pitfalls.md`](docs/02-pitfalls.md) | **Fourteen traps and how to detect each one** |
| [`docs/03-reference.md`](docs/03-reference.md) | INCAR tags, phonopy flags, command cheat-sheet |
| [`docs/04-benchmark-design.md`](docs/04-benchmark-design.md) | The three-branch scheme for separating curvature from geometry error |
| [`scripts/`](scripts/) | Setup, convergence, validation and plotting tools |
| [`examples/hBN/`](examples/hBN/) | Cells, INCARs, KPOINTS, conf files and the results we obtained |
| [`figures/`](figures/) | The published benchmark figure |

---

## Requirements

- **VASP 6.3+** — `IVDW`, `LCALCEPS`/`LPEAD`, and `ML_LMLFF` for the force-field branch
- **phonopy 2.20+** — `conda install -c conda-forge phonopy`
- **Python 3.9+** — `pip install -r requirements.txt`

A VASP licence and your own `POTCAR` files are assumed; neither is distributed
here.

---

## The workflow in one page

Phonons are second derivatives of the total energy, so every error in the
geometry or the force engine is amplified before you see it. The protocol is
ordered so that each stage is validated before the next one can inherit its
mistakes.

```
Stage 0   Cells            hexagonal primitive, symmetry verified
Stage 1   Geometry         |sigma_xx| below a q-dependent threshold
Stage 2   Supercell        force-constant decay, commensurate comparisons
Stage 3   Displacements    three branches sharing one displacement pattern
Stage 4   Born charges     bulk only; LO-TO verified per direction
Stage 5   Post-process     ASR, degeneracies, branch count, ZA exponent
Stage 6   Analysis         error decomposition, Gamma table, figures
```

### Stage 0 — Cells

Use the **hexagonal primitive cell**, not a Cartesian-friendly orthorhombic one.
The 3-fold rotation is not an integer matrix in a C-centred orthorhombic basis,
so spglib finds 8 operations instead of 24 and the shear doublets split by a few
tenths of a meV — a purely numerical artefact that is easy to mistake for
physics. ([Pitfall 1](docs/02-pitfalls.md#1-orthorhombic-cells-lose-the-3-fold))

```bash
for s in bulk monolayer bilayer trilayer; do
    phonopy --symmetry --tolerance=1e-4 -c examples/hBN/cells/POSCAR_$s.vasp | head -3
done
```

Expect `P6_3/mmc` (#194, 24 ops) for bulk, `P-6m2` (#187) for monolayer and
trilayer, `P-3m1` (#164) for bilayer. **Never set `--tolerance` or `SYMPREC`
below `1e-4`** — the boron sublattice sits ~0.01 mÅ off ideal, enough to drop
spglib into a lower space group.

### Stage 1 — Geometry

```bash
python scripts/latscan.py gen monolayer     # then bilayer, trilayer
# run VASP in every generated directory
python scripts/latscan.py fit monolayer
python scripts/stressfit.py
```

**Do not use `ISIF = 3` on a slab with vacuum.** σ_zz is referenced to a volume
dominated by empty space, the line search cannot bracket, and VASP aborts with
`ZBRENT: fatal error in bracketing`. Scan E(a) instead, fit a **quartic** (a
parabola over a ±0.04 Å window biases the minimum by ~0.6 mÅ), then correct to
σ_xx = 0 with three statics per system.
([Pitfalls 2, 3](docs/02-pitfalls.md))

Why the stress criterion is tight: in-plane tension turns the flexural branch
from ω ∝ q² into ω ∝ q below q_c = √(σ_2D/κ). With κ ≈ 0.9 eV and a smallest
sampled q of 0.037 Å⁻¹, anything above ~0.03 kBar puts your first q-points in
the wrong regime. Targets: **0.014 / 0.024 / 0.030 kBar** for mono/bi/trilayer.

### Stage 2 — Supercell convergence

```bash
python scripts/supercell_conv.py status trilayer
python scripts/supercell_conv.py gen    trilayer
python scripts/supercell_conv.py setup  trilayer
# run VASP in every supercell_conv/<system>/*/disp-*/
python scripts/supercell_conv.py fit    trilayer
```

Nothing on disk is overwritten: `gen` skips built directories, `setup` skips
finished runs, `fit` discovers every DIM directory including hand-made ones.

Two things this script exists to enforce. **Commensurability:** M and K are
exactly representable only in supercells whose DIM makes them so, and comparing
a frequency at an incommensurate q against an interpolated one produced apparent
53 cm⁻¹ errors where the truth was 0.1 cm⁻¹. **Minimum image in a hexagonal
lattice:** rounding fractional coordinates overestimates in-plane separations by
up to 33%, so the decay analysis does an explicit translation search.
([Pitfalls 4, 5](docs/02-pitfalls.md))

Do not expect force constants to decay below ~10⁻³ in a polar material. The
r^−2.5 tail is physical, not noise.

### Stage 3 — Displacements

Generate the displacement pattern **once per geometry** and share the
`POSCAR-{number}` files between the ab initio and MLFF runs. That is what makes
the comparison conclusive: the only remaining variable is the force engine.

```bash
phonopy -d --dim="4 4 1" -c CONTCAR      # phonopy-init in v4+
phonopy -f disp-*/vasprun.xml
```

`IBRION = -1` with `NSW = 0`: phonopy applies the displacement, VASP must not
move the ions. Any ionic step in these OUTCARs means the force constants are
wrong. Use `EDIFF = 1E-8`, not `1E-6` — with a 0.01 Å displacement a force error
δF becomes a force-constant error of δF/0.01, a 100× amplification that lands
hardest on the weak interlayer terms.

### Stage 4 — Born charges (3D systems only)

```bash
phonopy-vasp-born > BORN      # then trim, or:
python scripts/make_born.py
python scripts/check_loto.py
```

**`LCALCEPS` + `LPEAD`, not `LEPSILON`.** Linear response (`IBRION = 7`/`8`)
does not permit a vdW correction, so with `IVDW` active DFPT is closed to you.

`phonopy-vasp-born` writes every atom; phonopy wants only the
symmetry-independent ones and otherwise reports *"Too many atoms in the BORN
file"*. `make_born.py` deduplicates and checks the acoustic sum rule on Z*.
Apply NAC to the bulk only — phonopy's NAC assumes 3D periodicity, and in a true
2D slab the splitting vanishes at Γ and grows linearly in q.
([Pitfalls 8, 9, 10](docs/02-pitfalls.md))

### Stage 5 — Validation

```bash
python scripts/stage5_checks.py                       # all systems, all branches
python scripts/stage5_checks.py --path <dir> <system>  # one band.yaml
python scripts/za_fit.py monolayer_ai/*/
python scripts/irreps.py
```

| Check | Criterion |
|---|---|
| Branch count | 3 × N_atoms |
| Acoustic modes at Γ | \|ω\| < 0.02 meV |
| Imaginary modes near Γ | none within 0.05 Å⁻¹ |
| Degenerate doublet splitting | < 0.01 meV |
| ZA exponent (monolayer) | 1.9 – 2.1 |
| Γ irrep labels | match the point group |
| Symmetry-tolerance sweep | Γ moves < 0.05 meV over 1e-4 → 1e-2 |

`FC_SYMMETRY = .TRUE.` in `band.conf` applies the permutation and translational
(acoustic) sum rules. **VASP applies neither.** It does *not* apply rotational
(Born–Huang) invariance, which is why a residual ZA offset can survive a
perfectly relaxed cell — see [Pitfall 11](docs/02-pitfalls.md).

### Stage 6 — Analysis and figures

```bash
python scripts/branch_analysis.py compare trilayer
python scripts/table1.py --latex
python scripts/setup_panels.py && python scripts/plot_dispersions.py
```

---

## Adapting this to another material

The protocol is material-agnostic; the h-BN numbers are not. What transfers
directly:

- The staged ordering and every check in Stage 5
- The `ISIF = 3` → E(a)-scan substitution for any cell containing vacuum
- Commensurability, minimum-image and ASR handling in the analysis scripts
- The three-branch benchmark design ([`docs/04-benchmark-design.md`](docs/04-benchmark-design.md))

What you must re-derive for your system: the stress threshold (it depends on κ
and your smallest sampled q), the supercell size, the expected point groups and
irrep labels, and whether NAC applies at all.

The scripts read system names from a dictionary at the top of each file.
`latscan.py`, `stressfit.py`, `stage5_checks.py` and `table1.py` are the four
that carry h-BN-specific constants; they are flagged in comments.

---

## Citation

If this protocol or its scripts are useful, please cite the repository
([`CITATION.cff`](CITATION.cff)) along with phonopy:

> A. Togo and I. Tanaka, *First principles phonon calculations in materials
> science*, Scripta Mater. **108**, 1 (2015).

## Licence

MIT, see [`LICENSE`](LICENSE). `POTCAR` files are not redistributable and are
excluded by `.gitignore`.
