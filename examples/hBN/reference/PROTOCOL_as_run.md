# Phonon benchmark — protocol

Living document. Completed work recorded inline. All INCAR blocks are one tag
per line; VASP's in-line separator is `;`, not whitespace.

**Everything runs at ENCUT = 650 eV**, the cutoff the MLFF was trained at.

---

## Settled during diagnosis

| Question | Answer |
|---|---|
| ENCUT consistency | 650 eV both branches — not a source of the softening |
| Vacuum | 15.157 Å (4.6× interlayer) — adequate |
| vdW | trained with `IVDW`; double-counting ruled out on sign grounds |
| Cell symmetry | orthorhombic cell finds 8 ops, **no 3-fold** — use hexagonal |
| Coordinate precision | B atoms 0.011 mÅ off ideal; spglib flips to `Amm2` below tol 1e-4 |
| Bulk c-drift in shear mode | **MLFF artefact** — DFT is flat (Stage 2d) |

| Symptom | Cause | Fixed by |
|---|---|---|
| Shear doublets split 0.42 / 0.24 meV | orthorhombic cell, 3-fold absent | Stage 0 ✅ |
| Imaginary acoustic modes near Γ | unenforced ASR | Stage 5 |
| ZA exponent 1.5 instead of 2 | residual in-plane strain + ASR | Stage 1 ✅ / 5 |
| LO/TO degenerate at Γ | no NAC — correct for 2D | Stage 4 |
| MLFF shear −29 %, breathing −5 % | force field: sliding response | Stage 7 |

---

## Stage 0 — Cells ✅ DONE

**Hexagonal primitive for all four pristine systems.** The orthorhombic cell
exists only to serve the embedding tooling; for a pristine band structure it
costs the 3-fold and doubles the atom count.

| System | Atoms | Branches | Space group | Ops | NAC |
|---|---|---|---|---|---|
| Bulk | 4 | 12 | `P6_3/mmc` (#194) | 24 | **yes** |
| Monolayer | 2 | 6 | `P-6m2` (#187) | 12 | no |
| Bilayer | 4 | 12 | `P-3m1` (#164) | 12 | no |
| Trilayer | 6 | 18 | `P-6m2` (#187) | 12 | no |

The bilayer is $D_{3d}$, not $D_{3h}$ — Γ labels `A_1g A_2u E_g E_u`.

```bash
for s in bulk monolayer bilayer trilayer; do
    phonopy --symmetry --tolerance=1e-4 -c cells_zerostress/POSCAR_$s.vasp | head -3
done
```

**Never set `--tolerance` or `SYMPREC` below 1e-4.** The boron sublattice sits
0.011 mÅ off ideal, enough to drop spglib into `Amm2` (#38).

### The trilayer is computed twice

| Purpose | Panel | Cell | Sampling |
|---|---|---|---|
| Dispersion + converged DOS | (d) | hexagonal primitive, 6 atoms | dense q-mesh |
| Pristine reference for defects | (h) | orthorhombic supercell, 180 atoms | Γ-only |

Panel (h) must match panels (e)–(g) in cell and sampling, with identical
`SIGMA`. Substituting C for B or N leaves the atom count unchanged, so all four
panels carry exactly 540 modes.

---

## Stage 1 — Relax, verify stress ✅ DONE

**Do not use `ISIF = 3` on the slabs.** For a cell containing vacuum, σ_zz is
numerical noise referenced to a volume dominated by empty space; the line search
cannot bracket and VASP aborts with `ZBRENT: fatal error in bracketing`.

| System | Free ionic DOF | Procedure |
|---|---|---|
| Bulk | 0 | `ISIF = 3` |
| Monolayer | **0** | static only |
| Bilayer | 2 | ions only, `ISIF = 2` |
| Trilayer | 2 | ions only, `ISIF = 2` |

### Results

| System | a (Å) | σ_xx (kBar) | Interlayer (Å) |
|---|---|---|---|
| Monolayer | 2.508960 | +0.00024 | — |
| Bilayer | 2.507926 | −0.00038 | 3.3597 |
| Trilayer | 2.507510 | −0.00107 | 3.3384 |
| Bulk | 2.506566 | +0.00782 | 3.2751 (c = 6.55028) |

Interlayer spacing is monotonic — bilayer > trilayer > bulk. Fewer neighbours,
weaker binding, larger spacing. Outer layers buckle, splitting B and N by
1.6 mÅ; the trilayer middle layer stays planar by symmetry. **These numbers
belong in the manuscript.**

Bulk a is +0.10 % against the experimental 2.504 Å; bulk c/2 is **−1.66 %**
against 3.3305 Å. The contraction matters for interpreting the shear mode (see
Stage 2d).

Three independent consistency checks passed, all scaling 1 : 2 : 3 with layer
count — the anharmonic residual, the 650-vs-850 energy shift, and the in-plane
curvature (38.48 / 77.17 / 115.54 eV/Å²). In-plane stiffness is exactly
additive per sheet. Interlayer binding came out at −93.4 meV per interface with
a trilayer/bilayer ratio of 2.088 against an expected 2.0, the 4 % excess being
the middle layer's second neighbour through the D3 three-body term.

### Why the stress criterion is so tight

In-plane tension changes ZA from $q^2$ to $q$ below
$q_c = \sqrt{\sigma_{2D}/\kappa}$. With κ ≈ 0.9 eV per sheet and a smallest
sampled q of 0.037 Å⁻¹, anything above ~0.03 kBar puts the first few q-points
into the tension-dominated regime — which produced the 1.5 exponent. Targets
were 0.014 / 0.024 / 0.030 kBar; achieved 0.001 or better.

E(a) and σ(a) do not share a zero: the scan minimises energy with the basis
varying between points, whereas VASP's stress is the derivative at fixed basis.
The gap is Pulay stress, worth 0.85 mÅ.

### 1a — Bulk relaxation

KPOINTS 12×12×4 Γ-centred.

```
SYSTEM  = hBN bulk relax
PREC    = Accurate
ENCUT   = 650
EDIFF   = 1E-8
EDIFFG  = -1E-3
NELM    = 200
NELMIN  = 5
ALGO    = Normal
IBRION  = 2
ISIF    = 3
NSW     = 200
POTIM   = 0.3
ISYM    = 2
SYMPREC = 1E-4
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.
IVDW    = 12
ADDGRID = .TRUE.
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

### 1b — Slab lattice scan

```bash
python latscan.py gen monolayer     # bilayer, trilayer
python latscan.py fit trilayer
```

KPOINTS 12×12×1 Γ-centred. **Monolayer** (zero ionic DOF, static):

```
SYSTEM  = hBN monolayer latscan
PREC    = Accurate
ENCUT   = 650
EDIFF   = 1E-8
NELM    = 200
NELMIN  = 5
ALGO    = Normal
IBRION  = -1
NSW     = 0
ISYM    = 2
SYMPREC = 1E-4
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.
IVDW    = 12
ADDGRID = .TRUE.
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

**Bilayer and trilayer** (two ionic DOF at fixed cell):

```
SYSTEM  = hBN slab latscan
PREC    = Accurate
ENCUT   = 650
EDIFF   = 1E-8
EDIFFG  = -1E-3
NELM    = 200
NELMIN  = 5
ALGO    = Normal
IBRION  = 1
ISIF    = 2
NSW     = 100
POTIM   = 0.2
ISYM    = 2
SYMPREC = 1E-4
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.
IVDW    = 12
ADDGRID = .TRUE.
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

**`IBRION = 1`, not `2`.** ZBRENT is the conjugate-gradient line minimiser;
RMM-DIIS avoids it by construction.

**Fit a quartic.** Over a ±0.04 Å window E(a) carries cubic anharmonicity — a
parabola biases the minimum by ~0.63 mÅ.

### 1c–1e — Stress scan and final check

Three statics per system at `a_opt ± 1 mÅ`, fit σ_xx(a) linearly, solve for
zero. Slope is $ak/2V$ — the factor 2 because scaling `a` strains the cell
**biaxially**. Predicted 935.5 / 1543.0 / 1961.5 kBar/Å against measured
936.1 / 1540.4 / 1961.3.

```bash
python stressfit.py
```

Stress-check INCAR (also used for the scan points):

```
SYSTEM  = stress check
PREC    = Accurate
ENCUT   = 650
EDIFF   = 1E-8
NELM    = 200
NELMIN  = 5
ALGO    = Normal
IBRION  = -1
NSW     = 0
ISYM    = 2
SYMPREC = 1E-4
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.
IVDW    = 12
ADDGRID = .TRUE.
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

Ignore σ_zz for the slabs — a basis-set artefact along c. Carry no
selective-dynamics flags into any phonon POSCAR.

---

## Stage 2 — Supercell convergence ✅ DONE

```bash
python supercell_conv.py status trilayer
python supercell_conv.py gen    trilayer
python supercell_conv.py setup  trilayer
# run VASP in every supercell_conv/<system>/*/disp-*/
python supercell_conv.py fit    trilayer
```

Nothing on disk is overwritten: `gen` skips built directories, `setup` skips
finished runs, `fit` discovers every DIM directory including hand-made ones.

### Three lessons the data forced

**Commensurability.** `DIM = N N M` gives exact frequencies only at
q = (i/N, j/N, k/M). M = (½,0,0) needs N even; K = (⅓,⅓,0) needs N divisible
by 3. Comparing exact against interpolated produced apparent 53 cm⁻¹ errors
where the true figure was 0.1.

**Single-axis series.** Changing the in-plane size and the c repetition
together is not a convergence test.

**In-plane separation.** `DIM = N N 1` truncates only in-plane. Using the full
3D distance for the decay confounds it with slab thickness — r_max comes out as
$\sqrt{\rho_{max}^2 + t^2}$, so a thick slab's outer shell holds weak
interlayer pairs and looks spuriously converged.

### Results

| System | Γ optical | rigid | M (4→6) | K (3→6) | Decision |
|---|---|---|---|---|---|
| Monolayer | 0.024 | — | 0.063 | 0.103 | **4×4×1** |
| Bilayer | 0.026 | 0.006 | 0.065 | 0.101 | **4×4×1** |
| Trilayer | 0.023 | 0.004 | 0.056 | 0.084 | **4×4×1** |
| Bulk | 0.027 | 0.132 | 0.077 | 0.127 | **6×6×2** |

All three slabs agree between 3×3×1 and 6×6×1 to ~0.1 cm⁻¹ at the zone
boundary. 4×4×1 gives margin at 96 atoms for the trilayer against 216 at 6×6×1.

The in-plane force constants show a long-ranged tail, |Φ| ∝ r^−2.5 to r^−3.3.
That is the polar dipole–dipole interaction, physical and not convergible in
real space — it is what NAC subtracts analytically. Not a truncation error.

### 2d — Bulk shear mode ✅ RESOLVED: force-field artefact

The bulk E₂g shear drifted with c repetition in the MLFF branch. Three
explanations were tested and two eliminated:

- **Real-space truncation along c** — force constants reach 1.2×10⁻⁶ of maximum
  by 13 Å, tail z^−7.4. Far too small to move a frequency by 0.6 cm⁻¹.
- **ASR redistribution** — Γ frequencies identical to three decimals with
  `FC_SYMMETRY` on and off, at every supercell size.
- **The force field itself** — confirmed by running the AI branch on matched
  supercells.

| 3×3×2 → 3×3×3 | E₂g shear |
|---|---|
| AI | 48.559 → 48.696 (**+0.137**) |
| MLFF | 35.304 → 34.469 (**−0.835**) |

DFT is flat in c. **The bulk supercell is fixed from the in-plane series at
6×6×2**, and the drift moves to Stage 7.

### Matched-supercell errors, bulk rigid modes

Same geometry, same displacements, same supercell — only the force engine
differs.

| Mode | AI | MLFF | Error |
|---|---|---|---|
| E₂g shear (3×3×2) | 48.559 | 35.304 | **−27.3 %** |
| E₂g shear (3×3×3) | 48.696 | 34.469 | **−29.2 %** |
| B₁g breathing | 126.22 | 119.71 | −5.2 % |

The force field fails specifically at the **sliding** response, not at
interlayer interaction generally.

### The DFT shear also sits below experiment

Converged AI estimate is 48.696 + 1.06 (in-plane correction from the 3×3→5×5
series) ≈ **49.76 cm⁻¹** against 51.8 measured, i.e. −3.9 %.

Bulk c/2 is 1.66 % *below* experiment, and a compressed interlayer spacing means
more registry corrugation, which should make the shear mode *stiffer*. The
functional therefore underestimates the sliding barrier by more than 3.9 %.
This is the expected signature of D3: a pairwise isotropic dispersion correction
reproduces binding energy and layer spacing but carries no registry dependence,
while E₂g is set by exchange and electrostatic corrugation that D3 does not
touch.

The two errors are independent and should be reported separately:

| Comparison | Error | Cause |
|---|---|---|
| DFT-D3 vs experiment | −3.9 % | functional: no registry dependence in D3 |
| MLFF vs its own DFT | −29.2 % | force field, Stage 7 |

**Do not reconcile against the old Table 1 value of 51.99.** That came from the
orthorhombic cell with mm2 symmetrisation, `EDIFF = 1E-6` and no ASR. Its
agreement with experiment to 0.2 cm⁻¹ was almost certainly fortuitous.
Regenerate the row.

### KPOINTS

Scaled inversely with `DIM`, holding ~0.24 Å⁻¹ spacing. Files in `kpoints/`.

| System | Unit cell | 3×3 | 4×4 | 5×5 | 6×6 |
|---|---|---|---|---|---|
| Slabs | 12×12×1 | 4×4×1 | 3×3×1 | 3×3×1 | 2×2×1 |
| Bulk | 12×12×4 | 4×4×2 | 3×3×2 | 3×3×2 | 2×2×2 |

Bulk 5×5×4 → 3×3×1; bulk 6×6×2 → 2×2×2. Use the same KPOINTS for both
branches — the MLFF does no SCF, and identical files remove one more difference
between `FC_ai` and `FC_ml`.

---

## Stage 3 — Displacements: three branches ⬜ IN PROGRESS

Directory convention:

```
<system>_ai       AI-relaxed geometry,   AI forces
<system>_ml_ai    AI-relaxed geometry,   MLFF forces
<system>_ml_ml    MLFF-relaxed geometry, MLFF forces
```

```bash
python branch_analysis.py status
python branch_analysis.py gamma   trilayer
python branch_analysis.py drift   bulk
python branch_analysis.py compare trilayer
```

Generate displacements **once** per geometry and share the `POSCAR-{number}`
files between `_ai` and `_ml_ai`. That is what makes the benchmark conclusive:
the only variable becomes the force engine.

```bash
phonopy -d --dim="4 4 1" -c CONTCAR      # 6 6 2 for bulk
```

**Ab initio branch:**

```
SYSTEM  = phonon displacement AI
PREC    = Accurate
ENCUT   = 650
EDIFF   = 1E-8
NELM    = 200
NELMIN  = 5
ALGO    = Normal
IBRION  = -1
NSW     = 0
ISYM    = 2
SYMPREC = 1E-4
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.
IVDW    = 12
ADDGRID = .TRUE.
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

**MLFF branch** — identical plus:

```
ML_LMLFF = .TRUE.
ML_MODE  = run
```

`IBRION = -1` with `NSW = 0`: phonopy applies the displacement, VASP must not
move the ions. Any ionic step in these OUTCARs means the force constants are
wrong.

`EDIFF = 1E-8`, not 1E-6. With a 0.01 Å displacement a force error δF becomes a
force-constant error of δF/0.01 — a 100× amplification, landing hardest on the
interlayer terms.

```bash
phonopy -f disp-*/vasprun.xml
```

### Status

| System | `_ai` | `_ml_ai` | `_ml_ml` |
|---|---|---|---|
| Monolayer | ✅ | ⬜ | geometry ✅ |
| Bilayer | ✅ | ⬜ | geometry ✅ |
| Trilayer | ✅ | ⬜ | geometry ✅ |
| Bulk | 3×3×2, 3×3×3 ✅ | from Stage 2 | geometry ✅ |

`_ml_ai` at the production DIM can be copied from the Stage 2 tree, which used
the same AI-relaxed geometries.

---

### 3b — MLFF geometries for the `_ml_ml` branch ✅ DONE

`_ml_ml` requires the force field's **own** equilibrium. Reusing the AI geometry
makes `FC_ml` identical to `FC_ml@ai`, and the geometric contribution — the
thing this branch exists to measure — becomes zero by construction.

**The key simplification: the MLFF has no plane-wave basis, so no Pulay
stress.** In Stage 1 the E(a) minimum and the σ_xx = 0 point differed by
0.85 mÅ, which is why the scan needed a separate stress correction. For the
force field they coincide, and the two-step procedure collapses to one.

Run both anyway as a cross-check: agreement to well under a mÅ confirms the
absence of Pulay stress in this branch. Disagreement means the MLFF stress
tensor is not what it should be, and `FC_ml` cannot be trusted until that is
understood.

**Slabs keep the scan.** `ISIF = 3` fails on a vacuum slab for reasons
unrelated to the force engine — VASP tries to relax c and collapses the vacuum.
**Bulk goes straight to `ISIF = 3`.**

Cost: 21 MLFF statics for the three scans, two ionic relaxations, one bulk
relaxation. Effectively free.

#### 3b-i — Bulk, MLFF relaxation

KPOINTS 12×12×4. Same as Stage 1a plus the two MLFF tags:

```
SYSTEM   = hBN bulk relax MLFF
PREC     = Accurate
ENCUT    = 650
EDIFF    = 1E-8
EDIFFG   = -1E-3
NELM     = 200
NELMIN   = 5
ALGO     = Normal
IBRION   = 2
ISIF     = 3
NSW      = 200
POTIM    = 0.3
ISYM     = 2
SYMPREC  = 1E-4
ISMEAR   = 0
SIGMA    = 0.01
LREAL    = .FALSE.
IVDW     = 12
ADDGRID  = .TRUE.
LWAVE    = .FALSE.
LCHARG   = .FALSE.
ML_LMLFF = .TRUE.
ML_MODE  = run
```

Compare the resulting c/2 against the AI value of 3.2751 Å. A significant
difference is another Stage 7 data point — the B₁g breathing error is −5 %, so
the force field's interlayer spacing is expected to be close but not identical.

#### 3b-ii — Slabs, MLFF lattice scan

```bash
python latscan.py gen monolayer     # bilayer, trilayer
# run with the MLFF INCAR below
python latscan.py fit monolayer
```

Widen `A_SPAN` in `latscan.py` from 0.04 to 0.06 Å first: the MLFF minimum is
not guaranteed to sit inside the AI bracket.

**Monolayer** (zero ionic DOF, static):

```
SYSTEM   = hBN monolayer latscan MLFF
PREC     = Accurate
ENCUT    = 650
EDIFF    = 1E-8
NELM     = 200
NELMIN   = 5
ALGO     = Normal
IBRION   = -1
NSW      = 0
ISYM     = 2
SYMPREC  = 1E-4
ISMEAR   = 0
SIGMA    = 0.01
LREAL    = .FALSE.
IVDW     = 12
ADDGRID  = .TRUE.
LWAVE    = .FALSE.
LCHARG   = .FALSE.
ML_LMLFF = .TRUE.
ML_MODE  = run
```

**Bilayer and trilayer** (two ionic DOF at fixed cell):

```
SYSTEM   = hBN slab latscan MLFF
PREC     = Accurate
ENCUT    = 650
EDIFF    = 1E-8
EDIFFG   = -1E-3
NELM     = 200
NELMIN   = 5
ALGO     = Normal
IBRION   = 1
ISIF     = 2
NSW      = 100
POTIM    = 0.2
ISYM     = 2
SYMPREC  = 1E-4
ISMEAR   = 0
SIGMA    = 0.01
LREAL    = .FALSE.
IVDW     = 12
ADDGRID  = .TRUE.
LWAVE    = .FALSE.
LCHARG   = .FALSE.
ML_LMLFF = .TRUE.
ML_MODE  = run
```

Still `IBRION = 1`: ZBRENT is the CG line minimiser and the failure mode is
independent of the force engine.

#### 3b-iii — Stress cross-check

Three statics per system at `a_ml ± 1 mÅ` using the Stage 1c stress INCAR plus
the MLFF tags, then:

```bash
python stressfit.py
```

**Criterion:** a(σ=0) agrees with the quartic E(a) minimum to < 0.1 mÅ. That is
the no-Pulay-stress confirmation. Predicted slopes are unchanged in form,
$ak/2V$, but use the MLFF's own fitted curvature rather than the AI values of
38.48 / 77.17 / 115.54 eV/Å².

#### 3b-iv — Final geometries ✅ DONE

The energy scan gave a_ml; a stress check then required a further −0.021 /
−0.030 / −0.006 mÅ, applied and re-verified. Final σ_xx is 0.0001–0.0008 kBar
against tolerances of 0.014–0.030, and the tension share of ω² at the smallest
sampled q fell from 15–22 % to under 0.4 %. **The ZA branch is now clean in
both branches.**

| System | a (Å) | Interlayer (Å) | Buckling | σ_xx (kBar) |
|---|---|---|---|---|
| Monolayer | 2.510369 | — | 0 | +0.00008 |
| Bilayer | 2.507180 | 3.3333 | 1.50 mÅ | +0.00079 |
| Trilayer | 2.507214 | 3.3163 | 1.47 mÅ | −0.00033 |
| Bulk | 2.507056 | 3.2842 | — | (`ISIF=3`) |

**σ_zz is the no-Pulay-stress confirmation.** AI branch: −0.743 / −1.246 /
−1.527 kBar. MLFF branch: 0.000 / −0.003 / −0.001. Three orders of magnitude
smaller. With no plane-wave basis, σ_zz for a vacuum slab is a clean derivative
and comes out at zero — direct evidence that the σ_zz you are told to ignore in
the AI branch really is an artefact.

Buckling matches the AI branch closely (1.50 / 1.47 mÅ against 1.6 mÅ), so the
force field gets the surface relaxation right.

#### Comparison with the AI branch

| Quantity | AI | MLFF | Difference |
|---|---|---|---|
| a monolayer | 2.508960 | 2.510369 | +1.41 mÅ |
| a bilayer | 2.507926 | 2.507180 | −0.75 mÅ |
| a trilayer | 2.507510 | 2.507214 | −0.30 mÅ |
| a bulk | 2.506566 | 2.507056 | +0.49 mÅ |
| interlayer bilayer | 3.3597 | 3.3333 | **−26.4 mÅ** |
| interlayer trilayer | 3.3384 | 3.3163 | **−22.1 mÅ** |
| interlayer bulk | 3.27514 | 3.28416 | +9.0 mÅ |
| in-plane curvature | 38.48 / 77.17 / 115.54 | 39.03 / 79.09 / 119.33 | +1.4 / +2.5 / +3.3 % |

**Three systematic findings, all pointing the same way.**

*The MLFF flattens the layer-count trend.* Successive differences in a:

| | mono→bi | bi→tri | tri→bulk |
|---|---|---|---|
| AI | −1.03 | −0.42 | −0.94 mÅ |
| MLFF | −3.18 | **+0.01** | −0.16 mÅ |

DFT contracts monotonically as layers are added. The force field takes one large
step from monolayer to bilayer and then goes flat — it registers "isolated sheet
versus stacked" but not the continued evolution with stacking.

*In-plane stiffness is no longer additive.*

| | AI | MLFF |
|---|---|---|
| curvature ratio | 1 : 2.005 : 3.003 | 1 : 2.026 : **3.057** |

Van der Waals contributes nothing to in-plane elasticity, so additivity is exact
by construction and any deviation is unambiguous. The MLFF is 1.9 % non-additive
at three layers — interlayer coupling leaking into the in-plane elastic response
where it should not appear at all.

*The interlayer spacing trend is compressed.* The sign of the interlayer error
flips between slabs (−22 to −26 mÅ) and bulk (+9 mÅ), so it is not a uniform
binding offset. The clearer statement is the spread across the sequence:

| | bilayer → bulk |
|---|---|
| AI | 84.6 mÅ |
| MLFF | 49.2 mÅ |

**The force field compresses the layer-count dependence to 58 % of the DFT
value.** It gets "few layers" versus "bulk" roughly right but flattens
everything in between.

Three independent measurements — in-plane lattice constants, in-plane curvature
additivity, and interlayer spacing — now say the same thing: the force field
does not resolve how the interlayer environment evolves with stacking. That is
the same defect as the −29 % shear error, seen from three different angles. All
belong in Stage 7b.

**Consequence for Stage 6a.** The trilayer interlayer error is 22 mÅ, roughly
8× the bulk offset, so the geometric contribution there will be larger than the
bulk numbers alone suggested. Whether it is large enough to matter against the
−29 % curvature error is exactly what `FC_ml` vs `FC_ml@ai` measures.

Final cells are in `cells_ml_corrected/`.

---

## Stage 4 — Born charges (bulk only) DONE

`LCALCEPS` + `LPEAD` + `IBRION = 6` on the bulk unit cell.

| Quantity | Value | Expected |
|---|---|---|
| eps_inf in-plane | 4.802 | 4.8-5.1 |
| eps_inf out-of-plane | 2.749 | 2.6-2.9 (electronic) |
| Z*(B) in-plane | +2.713 | ~+2.7 |
| charge neutrality | 0.00000 | < 0.01 |

The out-of-plane value looks low against Geick's 4.10, but that is the STATIC
constant; the electronic one is smaller and the modern static figure is
3.4 +/- 0.2. The 4.802 / 2.749 anisotropy is h-BN's hyperbolic behaviour.

**`phonopy-vasp-born` writes all four atoms; phonopy wants only the two
symmetry-independent ones.** Trim to one B and one N, or phonopy reports "Too
many atoms in the BORN file". `make_born.py` deduplicates automatically.

**The phonopy API path for NAC was unreliable** -- it reported "NAC on" while
returning unsplit frequencies. Use a conf file with `NAC = .TRUE.` and a short
band segment off Gamma; that worked every time.

### LO-TO, verified by direction

`check_loto.py` runs Gamma three ways -- no NAC, q in-plane, q along c -- and
confirms each mode splits only when its dipole lies along q:

| Mode | q in-plane | q along c |
|---|---|---|
| E_1u in-plane optical | 1356.03 -> **1591.36** | unmoved |
| A_2u ZO | unmoved at 741.42 | -> **805.01** |
| E_2g, B_1g | unmoved | unmoved |

E_2g is the out-of-phase Davydov partner: the two layers' dipoles cancel, so
there is no macroscopic polarisation and it never splits. On the Gamma-M-K-Gamma
path the A_2u dipole is always perpendicular to q, so panel (a) shows its TO
throughout.

In-plane splitting 235.33 cm^-1 against 243 measured (-3.2 %). Out-of-plane
63.6 against Geick's 45 (+40 %), consistent with eps_zz being the electronic
value. Table 1 reports the TO, so this does not affect it -- but do not quote
the A_2u LO without a caveat.

---

## Stage 4 — original instructions 

**`LCALCEPS`, not `LEPSILON`.** Linear response (`IBRION = 7` or `8`) does not
permit a vdW correction; with `IVDW = 12` active, DFPT is closed to you.

```
SYSTEM   = Born charges
PREC     = Accurate
ENCUT    = 650
EDIFF    = 1E-8
NELM     = 200
NELMIN   = 5
ALGO     = Normal
LCALCEPS = .TRUE.
LPEAD    = .TRUE.
IBRION   = 6
NFREE    = 2
POTIM    = 0.015
ISYM     = 2
SYMPREC  = 1E-4
ISMEAR   = 0
SIGMA    = 0.01
LREAL    = .FALSE.
IVDW     = 12
LWAVE    = .FALSE.
LCHARG   = .FALSE.
```

```bash
phonopy-vasp-born > BORN
```

---

## Stage 5 — Post-processing ✅ MOSTLY DONE

```bash
python stage5_checks.py                       # all systems and branches
python stage5_checks.py --path <dir> <system> # any single band.yaml
python za_fit.py monolayer_ai/*/ monolayer_ml_ml/*/
```

### Results

**Doublet splitting is zero everywhere**, all systems and branches. The Stage 0
problem (0.42 / 0.24 meV of spurious splitting from the orthorhombic cell) is
gone.

**Acoustic modes at Γ are 0.0000 meV** everywhere — `FC_SYMMETRY` imposes the
ASR exactly.

**ZA needed a bigger supercell than everything else.** At 4×4×1, where every
Stage 2 criterion passed, the flexural branch was imaginary out to ~0.2 Å⁻¹ in
all systems. Commensurate-point convergence (Γ, M, K — exact by construction)
does not guarantee correct Fourier interpolation at incommensurate small q, and
ZA is the branch where that distinction bites.

| Monolayer first q-point | `_ai` | `_ml_ml` |
|---|---|---|
| 4×4×1 | −0.368 | −0.246 |
| 6×6×1 | −0.093 | −0.099 |
| 8×8×1 | −0.103 | +0.012 |
| 10×10×1 | — | −0.011 |

Bulk showed the same threshold (imaginary at 3×3 and 4×4, clean at 5×5 and
6×6) — and bulk has no flexural branch at all, which is what identified this as
interpolation rather than a sum-rule problem.

**Fit to ω² = A q⁴ − c**, first 8 q-points:

| Branch | A converged | κ | c |
|---|---|---|---|
| `_ai` | 2 % spread | **0.82 eV** | holds at ~0.10 |
| `_ml_ml` | 4 % spread | 0.955 eV | collapses 10× |

The MLFF is rotationally invariant by construction (its energy depends only on
distances and angles), so its force constants satisfy Born–Huang automatically
and c converges away. The DFT branch does not. Remaining candidates: a genuine
rotational-invariance violation, or integration noise mimicking one. **Test:**
one pair of displacements at 6×6×1 with `NGXF/NGYF/NGZF` doubled. If c and the
residual both fall toward the MLFF values it was the augmentation grid; if both
hold, it is the sum rule.

κ = 0.82 eV sits in the accepted 0.8–0.9 eV range for monolayer h-BN, so the
branch is physically right despite the small-q residual — 0.1 meV² on a
spectrum reaching 180 meV.

### The ZA exponent criterion applies to the monolayer only

A stack of coupled layers is 2D-flexural only below a crossover where the
bending energy κq⁴ matches the interlayer coupling; above it the branch becomes
an ordinary 3D acoustic mode, linear in q. So the exponent must fall toward 1
as layers are added:

| System | `_ai` | `_ml_ml` |
|---|---|---|
| Monolayer | pending | 2.10 |
| Bilayer | **1.33** | **1.34** |
| Trilayer | pending | 1.19 |
| Bulk | 1 | 1 (3D acoustic) |

Bilayer agrees between branches to 0.8 %, so the crossover is physical, not a
method artefact. **This is a result worth stating in the manuscript**, not a
failed check.

### MLFF out-of-plane stiffness

| Observable | MLFF error |
|---|---|
| κ, monolayer | **+17 %** |
| ZA magnitude, bilayer | +9.6 % |
| E₂g shear, bilayer | +11.9 % |
| E₂g shear, trilayer | +21.6 / +17.5 % |
| E₂g shear, bulk | **−26.6 %** |

κ and ZA magnitude are consistent: ω ∝ √A, so +17 % in κ predicts +8.2 % in
frequency against +9.6 % measured. Same overstiff-out-of-plane defect as the
shear modes in slabs, flipping sign in bulk.

---

## Stage 5 — remaining items DONE

**Irreps** confirmed all four point groups: bulk D_6h, monolayer D_3h, bilayer
D_3d, trilayer D_3h. The monolayer returns D_3h, so the draft's D_6h labels
were wrong labels on a correct calculation.

The invalid label in the draft does not exist in any point group. The trilayer
shear doublets are E'' (22.05) and E' (37.97) -- two different irreps -- and the
TO/LO manifold is E', E'', E'. The trilayer ZO is a Davydov triplet
(791.56 / 798.80 / 801.74), not the single entry the draft carried.

Gamma decompositions:

| System | |
|---|---|
| Bulk | 2A_2u + 2B_1g + 2E_1u + 2E_2g |
| Monolayer | 2A_2'' + 2E' |
| Bilayer | 2A_1g + 2A_2u + 2E_g + 2E_u |
| Trilayer | 2A_1' + 4A_2'' + 4E' + 2E'' |

phonopy reports "Not found" for the acoustic triplet at 0 -- expected, the three
translations are degenerate and cannot be split. Assign by symmetry:
D_6h -> A_2u + E_1u; D_3h -> A_2'' + E'; D_3d -> A_2u + E_u.

### The ZA offset is the rotational sum rule -- RESOLVED

Three candidates, two eliminated by measurement:

| Candidate | Test | Result |
|---|---|---|
| supercell truncation | 6x6x1 -> 8x8x1 | MLFF c collapses 92 %, **DFT holds** |
| integration noise | FFT grid doubled at 6x6x1 | A stable 0.9 %, **c falls 10 %** |
| rotational invariance | -- | by elimination |

Fit to omega^2 = A q^4 - c at 6x6x1:

| Grid | A | kappa | c | rms |
|---|---|---|---|---|
| standard | 761.2 | 0.828 eV | 0.0882 | 5.20e-2 |
| doubled | 768.0 | 0.835 eV | 0.0792 | 4.65e-2 |

A stable to 0.9 %, so the two are comparable; c moves 10 %, where the MLFF moved
92 % when its truncation was removed. The offset does not converge away.

**kappa = 0.83 eV sits in the accepted 0.8-0.9 eV range**, so the branch is
physically right despite the small-q residual -- 0.08 meV^2 on a spectrum
reaching 180 meV.

The MLFF energy depends only on distances and angles, so it is rotationally
invariant by construction and gives a clean exponent of 2.10. `hiphive` would
remove the DFT offset at the cost of ~15 rattled supercells and would change no
conclusion. **State the limitation instead.**

### The ZA exponent criterion applies to the monolayer only

| System | _ai | _ml_ml |
|---|---|---|
| Monolayer | 2.10 | 2.10 |
| Bilayer | **1.33** | **1.34** |
| Trilayer | -- | 1.19 |
| Bulk | 1 | 1 |

A stack of coupled layers is 2D-flexural only below a crossover where kappa*q^4
matches the interlayer coupling; above it the branch becomes an ordinary 3D
acoustic mode. Bilayer agrees between branches to 0.8 %, so the crossover is
physical. **A result to state, not a failed check.**

---

## Stage 5 — original instructions 

`band.conf`:

```
DIM = 4 4 1
PRIMITIVE_AXES = 1 0 0  0 1 0  0 0 1
BAND = 0.0 0.0 0.0  0.5 0.0 0.0  0.3333 0.3333 0.0  0.0 0.0 0.0
BAND_LABELS = $\Gamma$ M K $\Gamma$
BAND_POINTS = 41
FC_SYMMETRY = .TRUE.
SYMMETRY = .TRUE.
SYMPREC = 1E-4
EIGENVECTORS = .TRUE.
```

`FC_SYMMETRY` applies the permutation and translational (acoustic sum rule)
conditions. VASP applies neither.

**Checks:**

| Quantity | Criterion | Was |
|---|---|---|
| Branch count | 18 (trilayer) | 18 |
| Acoustic modes at Γ | \|ω\| < 0.02 meV | −0.248 |
| Imaginary near Γ | none within 0.05 Å⁻¹ | −1.76 at q₂ |
| Shear doublet splitting | < 0.01 meV | 0.42 / 0.24 |
| ZA exponent | 1.9 – 2.1 | 1.5 |

With the stress now at 0.001 kBar the ZA exponent should recover on ASR alone.
If not, the remaining cause is unenforced rotational (Born–Huang) invariance,
needing `hiphive` and ~10–20 rattled supercells — only worth it if the
quadratic ZA claim is load-bearing, which it is not for the MLFF benchmark since
both branches are affected identically.

**Tolerance sweep:** 1e-4, 1e-3, 1e-2. Γ frequencies should move < 0.05 meV.

**Irreps:**

```bash
phonopy -c POSCAR-unitcell --irreps="0 0 0" --tolerance=1e-4 \
        --dim="4 4 1" --readfc
```

| System | Point group | Valid Γ labels |
|---|---|---|
| Bulk | `D_6h` | `A_1g A_2u B_1g B_2g E_1g E_1u E_2g E_2u` |
| Monolayer | `D_3h` | `A_1' A_2' E' A_1'' A_2'' E''` |
| Bilayer | `D_3d` | `A_1g A_2u E_g E_u` |
| Trilayer | `D_3h` | `A_1' A_2' E' A_1'' A_2'' E''` |

`E'''` does not exist in any point group. If the monolayer returns `D_6h`
labels, that calculation is not a monolayer — which would also explain the
byte-identical monolayer/trilayer ZO entry.

**NAC:** `NAC = .TRUE.` with `BORN` present.

- Bulk: correct and necessary; fixes the LO row of Table 1.
- Slabs: phonopy's NAC assumes 3D periodicity. In true 2D the splitting vanishes
  at Γ and grows linearly in q. A Γ-point splitting in the trilayer would be
  slab–slab interaction. The 0.085 meV manifold spread is the correct 2D answer
  — report slabs without NAC and say why.

**DOS** for the defect supercells and panel (h), `mesh.conf`:

```
DIM = 1 1 1
MP = 1 1 1
FORCE_CONSTANTS = READ
DOS = .TRUE.
DOS_RANGE = 0 200 0.1
SIGMA = 0.5
EIGENVECTORS = .TRUE.
```

For the in-plane/out-of-plane split, weight eigenvectors by Cartesian component
— `PDOS` gives atom projection, not direction projection.

---

## Stage 6 — Diagnostics ⬜

**6a. Error decomposition.**

```bash
python branch_analysis.py compare trilayer
```

| Comparison | Meaning |
|---|---|
| `ml_ai` vs `ai` | pure curvature error |
| `ml_ml` vs `ml_ai` | geometric contribution |
| `ml_ml` vs `ai` | what the manuscript reports |

If the rigid-mode error is mostly geometric, the fix is the relaxation, not
retraining.

This stage has no independent cross-check built in. Several metrics in this
protocol produced confident wrong answers that the physics then contradicted —
the ZBRENT recipe, the biaxial factor of two, the 3D decay distance. If a
comparison looks surprising, ask first whether it measures what it claims.

**6b. vdW consistency.** One displaced supercell, forces three ways: AI with
`IVDW = 12`; MLFF with `IVDW = 12`; MLFF without. Expect the first two to agree.

**6c. S(E) sensitivity — run before any retraining.** Recompute S(E) and S_tot
per defect with `FC_ai` and `FC_ml`, identical ΔQ.

| Quantity | Criterion |
|---|---|
| ΔS_tot | < 2 % |
| fraction of S_tot below 20 meV | quote it |
| Bhattacharyya distance between S(E) | < 0.01 |
| peak shifts in L(E) | < 1 meV |

Below 5 % of S_tot under 20 meV → report the discrepancy, skip retraining.
Above 5 % → retraining is mandatory. The Bhattacharyya distance is the metric
Turiansky *et al.* use, making the benchmark directly comparable.

---

## Stage 7 — MLFF interlayer accuracy ⬜

A **single** force field covers all four systems — a stronger result than four
separately fitted models, and worth stating. It also means the interlayer
discrepancy has a candidate cause a per-system model would not have: the model
must represent interlayer coupling across environments in which each layer has
two neighbours (bulk), one or two (trilayer), one (bilayer), or none at all
(monolayer).

**The failure is specific to sliding.** Matched-supercell bulk errors are
−29 % for E₂g shear against −5 % for B₁g breathing. Whatever the mechanism, it
affects the registry-dependent response and not interlayer binding generally.

**A complication worth taking seriously.** The DFT reference itself
underestimates the shear mode by ~3.9 % against experiment, and for a structural
reason: D3 is pairwise and isotropic, so it carries no registry dependence. The
MLFF is therefore trying to learn a sliding response that is weakly determined
in its own training data. That may limit how far 7c and 7d can go — and it is
worth stating in the manuscript regardless of the outcome.

**7a — training set composition.** Audit `ML_AB`: how many configurations from
each system, and what fraction? A monolayer configuration carries no
information about interlayer coupling. More sharply: does any configuration
sample **rigid-layer sliding**? MD at elevated temperature is dominated by
intralayer motion and will not.

**7b — cross-system error analysis.** Both a diagnostic and a publishable
result:

| | bulk | bilayer | trilayer |
|---|---|---|---|
| shear | **−29 %** | | |
| breathing | **−5 %** | | |

Already established from the geometries alone, before any phonon comparison:

| Observable | AI | MLFF | Reading |
|---|---|---|---|
| in-plane curvature ratio | 1 : 2.005 : 3.003 | 1 : 2.026 : 3.057 | 1.9 % non-additive |
| a, bilayer→trilayer step | −0.42 mÅ | +0.01 mÅ | trend flattened |
| interlayer spread, bilayer→bulk | 84.6 mÅ | 49.2 mÅ | compressed to 58 % |

- Growing with layer count → descriptor reach (7c).
- Tracking training fraction → coverage (7d).
- Uniform → systematic property; report as such.

The current Table 1 gives bulk shear at −18.8 % but bilayer shear at −4.1 %,
while bulk breathing is −2.2 % and bilayer breathing −18.7 % — inverted between
systems. Not a natural failure mode; may indicate MLFF values paired against the
wrong DFT mode. Regenerate before concluding anything from it.

**7c — descriptor cutoff.** `ML_RCUT2` (angular, default 5 Å) against
c/3 ≈ 3.3 Å — reaches the adjacent layer but not the next. Test 6 and 7 Å.
Also check `ML_RCUT1` (radial, default 8 Å) against supercell dimensions: the
finite-displacement requirement is that every dimension exceed twice the
interaction range, and 6×6 in-plane is only 15.0 Å.

**7d — targeted continuation training.** Rigid-layer **shear** displacements
±0.05–0.20 Å are the priority given the −29 % result; layer-separation changes
±1–4 % of c/3; low-T MD at 50–100 K. Sample in **bulk and trilayer**, not
monolayer — the coordinates do not exist there.

**Criterion:** breathing within 5 %, shear within 3 %, every system.

---

## Stage 8 — Rebuild Table 1 ⬜

- uniform precision, 2 decimals in cm⁻¹
- trilayer shear and breathing rows populated
- monolayer irreps in `D_3h`
- experimental 1366.2 / 1370 / 1364 cm⁻¹ moved from the LO row to the TO row
- byte-identical monolayer/trilayer ZO entries resolved
- bulk shear regenerated; do not reconcile against the old 51.99

Also add to the manuscript: the interlayer spacing sequence (3.3597 / 3.3384 /
3.2751 Å), the 1.6 mÅ surface buckling, and the separation of DFT-vs-experiment
from MLFF-vs-DFT error for the shear mode.

---

## Order

| Stage | Status | Note |
|---|---|---|
| 0 cells | ✅ | hexagonal primitive, tol ≥ 1e-4 |
| 1 relax + stress | ✅ | σ_xx < 0.0011 kBar |
| 2 supercell | ✅ | 4×4×1 slabs, 6×6×2 bulk |
| 2d bulk shear | ✅ | MLFF artefact, moved to Stage 7 |
| 3 displacements | ⬜ **in progress** | `_ai` done for the slabs; 3b geometries all ✅ |
| 4 Born | done | eps 4.802/2.749, Z*(B) +2.713 |
| 5 checks | ✅ | doublets 0, ASR exact, ZA resolved, kappa 0.82 eV |
| 5 irreps + NAC | done | Table 1 complete |
| 6b vdW | ⬜ | one force eval, do early |
| 6a decomposition | ✅ | curvature dominates; retraining, not relaxation |
| 6c S(E) | ⬜ | **before** Stage 7 |
| 7 MLFF | ⬜ | conditional on 6c |
| 8 table | done | table1.tex, all irreps verified |
| figures (a)-(d) | done | plot_dispersions.py |
| figures (e)-(h) | done | plot_defect_dos.py |

## Files

| File | Purpose |
|---|---|
| `cells_zerostress/` | AI-branch geometries, σ_xx ≈ 0 |
| `cells_ml_corrected/` | MLFF-branch geometries, σ_xx ≈ 0 |
| `kpoints/` | meshes for every cell and DIM |
| `latscan.py` | lattice scan, quartic fit |
| `stressfit.py` | σ_xx(a) → zero-stress lattice constant |
| `supercell_conv.py` | Stage 2 driver and analysis |
| `branch_analysis.py` | three-branch tree: status, gamma, drift, compare |
| `stage5_checks.py` | ASR, imaginary modes, doublets, ZA exponent |
| `za_fit.py` | omega^2 = A q^4 - c, rigidity, offset diagnosis |
| `make_born.py` | build and validate BORN, deduplicated |
| `check_loto.py` | which modes split, and along which direction |
| `irreps.py` | Gamma irreps across the tree |
| `table1.py`, `table1.tex` | Table 1 generator and output |
| `plot_dispersions.py` | panels (a)-(d) |
| `plot_defect_dos.py` | panels (e)-(h) |
| `setup_panels.py` | writes band.conf / mesh.conf for the figures |
| `phonon_benchmark.py` | MLFF vs AI error statistics over the full path |
| `stage5_checks.py` | branch count, ASR, imaginary modes, doublets, ZA |
| `za_fit.py` | omega^2 = A q^4 - c, flexural rigidity, offset diagnosis |
| `make_branch.py` | build one <system>_<branch>/<DIM>/ end to end |
| `ref_sup_additions.bib` | 15 references for the supplementary |
| `phonons_revised.tex` | revised supplementary section |
