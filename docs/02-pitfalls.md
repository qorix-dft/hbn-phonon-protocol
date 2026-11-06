# Twenty-three ways a phonon calculation goes quietly wrong

Every entry here cost real time. They are grouped by the stage at which they
bite, and each gives the symptom you actually see, the cause, how to detect it,
and the fix.

The common feature is that **none of them crash**. The calculation finishes, the
bands look like bands, and the error sits in a number you were not looking at.

---

## Cells and symmetry

### 1. Orthorhombic cells lose the 3-fold

**Symptom.** Modes that must be exactly degenerate are split by 0.1–0.5 meV.
In h-BN the interlayer shear doublets came out 0.42 and 0.24 meV apart.

**Cause.** A 3-fold rotation is not an integer matrix in a C-centred
orthorhombic basis. spglib therefore finds 8 symmetry operations instead of 24
and reports point group `mm2`. Everything downstream — the symmetrisation of the
force constants, the irrep assignment, the degeneracies — is then done in the
wrong group.

It is tempting to blame the electronic convergence; we did, and were wrong. The
splitting is geometric.

**Detect.**

```bash
phonopy --symmetry --tolerance=1e-4 -c POSCAR | head -3
```

Count the operations. If you expected 24 and got 8, this is why.

**Fix.** Use the hexagonal primitive cell for the pristine band structure.
Orthorhombic supercells are often the right choice for defect embedding — just
do not use one for a dispersion.

### 2. `ISIF = 3` aborts on a slab with vacuum

**Symptom.**

```
ZBRENT: fatal error in bracketing
  please rerun with smaller EDIFF, or copy CONTCAR to POSCAR and continue
```

Following that advice does not help.

**Cause.** σ_zz for a cell that is mostly vacuum is numerical noise referenced to
a volume dominated by empty space. The conjugate-gradient line minimiser cannot
bracket a minimum along a direction where the gradient is noise, so it gives up.
Sometimes it does something worse first and collapses the vacuum gap.

**Fix.** Never relax the cell of a vacuum-containing slab. Scan the in-plane
lattice constant instead (`scripts/latscan.py`), relaxing ions only at each
point with `ISIF = 2`. For the ionic relaxation use **`IBRION = 1`** (RMM-DIIS),
not `IBRION = 2` — ZBRENT *is* the CG line minimiser, and RMM-DIIS avoids it by
construction.

Also: ignore σ_zz for slabs entirely. It is a basis-set artefact along c.

### 3. A parabola is not good enough for E(a)

**Symptom.** The E(a) minimum and the σ_xx = 0 point disagree by ~1.5 mÅ, more
than you can explain.

**Cause.** Two separate effects, both real.

Over a ±0.04 Å window, E(a) carries enough cubic anharmonicity that a parabolic
fit biases the minimum by ~0.63 mÅ. Fit a **quartic** with a centred abscissa.

The remaining ~0.85 mÅ is **Pulay stress** and is not a fitting error: the scan
minimises energy with the plane-wave basis varying between points, whereas
VASP's stress tensor is the derivative at fixed basis. They genuinely do not
share a zero. Correct to σ_xx = 0 with three statics per system
(`scripts/stressfit.py`), and take that as the production geometry.

A machine-learned force field has no plane-wave basis and therefore no Pulay
stress; for an MLFF branch the two criteria coincide. Running both anyway is a
cheap and informative check — disagreement means the MLFF stress tensor is not
what it should be.

### 4. Incommensurate q-points

**Symptom.** Two calculations that should agree differ by tens of cm⁻¹ at M or
K, with no physical explanation.

**Cause.** A finite-displacement supercell of dimension `DIM` represents exactly
those q-points that are commensurate with it; everything else is Fourier
interpolation. Comparing an exact frequency in one supercell against an
interpolated one in another compares a number against an approximation to a
different number.

In our case this produced apparent 53 cm⁻¹ errors where the true difference was
0.1 cm⁻¹.

**Detect.** `scripts/supercell_conv.py` has a `commensurate(dim, point)` helper;
every comparison it makes is gated on it. M = (½,0,0) needs an even first index;
K = (⅓,⅓,0) needs DIM divisible by 3.

**Fix.** Either compare only at commensurate q, or compare full interpolated
paths against each other (`scripts/phonon_benchmark.py`) — never one of each.

### 5. Minimum image by rounding fails in a hexagonal lattice

**Symptom.** Force-constant decay looks slower than it is; the required
supercell looks bigger than it is.

**Cause.** The usual trick for the minimum-image distance — round the fractional
coordinate difference to the nearest integer and convert — is valid only for an
orthogonal cell. In a hexagonal lattice, with a 120° angle between **a** and
**b**, it overestimates in-plane separations by up to 33%.

**Fix.** Search explicitly over the 27 (or 9, in-plane) neighbouring lattice
translations and take the minimum. `supercell_conv.py` does this.

A related trap: in a slab, the 3D pair distance is √(ρ² + t²) where t is the
layer separation, so a "decay with distance" plot is confounded by slab
thickness. Decompose, and plot against **in-plane** separation ρ.

### 6. phonopy's `distance` field has no 2π

**Symptom.** A fit restricted to "small q" silently includes q-points far
outside the intended window; a ZA exponent comes out wrong because the fit
range was six times too wide.

**Cause.** The `distance` entry in `band.yaml` is a reciprocal-space path length
computed **without** the 2π factor. A cut at `distance < 0.05` is not
q < 0.05 Å⁻¹; it is q < 0.31 Å⁻¹.

**Fix.** Multiply by 2π before interpreting `distance` as a physical wavevector.
`za_fit.py` and `stage5_checks.py` do.

### 7. Mismatched `ENCUT` between stages

**Symptom.** A finished, well-converged relaxation reports a residual stress of
several kBar that will not go away.

**Cause.** Relaxing at one cutoff and checking the stress at another builds in a
Pulay stress offset. Relaxing at 850 eV and checking at 650 eV gave us a
built-in −3.17 kBar.

**Fix.** One cutoff for everything, and if you are benchmarking a force field,
make it **the cutoff the force field was trained at**. Identical `KPOINTS` files
across branches for the same reason — the MLFF does no SCF, so sharing them
costs nothing and removes a difference.

---

## Born charges and LO–TO splitting

### 8. `LCALCEPS`, not `LEPSILON`, when `IVDW` is on

**Symptom.** VASP refuses to run, or runs and silently drops the vdW
correction.

**Cause.** Density-functional perturbation theory (`LEPSILON`, `IBRION = 7`/`8`)
is incompatible with a vdW correction. With `IVDW = 12` active, DFPT is closed
to you.

**Fix.** Finite differences: `LCALCEPS = .TRUE.` plus `LPEAD = .TRUE.` with
`IBRION = 6`, `NFREE = 2`, `POTIM = 0.015`. See
[`examples/hBN/incar/INCAR.born`](../examples/hBN/incar/INCAR.born).

### 9. `BORN` wants only the symmetry-independent atoms

**Symptom.**

```
Too many atoms in the BORN file
```

**Cause.** `phonopy-vasp-born` writes a Z* block for every atom in the cell.
phonopy expects one per symmetry-independent site.

**Fix.** Trim to one representative per site — for bulk h-BN, one B and one N,
not four atoms. `scripts/make_born.py` deduplicates automatically and also
checks the acoustic sum rule on the charges (Σ Z* should be 0 to <0.01) and
validates ε against expectation.

### 10. `Q_DIRECTION` is a conf-file tag, not a command-line flag

**Symptom.** `--nac-q-direction=1 0 0` is rejected as an unrecognised option.

**Cause.** It does not exist as a CLI flag.

**Fix.** Put `Q_DIRECTION = 1 0 0` in the conf file. And while you are there:
**the phonopy Python API path for NAC proved unreliable** — it reported "NAC on"
while returning unsplit frequencies. A conf file with `NAC = .TRUE.` and a short
band segment off Γ worked every time. See
[`examples/hBN/conf/loto_inplane.conf`](../examples/hBN/conf/loto_inplane.conf)
and `scripts/check_loto.py`.

Two further things about LO–TO that are easy to get backwards:

- **Only modes with a macroscopic dipole along q split.** Run Γ three ways — no
  NAC, q in-plane, q along c — and confirm each mode moves only in the direction
  where its dipole lies. In bulk h-BN E₂g never splits in either direction,
  because it is the out-of-phase Davydov partner and the two layers' dipoles
  cancel.
- **Do not apply NAC to a 2D slab.** phonopy's NAC assumes 3D periodicity. In a
  true 2D system the splitting vanishes at Γ and grows linearly in q; a Γ-point
  splitting in a slab calculation is slab–slab interaction, i.e. an artefact of
  your periodic images.

### 11. The acoustic sum rule is not rotational invariance

**Symptom.** A perfectly relaxed, zero-stress monolayer still gives a flexural
branch with ω ∝ q^1.5 instead of q², or a small negative offset in a ω² vs q⁴
fit that will not go to zero.

**Cause.** `FC_SYMMETRY = .TRUE.` enforces the permutation and **translational**
(acoustic) sum rules. It does not enforce **rotational** (Born–Huang)
invariance, which is the condition the quadratic flexural dispersion actually
depends on.

Three contributions have to be separated before blaming the force engine:

| Contribution | Diagnostic | Our result |
|---|---|---|
| Residual in-plane strain | reduce σ_xx, refit | fixed: 1.7 kBar → 0.001 |
| Supercell truncation | enlarge DIM, refit | fixed: 4×4×1 → 8×8×1 |
| Rotational invariance | double the FFT grid, refit | irreducible, c ≈ 0.08 meV² |

The FFT-grid test is the decisive one: if the offset is numerical it falls with
grid density. Ours fell by 10% while the MLFF's fell by 92% — so the MLFF's
offset was numerical and the DFT's was not.

**Fix, or rather decision.** Enforcing rotational invariance needs a
constrained fit (`hiphive` with ~10–20 rattled supercells). Whether that is
worth it depends on whether the quadratic claim is load-bearing. For a force
field benchmark it usually is not — both branches are affected the same way — so
state the limitation and move on. Our κ = 0.83 eV agrees with the literature
regardless.

### 12. Force constants do not decay to zero in a polar material

**Symptom.** The supercell convergence test never passes at a 10⁻³ threshold,
however large the cell.

**Cause.** A polar material has a long-range dipole–dipole contribution to the
force constants with an r^−2.5 tail. That is physics, not truncation error.

**Fix.** Judge convergence on the frequencies you care about, at commensurate
q-points, not on an absolute force-constant threshold. If you do want a
threshold, set it from the frequency sensitivity, not from a round number.

---

## Post-processing and analysis

### 13. Pairing degenerate modes by index

**Symptom.** Enormous "splittings" — we briefly reported 76–89 meV.

**Cause.** Sorting modes by frequency and pairing adjacent ones assumes the
degenerate partners are adjacent. They are not. Our pairing matched a breathing
mode with a ZO mode.

**Fix.** Select the modes to pair by **irrep or eigenvector character**, not by
index. For an N-layer stack the rigid-layer manifold has N−1 shear doublets and
N−1 breathing singlets; pair only within the shear set. `stage5_checks.py` sets
`n_doublet = 1 if system == "bulk" else nl - 1`.

### 14. Comparing two branches at different supercell sizes

**Symptom.** An error metric that shifts by ~10% depending on which runs
happened to finish.

**Cause.** A script that independently picks the largest available `DIM` per
branch will happily compare a 4×4×2 reference against a 5×5×4 test. The
difference then includes a supercell effect you are attributing to the force
field.

**Fix.** Pin the comparison branch to the reference branch's `DIM`.
`table1.py` does this, and prints which `DIM` it used.

### 15. The stress slope is a·k/2V, not a·k/V

**Symptom.** A perfectly good σ_xx(a) dataset appears to disagree with the
elastic constant from E(a) by a factor of two.

**Cause.** Scaling the in-plane lattice constant `a` strains the cell
**biaxially**. The extra factor of 2 follows.

Our predicted slopes of 935.5 / 1543.0 / 1961.5 kBar/Å matched the measured
936.1 / 1540.4 / 1961.3 to 0.2% once the factor was right. If your data and your
model disagree by exactly 2, suspect the model.

### 16. `EDIFF = 1E-6` is too loose for force constants

**Cause.** With a 0.01 Å displacement, a force error δF becomes a
force-constant error of δF/0.01 — a 100× amplification. It lands hardest on the
weak interlayer terms, which are exactly the ones a layered-material benchmark
is about.

**Fix.** `EDIFF = 1E-8` for all displacement runs and all geometry work.

### 17. VASP's in-INCAR separator is `;`, not whitespace

**Symptom.** Tags are silently ignored, and the run uses defaults you did not
ask for.

**Fix.** One tag per line, or separate with `;`. Every INCAR in
[`examples/hBN/incar/`](../examples/hBN/incar/) is one tag per line for this
reason.

### 18. `LPHON_READ_FORCE_CONSTANTS` makes every other tag inert

**Symptom.** Changing INCAR settings has no effect on the output at all.

**Cause.** When VASP reads force constants from file, it is not doing a
calculation; the rest of the INCAR is decoration.

**Fix.** Know which mode you are in. More generally, this is a reason to prefer
phonopy over VASP's internal `LPHON_DISPERSION` machinery: with phonopy, the
symmetrisation, the sum rules and the interpolation are all explicit and
inspectable, and `FC_SYMMETRY` gives you sum rules VASP does not apply at all.

### 19. An MLFF run writes no `TOTEN` line

**Symptom.** An energy-parsing script returns nothing for the force-field runs
while working fine for the ab initio ones.

**Cause.** With `ML_MODE = run` there is no SCF cycle, so `OUTCAR` never gets
the `free  energy   TOTEN` line.

**Fix.** Parse `e_fr_energy` from `vasprun.xml`, which both branches write
identically. `latscan.py` reads `vasprun.xml` first and falls back to `OUTCAR`.

### 20. Mismatched `SIGMA` between DOS panels you intend to compare

**Symptom.** Two phonon DOS curves differ in a way that looks like physics and
is smearing.

**Fix.** Fix `SIGMA` once in the conf file and reuse it for every panel. If the
panels also need to match in cell and q-sampling — e.g. a pristine reference
against defect supercells — match those too, and say so in the caption.

### 21. ε_zz: electronic versus static

**Symptom.** Your out-of-plane dielectric constant (2.75 for h-BN) looks far
too low against the textbook 4.10, and you start hunting for a bug.

**Cause.** 4.10 is the **static** constant, which includes the ionic
contribution. `LCALCEPS` gives you the **electronic** one, which is smaller. The
modern static figure is 3.4 ± 0.2.

**Fix.** Compare like with like. And when quoting an out-of-plane LO–TO
splitting computed from the electronic ε_zz, caveat it — ours came out +40%
against experiment for exactly this reason, while the in-plane splitting was
−3.2%.

### 22. `E'''` does not exist

**Symptom.** An irrep label that does not appear in any character table.

**Fix.** Check the returned labels against the point group you expect. Valid Γ
label sets:

| System | Point group | Labels |
|---|---|---|
| Bulk | `D_6h` | `A_1g A_2u B_1g B_2g E_1g E_1u E_2g E_2u` |
| Monolayer | `D_3h` | `A_1' A_2' E' A_1'' A_2'' E''` |
| Bilayer | `D_3d` | `A_1g A_2u E_g E_u` |
| Trilayer | `D_3h` | `A_1' A_2' E' A_1'' A_2'' E''` |

A bilayer is **D_3d**, not D_3h — the inversion centre between the layers
changes the group, and the Γ labels with it. If a monolayer returns `D_6h`
labels, it is not a monolayer.

`scripts/irreps.py` checks the returned set against the expected one for every
system in the tree.

### 23. Retracting a conclusion drawn from one branch

Not a software trap. We twice declared a question settled on the basis of the
MLFF data alone and had the ab initio data contradict it a day later. The
three-branch design exists precisely so that neither engine's quirks are
mistaken for physics — use all three before concluding anything, and prefer a
test whose outcome discriminates (the FFT-grid test in
[Pitfall 11](#11-the-acoustic-sum-rule-is-not-rotational-invariance)) over one
that merely looks consistent.
