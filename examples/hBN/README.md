# Worked example: hexagonal boron nitride

Four systems — bulk, monolayer, bilayer, trilayer — carried through the full
protocol with a machine-learned force field benchmarked against its own ab
initio reference.

**Setup.** VASP 6, PBE + DFT-D3(BJ) (`IVDW = 12`), ENCUT 650 eV (the cutoff the
force field was trained at), finite displacements of 0.01 Å via phonopy. Slabs
carry 15.157 Å of vacuum — 4.6× the interlayer spacing.

Every number below is measured, not quoted from the literature, unless marked
as experimental.

---

## Contents

| Path | What |
|---|---|
| `cells/` | Zero-stress ab initio geometries, hexagonal primitive |
| `cells_mlff/` | The force field's own equilibrium geometries |
| `incar/` | INCAR templates, one tag per line |
| `kpoints/` | KPOINTS scaled with supercell size (~0.24 Å⁻¹) |
| `conf/` | phonopy conf files: band, DOS, NAC, LO–TO tests |
| `reference/table1.tex` | The Γ-frequency table as published |
| `reference/PROTOCOL_as_run.md` | The working log, with dead ends intact |
| `../../figures/` | The benchmark figure |

`reference/PROTOCOL_as_run.md` is the unedited working document. It is messier
than `docs/01-workflow.md` and more useful if you hit something the clean
protocol does not cover — it records what was tried and what the data said at
each point, including the retractions.

---

## Stage 0 — Cells

```bash
for s in bulk monolayer bilayer trilayer; do
    phonopy --symmetry --tolerance=1e-4 -c cells/POSCAR_$s.vasp | head -3
done
```

| System | Atoms | Modes | Space group | Ops | Point group |
|---|---|---|---|---|---|
| Bulk | 4 | 12 | `P6_3/mmc` (#194) | 24 | D_6h |
| Monolayer | 2 | 6 | `P-6m2` (#187) | 12 | D_3h |
| Bilayer | 4 | 12 | `P-3m1` (#164) | 12 | D_3d |
| Trilayer | 6 | 18 | `P-6m2` (#187) | 12 | D_3h |

The bilayer is D_3d, not D_3h — Γ labels are `A_1g A_2u E_g E_u`.

Vacuum convergence was checked against 12, 16, 20, 24 and 28 Å; 15.157 Å is
converged.

---

## Stage 1 — Geometry

```bash
python ../../scripts/latscan.py gen monolayer    # then bilayer, trilayer
python ../../scripts/latscan.py fit monolayer
python ../../scripts/stressfit.py
```

| System | a (Å) | σ_xx (kBar) | Interlayer (Å) |
|---|---|---|---|
| Monolayer | 2.508960 | +0.00024 | — |
| Bilayer | 2.507926 | −0.00038 | 3.3597 |
| Trilayer | 2.507510 | −0.00107 | 3.3384 |
| Bulk | 2.506566 | +0.00782 | 3.2751 (c = 6.55028) |

Interlayer spacing is monotonic — bilayer > trilayer > bulk. Fewer neighbours,
weaker binding, larger spacing. The outer layers buckle, splitting B and N by
1.6 mÅ; the trilayer middle layer stays planar by symmetry.

Against experiment: bulk a is +0.10% (expt. 2.504 Å), bulk c/2 is **−1.66%**
(expt. 3.3305 Å). The c contraction matters for interpreting the shear mode.

**Consistency checks, all passing:**

| Check | Result |
|---|---|
| In-plane curvature, mono : bi : tri | 38.48 : 77.17 : 115.54 eV/Å² = 1 : 2.005 : 3.003 |
| Interlayer binding per interface | −93.4 meV |
| Trilayer/bilayer binding ratio | 2.088 (expected 2.0; 4% excess from the D3 three-body term) |
| Stress slope vs predicted a·k/2V | 936.1 / 1540.4 / 1961.3 vs 935.5 / 1543.0 / 1961.5 kBar/Å |

Stress targets were 0.014 / 0.024 / 0.030 kBar; achieved 0.001 or better.

---

## Stage 2 — Supercell convergence

Production sizes:

| System | DIM | Displacements |
|---|---|---|
| Bulk | 4×4×2 | 2 |
| Monolayer | 8×8×1 | 2 |
| Bilayer | 6×6×1 | 4 |
| Trilayer | 6×6×1 | 6 |

The bulk shear mode showed an apparent c-axis drift with supercell size that
turned out to be a **force-field artefact**: on matched supercells the ab initio
result is flat at +0.137 while the MLFF drifts to −0.835.

---

## Stage 3 — Branches

```
bulk_ai       monolayer_ai       bilayer_ai       trilayer_ai
bulk_ml_ai    monolayer_ml_ai    bilayer_ml_ai    trilayer_ml_ai
bulk_ml_ml    monolayer_ml_ml    bilayer_ml_ml    trilayer_ml_ml
```

The `_ai` and `_ml_ai` trees share displacement patterns, so the only variable
between them is the force engine. `cells_mlff/` holds the geometries for the
`_ml_ml` branch, obtained by repeating Stage 1 with the force field. As
expected, the MLFF's E(a) minimum and σ_xx = 0 point coincide — no plane-wave
basis, no Pulay stress.

---

## Stage 4 — Born charges (bulk only)

```bash
python ../../scripts/make_born.py
python ../../scripts/check_loto.py
```

| Quantity | Value | Expected |
|---|---|---|
| ε_∞ in-plane | 4.802 | 4.8 – 5.1 |
| ε_∞ out-of-plane | 2.749 | 2.6 – 2.9 (electronic) |
| Z*(B) in-plane | +2.713 | ~+2.7 |
| Charge neutrality | 0.00000 | < 0.01 |

The out-of-plane value looks low against Geick's 4.10 — but that is the
**static** constant. The electronic one is smaller and the modern static figure
is 3.4 ± 0.2. The 4.802 / 2.749 anisotropy is h-BN's hyperbolic behaviour.

Direction-resolved LO–TO:

| Mode | q in-plane | q along c |
|---|---|---|
| E_1u (in-plane optical) | 1356.03 → **1591.36** cm⁻¹ | unmoved |
| A_2u (ZO) | unmoved at 741.42 | → **805.01** |
| E_2g, B_1g | unmoved | unmoved |

In-plane splitting **235.33 cm⁻¹** against 243 measured (−3.2%). Out-of-plane
63.6 against Geick's 45 (+40%) — consistent with ε_zz being the electronic
value, so do not quote the A_2u LO without that caveat. Table 1 reports the TO,
which is unaffected.

E_2g is the out-of-phase Davydov partner: the two layers' dipoles cancel, so it
has no macroscopic polarisation and never splits. On the Γ–M–K–Γ path the A_2u
dipole is always perpendicular to **q**, so the dispersion panel shows its TO
throughout.

---

## Stage 5 — Validation

| Check | Result |
|---|---|
| Branch count | 12 / 6 / 12 / 18 as required |
| Acoustic modes at Γ | within 1e-6 THz after `FC_SYMMETRY` |
| Degenerate doublet splitting | **exactly 0** (was 0.42 / 0.24 meV in the orthorhombic cell) |
| ZA exponent, monolayer | 2.10 |
| Γ irreps | match D_6h / D_3h / D_3d / D_3h |

Flexural rigidity from ω² = A q⁴ − c: **κ = 0.83 eV** (DFT), 0.955 eV (MLFF,
+17%). The residual offset c = 0.08 meV² is the rotational-invariance limitation
discussed in [Pitfall 11](../../docs/02-pitfalls.md#11-the-acoustic-sum-rule-is-not-rotational-invariance);
doubling the FFT grid moved the DFT offset by 10% and the MLFF offset by 92%,
confirming the DFT residual is physical to the method and the MLFF one was
numerical.

**2D → 3D flexural crossover**, exponents by system:

| Monolayer | Bilayer | Trilayer | Bulk |
|---|---|---|---|
| 2.10 | 1.33 | 1.24 | 1 |

Both branches agree on this, which is the point — it is a real dimensional
crossover, not a force-field artefact.

---

## Stage 6 — Benchmark results

Full-path trilayer RMSE: **0.475 meV** (0.690 meV in the pre-protocol setup).

Where the force field is accurate and where it is not:

| Quantity | DFT | MLFF | Error |
|---|---|---|---|
| Optical modes, all systems | — | — | within ~3% |
| Shear, bilayer | — | — | **+16.6%** |
| Shear, trilayer | — | — | **+26.2% / +21.8%** |
| Shear, bulk | — | — | **−18.6%** |
| Shear stiffening, bilayer → bulk | +63% | +7% | — |
| In-plane elastic additivity | 1 : 2.005 : 3.003 | 1 : 2.026 : 3.057 | 1.9% |
| Interlayer spacing spread | 84.6 mÅ | 49.2 mÅ | 58% of DFT |
| Flexural rigidity κ | 0.83 eV | 0.955 eV | +17% |

The diagnosis is consistent across four independent measures: the force field
under-resolves how interlayer binding changes with coordination. Optical modes,
which are intra-layer, are fine.

**Error cancellation.** In the bulk the curvature error is 9.1 cm⁻¹ while the
end-to-end error is 13.0 cm⁻¹, and for several modes the geometric and curvature
contributions carry opposite signs — so a two-set benchmark would have reported
a smaller number than the force field deserves. See
[`docs/04-benchmark-design.md`](../../docs/04-benchmark-design.md).

**Separating the functional from the force field:**

| Comparison | Error | Attributable to |
|---|---|---|
| DFT-D3 vs experiment (bulk shear) | −3.9% | the functional — pairwise dispersion has no registry dependence, so the sliding barrier is too soft |
| MLFF vs its own DFT (bulk shear) | −29.2% | the force field |

---

## Reproducing the figure

```bash
python ../../scripts/setup_panels.py
python ../../scripts/plot_dispersions.py     # panels (a)-(d)
python ../../scripts/plot_defect_dos.py      # panels (e)-(h)
```

`setup_panels.py` writes `band.conf` (BAND_POINTS = 101, NAC for bulk only) and
`mesh.conf` (MP 40×40×1 slabs, 30×30×10 bulk) into all eight panel directories.
Both plotting scripts use the `_ml_ml` branch as the comparison; keep that
consistent with whatever `table1.py --cmp` you report.

---

## Known limitations of this example

Stated plainly, because they bound what the numbers support:

- **Rotational invariance is not enforced.** The residual flexural offset
  survives convergence. Both branches are affected identically, so the benchmark
  is unaffected, but the absolute quadratic claim is approximate.
- **Slabs are reported without NAC**, which is correct for a true 2D system but
  means slab LO–TO cannot be read off these panels.
- **The defect DOS panels shade fixed energy windows**, not an eigenvector
  projection. The in-plane/out-of-plane character labels are indicative.
- **The out-of-plane LO–TO splitting uses the electronic ε_zz** and is
  consequently ~40% high against experiment.
