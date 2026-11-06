# Benchmark design: three force-constant sets, not two

The obvious way to benchmark a machine-learned force field against phonons is to
relax with the force field, compute its phonons, and compare against a DFT
dispersion. That number is what a user of the force field experiences, so it is
worth reporting — but on its own it cannot tell you *what is wrong*, and it can
be small for the wrong reason.

## The decomposition

A phonon frequency depends on two things the force field can get wrong: the
curvature of the potential energy surface, and the location of the minimum it is
evaluated at. Three force-constant sets separate them.

| Set | Geometry | Forces | Directory |
|---|---|---|---|
| `FC_ai` | ab initio | ab initio | `<system>_ai/` |
| `FC_ml@ai` | **ab initio** | MLFF | `<system>_ml_ai/` |
| `FC_ml` | MLFF | MLFF | `<system>_ml_ml/` |

- **`FC_ml@ai` vs `FC_ai`** — same geometry, different force engine. Pure
  curvature error.
- **`FC_ml` vs `FC_ml@ai`** — same force engine, different geometry. The
  geometric contribution.
- **`FC_ml` vs `FC_ai`** — the end-to-end error. Report this, but do not
  interpret it alone.

For the middle branch to mean anything, `_ai` and `_ml_ai` must share the
**same displacement pattern**, not merely the same geometry. Generate the
displacements once with `phonopy -d` and copy the `POSCAR-{number}` files into
both trees. Then the only variable is which code computes the forces.

## Why two numbers are not enough: errors cancel

In our bulk h-BN case the curvature error was 9.1 cm⁻¹ and the end-to-end error
13.0 cm⁻¹ — but for several modes the geometric and curvature contributions
carried **opposite signs**, so the end-to-end error was smaller than the
curvature error alone.

A force field can therefore look better than it is. The force field's
equilibrium geometry is slightly wrong in a direction that partly compensates
for its slightly wrong curvature, and a two-set benchmark reports the
compensated number as if it were the model's accuracy. That cancellation is not
stable: it depends on the system, and in our data it held for the bulk while the
shear modes flipped sign between bilayer (+16.6%), trilayer (+26.2%) and bulk
(−18.6%).

This is the methodological point the comparable literature does not generally
make, and it is the main reason this repository exists in the form it does.

## Diagnostics that discriminate

Once the three sets are in hand, a few derived quantities localise the error far
better than an RMSE:

**Sign-resolved per-mode errors.** A uniform sign across a manifold suggests a
systematic curvature problem; sign changes with layer count suggest the
interlayer term specifically.

**Elastic additivity.** In-plane stiffness should scale exactly with layer
count. Ours did in DFT (1 : 2.005 : 3.003) and was violated by 1.9% in the MLFF
(1 : 2.026 : 3.057) — a clean, geometry-independent signal that the force field
is not treating the layers as independent sheets.

**Interlayer spacing spread.** DFT gave an 84.6 mÅ spread across bilayer →
trilayer → bulk; the MLFF compressed it to 49.2 mÅ (58%). The force field is
under-resolving how binding changes with coordination.

**Stiffening with layer count.** DFT stiffens the shear mode by 63% from bilayer
to bulk; the MLFF by only 7%. This is the same deficiency seen from a third
angle, and it is the one that explains the shear errors.

**Rigid-layer manifold.** For an N-layer stack there are N−1 shear doublets and
N−1 breathing singlets. These are the softest modes and the most sensitive to
the interlayer description, so they are where a layered-material force field
fails first — and where you should look before the optical modes.

## A note on what the reference can and cannot settle

The ab initio branch is the reference for the *force field*, not for nature. In
our case DFT-D3 itself underestimates the bulk shear mode by 3.9% against
experiment, because pairwise dispersion carries no registry dependence and so
gives too soft a sliding barrier. A force field trained on DFT-D3 inherits that
error and cannot be blamed for it.

Keep the two comparisons separate and label them:

| Comparison | Our result | Attributable to |
|---|---|---|
| DFT-D3 vs experiment | −3.9% | the functional |
| MLFF vs its own DFT | −29.2% | the force field |

Collapsing these into "the force field is off by 33%" is wrong, and so is
quoting the −29.2% as if the reference were experiment.
