#!/usr/bin/env python3
"""
Build 6-atom hexagonal primitive trilayer h-BN cells with a controlled
vacuum gap, for the Stage 00c vacuum-convergence test.

The 12-atom orthorhombic cell cannot express the 3-fold rotation as an
integer matrix, so spglib/phonopy/VASP find only mm2 (8 operations) and
the E' / E'' doublet degeneracies are never enforced. The hexagonal
primitive cell restores all 12 operations of D_3h.
"""
import numpy as np, spglib

A       = 2.5065290928      # in-plane lattice constant (from bulk relaxation)
D_INTER = 3.28229           # interlayer spacing
GAPS    = [12, 16, 20, 24, 28]   # vacuum gap in Angstrom
OUTDIR  = "/mnt/user-data/outputs/cells"

# hexagonal primitive in-plane vectors
a1 = np.array([A, 0.0, 0.0])
a2 = np.array([-A/2, A*np.sqrt(3)/2, 0.0])

# AA'A stacking, fractional in-plane sites
B_SITE = np.array([1/3, 2/3])
N_SITE = np.array([2/3, 1/3])

for gap in GAPS:
    slab = 2*D_INTER
    c    = slab + gap
    a3   = np.array([0.0, 0.0, c])
    cell = np.array([a1, a2, a3])
    z0   = (c - slab)/2.0            # centre the slab in the cell

    frac, nums, labels = [], [], []
    for L, zc in enumerate([z0, z0 + D_INTER, z0 + 2*D_INTER]):
        # AA' : layers 0 and 2 identical, layer 1 has B/N swapped
        b_xy, n_xy = (B_SITE, N_SITE) if L % 2 == 0 else (N_SITE, B_SITE)
        frac.append([b_xy[0], b_xy[1], zc/c]); nums.append(5);  labels.append('B')
        frac.append([n_xy[0], n_xy[1], zc/c]); nums.append(7);  labels.append('N')

    frac = np.array(frac)
    order = np.argsort(nums, kind='stable')
    frac, nums = frac[order], [nums[i] for i in order]

    ds = spglib.get_symmetry_dataset((cell, frac, nums), symprec=1e-4)
    has3 = any((not np.allclose(R, np.eye(3)))
               and np.allclose(np.array(R) @ R @ R, np.eye(3))
               for R in ds.rotations)

    fn = f"{OUTDIR}/POSCAR_trilayer_vac{gap}"
    with open(fn, "w") as f:
        f.write(f"trilayer hBN AA'A | 6-atom hexagonal primitive | vacuum {gap} A\n1.0\n")
        for v in cell:
            f.write(f"  {v[0]:22.16f}{v[1]:22.16f}{v[2]:22.16f}\n")
        f.write("  B  N\n   3   3\nDirect\n")
        for p in frac:
            f.write(f"  {p[0]%1:22.16f}{p[1]%1:22.16f}{p[2]%1:22.16f}\n")

    print(f"vac={gap:3d} A  c={c:7.3f}  {ds.international} (#{ds.number})  "
          f"ops={len(ds.rotations):3d}  3-fold={has3}  -> {fn.split('/')[-1]}")
