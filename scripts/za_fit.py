#!/usr/bin/env python3
"""
Fit the flexural branch to  omega^2 = A q^4 - c  and test whether the offset c
depends on supercell size.

    python za_fit.py monolayer_ai/6x6x1 monolayer_ai/8x8x1
    python za_fit.py monolayer_*/*/                       # everything at once

Why this distinguishes the two candidates
-----------------------------------------
Harmonic 2D theory gives omega^2 = (kappa/rho) q^4 for an isolated sheet. A
violation of the rotational (Born-Huang) sum rules adds a q-independent shift,
so the branch becomes omega^2 = A q^4 - c with c > 0 fixed by the size of the
violation -- NOT by the supercell.

  c constant across supercells  -> a genuine sum-rule violation. It will not
                                   converge away; hiphive is the fix.
  c falling with supercell size -> real-space truncation. A larger supercell
                                   fixes it and hiphive is unnecessary.
  c scattered                   -> numerical noise. Tighter integration
                                   (denser FFT grid, tighter EDIFF) may help.

Reported alongside:
  A     the q^4 coefficient. It must be supercell-independent; if it is not,
        the branch itself is unconverged and c cannot be interpreted.
  q_c   (c/A)^(1/4), the wavevector below which the branch is imaginary.
        Compare it against your smallest sampled q (~0.036 1/A).
"""
import sys, os, glob
import numpy as np

THZ_TO_CM = 33.35641
CM_TO_MEV = 1.0 / 8.06554
N_FIT = 8          # q-points from Gamma to include


def load_za(path, a_inplane=2.5104):
    """Return q in 1/A and the lowest branch in meV, from band.yaml."""
    import yaml
    d = yaml.safe_load(open(os.path.join(path, "band.yaml")))
    f = np.array([[b["frequency"] for b in p["band"]] for p in d["phonon"]])
    q = np.array([p["distance"] for p in d["phonon"]]) * 2 * np.pi
    return q, f[:, 0] * THZ_TO_CM * CM_TO_MEV


def fit(path):
    q, w = load_za(path)
    sl = slice(1, 1 + N_FIT)
    qq, ww = q[sl], w[sl]
    y = np.sign(ww) * ww ** 2          # signed omega^2, so imaginary points count
    x = qq ** 4
    A, negc = np.polyfit(x, y, 1)
    c = -negc
    resid = y - (A * x + negc)
    rms = np.sqrt((resid ** 2).mean())

    qc = (c / A) ** 0.25 if (c > 0 and A > 0) else 0.0

    # flexural rigidity from A = (kappa/rho) q^4, with omega in meV, q in 1/A
    a0 = 2.5104e-10
    rho = (10.811 + 14.007) * 1.66053907e-27 / (a0 ** 2 * np.sqrt(3) / 2)
    hbar, mev = 1.0545718e-34, 1.602176634e-22
    kappa = A * rho * mev ** 2 / (hbar ** 2 * 1e40) / 1.602176634e-19

    return dict(A=A, c=c, rms=rms, qc=qc, kappa=kappa,
                npos=int((ww > 0).sum()), first=w[1], q1=q[1])


if __name__ == "__main__":
    paths = []
    for p in sys.argv[1:]:
        paths.extend(sorted(glob.glob(p)))
    paths = [p for p in paths
             if os.path.isdir(p) and os.path.exists(os.path.join(p, "band.yaml"))]
    if not paths:
        print(__doc__)
        sys.exit(1)

    print(f"fitting omega^2 = A q^4 - c over the first {N_FIT} q-points\n")
    # Group by branch. Pooling branches is meaningless: each force engine has
    # its own convergence behaviour, and the whole point is to compare them.
    groups = {}
    for p in paths:
        try:
            r = fit(p)
        except Exception as e:
            print(f"{p:<34}  {type(e).__name__}: {e}")
            continue
        key = os.path.dirname(p.rstrip("/")) or p
        groups.setdefault(key, []).append((os.path.basename(p.rstrip("/")), r))

    print(f"{'branch / DIM':<34}{'A':>10}{'kappa (eV)':>12}{'c (meV^2)':>12}"
          f"{'q_c (1/A)':>11}{'rms resid':>12}")
    for key, rows in groups.items():
        print(f"\n{key}")
        for tag, r in rows:
            print(f"  {tag:<32}{r['A']:>10.1f}{r['kappa']:>12.3f}"
                  f"{r['c']:>12.5f}{r['qc']:>11.4f}{r['rms']:>12.2e}")
        if len(rows) < 2:
            continue
        cs = np.array([r['c'] for _, r in rows])
        As = np.array([r['A'] for _, r in rows])
        drop = cs[-1] / cs[0] if cs[0] != 0 else np.inf
        spread = (cs.max() - cs.min()) / max(abs(cs.mean()), 1e-12)
        print(f"    c: {cs[0]:.5f} -> {cs[-1]:.5f}   ratio {drop:.2f}")
        print(f"    A spread {100*(As.max()-As.min())/abs(As.mean()):.0f} %")
        if drop < 0.5:
            print("    -> c collapses with supercell size: real-space")
            print("       truncation, converged. No sum-rule problem here.")
        elif spread < 0.30:
            print("    -> c holds steady: a fixed offset. Either a rotational")
            print("       (Born-Huang) violation or numerical noise mimicking")
            print("       one. The rms residual tells you which: if it is a")
            print("       large fraction of c, the fit is noise-dominated and")
            print("       a denser FFT grid should be tried before hiphive.")
        else:
            print("    -> c neither collapses nor holds steady; the branch is")
            print("       probably not converged. Check A first.")

    print("\n  A must be supercell-independent whatever c does -- it is the")
    print("  physical bending rigidity. Literature kappa for monolayer h-BN")
    print("  is ~0.8-0.9 eV.")
