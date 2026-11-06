#!/usr/bin/env bash
# ---------------------------------------------------------------------
# phonon_triage.sh
#
# Finds and classifies existing phonon calculation directories so you can
# work out what you already have before rerunning anything.
#
# Usage:   bash phonon_triage.sh /path/to/search/root
# Example: bash phonon_triage.sh ~/wrk/UQ/PhD/first_calc/MLFF
# ---------------------------------------------------------------------

ROOT="${1:-.}"

echo "==================================================================="
echo " Searching under: $ROOT"
echo "==================================================================="

# ---------------------------------------------------------------------
# 1. Every vasprun.xml that actually contains a Hessian.
#    These are the force-constant runs -- the valuable ones.
# ---------------------------------------------------------------------
echo
echo "--- [1] vasprun.xml files containing <hessian> --------------------"
echo "    (these can be converted directly with: phonopy --fc vasprun.xml)"
echo
find "$ROOT" -name vasprun.xml -type f 2>/dev/null | while read -r f; do
    if grep -q "hessian" "$f" 2>/dev/null; then
        sz=$(du -h "$f" | cut -f1)
        printf "  HESSIAN  %6s  %s\n" "$sz" "$f"
    fi
done

# ---------------------------------------------------------------------
# 2. INCARs that generated force constants
# ---------------------------------------------------------------------
echo
echo "--- [2] INCARs with IBRION = 5 or 6 (force-constant runs) ---------"
echo
find "$ROOT" -name INCAR -type f 2>/dev/null | while read -r f; do
    if grep -qE "^[[:space:]]*IBRION[[:space:]]*=[[:space:]]*[56]" "$f"; then
        echo "  $f"
        grep -HE "^[[:space:]]*(ENCUT|EDIFF|POTIM|NFREE|ISYM|ISMEAR|SIGMA|LREAL|IVDW|PREC|ADDGRID|ML_)" "$f" \
            | sed 's|^.*INCAR:|      |'
        echo
    fi
done

# ---------------------------------------------------------------------
# 3. INCARs that only read force constants (dispersion-only runs)
# ---------------------------------------------------------------------
echo "--- [3] INCARs with LPHON_READ_FORCE_CONSTANTS (dispersion only) --"
echo "    NOTE: in these runs every other INCAR tag is INERT."
echo
find "$ROOT" -name INCAR -type f 2>/dev/null \
    -exec grep -lie "LPHON_READ_FORCE_CONSTANTS" {} \; 2>/dev/null \
    | sed 's|^|  |'

# ---------------------------------------------------------------------
# 4. MLFF artefacts -- provenance of the trained model
# ---------------------------------------------------------------------
echo
echo "--- [4] MLFF files (ML_AB / ML_FF / ML_ABN / ML_FFN) --------------"
echo
find "$ROOT" \( -name "ML_AB*" -o -name "ML_FF*" \) -type f 2>/dev/null \
    | while read -r f; do
        printf "  %8s  %s\n" "$(du -h "$f" | cut -f1)" "$f"
      done

# ---------------------------------------------------------------------
# 5. Existing phonopy artefacts
# ---------------------------------------------------------------------
echo
echo "--- [5] phonopy files already present -----------------------------"
echo
find "$ROOT" \( -name "FORCE_CONSTANTS*" -o -name "FORCE_SETS" \
               -o -name "phonopy*.yaml" -o -name "BORN" \
               -o -name "band.conf" -o -name "mesh.conf" \
               -o -name "*.hdf5" \) -type f 2>/dev/null \
    | sed 's|^|  |'

# ---------------------------------------------------------------------
# 6. ENCUT audit -- the check that gates interpretation of the benchmark
# ---------------------------------------------------------------------
echo
echo "--- [6] ENCUT AUDIT ------------------------------------------------"
echo "    If the ab initio and MLFF branches used different ENCUT, the"
echo "    systematic -0.40 meV softening is partly a basis-set offset,"
echo "    not a force-field property."
echo
find "$ROOT" -name INCAR -type f 2>/dev/null | while read -r f; do
    e=$(grep -iE "^[[:space:]]*ENCUT" "$f" | head -1 | tr -d ' ')
    m=$(grep -icE "^[[:space:]]*ML_(MODE|LMLFF)" "$f")
    i=$(grep -iE "^[[:space:]]*IBRION" "$f" | head -1 | tr -d ' ')
    [ -n "$e" ] && printf "  %-12s %-12s mlff=%s  %s\n" "$e" "$i" "$m" "$f"
done | sort

# ---------------------------------------------------------------------
# 7. Detected symmetry in the force-constant runs
# ---------------------------------------------------------------------
echo
echo "--- [7] Space group detected by VASP in OUTCARs -------------------"
echo "    Trilayer AA'A h-BN should be P-6m2 (#187), point group D_3h."
echo
find "$ROOT" -name OUTCAR -type f 2>/dev/null | while read -r f; do
    sg=$(grep -m1 -iE "space group|Routine SETGRP|LATTYP" "$f" 2>/dev/null)
    nd=$(grep -c "Degrees of freedom" "$f" 2>/dev/null)
    if [ -n "$sg" ]; then
        echo "  $f"
        echo "      $sg"
        grep -m1 -i "Degrees of freedom" "$f" | sed 's|^|      |'
    fi
done

# ---------------------------------------------------------------------
# 8. Residual stress
# ---------------------------------------------------------------------
echo
echo "--- [8] Residual stress in OUTCARs (want |in-plane| < 0.5 kBar) ---"
echo
find "$ROOT" -name OUTCAR -type f 2>/dev/null | while read -r f; do
    s=$(grep "in kB" "$f" 2>/dev/null | tail -1)
    [ -n "$s" ] && printf "  %s\n      %s\n" "$f" "$s"
done

echo
echo "==================================================================="
echo " Done."
echo "==================================================================="
