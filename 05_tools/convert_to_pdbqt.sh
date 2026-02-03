#!/bin/bash

module load openbabel 2>/dev/null || true
module load parallel 2>/dev/null || true

SDFDIR="ligands"
PDBQTDIR="ligands_pdbqt"
mkdir -p "$PDBQTDIR"

find "$SDFDIR" -type f -name "*.sdf" -print0 |
  parallel -0 -j "${SLURM_CPUS_PER_TASK:-12}" --bar '
    in="{}"
    base=$(basename "$in" .sdf)
    out="'"$PDBQTDIR"'/${base}.pdbqt"

    [[ -s "$out" ]] && exit 0

    obabel -isdf "$in" -opdbqt -O "$out" \
      --partialcharge gasteiger \
      -xh
  '
