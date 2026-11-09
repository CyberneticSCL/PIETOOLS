#!/bin/bash
# sol_pull.sh -- fetch Sol harness results into a local folder.   CC, 09/27/2026
#
#   bash sol/sol_pull.sh 'C:\Users\mpeet\cuadmm_out\sol_0927'
#
# Copies regime/ (manifest, per-item .mat, block logs, probe TSVs, env.txt, the
# Slurm .out files) but not the kept iterates (run_NNN_X/y.txt), which are
# large; fetch those by hand for re-certification.
set -euo pipefail
DSTW="${1:?local destination}"; DST="$(cygpath -u "$DSTW")"; mkdir -p "$DST"
SOL=/scratch/mpeet/pietools
ssh sol "cd $SOL && tar czf - --exclude='run_*_X.txt' --exclude='run_*_y.txt' --exclude='At.txt' \
         harness_out/regime harness_out/slurm harness/VERSION" | tar xzf - -C "$DST"
echo "pulled into $DST:"; ls "$DST"
