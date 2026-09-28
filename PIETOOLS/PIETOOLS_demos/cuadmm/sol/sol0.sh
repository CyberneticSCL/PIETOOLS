#!/bin/bash
# sol0.sh -- the Sol-0 smoke job (approved by the maintainer 2026-09-27).   CC, 09/27/2026
#
#   bash sol/sol0.sh            (needs the ASU VPN)
#
# Pushes the tree and the three X0 dumps from the workstation's x0_desktop run,
# then submits bl_regime blocks N1 (driver checks) and X0 (warm-up, then
# stab_tr1, stab_rd1_hv, scale_rot_f0p50_n04 at the certified standard) on one
# pinned A100.  Desktop reference, same dumps, measured 09-27 04:21:
#   stab_tr1 1,973 it 1.57 s F;  stab_rd1_hv 14,092 it 19.0 s F;
#   scale_rot_f0p50_n04 10,894 it 52.2 s F (clipped eta 1.0e-8); N1 all PASS.
# Estimated cost: ~2-4 min of a 4-cpu / 32G / a100 job at ~37 CHE/h, i.e.
# ~2-3 CHE (billing is on elapsed time; the 15 min walltime is a ceiling).
set -euo pipefail
PKG="$(cd "$(dirname "$0")/.." && pwd)"
bash "$PKG/sol/sol_push.sh" 'C:\Users\mpeet\cuadmm_out\x0_desktop' stab_tr1 stab_rd1_hv scale_rot_f0p50_n04
ssh sol "cd /scratch/mpeet/pietools/harness/PIETOOLS/PIETOOLS_demos/cuadmm && \
         sbatch --time=00:15:00 --export=ALL,BLOCKS='N1 X0' sol/sol_harness.slurm"
echo "watch: ssh sol 'squeue -u mpeet; tail -20 /scratch/mpeet/pietools/harness_out/slurm/*.out'"
echo "fetch: bash $PKG/sol/sol_pull.sh 'C:\\Users\\mpeet\\cuadmm_out\\sol0'"
