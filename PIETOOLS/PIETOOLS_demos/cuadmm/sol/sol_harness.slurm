#!/bin/bash
# sol_harness.slurm -- run bl_regime block(s) of the cuADMM harness on one Sol A100.   CC, 09/27/2026
#
#   sbatch --time=00:20:00 --mem=32G --export=ALL,BLOCKS="N1 S1" sol/sol_harness.slurm
#
# Size every submission from the desktop's measured phase costs (bl_bisect probe
# rows: t_init, t_s, t_cert, mem_mb) and ask the maintainer before any job above
# ~10 CHE (a 4-cpu / 32G / a100 job bills ~37 CHE per hour).  Defaults below are
# for a smoke job only.
#
# Sol facts this script encodes (memory: asu-sol-access-and-environment,
# sol-gpu-job-gotchas):
#   --gres=gpu:a100:1   explicit type: gpu:1 can land on a MIG slice (~2x slower,
#                       bills 15), and nvidia-smi names the physical A100 anyway;
#                       CUADMM_REQUIRE_FULL_GPU=1 makes bl_bisect refuse a slice
#   --mem >= 24G        Sol rejects GPU jobs asking for less than 24000 MB
#   gcc/13.4.0          the Sol cuadmm_exe needs its libstdc++ AT RUN TIME
#   matlab -batch "run('file.m')"   a multi-line -batch string can fail
#   $SLURM_SUBMIT_DIR   Slurm runs a spooled copy, so BASH_SOURCE is useless
#   no Mosek            only cuADMM (and SeDuMi) blocks can run here; Mosek
#                       brackets come from the desktop's regime/*.mat, copied in
#                       by sol_push.sh
#SBATCH --job-name=cuadmm_harness
#SBATCH --partition=htc
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/mpeet/pietools/harness_out/slurm/%x-%j.out
# (#SBATCH lines are not variable-expanded; sol_push.sh creates that folder)
set -euo pipefail
: "${BLOCKS:?pass --export=ALL,BLOCKS='<bl_regime block ids>'}"
SOL=/scratch/mpeet/pietools
HARNESS_ROOT="${HARNESS_ROOT:-$SOL/harness/PIETOOLS}"        # sol_push.sh puts the tree here
export CUADMM_OUT="${CUADMM_OUT:-$SOL/harness_out}"          # dumps + results (sol_push/sol_pull)
export CUADMM_EXE="${CUADMM_EXE:-$SOL/solvers/cuADMM/build_sol/cuadmm_exe}"
export CUADMM_SEDUMI="$SOL/solvers/SeDuMi_1_3"
export CUADMM_REQUIRE_FULL_GPU=1

module purge
module load gcc/13.4.0 nvhpc/25.7-cuda-12.9-nompi matlab/2025b
# MATLAB's system() prepends MATLAB's own libstdc++ (no GLIBCXX_3.4.32), so
# cuadmm_exe launched from MATLAB fails to load (Sol-0, job 64075050).  Save the
# module environment's path now; bl_bisect restores it for cuadmm_exe only.
export CUADMM_LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# a missing library would make cuadmm_exe fail with exit 0 (it does on CUDA
# failure): stop here instead of spending the walltime
nf=$(ldd "$CUADMM_EXE" | grep -c 'not found' || true)
[ "$nf" = "0" ] || { echo "cuadmm_exe: $nf unresolved libraries"; ldd "$CUADMM_EXE"; exit 2; }
nvidia-smi -L
nvidia-smi -L | grep -q MIG && { echo "MIG slice allocated; refusing"; exit 3; }

# deadline: the job's own end time minus 5 minutes, so bl_regime refuses items
# that cannot finish inside the walltime
end=$(scontrol show job "$SLURM_JOB_ID" | grep -oP 'EndTime=\K\S+')
export REGIME_DEADLINE=$(( $(date -d "$end" +%s) - 300 ))
echo "job $SLURM_JOB_ID on $(hostname): blocks [$BLOCKS], deadline $(date -d @$REGIME_DEADLINE)"

mkdir -p "$CUADMM_OUT/regime"
run="$CUADMM_OUT/regime/run_${SLURM_JOB_ID}.m"
{
  echo "addpath('$HARNESS_ROOT/PIETOOLS_demos/cuadmm'); cuadmm_path;"
  for b in $BLOCKS; do echo "bl_regime('$b');"; done
} > "$run"
matlab -nodisplay -nosplash -batch "run('$run')"
sacct -j "$SLURM_JOB_ID" --format=JobID,Elapsed,MaxRSS,TotalCPU,AllocTRES%60 || true
