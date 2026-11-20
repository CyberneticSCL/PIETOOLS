#!/bin/bash
#SBATCH --job-name=pde55
#SBATCH --partition=public
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/mpeet/pietools/pde55/jobs/%x_%j.out
# PDE 55 in 'P' stripped, one build and one solve (sol55.m). Projected 3-4 h
# and 100-200 GB; the 6 h limit caps the charge at 768 CHE (64 + 0.25*256 per h).
set -u
echo "=== node $(hostname) job $SLURM_JOB_ID cpus $SLURM_CPUS_PER_TASK mem 256G started $(date -Is)"
module purge
module load matlab/2025b || echo "module load failed; using the install path"
command -v matlab >/dev/null 2>&1 || export PATH=/packages/apps/matlab/r2025b/bin:$PATH
echo "=== matlab: $(command -v matlab)"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MATLAB_PREFDIR=${TMPDIR:-/tmp}/mlprefs_$SLURM_JOB_ID; mkdir -p "$MATLAB_PREFDIR"
export MOSEKLM_LICENSE_FILE=$HOME/mosek/mosek.lic
cd /scratch/mpeet/pietools/pde55/jobs
matlab -nodisplay -nosplash -batch "maxNumCompThreads($SLURM_CPUS_PER_TASK); run('/scratch/mpeet/pietools/pde55/repo/testing_codes/private/baselines/2026-11-18_studies/private/sol55.m')"
echo "=== matlab exit $? finished $(date -Is)"
