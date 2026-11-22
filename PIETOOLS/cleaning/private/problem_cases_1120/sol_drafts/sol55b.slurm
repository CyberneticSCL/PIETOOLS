#!/bin/bash
#SBATCH --job-name=pde55b
#SBATCH --partition=general
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/mpeet/pietools/pde55/jobs/%x_%j.out
# DRAFT, NOT SUBMITTED. PDE 55 'P' stripped as in sol55 (job 2067511, 6 h limit reached inside the solve), on the
# general partition with 48 h; cap 48 * (64 + 0.25*256) = 6,144 CHE.
set -u
echo "=== node $(hostname) job $SLURM_JOB_ID started $(date -Is)"
module purge
module load matlab/2025b || export PATH=/packages/apps/matlab/r2025b/bin:$PATH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MATLAB_PREFDIR=${TMPDIR:-/tmp}/mlprefs_$SLURM_JOB_ID; mkdir -p "$MATLAB_PREFDIR"
export MOSEKLM_LICENSE_FILE=$HOME/mosek/mosek.lic
cd /scratch/mpeet/pietools/pde55/jobs
matlab -nodisplay -nosplash -batch "maxNumCompThreads($SLURM_CPUS_PER_TASK); run('/scratch/mpeet/pietools/pde55/repo/testing_codes/private/baselines/2026-11-18_studies/private/sol55.m')"
echo "=== matlab exit $? finished $(date -Is)"
