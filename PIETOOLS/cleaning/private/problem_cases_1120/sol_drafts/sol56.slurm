#!/bin/bash
#SBATCH --job-name=pde56
#SBATCH --partition=general
#SBATCH --cpus-per-task=64
#SBATCH --mem=384G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/mpeet/pietools/pde55/jobs/%x_%j.out
# DRAFT, NOT SUBMITTED. PDE 56 'Q' light dD 0 with the cap at deg_Q{3} = [1 2 2] (the degree that certifies its 1-D
# twin, PDE 43): 7,857,280 variables, m 135,522; cap 48 * (64 + 0.25*384) = 7,680 CHE.
set -u
echo "=== node $(hostname) job $SLURM_JOB_ID started $(date -Is)"
module purge
module load matlab/2025b || export PATH=/packages/apps/matlab/r2025b/bin:$PATH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MATLAB_PREFDIR=${TMPDIR:-/tmp}/mlprefs_$SLURM_JOB_ID; mkdir -p "$MATLAB_PREFDIR"
export MOSEKLM_LICENSE_FILE=$HOME/mosek/mosek.lic
cd /scratch/mpeet/pietools/pde55/jobs
matlab -nodisplay -nosplash -batch "maxNumCompThreads($SLURM_CPUS_PER_TASK); run('/scratch/mpeet/pietools/pde55/repo/testing_codes/private/baselines/2026-11-18_studies/private/sol56.m')"
echo "=== matlab exit $? finished $(date -Is)"
