#!/bin/bash
#SBATCH --job-name=c5c_box3d
#SBATCH --partition=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=320G
#SBATCH --time=10:00:00
#SBATCH --output=/scratch/mpeet/pietools/psatz_faces/jobs/%x_%j.out
# The 3-D box set at the next cap (psatz_faces C5c), one MOSEK solve from
# MATLAB; billed 64 + 0.25*320 = 144 CHE/h. SCRATCH, MMP, 11/15/2026.
set -u
B=/scratch/mpeet/pietools/psatz_faces
echo "=== node $(hostname) job $SLURM_JOB_ID cpus $SLURM_CPUS_PER_TASK mem ${SLURM_MEM_PER_NODE:-?}M started $(date -Is)"
module purge
module load matlab/2025b 2>/dev/null || echo "module load failed; using the install path"
command -v matlab >/dev/null 2>&1 || export PATH=/packages/apps/matlab/r2025b/bin:$PATH
echo "=== matlab: $(command -v matlab)"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MATLAB_PREFDIR=${TMPDIR:-/tmp}/mlprefs_$SLURM_JOB_ID; mkdir -p "$MATLAB_PREFDIR"
export MOSEKLM_LICENSE_FILE=$HOME/mosek/mosek.lic
cd $B/psatz_faces
matlab -nodisplay -nosplash -batch "maxNumCompThreads($SLURM_CPUS_PER_TASK); run('$B/psatz_faces/pf_sol_c5c.m')"
rc=$?
echo "=== matlab exit $rc finished $(date -Is)"
exit $rc
