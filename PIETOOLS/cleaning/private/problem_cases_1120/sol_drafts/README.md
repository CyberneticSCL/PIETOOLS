# Sol job drafts – NOT submitted (maintainer, 11/20 16:34: "Not now")

Both exceed 10 CHE and wait for the maintainer's word (standing rule). CHE per hour = cores + 0.25 * GB.

| Job | Program | Resources | Cap |
|---|---|---|---|
| sol55b | PDE 55, 'P' stripped (dD 1): 7,340,048 variables, m 105,717 | general, 64 cores, 256 GB, 48 h | 6,144 CHE |
| sol56 | PDE 56, 'Q' light dD 0 with the cap, deg_Q{3} = [1 2 2]: 7,857,280 variables, m 135,522 | general, 64 cores, 384 GB, 48 h | 7,680 CHE |

Basis: job 2067511 (11/19-11/20, sol55) ran the 55 program for 6 h at 64 threads (191 GiB peak, about 327
core-hours of CPU) and reported no interior-point iterate; the run time needed is unknown. For m = 105,717 one
dense Cholesky of the Schur complement is about m^3/3 = 3.9e14 flops; 56 (m = 135,522) about 8.3e14 per
iteration. The job log of sol55 stayed empty until the limit: add a MOSEK log file (or flush the MATLAB output)
before resubmitting, so a partial run reports its iterations.

Setup: the repository copy /scratch/mpeet/pietools/pde55/repo of 11/19 (sol55); copy sol56.m and sol56_starts.csv
into its harness folder testing_codes/private/baselines/2026-11-18_studies/private, and the patched ex_search.m of
this run (results/aad9606c/ex_search.m) if the copy predates it.
