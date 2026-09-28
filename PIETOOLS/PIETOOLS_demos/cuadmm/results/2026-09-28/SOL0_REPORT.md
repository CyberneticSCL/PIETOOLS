# Sol-0: the harness on Sol, first run (2026-09-28)

Blocks N1 (driver checks) and X0 (a discarded warm-up, then three feasibility runs at the certified
standard: tol 1e-7, psd_clip), on the dumps pushed from the workstation. The desktop reference is
the same blocks on the same bytes: `results/2026-09-27/x0_desktop_*`.

## Result

| case | m | desktop RTX 4090: it / solve | Sol A100-80GB: it / solve | speedup | verdict (both) |
|---|---|---|---|---|---|
| `stab_tr1` | 96 | 1,973 / 1.57 s | 1,973 / 0.97 s | 1.62× | F |
| `stab_rd1_hv` | 135 | 14,092 / 19.02 s | 14,092 / 8.55 s | 2.22× | F |
| `scale_rot_f0p50_n04` (coupled) | 2,160 | 10,894 / 52.17 s | 10,894 / 27.01 s | 1.93× | F |

- **Iteration counts are identical to the digit, and so are the verdicts.** The repaired
  λ_min is -5.0e-8 on `stab_tr1` on both machines.
- **N1 on Linux: 4/4 PASS**, including STOP inside `bl_bisect`.
- **The GPU was a full A100-SXM4-80GB, not a MIG slice** (`sol0_env.txt`; the MIG check passed).
- The Sol binary is a different build from the desktop's: different sha256, fat for
  sm_80/89/90. Identical iteration counts on all three cases say it is numerically the same solver
  at these sizes.
- **Warm-up matters:** the first `stab_tr1` solve took 1.57 s, the second 0.97 s.
- **GPU memory:** about 810 MiB for all three (`CUADMM_MEM`).
- **Cost:** job 64075208, 1 min 39 s at billing weight 40, about 1.1 CHE. The whole Sol-0 exercise
  came to about 2 CHE (below).

## The failure it found first (job 64075050)

Every cuADMM call died at load:

    .../matlab/r2025b/sys/os/glnxa64/libstdc++.so.6: version `GLIBCXX_3.4.32' not found

MATLAB's `system()` prepends MATLAB's own libstdc++ for child processes, which shadows
gcc/13.4.0's.
- **Fix:** `sol_harness.slurm` saves the module environment's path as `CUADMM_LD_LIBRARY_PATH`
  before MATLAB starts, and the Linux launcher in `bl_bisect` restores it for `cuadmm_exe` only.
- **Caught because** the harness gates on the `CUADMM_TIMING` line: `cuadmm_exe` exits 0 on
  failure. Each failed call was read as i, never as F.
- **Resume trap:** a failed solver call still completes as an `ok` item, so the Sol manifest was
  cleared before the rerun. The first attempt's files are in `harness_out/attic_sol0_try1` on Sol.
- **Cost:** about 0.8 CHE (1 min 14 s), plus a cancelled resubmission at about 0 CHE.

## What it establishes, and what it does not

**Established:**
- The harness runs unmodified on Sol from a pushed tree: launcher, environment record, MIG check,
  deadline from the job's EndTime, results pulled back.
- Desktop and Sol agree exactly at small m.
- The A100 is 1.6–2.2× faster per solve here, consistent with the 09-21 raw-binary measurement.

**Not established:**
- Behaviour at large m, where iteration counts drifted about 1% on Th_n3 in 09-21's raw runs.
- The cost of the host phases on Sol.
- Anything about the coupled-ladder fill problem (see 2026-09-27 §4).
