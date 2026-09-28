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

## X1: `scale_stab_n32` on Sol (job 64075580, 5 min 13 s, about 3.5 CHE)

A warm-up, then the decoupled rung n = 32 (m = 138,240, 1.42M nnz), run exactly as the
workstation's 09-27 S1 did it (`results/2026-09-27/probes/scale_stab_n32_S1_cuadmm_probes.tsv`):

| | desktop RTX 4090 (09-27) | Sol A100-80GB (09-28) | ratio |
|---|---|---|---|
| iterations | 14,092 | **14,092** | identical |
| verdict, clipped η | F, 6.67e-9 | **F, 6.67e-9** | identical |
| solve | 354.8 s (25.2 ms/it) | **253.0 s (18.0 ms/it)** | **1.40× faster** |
| setup (`t_init`) | 0.4 s | 0.59 s | 1.5× slower (host) |
| certification (`t_cert`) | 4.2 s | 4.63 s | 1.1× slower (host) |
| GPU memory | — | 988 MiB | |
| MATLAB peak memory | 2.3 GB (current use, Windows) | 1.26 GB (VmHWM); sacct MaxRSS 2.23 GB | |

**The speedup shrinks with size: 1.6–2.2× at m ≤ 2,160, 1.40× at m = 138,240.** This matches the
09-22 raw-binary measurement on Th_n3 (1.45× at m = 80,550) and the memory-bandwidth reading
(asu-sol-access-and-environment): the A100's advantage does not grow toward its FP64 ratio as the
problem grows. The host phases run 1.1–1.5× slower on Sol, as measured before. At this size they
are seconds, so they don't matter here, but for Th_n3 they are about 13 minutes per probe
(2026-09-27 §3).

## What it establishes, and what it does not

**Established:**
- The harness runs unmodified on Sol from a pushed tree: launcher, environment record, MIG check,
  deadline from the job's EndTime, results pulled back.
- Desktop and Sol agree exactly at small m.
- The A100 is 1.6–2.2× faster per solve here, consistent with the 09-21 raw-binary measurement.

**Not established:**
- ~~Behaviour at large m~~: X1 shows identical iterations and verdict at m = 138,240 on the
  decoupled ladder. Only a coupled or large-nnz case (Th_n3 drifted about 1% in 09-21's raw
  runs) can still differ.
- Host-phase cost on Sol for a large-nnz dump.
- Anything about the coupled-ladder fill problem (see 2026-09-27 §4).
