# results/2026-09-24

The measured tables and raw logs behind `../../BASELINE_REPORT.md`. All of it was produced on
2026-09-24/25 by the harness in its original scratchpad layout, on one workstation
(i9-14900KF, 64 GB, RTX 4090). Frozen: nothing re-runs into this folder.

## Tables

| file | produced by | notes |
|---|---|---|
| `bl_mosek.tsv` | `bl_run` | Mosek arm, 35 cases. **The pre-`lpi_eq`-fix reference** that `experiments/regress_lpieq.m` compares against |
| `bl_mosek_v1.tsv` | `bl_run`, first pass | superseded; its `psd_*` columns sliced the wrong vector (`RRx`, not the cone vector `x`) |
| `bl_expect.tsv` | derived from `bl_mosek.tsv` | snapshot of the package-root `bl_expect.tsv`; **edit that one, not this** |
| `bl_scale.tsv` | `bl_scale` | the ladder's shapes up to n=32 (shape only above n=16) |
| `bl_big.tsv` | `bl_big` | n=24 and n=32, built and dumped after the `lpi_eq` fix |
| `bl_gpu.tsv` | `gpu/bl_gpu.sh` | **raw and contaminated**: two instances ran concurrently over part of the arm, so 19 of 70 (case, tolerance) keys appear twice. The first line is a data row, not a header |
| `bl_gpu_clean.tsv` | a shell loop over the per-case cuADMM logs | one row per (case, tolerance), with convergence decided by the presence of a `CUADMM_TIMING` line. The master table's GPU columns come from here. Rebuild it from `logs/cuadmm/` |
| `bl_gpu_last.tsv` | `gpu/bl_gpu_post.sh` | last-iteration residuals and which one was binding (`limiter`) |
| `bl_verify.tsv` | `bl_verify` | cuADMM certificates re-scored against the pre-dump data |
| `bl_bisect.tsv` | `bl_bisdump` | Mosek on the pinned-γ feasibility programs. **The `h2oco_rd1_*` rows are invalid** (γ pinned instead of γ²; see the report's §7 retraction) |
| `bl_bisgpu.tsv` | `gpu/bl_bisgpu.sh` | per-step times are valid; **every residual column is empty**, because of an `awk` bug that's now fixed. The pinf values the report cites were read from `logs/bisect/` |
| `check_smoke.tsv`, `check_structure.tsv` | `bl_check` | 4/4 and 22/22 PASS; `structure` is what cleared the `lpi_eq` fix |

## Logs

| folder | contents |
|---|---|
| `logs/*.log` | MATLAB output of each script and runner, by script name. `hinf_gam.log` is the output of what is now `experiments/cu_hinf_gam.m`, renamed to avoid a Robust Control Toolbox clash |
| `logs/cuadmm/<id>_<tol>.log` | full cuADMM output per case and tolerance: iteration trace, final residuals, `CUADMM_TIMING`, `CUADMM_MEM` |
| `logs/bisect/<lab>_step5000.log` | cuADMM traces on the pinned-γ programs at a 5000-iteration budget, the source of §7's plateau and 2-D control tables. `<lab>` is `<case>_g<pct>` for fractions of γ\*, `hinf2_rd_a<pct>` for absolute γ = pct/100, and `hinf2_rd_b<pct>` for the below-analytic controls |

The SDP dumps themselves (~1.9 GB) were not kept. Re-running the harness regenerates them
under `cuadmm_outdir()`.
