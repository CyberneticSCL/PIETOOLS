# Stress-test inventory A: tests and harnesses under `testing_codes/private/` (11/22/2026)

Queue for the stress test of 11/22/2026 (maintainer: "everything that can be run, should be. This
is a stress test."), scope `testing_codes/private/`: every `Test_*`, `test_*`, `*_validate`, `e2e_*`
file, the baseline harnesses of `testing_codes/private/baselines/` and the cuADMM suite. Nothing was
run to build it: two MATLAB probes of under 6 s checked entry points (`cuadmm_path`, `bl_suites`,
`bl_cases`, `heatNd_path`, `pielr_path`, the three `lb_setup` copies). Rows are in `inventory_A.tsv`
(415 rows, 12 columns); `runnable` is `yes` for 382 rows and `no` for 33.

## How a row runs

One batch MATLAB process per row, on a clean worktree of HEAD, as `rel_run2` of the release run of
11/19 (`C:\Users\mpeet\cuadmm_out\release_1119\runner`): `restoredefaultpath`; `cd` to the
worktree `PIETOOLS/`; `pietools_path_update`; `addpath(genpath(fullfile(pwd,'testing_codes','private')))`
(one warning for `private/`); SeDuMi
(`C:\YALMIP\yalmip_install\tbxmanager\toolboxes\sedumi\1.3\pcwin64`) and MOSEK 11
(`C:\Program Files\Mosek\11.0\toolbox\r2019b`); figures invisible; `rng(1)`; `cd(folder)`; evaluate
`cmd` in the base workspace. `folder` is relative to the worktree `PIETOOLS/`. Harness rows reset the
path themselves (`lb_setup` through `lb_run`, `p0_path`, `cuadmm_path`, `pielr_path`, `heatNd_path`).
A pass is no error, unless `check` names a comparison; several `cmd` add an `assert` to scripts that
only print (`check` says so).

`est_s` is the run time of the row without MATLAB start and path setup (about 20-40 s per process:
about 3 h over 382 processes). `cap_s` = max(300, 3 est_s). `mem_gb` is a recorded peak where one
exists, otherwise an estimate marked in `source`, otherwise empty. `source` names the timing record.

## Counts and time

| group | rows | runnable | est (h) |
|---|---|---|---|
| PIESOS (7 groups: top, intop, polyopvar, tensopmat, tensopvar, sosprog, quadPoly) | 118 | 118 | 1.2 |
| PIESOS_perf | 1 | 0 | – |
| SOPVAR | 16 | 16 | 0.03 |
| LPI_SOP | 15 | 14 | 0.46 |
| EXEC_SOP | 2 | 2 | 0.01 |
| CONVERTERS | 6 | 6 | 0.74 |
| HEATND | 16 | 15 | 1.7 |
| FI | 1 | 1 | 0.01 |
| LOWRANK | 4 | 3 | 3.1 |
| CUADMM | 30 | 22 | 0.82 |
| BL_1028 (example library 10/28) | 4 | 0 | – |
| BL_P0 (opvar Phase 0, 11/02) | 18 | 16 | 1.0 |
| BL_CLASSOPS (11/05) | 9 | 6 | 0.33 |
| BL_MIXED (11/12) | 1 | 0 | – |
| BL_TIERS (11/17) | 12 | 8 | 8.7 |
| BL_STUDIES (11/18) | 162 | 155 | 37.3 |
| **total** | **415** | **382** | **55.3** |

14 runnable rows of BL_STUDIES are duplicates of BL1117_tiers_A/B rows (control tags L1ctl, M1L,
H1ctl, B1c0, M2L, H2ctl, B4ctl, C1ext, C1ext2, E1str, E2l, E2h, E3vh) or redundant (B1c1g); `reason`
marks them. Without them: 168,593 s = **46.8 h**, plus process overhead.

Longest rows: BL1117_tiers_A 21,909 s (light, heavy, veryheavy library pass; final_summary_1118.md);
ST_A3cap1, ST_A3cap1f09, ST_A3cap1f11 9,000 s each (D2 cap-1 rung, 2.72 M variables, m 82,551;
estimate 2-3 h near 60 GB, never completed); LR_pielr_bench_t3_2d 7,200 s (time to the cell cap of
09/2026, not to completion); BL1117_tiers_B 7,166 s; then ST_H2ctl 5,877 s (duplicate), ST_H2T3
5,692 s, ST_EXS2 5,484 s, ST_B4q866 5,422 s, ST_E2l 5,264 s (duplicate).

Memory at or above 10 GB: ST_A3cap1* about 60 (estimate; run alone); BL1117_tiers_A/B up to 40
(size gate 2e6 variables, m 45,000); ST_A3cap1ci 38.0; ST_H2T3 22.0; ST_H2ctl 21.7; ST_C1ext2 18.5;
ST_C2ext2 14.9; LPI_test_endtoend_E4 12 (killed 11/19 at free memory below 3 GB); HND_pv_3d_Edir1 12;
CU_check_scaling 12; HND_val_3d 11.5; HND_pv_3d_L2dir1 10.

## Setup the rows assume

- **Phase 0 (BL_P0).** `git archive HEAD:PIETOOLS` into `C:\Users\mpeet\wt\p0_stress1122` and
  `git rev-parse HEAD > ...\P0_COMMIT.txt` (`p0_path` refuses a root without it). The pass-3 data
  of F1, F2 and F4 are git-ignored and exist only in the main checkout: copy
  `F1/pass3`, `F2/pass3`, `F4/pass3` into the worktree harness folder before P0_CMP. Rows record
  pass 5; P0_CMP runs last.
- **Studies (BL_STUDIES).** Each `cmd` writes 1440 into `ST_STOPMIN` (default 15: no row starts after
  15 min) and moves the banked `results_<TAG>.csv`, `exs_<TAG>*.csv` and `exs_state*/` to `banked_*`,
  since `lb_run` and `ex_search` resume and would skip every banked row. Compare with
  `bash cmp_ident.sh banked_results_<TAG>.csv results_<TAG>.csv` and `perl cmp_runs.pl` in the same
  folder. ST_libcheck writes a fresh `libcheck_stress1122.csv`.
- **cuADMM.** `CUADMM_OUT` points at `tempdir/stress_1122/cuadmm`: `bl_run` and `bl_scale` resume from
  TSVs there. `bl_check` never re-banks here (mode check).
- **heatNd.** Every validation writes under `tempdir/stress_1122/heatNd/<stage>`; the '3d' stages
  resume from a trace in their folder, so the folder must be new.
- **Order.** CU_blrun_hinf/_estimator/_controller before CU_bl_bisdump; CU_regime_B0 before
  CU_regime_B1a; ST_BLOCKCHECK before ST_BLOCKCHECK2-5; ST_B1c1 before ST_B1c1g; P0_F1-P0_F4 before
  P0_CMP. ST_A3cap1*, ST_A3cap1ci, BL1117_tiers_A/B and LPI_test_endtoend_E4 run alone.

## Known failures recorded before this run

- LPI_SOP (release run 11/19, origin 7f7a0300): `Unrecognized field name "int"` in
  test_translators_sop, test_dpvar_hinf_chain_sop, test_getsol_sop, test_settings_prog_sop and
  test_endtoend E1/E2 (settings read `lpi.pos.int`; the tiers hold `len_*` storage rows since
  bb367fe1); test_dpvar_ops_sop 3 of 575 checks (stale, closed 11/19); E4 OK 10/27 (723 s), killed
  11/19 for memory.
- HEATND: test_heatNd_pie, 8 rows "(g) vs stock convert" fail (stale, closed 11/19).
- CUADMM: `bl_expect.tsv` last banked 11/08 (2c8f7143); on the 11/18 tree (= bb367fe1)
  STRUCTURE CHANGED in structure 5, twod 4, objective2d 3, scaling 4 (so smoke 1); executives
  expects ctrl_rd1 FAIL (no point); nonlinear 3 ERR (`fisher_prog` broken by 6fd03c15, no port).
- PIESOS: no recorded run since 10/23 ("Skip PIESOS testing for now"). Recorded: Test_polyopvar_subs
  2 of 12 cases, Test_convert_PDE_nonlinear_varnames case 1 (triage 11/12, dropped 11/19);
  inferred: Test_piesos_program line 254 calls the removed `ndopvar`.
- BL_CLASSOPS: `bench_lpi_degree` builds `lpi.pos` in the 11/12 layout; likely refused since
  bb367fe1 (inferred, not run).
- BL_TIERS: 11/17 error rows per tier in the `known_fail` column (PDE 50, 53-56, 61, 65, 69, 70).
- BL_STUDIES: banked files; libcheck 11/19 cases 7 and 11 dual infeasible, case 50 error;
  blockcheck DDE19 Hinf_control unknown, DDF03 not found or refused; every A3 rung infeasible.
- LOWRANK: tier 1 missed 3 of 13 (search failures); the 2-D low-rank arm certifies nothing.
- BL_P0: F1 1-D example 32 errors in every pass; F3 skips PDE 36, 38, 50 and its dual.

## Not runnable (33 rows)

- **Retired classes or code** (opvar, dopvar, opvar2d, `lpi_programming`, stock executives,
  `lpisettings`) – 11: LR_run_regression_1d (low-rank `tests_1d`), LPI_test_endtoend_E3,
  CO_bench_class_ops, CO_analyze_class_ops, CO_lpideg_legacy, BL1028_* (4: `lb_setup` fails at
  `lpisettings`, measured), BL1117_hd_diag, BL1117_lb_full.
- **GPU, WSL or Sol** – 8: bl_check lists sol_compare, sol_large, problematic; bl_regime GPU
  blocks; `gpu/*.sh`; `sol/*`; bl_verify/cuimport; ST_A3unc (57.2 GB at the memory stop, a Sol job).
- **Data or trees that no longer exist** – 4: BL1112_runners (deleted scratch exports),
  BL1117_ab_head (worktree of 33d30877), BL_P0_run_phase0_sh (deleted exports), HND_pv_3dsum
  (needs the certificate files of a previous '3d' run by name).
- **Main checkout by absolute path** – 2: ST_VDA, ST_st_merge_tests (its items run as separate rows).
- **Other** – 8: PERF_mtimes_tensopvar (driver launches R2026b, not installed; no spec file),
  CU_bl_reduce (library function, no driver), ST_NEWCAT and ST_CATFIX (catalog edits already
  applied), ST_LIBCHECK_subsets, ST_lb_validate_copy and BL_P0_lb_validate (subsets or copies),
  BL1117_lb_rerun_5356 (superseded).

No item needs an interactive prompt: no file under `testing_codes/private/` calls `input`,
`keyboard` or `pause`.

## Outside this inventory

The settings-rerun queue of 11/21 (`cleaning/private/settings_rerun_1121/rr_run.sh`, tags RR4_*,
shift checks, libcheck_rr) drives this harness from `cleaning/` and hardcodes `C:/wt1121r`; it is not
listed. `PIETOOLS_demos/cuadmm/heat2d_compare.m` and DEMO11 (from `st_merge_tests`) belong to the
demos inventory.
