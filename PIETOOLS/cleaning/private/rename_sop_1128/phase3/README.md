# Phase 3 – verification, 11/28/2026

Identity before and after the rename. A rename changes no arithmetic, so every number must match.

## Method

- After: export of 74177d39 (Phases 1, 1b and 2) at `C:/wt1128sop/after` (`git archive` without
  `cleaning/`, `orphaned_file/`, `.claude/`).
- Reference: the run of 11/27 on abfb0e83 (scratchpad `p3/out_new`, Phase 3 of the container
  rename). Between abfb0e83 and 91ea628e, the parent of the rename, only documents and records
  changed (the agent file, the two `sopvar` notes and their PDFs, `testing_codes/README.md`), so
  its toolbox code is that of 91ea628e. Control: two rows of that run, `lps_test_lpi_ineq_sop` and
  `demo01`, rerun on the export of 91ea628e (`C:/wt1128sop/base`), must reproduce it, so that the
  environment (MATLAB, MOSEK, SeDuMi) has not moved since 11/27.
- Harness: that of 11/27 (`../../rename_copvar_opvar_1125/phase3/`): `d3run.m`, whose resolution
  check now picks the names of each tree (`PIETOOLS_stability_sop`, `lpisettings_sop` before;
  `PIETOOLS_stability`, `lpisettings` after); the capture shadows of `Sedumi2Mosek`,
  `MosekSol2SedumiSol` and `sedumi` (SHA-256 of each array); `q_after.tsv` from its `q_new.tsv` by
  `mkq_after.pl` (folders and commands take the new names, the row ids stay); `p3s_drive.sh`,
  `p3s_finish.sh`; one MATLAB process per batch, none over 17 min.
- Order: T1 and D on the tree after the rename (03:14–03:24) and the control (03:24–03:25), within
  20 minutes; T2, P and X, about 35 minutes together, were queued for 04:30 and started at 03:51
  (maintainer: "start now"). The two P rows stopped at once (`p0_path` wants the marker
  `P0_COMMIT.txt` at the export root, which the export of 11/27 had); with the marker, rerun
  04:15–04:21.

## Results

- Rows (`results_after.tsv`, `cmp_sdp_all.txt`): 56 of 56 OK, each with the status of the
  reference (T1 19, T2 19, D 15, P 2, X 1).
- SDPs: every row has the same sequence of SDP records – 110 SDPs equal bit for bit in A, b, c and
  K, and in x and y wherever captured (219 records): the converter and LPI tests, the executives
  test, the manual snippet Ch7, the FI library, PIESOS, DEMO1–DEMO15 (DEMO4: 8 bisection steps) and
  the 3-D heat LPI (853 s; 1,090 s on 11/27).
- Control (`cmp_sdp_ctrl.txt`): both rows equal to the reference, 5 SDPs bit for bit in A, b, c,
  K, x, y.
- DEMO workspaces (`cmp_ws.txt`): 15 files, no differing leaf except three wall-clock times in
  DEMO13 (`info_cl.hist.t_build`, `t_solve`, `witness.t`), as on 11/27.
- F1, F2 (`cmp_f12.txt`, `p0_compare` of pass 7 of this export against pass 7 of 11/27): PIESIM on
  73 library examples, largest relative difference 0; the operators of 172 library PIEs applied to
  random polynomials at random points, largest relative difference 0.
- The guards of Phase 1b are not met by any of the 56 rows: no current call in the tests, the
  demos or the library passes a numeric third input to `getController`/`getObserver` or
  `opts.pure`/`opts.sep` to `lpi_ineq`.
