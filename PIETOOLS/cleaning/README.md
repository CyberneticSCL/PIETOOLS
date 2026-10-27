# cleaning/

Files staged for removal before the public release (maintainer decision, 10/23/2026). Nothing has
been deleted: every file below is a `git mv` of a tracked file, and `git mv <moved path> <original
path>` restores it.

The moved files sit in `private/`, which `genpath` (and so `pietools_path_update`) skips with
everything below it; MATLAB never searches the subfolders of a private folder. None of them is on
the MATLAB path, which matters: `dead_code/` holds an older `rand_sdopvar.m` that would otherwise
shadow the live one in `sopvar/Testfolder/sdopvar/`. Do not add these folders by hand. Above each
moved folder, `@X` became `X_class` and `private` became `private_files`, so adding a subfolder to
the path cannot create a second class or private definition.

- `private/definitely/` – the Claude test folders (`claude_tests`, `claude`, `Testfolder_Claude`),
  a backup of unpushed commits, session notes, to-do lists and `CLAUDE.md`.
- `private/dead_code/` – obsolete code with no caller anywhere in the tree on 10/23/2026: a word search for
  each function name over every tracked file, comments excluded, plus the files of
  `PIETOOLS_PDE_GUI.mlapp`. It includes the earlier `sopvar/private/dead_code/` (09/29/2026), whose
  own README lists the original paths of its 31 files.

## Not moved

| Files | Reason kept |
|---|---|
| 12 helpers in `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/` (`cx_exec/cx_PIE2PDEstability`, `cx_plant`, `cx_resid`, `cx_settings`, `cx_shape`, `cx_solve`, `cx_space_list`, `cx_stability_2D`, `cx_stability_eq_opts_2D`, `cx_stability_pos2d`; `pi_dpvar_eval`, `pi_poly_grid`) | called by 9 tests in `lpi_programming_sopvar/tests/` and by `PIETOOLS_demos/sopvar_demos/heatNd_solve.m` and `heatNd_path.m` |
| `PIETOOLS/converters/internal/legacy/` (8 files) | called by `PIETOOLS_DDE.m`, `convert_PIETOOLS_PDE.m` and `initialize_PIETOOLS_PDE.m` |
| `PIETOOLS/opvar/2D/@opvar2d/*_legacy.m`, `@dopvar2d/*_legacy.m` (6 files) | called by the 2-D H-infinity executives |
| `PIETOOLS/opvar/inv_opvar_old.m` | called by `inv_opvar.m` |
| `PIETOOLS/opvar/inv_opvar_2.m` | the reference of `sopvar/Testfolder/Test_copvar_inv.m` |
| `PIETOOLS/SOSTOOLS400/multipoly/s2p_old.m`, `@polynomial/p2s_old.m`, both copies of `peval.mexw32` and `peval.c` | SOSTOOLS is a separate toolbox |
| `PIETOOLS/PIETOOLS_demos/cuadmm/`, `PIETOOLS/PIETOOLS_demos/lowrank_2d_stability/` | left alone for now (maintainer, 10/23/2026) |

Not yet decided (candidates for `maybe/`): `sopvar/sopvar_lift_notes.tex` and `.pdf`,
`sopvar/sopvar_implementation_notes.pdf`, `sopvar/sopvar_duality_pruning_notes.tex`,
`ndopvar/Notes_representation_ndopvar_v5.pdf`, and `PIETOOLS/PIETOOLS_2025_User_Manual.pdf`, which
master added on 04/29/2026 and which is an older copy of
`PIETOOLS_documentation/PIETOOLS_2025_User_Manual.pdf`.

## Moved files

First commit is the earliest commit of the original path, renames followed.

| Original path | Moved to (under `PIETOOLS/cleaning/private/`) | First commit |
|---|---|---|
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_H2_control.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_H2_control.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_H2_estimator.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_H2_estimator.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_H2_norm_2D_c.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_H2_norm_2D_c.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_H2_norm_2D_o.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_H2_norm_2D_o.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_H2_norm_c.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_H2_norm_c.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_H2_norm_c_coercive.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_H2_norm_c_coercive.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_H2_norm_o.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_H2_norm_o.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_H2_norm_o_coercive.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_H2_norm_o_coercive.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_control.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_control.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_estimator.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_estimator.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_gain.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_gain.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_gain_2D.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_gain_2D.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_gain_2D_non_coercive.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_gain_2D_non_coercive.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_gain_coercive.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_gain_coercive.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_gain_dual.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_gain_dual.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_gain_dual_2D.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_gain_dual_2D.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_gain_dual_coercive.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_Hinf_gain_dual_coercive.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_PDEstability.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_PDEstability.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_PDEstability_dual.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_PDEstability_dual.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_PIE2PDEstability_dual.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_PIE2PDEstability_dual.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_bisect.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_bisect.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_cases_h2.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_cases_h2.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_cases_hinf.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_cases_hinf.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_cases_stability.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_cases_stability.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_compare.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_compare.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_h2_pos2d.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_h2_pos2d.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_h2_set2d.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_h2_set2d.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_h2_struct2d.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_h2_struct2d.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_h2_trace.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_h2_trace.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_2dpos.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_2dpos.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_capture.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_capture.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_eqopts2d.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_eqopts2d.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_lf2d.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_lf2d.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_op2d.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_op2d.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_qdeg2d.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_qdeg2d.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_set2d.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_set2d.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_slack2d.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_slack2d.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_struct2d.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_hinf_struct2d.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_on_registry.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_on_registry.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_run_1d.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_run_1d.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_stability_check.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_stability_check.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_stability_dual_2D.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_stability_dual_2D.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_stability_thresh.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_stability_thresh.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_stock_capture.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_stock_capture.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_stock_pinned.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_stock_pinned.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_well_posedness.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_well_posedness.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/mincap_dual.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/mincap_dual.m` | MATTHEW M PEET, 2026-10-08, 96ef81f6 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/pi_blk_kernels.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/pi_blk_kernels.m` | MATTHEW M PEET, 2026-09-22, a252e720 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/pi_gamma_index.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/pi_gamma_index.m` | MATTHEW M PEET, 2026-08-23, e928b130 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/pi_kernel_pairing.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/pi_kernel_pairing.m` | MATTHEW M PEET, 2026-08-23, e928b130 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/pi_monom_vector.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/pi_monom_vector.m` | MATTHEW M PEET, 2026-08-23, e928b130 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/pi_sdopvar_kernels.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/pi_sdopvar_kernels.m` | MATTHEW M PEET, 2026-08-23, e928b130 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/slack_sizing/hinf_tailor_1d.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/slack_sizing/hinf_tailor_1d.m` | MATTHEW M PEET, 2026-10-08, bc3be32f |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/slack_sizing/hinf_tailor_2d.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/slack_sizing/hinf_tailor_2d.m` | MATTHEW M PEET, 2026-10-08, bc3be32f |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/slack_sizing/volterra_tailor.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/slack_sizing/volterra_tailor.m` | MATTHEW M PEET, 2026-10-08, bc3be32f |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_canonical_multiplier.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_canonical_multiplier.m` | MATTHEW M PEET, 2026-08-30, 1b1be818 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_container_method_parity.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_container_method_parity.m` | MATTHEW M PEET, 2026-09-30, 3e61d8dc |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_container_plus_batch.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_container_plus_batch.m` | MATTHEW M PEET, 2026-10-02, c1e1f2ec |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_copvar2opvar2d_order.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_copvar2opvar2d_order.m` | MATTHEW M PEET, 2026-09-28, e6c78800 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_copvar_blockops.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_copvar_blockops.m` | MATTHEW M PEET, 2026-09-25, 2fac5f08 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_copvar_kyp.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_copvar_kyp.m` | MATTHEW M PEET, 2026-09-25, 2fac5f08 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_copvar_registry.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_copvar_registry.m` | MATTHEW M PEET, 2026-09-26, bf1fd5f3 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_copvar_silent_fixes.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_copvar_silent_fixes.m` | MATTHEW M PEET, 2026-09-29, bd8f4272 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_dopvar2sdopvar.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_dopvar2sdopvar.m` | Declan, 2026-10-05, 6d39aee9 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_dpvar2sdvar.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_dpvar2sdvar.m` | Declan, 2026-10-05, 6d39aee9 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_dpvars2sdvars.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_dpvars2sdvars.m` | Declan, 2026-10-05, 6d39aee9 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_index_maps.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_index_maps.m` | MATTHEW M PEET, 2026-09-30, 3e61d8dc |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_lpi_eq_sdopvar.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_lpi_eq_sdopvar.m` | MATTHEW M PEET, 2026-08-29, 725e5f23 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_lpi_eq_sdopvar_endtoend.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_lpi_eq_sdopvar_endtoend.m` | MATTHEW M PEET, 2026-08-29, 725e5f23 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_lpi_soseq.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_lpi_soseq.m` | MATTHEW M PEET, 2026-10-04, f08854c5 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_lpivar_cdopvar.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_lpivar_cdopvar.m` | MATTHEW M PEET, 2026-09-25, fbfebe30 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_lpivar_cdopvar_stability.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_lpivar_cdopvar_stability.m` | MATTHEW M PEET, 2026-09-25, fbfebe30 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_mincap_dual.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_mincap_dual.m` | MATTHEW M PEET, 2026-10-08, 96ef81f6 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_opvar2copvar_rect.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_opvar2copvar_rect.m` | MATTHEW M PEET, 2026-09-25, 295085df |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_opvar2d2copvar_zero.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_opvar2d2copvar_zero.m` | MATTHEW M PEET, 2026-09-26, bf1fd5f3 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_poscopvar.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_poscopvar.m` | MATTHEW M PEET, 2026-09-22, a252e720 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_poscopvar_direct.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_poscopvar_direct.m` | MATTHEW M PEET, 2026-10-08, bc3be32f |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_poscopvar_lift.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_poscopvar_lift.m` | MATTHEW M PEET, 2026-10-08, bc3be32f |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_poscopvar_stability.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_poscopvar_stability.m` | MATTHEW M PEET, 2026-09-22, a252e720 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_poscopvar_vs_poslpivar.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_poscopvar_vs_poslpivar.m` | MATTHEW M PEET, 2026-09-22, a252e720 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_poscopvar_vs_poslpivar_2d.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_poscopvar_vs_poslpivar_2d.m` | MATTHEW M PEET, 2026-09-22, a252e720 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_possopvar.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_possopvar.m` | MATTHEW M PEET, 2026-08-23, e928b130 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_possopvar_integration.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_possopvar_integration.m` | MATTHEW M PEET, 2026-08-23, e928b130 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_possopvar_vs_poslpivar.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_possopvar_vs_poslpivar.m` | MATTHEW M PEET, 2026-08-23, e928b130 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_sdopvar2dopvar.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_sdopvar2dopvar.m` | Declan, 2026-10-05, 6d39aee9 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_sdopvar_tier2.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_sdopvar_tier2.m` | MATTHEW M PEET, 2026-09-29, c1ea2161 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_sdvar2dpvar.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_sdvar2dpvar.m` | MATTHEW M PEET, 2026-09-26, bf1fd5f3 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_sopquadvar_faces.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_sopquadvar_faces.m` | MATTHEW M PEET, 2026-09-28, e6c78800 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_sopvar2nopvar_guard.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_sopvar2nopvar_guard.m` | MATTHEW M PEET, 2026-09-29, a8ca3dc0 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_sopvar_mtimes_ctranspose.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_sopvar_mtimes_ctranspose.m` | MATTHEW M PEET, 2026-08-29, 725e5f23 |
| `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_space_parser.m` | `definitely/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_space_parser.m` | MATTHEW M PEET, 2026-09-30, 3e61d8dc |
| `PIETOOLS/ndopvar/Testfolder/ndopvar/claude_tests/Test_mtimes_nopvar_regression.m` | `definitely/PIETOOLS/ndopvar/Testfolder/ndopvar/claude_tests/Test_mtimes_nopvar_regression.m` | MATTHEW M PEET, 2026-08-30, 1b1be818 |
| `PIETOOLS/ndopvar/Testfolder/ndopvar/claude_tests/apply_nopvar_poly.m` | `definitely/PIETOOLS/ndopvar/Testfolder/ndopvar/claude_tests/apply_nopvar_poly.m` | MATTHEW M PEET, 2026-08-30, 1b1be818 |
| `PIETOOLS/ndopvar/Testfolder/ndopvar/claude_tests/apply_nopvar_quad.m` | `definitely/PIETOOLS/ndopvar/Testfolder/ndopvar/claude_tests/apply_nopvar_quad.m` | MATTHEW M PEET, 2026-08-30, 1b1be818 |
| `PIETOOLS/ndopvar/Testfolder/ndopvar/claude_tests/test_lpi_eq_ndopvar.m` | `definitely/PIETOOLS/ndopvar/Testfolder/ndopvar/claude_tests/test_lpi_eq_ndopvar.m` | MATTHEW M PEET, 2026-08-29, 725e5f23 |
| `PIETOOLS/ndopvar/Testfolder/ndopvar/claude_tests/test_lpi_eq_ndopvar_endtoend.m` | `definitely/PIETOOLS/ndopvar/Testfolder/ndopvar/claude_tests/test_lpi_eq_ndopvar_endtoend.m` | MATTHEW M PEET, 2026-08-29, 725e5f23 |
| `PIETOOLS/sopvar/Testfolder/converters/claude/pi_ndopvar_kernels.m` | `definitely/PIETOOLS/sopvar/Testfolder/converters/claude/pi_ndopvar_kernels.m` | MATTHEW M PEET, 2026-08-30, 1b1be818 |
| `PIETOOLS/sopvar/Testfolder/converters/claude/test_ndopvar_sdopvar_converters.m` | `definitely/PIETOOLS/sopvar/Testfolder/converters/claude/test_ndopvar_sdopvar_converters.m` | MATTHEW M PEET, 2026-08-30, 1b1be818 |
| `PIETOOLS/PIESOS/Examples/LocalStability/Testfolder_Claude/Test_LocalStability_v2.m` | `definitely/PIETOOLS/PIESOS/Examples/LocalStability/Testfolder_Claude/Test_LocalStability_v2.m` | Declan, 2026-10-15, 6fd03c15 |
| `PIETOOLS/PIESOS/Examples/LocalStability/Testfolder_Claude/Test_SOS_DP_v2.m` | `definitely/PIETOOLS/PIESOS/Examples/LocalStability/Testfolder_Claude/Test_SOS_DP_v2.m` | Declan, 2026-10-15, 6fd03c15 |
| `PIETOOLS/PIESOS/Examples/LocalStability/Testfolder_Claude/Test_Weighted_Sobolev_Norm.m` | `definitely/PIETOOLS/PIESOS/Examples/LocalStability/Testfolder_Claude/Test_Weighted_Sobolev_Norm.m` | Declan, 2026-10-15, 6fd03c15 |
| `PIETOOLS/PIESOS/Examples/LocalStability/Testfolder_Claude/Test_polyopvar_equal.m` | `definitely/PIETOOLS/PIESOS/Examples/LocalStability/Testfolder_Claude/Test_polyopvar_equal.m` | Declan, 2026-10-15, 6fd03c15 |
| `PIETOOLS/ndopvar/misc/ndopvar/claude/fill_empty_C.m` | `definitely/PIETOOLS/ndopvar/misc/ndopvar/claude/fill_empty_C.m` | MATTHEW M PEET, 2026-08-30, 1b1be818 |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/README.txt` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/README.txt` | MATTHEW M PEET, 2026-09-11, bb526674 |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/SHA256SUMS-bundle.txt` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/SHA256SUMS-bundle.txt` | MATTHEW M PEET, 2026-09-11, bb526674 |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/filelist.txt` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/filelist.txt` | MATTHEW M PEET, 2026-09-11, bb526674 |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/files_before_e3a29ef0/NEW_FILES_NO_BASELINE.txt` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/files_before_e3a29ef0/NEW_FILES_NO_BASELINE.txt` | MATTHEW M PEET, 2026-09-11, bb526674 |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/@mopvar/mopvar.m.bak` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/@mopvar/mopvar.m.bak` | Declan, 2026-01-10, bb643b07 |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/@sdopvar/mtimes.m.bak` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/@sdopvar/mtimes.m.bak` | DanBraghini, 2026-08-17, 04a486b8 |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/@sdopvar/mtimes_AT.m.bak` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/@sdopvar/mtimes_AT.m.bak` | DanBraghini, 2026-08-17, 04a486b8 |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/@sdopvar/sdopvar.m.bak` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/@sdopvar/sdopvar.m.bak` | Sachin Shivakumar, 2026-08-01, 0acdf570 |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/@sopvar/mtimes.m.bak` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/@sopvar/mtimes.m.bak` | Sachin Shivakumar, 2026-01-07, 5a655a1a |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/@sopvar/private/leftShiftMonomials_SS.m.bak` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/@sopvar/private/leftShiftMonomials_SS.m.bak` | Sachin Shivakumar, 2026-04-15, 47f680dc |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/@sopvar/sopvar.m.bak` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/@sopvar/sopvar.m.bak` | Sachin Shivakumar, 2026-01-19, 1849fb19 |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/@sopvar/sopvar2opvar2d.m.bak` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/@sopvar/sopvar2opvar2d.m.bak` | Declan, 2026-05-27, efd66d8f |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_canonical_multiplier.m.bak` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_canonical_multiplier.m.bak` | MATTHEW M PEET, 2026-08-30, 1b1be818 |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/Testfolder/sdopvar/rand_sopvar.m.bak` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/Testfolder/sdopvar/rand_sopvar.m.bak` | Declan, 2026-06-07, 0ca269c9 |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/int_semisep.m.bak` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/files_before_e3a29ef0/PIETOOLS/sopvar/int_semisep.m.bak` | Talitsky, 2026-08-17, 08e23e4e |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/files_before_e3a29ef0/SHA256SUMS.txt` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/files_before_e3a29ef0/SHA256SUMS.txt` | MATTHEW M PEET, 2026-09-11, bb526674 |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/ndopvar-unpushed.bundle` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/ndopvar-unpushed.bundle` | MATTHEW M PEET, 2026-09-11, bb526674 |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/uncommitted-worktree.patch` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/uncommitted-worktree.patch` | MATTHEW M PEET, 2026-09-11, bb526674 |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/unpushed-commits.log` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/unpushed-commits.log` | MATTHEW M PEET, 2026-09-11, bb526674 |
| `PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/unpushed-commits.patch` | `definitely/PIETOOLS/sopvar/private_files/backup_unpushed_2026-09-11/unpushed-commits.patch` | MATTHEW M PEET, 2026-09-11, bb526674 |
| `PIETOOLS/sopvar/session_summary_2026_10_08.md` | `definitely/PIETOOLS/sopvar/session_summary_2026_10_08.md` | MATTHEW M PEET, 2026-10-08, 96ef81f6 |
| `PIETOOLS/opvar/inverse_dependency_map_2026_10_09.md` | `definitely/PIETOOLS/opvar/inverse_dependency_map_2026_10_09.md` | MATTHEW M PEET, 2026-10-09, b5912704 |
| `PIETOOLS/lpi_programming_sopvar/NOTES.md` | `definitely/PIETOOLS/lpi_programming_sopvar/NOTES.md` | MATTHEW M PEET, 2026-09-29, 0a66f793 |
| `PIETOOLS/sopvar/@sdopvar/to-do-list.txt` | `definitely/PIETOOLS/sopvar/sdopvar_class/to-do-list.txt` | Sachin Shivakumar, 2026-08-10, 5234ed5b |
| `CLAUDE.md` | `definitely/CLAUDE.md` | MATTHEW M PEET, 2026-08-31, e7ed3f1a |
| `PIETOOLS/GUI/getPDEstructure_old.m` | `dead_code/PIETOOLS/GUI/getPDEstructure_old.m` | Declan, 2021-12-31, e15f2543 |
| `PIETOOLS/GUI/parsePDEParams_old.m` | `dead_code/PIETOOLS/GUI/parsePDEParams_old.m` | Declan, 2021-12-31, e15f2543 |
| `PIETOOLS/opvar/inv_opvar_new.m` | `dead_code/PIETOOLS/opvar/inv_opvar_new.m` | Declan, 2021-12-31, 4d994ea2 |
| `PIETOOLS/opvar/@opvar/subsref_obsolete.m` | `dead_code/PIETOOLS/opvar/opvar_class/subsref_obsolete.m` | Declan, 2021-12-31, 4d994ea2 |
| `PIETOOLS/ndopvar/@quadPoly/mtimes_old.m` | `dead_code/PIETOOLS/ndopvar/quadPoly_class/mtimes_old.m` | Sachin Shivakumar, 2026-01-16, 5bba3db0 |
| `PIETOOLS/ndopvar/@quadPoly/subs_OLD.m` | `dead_code/PIETOOLS/ndopvar/quadPoly_class/subs_OLD.m` | Sachin Shivakumar, 2026-01-16, 5bba3db0 |
| `PIETOOLS/ndopvar/@quadPoly/incompleteFunctions/clean.m` | `dead_code/PIETOOLS/ndopvar/quadPoly_class/incompleteFunctions/clean.m` | Sachin Shivakumar, 2026-02-11, 0992dfd0 |
| `PIETOOLS/ndopvar/@quadPoly/incompleteFunctions/todo.txt` | `dead_code/PIETOOLS/ndopvar/quadPoly_class/incompleteFunctions/todo.txt` | Sachin Shivakumar, 2026-02-15, e67d118e |
| `PIETOOLS/sopvar/@sopvar/retired/README.md` | `dead_code/PIETOOLS/sopvar/sopvar_class/retired/README.md` | MATTHEW M PEET, 2026-08-31, e7ed3f1a |
| `PIETOOLS/sopvar/@sopvar/retired/int_semisep_AT_v1.m` | `dead_code/PIETOOLS/sopvar/sopvar_class/retired/int_semisep_AT_v1.m` | Talitsky, 2026-08-17, 08e23e4e |
| `PIETOOLS/sopvar/@sopvar/retired/int_semisep_AT_v2_wrapper.m` | `dead_code/PIETOOLS/sopvar/sopvar_class/retired/int_semisep_AT_v2_wrapper.m` | MATTHEW M PEET, 2026-08-31, a8c4f802 |
| `PIETOOLS/sopvar/@sopvar/retired/int_semisep_local_mtimes.m` | `dead_code/PIETOOLS/sopvar/sopvar_class/retired/int_semisep_local_mtimes.m` | MATTHEW M PEET, 2026-08-31, a8c4f802 |
| `PIETOOLS/sopvar/@sopvar/retired/mtimes_AT.m` | `dead_code/PIETOOLS/sopvar/sopvar_class/retired/mtimes_AT.m` | Sachin Shivakumar, 2026-01-07, 5a655a1a |
| `PIETOOLS/sopvar/@sopvar/retired/rearrangeCoef_old.m` | `dead_code/PIETOOLS/sopvar/sopvar_class/retired/rearrangeCoef_old.m` | MATTHEW M PEET, 2026-08-31, a8c4f802 |
| `PIETOOLS/PIESOS/polyopvar_old.pdf` | `dead_code/PIETOOLS/PIESOS/polyopvar_old.pdf` | MATTHEW M PEET, 2026-09-06, ebc302c3 |
| `PIETOOLS/sopvar/private/dead_code/README.md` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/README.md` | MATTHEW M PEET, 2026-09-29, b3808f0f |
| `PIETOOLS/sopvar/private/dead_code/Testfolder/example_of_leftshiftmonomials_AT.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/Testfolder/example_of_leftshiftmonomials_AT.m` | Talitsky, 2026-04-01, 8223bfd7 |
| `PIETOOLS/sopvar/private/dead_code/Testfolder/leftshiftMonomoials_ATv2.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/Testfolder/leftshiftMonomoials_ATv2.m` | Talitsky, 2026-03-17, 56f9b636 |
| `PIETOOLS/sopvar/private/dead_code/Testfolder/rightshiftMonomials_ATv2.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/Testfolder/rightshiftMonomials_ATv2.m` | Talitsky, 2026-04-14, 80a4c164 |
| `PIETOOLS/sopvar/private/dead_code/Testfolder/test_leftshift_monomials_sopvar.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/Testfolder/test_leftshift_monomials_sopvar.m` | Talitsky, 2026-07-31, b92f04ba |
| `PIETOOLS/sopvar/private/dead_code/sdopvar_class/lrmultiply_batch.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sdopvar_class/lrmultiply_batch.m` | Sachin Shivakumar, 2026-08-17, a758cc6d |
| `PIETOOLS/sopvar/private/dead_code/sdopvar_class/minus_decparam.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sdopvar_class/minus_decparam.m` | Sachin Shivakumar, 2026-08-10, bc779a38 |
| `PIETOOLS/sopvar/private/dead_code/sdopvar_class/mtimes_AT.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sdopvar_class/mtimes_AT.m` | DanBraghini, 2026-08-17, 04a486b8 |
| `PIETOOLS/sopvar/private/dead_code/sdopvar_class/plus_decparam.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sdopvar_class/plus_decparam.m` | Sachin Shivakumar, 2026-08-10, bc779a38 |
| `PIETOOLS/sopvar/private/dead_code/sdopvar_class/plus_decparam_batch.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sdopvar_class/plus_decparam_batch.m` | Sachin Shivakumar, 2026-08-17, a758cc6d |
| `PIETOOLS/sopvar/private/dead_code/sdopvar_class/rand_sdopvar.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sdopvar_class/rand_sdopvar.m` | Declan, 2026-06-07, 0ca269c9 |
| `PIETOOLS/sopvar/private/dead_code/sdopvar_class/test_script.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sdopvar_class/test_script.m` | DanBraghini, 2026-08-10, 54689a9f |
| `PIETOOLS/sopvar/private/dead_code/sdopvar_class_private/CombineDecisionBasis.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sdopvar_class_private/CombineDecisionBasis.m` | DanBraghini, 2026-08-10, 07ffd3d1 |
| `PIETOOLS/sopvar/private/dead_code/sdopvar_class_private/MatrixMultiply.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sdopvar_class_private/MatrixMultiply.m` | DanBraghini, 2026-08-10, d377d14d |
| `PIETOOLS/sopvar/private/dead_code/sdopvar_class_private/UnionBasisMonomials.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sdopvar_class_private/UnionBasisMonomials.m` | Declan, 2026-06-08, bc13c581 |
| `PIETOOLS/sopvar/private/dead_code/sopvar_class/ctranspose_old.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sopvar_class/ctranspose_old.m` | Sachin Shivakumar, 2026-01-19, 1849fb19 |
| `PIETOOLS/sopvar/private/dead_code/sopvar_class/mtimestest.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sopvar_class/mtimestest.m` | Sachin Shivakumar, 2026-05-10, ecb5a76f |
| `PIETOOLS/sopvar/private/dead_code/sopvar_class/randOpvar.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sopvar_class/randOpvar.m` | Sachin Shivakumar, 2026-01-20, 3fedd01a |
| `PIETOOLS/sopvar/private/dead_code/sopvar_class/sopvar_old.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sopvar_class/sopvar_old.m` | Sachin Shivakumar, 2026-01-19, 1849fb19 |
| `PIETOOLS/sopvar/private/dead_code/sopvar_class/test_script.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sopvar_class/test_script.m` | Sachin Shivakumar, 2026-01-20, 3fedd01a |
| `PIETOOLS/sopvar/private/dead_code/sopvar_class_private/int_2b.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sopvar_class_private/int_2b.m` | Sachin Shivakumar, 2026-03-18, 4da55ecb |
| `PIETOOLS/sopvar/private/dead_code/sopvar_class_private/leftPermuteVec.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sopvar_class_private/leftPermuteVec.m` | Sachin Shivakumar, 2026-02-11, 0992dfd0 |
| `PIETOOLS/sopvar/private/dead_code/sopvar_class_private/leftShiftMonomials_old.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sopvar_class_private/leftShiftMonomials_old.m` | Sachin Shivakumar, 2026-04-15, 47f680dc |
| `PIETOOLS/sopvar/private/dead_code/sopvar_class_private/leftshiftMonomoials_AT.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sopvar_class_private/leftshiftMonomoials_AT.m` | Talitsky, 2026-03-17, 56f9b636 |
| `PIETOOLS/sopvar/private/dead_code/sopvar_class_private/mapAlphaBetaToGamma.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sopvar_class_private/mapAlphaBetaToGamma.m` | Sachin Shivakumar, 2026-02-11, 11348479 |
| `PIETOOLS/sopvar/private/dead_code/sopvar_class_private/monomial_outer.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sopvar_class_private/monomial_outer.m` | Sachin Shivakumar, 2026-01-29, 23edd836 |
| `PIETOOLS/sopvar/private/dead_code/sopvar_class_private/monomial_shift_left.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sopvar_class_private/monomial_shift_left.m` | Sachin Shivakumar, 2026-06-26, 8a95a2f7 |
| `PIETOOLS/sopvar/private/dead_code/sopvar_class_private/monomial_shift_right.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sopvar_class_private/monomial_shift_right.m` | Sachin Shivakumar, 2026-06-26, 8a95a2f7 |
| `PIETOOLS/sopvar/private/dead_code/sopvar_class_private/rightshiftMonomials_AT.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sopvar_class_private/rightshiftMonomials_AT.m` | Talitsky, 2026-03-17, 56f9b636 |
| `PIETOOLS/sopvar/private/dead_code/sopvar_class_private/termCompose.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sopvar_class_private/termCompose.m` | Sachin Shivakumar, 2026-02-11, 11348479 |
| `PIETOOLS/sopvar/private/dead_code/sopvar_class_private/termCompose_old.m` | `dead_code/PIETOOLS/sopvar/private_files/dead_code/sopvar_class_private/termCompose_old.m` | Sachin Shivakumar, 2026-02-18, dc7378f0 |
