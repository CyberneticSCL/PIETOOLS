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
  a backup of unpushed commits, session notes, to-do lists and `CLAUDE.md`; since 11/07/2026 also
  the one-off study scripts of the cuADMM harness (`PIETOOLS_demos/cuadmm/experiments/`) and the
  helpers only they called.
- `private/dead_code/` – obsolete code with no caller anywhere in the tree on 10/23/2026: a word search for
  each function name over every tracked file, comments excluded, plus the files of
  `PIETOOLS_PDE_GUI.mlapp`. It includes the earlier `sopvar/private/dead_code/` (09/29/2026), whose
  own README lists the original paths of its 31 files.
- `private/lpi_programming/`, `private/opvar/`, `private/executives/` (11/08/2026, maintainer: "move
  lpi_programming, opvar, executives, and settings folders into cleaning/private") – the stock LPI
  layer, the opvar classes with their helpers, and the stock executives, each moved whole from
  `PIETOOLS/<folder>/`; `@opvar`, `@dopvar`, `2D/@opvar2d`, `2D/@dopvar2d` became `*_class`, and the
  empty `@quadPoly`, `Testing` and `legacy` folders were dropped. `opvar/2D/int_simple.m` moved to
  `misc/` first. The 13 files of `settings/`, which code on the path still calls, moved to
  `orphaned_file/needed/` instead (maintainer: "anything that is needed can move into
  orphaned_files/needed"), which lists their callers. Four more went there at first and left later
  on 11/08/2026: the stock `lpi_ineq` and `lpigetsol` rejoined `private/lpi_programming/` once three
  container tests and `lpi_ineq_sop` stopped calling them (maintainer: "find a way to get rid of
  it"); `closedLoopPIE` moved to `executives_sopvar/utility_functions/`, and `get_eq_opts_2D`, which
  nothing on the path calls, to `private/executives/utility_functions/` (maintainer: "now address
  closedLoopPIE and get_eq_opts_2D").
- `private/settings/` (11/08/2026, maintainer: "Move it now") – the 13 files of `settings/`
  (`lpisettings`, the 1-D presets and the `2D/` presets), from `orphaned_file/needed/settings/`
  once every caller on the path had switched to `lpisettings_sop`; the tiers of `settings_sopvar`
  set the values of the presets they use themselves (copies in `settings_sopvar/private/storage_*.m`
  on 11/08/2026, merged into the tier files 11/09/2026), and an `extreme` tier stands in for the
  stock extreme level. The only remaining call is in
  `PIETOOLS_demos/lowrank_2d_stability/tests_1d/`, which already calls retired functions.

- `ndopvar/` (11/08/2026, maintainer: "move the ndopvar folder into the orphan_code directory") â€“
  the whole folder moved to `orphaned_file/private/ndopvar/` (names kept; `private/` keeps it off the
  path). Measured first with `ndopvar/` off the path: 65 demos and examples gave no new failure.
  Taken out of it beforehand: `rand_poly.m` and `rand_dpvar.m` to `PIETOOLS/misc/` (about 100 tests
  call them); to `private/dead_code/PIETOOLS/ndopvar/`: `Test_quadPoly_common_vars.m`,
  `Test_quadPoly_mtimes.m`, `Test_quadPoly_set_vars.m`, `coeff2poly.m` and `ndopvar2tensopvar.m`
  (untracked since 11/08/2026; committed at `PIESOS/misc/` in `HEAD`); and the six converters that
  only create or read an `nopvar`/`ndopvar` (no caller outside tests): `nopvar2sopvar`,
  `nopvar2copvar`, `ndopvar2sdopvar` (`converters/operator_converters/`), `sopvar2nopvar`,
  `copvar2nopvar`, `sdopvar2ndopvar` (`sopvar/{sopvar,copvar,sdopvar}_class/`). The tests that still
  call them (`Test_ndopvar_to_sdopvar_conversion`, section 6 of `Test_copvar_converters`,
  `Bench_copvar_vs_legacy`) already failed on `opvar2d`/`rand_opvar`, and
  `copvar2opvar`/`copvar2opvar2d` still name `copvar2nopvar` in an error message.

- `mopvar` and six files that only retired code calls (11/09/2026, maintainer: "move this into
  orphaned code", for the `mopvar` constructor, and "move to orphaned code") – to
  `orphaned_file/private/`: the class folder `sopvar/@mopvar/` (`mopvar`, `canonicalize`, `verify`)
  as `@mopvar/`; `sopvar/@copvar/canonicalize.m` as `@copvar/canonicalize.m`; `sdopvar2dopvar`
  (`converters/operator_converters/`); `possopvar` and `poscopvar_lift`
  (`lpi_programming_sopvar/misc/`), then `sopquadvar` and `posmult_cdopvar` from the same folder,
  whose last live callers were those two. A live caller is code outside `testing_codes/`,
  `cleaning/private/` and `orphaned_file/private/`; none had one. `SOSTOOLS400/mpmonomials.m`, also
  called only by retired code, stays: SOSTOOLS is a separate toolbox (maintainer: leave it).

- `poly_decomp_1D`, `initialize_PIETOOLS_PIE` and seven demo tools that nothing calls (11/09/2026,
  maintainer: "move to orphaned code") – to `orphaned_file/private/`: `misc/poly_decomp_1D`,
  `converters/initialize_PIETOOLS_PIE`, `cuadmm_linear_psatz` (`PIETOOLS_demos/cuadmm/`),
  `pielr_face_encode`, `regen_face`, `set2d_psatz4` (`PIETOOLS_demos/lowrank_2d_stability/`),
  `jp_wave_plate_rates` (`PIETOOLS_demos/multivariate_PIE_2508_14840/`), `heatNd_ladder` and
  `heatNd_tailor` (`PIETOOLS_demos/sopvar_demos/`). None had a caller, tests and retired code
  included. `pielr_ipm_ref`, uncalled by name too, stays: the low-rank benches reach it through
  `pielr_private(name)`.

- Tests that need the retired `opvar2d` (11/09/2026, maintainer: "move to cleaning") – from
  `testing_codes/private/` to `private/dead_code/PIETOOLS/testing_codes/private_files/`:
  `sopvar/Testfolder/sdopvar/Test_apply_sopvar_2d`, `Test_ctranspose_sopvar_2d` and
  `Test_mtimes_sopvar_2d`, which stop at `rand_opvar2d`, with `sopvar2opvar2d`, which only the last
  two call; `lpi_programming_sopvar/tests/e2e_stability_2d` (E3), which stops in its reference
  `cx_stability_2D`, with `converters/operator_converters/opvar2d2copvar`; and `cx_stability_2D`
  itself, to `private/definitely/PIETOOLS/testing_codes/private_files/sopvar/Testfolder/sdopvar/
  claude_tests/cx_exec/`. Measured on 11/09/2026. `test_endtoend_sop('all')` now runs E1, E2 and
  E4. The Phase 0 and class-ops baselines, which name `opvar2d` too, stay: they run on exports of
  older trees.

- A test that needs the retired `opvar` (11/13/2026, maintainer: "Yes, apply both", the move and the
  update of `Test_sopvar_diff`, after a triage of the tests failing on 019b383e) – from
  `testing_codes/private/` to `private/dead_code/PIETOOLS/testing_codes/private_files/`:
  `sopvar/Testfolder/sdopvar/Test_apply_sopvar_1d`, which stops at `rand_opvar` (line 29); its
  subject was `apply_sopvar` against `apply_opvar` on a random `opvar`. `Test_int2sopvar` checks
  `apply_sopvar` against `int` of the kernel in 1 to 3 variables, but stops in its first
  configuration until the `int2sopvar` defect of 11/08/2026 (791da3fb) is fixed. Measured on an
  export of 019b383e, 11/12/2026. The last section of `Test_sopvar_diff`, the one other test that
  used `opvar`, was rewritten in place on 11/13/2026.

- `copvar2opvar`, `copvar2opvar2d`, `sopvar2opvar` and 7 failing tests (11/08/2026, maintainer: "Move
  copvar2opvar, copvar2opvar2d and sopvar2opvar, as well as the failing tests into the clean folder
  for deletion") â€“ to `private/dead_code/PIETOOLS/sopvar/{copvar,sopvar}_class/`. All three build an
  `opvar`/`opvar2d`, which is no longer on the path, so each already failed when called. Tests, all
  failing on 11/08/2026 and all calling one of the converters retired above: `Test_copvar_converters`,
  `Bench_copvar_vs_legacy`, `Test_copvar_inv`, `converters/Test_ndopvar_to_sdopvar_conversion`,
  `sdopvar/Test_mtimes_sopvar_1d`, `sdopvar/test_sopvar_AT_bug` (all `sopvar/Testfolder/`), and
  `PIESOS/Examples/LocalStability/Tests/WSB_test`. Callers left behind, none reachable without an
  opvar class: `pie2opvar` (lines 198, 200; no live caller), `lift_operator` in
  `converters_sopvar/tmap_separable_sopvar.m` (line 271), and the handles `tl.toop` in
  `exec_tools_sop` (lines 193, 204; every executive call is commented out).

- `converters_sopvar/` and `pie2opvar` (11/08/2026, maintainer: "move the old converter
  `convert_PIETOOLS_PDE_sopvar` and associated functions to the clean folder ... Do the PIESIM codes
  work without pie2opvar? If so, move that function as well") â€“ the 12 files of `converters_sopvar/`
  (`convert_PIETOOLS_PDE_sopvar`, `tmap_sopvar`, `tmap_separable_sopvar`, `periodic_map_sopvar`,
  `bc_operators_sopvar`, `pde_operators_sopvar`, `ftc_expansion_sopvar` and the five of `private/`,
  now `private_files/`) to `private/dead_code/PIETOOLS/converters_sopvar/`, and
  `converters/operator_converters/pie2opvar.m` to `private/dead_code/PIETOOLS/converters/operator_converters/`.
  No caller outside the folder on the path (word search, comments excluded, 11/08/2026); nothing the
  folder calls is left without another caller. `PIESIM` converts with `pie2piesim`
  (11/04/2026), its three `pie2opvar` lines are commented out. Remaining mentions: `p0_path.m` and
  `bench_lpi_degree.m` (legacy series only) in `baselines/`, and the comments of `pie_witness_sop`.

- `piess_opvar` (11/08/2026, maintainer: "Yes, move it") – from `cleaning/` itself, where c821eb72
  put it on 11/06/2026 (from `converters/`), to `private/dead_code/PIETOOLS/converters/`. In
  `cleaning/` it was on the MATLAB path, since `genpath` skips only `private/`. It builds `opvar`
  and `opvar2d` objects, so it failed when called once `opvar/` left the path ("Unrecognized
  function or variable 'opvar'", 11/08/2026). Its one caller, lines 316–323 of `converters/piess.m`,
  is reached only by an opvar-class input (see the dead branches below).

- `getController_direct_sop` (11/08/2026, maintainer: "move the unused files to cleaning and
  commit") – from `executives_sopvar/utility_functions/` to
  `private/dead_code/PIETOOLS/executives_sopvar/utility_functions/`. The controller gains without
  forming the inverse (10/09/2026). No code caller outside `cleaning/` (word search over the tree,
  comments excluded, 11/08/2026); its one caller, `Test_copvar_inv`, is in `dead_code/` since
  11/08/2026 (above). It builds the 'opvar' form of the gain with the `opvar` constructor, off the
  path since 11/08/2026, so `opvar Kop;` (line 180) fails on the current path. Its helpers
  `gk_grid` and `chebfit_mono` (`sopvar/misc/`) stay: `sopvar/@sopvar/inv.m` calls both. Remaining
  mentions: the comments of `@sopvar/inv.m` and `gk_grid.m`, and
  `baselines/2026-11-02_opvar_phase0/opvar_removal_plan.md`.

- `PIESOS/Examples/` (11/09/2026, maintainer: "Move all the test functions and the example scripts from
  PIESOS/Examples (linear and nonlinear) to the cleaning folder for deletion") â€“ the 42 example
  scripts (`Linear_Examples/`, `LocalStability/` with `Examples/`, `Examples/Sim/`, `Interpret_LF/`,
  and the 10 in the root of `Examples/`) and the 13 tests that came from it, from
  `testing_codes/private/PIESOS/Examples/` (`Test_Linear_Examples*.m`, `examples_PDE_library_PIETOOLS_nouinput.m`,
  `LocalStability/Testfolder/`, `LocalStability/Tests/`), to `private/dead_code/PIETOOLS/PIESOS/Examples/`, paths
  kept. No file outside `PIESOS/Examples/` and `testing_codes/` names any of them (word search,
  comments excluded; no hard-coded path to the folder). `Weighted_Sobolev_Norm.m` no longer shares
  its name, up to case, with `PIESOS/executives/weighted_sobolev_norm.m`. The other PIESOS tests
  (`testing_codes/private/PIESOS/Testfolder/`, `sosprogramming/Testfolder/`) stay in `testing_codes/`.
  Six of the 42 came back the same day (maintainer: "Those are not really example files, and could
  potentially be useful"): `LocalStability/polys2sopvar.m` to `PIESOS/converters/`,
  `Interpret_LF/polyopvar_equal.m`, `Interpret_LF/filter_Pcell.m` and `DMB.m` to `PIESOS/custom/`, and
  `SOS_DP_v2.m` and `kernel_degree2.m` to `PIESOS/custom/` as `SOS_DP.m` and `kernel_degree.m`. The
  earlier `SOS_DP.m` and `kernel_degree.m` stay here. `polys2sopvar.m` returned here on 11/09/2026
  (maintainer: no caller left; only the retired `SOS_DP.m` and `V3_Liediff.m` here call it), to
  `private/dead_code/PIETOOLS/PIESOS/polyopvar/`. `Tests/SOS_DP_test.m` expects the v1 `SOS_DP`
  (`Zs` a cell; the renamed one returns the `tensopmat` `Zop`), `Tests/kernel_degree_test.m` calls
  `kernel_degree2`, `Tests/DMB_test.m` builds an `opvar()`; none passes against the live functions as it is.
- `lpi_programming_sopvar/examples/` (11/10/2026, maintainer: "move the lpi_programming/examples
  folder to the cleaning folder") – its four scripts (`hinf_gain_1d_sop`, `stability_1d_sop`,
  `stability_2d_sop`, `volterra_norm_sop`), which use the user functions of
  `lpi_programming_sopvar/`, to `private/dead_code/PIETOOLS/lpi_programming_sopvar/examples/`; the
  folder is gone. No code or test calls them (word search, comments excluded); tests and the opvar
  removal plan name them in comments. The rest of `lpi_programming_sopvar/` stays: each file has a
  caller in the library or the DEMOs, directly or through another file of the folder, except
  `misc/copquadvar` and the four `misc/private/` helpers only it calls, which stay (maintainer,
  11/09/2026: in development).

## Not moved

| Files | Reason kept |
|---|---|
| 12 helpers in `PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/` (`cx_exec/cx_PIE2PDEstability`, `cx_plant`, `cx_resid`, `cx_settings`, `cx_shape`, `cx_solve`, `cx_space_list`, `cx_stability_2D`, `cx_stability_eq_opts_2D`, `cx_stability_pos2d`; `pi_dpvar_eval`, `pi_poly_grid`) | called by 9 tests in `lpi_programming_sopvar/tests/` and by `PIETOOLS_demos/sopvar_demos/heatNd_solve.m` and `heatNd_path.m`; since 11/08/2026 ten of them are with the tests, now in `testing_codes/private/sopvar/Testfolder/sdopvar/claude_tests/`, and `cx_resid` and `cx_shape`, which the two demos call, stay (`testing_codes/README.md`); `cx_stability_2D` went to `definitely/` on 11/09/2026 (it needs the retired `opvar2d`) |
| `PIETOOLS/converters/internal/legacy/` (8 files) | called by `PIETOOLS_DDE.m`, `convert_PIETOOLS_PDE.m` and `initialize_PIETOOLS_PDE.m` |
| `PIETOOLS/opvar/2D/@opvar2d/*_legacy.m`, `@dopvar2d/*_legacy.m` (6 files) | called by the 2-D H-infinity executives; moved with `opvar/` on 11/08/2026 |
| `PIETOOLS/opvar/inv_opvar_old.m` | called by `inv_opvar.m`; moved with `opvar/` on 11/08/2026 |
| `PIETOOLS/opvar/inv_opvar_2.m` | the reference of `sopvar/Testfolder/Test_copvar_inv.m`; moved with `opvar/` on 11/08/2026 |
| `PIETOOLS/SOSTOOLS400/multipoly/s2p_old.m`, `@polynomial/p2s_old.m`, both copies of `peval.mexw32` and `peval.c` | SOSTOOLS is a separate toolbox |
| `PIETOOLS/PIETOOLS_demos/cuadmm/` (apart from `experiments/` and its helpers, moved 11/07/2026), `PIETOOLS/PIETOOLS_demos/lowrank_2d_stability/` | left alone for now (maintainer, 10/23/2026) |

Not yet decided (candidates for `maybe/`): `sopvar/sopvar_lift_notes.tex` and `.pdf`,
`sopvar/sopvar_implementation_notes.pdf`, `sopvar/sopvar_duality_pruning_notes.tex`,
`ndopvar/Notes_representation_ndopvar_v5.pdf`, and `PIETOOLS/PIETOOLS_2025_User_Manual.pdf`, which
master added on 04/29/2026 and which is an older copy of
`PIETOOLS_documentation/PIETOOLS_2025_User_Manual.pdf`.

## Code on the path that still calls a retired function (11/08/2026)

Measured on 11/07/2026 with `opvar/` and `lpi_programming/` off the path, unless marked static; the
`executives/` and `settings/` names by a word search on 11/08/2026.

- Tests whose reference or subject is an opvar class or the stock LPI layer (decision 5 of the opvar
  removal plan): 13 Phase 6 tests, 10 in `sopvar/Testfolder` and 3 in `ndopvar/Testfolder` (2 of the
  13 only through the dead line `rand_ndopvar.m:5`), and 4 in
  `PIESOS/Examples/LocalStability/Testfolder`; `Test_copvar_inv` also calls `getController`,
  `getObserver`. The legacy sections of `test_dpvar_ops_sop`, `test_lpi_eq_sop`, `test_lpi_ineq_sop`,
  `test_tier1_fixes_sop`, `test_translators_sop` and `test_getsol_sop` (iii) were deleted later on
  11/08/2026 (Phase 3 cleanup).
- `converters/internal/2D/PIETOOLS_FTC_Expansion_2D.m:225`, the 2-D nonlinear converter, builds
  opvars; so do the legacy 2-D converter chain and `convert_PIETOOLS_PDE_sopvar` (static).
- The 2022 declaration format: the `parser/@state` and `@terms` methods build opvars (static;
  decision 2), and the GUI follows it.
- `PIETOOLS_examples/examples_PDE_library_PIETOOLS.m:1060, 1069` called `PIETOOLS_auto_execute`, which
  already failed on every copvar PIE, until 11/08/2026 (96529608); they now call its container
  counterpart `PIETOOLS_auto_execute_sop`.
- `PIETOOLS_demos/multivariate_PIE_2508_14840/jp_wave_plate_rates.m` (the stock LPI layer),
  `PIETOOLS_demos/lowrank_2d_stability/` (the stock executives and LPI layer) and
  `PIETOOLS_demos/cuadmm/heat2d_compare.m` (through `lowrank_2d_stability`) (static).
  `jp_wave_plate_rates` is in `orphaned_file/private/` since 11/09/2026.
- The manual snippets of Ch3 and Ch5 build opvars (static). `lpi_programming_sopvar/examples/
  volterra_norm_sop.m` did too, until 11/08/2026.
- `sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_stability_2D.m:59, 72–74`, the reference of
  E3 (`lpi_programming_sopvar/tests/e2e_stability_2d.m`, check (0)), passes the container `PIE.T`
  and `PIE.A` to `opvar2d2copvar`, which accepts only an `opvar2d`, and builds the `eppos` margin as
  an `opvar2d` (a 2-D transcription, left alone). Measured on 11/08/2026: E3 stops there at line 59
  (`opvar2d2copvar:badInput`) with the `extreme` settings, and at line 52, its `use_sosineq = 1`
  error, with the default `light`. The transcription in E3 (`build_2d`) and part 4 of
  `lpi_programming_sopvar/tests/bench_getsol_sop.m` (`lpivar`, `getsol_lpivar`, `lpigetsol`) called
  retired functions until 11/08/2026. E3, `cx_stability_2D` and `opvar2d2copvar` are in `private/`
  since 11/09/2026 (above).
- Branches that only an opvar-class input reaches, now dead: `piesos_getsol`, `piesos_ineq`,
  `dopvar2sdopvar` (static), and lines 316–323 of `converters/piess.m`, which call `piess_opvar`
  (off the path since 11/08/2026). `lpi_ineq_sop`, `lpi_eq_sop`, `getsol_lpivar_sop`, `lpigetsol_sop`,
  `poslpivar_sop` and `get_lpivar_degs_sop` raise their own errors for such an input since
  11/08/2026.

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
| `PIETOOLS/PIESOS/Examples/Linear_Examples/Test_Conversion_Linear.m` | `definitely/PIETOOLS/PIESOS/Examples/Linear_Examples/Test_Conversion_Linear.m` | Declan, 2026-10-15, 6fd03c15 |
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
| `PIETOOLS/converters/@pie_struct_opvar/` (10 files: `ctranspose.m`, `display.m`, `eq.m`, `expand_tderivatives.m`, `horzcat.m`, `initialize.m`, `pie_struct.m`, `subsasgn.m`, `subsref.m`, `vertcat.m`) | `dead_code/PIETOOLS/converters/pie_struct_opvar_class/` | never committed: an untracked copy of the `opvar`-based `@pie_struct` as it was before the 10/29/2026 update (now `copvar`-based, with `npie_struct` for nonlinear PIEs); moved with `Move-Item`, not `git mv`; nothing in the tree refers to it by name |
| `PIETOOLS/lpi_programming_sopvar/lpiprogram_sop.m` (11/04/2026: no container routine reads the program variables its extra input forms filled; callers use `lpiprogram`) | `dead_code/PIETOOLS/lpi_programming_sopvar/lpiprogram_sop.m` | MATTHEW M PEET, 2026-09-29, 0a66f793 |
| `PIETOOLS/lpi_programming_sopvar/tests/test_lpiprogram_sop.m` (11/04/2026: the test of `lpiprogram_sop`) | `dead_code/PIETOOLS/lpi_programming_sopvar/tests/test_lpiprogram_sop.m` | MATTHEW M PEET, 2026-09-29, 0a66f793 |
| `PIETOOLS/PIETOOLS_demos/cuadmm/experiments/` (74 files: 73 scripts and `README.md`; 11/07/2026, maintainer: "move the experiments folder to cleaning"; they call the stock LPI functions and executives, which the harness no longer uses, and no file outside the folder calls them) | `definitely/PIETOOLS/PIETOOLS_demos/cuadmm/experiments/` | MATTHEW M PEET, 2026-09-25, db49553e |
| `PIETOOLS/PIETOOLS_demos/cuadmm/cuadmm_private.m` (11/07/2026, maintainer: "move them too"; called only by `experiments/`) | `definitely/PIETOOLS/PIETOOLS_demos/cuadmm/cuadmm_private.m` | MATTHEW M PEET, 2026-09-25, db49553e |
| `PIETOOLS/PIETOOLS_demos/cuadmm/cuadmm_shadow.m` (as above) | `definitely/PIETOOLS/PIETOOLS_demos/cuadmm/cuadmm_shadow.m` | MATTHEW M PEET, 2026-09-25, db49553e |
| `PIETOOLS/PIETOOLS_demos/cuadmm/shadows/lpisolve.m.txt`, `poslpivar.m.txt` (read only by `cuadmm_shadow`) | `definitely/PIETOOLS/PIETOOLS_demos/cuadmm/shadows/` | MATTHEW M PEET, 2026-09-25, db49553e |
| `PIETOOLS/PIETOOLS_demos/cuadmm/private/lpi_resid.m`, `opnorm_pi.m`, `opgram.m`, `pgridval.m`, `gaussleg.m`, `legval.m` (called only by `experiments/` and by each other) | `definitely/PIETOOLS/PIETOOLS_demos/cuadmm/private_files/` | MATTHEW M PEET, 2026-09-25, db49553e |
| `PIETOOLS/converters/piess_opvar.m` (split from `converters/piess.m` on 10/30/2026; in `cleaning/` itself, on the MATLAB path, from 11/06/2026, c821eb72; 11/08/2026, maintainer: "Yes, move it"; its one caller, `converters/piess.m`, reaches it only for an opvar-class input) | `dead_code/PIETOOLS/converters/piess_opvar.m` | Declan, 2026-10-30, e504f9ba |
| `PIETOOLS/executives_sopvar/utility_functions/getController_direct_sop.m` (11/08/2026, maintainer: "move the unused files to cleaning and commit"; no code caller outside `cleaning/`; builds an `opvar`, off the path) | `dead_code/PIETOOLS/executives_sopvar/utility_functions/getController_direct_sop.m` | MATTHEW M PEET, 2026-10-09, 368097d3 |
| `PIETOOLS/testing_codes/private/sopvar/Testfolder/sdopvar/Test_apply_sopvar_2d.m` (11/09/2026, maintainer: tests that need the retired `opvar2d`, "move to cleaning") | `dead_code/PIETOOLS/testing_codes/private_files/sopvar/Testfolder/sdopvar/Test_apply_sopvar_2d.m` | Declan, 2026-06-07, 46aed89c |
| `PIETOOLS/testing_codes/private/sopvar/Testfolder/sdopvar/Test_ctranspose_sopvar_2d.m` | `dead_code/PIETOOLS/testing_codes/private_files/sopvar/Testfolder/sdopvar/Test_ctranspose_sopvar_2d.m` | Declan, 2026-05-27, 4da69e4f |
| `PIETOOLS/testing_codes/private/sopvar/Testfolder/sdopvar/Test_mtimes_sopvar_2d.m` | `dead_code/PIETOOLS/testing_codes/private_files/sopvar/Testfolder/sdopvar/Test_mtimes_sopvar_2d.m` | Declan, 2026-05-27, 5135d8f5 |
| `PIETOOLS/testing_codes/private/sopvar/Testfolder/sdopvar/sopvar2opvar2d.m` (a `sopvar` method until 11/09/2026) | `dead_code/PIETOOLS/testing_codes/private_files/sopvar/Testfolder/sdopvar/sopvar2opvar2d.m` | Declan, 2026-05-27, efd66d8f |
| `PIETOOLS/testing_codes/private/sopvar/Testfolder/sdopvar/Test_apply_sopvar_1d.m` (11/13/2026, maintainer: "Yes, apply both"; needs the retired `opvar`) | `dead_code/PIETOOLS/testing_codes/private_files/sopvar/Testfolder/sdopvar/Test_apply_sopvar_1d.m` | Declan, 2026-06-07, 46aed89c |
| `PIETOOLS/testing_codes/private/lpi_programming_sopvar/tests/e2e_stability_2d.m` | `dead_code/PIETOOLS/testing_codes/private_files/lpi_programming_sopvar/tests/e2e_stability_2d.m` | MATTHEW M PEET, 2026-09-29, 0a66f793 |
| `PIETOOLS/testing_codes/private/converters/operator_converters/opvar2d2copvar.m` | `dead_code/PIETOOLS/testing_codes/private_files/converters/operator_converters/opvar2d2copvar.m` | MATTHEW M PEET, 2026-09-19, 2ba10689 |
| `PIETOOLS/testing_codes/private/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_stability_2D.m` | `definitely/PIETOOLS/testing_codes/private_files/sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_stability_2D.m` | MATTHEW M PEET, 2026-09-25, 2be2e6b5 |
| `PIETOOLS/lpi_programming_sopvar/examples/hinf_gain_1d_sop.m` | `dead_code/PIETOOLS/lpi_programming_sopvar/examples/hinf_gain_1d_sop.m` | MATTHEW M PEET, 2026-09-29, 0a66f793 |
| `PIETOOLS/lpi_programming_sopvar/examples/stability_1d_sop.m` | `dead_code/PIETOOLS/lpi_programming_sopvar/examples/stability_1d_sop.m` | MATTHEW M PEET, 2026-09-29, 0a66f793 |
| `PIETOOLS/lpi_programming_sopvar/examples/stability_2d_sop.m` | `dead_code/PIETOOLS/lpi_programming_sopvar/examples/stability_2d_sop.m` | MATTHEW M PEET, 2026-09-29, 0a66f793 |
| `PIETOOLS/lpi_programming_sopvar/examples/volterra_norm_sop.m` | `dead_code/PIETOOLS/lpi_programming_sopvar/examples/volterra_norm_sop.m` | MATTHEW M PEET, 2026-09-29, 0a66f793 |
