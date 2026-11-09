# testing_codes/

Test code staged for removal (maintainer, 11/08/2026: "Move all the testing codes and files which
exist just for those codes into a testing codes folder in the repo root to be eventually
cleaned"; the maintainer then placed the folder in `PIETOOLS/` rather than in the git root). Nothing
was deleted: every tracked file here is a `git mv`, and `git mv testing_codes/<path> <path>`
restores one.

## On the path

`pietools_path_update` adds `genpath` of `PIETOOLS/`, which includes this folder, so every test,
harness and helper here resolves by name as before the move. Paths mirror the old ones:
`testing_codes/<old path under PIETOOLS/>`. The one exception is `PIESOS/polyopvar/@quadPoly/test/`,
now `PIESOS/polyopvar/quadPoly_class/test/`, so that no second `@quadPoly` folder lands on the path.
Each `private/` folder keeps its name and moved with the files that call it. No name in this
folder has a second copy elsewhere on the path (word search over the tracked files, 11/08/2026).

## What moved (897 tracked files, 11/08/2026)

| From | Files | What |
|---|---|---|
| `PIESOS/Testfolder/` | 114 | the PIESOS tests |
| `PIESOS/sosprogramming/Testfolder/` | 13 | |
| `PIESOS/Examples/LocalStability/Testfolder/`, `Tests/` | 8 | |
| `PIESOS/Examples/Linear_Examples/` | 5 | `Test_Linear_Examples*.m` (4) and `examples_PDE_library_PIETOOLS_nouinput.m`, which only they call |
| `PIESOS/polyopvar/@quadPoly/test/` | 2 | now `quadPoly_class/test/` |
| `sopvar/Testfolder/` | 34 | all but the three files listed under "What stayed" |
| `lpi_programming_sopvar/tests/` | 17 | `tests/private/` included |
| `executives_sopvar/tests/` | 1 | `test_executives_sop.m` |
| `converters/tests/` | 4 | |
| `misc/rand_poly.m`, `misc/rand_dpvar.m` | 2 | random test inputs; only tests call them |
| `PIETOOLS_demos/sopvar_demos/` | 6 | `test_heatNd_pie`, `test_heatNd_lpi`, `test_heatNd_poincare`, `test_copquadvar_faces`, `heatNd_validate`, `heatNd_poincare_validate` |
| `baselines/` | 219 | the three regression and benchmark harnesses with their reference data and reports, and `.gitignore` |
| `PIETOOLS_demos/cuadmm/` | 442 | the baseline suite: 12 `bl_*` files (`bl_expect.tsv` among them), `private/` (11), `results/` (409), the `gpu/` scripts `bl_bisgpu.sh`, `bl_gpu.sh`, `bl_gpu_post.sh` and `regime.sh`, `sol/` (4), `SUITES.md`, `BASELINE_REPORT.md` |
| `PIETOOLS_demos/lowrank_2d_stability/` | 30 | `tests_1d/` (26), `pielr_bench`, `pielr_bench_cases`, `pielr_bench_diff`, `BENCH.md` |

Also moved, with `mv`: 1,042 untracked files (55 MB), the fixture data of
`baselines/2026-11-02_opvar_phase0/` (`F1/` 292, `F2/` 688, `F4/` 22, `logs/` 40). The harness reads
them from its own folder; `testing_codes/baselines/.gitignore` keeps them out of git, as the old
`baselines/.gitignore` did.

New: `PIETOOLS_demos/cuadmm/.gitattributes` and `.gitignore`, copies of those of
`PIETOOLS_demos/cuadmm/`. The rule `*.sh text eol=lf` keeps the moved shell scripts LF on a
Windows checkout, which WSL bash needs.

## How the files were chosen

Seed: every tracked file under a folder named `Testfolder*`, `Testing`, `test`, `tests`, `tests_*`
or `claude_tests`; every `.m` file named `Test_*`, `test_*`, `*_test`, `*_tests`, `e2e_*` or
`*_validate`; and the three harness groups the maintainer chose. Not searched: `SOSTOOLS400/` (a
separate toolbox), `orphaned_file/` and `cleaning/private/`. Two rules then ran until nothing
changed:

1) A helper joins the set when every caller is in the set. A caller of a name is a live code file
   (`.m`, `.sh`, `.slurm`) holding the name outside comments; a file in `D/private/` is callable
   only from `D/` and `D/private/`, and a caller with its own private copy of the name calls that
   copy. Helpers were taken only from places that hold test scaffolding (`misc/rand_*`,
   `PIESOS/Examples/Linear_Examples/`, the cuADMM `private/`, `gpu/` and `sol/`): a function whose
   only callers are its tests is the code under test, not a helper of the tests.
2) A file of the set that code outside the set calls stays, and so does what it calls.

## What stayed, and why

| File | Called by |
|---|---|
| `PIETOOLS_demos/cuadmm/bl_bisect.m` | `PIETOOLS_demos/cuadmm/heat2d_compare.m` and `PIETOOLS_demos/multivariate_PIE_2508_14840/jp_wave_plate_rates.m` |
| `PIETOOLS_demos/cuadmm/bl_face.m`, `bl_fix.m`, `gpu/cuadmm_pinf_inf.patch` | `bl_bisect` |
| `sopvar/Testfolder/sdopvar/apply_sopvar.m` | `PIESOS/@polyopvar/combine_factors`, `PIESOS/@tensopmat/apply_op`, `PIESOS/@tensopvar/apply_op`, `PIESOS/misc/tensopvar2polynomial` |
| `sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_resid.m`, `cx_shape.m` | `PIETOOLS_demos/sopvar_demos/heatNd_solve.m` and `heatNd_path.m` (`cx_shape` through `cx_resid`) |
| `PIESOS/Examples/LocalStability/Examples/Sim/BCs_test.m` | `Sim_PDE.m`: the boundary conditions of a simulation example, not a test |

Functions that only tests call, kept as code under test or user functions: `Liediff`,
`piesos_lpivar` and `piesos_poslpivar` (`PIESOS/sosprogramming/`); `kernel_degree` and
`polyopvar_equal` (`PIESOS/Examples/LocalStability/`); `subs_dvar_sop` and the `examples/` of
`lpi_programming_sopvar`; `opvar2d2copvar`, `pie2copvar` and `int2copvar`; the heatNd demo routines;
the DEMOs that `run_phase0.sh` runs; and `pielr_ipm_ref` in the low-rank `private/`, which
`pielr_bench` reaches through `pielr_private`.

## Code changed for the move

Six files built a path from their own location; each change is marked "MMP, 11/08/2026".

- `PIETOOLS_demos/lowrank_2d_stability/pielr_path.m` (stays): the three `pielr_bench*` names move from
  `own`, which requires that folder, to `need`, which requires one copy. `pielr_path` errored
  otherwise.
- `PIETOOLS_demos/lowrank_2d_stability/tests_1d/t1d_path.m`: `D.root` is four levels up and `D.pkg`
  the package in `PIETOOLS_demos/`; `D.pkg` was the parent folder.
- `PIETOOLS_demos/cuadmm/gpu/regime.sh`: MATLAB also gets `PIETOOLS_demos/cuadmm`, which holds
  `cuadmm_path`, the first call of `bl_regime`.
- `PIETOOLS_demos/cuadmm/sol/sol_push.sh`: `TREE` climbs three levels. With two it would have
  deleted the tree on Sol and sent only `testing_codes/`.
- `PIETOOLS_demos/cuadmm/sol/sol0.sh`: the job is submitted from
  `testing_codes/PIETOOLS_demos/cuadmm` of the pushed tree.
- `baselines/2026-11-02_opvar_phase0/run_phase0.sh`: `P0W` names the new folder.

Valid without a change: `cuadmm_path` requires one copy of `bl_run`, `bl_check` and `bl_fix` inside
the PIETOOLS tree, which holds this folder; `heatNd_path` asserts the `heatNd_*` files of its own
folder; `sol_harness.slurm` adds `PIETOOLS_demos/cuadmm` of the tree on Sol. Not changed:
`baselines/2026-11-02_opvar_phase0/private/p0_path.m` still takes the old `PIETOOLS_demos/cuadmm` and
`PIETOOLS_demos/lowrank_2d_stability` off the path of an export; in an export of this layout, their
moved parts stay on it. Their 39 path-visible names match none of the 757 other path-visible names
of the tree, the 36 harness names or the 22 names `p0_path` asserts (measured 11/08/2026 by the
harness owner, who keeps the list as it is).

READMEs updated: `executives_sopvar/`, `lpi_programming_sopvar/`, `PIETOOLS_demos/sopvar_demos/`
(run the heatNd checks after `heatNd_path`), `PIETOOLS_demos/cuadmm/`, and the "Not moved" table of
`cleaning/README.md`. Comments elsewhere that name an old path (the low-rank documentation,
`cx_shape`, converter headers, the sopvar notes) are left as written.

## Checked (11/08/2026, on the default path)

- Resolution: each of the 240 moved functions resolves first, and only once in the tree, at its
  new path, except the seven `tests_1d` copies of the low-rank `private/` core (`bm_dr`, `bm_lm2`,
  `bm_proj`, `bm_report`, `bm_resid`, `bm_setup`, `raw_data`), which the suite keeps on purpose
  and which the package reaches as private functions first. `test` and `transfer_function_test`
  (`quadPoly_class/test/`) are on the path for the first time; nothing else on it has either name.
- Run by name: the eight container LPI tests, `test_executives_sop`, `e2e_stability_1d`,
  `Test_int2sopvar` and `Test_mtimes_id_bypass` – no error. `test_convert_PIETOOLS_PDE_periodic`
  fails 1 of 18 cases (`integral_boundaries`) and `test_heatNd_pie` stops at the `copvar` field
  `var1`, both identically on an export of 5553014c, before the move.
- `t1d_path` (tree, package and private folder found), `pielr_path`, `heatNd_path` and
  `cuadmm_path` – no error; `bl_check('structure')` – 16 PASS.
- The shell scripts need the GPU, WSL or Sol and were not run: `bash -n` parses all four, and
  the folders they compute exist (`TREE` is the PIETOOLS folder, the bridge holds `cuadmm_path`).
- Not run: the Phase 0 harness, the DEMOs and the example library. No code outside this folder
  changed but `pielr_path`.
