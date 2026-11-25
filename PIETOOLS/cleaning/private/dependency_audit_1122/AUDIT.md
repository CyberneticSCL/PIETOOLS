# Toolbox references to `cleaning/`, `orphaned_file/` and `testing_codes/` – audit of 11/22/2026

Maintainer, 11/22/2026: "make sure nothing in the repo depends on or references anything in cleaning,
orphaned code, or testing folders."

**Scope.** The 1,255 tracked files outside the three folders at HEAD 27513f00 (1,140 `.m` files, the GUI
app `PIETOOLS_PDE_GUI.mlapp`, 11 READMEs, tex notes, html, `demos.xml`, git configuration).

**Methods.**
1. Every text file searched for the folder names (`cleaning`, `orphaned`, `testing_code`, `Testfolder`,
   `testing folder`); hits in `data/text_hits.txt`.
2. An index of the names defined only in the three folders: 998 of their 1,047 function, class and
   package names (`data/zonly.txt`). Every toolbox code token, comments stripped and strings kept,
   checked against it; comments checked separately (`data/code_hits.tsv`, tool `tools/deps.pl`).
3. MATLAB's dependency analysis, `matlab.codetools.requiredFilesAndProducts(f,'toponly')`, of all 1,140
   `.m` files and the GUI app. The toolbox goes first on the path and the three folders after it, so that
   a name resolves into them only where the toolbox defines no such name (`data/dep_out.tsv`, tool
   `tools/dep_run.m`, 229 s).
4. The callers of every function found in 2 and 3, by name.

The runtime counterpart is the stress test of 11/22/2026, which runs the release on an export without
`cleaning/` and `orphaned_file/` (`cleaning/private/stress_test_1122/`).

## 1. No file-level dependency

- No toolbox code names the three folders: there is no `addpath`, `fullfile`, `load`, `fopen`, `run` or `cd`
  into them. Every line that names them is a comment, a README line or a tex line (§3).
- On the release path (`pietools_path_update` = `genpath`, which skips `private/`) nothing in the three
  folders is visible: each holds only a README at its top level, and everything else lies under `private/`.
- MATLAB's analysis resolves no toolbox file into `orphaned_file/` or `testing_codes/`. Its other hits fall
  into two groups:
  - 190 of its 240 hits resolve a method name (`initialize` 145, `isvalid` 35, `op_slice` 6,
    `expand_tderivatives` 3, `clean` 1) to a retired class folder in `cleaning/`. That folder is a plain
    folder (`@X` renamed `X_class`), so it appears only when `cleaning/private` is put on the path by
    hand. They are artifacts of method 3.
  - The other 50 hits are §2.
- The GUI app names no function of the three folders. Its one converter call is
  `convert_PIETOOLS_PDE_terms`, which has run on containers since 11/05/2026. The only token that matches
  is the word `test`.

## 2. Calls to functions that exist only in `cleaning/` (retired opvar classes and LPI routines)

These lines stop with "Undefined function" if reached.

| Group | Files | Retired calls | Reached from |
|---|---|---|---|
| Low-rank 2-D certify, unported since 11/08/2026 | `PIETOOLS_demos/lowrank_2d_stability/`: `pielr_certify`, `pielr_certify_pos`, `build_stab_2d_st2`, `demo_poincare_2d`, `private/opcheck_2d`, `private/opcheck_pos2d`, `private/pielr_evalgram`, `private/pielr_evalop` | `poslpivar_2d`, `lpi_eq_2d`, `lpi_eq`, `get_eq_opts_2D`, `lpigetsol`, `clean_opvar`, `opvar2d()`, `opvar()` | `demo_lowrank_2d`, `demo_your_system`, `demo_poincare_2d` (user-facing files of that folder); `pielr_certify` and `pielr_certify_pos` test `which('poslpivar_2d')` and stop with a message |
| Legacy 2-D converter, a closed group | `converters/legacy/convert_PIETOOLS_PDE_2D`, `convert_PIETOOLS_PDE_nonlinear_2D`, `legacy/2D/Compute_Tmap_PDE_2D_Separable`, `PIETOOLS_FTC_Expansion_2D`, `compute_PIETOOLS_PDE_periodic_map_2D`, and `converters/internal/compute_PIETOOLS_PDE_periodic_map` (DB, 10/28/2026; its only caller is `convert_PIETOOLS_PDE_2D`) | `opvar()`, `opvar2d()`, `clean_opvar`, `mat2opvar` | nothing outside the group |
| Branches for opvar inputs in live functions | `converters/operator_converters/op2copvar_sop` (`opvar2d()`), `opvar2copvar`, `opvar2sopvar`, `opvar2d2sopvar`, `dopvar2sdopvar` (`opvar2dopvar`); `parser/piess` (`piess_opvar`); `PIESOS/piesos_getsol` (`getsol_lpivar`), `PIESOS/piesos_ineq` (`lpi_ineq`) | as listed | an opvar argument only, which cannot be built without the retired class |
| Class attributes | `parser/@pde_struct/pde_struct.m` (`InferiorClasses` `?opvar`, `?opvar2d`), `parser/@pie_struct/pie_struct.m` (`?opvar2d`) | – | MATLAB loads both classes without the named ones (every run since 11/08/2026) |

## 3. Text that points into the three folders (71 lines in 11 files, plus 2 false hits)

| File | Lines | Content |
|---|---|---|
| `PIETOOLS_demos/sopvar_demos/README.md` | 21 | the heatNd benchmark, whose routines and tests moved to `testing_codes/private/PIETOOLS_demos/sopvar_demos/` (only `heatNd_apply` stays; only tests call it) |
| `lpi_programming_sopvar/README.md` | 20 | history of moves to `cleaning/`, `orphaned_file/`, `testing_codes/`; the sections on `examples/` and `tests/`. Stale: `poscopvar_lift` and `posmult_cdopvar` are listed in `orphaned_file/private/`, but they returned to `misc/` in bb367fe1 (11/18/2026) |
| `PIETOOLS_demos/cuadmm/README.md` | 9 | the baseline suite in `testing_codes/`, `experiments/` and shadows in `cleaning/`, `cuadmm_linear_psatz` in `orphaned_file/` |
| `PIETOOLS_demos/lowrank_2d_stability/pielr_path.m` | 6 | header and comments: the benches in `testing_codes/` |
| `executives_sopvar/README.md` | 4 | the cx_exec transcriptions and `Test_copvar_inv` under the old `sopvar/Testfolder`; the stock helpers retired to `cleaning/`; the location of `test_executives_sop` |
| `PIETOOLS_demos/cuadmm/cuadmm_path.m`, `cuadmm_outdir.m` | 3, 3 | header comments: the suite in `testing_codes/` |
| `settings_sopvar/README.md` | 2 | the heatNd PIE in `testing_codes/`; a measurement file in `cleaning/private/settings_1116/` |
| `sopvar/sopvar_lift_notes.tex`, `sopvar_duality_pruning_notes.tex` | 1, 1 | tests under the old `sopvar/Testfolder` (now in `testing_codes/` or `cleaning/`) |
| `PIETOOLS_demos/lowrank_2d_stability/README.md` | 1 | `regen_face` in `orphaned_file/private/` |

False hits: "orphaned row" (`PIESOS/polyopvar/@polyopvar/clean_poly.m`) and "cleaning process"
(`converters/legacy/convert_PIETOOLS_PDE_terms_internal.m`).

## 4. Names of retired functions in code and comments (no dependency)

- Code: 252 tokens name the opvar classes (`opvar`, `opvar2d`, `dopvar`, `dopvar2d`). By folder:
  converters 75, parser 69, sopvar 50, PIETOOLS_demos 38 (low-rank), PIESOS 17 and PIESIM 3. They
  appear in `isa` checks, in messages, and in the opvar-style field access of the containers. Message
  strings also name retired functions (`pielr_lpi`, `settings2possopvar`, `cuadmm_path`: `cuadmm_shadow`).
- Comments: 681 lines in 113 files name a retired function or a test that now lives in the three folders
  (`opvar` 183, `opvar2d` 154, `dopvar` 42, `dopvar2d` 40, `poslpivar` 22, `poslpivar_2d` 22, `lpiscript`
  16, `lpi_ineq` 16, `lpisettings` 15, ...). The counts leave out `scalar` and `test`, which are plain
  words in those comments.

## 5. The other direction (allowed, listed for completeness)

`PIETOOLS_demos/sopvar_demos/heatNd_apply.m` lies in the toolbox, and only tests call it.

## Actions of 11/22/2026 (the first four committed before the stress run of 04:30)

| Maintainer | Done |
|---|---|
| "move all the .md files into orphaned_code/claude" | The 8 markdown files written by Claude sessions moved to `orphaned_file/claude/<path>` (list in its README). Three stay: the repository README, the SOSTOOLS README and the agent definition in `.claude/agents/`. |
| low-rank certify: "ask the low rank agent to do this" | Retired by the low-rank session: 15 files moved to `orphaned_file/private/PIETOOLS_demos/lowrank_2d_stability/`, and `pielr_path.m`, `private/restrict_solve.m` and `private/pielr_evalop.m` edited. That session verified it with `pielr_bench` on 5 cases. |
| `compute_PIETOOLS_PDE_periodic_map`: "move to orphaned code" | Moved to `orphaned_file/private/`. |
| opvar-input branches: "move these to orphaned code" | `opvar2copvar`, `opvar2sopvar`, `opvar2d2sopvar` and `dopvar2sdopvar` moved to `orphaned_file/private/`. The branches were removed from `op2copvar_sop`, `piess`, `piesos_getsol`, `piesos_ineq`, `tensopmat`, `@tensopmat/initialize` and `@tensopvar/initialize`, and the originals kept in `orphaned_file/private/opvar_branches_1122/`. 13 demos and tests give the same output on HEAD and the edited tree, apart from solver times. |
| "move the legacy 2-D converters to orphaned code" | The 7 files of the closed group moved to `orphaned_file/private/converters/legacy/` with their paths: `convert_PIETOOLS_PDE_2D`, `convert_PIETOOLS_PDE_nonlinear_2D`, `pie_pde_maps_2D`, and `2D/` `Compute_Tmap_PDE_2D_Separable`, `PIETOOLS_FTC_Expansion_2D`, `expand_PIETOOLS_PDE_BCs`, `compute_PIETOOLS_PDE_periodic_map_2D`. The last three were found by caller search: nothing else calls them. Three 2-D library examples convert the same as before. |

After these moves, no toolbox code calls a function of the three folders. The names of the retired
classes remain in `isa` checks, messages and comments (§4).

## Decisions for the maintainer

1. Low-rank 2-D certify (§2, row 1): retired 11/22/2026 (Actions).
2. Legacy 2-D converter group (§2, row 2; DJ and DB files): moved to orphaned_file/private 11/22/2026 (Actions).
3. Opvar branches and opvar names (§2 rows 3–4, §4): these belong to the opvar-removal plan. Our files are
   the converters and the container messages. `parser/` and `PIESOS/` belong to DJ, `PIESIM/` to YP.
4. Text pointers (§3), core folders. Done (11/22/2026):
   `settings_sopvar/README.md` (2 lines), `executives_sopvar/README.md` (4), `sopvar/sopvar_lift_notes.tex`
   (2, among them a `claude_tests` pointer that the folder search missed) and
   `sopvar/sopvar_duality_pruning_notes.tex` (1). Not done: `lpi_programming_sopvar/README.md` needs
   revision rather than a scrub. Its introduction and dispatch tables describe a folder beside
   `lpi_programming/` that sends opvar objects to the stock routines. Those routines are in `cleaning/`,
   and since 11/08/2026 the `_sop` routines raise their own errors for those classes. The README also
   documents `tests/`, `examples/` and `subs_dvar_sop`, which live in the three folders. The lift notes
   cite `cuadmm-hinf-objective-class`, which is not a file of the repository.
5. Text pointers (§3), research folders (`PIETOOLS_demos/cuadmm`, `lowrank_2d_stability`, `sopvar_demos`):
   their READMEs document code that lives in `testing_codes/`. The documentation could move there with
   the code, or the folders could leave the release.
6. Comments naming retired functions (§4): a comment pass, if wanted.
