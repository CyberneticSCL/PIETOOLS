# New internal helpers outside PIESOS: type, use and need

20 November 2026, HEAD 6c2cd695 (working tree). Full table: `helpers_review.tsv` (this folder).

## Scope

The 107 files classed *New, internal* outside PIESOS in `CLASSIFICATION.md`: new since the 2025
release and located under `private/`, `misc/`, `internal/`, `utility_functions/`, `PIESIM_routines/`,
`converters/legacy/` or the third-party folder. Not included: the 9 internal replacements of 2025 files
(`getController_sop`, `getObserver_sop`, `get_lpivar_degs_sop`, `poslpivar_sop`, `lpi_eq_cdopvar`,
`lpi_eq_sdopvar`, `apply_sopvar`, `degbalance_core`, `eq_copvar`) and the 31 internal files of PIESOS.

## Method

1. Code: every `.m` file with its comments removed by the stripper of the release comment pass
   (`cleaning/private/release_comment_pass/codeident.awk`).
2. References: an identifier in the code, not preceded by `.`, equal to the name of a release file.
   Names inside strings count, so `feval('name')` is seen. MATLAB visibility applies: a file in
   `X/private/` is visible from the files under `X/` only, and there it shadows a public file of the
   same name; class methods are reached by dispatch and need no reference.
3. Reached: referenced, directly or through other reached files, from a user-level release file
   (any release file outside the internal folders above). A caller counts only if it is reached.
4. Copies: similarity of the code lines (whitespace collapsed, the own name replaced, multiset
   Jaccard) between the helpers and every release file of the same name; 0.8 or more is a copy.
5. Size: code lines, comments and blank lines excluded.
6. Tests: files under `testing_codes/`; *other*: `cleaning/`, `orphaned_file/` and the demo folders
   outside the release.

Limit of the method: a name assembled at run time (`feval` of a concatenated string) is not seen.
For the 12 unreached helpers, a search of all code, strings included, found no such construction for
their names, for the legacy converter, or for `copquadvar`, `poscopvar_lift` and `posmult_cdopvar`.

## Rule for an independent file

Section 1 of the project CLAUDE.md admits a separate helper for duplicated logic across two or more
real callers, or for a block too large to sit inline. Applied to each helper, in this order:

1. **Not needed**: not reached from the release code.
2. **Merge**: a copy (similarity 0.8 or more) of the same helper in the other class; one shared file
   can serve both classes.
3. **Keep, shared**: two or more reached callers.
4. **Keep, tests**: one reached caller, and test files call the helper directly.
5. **Keep, large**: one reached caller, 100 or more code lines.
6. **Fold**: one reached caller, at most 40 code lines: a local function of the caller, or inline.
7. **Optional fold**: one reached caller, 41 to 99 code lines.

The size bounds 40 and 100 are a choice of this review; the table lists every size.

## Result

| verdict | files |
|---|---|
| Keep: shared | 57 |
| Fold into its caller | 10 |
| Merge: copy of the same helper in the other class | 8 |
| Optional fold | 7 |
| Not needed in the release: reached only by tests and benchmarks | 7 |
| Keep: one caller, large | 7 |
| Keep: tests call it | 6 |
| Not needed: reached only from the legacy 2-D converter, which no release file calls | 4 |
| Not needed: no reference anywhere | 1 |

64 of the 107 helpers stand as independent files under the rule (57 shared, 7 large); 6 more stand
because tests call them. The other 37 are not needed in the release (12), copies (8), or single-caller
files of at most 99 code lines (17).

## By type and function

| category | files | shared (kept) | files (names) |
|---|---|---|---|
| A1 index and space conventions | 10 | 9 | cell_of_gamma, check_reserved_names, copvar_space_list, gamma_of_cell, kron_split, kron_strides, monomial_position, monomial_position_map, multiindex_grid, parse_copvar_spaces |
| A2 canonical form of multipliers | 6 | 5 | canonical_adjoint_map, canonical_var_order, canonicalize_multiplier, is_canonical_multiplier, monomial_gather, prune_zero_monomials |
| A3 container assembly and metadata | 9 | 8 | merge_dvar_lists, merge_dvar_pair, cat_copvar_grid, derive_copvar_meta, mat2copvar_grid, merge_copvar_registry, on_registry_sop, split_copvar, verify_copvar_meta |
| A4 block algebra | 12 | 9 | ChangeDecVar, apply_basis_map, leftShiftMonomials_SD, leftShiftMonomials_SS, lrmultiply, sync_basis, UnionBasisMonomials, leftShiftMonomials_SS, compose_sopvar, dpvar_op_copvar, int_semisep, spIkron |
| A5 kernel evaluation and display | 6 | 2 | decode_params, display_term, params2kernel, decode_params, display_term, params2kernel |
| A6 argument parsing | 3 | 3 | parse_diff_inputs, parse_int_inputs, pvar2varname |
| A7 operator inverse | 2 | 0 | chebfit_mono, gk_grid |
| B1 positive-operator construction | 12 | 3 | copquadvar, get_lift_degs, lift_copvar, poscopvar_direct, poscopvar_lift, posmult_cdopvar, build_exponent_grid, dvar_rows, expand_full, merge_monomial_product, process_degrees_one, unpack_sheets |
| B2 equality constraints | 4 | 2 | eq_opts_sopvar, lpi_soseq, collect_eq_rows, impose_eq_rows |
| B3 solve diagnostics | 3 | 0 | lpi_classify_sop, lpi_shape_sop, lpigetdual_sop |
| C1 executive setup | 7 | 4 | exec_pie_sop, exec_sop_settings, exec_tools_sop, poslpivar_settings_2d_sop, poslpivar_settings_sop, settings2possopvar, settings_2d_sop |
| C2 LPI assembly and solve loop | 3 | 3 | exec_kmult_sop, lpi_solve_loop_sop, trace_rn_sop |
| C3 results | 2 | 2 | exec_info_sop, qp_copvar |
| C4 failure witness (1-D Galerkin) | 2 | 1 | pie_disc_sop, pie_witness_sop |
| D1 BVP converter pipeline | 9 | 2 | compose_Tmaps, compute_Tmap_1D, compute_canonical_BCs, decouple_BVPs, eliminate_BVP_sources, fund_inputs, lift_BVP_sources, split_BVPs, convert_PDE_term |
| D2 legacy 2-D conversion | 4 | 0 | compute_PIETOOLS_PDE_periodic_map, compute_PIETOOLS_PDE_periodic_map_2D, convert_PIETOOLS_PDE_nonlinear_2D, pie_pde_maps_2D |
| E1 PIE structure compatibility | 7 | 0 | pie_compat_tab, pie_compat_vars, pie_space_canon, pie_compat_tab, pie_compat_vars, pie_rename_dummy, pie_space_canon |
| E2 equation display | 3 | 3 | print_terms, use_unicode, varname2display |
| F1 PIESIM setup and discretization | 2 | 1 | PIESIM_Mat_Truncate, weightPIE |
| F2 third-party | 1 | 0 | clencurt |

- **A, operator classes (48)**: the index maps between multi-indices, parameter cells, Kronecker
  positions and monomials (A1, one definition each, used across the classes and the LPI functions);
  the canonical storage of multiplier directions and its adjoint map (A2); the block grids of
  `copvar`/`cdopvar` (A3); sums and products of `sopvar`/`sdopvar` blocks (A4); kernels and display
  terms from the stored parameters (A5); argument checks of `diff` and `int` (A6); the inverse (A7).
- **B, LPI programming (19)**: the construction behind `poscopvar` (`poscopvar_direct`, `lift_copvar`,
  `get_lift_degs` and two private helpers) and the earlier quadratic-form and lift constructions
  (`copquadvar`, `poscopvar_lift`, `posmult_cdopvar` and four private helpers) (B1); the equality path
  of `lpi_eq_sop` (B2); solve diagnostics (B3).
- **C, executives (14)**: preparation of the PIE and the settings (C1), the LPI terms and the solve
  loop (C2), the summary and the least-squares reconstruction of gains (C3), and the 1-D Galerkin
  witness of a failed certificate (C4).
- **D, converters (13)**: the stages of `convert_PIETOOLS_BVP` (D1) and the legacy 2-D conversion (D2).
- **E, parser (10)**: the compatibility tables of `pie_struct`/`npie_struct` (E1) and the equation
  display (E2).
- **F, PIESIM (3)**: two setup and discretization steps and one third-party routine.

## Not needed in the release (12)

| file | folder | reached callers | code lines | test files | added by |
|---|---|---|---|---|---|
| `clencurt.m` | PIESIM/PIESIM_routines/Third_Party_Supporting_Files | none | 13 | 0 | Yulia |
| `compute_PIETOOLS_PDE_periodic_map.m` | converters/internal | none | 108 | 0 | DanBraghini |
| `compute_PIETOOLS_PDE_periodic_map_2D.m` | converters/legacy/2D | none | 351 | 0 | DanBraghini |
| `convert_PIETOOLS_PDE_nonlinear_2D.m` | converters/legacy | none | 989 | 1 | Declan |
| `copquadvar.m` | lpi_programming_sopvar/misc | none | 425 | 2 | MATTHEW M PEET |
| `dvar_rows.m` | lpi_programming_sopvar/misc/private | none | 11 | 0 | MATTHEW M PEET |
| `expand_full.m` | lpi_programming_sopvar/misc/private | none | 13 | 0 | MATTHEW M PEET |
| `merge_monomial_product.m` | lpi_programming_sopvar/misc/private | none | 38 | 0 | MATTHEW M PEET |
| `pie_pde_maps_2D.m` | converters/legacy | none | 153 | 0 | Declan |
| `poscopvar_lift.m` | lpi_programming_sopvar/misc | none | 131 | 0 | MATTHEW M PEET |
| `posmult_cdopvar.m` | lpi_programming_sopvar/misc | none | 100 | 0 | MATTHEW M PEET |
| `unpack_sheets.m` | lpi_programming_sopvar/misc/private | none | 40 | 0 | MATTHEW M PEET |

- `clencurt`: no code refers to it, in the release, the tests or the cleanup folders.
- The quadratic-form and lift constructions: `poscopvar` calls `poscopvar_direct`, and no release
  file calls `copquadvar`, `poscopvar_lift` or `posmult_cdopvar`. The private helpers `dvar_rows`,
  `expand_full`, `merge_monomial_product` and `unpack_sheets` serve `copquadvar` only. Two test files
  call `copquadvar`; files outside the release (`cleaning/private/`, and
  `orphaned_file/private/sopquadvar.m`) call `copquadvar` (7 files), `poscopvar_lift` (1) and
  `posmult_cdopvar` (1). The private helpers `build_exponent_grid` and `process_degrees_one` serve
  `copquadvar` and `poscopvar_direct`; without `copquadvar` each has one caller. A copy of
  `copquadvar` under `testing_codes/` would no longer see the private helpers of
  `lpi_programming_sopvar/misc/`.
- The legacy 2-D conversion: the N-D converter `convert_PIETOOLS_PDE` dispatches to
  `convert_PIETOOLS_BVP`, `convert_PIETOOLS_PDE_batch` and `convert_PIETOOLS_PDE_terms_legacy`. No
  release file calls `converters/legacy/convert_PIETOOLS_PDE_2D` (the 2025 converter), and that
  converter is the only caller of the four files above, directly or through
  `convert_PIETOOLS_PDE_nonlinear_2D`. One test file calls `convert_PIETOOLS_PDE_nonlinear_2D`
  directly. The same holds for the 2025 files of `converters/legacy/2D/` that call the legacy
  converter (`Compute_Tmap_PDE_2D_Separable`); they are outside this list.

## Copies (8)

| file | folder | reached callers | code lines | test files | added by |
|---|---|---|---|---|---|
| `display_term.m` | sopvar/@sdopvar/private | get_display | 146 | 0 | Declan |
| `display_term.m` | sopvar/@sopvar/private | get_display | 142 | 0 | Declan |
| `pie_compat_tab.m` | parser/@npie_struct/private | npie_struct | 66 | 0 | Declan |
| `pie_compat_tab.m` | parser/@pie_struct/private | pie_struct | 63 | 0 | Declan |
| `pie_compat_vars.m` | parser/@npie_struct/private | npie_struct | 30 | 0 | Declan |
| `pie_compat_vars.m` | parser/@pie_struct/private | pie_struct | 30 | 0 | Declan |
| `pie_space_canon.m` | parser/@npie_struct/private | vertcat;initialize;horzcat | 16 | 0 | Declan |
| `pie_space_canon.m` | parser/@pie_struct/private | vertcat;horzcat;initialize | 16 | 0 | Declan |

`pie_compat_vars` and `pie_space_canon` are identical in `@pie_struct/private` and
`@npie_struct/private`; `pie_compat_tab` agrees on 92 % of its code lines, `display_term` on 82 %
between `@sopvar/private` and `@sdopvar/private`. One file per pair, on a shared path
(`parser/internal/`, `sopvar/misc/`), serves both classes. Two further pairs share part of their
code and keep their verdicts:

| file | overlap | verdict |
|---|---|---|
| `sopvar/@sdopvar/private/leftShiftMonomials_SS.m` | partial copy of sopvar/@sopvar/private/leftShiftMonomials_SS.m (0.65) | Keep: shared |
| `sopvar/@sopvar/private/leftShiftMonomials_SS.m` | partial copy of sopvar/@sdopvar/private/leftShiftMonomials_SS.m (0.65) | Optional fold |
| `sopvar/@sdopvar/private/decode_params.m` | partial copy of sopvar/@sopvar/private/decode_params.m (0.68) | Keep: shared |
| `sopvar/@sopvar/private/decode_params.m` | partial copy of sopvar/@sdopvar/private/decode_params.m (0.68) | Keep: shared |

## Fold into the caller (10)

| file | folder | caller | code lines | test files | added by |
|---|---|---|---|---|---|
| `build_exponent_grid.m` | lpi_programming_sopvar/misc/private | poscopvar_direct | 40 | 0 | MATTHEW M PEET |
| `check_reserved_names.m` | sopvar/misc/conventions | parse_copvar_spaces | 7 | 0 | Declan |
| `lpi_classify_sop.m` | lpi_programming_sopvar/misc | lpi_solve_loop_sop | 39 | 0 | MATTHEW M PEET |
| `lpigetdual_sop.m` | lpi_programming_sopvar/misc | exec_info_sop | 39 | 0 | MATTHEW M PEET |
| `lrmultiply.m` | sopvar/@sdopvar/private | mtimes | 11 | 0 | Sachin Shivakumar |
| `params2kernel.m` | sopvar/@sdopvar/private | sdopvar | 26 | 0 | Declan |
| `params2kernel.m` | sopvar/@sopvar/private | sopvar | 13 | 0 | Declan |
| `pie_rename_dummy.m` | parser/@pie_struct/private | initialize | 37 | 0 | Declan |
| `process_degrees_one.m` | lpi_programming_sopvar/misc/private | poscopvar_direct | 39 | 0 | MATTHEW M PEET |
| `weightPIE.m` | PIESIM/PIESIM_routines/PIESIM_Setup | PIESIM | 18 | 0 | Yulia |

`lpi_classify_sop` is also called by scripts under `cleaning/` (7 files); folding it removes their
access. `weightPIE` follows the layout of `PIESIM_Setup`, one file per setup step.

## Optional fold (7)

| file | folder | caller | code lines | test files | added by |
|---|---|---|---|---|---|
| `chebfit_mono.m` | sopvar/misc | inv | 57 | 0 | MATTHEW M PEET |
| `compute_canonical_BCs.m` | converters/internal/BVP | compute_Tmap_1D | 87 | 0 | Declan |
| `decouple_BVPs.m` | converters/internal/BVP | convert_PIETOOLS_BVP | 88 | 0 | Declan |
| `leftShiftMonomials_SD.m` | sopvar/@sdopvar/private | mtimes | 48 | 0 | Sachin Shivakumar |
| `leftShiftMonomials_SS.m` | sopvar/@sopvar/private | mtimes | 99 | 0 | Sachin Shivakumar |
| `lpi_soseq.m` | lpi_programming_sopvar/misc | impose_eq_rows | 57 | 0 | MATTHEW M PEET |
| `split_copvar.m` | sopvar/misc | convert_PIETOOLS_BVP | 46 | 0 | Declan |

## Kept with one caller (13)

Large (100 or more code lines):

| file | folder | caller | code lines | test files | added by |
|---|---|---|---|---|---|
| `compute_Tmap_1D.m` | converters/internal/BVP | convert_PIETOOLS_BVP | 210 | 0 | Declan |
| `eliminate_BVP_sources.m` | converters/internal/BVP | convert_PIETOOLS_BVP | 243 | 0 | Declan |
| `fund_inputs.m` | converters/internal/BVP | compose_Tmaps | 313 | 0 | Declan |
| `gk_grid.m` | sopvar/misc | inv | 115 | 0 | MATTHEW M PEET |
| `lift_BVP_sources.m` | converters/internal/BVP | convert_PIETOOLS_BVP | 138 | 0 | Declan |
| `pie_disc_sop.m` | executives_sopvar/utility_functions | pie_witness_sop | 140 | 0 | MATTHEW M PEET |
| `split_BVPs.m` | converters/internal/BVP | convert_PIETOOLS_BVP | 173 | 0 | Declan |

Called directly by tests:

| file | folder | caller | code lines | test files | added by |
|---|---|---|---|---|---|
| `eq_opts_sopvar.m` | lpi_programming_sopvar/misc | lpi_ineq_sop | 118 | 6 | MATTHEW M PEET |
| `is_canonical_multiplier.m` | sopvar/misc/canonical_form | collect_eq_rows | 128 | 1 | MATTHEW M PEET |
| `lpi_shape_sop.m` | lpi_programming_sopvar/misc | exec_info_sop | 21 | 1 | MATTHEW M PEET |
| `poslpivar_settings_2d_sop.m` | executives_sopvar/utility_functions | exec_tools_sop | 198 | 3 | MATTHEW M PEET |
| `poslpivar_settings_sop.m` | executives_sopvar/utility_functions | exec_tools_sop | 23 | 5 | MATTHEW M PEET |
| `settings2possopvar.m` | executives_sopvar/utility_functions | poslpivar_settings_2d_sop | 245 | 1 | MATTHEW M PEET |

## Owners

The 37 verdicts that change a file fall on the code of five authors (first commit of each file):

- Declan Jagt (17): the 8 copies (`pie_compat_tab`, `pie_compat_vars`, `pie_space_canon`,
  `display_term`); the folds `check_reserved_names`, `params2kernel` (2), `pie_rename_dummy`; the
  optional folds `compute_canonical_BCs`, `decouple_BVPs`, `split_copvar`; and
  `convert_PIETOOLS_PDE_nonlinear_2D`, `pie_pde_maps_2D` (not needed).
- Commits under M. Peet, sessions of this repository included (13): the lift chain (7, not needed);
  the folds `build_exponent_grid`, `process_degrees_one`, `lpi_classify_sop`, `lpigetdual_sop`; the
  optional folds `chebfit_mono`, `lpi_soseq`.
- Sachin Shivakumar (3): `lrmultiply` (fold), `leftShiftMonomials_SD`, `leftShiftMonomials_SS`
  (optional folds).
- Dan Braghini (2): the two periodic maps of the legacy converter (not needed).
- Yulia Peet (2): `clencurt` (not needed), `weightPIE` (fold).

The project rules send changes to the code of another developer to that developer as a report.
Nothing is changed by this review.

## Decisions and changes, 20 November 2026

Maintainer decisions on the ten folds: no change to the code of Declan Jagt (`check_reserved_names`,
`params2kernel` (2), `pie_rename_dummy` stay files); `lpigetdual_sop` and `lpi_classify_sop` stay
files; `build_exponent_grid` and `process_degrees_one` stay files while `copquadvar` calls them;
`clencurt` moves to the cleanup folder. The changes waited for the release comment pass to commit
the files (ece48c96, 98790f2e, 05:31-05:32).

Made (uncommitted):

- `lrmultiply` is a local function of `sopvar/@sdopvar/mtimes.m`; `@sdopvar/private/lrmultiply.m`
  deleted.
- `weightPIE` is a local function of `PIESIM/PIESIM_routines/PIESIM.m`, with a log entry and
  BEGIN/END markers (a 2025 file); `PIESIM_Setup/weightPIE.m` deleted. `mtimes.m` is new in 2026 and
  carries no change log under the release rules.
- `clencurt.m` moved to `cleaning/private/dead_code/PIETOOLS/PIESIM/PIESIM_routines/
  Third_Party_Supporting_Files/`.

Checks: the moved code equals the deleted files line for line (comments stripped), and each caller
keeps its HEAD code with the function appended. Products of decision and fixed containers in 1-D
and 2-D, and PIESIM on examples 39, 40 and 41 of `examples_pde_library_PIESIM_1D` (the three with
`uinput.weight`) give identical results before and after (isequaln); the profiler shows
`sdopvar.mtimes>lrmultiply` and `PIESIM>weightPIE` running. `which` finds no `lrmultiply`,
`weightPIE` or `clencurt` on the path. checkcode: `PIESIM.m` 0 messages; `mtimes.m` the 14 of HEAD
plus one inherited with `lrmultiply` (a preallocation of `Bout` overwritten two lines later).
