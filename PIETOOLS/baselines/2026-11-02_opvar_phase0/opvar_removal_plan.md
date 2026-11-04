# Removing opvar, dopvar, opvar2d and dopvar2d from PIETOOLS – inventory, PIESIM routes, work plan

Tree: `PIETOOLS/` on branch `ndopvar`, HEAD `eb5d0600` (merge of 10/30), with the uncommitted working-tree
edits listed in §1.3. Date: 2026-11-01. The investigation read the repository only: no file in it was
edited, created or deleted, and nothing was committed.

Labels: **[M]** measured – read from the code, or produced by a run reported in §1.2; **[I]** inferred –
a conclusion or an estimate that was not run.

---

## Summary

- **[M] Size.** 1,943 `.m` files on disk (1,874 tracked, 69 untracked). The legacy core is 156 files and
  19,808 code lines: the four classes (115 files, 14,260 lines), the helpers in `opvar/` and `opvar/2D/`
  (24 files, 1,901 lines) and `lpi_programming/` (17 files, 3,647 lines). Of the 1,215 files outside the
  core and the excluded folders, 193 reference the classes, their converters, their helpers or the stock
  LPI functions once comments are stripped; 114 more contain only field names that other classes share
  (`.P`, `.Q1`, `.dim`, `.I`, `.var1`).
- **[M] Excluded folders**, files with references: `PIETOOLS_demos/cuadmm` 36 of 112,
  `PIETOOLS_demos/lowrank_2d_stability` 20 of 80, `cleaning` 54 of 162, `SOSTOOLS400` 0 of 218.
- **[M] Root.** Every PDE→PIE path builds opvar or opvar2d objects; `pie_struct/initialize` converts them
  to copvar at the end, so a copvar PIE exists only after an opvar PIE has been built. With `opvar/` off
  the path, `convert_PIETOOLS_PDE` fails on all three input formats tried, and the native converter
  `convert_PIETOOLS_PDE_sopvar` fails as well. The legacy input formats (batch, legacy terms) still return
  PIEs whose operators are opvar: 38 of the 41 PIESIM 1-D examples and the GUI use them.
- **[M] Downstream, with `opvar/` and `lpi_programming/` off the path:** PIESIM fails at `pie2opvar`;
  `PIETOOLS_Hinf_gain_sop` and `PIETOOLS_H2_norm_sop` fail at the stock scalar `lpi_ineq`; `lpigetsol_sop`
  fails on a `dpvar`; `PIETOOLS_stability_sop` and `PIETOOLS_well_posedness_sop` still certify, but return
  the numerical witness as `NaN` and lose their opvar outputs without raising an error.
- **PIESIM: route (b)**, one converter from the copvar PIE to the kernel structs PIESIM already consumes,
  sharing a kernel extractor with the native converter and the executives' witness. [I] Cost: one new
  converter (≈200 lines), the shared extractor (≈250 lines), ≈25 changed lines in 4 PIESIM files; route (a)
  rewrites ≈600 code lines in 9 PIESIM files.
- **Phases:** 0 baseline fixtures; 1 kernel extractor and constructor; 2 PIE construction without opvar
  (critical path, the converters' owner); 3 container LPI layer without stock dispatch; 4 executives
  without opvar exits; 5 PIESIM; 6 PIESOS, ndopvar, tests, converters; 7 user-facing material; 8
  retirement of the classes. Phases 2 and 5 run in parallel after Phase 1; Phase 3 depends on no other
  phase and can start now.

---

## 1. Method

### 1.1 Static scan [M]

- **Scope.** Every `.m` file under `PIETOOLS/` except `.git`: 1,943 files. Script `scan.pl` (Perl 5.42).
- **Comment stripping.** Lines inside `%{ … %}` blocks; text after `%` outside a string literal; text after
  a `...` continuation. A `'` opens a string unless it follows an identifier character, `)`, `]`, `}`, `.`
  or `'`. String literals are kept and marked, so `isa(X,'opvar')` counts as a type test. "Code lines"
  below are the non-blank lines left after stripping.
- **Categories.** Constructions (`opvar(…)`, `opvar2d(…)`, command syntax `opvar P Q;`); class names in
  strings (type tests, messages); `InferiorClasses` entries; calls to 28 converter names in both
  directions; calls to the stock LPI functions (`poslpivar`, `lpivar`, `lpi_ineq`, `lpi_eq`,
  `getsol_lpivar`, `lpigetsol`, `lpiscript` and their `_2d` forms, `checkdeg_lpi_eq_2d`); calls to opvar
  helpers (`clean_opvar`, `apply_opvar`, `inv_opvar*`, `rand_opvar*`, `int_simple`, `getController`,
  `getObserver`, …); field reads `.R.R0`, `.R.R1`, `.R.R2` and the 16 opvar2d fields `.R00 … .R22`.
- **Dependent file:** at least one hit in these categories, after removing the hit on the file's own name
  (its function line) and calls to `piess`, `pielft` and `closedLoopPIE`, which accept copvar operators
  (session 2 ran `piess` on copvars; `closedLoopPIE` converts an opvar gain with `op2copvar_sop` and works
  on copvar otherwise).
- **Ambiguous only:** `.P`, `.Q1`, `.Q2`, `.dim`, `.I`, `.var1`, `.var2`. `pde_struct` and `pie_struct`
  also have `.dim`, so these names alone do not mark a file.
- **Exclusions:** `PIETOOLS_demos/cuadmm`, `PIETOOLS_demos/lowrank_2d_stability`, `cleaning`,
  `SOSTOOLS400`.
- **Limit.** Names built at run time (`feval`, `str2func`, `evalin` strings) are not followed, except the
  stock executive names, which were also searched as strings (§2.6).

### 1.2 MATLAB runs [M]

MATLAB R2025b, MOSEK 11.0.30, started with `-sd` in the scratchpad. Path: `restoredefaultpath`,
`addpath(genpath(PIETOOLS))`, then `cuadmm`, `lowrank_2d_stability`, `baselines` and `cleaning` removed.
Runs were serial, each under a timeout of 10–18 min that none reached.

- **Session 1, full path.** PIESIM 1-D examples 1, 10, 17, 23, 31, 38 and 2-D examples 1, 10:
  `convert_PIETOOLS_PDE`, then `PIESIM` (N 16 in 1-D, 8 in 2-D; tf 0.1; dt 1e-3);
  `convert_PIETOOLS_PDE_sopvar` against `convert_PIETOOLS_PDE`, operator by operator. On the 1-D PDE
  x_t = x_ss + 5x + w, z = ∫₀¹ x ds, x(0) = x(1) = 0: `PIETOOLS_stability_sop`, `PIETOOLS_Hinf_gain_sop`
  and the KYP LPI of DEMO9 written with container calls.
- **Session 2, `opvar/` and `lpi_programming/` removed from the path.** `lpiprogram`, `lpidecvar`,
  `lpisetobj` and `lpisolve` were re-added from a scratch copy, since none of them references the four
  classes. First failure of 22 flows. Two of its probes were invalid – the check script overwrote the
  settings variable with `syms st` – and session 3b repeats them.
- **Session 3a, full path.** Every `.mat` outside `.git` and `cleaning` (51 files) loaded and searched
  recursively for objects of the four classes.
- **Session 3b, as session 2.** Four container executives on the PIE of session 1, outputs recorded.
- **Session 4, full path.** Native-converter coverage on 96 PDE sources: PIESIM 1-D 39–41 (the `pde_var`
  ones), PIESIM 2-D 1–32, and 61 function-form files of `PIETOOLS_examples/Examples_Library` called as
  `[PDE_t,PDE_b] = f(0,{})`. The first run stopped at the second PDE with a MATLAB "Illegal instruction"
  fault in `m_lxe.dll` (MATLAB's JIT); the rerun from that PDE completed without a fault. 33 of the 61
  library files did not accept that calling convention (27 "Too many output arguments", 6 others); they
  are not counted below.

### 1.3 Tree state [M]

Uncommitted edits by other sessions while the runs took place: `SOSTOOLS400/sossolve.m`,
`executives_sopvar/utility_functions/exec_info_sop.m`, `.../lpi_solve_loop_sop.m`,
`lpi_programming_sopvar/misc/lpi_classify_sop.m`, `.../misc/lpi_eq_cdopvar.m`, `.../misc/lpi_eq_sdopvar.m`,
`.../misc/private/collect_eq_rows.m`; untracked `baselines/` and
`SOSTOOLS400/internal/processing/cuadmm_sedumi.m`. The executive runs used these files as they were. By the
end of the investigation `git status` also listed `lpi_programming_sopvar/README.md` and
`lpi_programming_sopvar/tests/test_lpi_eq_sop.m` as modified, by another session.

---

## 2. Inventory

### 2.1 Totals [M]

| Group | Files | Code lines | Files with references |
|---|---:|---:|---:|
| `opvar/@opvar` | 30 | 1,677 | – |
| `opvar/@dopvar` | 25 | 1,213 | – |
| `opvar/2D/@opvar2d` | 33 | 6,413 | – |
| `opvar/2D/@dopvar2d` | 27 | 4,957 | – |
| helpers `opvar/*.m` (14), `opvar/2D/*.m` (10) | 24 | 1,901 | – |
| `lpi_programming/` (11) and `lpi_programming/2D` (6) | 17 | 3,647 | – |
| **Core** | **156** | **19,808** | |
| Outside the core and the exclusions | 1,215 | | **193** (+114 ambiguous only) |
| Excluded | 572 | | 110 |

- 4 of the 17 `lpi_programming` files – `lpiprogram`, `lpidecvar`, `lpisetobj`, `lpisolve` – contain no
  reference to the four classes and are called from 60, 35, 40 and 56 files outside the core, the
  container DEMOs among them.
- `opvar/@quadPoly` contains no `.m` file; the `quadPoly` on the path is `ndopvar/@quadPoly`
  (`which -all`). `objs/` contains no file, `polyopvar/` one `.asv`, `testfolder/@testobj` no file.
- Stored data: 0 of 51 `.mat` files hold an object of the four classes (session 3a). A `.mat` that holds
  one loads it as `double` when the class is absent (session 2, a saved 1-D PIE).

### 2.2 By subsystem [M]

Columns: `.m` files scanned; files with references; then call sites – constructions (function or command
syntax), class names in strings, `InferiorClasses` entries, converter calls, stock-LPI calls, helper calls,
1-D field reads (`.R.R*`), 2-D field reads (`.R00 … .R22`).

| Subsystem | Files | Dep. | New | Str | Inf | Conv | LPI | Help | Fld 1-D | Fld 2-D | Note |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `converters/@pde_struct`, `@pie_struct`, `@npie_struct` | 60 | 8 | 6 | 19 | 3 | 3 | 0 | 0 | 3 | 0 | class definitions |
| `converters/operator_converters` | 16 | 8 | 7 | 22 | 0 | 10 | 0 | 0 | 9 | 0 | |
| `converters/` other (PDE→PIE, `internal`, `tests`) | 36 | 18 | 112 | 56 | 0 | 10 | 0 | 39 | 158 | 84 | DJ |
| `converters_sopvar/` | 12 | 2 | 17 | 1 | 0 | 7 | 0 | 12 | 24 | 0 | DB, 10/28–10/29 |
| `executives/` (stock 1-D) | 20 | 19 | 0 | 0 | 0 | 22 | 223 | 4 | 4 | 0 | assigned to Declan |
| `executives/2D` | 14 | 14 | 10 | 0 | 0 | 12 | 179 | 10 | 0 | 4 | report-only |
| `executives/utility_functions` | 6 | 6 | 0 | 26 | 0 | 1 | 0 | 4 | 15 | 83 | 2-D helpers report-only |
| `executives_sopvar/` | 28 | 17 | 2 | 22 | 0 | 62 | 10 | 0 | 3 | 1 | |
| `lpi_programming_sopvar/` | 56 | 19 | 5 | 42 | 0 | 45 | 51 | 8 | 8 | 0 | 11 of 19 are tests |
| `sopvar/@*` (container classes) | 109 | 4 | 4 | 4 | 0 | 3 | 0 | 0 | 3 | 0 | `copvar2opvar(2d)`, `sopvar2opvar(2d)` |
| `sopvar/misc`, `sopvar/constructors` | 36 | 0 | | | | | | | | | |
| `sopvar/Testfolder` | 40 | 12 | 14 | 8 | 0 | 50 | 0 | 29 | 40 | 0 | opvar as test oracle |
| `ndopvar/` | 98 | 14 | 6 | 13 | 0 | 42 | 5 | 9 | 6 | 42 | 5 are error-message strings |
| `PIESIM/` | 62 | 7 | 1 | 3 | 0 | 3 | 0 | 0 | 21 | 160 | §3 |
| `PIESOS/` | 334 | 20 | 11 | 20 | 0 | 21 | 2 | 14 | 44 | 0 | testing deferred (10/23) |
| `parser/` (legacy `sys`, `state`, `terms`) | 50 | 13 | 10 | 0 | 0 | 0 | 0 | 0 | 42 | 0 | |
| `GUI/` (`.m` files) | 4 | 0 | | | | | | | | | `.mlapp`: §2.7 |
| `PIETOOLS_demos/` DEMO files | 11 | 3 | 1 | 0 | 0 | 4 | 0 | 0 | 1 | 0 | DEMO2, 7, 11 |
| `PIETOOLS_demos/snippets_from_manual` | 9 | 3 | 12 | 0 | 0 | 0 | 9 | 1 | 18 | 12 | Ch3, Ch5, Ch7 |
| `PIETOOLS_demos/sopvar_demos` | 18 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 3 | 1 | |
| `PIETOOLS_demos/multivariate_PIE_2508_14840` | 1 | 1 | 0 | 0 | 0 | 0 | 5 | 0 | 0 | 0 | |
| `PIETOOLS_examples/` | 99 | 0 | | | | | | | | | 1 runs a stock executive by string (§2.6) |
| `settings/`, `settings_sopvar/`, `misc/` | 23 | 0 | | | | | | | | | |
| root scripts | 5 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 4 | `PIETOOLS_PDE.m`; 2 more call stock executives |
| `baselines/` (untracked) | 68 | 3 | 0 | 0 | 0 | 5 | 1 | 0 | 3 | 2 | |

The per-file list is Appendix A.

### 2.3 Call sites that matter [M]

**(a) Constructions.**

- `converters/convert_PIETOOLS_PDE.m`: 17 (1 `opvar()`, 10 `opvar2d(…)`, 6 command syntax), and in its
  internals `internal/2D/PIETOOLS_FTC_Expansion_2D.m` (3; lines 225 and 238 are the first failures without
  opvar for 1-D and 2-D `pde_struct` input), `internal/2D/Compute_Tmap_PDE_2D_Separable.m` (2),
  `internal/compute_PIETOOLS_PDE_periodic_map.m` (3), `internal/2D/compute_PIETOOLS_PDE_periodic_map_2D.m` (6).
- Legacy converters: `internal/legacy/convert_PIETOOLS_PDE_batch.m` (3), `…_terms_internal.m` (8),
  `…_terms_legacy.m` (9), `convert_PIETOOLS_DDE2PIE_legacy.m` (1).
- `converters/convert_PIETOOLS_DDF.m`: one command-syntax line declaring 13 opvars (line 152) and 36
  `.R.R*` field accesses; DDE and NDS inputs pass through it (`convert_PIETOOLS_DDE:69`,
  `convert_PIETOOLS_NDS:88`).
- `converters/convert_PIETOOLS_PDE_nonlinear.m`: 15, converted with `opvar2sopvar` into `tensopvar`.
- `converters_sopvar/convert_PIETOOLS_PDE_sopvar.m`: 17, in its fallback (a copy of `convert_PIETOOLS_PDE`).
- `converters/initialize_PIETOOLS_PIE.m`: 40 (20 `opvar`, 20 `opvar2d`, command syntax); no reference to
  it outside its own file, comments included, `cleaning/` excluded.
- `converters/@pie_struct/expand_tderivatives.m`: 6; it converts the PIE with `copvar2opvar(2d)` (lines
  107, 109) and expands higher-order time derivatives in opvar.
- `PIESIM/PIESIM_routines/PIESIM_Setup/weightPIE.m`: `opvar W` (35) and opvar products `W*PIE.T`,
  `W*PIE.Tu`, `W*PIE.Tw` (44–46).
- `executives_sopvar/utility_functions`: `exec_tools_sop.m` (an `opvar2d` and two `mat2opvar`, identities
  for opvar PIEs, 277–281 and 335), `getController_direct_sop.m` (`opvar Kop`, 180),
  `getController_sop.m`/`getObserver_sop.m` (`mat2opvar` unit vectors, 125/124).
- Converters to opvar: `sopvar/@copvar/copvar2opvar.m`, `copvar2opvar2d.m`, `sopvar/@sopvar/sopvar2opvar.m`,
  `sopvar2opvar2d.m`, `converters/operator_converters/pie2opvar.m` (3).
- User-facing: DEMO2 (`opvar Top`, 53); snippets Ch3 (6), Ch5 (6); `parser/` (10).

**(b) Field reads.** PIESIM 181 in 6 files (§3.2); converters 242; `executives/utility_functions` 98
(`get_eq_opts_2D`, `get_lpivar_degs`, stock); `PIETOOLS_PDE.m:100–106` (`PIE.T.P`, `.R.R0`, `.R00`,
`.Rxx{1}`, `.Ryy{1}`, `.R22{1}`); `executives_sopvar/utility_functions/pie_disc_sop.m` (`R.R0–R2`, opvar
input only, 94–108); `converters_sopvar/tmap_separable_sopvar.m:271–288` (`copvar2opvar`, then
`.P`, `.Q1`, `.Q2`, `.R.R*`); `converters/@pde_struct/mtimes.m:144–147` (`C.R.R0–R2`).

**(c) Type tests and class names in strings.** `pie_struct/initialize` (98, 136; messages 115, 151),
`pie_struct/display` (104–107), `pie_struct/private/pie_rename_dummy` (55), `pde_struct/mtimes` (62–113),
`pde_struct/times` (54–56), `closedLoopPIE`, `piess` (38, in messages and one type test), `op2copvar_sop`, `opvar2sopvar`,
`getsol_lpivar_sop` (92), `lpigetsol_sop` (96–102), `trace_rn_sop` (72–77), PIESOS `tensopmat` (153, 175),
`tensopmat/initialize` (161), `tensopvar/initialize` (61), `piesos_ineq` (132), `piesos_getsol` (85, 92),
`PIESIM_input_check` (469), `transl_mod` (58); "opvar objects" in 5 error messages of `@ndopvar` and
`@nopvar`.

**(d) Converters in either direction** (call sites by folder, the converter files' own lines included).

| Converter | Direction | Call sites |
|---|---|---|
| `op2copvar_sop` | opvar, opvar2d, dopvar → copvar; a copvar passes through | `executives_sopvar` 44, `lpi_programming_sopvar` 10, `converters` 3 (`pie_struct/initialize`, `pie_pde_maps` ×2), `converters_sopvar` 3, DEMO2/DEMO7 3, `closedLoopPIE` 1 |
| `opvar2copvar` | opvar → copvar | `lpi_programming_sopvar` 29 (tests), `executives_sopvar` 9, `sopvar` 8 |
| `opvar2sopvar` | opvar → sopvar | PIESOS 18, `converters` 5 (nonlinear converter 3), `sopvar` 6 |
| `opvar2d2sopvar`, `opvar2d2copvar` | opvar2d → containers | `converters` 3, `sopvar` 3 + 8, `lpi_programming_sopvar` 3 |
| `copvar2opvar`, `copvar2opvar2d` | copvar → opvar, opvar2d | `executives_sopvar` 3 + 1, `converters` 2 + 2, `converters_sopvar` 1, `sopvar` 7 + 4 |
| `sopvar2opvar`, `sopvar2opvar2d` | sopvar → opvar, opvar2d | PIESOS 2, `sopvar` 8 + 6 |
| `pie2opvar` | copvar PIE → struct of opvars | PIESIM 3, `pie_witness_sop` 1, PIESOS 1 |
| `pie2copvar` | opvar PIE → copvar PIE | DEMO11 1, `converters_sopvar` 1 |
| `mat2opvar` | matrix → opvar | stock executives 34, `converters` 5, `executives_sopvar` 4, `lpi_programming_sopvar` 3, `sopvar` 3, `converters_sopvar` 2 |
| `dopvar2ndopvar`, `dopvar2d2ndopvar`, `ndopvar2dopvar` | | `ndopvar/Testfolder` 41 |

**(e) `InferiorClasses`.** `pde_struct`: `{?state, ?terms, ?opvar, ?opvar2d, ?dpvar, ?polynomial}`;
`pie_struct`: `{?opvar2d, ?dpvar, ?polynomial}` (copvar is not listed). `parser/@terms` declares the
property `operator opvar;` (line 35). `sopvar`, `sdopvar`, `copvar`, `cdopvar`, `nopvar`, `ndopvar` and the
PIESOS classes list none of the four.

**(f) Stock LPI functions and stock utilities reached from container code.**

| Caller | Stock function | When | Effect without `lpi_programming/` |
|---|---|---|---|
| 6 container executives (H2 control, H2 estimator, H2 norm, H∞ control, H∞ estimator, H∞ gain), 9 sites | `lpi_ineq(prog,gam…)` on a `dpvar` | the γ constraint | `PIETOOLS_Hinf_gain_sop:251`, `PIETOOLS_H2_norm_sop:262`: "Undefined function 'lpi_ineq'" |
| `lpi_ineq_sop:202–204` | `lpi_ineq` | P not a container | – |
| `lpi_eq_sop:104–105` | `lpi_eq` | P not a container (incl. `dpvar`) | – |
| `getsol_lpivar_sop:92` | `getsol_lpivar` | opvar classes | – |
| `lpigetsol_sop:96–103` | `lpigetsol` (→ `sosgetsol` for `dpvar`, char, cell) | `dpvar`, double, polynomial, char, cell, opvar classes | "Undefined function 'lpigetsol' for input arguments of type 'dpvar'"; DEMO2, 3, 5, 6, 8, 9 call it on the `dpvar` γ |
| `poslpivar_sop:104–106` | `poslpivar` | n not a container | – |
| `get_lpivar_degs_sop:97–98` | `get_lpivar_degs` | R not a container | – |
| `trace_rn_sop:72–75` | field trace | X legacy | – |
| `PIETOOLS_Hinf_estimator_sop:163` | `getObserver_2D` (opvar2d inverse) | 2-D | – (`@copvar/inv` is 1-D: error `copvar:inv:nd`) |
| `pie_witness_sop:145` | `pie2opvar`, then `pie_disc_sop` (opvar only) | witness on | certificate kept, `info.maxre = NaN`, no error |
| `exec_tools_sop:186/197` (`tl.toop`) | `copvar2opvar(2d)` | the outputs P, R, Q, W, Z | `stability_sop` returns the copvar; `well_posedness_sop` returns `[]` for P and R |
| `getController_sop`, `getObserver_sop` | `copvar2opvar`, `mat2opvar` | `Kop`/`Lop` output; `qp` unit vectors | – |

### 2.4 Run-time dependency [M]

`opvar/` and `lpi_programming/` off the path (sessions 2 and 3b):

| Flow | Result | First failure |
|---|---|---|
| `convert_PIETOOLS_PDE`, batch format (PIESIM 1-D example 1) | fails | `convert_PIETOOLS_PDE_batch:191` |
| `convert`, `pde_struct`, 1-D (the PDE of §1.2) | fails | `PIETOOLS_FTC_Expansion_2D:225` ← `convert_PIETOOLS_PDE:453` |
| `convert_PIETOOLS_PDE`, 2-D (PIESIM 2-D example 1) | fails | `PIETOOLS_FTC_Expansion_2D:238` |
| `convert_PIETOOLS_PDE_sopvar`, legacy formats (6 PDEs) | fails | `convert_PIETOOLS_PDE_batch:191`, `…_terms_legacy:121` |
| `convert_PIETOOLS_PDE_sopvar`, 2-D (examples 1, 10) | fails | `PIETOOLS_FTC_Expansion_2D:238` ← `convert_PIETOOLS_PDE_sopvar:454`, the fallback after the native attempt |
| `initialize`, `pie2copvar`, `piess` on a copvar PIE | run | – |
| `PIESIM(PIE)` with a copvar PIE, 1-D and 2-D | fails | `pie2opvar:155`, `:158` ← `PIESIM:190` |
| `PIESIM(PDE)`, 2-D | fails | in `convert` (above) |
| KYP LPI with `poscopvar`, `lpi_ineq_sop`, `lpisolve` | solves; extraction fails | `lpigetsol_sop:103` |
| `PIETOOLS_stability_sop`, form Q | certified; `maxre = NaN`; P returned as copvar | errors caught at `PIETOOLS_stability_sop:205–210` and in `pie_witness_sop:144–147` |
| `PIETOOLS_well_posedness_sop` | certified; `maxre = NaN`; P, R returned as `[]` | caught |
| `PIETOOLS_Hinf_gain_sop`, form Q | fails | `PIETOOLS_Hinf_gain_sop:251` (`lpi_ineq`) |
| `PIETOOLS_H2_norm_sop`, form o | fails | `PIETOOLS_H2_norm_sop:262` (`lpi_ineq`) |
| `lpiscript_sop(PIE,'stability','light')` | certified; "numerical spectrum max Re NaN" | caught |

Full path (session 1): `PIETOOLS_stability_sop` certified in 1.3 s; `PIETOOLS_Hinf_gain_sop` γ = 0.16772443
in 1.5 s; the KYP LPI γ = 0.16772409 in 0.5 s. `class(PIE.T)` from `convert_PIETOOLS_PDE` is `opvar` for
the 6 PIESIM 1-D examples (legacy formats) and `copvar` for the 2 PIESIM 2-D examples. PIESIM times:
1-D 0.31–2.15 s, 2-D 3.71–4.50 s.

**Native converter coverage (sessions 1 and 4, full path).** On 63 `pde_struct` PDEs that convert with
`convert_PIETOOLS_PDE`:

| Source | PDEs | Native | Max coefficient difference vs `convert`, 12 operators | Fallback |
|---|---:|---:|---|---|
| PIESIM 1-D 39–41 | 3 | 3 | ≤ 1.1e-14, all dimensions equal | – |
| Example library, `PDE_t` (1-D) | 28 | 28 | ≤ 1.1e-14, all dimensions equal | – |
| PIESIM 2-D 1–32 | 32 | 14 | ≤ 1.1e-14, all dimensions equal | 18: BCs that need derivatives of inputs (14), composition with a derivative not expressible as a PI operator (3), mixed BCs (1); each fallback then fails with "Unrecognized … field 'F' for class 'pie_struct'" |

The 12 operators are T, Tw, Tu, A, Bw, Bu, Cz, Dzw, Dzu, Cy, Dyw, Dyu. `convert_PIETOOLS_PDE_sopvar` has no
caller in the repository. The PIESIM legacy-format examples (1-D 1–38) fail in its legacy branch on the
same removed field `PIE.F` (6 of 6 tried, session 1).

### 2.5 Changes since the inventory of 10/27 (`2901956c`) [M]

- DEMO1–DEMO11 call the container API (`lpiprogram_sop`, `poscopvar`, `lpivar_cdopvar`, `lpi_ineq_sop`,
  `lpigetsol_sop`, `lpiscript_sop`). Three still use legacy names: DEMO2 (opvar construction and
  `op2copvar_sop`), DEMO7 (`op2copvar_sop` on the gains `lpiscript_sop` returns), DEMO11 (`pie2copvar` on a
  copvar PIE).
- `pie_struct` holds copvar operators for `pde_struct`, DDE and DDF input (`initialize` converts them);
  19 of the 20 files in `executives/` read `<operator>.dim`, which copvar does not have (the 20th is the
  script `PIETOOLS_auto_execute`).
- `closedLoopPIE` and `piess` accept copvar operators (from the code; `piess` ran on copvars in session 2);
  `pielft` names none of the four classes, and DEMO5 and DEMO7 call it on copvar PIEs (not run here).
- `converters_sopvar/` (DB, 10/28–10/29) adds the native copvar converter, not called by `convert`.
- PIESIM and `pie_witness_sop` convert copvar PIEs back to opvar (`pie2opvar`).
- Unchanged: the legacy formats return opvar PIEs; the container executives return opvar outputs, call
  the stock scalar `lpi_ineq` and dispatch to stock routines in the four `_sop` LPI functions.

### 2.6 Stock executives and their callers [M]

Assigned to Declan; no edits are proposed for them.

- `executives/` (20): `PIETOOLS_FEPDEstability`, `_FEPIE2PDEstability`, `_H2_control`, `_H2_estimator`,
  `_H2_norm_c`, `_H2_norm_c_coercive`, `_H2_norm_o`, `_H2_norm_o_coercive`, `_Hinf_control`,
  `_Hinf_estimator`, `_Hinf_gain`, `_Hinf_gain_coercive`, `_Hinf_gain_dual`, `_Hinf_gain_dual_coercive`,
  `_PDEstability`, `_PDEstability_dual`, `_PIE2PDEstability`, `_PIE2PDEstability_dual`,
  `_well_posedness`, and the script `_auto_execute`.
- `executives/2D` (14, report-only): `_FEPDEstability_2D`, `_FEPIE2PDEstability_2D`, `_H2_norm_2D_c`,
  `_H2_norm_2D_c_non_coercive`, `_H2_norm_2D_o`, `_H2_norm_2D_o_non_coercive`, `_Hinf_estimator_2D`,
  `_Hinf_gain_2D`, `_Hinf_gain_2D_non_coercive`, `_Hinf_gain_dual_2D`, `_PDEstability_2D`,
  `_PDEstability_dual_2D`, `_PIE2PDEstability_2D`, `_PIE2PDEstability_dual_2D`.
- `executives/utility_functions`: `getController`, `getObserver` (stock 1-D, `inv_opvar`), `getObserver_2D`
  and `get_eq_opts_2D` (2-D, report-only), `get_lpivar_degs`, `closedLoopPIE` (copvar-capable).
- Together they hold 402 stock-LPI call sites (223 in 1-D, 179 in 2-D).
- Callers outside `executives/`: `PIETOOLS_DDE.m` (10 calls, lines 58–85) and `PIETOOLS_DDF.m`. Neither
  root script runs as written: `PIETOOLS_DDF.m` calls `PIETOOLS_stability`, `PIETOOLS_stability_dual` and
  `convert_PIETOOLS_DDE2PIE`, and `PIETOOLS_DDE.m` calls `PIETOOLS_H2_controller`; none of the four exists
  in the tree. `PIETOOLS_examples/examples_PDE_library_PIETOOLS.m:1006, 1015`
  (`evalin('base','PIETOOLS_auto_execute')`); `lpi_programming_sopvar/tests/e2e_hinf_gain_1d.m` and
  `test_dpvar_hinf_chain_sop.m` (reference values). No file outside `executives/` and `lpi_programming/`
  calls the stock `lpiscript`.

### 2.7 Material outside `.m` files [M]

- **GUI.** `PIETOOLS_PDE_GUI.mlapp` (`matlab/document.xml`, line 1391) calls `convert_PIETOOLS_PDE_terms`,
  which calls `convert_PIETOOLS_PDE_terms_legacy` and returns an opvar PIE.
- **Manual.** `PIETOOLS_documentation/PIETOOLS_2025_User_Manual.pdf` (213 pages, last changed 06/02/2026):
  opvar 275, dopvar 41, opvar2d 41, dopvar2d 15; sopvar, sdopvar, copvar, cdopvar 0; poslpivar 29,
  lpi_ineq 25, lpiscript 25. The 227-page copy at the repository root differs as a file and gives the same
  class counts. The TeX source is not in the repository. By chapter (pdftotext):

| Chapter | Pages | opvar | dopvar | opvar2d | dopvar2d | poslpivar | lpi_ineq | lpiscript |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| front matter | 1–8 | 8 | 3 | 2 | 0 | 1 | 1 | 0 |
| 2 Scope of PIETOOLS | 12–18 | 0 | 0 | 0 | 0 | 0 | 0 | 3 |
| 3 PI Operators in PIETOOLS | 19–30 | 23 | 1 | 19 | 1 | 0 | 0 | 0 |
| 4 Setup and Representation of PDEs and DDEs | 31–45 | 25 | 0 | 0 | 0 | 0 | 0 | 0 |
| 5 Conversion of PDEs and DDEs to PIEs and Closing the Loop | 46–62 | 126 | 0 | 0 | 0 | 0 | 0 | 0 |
| 6 PIESIM: A General-Purpose Simulation Tool | 63–82 | 1 | 0 | 0 | 0 | 0 | 0 | 1 |
| 7 Declaring and Solving Convex Optimization Programs on PI Operators | 83–105 | 22 | 11 | 9 | 10 | 21 | 14 | 9 |
| 8 Input Formats for ODE-PDE Systems | 106–141 | 1 | 0 | 1 | 0 | 0 | 0 | 0 |
| 10 Operations on PI Operators in PIETOOLS: opvar and dopvar | 149–160 | 59 | 18 | 10 | 4 | 1 | 0 | 0 |
| 11 PIETOOLS Demonstrations | 161–186 | 2 | 0 | 0 | 0 | 6 | 10 | 1 |
| 13 Standard Applications of LPI Programming | 191–201 | 8 | 8 | 0 | 0 | 0 | 0 | 11 |

  Chapters 1, 9, 12 and Appendix A: none. The mentions in chapter 4 are `pie_struct` listings such as
  `T: [2×2 opvar]`.

---

## 3. PIESIM

### 3.1 Data flow today [M]

`PIESIM(PDE | DDE | PIE, opts, uinput)` → for a PDE, `convert_PIETOOLS_PDE` → for a `pie_struct`,
`pie2opvar` (`PIESIM.m:172`; DDE `:182`; PIE `:190`) → `PIESIM_input_check` → `weightPIE` (1-D, optional)
→ `rescalePIE` / `rescalePIE_2D` (affine map of the domain onto [−1, 1] by `transl_mod` / `transl_2D`) →
`PIESIM_discretize_all(_2D)` → `PIESIM_discretize_ops(_2D)` → one call per operator to
`PIESIM_4PI2Mat_cheb` (1-D) or `PIESIM_fullPI2Mat_cheb_2D` (2-D), which pass polynomial pieces to
`PIESIM_Poly2Mat_cheb`, `PIESIM_PI2Mat_opint_cheb`, `PIESIM_3PI2Mat_cheb` and the 9-PI and 1-D→2-D
routines → time integration on the matrices `Dop` → `PIESIM_transform_to_solution(_2D)`.

- `pie2opvar` (113 code lines) returns a plain struct whose fields T, Tw, Tu, A, B1, B2, C1, D11, D12, C2,
  D21, D22 (and K) are opvar or opvar2d objects, assembled block by block with `copvar2opvar` /
  `copvar2opvar2d`, plus `dom`, `vars`, `dim` and the tables `x_tab … z_tab` sorted by space kind.
- Legacy-format PIEs arrive with opvar operators and skip `pie2opvar` (session 1: examples 1, 10, 17, 23,
  31, 38). Of the 41 PIESIM 1-D examples, 16 take the batch branch of `convert_PIETOOLS_PDE` (fields
  `n0`, `n1` or `n2`), 22 the legacy-terms branch and 3 are
  `pde_var`; all 32 2-D examples use `pde_var`.
- Without opvar, every PIESIM call fails (§2.4).

### 3.2 Routines that read opvar or opvar2d data [M]

| File | Code lines | Data read | Class-specific code | Role |
|---|---:|---|---|---|
| `PIESIM.m` | 54 | – | `pie2opvar` ×3 | entry |
| `PIESIM_Setup/PIESIM_input_check.m` | 714 | `PIE.T.dim` ×4, `C1.Q1`, `C2.Q1` ×2 each | `isa(PIE.T,'opvar2d')` (469); `size(PIE.Tu/Tw/C1/C2,k)` with the opvar overload (352–366, 497–501) | sizes for type PIE and DDE |
| `PIESIM_Setup/weightPIE.m` | 18 | `T.R.R0`, `T.I` | `opvar W` (35); opvar products (44–46) | weighted formulation |
| `PIESIM_Setup/transl_mod.m` | 31 | `Q1`, `Q2`, `R.R0–R2`, `I`, `var1`, `var2` | `isa(T,'opvar')` (58) | 1-D rescaling |
| `PIESIM_Setup/transl_2D.m` | 107 | the 16 fields, `I`, `var1`, `var2` (58 reads) | – | 2-D rescaling |
| `Discretization_Of_Operators/PIESIM_4PI2Mat_cheb.m` | 53 | `P`, `Q1`, `Q2`, `R` (17 reads) | – | 1-D Chebyshev matrix of a 4-PI operator |
| `2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m` | 217 | the 16 fields, `var1` (65 reads) | – | 2-D Chebyshev matrix |

The other 55 PIESIM files take polynomial matrices or arrays; `rescalePIE(_2D)` and
`PIESIM_discretize_ops(_2D)` read PIE fields and hand them to the routines above.

### 3.3 Field meaning and the same data in a container [M: class headers, `sopvar2opvar`, `apply_sopvar`]

A copvar block `C{i,j}` is a sopvar with `vars.out = [S2, S3]`, `vars.in = [S3, S1]` (S1: variables of the
input only, integrated; S2: variables of the output only; S3: shared), exponent cells `ZL` (output
variables) and `ZR` (input variables), `dims = [q p]`, `dom`, and `params`, a cell of 3^|S3| sparse
matrices. For a multi-index α the kernel is

  K_α(S2, S3; S3′, S1) = (I_q ⊗ Z_L(S2, S3))ᵀ · C_α · (I_p ⊗ Z_R(S3′, S1)),

where index 1 along a variable s of S3 is the multiplier (s′ = s), 2 the integral from the lower bound to
s, 3 the integral from s to the upper bound. In canonical multiplier form a multiplier index carries no
`ZR` degree in that variable. The container itself holds `space_out`, `space_in` (which variables each
component depends on), `dim_out`, `dim_in`, `vars` and `dom`.

**1-D.** opvar: R^{m0} × L2^{m1}[a,b] ← R^{n0} × L2^{n1}[a,b]; `dim = [m0 n0; m1 n1]`, `I = [a b]`,
`var1 = s`, `var2 = θ`.

| opvar field | Operator piece | copvar block | sopvar data |
|---|---|---|---|
| `P` (m0×n0) | x0 ↦ P x0 | output and input space R | S1 = S2 = S3 = ∅; `params{1}` |
| `Q1(s)` (m0×n1) | x1 ↦ ∫ₐᵇ Q1(s) x1(s) ds | output R, input L2[s] | S1 = {s}; `params{1}`, `ZR` |
| `Q2(s)` (m1×n0) | x0 ↦ Q2(s) x0 | output L2[s], input R | S2 = {s}; `params{1}`, `ZL` |
| `R.R0(s)` | multiplier | L2[s] → L2[s] | S3 = {s}; `params{1}` |
| `R.R1(s,θ)` | ∫ₐˢ | | `params{2}` |
| `R.R2(s,θ)` | ∫ₛᵇ | | `params{3}` |
| `dim` | component counts | `dim_out`, `dim_in`, `space_out`, `space_in` | `dims` per block |
| `I`, `var1`, `var2` | domain, variables | `dom`, `vars` | `dom`, `vars` (dummy names are implicit) |

An opvar stacks all R^n components before all L2 components; a copvar keeps the PIE's component order,
in which the two kinds may alternate. `pie2opvar` sorts the components by space kind and keeps their
order within a kind.

**2-D.** opvar2d: spaces R, L2[x], L2[y], L2[x,y]; `dim` 4×2; field `R_ab`, a the output space, b the
input space, a, b ∈ {0, x, y, 2}.

| Field | S2 | S3 | S1 | `params` | opvar2d cell |
|---|---|---|---|---:|---|
| `R00` | ∅ | ∅ | ∅ | 1 | – |
| `R0x`, `R0y`, `R02` | ∅ | ∅ | {x}, {y}, {x,y} | 1 | – |
| `Rx0`, `Ry0`, `R20` | {x}, {y}, {x,y} | ∅ | ∅ | 1 | – |
| `Rxx`, `Ryy` | ∅ | {x}, {y} | ∅ | 3 | `{1:3}`: multiplier, lower, upper |
| `Rxy`, `Ryx` | {x}, {y} | ∅ | {y}, {x} | 1 | – |
| `Rx2`, `Ry2` | ∅ | {x}, {y} | {y}, {x} | 3 | `{1:3}` |
| `R2x`, `R2y` | {y}, {x} | {x}, {y} | ∅ | 3 | `{1:3}` |
| `R22` | ∅ | {x,y} | ∅ | 9 | `{i,j}`: i along x, j along y |

Orientation: for |S3| = 2, `params` is indexed along the variables in order of their **names**
(`sopvar_from_kernels` header; the `possopvar` sorted-S3 convention), whereas opvar2d indexes `R22{i,j}`
by the PIE's (x, y) (`copvar2opvar2d(Pc, vars_xy)`). When the names sort as (y, x) the 3×3 cell is
transposed.

### 3.4 What a conversion takes [M: `sopvar2opvar`, `sopvar2opvar2d`, `pie2opvar`]

For each block: build Z_L and Z_R from the exponent cells, form K_α for each α (`quadPoly` →
`polynomial`), and place it in the slot of the block's (output, input) space kinds, at the row and column
offset of its components; in 2-D, map the cell index by the variable order. PIE operators carry no
decision variables, so the cost grows with blocks × monomials and not with the decision-variable count
[I].

### 3.5 Routes

**(a) PIESIM reads copvar PIEs directly.** The discretizers loop over the blocks of each copvar operator,
extract each block's kernels and place the Chebyshev matrices at the block offsets, in PIESIM's state order
(R^n components first, PDE components by differentiability order, read from `x_tab`). Rescaling acts on the
extracted kernels; a copvar-level affine change of variables does not exist (no such method in `@copvar`
or `@sopvar`). `weightPIE` uses a copvar multiplier and copvar `mtimes`; `PIESIM_input_check` reads sizes
from `dim_out`, `dim_in`, `space_out`.

**(b) Converter to PIESIM's internal arrays.** A new converter (name to choose, e.g. `pie2piesim`)
returns the operators as plain structs with the field layout PIESIM consumes today – 1-D `P, Q1, Q2,
R.R0, R.R1, R.R2, dim, I, var1, var2`; 2-D `R00 … R22` with opvar2d's cells, `dim, I, var1, var2` –
holding `polynomial` and `double` values, with no class. The blocks are read with the Phase-1 extractor;
the placement logic is that of `pie2opvar`'s `assemble_op`. PIESIM then changes only where it relies on
the class: `PIESIM.m` (3 calls), `transl_mod.m:58` (class check), `weightPIE.m:35–46` (constructor and
products: apply the weight as a copvar multiplier before the conversion), `PIESIM_input_check.m`
(class test at 469, overloaded `size` at 352–366 and 497–501).

| | (a) direct | (b) converter |
|---|---|---|
| Files edited, 1-D [I] | `PIESIM.m`, `PIESIM_input_check`, `weightPIE`, `transl_mod`, `PIESIM_discretize_ops`, `PIESIM_4PI2Mat_cheb` (6) | `PIESIM.m`, `PIESIM_input_check`, `weightPIE`, `transl_mod` (4) |
| Files edited, 2-D [I] | + `PIESIM_discretize_ops_2D`, `PIESIM_fullPI2Mat_cheb_2D`, `transl_2D` (3) | none beyond the 1-D list |
| New files [I] | the extractor (Phase 1) | the extractor (Phase 1) and the converter |
| Lines [I] | ≈600 code lines rewritten (53 + 217 + 26 + 18 + 31 + 107 + 18 + ≈150 of `PIESIM_input_check`) | ≈25 lines changed, ≈200 new |
| Files unchanged | the polynomial-level routines | all 55 polynomial-level, integration and post-processing files, and `transl_2D`, `rescalePIE(_2D)` |
| Risk, 1-D [I] | block placement in PIESIM's state order, spread over the discretizers | the same placement in one function, whose logic `pie2opvar` has run since 10/29 |
| Risk, 2-D [I] | 16 slots, cell orientation and placement spread over 3 files | 16 slots and orientation in one function; orientation test (§3.6) |
| Use for N-D [I] | a block-based discretizer is the starting point of an N-D PIESIM | none; PIESIM stays 1-D/2-D |

**Recommendation: route (b).** Metric: ≈225 changed or new lines in 5 files against ≈600 rewritten lines
in 9 files, and the 55 polynomial-level files keep their inputs; the extractor of Phase 1 also serves the
native converter (Phase 2) and the witness (Phase 4). Limit: PIESIM keeps its internal field names, and an
N-D PIESIM would need the block-based discretizer of route (a); neither route gives PIESIM a third
spatial variable.

### 3.6 Verification

- **Before and after, trajectories.** Phase 0 records, for all 41 1-D and 32 2-D library examples at fixed
  options, the time histories in `solution`, the final states, and every discretized operator in `Dop`
  (`Acheb`, `Tcheb`, `B1cheb`, …). After the change: max over time steps and grid points of
  |x_after − x_before| / max |x_before| for states and outputs, and the entrywise difference of each `Dop`
  matrix. [I] Acceptance at 1e-10 relative: the change reorganizes the same polynomial kernels.
- **Independent references.** (i) The error against `uinput.exact` must not increase: 35 of the 41 1-D
  examples (1–32, 39–41) and all 32 2-D examples define an exact solution [M]. (ii) The 38 legacy-format 1-D
  examples were built in opvar without any container code [M], so their Phase-0 trajectories do not share
  code with the new converter.
- **Orientation (2-D).** Rename the variables of a 2-D example with a non-symmetric `R22` so that their
  names sort opposite to (x, y); the trajectories must equal the original ones under the renaming.
- **Negative control.** In a scratch copy, exchange R1 and R2 in the converter; the trajectory check must
  fail.
- **Not used.** A reverse converter applied to the output, or a comparison of the new converter with
  `pie2opvar` alone: if both read the blocks with the same extraction code, a shared index error passes.
- **Time** [I, from the session-1 unit costs]: ≈4 min per pass (41 × ≈1 s in 1-D, 32 × ≈4–5 s in 2-D).

---

## 4. Work plan

Each phase lists its files, what replaces opvar there, its dependencies, its verification down the chain
pvar/dpvar → operator classes → `lpi_programming_sopvar` → executives → examples, and its size [I]. Every
edit follows the modification protocol of `CLAUDE.md` §1 (minimal change, header entry, column-77
markers, old code kept). Dead code of other developers moves to `cleaning/private/dead_code/` with a
README rather than being deleted.

### Phase 0 – Baseline fixtures (before any edit)

- **Files.** None in the repository: scripts and fixtures in a test location chosen by the maintainer.
- **Content.** F1, the PIESIM library (§3.6), ≈4 min per pass. F2, every library PIE as built today
  (copvar `.mat`) and its action on 5 random polynomial test functions per operator, evaluated with
  `apply_sopvar`; for the 38 legacy-format 1-D examples also with `apply_opvar` on their opvar operators,
  which is independent of the container code. F3, the container rows of the stored example-library
  baseline (`baselines/2026-10-28_example_library`) re-run on the current tree: status, γ, m, `K.s`, nnz,
  witness. F4, DEMO1–DEMO11 outputs (γ values, closed-loop trajectories).
- **Time.** [M, `baselines/2026-10-28_example_library/runs.log`] consecutive container runs of the 130
  library rows started 3.4 min (light → heavy) and 14.7 min (heavy → veryheavy) apart; the baseline skips
  2-D rows whose projected time per solve exceeds 600 s (its README). The stored rows come from trees before
  `3df33dcc` (`b515277f`, `8354ecb0`), so F3 re-runs them on the current tree.
- **Verification.** Two passes give identical fixtures, apart from timings.
- **Tree.** Record on a named commit. The main checkout carried, on 11/01–11/02, another session's
  uncommitted equality-row change (`collect_eq_rows`), which lowers the row count m and flips some
  library verdicts through row order (its 11/01 note). F3 and F4 depend on it; F1 and F2 do not.
  A later comparison of m on a tree with a different LPI layer mixes the two effects, so F3 and F4
  are re-recorded whenever the LPI layer underneath changes.
- **Size** [I]: 3 scripts, ≈300 lines. **Parallel:** gates every edit below; Phase 1 design can start.
- **Done 11/02**: `PIETOOLS/baselines/2026-11-02_opvar_phase0/` (untracked, beside the 10/28
  baseline; README there), on a frozen export of `eb5d0600` without the uncommitted edits. Recorded:
  F1, 72 of 73 examples (1-D example 32 errors inside PIESIM); F2, 172 PIEs (on the 38 opvar PIEs
  `apply_opvar` and `apply_sopvar` agree to a relative difference of 2.6e-14); F3 at light, 135 rows
  (PDE 36 and 38 skipped; PDE 50, 2,421 s for one solve, in pass 1 only); F4, the 11 DEMOs. Pass 2
  equals pass 1 on F1, F2 and the 134 F3 rows of both passes. Not recorded: F3 at heavy (estimated
  29 min per pass) and the PDE 50 dual. No fixture file stores an object of the four classes.
- **Re-recorded 11/04** on `c81c3d04` (passes 3 and 4, the reference for every later phase; passes 1
  and 2 kept). Pass 4 equals pass 3 on F1, F2 and the 134 F3 rows. Against pass 1: the converter of
  `b4256de6` splits PIE states into more blocks and reorders the state components of 10 PIEs; by
  components and up to that order, every F2 operator agrees to at most 3.5e-12. F3: m lower on 116
  rows (−3.9% in total), 5 status changes, 17 further result changes, all on rows whose PIE is
  unchanged, so from the LPI layer (`e3a1f667`) and the `parse_copvar_spaces` fix. Phases 1–2 change
  the block partition again: compare F2 with `p0_f2cmp`, not with the stored outputs.

### Phase 1 – Kernel extraction and construction for containers

- **Files.** New `sopvar/misc/copvar_kernels.m` (copvar → plain struct of kernel polynomials: the 1-D and
  2-D layouts of §3.3, and per-block kernel cells for any number of variables), taken from the
  extraction in `sopvar2opvar.m` / `sopvar2opvar2d.m` without the class constructors. Promote
  `converters_sopvar/private/sopvar_from_kernels.m` and `assemble_copvar_blocks.m` to a public
  constructor (kernels → copvar), e.g. in `sopvar/constructors/`.
- **Replaces.** `copvar2opvar(2d)` and `sopvar2opvar(2d)` as kernel readers (callers:
  `tmap_separable_sopvar`, `pie2opvar`, `expand_tderivatives`, `exec_tools_sop`, `getController_sop`,
  `getObserver_sop`); opvar construction with field assignment as a kernel writer (DEMO2, snippets Ch3 and
  Ch5, `getController_direct_sop`, `convert_PIETOOLS_DDF`).
- **Depends on** Phase 0 (F2).
- **Verification.** Each direction against a class definition, never a round trip (`CLAUDE.md` §4).
  (i) Extractor: on random copvars with 1, 2 and 3 variables, the extracted kernels, applied to random
  polynomial test functions by exact integration (with `apply_opvar` while the class exists), equal
  `apply_sopvar` of the blocks (the sopvar definition) to 1e-12. (ii) The forward converter on which F2
  relies: on the 38 legacy-format 1-D PIEs, `apply_opvar` of the original opvar operators equals
  `apply_sopvar` of `op2copvar_sop(PIE.X)` – two evaluators, each from its own class definition. Comparing
  kernels extracted from `op2copvar_sop(PIE.X)` with the original fields is a round trip of two mutually
  inverse routines and is not used. (iii) Negative controls that must fail: R1 and R2 exchanged; ZL and ZR
  exchanged on blocks whose bases differ in size (equal sizes hide the exchange, `CLAUDE.md` §5); the 2-D
  cell transposed for names that sort as (y, x). (iv) Constructor: `apply_sopvar` of the constructed
  copvar equals direct evaluation of the input kernels.
- **Performance** [I]: fixed operators only; cost in blocks × monomials, no array scaled by the number of
  decision variables.
- **Size** [I]: 2–3 new files, ≈400 lines; one test, ≈200 lines.

### Phase 2 – PIE construction without opvar (the root; critical path)

- **Owner.** The converters are DJ's; the message of `3df33dcc` says a new converter that supports N-D PDEs
  is in progress [M]. Phase 2 is that work or the steps below; the maintainer assigns it.
- **2a – `pde_struct` input.** Make the native construction of `converters_sopvar/` the path of
  `convert_PIETOOLS_PDE` for `pde_struct` input, after three repairs: `tmap_separable_sopvar.m:271` reads
  kernels through `copvar2opvar` (→ Phase-1 extractor); `pie_ops2copvar` reads the removed `PIE.F`; field
  names of the 10/29 `pie_struct` (`B1` → `Bw`, `F` → `pde2pie.x.C`, writes to `x_tab`). Then close the
  gaps: [M] 18 of the 32 PIESIM 2-D PDEs (BCs that need derivatives of inputs 14, compositions with a
  derivative 3, mixed BCs 1); [M, from the code] higher-order time derivatives (port
  `pie_struct/expand_tderivatives`, which runs in opvar), PDEs without spatial variables, generalized
  periodic BCs with ODE coupling or boundary inputs. Coverage today [M]: 45 of 63 `pde_struct` PDEs,
  operators equal to `convert`'s to ≤ 1.1e-14.
- **2b – Legacy formats.** Batch (`n0`, `n1`, `n2` → `convert_PIETOOLS_PDE_batch`) and legacy terms
  (`PDE.PDE`, `PDE.ODE` → `convert_PIETOOLS_PDE_terms_legacy`) return opvar PIEs without `initialize` [M].
  Users: 38 of 41 PIESIM 1-D examples, the `PDE_b` output of the example-library functions, the GUI.
  Options: translators to `pde_struct` followed by 2a, or retirement of the formats with the 38 examples
  and the GUI moved to `pde_var`.
- **2c – DDE, DDF, NDS.** Rebuild the operator assembly of `convert_PIETOOLS_DDF.m` (13 opvars, 36 `.R.R*`
  field accesses) with `mat2copvar_sop` and the Phase-1 constructor. `convert_PIETOOLS_DDE2PIE_legacy`, called only
  by `PIETOOLS_DDE.m:53`, retires.
- **2d – Nonlinear PDEs.** `convert_PIETOOLS_PDE_nonlinear.m` (15 opvars, `opvar2sopvar` → `tensopvar`):
  build the sopvar blocks with the helpers of 2a.
- **2e – Class definitions.** Remove `?opvar` and `?opvar2d` from `pde_struct` and `?opvar2d` from
  `pie_struct`; `pie_struct/initialize` raises an error for opvar operators (98–151); `display`
  (104–107); `pie_rename_dummy` (55); `pde_struct/mtimes` and `times` (opvar multipliers of PDE terms:
  accept copvar through the extractor, or drop); `pie_pde_maps` (2 `op2copvar_sop`); `piess` (the
  `piess_opvar` branch); `initialize_PIETOOLS_PIE.m` (no caller) to `cleaning/`.
- **Depends on** Phase 1 and the decision on 2b.
- **Verification.** (i) Operators of every library PIE against F2: coefficient difference and action on
  the test functions. (ii) `lpi_programming_sopvar` tests and end-to-end tests. (iii) Executives: the 23
  library builds give identical SDPs (m, `K.s`, nnz) and values, and F3 holds. (iv) PIESIM trajectories
  against F1. (v) Every library file converts; DEMO1–DEMO11 run. (vi) The configuration of session 2
  (`opvar/` off the path) converts all inputs.
- **Size** [I]: ≈20 files. Largest items: DDF ≈300 lines, nonlinear ≈500, `expand_tderivatives` ≈450,
  two legacy translators ≈300 each (if kept), the 2a gaps ≈300–600.
- **Parallel:** 2c, 2d and 2e edit disjoint files.

### Phase 3 – Container LPI layer without stock dispatch

- **Files.** `lpi_programming_sopvar/lpi_ineq_sop.m:202–204` (a `dpvar` or `polynomial` goes to `sosineq`,
  as in the stock `lpi_ineq`; legacy classes raise an error); `lpi_eq_sop.m:104–105` (`dpvar` → `soseq`);
  `lpigetsol_sop.m:96–103` (`dpvar`, char, cell → `sosgetsol`; double and polynomial returned unchanged,
  as in the stock `lpigetsol`); `getsol_lpivar_sop.m:92` (drop the legacy branch). Move `lpiprogram`,
  `lpidecvar`, `lpisetobj`, `lpisolve` out of `lpi_programming/` (into `lpi_programming_sopvar/`, for
  example). `lpiscript_sop`: drop the `sys` branch; keep the call to the stock `lpisettings` for
  'extreme', which references no class.
- **Tests.** The 11 test files of `lpi_programming_sopvar/tests` that compare with stock routines compare
  with F2 and F3 instead.
- **Depends on** nothing; it can start now.
- **Verification.** Unit tests; DEMO2, 3, 5, 6, 8, 9 (each calls `lpigetsol_sop` on the `dpvar` γ) with
  `lpi_programming/` off the path; the KYP LPI of session 1 returns γ = 0.16772409 [M with the full path].
- **Performance** [I]: scalar branches; nothing on the decision-variable axis changes.
- **Size** [I]: 4 files, ≈40 changed lines; 4 files moved; 11 tests adjusted.

### Phase 4 – Container executives without opvar exits

- **γ constraints.** Six executives, 9 sites: `lpi_ineq(prog,gam…)` → `lpi_ineq_sop` (after Phase 3).
- **Outputs.** Drop `tl.toop` (`exec_tools_sop:186, 197`); the executives return the containers P, R, Q, W,
  Z, K, L. `getController_sop` / `getObserver_sop`: drop the opvar outputs `Kop` / `Lop` (146 / 145) and
  build the unit vectors of the `qp` method with `mat2copvar_sop` (125 / 124). `getController_direct_sop`:
  build K with the Phase-1 constructor (180). The return type is user-facing: DEMO7:124 converts the gains
  with `op2copvar_sop`; snippets Ch2 and Ch6 pass them to `closedLoopPIE`, which takes copvar (Phase 7).
- **Witness.** `pie_witness_sop:145` (`pie2opvar`) and `pie_disc_sop` (opvar only, 94) read the kernels
  with the Phase-1 extractor. Without opvar the witness returns `NaN` without an error [M].
- **Opvar PIE support.** `exec_pie_sop` (the `initialize` branch), the identities of `exec_tools_sop`
  (277–281, 335), and the 44 `op2copvar_sop(PIE.X)` pass-through calls, once Phase 2 guarantees copvar PIEs.
- **Legacy fallbacks.** `poslpivar_sop:104–106`, `get_lpivar_degs_sop:97–98`, `trace_rn_sop:72–75` raise
  an error instead.
- **2-D H∞ estimator gain.** `PIETOOLS_Hinf_estimator_sop:163` calls the stock `getObserver_2D` (opvar2d
  inverse); `@copvar/inv` is 1-D [M]. Options: no 2-D gain (`[]` with a message, as in N-D), or a container
  2-D inverse.
- **Swallowed errors.** The `try … catch` blocks around extraction and witness (e.g.
  `PIETOOLS_stability_sop:205–210`, `pie_witness_sop:144–147`) hide a missing class. The tests assert that
  `info.witness.ok` is true in 1-D and that the outputs are non-empty containers.
- **Tests.** `executives_sopvar/tests/test_executives_sop.m` (the stock-slack comparison) compares with F3.
- **Depends on** Phases 1 and 3; the removal of opvar-PIE support on Phase 2.
- **Verification.** `test_executives_sop(false)` (20 checks); F3: same SDPs (m, `K.s`, nnz), same values,
  same witness values; the configuration of session 2 runs every executive with no error and no `NaN`
  witness in 1-D.
- **Performance** [I]: removing `tl.toop` removes a conversion per output; the LPIs are unchanged.
- **Size** [I]: 9 executives and 8 utilities, ≈150 changed lines.

### Phase 5 – PIESIM on container PIEs (route b of §3)

- **Files.** One new converter in `PIESIM/PIESIM_routines/PIESIM_Setup/`; edits in `PIESIM.m` (172, 182,
  190), `transl_mod.m` (58), `weightPIE.m` (35–46), `PIESIM_input_check.m` (352–366, 469, 497–501).
- **Depends on** Phase 1; the 38 legacy-format 1-D examples on Phase 2b, since their PIEs are opvar today.
- **Verification.** §3.6.
- **Size** [I]: ≈200 new lines, ≈25 changed lines.

### Phase 6 – Remaining code

- **PIESOS (20 files).** Drop the opvar branches of `tensopmat` (153, 175), `tensopmat/initialize` (161),
  `tensopvar/initialize` (61), `piesos_ineq` (132–137), `piesos_getsol` (85–94).
  `@intop/apply_op.m:277` calls `opvar/2D/int_simple.m`, a `dpvar`/`polynomial` integrator that references
  no class: move it out of `opvar/` before Phase 8, whatever the PIESOS schedule. `PIE2polyPIE.m`
  (`pie2opvar`, `opvar2sopvar`) and 12 example or test files that build opvars use the Phase-1 constructor.
  PIESOS testing is deferred (maintainer, 10/23).
- **ndopvar (14).** 5 error messages name "opvar objects" (text only). `Testfolder/converters`
  (`dopvar2ndopvar`, `dopvar2d2ndopvar`, `ndopvar2dopvar`) and 6 tests use opvar as the reference:
  fixtures, or retirement with the classes.
- **Container classes (4).** `@copvar/copvar2opvar`, `copvar2opvar2d`, `@sopvar/sopvar2opvar`,
  `sopvar2opvar2d` retire when no caller remains.
- **`sopvar/Testfolder` (12).** The opvar classes are the independent oracle of the container tests
  (`apply_opvar` against `apply_sopvar`, `Test_copvar_converters`, `Bench_copvar_vs_legacy`,
  `Test_copvar_inv`). Replace with the F2 fixtures, or keep a copy of the classes off the default path that
  only these tests add.
- **Operator converters (8 of 16).** `op2copvar_sop`, `opvar2copvar`, `opvar2sopvar`, `opvar2d2sopvar`,
  `opvar2d2copvar`, `dopvar2sdopvar`, `sdopvar2dopvar`, `pie2copvar`, `pie2opvar` retire when their callers
  are gone.
- **Legacy parser (13 of 50).** `sys`, `state`, `terms` build opvars; retire them together with
  `state2pde_struct`, `terms2pde_struct` and the `sys` branches of `convert_PIETOOLS_PDE` and
  `lpiscript_sop`.
- **`converters/tests` (2).** Class assertions `'opvar'`, `'opvar2d'` → copvar.
- **`lpi_programming_sopvar/examples`.** `volterra_norm_sop` (opvar construction) and 3 files with
  `op2copvar_sop` on PIE fields.
- **Parallel:** disjoint folders. **Size** [I]: ≈50 files, most of them tests.

### Phase 7 – User-facing material

- **DEMOs (3 of 11).** DEMO2 (53–60: the Volterra operator through the Phase-1 constructor); DEMO7 (124:
  remove `op2copvar_sop` once the executives return containers); DEMO11 (17: remove `pie2copvar`).
- **Snippets.** Ch3 (3 opvar, 3 opvar2d, field assignments); Ch5 (6 opvar, field access, `piess`,
  `closedLoopPIE`); Ch7 (stock LPI API → container API); Ch2 and Ch6 need no edit once the executives return
  containers.
- **Root scripts.** `PIETOOLS_PDE.m:100–106` (opvar fields size the initial condition → `PIE.space.x` or
  `T.dim_in`); `PIETOOLS_DDE.m`, `PIETOOLS_DDF.m` (stock executives → `_sop` executives or `lpiscript_sop`,
  and the four missing functions of §2.6).
- **Example library.** `examples_PDE_library_PIETOOLS.m:1006, 1015` → `PIETOOLS_auto_execute_sop`.
- **GUI.** The `.mlapp` follows the decision on 2b; it is edited in App Designer.
- **Demo folders.** `multivariate_PIE_2508_14840/jp_wave_plate_rates.m` (stock `lpi_eq`, `lpi_ineq`,
  `lpivar`, `lpivar_2d`, `sys`) and `sopvar_demos/test_heatNd_pie.m` (class string, fields).
- **Manual.** Chapters 3, 4, 5, 7, 10, 11, 13 and the front matter (§2.7), written as if the four classes
  never existed: no conversion sections, correspondence tables or "earlier releases", and none of
  `poslpivar`, `lpivar`, `op2copvar_sop`, `copvar2opvar`, `pie2opvar` named (rule of 11/02). Chapter 10
  needs the container counterpart of each operation it lists. The TeX source is outside the repository;
  chapters cite only the manual (rule of 11/01).
- **Depends on** Phases 2–5 (final API and return types).
- **Verification.** Every DEMO and snippet runs with `opvar/` and `lpi_programming/` off the path, with
  values equal to F4; the rebuilt PDF contains none of the four class names and none of the five function
  names above.
- **Size** [I]: ≈12 code files; 7 chapters.

### Phase 8 – Retirement of the classes

- **Files.** Move `opvar/` (classes and helpers) and the 13 `lpi_programming/` files other than the four
  moved in Phase 3 (12 reference the classes; `lpiscript` dispatches to the stock executives) to
  `cleaning/private/dead_code/` with a README; remove the empty `objs/opvar` and `opvar/@quadPoly`
  folders.
- **Precondition.** The stock executives (34 files, 402 stock-LPI call sites) are ported or retired by
  their owner; otherwise they stop resolving `poslpivar`, `lpi_ineq`, … when `lpi_programming/` leaves the
  path.
- **Verification.** The full suite on the default path: tests, DEMOs, snippets, F1, F3; `which('opvar')`
  is empty; the scan of §1.1 finds none of the four names in code outside `cleaning/`.

### Order and parallelism

| Phase | Folders edited | Needs | Runs in parallel with |
|---|---|---|---|
| 0 Baseline fixtures | none | – | design of 1 and 3 |
| 1 Kernel extractor and constructor | `sopvar/` | 0 | 3 |
| 2 PIE construction | `converters/`, `converters_sopvar/` | 1; the decision on 2b | 3, 4, 5, 6 |
| 3 LPI layer | `lpi_programming_sopvar/` (root), `lpi_programming/` (4 moves) | – | 0, 1, 2, 5, 6 |
| 4 Executives | `executives_sopvar/` | 1, 3; 2 for removing opvar-PIE support | 2, 5, 6 |
| 5 PIESIM | `PIESIM/` | 1; 2b for the legacy-format examples | 2, 3, 4, 6 |
| 6 Remaining code | `PIESOS/`, `ndopvar/Testfolder`, `sopvar/Testfolder`, `converters/operator_converters`, `parser/` | 1; 4 for the executive tests | 2, 3, 4, 5 |
| 7 User-facing | `PIETOOLS_demos/`, root scripts, `PIETOOLS_examples/` (1 file), GUI, manual | 2, 4, 5 | – |
| 8 Retirement | `opvar/`, `lpi_programming/`, `cleaning/` | all; the stock executives' owner | – |

[I, from the file lists] No two phases that run in parallel edit the same file; Phase 4 removes the
`op2copvar_sop` call sites in `executives_sopvar/`, and Phase 6 removes the converter file afterwards.

### Decisions for the maintainer

1. Who carries Phase 2: DJ's N-D converter or the native path of `converters_sopvar/`.
2. Legacy input formats: translate to `pde_struct`, or retire (38 PIESIM 1-D examples, the `PDE_b`
   outputs of the example library, the GUI).
3. Executive outputs as containers instead of opvar (a user-facing change).
4. 2-D H∞ estimator gain: none, or a container 2-D inverse.
5. Test oracle after the removal: fixtures, or an off-path copy of the classes for tests only.
6. The stock executives (Declan): ported or retired before Phase 8.
7. Access to the manual's TeX source.

---

## Appendix A – Files with references, by subsystem [M]

Comment-stripped; hits on a file's own name excluded. Tags: `new` constructions; `str` class names in
strings; `inf` `InferiorClasses` entries; `cls` other bare uses; `conv(x)`, `lpi(x)`, `help(x)` calls;
`fld` 1-D field reads; `fld2` 2-D field reads.

**PDE/PIE structs (converters/@*)** (8 files)

- `converters/@pde_struct/mtimes.m`: fld=3 str=4
- `converters/@pde_struct/pde_struct.m`: inf=2
- `converters/@pde_struct/times.m`: str=2
- `converters/@pie_struct/display.m`: str=4
- `converters/@pie_struct/expand_tderivatives.m`: conv(copvar2opvar)=1 conv(copvar2opvar2d)=1 new=6
- `converters/@pie_struct/initialize.m`: conv(op2copvar_sop)=1 str=8
- `converters/@pie_struct/pie_struct.m`: inf=1
- `converters/@pie_struct/private/pie_rename_dummy.m`: str=1

**operator_converters** (8 files)

- `converters/operator_converters/dopvar2sdopvar.m`: conv(opvar2dopvar)=1 fld=3 str=5
- `converters/operator_converters/op2copvar_sop.m`: conv(opvar2copvar)=1 conv(opvar2d2sopvar)=1 new=1 str=5
- `converters/operator_converters/opvar2copvar.m`: conv(opvar2sopvar)=1 new=1 str=2
- `converters/operator_converters/opvar2d2copvar.m`: conv(opvar2d2sopvar)=1 new=1 str=3
- `converters/operator_converters/opvar2d2sopvar.m`: conv(opvar2sopvar)=1 str=3
- `converters/operator_converters/opvar2sopvar.m`: conv(dopvar2sdopvar)=1 conv(opvar2d2sopvar)=1 fld=3 str=4
- `converters/operator_converters/pie2opvar.m`: conv(copvar2opvar)=1 conv(copvar2opvar2d)=1 new=3
- `converters/operator_converters/sdopvar2dopvar.m`: fld=3 new=1

**converters (PDE->PIE, internal, tests)** (18 files)

- `converters/convert_PIETOOLS_DDF.m`: fld=36 new=1
- `converters/convert_PIETOOLS_PDE.m`: conv(mat2opvar)=2 fld=21 help(clean_opvar)=12 new=17
- `converters/convert_PIETOOLS_PDE_nonlinear.m`: conv(opvar2sopvar)=3 fld=15 help(clean_opvar)=10 new=15
- `converters/initialize_PIETOOLS_PIE.m`: new=40
- `converters/internal/2D/Compute_Tmap_PDE_2D_Separable.m`: conv(mat2opvar)=1 fld2=60 help(clean_opvar)=2 new=2
- `converters/internal/2D/PIETOOLS_FTC_Expansion_2D.m`: fld=6 fld2=18 new=3
- `converters/internal/2D/compute_PIETOOLS_PDE_periodic_map_2D.m`: fld2=5 new=6
- `converters/internal/compute_PIETOOLS_PDE_periodic_map.m`: fld=6 new=3
- `converters/internal/legacy/convert_PIETOOLS_DDE2PIE_legacy.m`: fld=36 new=1
- `converters/internal/legacy/convert_PIETOOLS_PDE_batch.m`: fld=9 help(op_slice)=12 new=3
- `converters/internal/legacy/convert_PIETOOLS_PDE_terms_internal.m`: fld=8 new=8
- `converters/internal/legacy/convert_PIETOOLS_PDE_terms_legacy.m`: fld=12 new=9
- `converters/internal/pie_pde_maps.m`: conv(op2copvar_sop)=2
- `converters/piess.m`: help(piess_opvar)=1 str=38
- `converters/piess_opvar.m`: conv(mat2opvar)=2 help(subs_vars_op)=2 new=4 str=12
- `converters/terms2pde_struct.m`: fld=3 str=2
- `converters/tests/test_convert_PIETOOLS_PDE_periodic.m`: fld=3 fld2=1 str=3
- `converters/tests/test_convert_PIETOOLS_PDE_periodic_higher_order.m`: fld=3 str=1

**converters_sopvar** (2 files)

- `converters_sopvar/convert_PIETOOLS_PDE_sopvar.m`: conv(mat2opvar)=2 conv(op2copvar_sop)=3 conv(pie2copvar)=1 fld=21 help(clean_opvar)=12 new=17 str=1
- `converters_sopvar/tmap_separable_sopvar.m`: conv(copvar2opvar)=1 fld=3

**executives (stock 1-D)** (19 files)

- `executives/PIETOOLS_FEPDEstability.m`: lpi(getsol_lpivar)=1 lpi(lpi_eq)=1 lpi(lpi_ineq)=1 lpi(poslpivar)=4
- `executives/PIETOOLS_FEPIE2PDEstability.m`: lpi(getsol_lpivar)=1 lpi(lpi_eq)=2 lpi(lpi_ineq)=1 lpi(lpivar)=1 lpi(poslpivar)=4
- `executives/PIETOOLS_H2_control.m`: conv(mat2opvar)=1 fld=2 help(getController)=1 lpi(lpi_eq)=2 lpi(lpi_ineq)=4 lpi(lpigetsol)=4 lpi(lpivar)=2 lpi(poslpivar)=6
- `executives/PIETOOLS_H2_estimator.m`: conv(mat2opvar)=1 fld=2 help(getObserver)=1 lpi(lpi_eq)=2 lpi(lpi_ineq)=4 lpi(lpigetsol)=4 lpi(lpivar)=2 lpi(poslpivar)=6
- `executives/PIETOOLS_H2_norm_c.m`: conv(mat2opvar)=1 lpi(lpi_eq)=3 lpi(lpi_ineq)=4 lpi(lpigetsol)=4 lpi(lpivar)=1 lpi(poslpivar)=7
- `executives/PIETOOLS_H2_norm_c_coercive.m`: lpi(lpi_eq)=1 lpi(lpi_ineq)=3 lpi(lpigetsol)=2 lpi(poslpivar)=4
- `executives/PIETOOLS_H2_norm_o.m`: conv(mat2opvar)=1 lpi(lpi_eq)=3 lpi(lpi_ineq)=3 lpi(lpigetsol)=4 lpi(lpivar)=1 lpi(poslpivar)=7
- `executives/PIETOOLS_H2_norm_o_coercive.m`: lpi(lpi_eq)=1 lpi(lpi_ineq)=3 lpi(lpigetsol)=2 lpi(poslpivar)=4
- `executives/PIETOOLS_Hinf_control.m`: conv(mat2opvar)=3 help(getController)=1 lpi(lpi_eq)=1 lpi(lpi_ineq)=2 lpi(lpigetsol)=4 lpi(lpivar)=1 lpi(poslpivar)=4
- `executives/PIETOOLS_Hinf_estimator.m`: conv(mat2opvar)=3 help(getObserver)=1 lpi(lpi_eq)=1 lpi(lpi_ineq)=2 lpi(lpigetsol)=4 lpi(lpivar)=1 lpi(poslpivar)=4
- `executives/PIETOOLS_Hinf_gain.m`: conv(mat2opvar)=2 lpi(lpi_eq)=2 lpi(lpi_ineq)=2 lpi(lpigetsol)=3 lpi(lpivar)=1 lpi(poslpivar)=4
- `executives/PIETOOLS_Hinf_gain_coercive.m`: conv(mat2opvar)=2 lpi(lpi_eq)=1 lpi(lpi_ineq)=2 lpi(lpigetsol)=3 lpi(poslpivar)=4
- `executives/PIETOOLS_Hinf_gain_dual.m`: conv(mat2opvar)=2 lpi(lpi_eq)=2 lpi(lpi_ineq)=2 lpi(lpigetsol)=3 lpi(lpivar)=1 lpi(poslpivar)=4
- `executives/PIETOOLS_Hinf_gain_dual_coercive.m`: conv(mat2opvar)=2 lpi(lpi_eq)=1 lpi(lpi_ineq)=2 lpi(lpigetsol)=3 lpi(poslpivar)=4
- `executives/PIETOOLS_PDEstability.m`: conv(mat2opvar)=1 lpi(getsol_lpivar)=1 lpi(lpi_eq)=1 lpi(lpi_ineq)=1 lpi(poslpivar)=4
- `executives/PIETOOLS_PDEstability_dual.m`: conv(mat2opvar)=1 lpi(getsol_lpivar)=1 lpi(lpi_eq)=1 lpi(lpi_ineq)=1 lpi(poslpivar)=4
- `executives/PIETOOLS_PIE2PDEstability.m`: lpi(getsol_lpivar)=1 lpi(lpi_eq)=2 lpi(lpi_ineq)=1 lpi(lpivar)=1 lpi(poslpivar)=4
- `executives/PIETOOLS_PIE2PDEstability_dual.m`: lpi(getsol_lpivar)=1 lpi(lpi_eq)=2 lpi(lpi_ineq)=1 lpi(lpivar)=1 lpi(poslpivar)=4
- `executives/PIETOOLS_well_posedness.m`: conv(mat2opvar)=2 lpi(getsol_lpivar)=2 lpi(lpi_eq)=2 lpi(lpi_ineq)=2 lpi(poslpivar)=8

**executives/2D (stock)** (14 files)

- `executives/2D/PIETOOLS_FEPDEstability_2D.m`: help(clean_opvar)=1 lpi(checkdeg_lpi_eq_2d)=2 lpi(lpi_eq_2d)=1 lpi(lpi_ineq_2d)=1 lpi(poslpivar_2d)=5
- `executives/2D/PIETOOLS_FEPIE2PDEstability_2D.m`: help(clean_opvar)=1 lpi(lpi_eq_2d)=2 lpi(lpi_ineq_2d)=2 lpi(lpivar_2d)=1 lpi(poslpivar_2d)=4
- `executives/2D/PIETOOLS_H2_norm_2D_c.m`: fld2=1 lpi(lpi_eq_2d)=1 lpi(lpi_ineq)=1 lpi(lpi_ineq_2d)=1 lpi(lpigetsol)=2 lpi(poslpivar_2d)=4 new=1
- `executives/2D/PIETOOLS_H2_norm_2D_c_non_coercive.m`: conv(mat2opvar)=1 fld2=1 lpi(lpi_eq_2d)=3 lpi(lpi_ineq)=1 lpi(lpi_ineq_2d)=2 lpi(lpigetsol)=4 lpi(lpivar_2d)=1 lpi(poslpivar_2d)=7
- `executives/2D/PIETOOLS_H2_norm_2D_o.m`: fld2=1 lpi(lpi_eq_2d)=1 lpi(lpi_ineq)=1 lpi(lpi_ineq_2d)=1 lpi(lpigetsol)=2 lpi(poslpivar_2d)=4 new=1
- `executives/2D/PIETOOLS_H2_norm_2D_o_non_coercive.m`: conv(mat2opvar)=1 fld2=1 lpi(lpi_eq_2d)=3 lpi(lpi_ineq)=1 lpi(lpi_ineq_2d)=2 lpi(lpigetsol)=4 lpi(lpivar_2d)=1 lpi(poslpivar_2d)=7
- `executives/2D/PIETOOLS_Hinf_estimator_2D.m`: conv(mat2opvar)=2 help(clean_opvar)=1 lpi(checkdeg_lpi_eq_2d)=2 lpi(lpi_eq_2d)=2 lpi(lpi_ineq_2d)=2 lpi(lpigetsol)=8 lpi(lpivar_2d)=1 lpi(poslpivar_2d)=7 new=2
- `executives/2D/PIETOOLS_Hinf_gain_2D.m`: conv(mat2opvar)=2 help(clean_opvar)=1 lpi(checkdeg_lpi_eq_2d)=2 lpi(lpi_eq_2d)=2 lpi(lpi_ineq_2d)=2 lpi(lpigetsol)=4 lpi(poslpivar_2d)=7 new=2
- `executives/2D/PIETOOLS_Hinf_gain_2D_non_coercive.m`: conv(mat2opvar)=2 help(clean_opvar)=1 lpi(checkdeg_lpi_eq_2d)=2 lpi(lpi_eq)=1 lpi(lpi_eq_2d)=2 lpi(lpi_ineq_2d)=2 lpi(lpigetsol)=4 lpi(lpivar)=1 lpi(poslpivar_2d)=7 new=1
- `executives/2D/PIETOOLS_Hinf_gain_dual_2D.m`: conv(mat2opvar)=2 help(clean_opvar)=1 lpi(checkdeg_lpi_eq_2d)=2 lpi(lpi_eq_2d)=1 lpi(lpi_ineq_2d)=1 lpi(lpigetsol)=4 lpi(poslpivar_2d)=5 new=1
- `executives/2D/PIETOOLS_PDEstability_2D.m`: help(clean_opvar)=1 lpi(checkdeg_lpi_eq_2d)=2 lpi(lpi_eq_2d)=1 lpi(lpi_ineq_2d)=1 lpi(poslpivar_2d)=5 new=1
- `executives/2D/PIETOOLS_PDEstability_dual_2D.m`: help(clean_opvar)=1 lpi(checkdeg_lpi_eq_2d)=2 lpi(lpi_eq_2d)=1 lpi(lpi_ineq_2d)=1 lpi(poslpivar_2d)=5 new=1
- `executives/2D/PIETOOLS_PIE2PDEstability_2D.m`: conv(mat2opvar)=1 help(clean_opvar)=1 lpi(lpi_eq_2d)=2 lpi(lpi_ineq_2d)=2 lpi(lpivar_2d)=1 lpi(poslpivar_2d)=4
- `executives/2D/PIETOOLS_PIE2PDEstability_dual_2D.m`: conv(mat2opvar)=1 help(clean_opvar)=1 lpi(lpi_eq_2d)=2 lpi(lpi_ineq_2d)=2 lpi(lpivar_2d)=1 lpi(poslpivar_2d)=4

**executives/utility_functions** (6 files)

- `executives/utility_functions/closedLoopPIE.m`: conv(op2copvar_sop)=1 str=4
- `executives/utility_functions/getController.m`: help(clean_opvar)=1 str=1
- `executives/utility_functions/getObserver.m`: str=6
- `executives/utility_functions/getObserver_2D.m`: help(clean_opvar)=1 help(is_separable)=1 str=3
- `executives/utility_functions/get_eq_opts_2D.m`: fld2=48
- `executives/utility_functions/get_lpivar_degs.m`: fld=15 fld2=35 help(getdeg)=1 str=12

**executives_sopvar** (17 files)

- `executives_sopvar/PIETOOLS_H2_control_sop.m`: conv(op2copvar_sop)=6 lpi(lpi_ineq)=2
- `executives_sopvar/PIETOOLS_H2_estimator_sop.m`: conv(op2copvar_sop)=6 lpi(lpi_ineq)=2
- `executives_sopvar/PIETOOLS_H2_norm_sop.m`: conv(op2copvar_sop)=4 lpi(lpi_ineq)=2
- `executives_sopvar/PIETOOLS_Hinf_control_sop.m`: conv(op2copvar_sop)=7 lpi(lpi_ineq)=1
- `executives_sopvar/PIETOOLS_Hinf_estimator_sop.m`: conv(op2copvar_sop)=8 lpi(lpi_ineq)=1
- `executives_sopvar/PIETOOLS_Hinf_gain_sop.m`: conv(op2copvar_sop)=6 lpi(lpi_ineq)=1
- `executives_sopvar/PIETOOLS_stability_sop.m`: conv(op2copvar_sop)=2
- `executives_sopvar/PIETOOLS_well_posedness_sop.m`: conv(op2copvar_sop)=2
- `executives_sopvar/tests/test_executives_sop.m`: str=5
- `executives_sopvar/utility_functions/exec_tools_sop.m`: conv(copvar2opvar)=1 conv(copvar2opvar2d)=1 conv(mat2opvar)=2 conv(op2copvar_sop)=3 new=1
- `executives_sopvar/utility_functions/getController_direct_sop.m`: conv(opvar2copvar)=3 new=1 str=3
- `executives_sopvar/utility_functions/getController_sop.m`: conv(copvar2opvar)=1 conv(mat2opvar)=1 conv(opvar2copvar)=3 str=3
- `executives_sopvar/utility_functions/getObserver_sop.m`: conv(copvar2opvar)=1 conv(mat2opvar)=1 conv(opvar2copvar)=3 str=3
- `executives_sopvar/utility_functions/pie_disc_sop.m`: fld=3 str=2
- `executives_sopvar/utility_functions/pie_witness_sop.m`: conv(pie2opvar)=1
- `executives_sopvar/utility_functions/poslpivar_sop.m`: lpi(poslpivar)=1
- `executives_sopvar/utility_functions/trace_rn_sop.m`: fld2=1 str=6

**lpi_programming_sopvar** (19 files)

- `lpi_programming_sopvar/examples/hinf_gain_1d_sop.m`: conv(op2copvar_sop)=5
- `lpi_programming_sopvar/examples/stability_1d_sop.m`: conv(op2copvar_sop)=2
- `lpi_programming_sopvar/examples/stability_2d_sop.m`: conv(op2copvar_sop)=2
- `lpi_programming_sopvar/examples/volterra_norm_sop.m`: conv(op2copvar_sop)=1 fld=1 new=1
- `lpi_programming_sopvar/getsol_lpivar_sop.m`: lpi(getsol_lpivar)=1 str=8
- `lpi_programming_sopvar/lpi_eq_sop.m`: lpi(lpi_eq)=2
- `lpi_programming_sopvar/lpi_ineq_sop.m`: lpi(lpi_ineq)=2
- `lpi_programming_sopvar/lpigetsol_sop.m`: lpi(lpigetsol)=2 str=14
- `lpi_programming_sopvar/tests/bench_getsol_sop.m`: lpi(getsol_lpivar)=1 lpi(lpigetsol)=1 lpi(lpivar)=1 str=1
- `lpi_programming_sopvar/tests/e2e_hinf_gain_1d.m`: conv(opvar2copvar)=5 lpi(lpi_ineq)=1
- `lpi_programming_sopvar/tests/e2e_stability_2d.m`: conv(opvar2d2copvar)=3 new=1
- `lpi_programming_sopvar/tests/test_dpvar_hinf_chain_sop.m`: conv(opvar2copvar)=5 lpi(lpi_ineq)=1 lpi(lpigetsol)=1
- `lpi_programming_sopvar/tests/test_dpvar_ops_sop.m`: conv(mat2opvar)=3 conv(opvar2copvar)=12 help(rand_opvar)=7 new=1 str=2
- `lpi_programming_sopvar/tests/test_getsol_sop.m`: conv(opvar2copvar)=5 lpi(getsol_lpivar)=2 lpi(lpi_ineq)=1 lpi(lpigetsol)=9 lpi(lpivar)=2 lpi(poslpivar)=1 str=9
- `lpi_programming_sopvar/tests/test_lpi_eq_sop.m`: lpi(lpi_eq)=9 lpi(lpivar)=2 str=2
- `lpi_programming_sopvar/tests/test_lpi_ineq_sop.m`: lpi(lpi_ineq)=1 lpi(poslpivar)=1 str=1
- `lpi_programming_sopvar/tests/test_lpiprogram_sop.m`: lpi(lpi_ineq)=1
- `lpi_programming_sopvar/tests/test_tier1_fixes_sop.m`: fld=6 help(rand_opvar)=1 lpi(lpigetsol)=1 str=5
- `lpi_programming_sopvar/tests/test_translators_sop.m`: conv(opvar2copvar)=2 fld=1 lpi(lpivar)=1 lpi(poslpivar)=7 new=2

**sopvar classes** (4 files)

- `sopvar/@copvar/copvar2opvar.m`: conv(sopvar2opvar)=1 new=1 str=2
- `sopvar/@copvar/copvar2opvar2d.m`: conv(sopvar2opvar2d)=1 new=1 str=2
- `sopvar/@sopvar/sopvar2opvar.m`: conv(sopvar2opvar2d)=1 fld=3 new=1
- `sopvar/@sopvar/sopvar2opvar2d.m`: new=1

**sopvar/Testfolder** (12 files)

- `sopvar/Testfolder/Bench_copvar_vs_legacy.m`: conv(opvar2copvar)=1 conv(opvar2d2copvar)=2 conv(sopvar2opvar)=1 conv(sopvar2opvar2d)=1 help(rand_opvar)=4 help(rand_opvar2d)=2 new=1 str=2
- `sopvar/Testfolder/Test_copvar_converters.m`: conv(copvar2opvar)=2 conv(copvar2opvar2d)=4 conv(opvar2copvar)=2 conv(opvar2d2copvar)=3 conv(opvar2sopvar)=1 conv(sopvar2opvar)=2 conv(sopvar2opvar2d)=1 fld=3 help(rand_opvar)=4 help(rand_opvar2d)=2 new=2 str=1
- `sopvar/Testfolder/Test_copvar_inv.m`: conv(copvar2opvar)=5 conv(mat2opvar)=3 conv(opvar2copvar)=3 fld=29 help(getController)=1 help(getObserver)=1 help(inv_opvar_2)=1 new=9 str=4
- `sopvar/Testfolder/Test_sopvar_diff.m`: conv(opvar2sopvar)=1 fld=2 help(apply_opvar)=1 new=1 str=1
- `sopvar/Testfolder/sdopvar/Test_apply_sopvar_1d.m`: conv(opvar2sopvar)=1 help(apply_opvar)=1 help(rand_opvar)=1
- `sopvar/Testfolder/sdopvar/Test_apply_sopvar_2d.m`: conv(opvar2sopvar)=1 help(apply_opvar)=1 help(rand_opvar2d)=1
- `sopvar/Testfolder/sdopvar/Test_ctranspose_sopvar_2d.m`: conv(opvar2d2sopvar)=1 conv(sopvar2opvar2d)=1 help(clean_opvar)=1 help(rand_opvar2d)=1
- `sopvar/Testfolder/sdopvar/Test_mtimes_sopvar_1d.m`: conv(opvar2sopvar)=2 conv(sopvar2opvar)=1 fld=3 help(clean_opvar)=1 help(rand_opvar)=2
- `sopvar/Testfolder/sdopvar/Test_mtimes_sopvar_2d.m`: conv(opvar2d2sopvar)=2 conv(sopvar2opvar2d)=1 help(clean_opvar)=1 help(rand_opvar2d)=2
- `sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_PIE2PDEstability.m`: conv(opvar2copvar)=2
- `sopvar/Testfolder/sdopvar/claude_tests/cx_exec/cx_stability_2D.m`: conv(opvar2d2copvar)=3 new=1
- `sopvar/Testfolder/sdopvar/test_sopvar_AT_bug.m`: conv(sopvar2opvar)=3 fld=3 help(clean_opvar)=1

**ndopvar** (14 files)

- `ndopvar/@ndopvar/blkdiag.m`: str=1
- `ndopvar/@ndopvar/subsasgn.m`: str=1
- `ndopvar/@ndopvar/vertcat.m`: str=1
- `ndopvar/@nopvar/subsasgn.m`: str=1
- `ndopvar/@nopvar/vertcat.m`: str=1
- `ndopvar/Testfolder/converters/dopvar2d2ndopvar.m`: conv(opvar2dopvar2d)=1 fld2=10 str=3
- `ndopvar/Testfolder/converters/dopvar2ndopvar.m`: conv(dopvar2d2ndopvar)=3 fld=3 str=5
- `ndopvar/Testfolder/converters/ndopvar2dopvar.m`: fld=3 fld2=2 new=4
- `ndopvar/Testfolder/ndopvar/Test_functions_nopvar_ndopvar_AT.m`: conv(dopvar2d2ndopvar)=24 help(rand_opvar2d)=4 lpi(lpivar)=3
- `ndopvar/Testfolder/ndopvar/Test_mtimes_ndopvar.m`: conv(dopvar2ndopvar)=4 conv(ndopvar2dopvar)=2 help(rand_opvar)=2 lpi(lpivar)=1
- `ndopvar/Testfolder/ndopvar/Test_mtimes_ndopvar2d.m`: conv(dopvar2d2ndopvar)=2 conv(dopvar2ndopvar)=2 conv(ndopvar2dopvar)=2 help(rand_opvar2d)=2 lpi(lpivar)=1
- `ndopvar/Testfolder/ndopvar/ctranspose_test.m`: conv(dopvar2d2ndopvar)=2 fld2=2 help(rand_opvar2d)=1
- `ndopvar/Testfolder/ndopvar/rand_ndopvar.m`: new=1
- `ndopvar/Testfolder/rand_opvar2d.m`: fld2=28 new=1

**PIESIM** (7 files)

- `PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m`: fld2=64
- `PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_4PI2Mat_cheb.m`: fld=5
- `PIESIM/PIESIM_routines/PIESIM.m`: conv(pie2opvar)=3
- `PIESIM/PIESIM_routines/PIESIM_Setup/PIESIM_input_check.m`: str=1
- `PIESIM/PIESIM_routines/PIESIM_Setup/transl_2D.m`: fld2=96
- `PIESIM/PIESIM_routines/PIESIM_Setup/transl_mod.m`: fld=12 str=2
- `PIESIM/PIESIM_routines/PIESIM_Setup/weightPIE.m`: fld=4 new=1

**PIESOS** (20 files)

- `PIESOS/@intop/apply_op.m`: help(int_simple)=1
- `PIESOS/@tensopmat/initialize.m`: conv(opvar2sopvar)=1 str=2
- `PIESOS/@tensopmat/tensopmat.m`: conv(opvar2sopvar)=2 str=4
- `PIESOS/@tensopvar/initialize.m`: conv(opvar2sopvar)=1 str=2
- `PIESOS/Examples/Linear_Examples/PIE2polyPIE.m`: conv(opvar2sopvar)=2 conv(pie2opvar)=1 help(clean_opvar)=2 str=3
- `PIESOS/Examples/Linear_Examples/Test_Conversion_Linear.m`: help(apply_opvar)=6
- `PIESOS/Examples/LocalStability/SOS_DP.m`: conv(opvar2sopvar)=3 fld=9 new=3
- `PIESOS/Examples/LocalStability/Testfolder/Test_innerprod_v2.m`: conv(opvar2sopvar)=1 fld=3 help(apply_opvar)=1 new=1
- `PIESOS/Examples/LocalStability/Testfolder/Test_polyopvar_sosquadvar.m`: conv(opvar2sopvar)=1 fld=3 new=1
- `PIESOS/Examples/LocalStability/Testfolder/Test_polyopvar_times_v2.m`: conv(opvar2sopvar)=2 fld=6 new=2
- `PIESOS/Examples/LocalStability/Testfolder/Test_quad2lin_v2.m`: conv(opvar2sopvar)=1 fld=3 help(apply_opvar)=1 new=1
- `PIESOS/Examples/LocalStability/Tests/DMB_test.m`: conv(opvar2sopvar)=1 fld=3 help(apply_opvar)=1 new=1
- `PIESOS/Examples/LocalStability/Tests/V3_Liediff_test.m`: conv(opvar2sopvar)=1 fld=6 new=1
- `PIESOS/Examples/LocalStability/Tests/WSB_test.m`: conv(opvar2sopvar)=1 conv(sopvar2opvar)=1
- `PIESOS/Examples/LocalStability/V3_Liediff.m`: conv(opvar2sopvar)=1 fld=6 new=1
- `PIESOS/Testfolder/Test_is_id_op.m`: conv(sopvar2opvar)=1 fld=5
- `PIESOS/Testfolder/polyopvar/Test_polyopvar_innerprod.m`: help(int_simple)=1
- `PIESOS/Testfolder/tensopvar/Test_tensopvar_int.m`: help(int_simple)=1
- `PIESOS/sosprogramming/piesos_getsol.m`: lpi(getsol_lpivar)=1 str=5
- `PIESOS/sosprogramming/piesos_ineq.m`: lpi(lpi_ineq)=1 str=4

**parser** (13 files)

- `parser/@state/int.m`: fld=2 new=1
- `parser/@state/minus.m`: fld=1 new=1
- `parser/@state/mtimes.m`: fld=2 new=1
- `parser/@state/plus.m`: fld=1 new=1
- `parser/@state/state2terms.m`: fld=1 new=1
- `parser/@state/uminus.m`: fld=1 new=1
- `parser/@state/uplus.m`: fld=1 new=1
- `parser/@sys/getPDEparams.m`: fld=24
- `parser/@terms/int.m`: fld=6 new=1
- `parser/@terms/mtimes.m`: fld=2 new=1
- `parser/@terms/plus.m`: fld=1
- `parser/@terms/terms.m`: cls=1
- `parser/@terms/vertcat.m`: new=1

**DEMOs (PIETOOLS_demos root)** (3 files)

- `PIETOOLS_demos/DEMO11_CdopvarControl.m`: conv(pie2copvar)=1
- `PIETOOLS_demos/DEMO2_volterra_operator_norm.m`: conv(op2copvar_sop)=1 fld=1 new=1
- `PIETOOLS_demos/DEMO7_observer_based_control.m`: conv(op2copvar_sop)=2

**snippets_from_manual** (3 files)

- `PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch3_PI_Operators.m`: fld=7 fld2=12 new=6
- `PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch5_PIEs.m`: fld=11 new=6
- `PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch7_LPI_Programming.m`: help(getObserver)=1 lpi(lpi_eq)=1 lpi(lpi_ineq)=1 lpi(lpigetsol)=3 lpi(lpivar)=1 lpi(poslpivar)=3

**sopvar_demos** (1 file)

- `PIETOOLS_demos/sopvar_demos/test_heatNd_pie.m`: fld=3 fld2=1 str=1

**multivariate_PIE demo** (1 file)

- `PIETOOLS_demos/multivariate_PIE_2508_14840/jp_wave_plate_rates.m`: lpi(lpi_eq)=1 lpi(lpi_ineq)=2 lpi(lpivar)=2

**root scripts** (1 file)

- `PIETOOLS_PDE.m`: fld=1 fld2=4

**baselines (untracked)** (3 files)

- `baselines/2026-10-28_example_library/private/hd_build.m (untracked)`: fld2=2
- `baselines/2026-10-28_example_library/private/hd_qform.m (untracked)`: conv(op2copvar_sop)=5 lpi(lpi_ineq)=1
- `baselines/2026-10-28_example_library/private/lb_fp.m (untracked)`: fld=3

## Appendix B – Scripts, logs and data [M]

All in the scratchpad folder `C:\Users\mpeet\AppData\Local\Temp\claude\C--Users-mpeet-OneDrive-Documents-GitHub-PIETOOLS-PIETOOLS\05ad1b64-cee9-4875-8809-890e152518e3\scratchpad\opvar_inv`:

- `scan.pl` (scanner), `agg.pl` (subsystem table), `idx.pl` (callers of a list of names), `appx.pl` (Appendix A);
  `hits.tsv` (raw hits, 17,581 rows), `hits_final.tsv` (hits used for the counts).
- `m/chk1.m` … `m/chk4.m` (MATLAB checks of §1.2); `gen/` (the four copied LPI functions of session 2).
- `out/session1_log.txt`, `session2_log.txt`, `session3a_log.txt`, `session3b_log.txt`, `session4_log.txt`,
  `session4_log_2.txt`, `session4_rows.txt` (one row per PDE); `out/session1.mat` (the PIEs and PIESIM
  solutions of session 1).
- `manual.txt`, `manual_root.txt` (pdftotext of the two manual copies); `gui_doc.xml` (the GUI code).
