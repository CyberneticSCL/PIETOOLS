# lpi_programming_sopvar — LPI programming for both operator families

Initial coding MMP, 09/29/2026 (Tier 1 of the container parity map).
MMP, 10/06/2026: the 1-D executive translators (Tier 2), see §2.

This folder sits beside `lpi_programming/` and lets one LPI body serve both operator
families:

- the legacy classes: `opvar`, `dopvar`, `opvar2d`, `dopvar2d`, `dpvar`, `polynomial`
  and `double`;
- the container family: `sopvar`, `sdopvar`, `copvar` and `cdopvar`, in any number of
  spatial variables.

Every `_sop` function sends a legacy class to the existing legacy routine, unchanged,
so its result is the legacy result. It sends a container-family class to container
code. No file in `lpi_programming/` was modified.

## 1. Naming rule

`pietools_path_update` adds every folder to the path with one `genpath`, so a file here
must never reuse a legacy name. Such a file would shadow the legacy one for every
caller. Every function here therefore ends in **`_sop`**: `lpiprogram_sop`,
`lpi_eq_sop`, `lpigetsol_sop`, and so on. Helpers follow the same rule:

- `tests/private/same_val_sop` and `examples/private/opcheck_sop` are reachable only
  from their own folders.

MMP, 09/30/2026: `private/spaces2meta_sop` (a verbatim copy of the space parser of
`lpivar_cdopvar`) is deleted. The constructors call the shared
`parse_copvar_spaces` in `sopvar/misc/conventions/` instead, a name that does not
collide with a legacy one.

`which -all` resolves each function to one file in the tree.

## 2. Dispatch tables

### `prog = lpiprogram_sop(...)`

This is `lpiprogram` with three more input forms (the last three rows below). Since
10/01/2026 `lpiprogram` itself takes any number of variables, so `lpiprogram_sop` only
converts those forms and calls it; before, it was a copy of lpiprogram without the
2-variable cap. Its outputs and error messages were checked unchanged on 93 input forms.

| input | result |
|---|---|
| every `lpiprogram` form: `(vartab,dumvartab,dom[,decvartab[,freevartab]])`, `(vartab,dom,...)`, an n×2 `vartab`, a 1×2 `dom` for all variables | as `lpiprogram`, for any n |
| `(names,dom,...)` with `names` a cellstr or string, e.g. a container's registry `P.vars`, `P.dom` | polynomial variables in that order; dummies `s_dum` |
| `(P[,decvartab[,freevartab]])` with `P` a `copvar`, `cdopvar`, `sopvar` or `sdopvar` | `P`'s sorted registry and domains; dummies `s_dum` |
| `(struct('vars',..,'dom',..),...)` | the same |

It has these fields: `vartable = [vars; dummies; free]` (a `polynomial`), `dom` (n×2),
`decvartable` (a cellstr), and the SOSTOOLS fields.

The former copy deviated from `lpiprogram` in three places. Since 10/01/2026 `lpiprogram`
has two of them itself (it checks the dummy variables with `ispvar`, and with no spatial
variable sets `dom = zeros(0,2)`); the cellstr and container forms stay here.

### `prog = lpi_eq_sop(prog,P[,opts])`

| class of `P` | routine |
|---|---|
| `cdopvar`, `copvar` | `lpi_eq_cdopvar` |
| `sdopvar`, `sopvar` | `lpi_eq_sdopvar`, which collects and imposes the rows (MMP, 09/30/2026: it has one output and no internal modes since then) |
| any other class (`dpvar`, `dopvar`, `dopvar2d`, `opvar`, `opvar2d`, `polynomial`, `double`, ...) | `lpi_eq`, unchanged |

`opts` (`'symmetric'`) is passed on only when given. Each routine keeps its own errors,
for example for a fixed `copvar`/`sopvar`, a polynomial, or an unknown option. The
program it returns equals the routine's own program field by field (test_lpi_eq_sop).

### `sol = lpigetsol_sop(prog,X)`, `Psol = getsol_lpivar_sop(prog,P)`, `Xsol = subs_dvar_sop(X,names,vals)`

These were written by the getsol agent. They are verified again here in the chain and
through E1–E4.

| class | `lpigetsol_sop` | `getsol_lpivar_sop` |
|---|---|---|
| `sdopvar` → `sopvar`, `cdopvar` → `copvar` | new code | new code |
| `sopvar`, `copvar` | returned as is | returned as is |
| a cell holding any of the four classes | taken apart, one element at a time | — |
| `[]`, `double`, `polynomial`, `opvar`, `opvar2d`, `dpvar`, char, another cell, `dopvar`, `dopvar2d` | `lpigetsol`, unchanged | `opvar`…`dopvar2d`: `getsol_lpivar` |
| anything else | error `lpigetsol_sop:badClass` | error `getsol_lpivar_sop:badClass` |

Values are read from `prog.solinfo.RRx`, never from `solinfo.x`: RRx is in decvartable
order, and x is the solver's cone vector.

- Each operator's decision list `Zd` is looked up **once** per container.
- The solved parameter is `A_γ + B_γ'd` for each γ cell. No `dpvar` of q names is
  formed.
- A name that is listed twice takes its value from its lowest row, the row where
  getequation places it.
- `subs_dvar_sop` is the same code applied to a given (names, values) pair, for
  certifying a candidate point.

### Constructors (Tier 1c)

The arguments `dims`, `spaces` and `dom` follow the conventions of `lpivar_cdopvar` and
`poscopvar`:

- `dims` is M×1, or `struct('out',..,'in',..)`;
- `spaces` is a 1×M cell of cellstr, where `{}` is R^q, or `struct('out',..,'in',..)`;
- `dom` is nv×2 in the order of the sorted registry, one 1×2 row for all variables, or
  `struct('vars',..,'dom',..)`.

| function | result |
|---|---|
| `mat2copvar_sop(Mat,dims,spaces,dom[,opts])` | the container of a constant matrix: `copvar` for a `double`, `cdopvar` for a `dpvar`. A block acts as a multiplier in shared variables, is constant in output-only variables and is integrated over input-only ones. `opts.mult_only` refuses an integral block, as `mat2opvar` does. |
| `eye_copvar_sop(dims,spaces,dom)` | the identity. Square only (`eye_copvar_sop:notSquare`). |
| `zeros_copvar_sop(dims,spaces,dom)` | the zero operator. Rectangular is allowed. |

All three pass `verify` and satisfy the canonical multiplier form. They replace the
`opvar2copvar(mat2opvar(...))` route, which has no N-D form.

### dpvar as an operator (Tier 1b, in the class files)

These changes are in the class files, not in this folder: `sopvar/@*/mtimes`, `plus`,
`horzcat` and `vertcat`, plus `sopvar/misc/dpvar_op_copvar.m`, `mat2copvar_grid.m` and
`dpvar2sdvar.m`. Each class method has one branch at the top of its body, and the
composition code below it is untouched.

They follow the legacy `opvar`/`dopvar` semantics:

- `gam*X` and `X*gam` for a scalar dpvar;
- `D*X` for a dpvar matrix;
- `gam ± X`, read as `gam*I ± X`;
- a dpvar entry inside `[ ]`, e.g. `[-gam, D'; D, -gam]`.

A product of two decision operators is refused (`cdopvar:decisionTimesDecision`).
Numeric operands keep their earlier behaviour: `P + 1e-3` and `[X, 0]` are still
refused. Use `dpvar(c)` or `mat2copvar_sop` instead.

### Legacy routines that need no `_sop` version

`lpidecvar`, `lpisetobj`, `lpisolve` and the scalar `lpi_ineq(prog,gam)` (sosineq) work
unchanged on container programs. **Measured** in two places:

- E1: 1-D, gamma declared, `gam >= 0`, and minimized;
- test_lpiprogram_sop (d): N = 1, 2, 3, 4 programs from `lpiprogram_sop`, with
  `lpidecvar`, `lpi_ineq`, `lpisetobj`, `lpisolve` and `lpigetsol_sop`.

Neither `prog.vartable` nor `prog.dom` is read by lpidecvar, lpisetobj, lpisolve or any
container routine.

### Executive translators (Tier 2, MMP, 10/06/2026)

The settings-to-operator steps the 1-D executives share. Each sends a legacy class to
the stock routine unchanged; the container branch mirrors it.

| function | legacy class | container |
|---|---|---|
| `poslpivar_sop(prog,n/X,d,options[,side,dom])` | `poslpivar` | `poscopvar` over X's `side` spaces, with poslpivar's degree filling, psatz (weight only for 1), exclude and sep translated; 1-D |
| `poslpivar_settings_sop(prog,n/X,settings,role[,side,dom])` | the `dd1`/`dd12` (`'lf'`) or `dd2`/`dd3` (`'slack'`) pair of `poslpivar` calls | the same pair through `poslpivar_sop`; P1 first, `P1 + P2`; imposes nothing |
| `get_lpivar_degs_sop(R[,T])` | `get_lpivar_degs` | the same 1-D rule read off the nonzero coefficients |
| `trace_rn_sop(X)` | `trace(X.P)`, `trace(X.R00)` | trace of the R^q diagonal blocks, a 1 x 1 dpvar over `X.Zd` |
| `copvar_space_list(X,side)` (sopvar/misc/conventions) | — | the inverse of `parse_copvar_spaces` |

They replace twelve test-folder translators of `cx_exec`, deleted 10/06/2026:
`cx_pl2pm`, `cx_hinf_posdeg`, `cx_h2_pos`, `cx_stability_pos`, `cx_hinf_lf`,
`cx_h2_lf`, `cx_hinf_slack`, `cx_h2_slack`, `cx_hinf_qdeg`, `cx_stability_lpivar_degs`,
`cx_hinf_dimsp`, `cx_h2_qdeg` (all MMP 09/25/2026). The three poslpivar translations
among them disagreed on psatz other than 0/1, on sep and on degree filling;
`poslpivar_sop` follows poslpivar on all three. `cx_space_list`, `cx_h2_trace` and
`cx_on_registry` stay unchanged for the frozen 2-D helpers.

**Measured:** the 17 1-D transcriptions at all six presets (126 programs; the 5 at
`extreme` fail identically, Qdeg(3) = -1) are bit-identical to before. Generator spans
equal poslpivar's on 16 cases, sep included (test_poscopvar_vs_poslpivar, part 2); the
container Q degrees equal stock `get_lpivar_degs` on 5 plants x 6 presets x 3 operator
forms.

## 3. Tier 1 scope

**Covered:**

- getsol for `sdopvar`/`cdopvar` (1a);
- dpvar operators (1b);
- constructors (1c);
- `lpiprogram_sop`, `lpi_eq_sop`, `lpigetsol_sop` and `getsol_lpivar_sop` (1d).

**Not covered:**

- **Container `lpi_ineq` is Tier 2.** A positive operator is declared with
  `poscopvar`/`copquadvar` and cancelled by `lpi_eq_sop`, as E1–E4 do. A dpvar scalar
  inequality goes through the legacy `lpi_ineq`.
- No `_sop` positive/indefinite variable dispatch: `lpivar`/`poslpivar` versus
  `lpivar_cdopvar`/`poscopvar` are still chosen by the caller. (MMP, 10/06/2026:
  `poslpivar_sop` now dispatches the positive variable, 1-D; there is still no
  `lpivar_sop`.)
- The PDE→PIE conversion and `pie_struct` still stop at 2 variables. heatNd builds its
  N-D PIE directly.
- There is no regrid, no container `inv`, no closed-loop PIE, no PIESIM.
- A dpvar that depends on a spatial variable is refused (`mat2copvar_grid:spatial`).
  Polynomial operands, `blkdiag`/`eq` with a dpvar, and numeric operands in `[ ]`/`+`
  are not handled.
- The 2-D executives, the cx_* 2-D helpers and `settings2possopvar` are frozen by the
  maintainer and were not touched. E3 uses `cx_stability_2D` unmodified as its
  reference.

## 4. Cost (measured unless marked)

**lpiprogram_sop.** O(n) in the spatial variables. It takes 2.90 ms against 2.91 ms for
`lpiprogram` at N = 2, and 1.7 / 3.2 / 3.2 / 3.7 ms at N = 1 / 3 / 5 / 9 (min of 20).

**lpi_eq_sop.** It adds 2–4 class tests to the routine's own cost, and passes its
arguments through with no copy. Measured with `tests/bench_lpi_eq_sop`: P = 0 for an
`lpivar_cdopvar` operator of q decision variables, one row per coefficient, min of 3,
alternating with the direct call.

| q | N = 1 | N = 2 | N = 3 |
|---|---|---|---|
| 1e4 | 22.9 / 17.8 ms | 23.7 / 17.1 ms | 19.7 / 16.5 ms |
| 1e5 | 0.248 / 0.269 s | 0.297 / 0.276 s | 0.234 / 0.203 s |
| 1e6 | 5.07 / 4.97 s | 5.29 / 5.34 s | 5.14 / 5.46 s |

Each cell is lpi_eq_sop / lpi_eq_cdopvar. The two programs are equal at every point.

- At q = 1e6 the ratio is 0.94–1.02. At 1e4 the dispatcher is 3–7 ms slower; it always
  ran first in each pair (*inferred*: warm-up order, since the branch is 2–4 `isa`
  calls).
- Profiler PeakMem at q = 1e6: 35.3 / 30.5 MB (N = 1), 30.6 / 30.6 MB (N = 2),
  31.3 / 31.3 MB (N = 3).
- The time is linear in q and flat in N.

**getsol, i9-14900KF** (numbers from the getsol agent's benchmark). The time is
dominated by the one name lookup, O(N log q). No q-sized dense array is formed.

- q = 3e6: 2.95 s. The time is flat over 1, 2 and 3 spatial variables (2.92–2.98 s).
- Peak memory is about 9 B per table entry: 77.5 MB at q = 3e6.
- 3-D heat (E4), 2.7e6 decision variables: the 9 extractions took 10.8 s in total,
  including one P*T composition.

**dpvar operators** (numbers from the operators agent):

| operation | cost |
|---|---|
| `gam*X` | flat in q |
| `gam + Xd` at q = 1e6 | 1.3–1.7 s over 1–3 variables: one `unique` in `merge_dvar_lists`, O(q log q) |
| hot path | the 5 assembly programs are bit-identical; 36 branch predicates per 2-D H∞ build cost ≤ 20 µs in total |

## 5. Running

Set the path with `run pietools_path_update`, then remove the `.claude`/`worktrees`
entries. A saved pathdef or a nested worktree can shadow the tree. The examples need
the cx_exec helpers (`sopvar/Testfolder/sdopvar/claude_tests/cx_exec`) and
`heatNd_apply` (`PIETOOLS_demos/sopvar_demos`), which `genpath` includes. The solvers
are MOSEK 11 and SeDuMi (E3).

| command | what | time |
|---|---|---|
| `test_lpiprogram_sop` | (a) 16 lpiprogram forms field by field and 8 error messages; (b) N = 3, 4 against the hand-built N-D program; (c) the 7 new input forms; (d) N-D solves (Volterra N = 1, 2; multiplier N = 1..4) | 92 s |
| `test_lpi_eq_sop` | dispatch of 16 cases field by field, the option pass-through, 7 errors | 2 s |
| `test_getsol_sop`, `bench_getsol_sop` | getsol (getsol agent) | 18 s, minutes |
| `test_dpvar_ops_sop`, `test_dpvar_hinf_chain_sop` | dpvar operators and constructors (operators agent) | 17 s, 10 s |
| `test_endtoend_sop` | E1, E2 by default; `'all'` adds E3 and E4 | 15 s; E3 5 min; E4 15–20 min |
| `examples/volterra_norm_sop` | DEMO2 on containers | 5 s |
| `test_translators_sop` | the executive translators: legacy dispatch, errors, the settings pair, Q degrees against stock, trace against its definition, space lists against opvar dims (MMP, 10/06/2026) | 9 s |

## 6. End-to-end results

MOSEK 11 and SeDuMi 1.3, R2025b, 09/29/2026. The machine was shared with other jobs,
so times may be high. A certified verdict needs all three of: MOSEK status, rel_b ≤ 1e-6
on the solinfo.RRx rows, and PSD blocks.

**E1 — `hinf_gain_1d_sop`.** 1-D H∞ gain, io1 plant, light settings, gamma a decision
variable, one SDP.

- The solve is certified: numerr 0, rel_b 4.88e-7, psd_relmin 2.4e-11.
- **γ = 0.18262633**, against the stock `PIETOOLS_Hinf_gain` value 0.18262634 (relative
  difference 4.7e-8). It lies in the container bisection bracket [0.182503, 0.182626]
  to within the bracket's 6-digit rounding.
- The shape is the stock one: 1274 decision variables, Ks [31 13 10 1], m 237 against
  242.
- Extracted against the rows: ‖coef(T'Q − R)‖/‖r‖ = 1.000000, and 1.4133 for K + N,
  which lies in [1, √2].
- By action:

| check | value |
|---|---|
| T'Q = R, weak form | 3.8e-7 |
| K + N | 1.0e-6 |
| getsol(K) against K rebuilt from γ and Q | 8.2e-17 |
| R ≥ 0, N ≥ 0, K ≤ 0 on test functions | min 0.68 / 0.43 / max −0.43 |
| control (solinfo.x) | 0.71 and 1.13 |

**E2 — `stability_1d_sop`.** 1-D stability, rd plant, 0.5·λ*.

- The solve is certified: rel_b 1.33e-8.
- Rows: 1.000000 for T'Q − P and 1.347 for D + N. The rebuilt row sets have exactly the
  program's row norm.
- By action: T'Q = P (weak form) 2.1e-8; D = A'Q + Q'A (weak form) 9e-17; D + N 4.9e-9.
- P and N ≥ 0 on the test functions.
- The control (solinfo.x) gives 0.97 and 0.87.

**E3 — `stability_2d_sop`.** 2-D heat, rd2 plant, 0.1·λ*, light settings, SeDuMi.
**This checks extraction only; it is not a certificate.**

- The transcription builds cx_stability_2D's program bit for bit: m 3456, 179840
  decision variables, Ks [424 8].
- SeDuMi: numerr 1, feasratio 0.925, rel_b 3.2e-6 (above the gate, as expected), in
  215–273 s.
- ‖coef(Qe + Q)‖/‖r‖ = 1.4142136, which lies in [1, √2]: the extracted operators
  reproduce the solver's own row residual exactly.
- By action:
  - getsol(Q) against Q(getsol(P)): 3e-14;
  - Qe + Q: 8.5e-5;
  - P ≥ 0 and Qe ≥ 0 on the test functions.
- The control is RRx plus 1e-3 rms noise, which gives 3.5e-2. solinfo.x cannot serve as
  a control here: it equals RRx, because the program has no free variable and no
  inequality.

**E4 — `heat3d_stability_sop`.** 3-D heat (DD×DN×DN), κ = 14.0 against κ* = 14.804407,
bench preset, d = 0. The SDP has 2,713,563 decision variables, m 31770, and Ks 288×7 plus
552×7.

- The run reused the row set of the stored κ = 14.60 dump of the same family, 22226
  rows. The check found m, K and b equal.
- **Certified +1**: MOSEK OPTIMAL, rel_b 5.81e-7 on all rows, psd_relmin 5.8e-14. MOSEK
  took 852 s and 20 iterations; the process peak was 12.97 GB.
- The build took 133 s, and extracting all 9 operators (P, R, Q, the constraint
  operators and P*T) took 10.8 s.
- Rows against coefficients: 1.4142136 for (15a) and for (15b). The two row sets
  reproduce ‖At'RRx − b‖ exactly (9.89e-9).
- By action (3-D quadrature, 57 s):
  - P*T − ε²T*T = R: 2.7e-8;
  - X1 + κX2 = −Q: 2.7e-9;
  - R ≥ 0, Q ≥ 0 and −(X1 + κX2) ≥ 0 on the test functions (0.74, 0.61, 0.61
    normalized).
- P*T = (P*T)' holds with residual exactly 0. At d = 0 this is structural, so it is not
  a discriminating check.
- The control (RRx plus 1e-3 noise) gives 2.4e-3 and 1.2e-3.
- Total run time was 18 min.

**test_lpiprogram_sop (d).** All solves use MOSEK and pass the certification gate.

| problem | result |
|---|---|
| Volterra N = 1 | ‖V‖ ≤ 0.63664911 (exact 2/π = 0.63661977) |
| Volterra N = 2, 4 linear face terms | 0.40528613 (exact 0.40528473) |
| multiplier [2 1; 1 2], N = 1–4 | γ = 3.00000000; N = 4 has 26245 decision variables |

The coefficient/row ratios lie in [1.23, 1.414].

## 7. Regression chain (09/29/2026, working tree, nothing committed)

| suite | result |
|---|---|
| `sg_battery` | 1 error, 0 failures; the error is the known `test_leftshift_monomials_sopvar`, which calls a class-private function (moved to sopvar/private/dead_code/ later on 09/29/2026 with the dead helpers it tested) |
| `cx_run_1d` (called directly) | 21 case lines identical to the operators agent's pre-change reference (timings excluded) |
| `test_heatNd_pie`, `test_heatNd_lpi` (with `test_copquadvar_faces`), `test_heatNd_poincare` | 70, 74 (97) and 44 checks pass |
| `test_copvar_silent_fixes` | 643 checks |
| `test_getsol_sop` | parts i–iv |
| `test_dpvar_ops_sop` | 592 passed, 0 failed |
| `test_dpvar_hinf_chain_sop` | 7 passed |
| `volterra_norm_sop` | passes |
| `test_lpi_eq_sop` | 32 checks |
| `test_lpiprogram_sop` | parts a–d |
| `test_endtoend_sop` | E1, E2 and E3 pass; E4 passes (run separately, with the stored row set) |
