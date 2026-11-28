# Phase 1b – 2025 executives and 2025 inputs (E1, E2), 11/28/2026

Maintainer, 11/28/2026: E1 "Error, but with specific guidance on the correct calling format for
the case which led to the specific error"; E2 the D6 rule (stop with an error that gives the
current call). Committed 0d02ac08.

## Change

- `executives/retired_2025/` (new): one file for each of the 21 executives of PIETOOLS 2025 whose
  LPI is now a form of one executive, written by `make_stubs.pl` (its table is the record). A call
  raises `<name>:retired` with the call that poses the LPI in this release; the help gives the
  2025 call and the call now. The forms: the LPI in the 2025 headers (`A'*P*T + T'*P*A` for
  `PDEstability`, `R = T'*Q` for `PIE2PDEstability` and the 2-D `_non_coercive` names), the map of
  `lpiscript` (2025 and now: `'stability'` → `PIE2PDEstability` = `'Q'`, `'l2gain-dual'` →
  `Hinf_gain_dual` = `'Qd'`, `'h2norm'` → `H2_norm_o` = `'o'`), the identification of 11/27
  (`Hinf_gain_dual_coercive` = `'Pd'`, `H2_norm_c_coercive` = `'cco'`, `H2_norm_o_coercive` =
  `'oco'`), and the table of `bl_sop.m` for the 2-D names. Notes in the message where 2025 differs:
  the 2025 `OPTIONS.h2` (a fixed H2 bound) has no counterpart; a 2-D name says that the executive
  takes a PIE in any number of variables; `GAIN` keeps its 2025 meaning.

| 2025 | now |
|---|---|
| `PIETOOLS_PDEstability`, `_dual` | `PIETOOLS_stability` `'P'`, `'Pd'` |
| `PIETOOLS_PIE2PDEstability`, `_dual` | `PIETOOLS_stability` `'Q'`, `'Qd'` |
| `PIETOOLS_Hinf_gain_coercive`, `_dual`, `_dual_coercive` | `PIETOOLS_Hinf_gain` `'P'`, `'Qd'`, `'Pd'` |
| `PIETOOLS_H2_norm_o`, `_c`, `_o_coercive`, `_c_coercive` | `PIETOOLS_H2_norm` `'o'`, `'c'`, `'oco'`, `'cco'` |
| `PIETOOLS_stability_2D`, `_dual_2D` | `PIETOOLS_stability` `'P'`, `'Pd'` |
| `PIETOOLS_Hinf_gain_2D`, `_2D_non_coercive`, `_dual_2D` | `PIETOOLS_Hinf_gain` `'P'`, `'Q'`, `'Pd'` (with `gain`) |
| `PIETOOLS_H2_norm_2D_c`, `_c_non_coercive`, `_o`, `_o_non_coercive` | `PIETOOLS_H2_norm` `'cco'`, `'c'`, `'oco'`, `'o'` |
| `PIETOOLS_Hinf_estimator_2D` | `PIETOOLS_Hinf_estimator` (with `gain`) |

- `getController`, `getObserver`: a non-empty numeric third input (the 2025 `TOL`, which `isfield`
  ignored) raises `getController:call2025` / `getObserver:call2025`; `[]` and a `struct` pass.
- `lpi_ineq`: `opts.pure` or `opts.sep` (2025 options, not read; `popts` is built anew,
  `lpi_ineq.m` near line 400) raise `lpi_ineq:call2025`. No live caller passes either field (the
  `sep` options of the toolbox go to `poslpivar`). The help states the reading of a numeric
  `psatz` (below).
- One line of help in each of the three headers; they are new in this release: no log line, no
  markers (D1). The stubs and the test are new files.

## Measured

`test_retired2025` (`testing_codes/private/executives/tests/`, log `test_retired2025.log`): 27
checks pass – each stub resolves once in `executives/retired_2025/` and raises its identifier
with the call of the table; the guards of `getController`, `getObserver`, `lpi_ineq`.

A numeric `opts.psatz` of `lpi_ineq` is a valid current call (term codes), so it cannot be caught
by type. On `T'*T` of the heat equation on [0,1]:

| `opts.psatz` | codes | Gram | ndec |
|---|---|---|---|
| `0` | `0` (plain) | 18 | 324 |
| `1` | `1` (box product alone) | 18 | 324 |
| `[0 1]` | `[0 1]` | [18 18] | 648 |
| `'none'` | `0` | 18 | 324 |
| `'auto'` | `[0 1]`, box product one degree lower | [18 12] | 468 |

The 2025 `psatz = 0` keeps its meaning; the 2025 `psatz = 1` (the plain and the psatz term) is
`'auto'` (or `[0 1]`) here, and `1` poses the box product without the plain term.

`checkcode` (`cc_1b.m`, `cc_1b.log`): no message and no mtree error in the 25 files; none before in
the three edited files.
