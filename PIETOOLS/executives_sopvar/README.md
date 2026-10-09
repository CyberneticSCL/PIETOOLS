# executives_sopvar: the executives on the container path (MMP, 10/08/2026)

Parallel versions of every file in `executives/` (1-D) and `executives/2D/`, one `PIETOOLS_<name>_sop`
per stock `PIETOOLS_<name>`, built on `lpi_programming_sopvar` and `sopvar/lpis_sopvar`. The stock files
are untouched; `pietools_path_update` is one `genpath`, so every name here ends in `_sop` and
`which -all` resolves each once. The stock output list is kept, in the stock order, and one `info`
struct is appended.

| stock executive | container version | builder (private/) | form |
|---|---|---|---|
| PIETOOLS_PDEstability, _dual | PIETOOLS_PDEstability_sop, _dual_sop | stability_build_sop | P, Pd |
| PIETOOLS_PIE2PDEstability, _dual | PIETOOLS_PIE2PDEstability_sop, _dual_sop | stability_build_sop | Q, Qd |
| PIETOOLS_Hinf_gain, _coercive, _dual, _dual_coercive | ..._sop | hinf_build_sop | Q, P, Qd, Pd |
| PIETOOLS_Hinf_control, _estimator | ..._sop | synth_build_sop | ctrl, est |
| PIETOOLS_H2_norm_c, _o, _c_coercive, _o_coercive | ..._sop | h2_build_sop | c, o, cco, oco |
| PIETOOLS_H2_control, _estimator | ..._sop | h2_build_sop | ctrl, est |
| PIETOOLS_well_posedness | PIETOOLS_well_posedness_sop | (in the file) | |
| PIETOOLS_auto_execute (script) | PIETOOLS_auto_execute_sop | | |
| 2D/PIETOOLS_stability_2D, _dual_2D | ..._sop | stability_build_sop | P, Pd |
| 2D/PIETOOLS_Hinf_gain_2D, _2D_non_coercive, _dual_2D | ..._sop | hinf_build_sop | P, Q, Pd |
| 2D/PIETOOLS_Hinf_estimator_2D | PIETOOLS_Hinf_estimator_2D_sop | synth_build_sop | est |
| 2D/PIETOOLS_H2_norm_2D_c, _o, _c_non_coercive, _o_non_coercive | ..._sop | h2_build_sop | cco, oco, c, o |

The 1-D wrappers route a 2-D PIE to the 2-D version with the same LPI, as the stock wrappers do
(`Hinf_gain` to `_2D_non_coercive`, `Hinf_gain_coercive` to `_2D`, both duals to `_dual_2D`,
`PDEstability` and `PIE2PDEstability` to `stability_2D`, the H2 norms to their 2-D forms).

## What is different from the stock executives

1. **Operators.** The PIE operators are converted once to containers (`op2copvar_sop`: `opvar2copvar`
   in 1-D, a per-component `opvar2d2sopvar` grid in 2-D); every product, adjoint and sum is a container
   operation. The LPIs are those of the stock files as the 10/06/2026 transcriptions
   (`sopvar/Testfolder/sdopvar/claude_tests/cx_exec`) wrote them, which produced programs identical to the
   stock on the test plants.
2. **Storage.** Declared as the stock executive declares it, through `exec_tools_sop`: in 1-D
   `poslpivar_settings_sop` 'lf' (dd1, options1, dd12, options12 of the settings), in 2-D
   `poslpivar_settings_2d_sop` 'lf' (LF_deg, LF_opts and the psatz terms of `settings_2d`, translated by
   `settings2possopvar`); eppos added where the stock adds it (`exec_tools_sop` eye: eppos, eppos2 in
   1-D, the four-entry eppos of `settings_2d` in 2-D). The user's level (light, heavy, ...) keeps its
   meaning.
3. **Negativity.** Every `Dop <= 0` is `lpi_ineq_sop(prog,-Dop,opts)`: the slack is sized by
   `get_lift_degs` on the storage operator (`opts.like`), with the weight raised by `settings.sop.dw_Q`
   (default 1) in the Q forms and `settings.sop.dw_P` (default 2) in the coercive forms (measured,
   `lpi_programming_sopvar/README.md`), the 1-D Markov-Lukacs pair or the N-D faces one degree below the
   plain term, through one `poscopvar_direct` call. `settings.sop.slack = 'stock'` restores the stock
   slack (`poslpivar_settings_sop` 'slack' in 1-D; in 2-D `poslpivar_settings_2d_sop` 'slack', eq_deg and
   eq_opts with the exclusions and separations `get_eq_opts_2D` derives from the operator) for
   comparison.
4. **Degree loop.** The program is built by a closure of the raise `dw`; `lpi_solve_loop_sop` solves,
   classifies the exit (`lpi_classify_sop`), and on MOSEK status UNKNOWN (the fixed-degree alternative of
   the duality note: no certificate at any finite value) or on primal infeasibility (a slack too small
   for the operator) raises `dw` by one and rebuilds, up to `settings.sop.max_raise` (default 2) times.
   For a feasibility LPI (stability, well-posedness) the raises cost solves and cannot change an
   infeasible verdict; set `max_raise = 0` to skip them.
5. **Dual.** The builders set `prog.sopeq = {}` after `lpiprogram` (when `settings.sop.keepdual`), which
   asks `lpi_eq_cdopvar` to record which row of the program is which coefficient of which cell (opt-in:
   programs built elsewhere gain no field); `lpigetdual_sop` assembles the solver's multipliers into the
   signed dual kernel of the negativity constraint (`info.dual`, a `copvar`; `info.dual_norm1`). The
   multipliers are in the original row coordinates with `sos_opts.simplify` on or off (sossolve restores
   removed rows with 0; checked: b'y equals the primal objective both ways). `pie_witness_sop` evaluates the numerical
   counterpart of the certificate on the Galerkin discretization of the 1-D PIE in the orthonormal
   Legendre basis (`pie_disc_sop`, degree `settings.sop.N_cheb`, exact pairings by quadrature): the
   frequency response for a gain (a numerical lower bound, `info.gain_lb`, at `info.omega`, and the
   relative gap `info.gap = gam/gain_lb - 1`), the spectrum of the pencil (A,T) for a rate
   (`info.maxre`), the Lyapunov H2 norm (`info.h2_num`), and the alignment of the worst input with the
   dual kernel (`info.dual_alignment`, 1 when the solver's dual points at the worst input). Not
   available in 2-D (`info.witness.ok` false).
6. **Outputs.** The stock outputs in the stock order, then `info`: status ('optimal', 'inaccurate',
   'infeasible', 'unbounded', 'unknown'), certified, the value, the witness fields, the degrees and
   terms of the slack, the SDP shape (`lpi_shape_sop`: decision count, cone sizes, rows), the loop
   history and the dual kernel. The synthesis executives return K or L as an opvar (opvar2d in 2-D)
   through the stock `getController` / `getObserver` / `getObserver_2D`.
7. **Margin.** `settings.sop.margin = eps` certifies `-Dop + eps*e_D >= 0` with `e_D` the lift's own
   energy (the regularizer of the proof program under which polynomial certificates exist at a fixed
   lift). Default 0.
8. **Fixed gamma.** `settings.sop.gam_fixed = g` (or the third argument of the 2-D gain and estimator
   executives, the stock's `gain`) poses the feasibility test at gamma = g; the value returned is g when
   certified and Inf otherwise. The stock bisection options are not reproduced.

## Settings

Stock `lpisettings` structs work unchanged; `exec_sop_settings` adds the `sop` field with defaults
(slack 'like', dw_Q 1, dw_P 2, max_raise 2, margin 0, witness true, N_cheb 24, nfreq 240, keepdual true,
verbose true) and selects MOSEK when installed; `sos_opts.simplify` is left as given. In 2-D the free
operators (Q of the Q forms, Z of the
estimator) take their degrees from the storage (`get_lpivar_degs_sop`, the per-role reading) and from
the largest entry of `settings_2d.Zop_deg` (one cap on every role of `lpivar_cdopvar`, a superset of the
stock `lpivar_2d` monomials); in 1-D from `get_lpivar_degs_sop` and `settings.ddZ`, as the stock.

## Library files added for this folder (lpi_programming_sopvar)

`exec_sop_settings`, `lpi_classify_sop`, `lpi_solve_loop_sop`, `lpi_shape_sop`, `lpigetdual_sop`,
`pie_disc_sop`, `pie_witness_sop`, `op2copvar_sop`, `on_registry_sop`, `settings_2d_sop`,
`poslpivar_settings_2d_sop`; `get_lpivar_degs_sop` gained the N-D reading;
`lpi_ineq_sop` the margin option and the equality tag; `lpi_eq_sop` returns the tag;
`sopvar/@copvar/copvar2opvar` converts non-square grids. `sopvar/lpis_sopvar/lpi_eq_cdopvar` and
`private/collect_eq_rows` record the row bookkeeping the dual read-back uses.

## Status (measured 10/08/2026, MOSEK 11, light settings)

`tests/test_executives_sop.m` asserts the 1-D part (19 checks; `test_executives_sop(false)` adds the
2-D stability executive). The numbers below are from the session drivers (`run_exec_1d`,
`run_exec_1d_b`, `run_exec_2d`, `run_exec_2d_b`); 'stock' is the stock executive, 'stock slack' this
folder with `settings.sop.slack = 'stock'`, 'like' the default. "numerical" is `pie_witness_sop`.

**Structural note.** Every dual-form LPI (Hinf_gain_dual, Hinf_gain_dual_coercive, Hinf_control,
H2_control, H2_norm_c_coercive) is infeasible or UNKNOWN for a distributed disturbance in both modes
and in the stock: the Lyapunov block (AP)T' + T(AP)' is compact (T' is purely integral) and the Schur
complement asks it to dominate B B' = I on L2. H2_norm_o and _o_coercive take the trace over the R^nw
block, empty for a distributed w, and return the solver floor, as the stock. Test those executives on a
finite-dimensional disturbance; the 1-D table below uses one (reaction-diffusion, w through s(1-s), u
constant, y = int x, z = [int x; u]; numerical gain 0.0333106, numerical H2 0.0523601).

| executive (1-D) | stock | stock slack | like | note |
|---|---|---|---|---|
| PDEstability, _dual, PIE2PDEstability, _dual (rd 0.5) | certified | certified | certified | max Re eig(A,T) -4.9348 = -pi^2/2 |
| well_posedness (rd 0.5) | certified | certified | certified | |
| Hinf_gain (io1) | 0.18262633 | 0.18262633 | 0.18248166 | numerical 0.18247908, gap 1.4e-5, dual alignment 0 (distributed w) |
| Hinf_gain | 0.033332351 | 0.033332353 | 0.033310666 | gap 1.5e-6, dual alignment 1.000 |
| Hinf_gain_coercive | UNKNOWN 0.0574 | UNKNOWN 0.0575 | 0.033311035 | stock light slack short in the lift |
| Hinf_gain_dual | 0.03333238 | 0.033332403 | 0.033310842 | |
| Hinf_gain_dual_coercive | UNKNOWN 0.0584 | UNKNOWN 0.0583 | 0.033313355 | |
| H2_norm_c (io1) | 0.28878602 | 0.28879091 | 0.28775623 | numerical 0.28768384, gap 2.5e-4 |
| H2_norm_c | 0.052715641 | 0.052715636 | 0.052428236 | gap 1.3e-3 |
| H2_norm_o | 0.052715396 | 0.052715396 | 0.052428398 | |
| H2_norm_c_coercive | 0.0825005 | 0.0824501 | 0.052396837 | stock 57% loose; gap 7e-4 |
| H2_norm_o_coercive | 0.0827630 | UNKNOWN 0.0823 | 0.052393641 | |
| Hinf_control | UNKNOWN 0.0569 | UNKNOWN 0.0569 | 0.032865781 | K opvar; open loop 0.0333 |
| Hinf_estimator | 2.11e-4 | 3.18e-4 | 1.39e-4 | at the eppos floor |
| H2_control | UNKNOWN 0.0829 | UNKNOWN 0.0829 | 0.052222382 | K opvar; open loop 0.0524 |
| H2_estimator | 1.71e-4 | 2.40e-4 | 1.63e-4 | at the eppos floor |
| auto_execute (io1: stability, Hinf_gain, H2_norm_dual, well_posedness) | | | runs | the stock names plus info_* |

| executive (2-D) | stock | stock slack | like | note |
|---|---|---|---|---|
| stability_2D (rd2 0.5) | UNKNOWN (23 s) | UNKNOWN (21 s) | certified (11 s) | like: nx 187984, m 3075, faces [2 2] at w [2 2] |
| stability_dual_2D (rd2 0.5) | UNKNOWN (23 s) | UNKNOWN (17 s) | certified (8 s) | |
| Hinf_gain_2D (io2) | UNKNOWN 1.2722 (413 s) | not run | 0.17172212 (26-33 s) | nx 200426, m 4096; matches the 10/08 KYP measurement |
| Hinf_gain_2D_non_coercive (io2) | | | UNKNOWN 277 (15 s) | the 2-D Q form, as in the 10/08 sizing study |
| Hinf_gain_dual_2D (io2) | | | INFEASIBLE (3 s) | distributed w: structural |
| H2_norm_2D_c (io2) | | | INFEASIBLE (23 s) | distributed w: structural |
| H2_norm_2D_o (io2) | | | 8.6e-5 (19 s) | empty R^nw trace: the floor |
| H2_norm_2D_c_non_coercive (io2) | | | INFEASIBLE at dw 0 (6 s) | nx 303234; raises disabled in the run |
| H2_norm_2D_o_non_coercive (io2) | | | UNKNOWN 11.0 at dw 0 (95 s) | nx 624536: W on the L2 input space through the LF storage |
| Hinf_estimator_2D (io2 + y = int x) | | | 5.1e-4 (24 s) | L opvar2d; nx 200451 |

The 2-D rows are 'like' only with `max_raise = 0` and no witness; of the stock 2-D gain and H2 programs
only Hinf_gain_2D was run to completion (UNKNOWN after 413 s); the 'stock' slack mode rebuilds those
programs and was not run for them within the daytime budget.
