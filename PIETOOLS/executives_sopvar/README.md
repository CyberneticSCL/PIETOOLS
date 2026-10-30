# executives_sopvar: the executives on the container path

Nine executives, one per analysis, built on `lpi_programming_sopvar` and the helpers in
`utility_functions/` (below). Each takes a `pie_struct` (1-D or 2-D, or a container PIE, below) and a settings structure, writes
its LPI in its own file, and returns the outputs of the stock executive followed by one `info` struct.
The LPI form is an argument; the dimension is read from the PIE. Every name ends in `_sop`, so
`which -all` resolves each once beside the stock executives, which are untouched.

**Container PIEs and 3 or more variables (since 10/29/2026).** The PIE may instead be a container
PIE: a `pie_struct` whose operators are `copvar` objects, as `pie2copvar` returns, in any number nv
of spatial variables (vars nv x 2, dom nv x 2).
- `exec_pie_sop` initializes an opvar PIE as before. For a container PIE it sets every absent
  operator ([]) to 0.
- In 1-D and 2-D a container PIE builds the same SDP as the opvar PIE: 23 executive/form builds
  matched, with identical fingerprints.
- In 3 or more variables, `exec_tools_sop` declares the storage from `settings.settings_nd`. The
  default follows the stock light 2-D storage: a multiplier and one separable full integral per
  space. Only the 'like' slack is available, and solved operators stay containers.
- Stability, H-infinity gain and H2 norm run in N-D. The H-infinity estimator solves its LPI but
  returns no gain, because the gain is extracted in 1-D and 2-D only. The controllers, the H2
  estimator and well-posedness reject nv >= 2.
- Size: the 'like' slack is a tensor product of the 1-D lift basis, 12^3 = 1,728 per term over 7
  terms in 3-D at weight 1 (m 145,123 with the first all-basis storage), beyond MOSEK on this
  workstation. In 3 or more variables the tiers bound it per direction (`cap_delta`, see
  `settings_sopvar/README.md`): 360 per term at light.
- Measured 10/29/2026 at light (MOSEK 11.0.30, this workstation, one pass each) on heat fields on
  [0,1]^3, s1 Dirichlet-Dirichlet and s2, s3 Dirichlet-Neumann, against exact modal values. One
  field: u_t = lap u + w, z = int u. ODE + field: x_t = -x + int u, u_t = lap u + x + w, z = x.

  | plant | LPI | rows | build | solve | result |
  |---|---|---|---|---|---|
  | one field | stability | 28,956 | 34 s | 313 s | certified |
  | one field | H-infinity gain | 31,154 | 25 s | 564 s | 0.0585827, exact 0.0414732 (+41%) |
  | one field, u_t = lap u + 5u | stability | 33,231 | 47 s | 312 s | certified (r = 5, 34% of the threshold; 10/30/2026) |
  | ODE + field | stability | 29,155 | 91 s | 353 s | certified |
  | ODE + field | H-infinity gain | 31,493 | 26 s | 442 s | 0.0605787, exact 0.0432677 (+40%) |
  | two coupled fields | stability, gain | 166,572, 170,772 | 280 s, 351 s | not run | beyond MOSEK here |

  Rows grow about as the square of the number of PDE fields. MOSEK memory grows about as rows^1.9-2
  and peaked at 13.7 GB on the r = 5 program (10/30/2026), so two coupled fields project to about
  290-340 GB. At heavy the one-field programs build in 163-167 s with 58,922 and 59,718 rows, and
  their solves project to about 41-44 GB (not run); `settings_sopvar/README.md` has the tier table.

| executive | FORM | stock executives it covers |
|---|---|---|
| `PIETOOLS_stability_sop(PIE,settings,form)` | 'P' (default), 'Pd', 'Q', 'Qd', and the finite-energy forms 'FP', 'FQ' (since 10/28/2026; epneg > 0); in 2-D 'Q', 'Qd', 'FQ' are posed as 'P', 'Pd', 'FP' | PDEstability, _dual, PIE2PDEstability, _dual, FEPDEstability, FEPIE2PDEstability, 2D/stability_2D, _dual_2D, 2D/FEPDEstability_2D, FEPIE2PDEstability_2D |
| `PIETOOLS_well_posedness_sop(PIE,settings)` | | well_posedness |
| `PIETOOLS_Hinf_gain_sop(PIE,settings,form,gain)` | 'Q' (default), 'Qd', 'P', 'Pd'; in 2-D 'Qd' is posed as 'Pd' | Hinf_gain, _coercive, _dual, _dual_coercive, 2D/Hinf_gain_2D, _2D_non_coercive, _dual_2D |
| `PIETOOLS_Hinf_control_sop(PIE,settings)` | | Hinf_control |
| `PIETOOLS_Hinf_estimator_sop(PIE,settings,gain)` | | Hinf_estimator, 2D/Hinf_estimator_2D |
| `PIETOOLS_H2_norm_sop(PIE,settings,form)` | 'c' (default), 'o', 'cco', 'oco' | H2_norm_c, _o, _c_coercive, _o_coercive, 2D/H2_norm_2D_c_non_coercive, _o_non_coercive, _2D_c, _2D_o |
| `PIETOOLS_H2_control_sop(PIE,settings)` | | H2_control |
| `PIETOOLS_H2_estimator_sop(PIE,settings)` | | H2_estimator |
| `PIETOOLS_auto_execute_sop` (script) | | auto_execute |

Without settings, each executive takes the default of the stock executive it covers: `lpisettings('light')`,
except stability in 1-D (`lpisettings('heavy')`, eppos 1e-4, eppos2 1e-6) and in 2-D (`lpisettings('light')`,
simplify 1, `settings_2d.eppos` 1e-2), and well-posedness (`lpisettings('heavy')`, simplify 1, eppos = eppos2 =
1e-2). `lpisettings_sop` (folder `settings_sopvar`) gives tiers tuned for this path. `lpiscript_sop`
(folder `lpi_programming_sopvar`, since 10/28/2026) runs them by the names of `lpiscript`.

On 10/28/2026 these nine replaced the 28 files of 10/08/2026 (one per stock executive, each calling a
builder in `private/`). Each form poses the program its former file posed: 31 cases compared, with status,
value, decision count, rows and Gram sizes identical.

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
   `lpi_programming_sopvar/NOTES.md`), the 1-D Markov-Lukacs pair or the N-D faces one degree below the
   plain term, through one `poscopvar_direct` call. `settings.sop.slack = 'stock'` restores the stock
   slack (`poslpivar_settings_sop` 'slack' in 1-D; in 2-D `poslpivar_settings_2d_sop` 'slack', eq_deg and
   eq_opts with the exclusions and separations `get_eq_opts_2D` derives from the operator) for
   comparison.
4. **Degree loop.** The program is built by a closure of the raise `dw`; `lpi_solve_loop_sop` solves,
   classifies the exit (`lpi_classify_sop`), and on MOSEK status UNKNOWN (the fixed-degree alternative of
   the duality note: no certificate at any finite value) or on primal infeasibility (a slack too small
   for the operator) raises `dw` by one and rebuilds, up to `settings.sop.max_raise` (default 2) times.
   An infeasible verdict is a statement about the slack degree, not about the LPI: on the 2-D stability
   executive the slack at `dw_P = 0` is infeasible (6.5e3 variables, 0.46 s) and one raise certifies
   (1.04e5 variables), 2.95 s in total against 6.35 s from the default `dw_P = 2` (measured 10/16/2026).
   What a raise cannot do is show that no degree certifies; set `max_raise = 0` to skip the raises.
5. **Dual.** The executives set `prog.sopeq = {}` after `lpiprogram` (when `settings.sop.keepdual`), which
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
   through `getController_sop` / `getObserver_sop` in 1-D (since 10/09/2026; `utility_functions/`,
   below) and the stock `getObserver_2D` in 2-D.
7. **Margin.** `settings.sop.margin = eps` certifies `-Dop + eps*e_D >= 0` with `e_D` the lift's own
   energy (the regularizer of the proof program under which polynomial certificates exist at a fixed
   lift). Default 0.
8. **Fixed gamma.** `settings.sop.gam_fixed = g` (or the argument `gain` of `PIETOOLS_Hinf_gain_sop` and
   `PIETOOLS_Hinf_estimator_sop`, the stock 2-D executives' `gain`) poses the feasibility test at gamma = g; the value returned is g when
   certified and Inf otherwise. The stock bisection options are not reproduced.

## Settings

Stock `lpisettings` structs work unchanged; `exec_sop_settings` adds the `sop` field with defaults
(slack 'like', dw_Q 1, dw_P 2, max_raise 2, margin 0, witness true, N_cheb 24, nfreq 240, keepdual true,
verbose true; `dw_feas`, `max_raise_feas`, `eppos_stab` empty, which leaves stability and well-posedness
on `dw_P`/`dw_Q`, `max_raise` and the stock eppos) and selects MOSEK when installed; `sos_opts.simplify` is
left as given. In 2-D the free
operators (Q of the Q forms, Z of the
estimator) take their degrees from the storage (`get_lpivar_degs_sop`, the per-role reading) and from
the largest entry of `settings_2d.Zop_deg` (one cap on every role of `lpivar_cdopvar`, a superset of the
stock `lpivar_2d` monomials); in 1-D from `get_lpivar_degs_sop` and `settings.ddZ`, as the stock.

## utility_functions/ (moved from lpi_programming_sopvar, 10/16/2026)

The helpers the executives call. They are on the path like every folder: the tests and examples of
`lpi_programming_sopvar`, the cx_exec tests, `sopvar/@sopvar/inv` and the stock `getController` call some
of them. The stock counterparts of `get_lpivar_degs_sop`, `getController_sop` and `getObserver_sop` are
in `executives/utility_functions/`.

| function | what it does |
|---|---|
| `exec_sop_settings` | completes a settings struct with the `sop` field (Settings, above) |
| `exec_tools_sop` | the storage, stock slack, margins, identities and opvar conversion the settings define for a 1-D or 2-D PIE, and for a container PIE in 3 or more variables (`settings_nd`, since 10/29/2026) |
| `exec_pie_sop` | the PIE an executive receives: an opvar PIE initialized, a container PIE with absent operators set to 0 (10/29/2026) |
| `exec_info_sop` | the part of `info` every executive shares: status, loop history, slack, SDP shape, dual kernel |
| `lpi_solve_loop_sop` | the degree loop: build, solve, classify the exit, raise `dw` |
| `settings_2d_sop` | the 2-D settings the stock 2-D executives read, with their defaults |
| `poslpivar_sop`, `poslpivar_settings_sop` | the 1-D storage and slack from the stock degrees and settings (Executive translators, below) |
| `poslpivar_settings_2d_sop` | the same in 2-D, for R^n, L2[x], L2[y] and L2[x,y] spaces (L2[x], L2[y] since 10/28/2026, from `dx`, `dy`) |
| `settings2possopvar` | a stock settings struct (`settings_PIETOOLS_*`, `*_2D`) as `possopvar` degrees and options, for `poslpivar_settings_2d_sop` and the cx_exec 2-D tests (from `sopvar/lpis_sopvar`, 10/16/2026) |
| `get_lpivar_degs_sop` | the degrees of the free operator, read off the storage |
| `trace_rn_sop` | the trace over the R^q blocks, for the H2 executives |
| `getController_sop`, `getObserver_sop` | the 1-D gains through the container inverse (below) |
| `getController_direct_sop` | the controller gains without forming the inverse; called by `Test_copvar_inv` only |
| `pie_disc_sop`, `pie_witness_sop` | the Legendre-Galerkin matrices of a 1-D PIE, and the numerical counterpart of a certificate |

Added for this folder on 10/08/2026 and kept in the library: `lpi_classify_sop`, `lpi_shape_sop`,
`lpigetdual_sop` (in `lpi_programming_sopvar/misc` since 10/21/2026), `op2copvar_sop` (now
`converters/operator_converters/`), `on_registry_sop` (now `sopvar/misc/`); `get_lpivar_degs_sop` gained
the N-D reading; `lpi_ineq_sop` the margin option and the equality tag; `lpi_eq_sop` returns the tag;
`sopvar/@copvar/copvar2opvar` converts non-square grids. `lpi_eq_cdopvar` and `private/collect_eq_rows`
(in `lpi_programming_sopvar/misc` since 10/21/2026; `sopvar/lpis_sopvar` until 10/16/2026) record the
row bookkeeping the dual read-back uses.

The two sections below moved here from `lpi_programming_sopvar/README.md` on 10/16/2026, verbatim.

### Operator inverse and gain reconstruction (MMP, 10/09/2026)

1-D only. The class methods are in `sopvar/@sopvar/inv.m` and `sopvar/@copvar/inv.m`;
the wrappers are in this folder. Measured on `sopvar/Testfolder/Test_copvar_inv.m` (17
checks, 10/09/2026): operator residuals `max|P*Pinv - I|` of 1e-9 to 1e-11 on seven
hand-built 3-PI and 4-PI operators (0.05 to 0.3 s each), and on the synthesis executives
`|K P - Z|/|Z|` = 2.4e-11 and `|P L - Z|/|Z|` = 7.0e-9 where the stock `getController`
and `getObserver` give 8.1e-5 and 7.4e-3 on the same solved P, Z.

| function | result |
|---|---|
| `[Rinv,info] = inv(R[,opts])` on a `sopvar` block L2^m[s] -> L2^m[s] | Gohberg-Krein (arXiv 2208.13104, Lemma 16, Cor. 17) on the stored coefficient matrices: the kernels `R_i = (I kron ZL') C_i (I kron ZR)` are split by one SVD each, `U` and `V` by RK4 on `opts.N` nodes (101), the three parameters fitted by Chebyshev least squares at one degree raised from `opts.deg0` (4) until the relative RMS residual is below `opts.tol` (1e-8) or `opts.degmax` (16). `info`: ranks, `rcondU22`, `condR0`, `d`, `relrms`. |
| `[Pinv,info] = inv(P[,opts])` on a `copvar` over R^k x L2^m[s] | the L2 block by the method above, the rest by the block inverse through the finite-dimensional Schur complement `T = Pm - Q1 Rh Q2` (Lemma 18): one operator inverse and one matrix inverse. `info.condT`. N-D is refused. |
| `[K,Kop,info] = getController_sop(P,Z[,opts])` | `K = Z P^{-1}` as a `copvar`; `Kop` its `opvar` for `closedLoopPIE`, `piess`, PIESIM. Accepts `opvar` inputs (`opvar2copvar`). No coefficient truncation (the stock routine zeroes coefficients of K below 1e-4). |
| `[L,Lop,info] = getObserver_sop(P,Z[,opts])` | `L = P^{-1} Z`, the same way. |
| `[K,info] = getController_direct_sop(P,Z[,opts])` | the gains `K1`, `K2(s)` of `u = Z P^{-1} x` WITHOUT forming `P^{-1}` (the analogue of Cor. 11 of arXiv 1806.08071): grid values of the Gohberg-Krein data and their cumulative integrals (`gk_grid`), Simpson for the Schur complement; `K.apply` evaluates `u` by quadrature, `K.op` is the one polynomial fit (of `K2`) for `closedLoopPIE`, `K.fit.relrms` its residual. Agrees with the inverse route to 1e-9 (Test_copvar_inv). The QT formulation is discussed in `opvar/inverse_dependency_map_2026_10_09.md`, Sec. 8.2. |

`PIETOOLS_Hinf_control_sop`, `_Hinf_estimator_sop`, `_H2_control_sop` and `_H2_estimator_sop` call
these in 1-D and return the `opvar` gain; `info.G` is the container gain and `info.inv` the inverse diagnostics. 2-D keeps
`getObserver_2D`. The defects of the stock routines are documented in
`opvar/inverse_dependency_map_2026_10_09.md`.

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

## Status (measured 10/08/2026, MOSEK 11, light settings)

The rows name the executives of 10/08/2026; the table at the top gives the executive and form that pose
each of these programs now.

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
