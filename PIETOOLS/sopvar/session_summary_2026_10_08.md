# Session summary, 10/07 to 10/08/2026: separated form, direct map, container inequality, slack sizing, the dual

Maintainer: M. Peet. Branch `ndopvar`. Everything below is measured on the workstation (i9-14900KF, MOSEK 11)
unless marked as inferred. Commit bc3be32f (10/08) holds the work up to the H-infinity slack sizing; the
examination of the dual, the tensor term set and the heat test are uncommitted on top of it (see Sec. 9).

## 1. What was built

| piece | file | what it does |
|---|---|---|
| separated form | `sopvar/lpis_sopvar/lift_copvar.m`, `posmult_cdopvar.m`, `poscopvar_lift.m` | the positive operator as `A'*M*A`: a lift `A` of basis operators and a polynomial matrix weight `M`; span-identical to `copquadvar` on 18 cases |
| direct map | `sopvar/lpis_sopvar/poscopvar_direct.m` | the positive `cdopvar` straight from the Gram blocks, no intermediate operator; identical to `copquadvar` on 24 cases (joint and subset caps, faces, 3-D); 'fast' path with Gram-side shifts for the box weights and a position-map cache; several Positivstellensatz terms in one call; single-direction quadratic codes `2nv+2+d` (10/08, late) |
| degree reader | `sopvar/lpis_sopvar/get_lift_degs.m` | lift `D = Dmin + dD`, weight `w = max(ceil((Dmax-D-1)/2), ceil(Mdeg/2), J-D-wR-1, 0) + dw` per direction from the lower-kernel support of the target; option `like` reads the support of another operator and lays the specification out over the target's spaces |
| container inequality | `lpi_programming_sopvar/lpi_ineq_sop.m` | `P >= 0` for `cdopvar`/`copvar` through the reader, one `poscopvar_direct` call with every term, and `lpi_eq_sop`; legacy classes go to `lpi_ineq` unchanged |
| tests | `lpi_programming_sopvar/tests/test_lpi_ineq_sop.m` (6 checks), `claude_tests/test_poscopvar_direct.m` (6 parts), `test_poscopvar_lift.m`, `test_mincap_dual.m` (2 checks) | all passing on 10/08 |
| benchmarks and harnesses | `PIETOOLS_demos/sopvar_demos/heatNd_tailor.m`, `heatNd_lpi.m` ('custom', `psatz_offset`, 'tensor'), `lpi_programming_sopvar/examples/hinf_tailor_1d.m` (forms Q, Qd, P, Pd), `hinf_tailor_2d.m` (forms Q, P), `volterra_tailor.m` | degree and term sweeps with certified bisection or MOSEK gains |
| minimum-cap dual | `claude_tests/mincap_dual.m` | the proof program's dual hierarchy at one level, a small SDP |
| notes | `sopvar/sopvar_lift_notes.tex`, `.pdf` (18 pages) | Secs. 1 to 6 the form, 7 degrees and every measurement, 8 cost, 9 the direct map, 10 the optimization problem and its dual, 11 what remains |

## 2. Cost of the direct map (3-D plain term, measured 10/07)

| route | time cold | time warm | peak memory |
|---|---|---|---|
| `copquadvar` | 1 | | 2303 MB |
| `poscopvar_direct` | 1.7 to 8x faster | 18 to 37x faster | 465 MB |

The multiplier padding of `copquadvar` (every block stored on the reachable bases, 64x column padding in 3-D)
was the real cost; pruning it in `posmult_cdopvar` gave 14x in time and 11x in memory and made 3-D finish.

## 3. Heat benchmark, degree tailoring (certified bisection, eps 0.1)

1-D Dirichlet, kappa* = 9.869604: `heavy` certifies 9.8695954 at nx 1854 to 1866 and m 87 to 96 for P
degrees 0 to 2; R at (D,w) = (2,1), Q at (2,2) and the product term at w-1 certify the same rate at nx 820
to 832 and m 78 to 87: 2.25 times fewer variables and 10% fewer rows. Q binds and needs weight 2, one above
the support rule, which the reader's default dw = 1 supplies. Dirichlet-Neumann: `bench` already exact and
minimal.

2-D DD x DN, kappa* = 12.337006, d = 0: `bench` 12.183764 at nx 40409; `heavy` 12.336731 at nx 565209,
m 4282. The product pair certifies nothing; faces are needed. R at the bench degrees with Q at the heavy
degrees, faces at the same degree, certifies heavy's 12.336731 at nx 392409, m 2884 (31% fewer variables).
Uncapped faces at w-1 with Q at (2,2) reach 12.335223 at nx 339209 to 360009. The tensor term set is being
tested in this configuration now (Sec. 8).

## 4. H-infinity gain of io1, 1-D: four forms, slack sizing

Plant io1: x_t = x_ss + (pi^2/2) x + w, Dirichlet, z = int x. Closed-form gain 0.182479080 (attained at
frequency 0). Storage at the stock light degrees unless stated; entries are the excess over the closed form,
nx in parentheses.

| slack | Q form (primal) | dual Q form | coercive P form |
|---|---|---|---|
| stock light | 1.5e-4 (1274) | 7.2e-5 (1274) | 1.0e-2, UNKNOWN (1231) |
| stock heavy, heavy storage | 2.5e-6 (2377) | 8.9e-8 (2377) | 6.6e-7 (2332) |
| reader on K | (5,3) 2.5e-6 (5242) | (5,5) 1.1e-7 (12410) | (4,3) 1.4e-6 (3282) |
| reader on K, dw 0 | (5,2) 2.6e-6 (2834) | (5,4) 1.0e-7 (8434) | (4,2) 1.1e-6 (1786) |
| reader on the storage (`like`) | (2,2) 2.6e-6 (1058) | (2,2) 2.7e-6 (1058) | (2,2) 1.1e-5 (810) |
| (2,3) | 2.6e-6 | 1.5e-7 (1858) | 1.6e-6 (1426) |
| (1,1) | 1.1e-2 | 1.1e-2 | diverges |
| (1,2) | 1.4e-4 | 3.6e-4 | diverges |
| (2,1) | 1.7e-3 | 1.2e-3 | 1.9e-1 |
| (2,2) cap 3 | 2.6e-6 (898) | | 4.1e-3 (670) |

Findings:
- In the Q forms the reader on K overshoots because K = A'Q + Q'A is composed before T'Q = R is imposed
  and carries the support of the whole lpivar family of Q (lower kernel to degree 9, Dmin 4). The reader on
  R gives the measured minimum (2,2). The stock light slack is one short in the lift.
- The 2.5e-6 floor belongs to the primal Q form alone: the dual Q form and the coercive form reach 1e-7.
- In the coercive form the reader on K is right (P is positive at its declared degrees), lift-1 slacks make
  the program infeasible under the margin, the stock light slack is 5.7% off with a solver failure, and the
  joint cap costs 4e-3.
- The dual coercive executive is unusable on io1 as shipped: the stock returns 8252, every tailored slack is
  primal infeasible. Not diagnosed.
- Rule: Q forms, `like` with the default raise; coercive form, `like` with the weight raised by 2.

## 5. H-infinity gain of io2, 2-D

Q form: no tailored slack up to (2,2) with faces admits a finite gain (MOSEK UNKNOWN, objective diverging);
the stock light slack is a 1.5M-variable program. Coercive form at light storage:

| slack | gain | excess | nx | rows | solve |
|---|---|---|---|---|---|
| reader on P, (2,1), four faces | 0.2665 | +56% | 1.1e5 | 2788 | 5 s |
| (3,1), four faces | 0.1974 | +15% | 3.4e5 | 5513 | 25 s |
| (2,2), four faces | 0.171718 | +0.36% | 557845 | 5420 | 57 s |
| (2,2), tensor set | 0.171720 | +0.36% | 233280 | 4116 | 17.6 s |
| (2,2), tensor set, cap 6 | 0.171775 | +0.39% | 197720 | 3801 | 16.4 s |
| (3,3), product pair | 0.171777 | +0.40% | 1427142 | 13659 | 491 s |
| (2,2), (3,2), (2,3), product pair | UNKNOWN | | | | |
| stock light, degrees + 3 | stopped at 459 s | | 755227 | 13978 | |

The baseline suite (09/24) had the stock coercive executive on this plant at 1.272 with status UNKNOWN. The
closed-form gain is 0.1711.

## 6. The optimization problem and its dual (notes Sec. 10)

- Rows: MOSEK's presolve found no dependent row in any of the twelve 2-D solves; the storage
  parametrisation of the null weights, (N Y N')', would remove exactly the dependent rows: nothing.
- The Gram has (w+1)^2/(2w+1) times the coefficients of the polynomial weight it certifies, 1.8 at w = 2;
  that is the certificate, not slack.
- The term set is the one lever. The Markov-Lukacs tensor set (plain at w, each box quadratic at w-1 in
  its own direction, their product at w-1; `poscopvar_direct` codes `2nv+2+d`, in 1-D equal to the product
  code) gives the faces' gain at 42% of the variables on the 2-D coercive KYP slack, but the faces at w-1
  do better still (36%, same gain) and win on the heat benchmark: Sec. 8.
- The dual MOSEK solves has the same Gram-sized slack; a dual-only interior-point method forms the same
  Schur complement; the first-order dual route was 265x slower on this class (09/26).
- The proof's pointwise dual, as the finite hierarchy of `mincap_dual`: cap 1 for 2 - max(s,t) at D = 0;
  unbounded at D = 0 and 46830 <= 65700 at D = 2 for the rank-one target; the critical Poincare operator 191
  at input degree 8 and unbounded at 12 for D = 1 (no bounded weight), 8.4 at degree 8 for D = 2 with no
  solver status at 12. Under a second per level. It detects a missing lift degree and cannot resolve a
  borderline one; it bounds the cap of a bounded, possibly nonpolynomial, weight from below.
- The UNKNOWN exits with diverging objective and residuals at 1e-9 are the fixed-degree alternative of the
  duality note: infeasible at every finite gain without an exact Farkas witness. Raising the weight by one
  fixed every such run.

## 7. Memory files written (C:/Users/mpeet/.claude/projects/.../memory)

`separated-form-poscopvar-lift`, `poscopvar-direct`, `copquadvar-reachable-bases-padding`, `lpi-ineq-sop`,
`heatNd-degree-tailoring`, `hinf-slack-sizing`, `certificate-dual-structure`, `degree-map-proof-to-code`,
`session-handoff-2026-10-07` (status updated), `bash-heredoc-collapses-backslashes` (the sed hazard).

## 8. Running and open

- The tensor set on the 2-D heat benchmark (R (1,1), Q (2,2) uncapped, d = 0, certified bisection):
  faces at w-1 certify 12.336731 (heavy's rate) at nx 339209 in 16 solves; the tensor set at w-1
  certifies 12.326178 at nx 395785 in 15 solves and twice the time; the tensor set at w has nx 762889.
  On the KYP slack faces at w-1 give 0.171722 at nx 200425, below the tensor set's 233280. Verdict: the
  faces one degree below the plain term keep the reach at a third of the variables on both 2-D problems;
  the tensor set is neither smaller nor stronger. The N-D default of `get_lift_degs` is now faces at w-1
  whenever every weight is at least 2 (faces at w otherwise); the 1-D pair is unchanged.
- Open: the 2.5e-6 floor of the primal Q form; a weight rule that reproduces the measured 2-D choice
  without the manual raise; P at the heavy 2-D degrees and (2,3) or (3,2) with faces in 2-D (beyond the
  daytime budget); the dual coercive executive; the stock coercive 2-D program to convergence (night window).

## 9. Repository state

- Committed: bc3be32f on `ndopvar` (18 files, 4345 lines): the separated form, the direct map, the reader,
  `lpi_ineq_sop`, the tests, the harnesses, the heat tailoring, the notes and READMEs. Not pushed.
- Uncommitted on top of it: `poscopvar_direct.m` (quadratic codes), `get_lift_degs.m` (N-D default:
  faces at w-1 when every weight is at least 2), `lpi_ineq_sop.m` (doc), `heatNd_lpi.m` ('tensor'),
  `hinf_tailor_1d.m` and `hinf_tailor_2d.m` ('quad', 'tensor', 'facesm1' options),
  `sopvar_lift_notes.tex/.pdf` (Sec. 10), `lpi_programming_sopvar/README.md`,
  `PIETOOLS_demos/sopvar_demos/README.md` (Sec. 12 addition), `claude_tests/mincap_dual.m`,
  `test_mincap_dual.m`, and this file. `test_lpi_ineq_sop` passes with the new default (6 checks).
- Left alone on purpose: `pielr_solve.m`, `pielr_opcheck.m`, `jp_wave_plate_rates.m` (another session's
  gate-rule work, CC 10/07).
- Committed later on 10/08: 96ef81f6 (the certificate as an optimization problem and its dual; faces at
  w-1). Uncommitted after it: everything of Sec. 10 below (`executives_sopvar/`, the new
  `lpi_programming_sopvar` files, and the edits to `get_lpivar_degs_sop`, `lpi_eq_sop`, `lpi_ineq_sop`,
  `copvar2opvar`, `lpi_eq_cdopvar`, `collect_eq_rows`), plus this file. `eye_copvar_sop` was
  overwritten by mistake during the work and restored from git (no change).

## 10. The executives on the container path (executives_sopvar, 10/08, UNCOMMITTED)

One `PIETOOLS_<name>_sop` per stock executive, 28 files (18 1-D including the auto_execute script, 10
2-D), on four private builders (`hinf_build_sop`, `stability_build_sop`, `h2_build_sop`,
`synth_build_sop`) with `exec_tools_sop` (per-dimension storage, stock slack, eppos identity, free-variable
degrees, container-to-opvar) and `exec_info_sop`. The LPIs are the stock ones as the 10/06 transcriptions
wrote them; what is new: the slack is `lpi_ineq_sop` 'like' the storage (dw_Q 1, dw_P 2), rebuilt with
the weight raised while MOSEK returns UNKNOWN or primal infeasible (`lpi_solve_loop_sop`, at most 2
raises); the dual kernel of the negativity constraint is read back (`lpigetdual_sop`); a Legendre
Galerkin discretization (`pie_disc_sop`, quadrature on the recurrence-evaluated basis, exact to degree 40)
gives the numerical counterpart (`pie_witness_sop`: frequency response, spectrum of (A,T), Lyapunov H2,
alignment of the worst input with the dual kernel); `settings.sop.slack = 'stock'` restores the stock
slack for comparison; `settings.sop.gam_fixed` poses the feasibility test at a fixed gamma. Library
changes outside the folder: `get_lpivar_degs_sop` (N-D per-role reading), `copvar2opvar` (non-square
grids; variable-free containers), `lpi_eq_sop` / `lpi_ineq_sop` (tag, margin), `lpi_eq_cdopvar` /
`collect_eq_rows` (row bookkeeping), and the new `settings_2d_sop`,
`poslpivar_settings_2d_sop`, `exec_sop_settings`, `lpi_classify_sop`, `lpi_solve_loop_sop`,
`lpi_shape_sop`, `op2copvar_sop`, `on_registry_sop`.

Measured (light, MOSEK 11; drivers `run_exec_1d`, `run_exec_1d_b`, `run_exec_2d` in the session
scratchpad; `executives_sopvar/tests/test_executives_sop.m` asserts the 1-D part):

- **Discretization.** The first `pie_disc_sop` expanded the Legendre basis in monomials and summed
  closed-form moments: exact at degree 8, spurious eigenvalues at 16, singular T at 24 (the default), so
  every witness of the first stability run was wrong (max Re 548 for a plant at -4.93). The quadrature
  version gives eig(A,T) = -4.9348, -34.5436, -83.8916 and the static gain 0.1824791 at every degree 8 to
  40, cond(T) 1.3e5 at 40; the io1 witness is 0.1824791 at omega 0 against the certified 0.18248166 (gap
  1.4e-5), the rate -4.93480, the H2 norm 0.287684.
- **1-D, distributed disturbance (io1, syn1, rd 0.5).** Stability x4 and well-posedness certify in both
  modes and agree with the stock. H2_norm_c: stock 0.28878602, 'stock' slack 0.28879091, 'like'
  0.28775623 against the numerical 0.28768384 (gap 2.5e-4 against the stock's 3.9e-3). Every DUAL-form
  LPI (Hinf_gain_dual, _dual_coercive, Hinf_control, H2_control, H2_norm_c_coercive) is infeasible or
  UNKNOWN in both modes and in the stock (numerr 2, gamma 1e4): the Lyapunov block (AP)T' + T(AP)' is
  compact (T' is purely integral) and the Schur complement asks it to dominate B B' = I on L2; not a
  defect of either implementation. H2_norm_o and _o_coercive take the trace over the R^nw block, which
  is empty for a distributed w: both return the solver floor (1e-4, numerical H2 0.2877), as the stock.
- **1-D, finite disturbance and input** (reaction-diffusion, w through s(1-s), u constant, y = int x,
  z = [int x; u]; numerical gain 0.0333106, numerical H2 0.0523601). All twelve gain, H2 and synthesis
  executives certify with the 'like' slack at dw 0; the 'stock' slack reproduces the stock in every
  case (value to 5-7 digits where the stock solves; UNKNOWN where it does not):

  | executive | stock (light) | 'stock' slack | 'like' slack | numerical |
  |---|---|---|---|---|
  | Hinf_gain | 0.033332351 | 0.033332353 | 0.033310666 (gap 1.5e-6) | 0.0333106 |
  | Hinf_gain_coercive | UNKNOWN 0.0574 | UNKNOWN 0.0575 | 0.033311035 | 0.0333106 |
  | Hinf_gain_dual | 0.03333238 | 0.033332403 | 0.033310842 | 0.0333106 |
  | Hinf_gain_dual_coercive | UNKNOWN 0.0584 | UNKNOWN 0.0583 | 0.033313355 | 0.0333106 |
  | H2_norm_c | 0.052715641 | 0.052715636 | 0.052428236 (gap 1.3e-3) | 0.0523601 |
  | H2_norm_o | 0.052715396 | 0.052715396 | 0.052428398 | 0.0523601 |
  | H2_norm_c_coercive | 0.0825005 | 0.0824501 | 0.052396837 (gap 7e-4) | 0.0523601 |
  | H2_norm_o_coercive | 0.0827630 | UNKNOWN 0.0823 | 0.052393641 | 0.0523601 |
  | Hinf_control | UNKNOWN 0.0569 | UNKNOWN 0.0569 | 0.032865781 | open loop 0.0333 |
  | Hinf_estimator | 2.11e-4 | 3.18e-4 | 1.39e-4 | floor |
  | H2_control | UNKNOWN 0.0829 | UNKNOWN 0.0829 | 0.052222382 | open loop 0.0524 |
  | H2_estimator | 1.71e-4 | 2.40e-4 | 1.63e-4 | floor |

  The stock light slack is short in the lift for every coercive form (UNKNOWN, or 57% loose for the
  coercive H2 norms); the 'like' slack at dw 2 certifies all of them within 1e-3 of the numerical value.
  The two estimator values are at the eppos floor (the error output is nearly unobservable from y).
  The dual alignment is 1.000 on the finite-disturbance plant for both primal gain forms.
- **2-D (rd2 at 0.5 lambda*, io2).** stability_2D and stability_dual_2D: the stock returns UNKNOWN
  (numerr 2, 23 s); 'like' certifies in 8-11 s (nx 187984, m 3075, faces [2 2] at w [2 2]); 'stock'
  slack reproduces the stock's UNKNOWN (nx 179840, m 3456). On io2 with 'like', no raises: Hinf_gain_2D
  0.17172212 in 33 s (nx 200426, the 10/08 KYP value), Hinf_gain_2D_non_coercive UNKNOWN (the 2-D Q
  form), Hinf_gain_dual_2D and H2_norm_2D_c INFEASIBLE (distributed w, structural), H2_norm_2D_o at the
  floor (empty R^nw trace), H2_norm_2D_c_non_coercive INFEASIBLE at dw 0 (6 s), H2_norm_2D_o_non_coercive
  UNKNOWN at dw 0 (nx 624536, 95 s), Hinf_estimator_2D 5.1e-4 in 24 s. The stock Hinf_gain_2D on io2 returns UNKNOWN (1.2722) after 413 s; the other stock 2-D gain
  and H2 programs were not run to completion within the daytime budget.
  Two defects surfaced by the 2-D H2 Q forms and fixed: `lpi_ineq_sop` passed the full-integral index 4
  of a separable direction to the unseparated lift (now expanded to 2 and 3, lossless), and
  `h2_build_sop` declared W on a distributed input space through the 1-D `poslpivar_sop` (now the LF
  storage in 2-D).
- **Dual read-back, checked on the Posopvar session's question.** With `sos_opts.simplify` 0 and 1
  alike, b'y equals the primal objective (0.033310669 against 0.033310666; 0.033310717 against
  0.033310712): sossolve restores the multipliers of removed rows as zeros in the original row
  coordinates, and the b normalization leaves y unscaled; the two dual kernels differ only in how the
  multipliers are spread over dependent rows. The row record `prog.sopeq` is opt-in (the builders set
  `prog.sopeq = {}` after `lpiprogram`): a program built by any other caller of `lpi_eq_cdopvar` gains
  no field.
- **Independent regression check (Posopvar session, 10/08, against b3e24923, clean copy vs this working
  tree, single resolution asserted).** Programs bit-identical: 91/91 SDP leaves of the w1-w3 benchmark
  builds, 126/126 of the 1-D container programs (cx_run_1d cases x 6 presets; the same 5 'extreme'
  cases fail in both). 3-D heatNd timing in alternating fresh processes: equality stages equal within
  noise (eq15a 1.51/1.63 vs 1.74/1.74 s, eq15b 2.11/2.06 vs 2.14/2.28 s); totals 46.64/47.07 vs
  47.62/47.90 s, the 0.8-1.0 s in the Q stage, which no edit touches (path length is one candidate;
  two pairs cannot separate it from noise). Memory not measured.
- **Tests.** `executives_sopvar/tests/test_executives_sop.m`: 20 checks pass (every 1-D executive on the
  finite-disturbance plant within tolerance of the numerical value, the stock-slack equalities, the
  four stability executives and well-posedness at -pi^2/2, and 2-D stability).
- **Write-up (10/09, second version at the maintainer's request).** `sopvar/sopvar_duality_pruning_notes.tex`
  / `.pdf` (19 pages, for PhD students familiar with PIEs): the canonical lift and the certificate map
  with the kernel formula; the necessity statement (proof_roadmap (2)) with its status, the perturbation
  theorem and the Poincare degree lower bound; the minimum-cap problem and the uniform-cap criterion
  (proved); the exact dual over signed covariance fields (Theorem 1 of the min-cap document, proof
  included), the capped criterion, the preannihilator form, the fixed-degree alternative, the finite
  hierarchy with the exact level dual and the lower-bound/monotonicity lemma (proved); the SOS weight
  cone and the primal SDP of lpi_ineq_sop with the cap it certifies; the SDP dual in coefficient
  coordinates (dual kernel, pairing lemma, Farkas form with the closedness caveat); the map to the
  dual problem: coefficient pairing = kernel pairing after the monomial Gram, the covariance field of a
  dual kernel as the diagonal of its lift (proved), and the theorem that the finite-degree dual is the
  exact dual with pointwise positivity of the field relaxed to weighted moment positivity of order w,
  with the positive-part normalization relaxed to a moment form (the cap SDP dual); the H-infinity
  instance (unit input-plus-output energy from the gamma row) and the UNKNOWN exit as weak
  infeasibility; sizing and pruning with proved/measured/assumed labels; the executives, witness,
  brackets and the compactness argument for the dual forms.
- **Open.** The 2-D Q forms (gain and H2) need a weight rule beyond dw 1 (the degree loop is the
  current answer); the dual alignment is a coefficient-space heuristic (0 for a distributed w on io1,
  1.000 on the finite-disturbance plant); the stock bisection options of the 2-D executives are not
  reproduced; the 2-D `Zop_deg` cap is one number on every role.
