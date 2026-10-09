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
