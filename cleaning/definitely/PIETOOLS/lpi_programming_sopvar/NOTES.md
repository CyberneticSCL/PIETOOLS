# lpi_programming_sopvar — measurements and history

Moved verbatim from `README.md` on 10/16/2026 (MMP). File locations are as of the dates given.
On 10/16/2026 the executive helpers moved to `executives_sopvar/utility_functions/`
(`poslpivar_sop`, `poslpivar_settings_sop`, `get_lpivar_degs_sop`, `trace_rn_sop` among them), the
constructors to `sopvar/constructors/`, `op2copvar_sop` and the operator converters of
`sopvar/Testfolder/converters/` (`opvar2copvar`, `opvar2d2sopvar`, `dpvar2sdvar`, ...) to
`converters/operator_converters/`, `on_registry_sop` to `sopvar/misc/`, and `examples/hinf_tailor_1d`,
`examples/hinf_tailor_2d` and `examples/volterra_tailor` to
`sopvar/Testfolder/sdopvar/claude_tests/slack_sizing/`.

## The slack sizing of `lpi_ineq_sop`, measured (10/08/2026)

**Measured on the H-infinity gain LPI, 1-D (`examples/hinf_tailor_1d`, 10/08/2026).** The primal
KYP LPI of the plant `io1` of `cx_plant` (reaction-diffusion, Dirichlet, `z = int(x)`, closed-form
gain `0.182479080`), `R` at the stock `light` or `heavy` degrees, `Q` from `get_lpivar_degs_sop`,
the slack `N = -K` sized three ways. Every slack with `D >= 2` and `w >= 2`, and `(1,3)`, gives
`gamma = 0.1824816 +- 1e-7` at both levels of `R` (2.5e-6 above the closed form, a residual neither
`R` nor the slack sets); the stock `light` slack `{2,[1 2 3],[1 2 3]}` is `(D,w) = (1,2)` with joint
cap 3 and gives `0.182626`; `(2,1)` gives `0.184185`. The smallest sufficient slack, `(2,2)` with
joint cap 3 and the product term at `w - 1`, has 396 decision variables and `nx = 898` against 587
and 1274 for the stock `light` slack and 1137 and 2377 for the stock `heavy` one
(`{3,[2 2 3],[2 2 3]}`: the same `(2,2)` and cap beside a separate degree-3 multiplier block). The
reader sizes the same slack at `(5,3)`, `nx = 5242` (`dw = 0`: `(5,2)`, 2834): `A'*Q + Q'*A` is
composed before `T'*Q = R` is imposed, so its support (lower kernel to degree 9, `Dmin = 4`) is
that of the whole `lpivar` family of `Q`, not of the optimal `Q`, whose unreachable part the
solver sets to zero at no cost in `gamma`. For a target linear in a free operator the support rule
measures the family; applied to `R` instead (`Dmin 1`, `Dmax 4`) the reader gives `(2,2)` at
light with its defaults: `opts.like` of `get_lift_degs` / `lpi_ineq_sop` (`nx` 1058, same `gamma`;
at heavy, where `R` has a degree-4 multiplier, `(2,3)` at `nx` 1881 against the stock 2377). At
`(1,1)` the gain stays finite, `0.193093` with the product pair and `0.186930` with two faces. On a fixed target the rule holds (`examples/volterra_tailor`, Volterra
norm `2/pi`): `D = 1`, excess `1.9e-4` at `w = 1` with the product pair (`nx` 126), `1.1e-6` with
two faces (301), `3.2e-8` at `w = 2` with the product pair (326), `9.9e-9` with two faces (676);
the stock-style pair at degree 2 gives `8.0e-9` at 883.


**The coercive form, 1-D (`hinf_tailor_1d(...,'P')`, `PIETOOLS_Hinf_gain_coercive`: `P >= 0`
at the storage degrees, `K33 = A'*P*T + T'*P*A`, no free operator).** A different picture. The
stock `light` slack is short: `gamma = 0.1928`, 5.7% above the closed form, MOSEK status UNKNOWN
(as the baseline suite recorded for this executive); the stock `heavy` slack is exact to 6.6e-7.
Slacks with lift 1 diverge (objective to 1e6: under the margin `P >= eps I` the restricted program
is infeasible); `(2,1)` doubles the gain. The reader on `K` is the right size here, because the
support of `A'*P*T + T'*P*A` is that of the positive `P` at its declared degrees under
composition: `(4,3)`, or `(4,2)` with `dw = 0`, reaches the floor of the light storage (1.1e-6 to
1.4e-6 above the closed form) at `nx` 1786 against the stock `heavy` 2332; `like` (the reader on
`P`, `(2,2)`) is 1.1e-5 above at `nx` 810; `(3,2)` and `(2,3)` are at the floor for 1246 and
1426. The joint cap 3 costs 4e-3 here. With `P` at `heavy` the floor moves to 1e-7, so in the
coercive form the storage binds at 1e-6, while the `Q` form stalls at 2.5e-6 with either `R`: that
residual belongs to the `Q` form (free `Q` under the `get_lpivar_degs` caps, `T'*Q = R`), not to
the slack. The two-component plant repeats every row.


**The dual forms, 1-D.** `'Qd'` (`PIETOOLS_Hinf_gain_dual`: `T*Q = R`, `K33 = Q'*A' + A*Q`)
behaves as the primal `Q` form but reaches the closed form: at light `R` the stock slack is 7.2e-5
above, `like` (`(2,2)`) 2.7e-6 at `nx` 1058, `(2,3)` 1.5e-7 at 1858, `(2,2)` faces 2.8e-7 at
2019, the reader on `K` (`(5,5)`) 1.1e-7 at 12410; at heavy every slack from `(2,2)` up is within
2e-7. The 2.5e-6 floor is specific to the primal `Q` form. `'Pd'`
(`PIETOOLS_Hinf_gain_dual_coercive`, `K33 = T*P*A' + A*P*T'`) is unusable on this plant as
shipped: the stock executive returns `gamma = 8252` at light (as the baseline suite recorded), the
container stock slack 10785 with status UNKNOWN, and every tailored slack from `(2,2)` to `(3,3)`
(pair or faces) is reported primal infeasible at both storage levels.


**2-D (`examples/hinf_tailor_2d`, io2 at light `R`, 10/08/2026).** No tailored slack admits a
finite gain: `(1,1)`, `(2,1)`, `(1,2)` with four faces, `(2,2)` with the product pair (with and
without a total-degree cap of 6) and `(2,2)` with four faces (`nx` 5.1e4 to 9.8e5) all end with
MOSEK status UNKNOWN, residuals at 1e-9 and the objective running to 1e1-3e3 (the uncapped
`(2,2)` faces run stalls at `gamma` 12.4). The stock light 2-D slack is `R`'s degrees + 3 in every
slot (`Dup = 3`), a 1.5M-variable program stopped after 9 minutes, and the baseline suite (09/24)
had found the stock 2-D executive on this plant 7.4x conservative with the same status. In the `Q` form the LPI gives no clean 2-D reference at any sizing tried; the heat benchmark
(`PIETOOLS_demos/sopvar_demos/README.md` Sec. 12) and the coercive form below carry the 2-D
evidence.


**Variables and the dual (notes Sec. 9, 10/08/2026).** Against the proof's minimum-cap problem and
its dual: the rows are independent in every solve (MOSEK presolve), so the storage parametrisation
of the null weights removes nothing; the Gram overhead `(w+1)^2/(2w+1)` is the price of the
certificate; the one lever is the term set. New `poscopvar_direct` codes `2nv+2+d` (the box
quadratic in direction `d` alone, offset lowering `int` there only) give the Markov-Lukacs tensor
set `[0, 2nv+2+(1:nv), 1]` at offsets `[0, 1..1, 1]`: on the 2-D coercive KYP slack at `(2,2)` it
reaches `gamma = 0.171720` at `nx` 233280, `m` 4116, 17.6 s against the four faces' 0.171718 at
557845, 5420, 57 s (`'deg 2 2 tensor'`); in 1-D the new code reproduces the product code exactly.
The 2-D heat benchmark decided against it (R (1,1), Q (2,2) uncapped, d = 0, certified bisection):
faces at `w - 1` certify 12.336731, the heavy rate, at `nx` 339209 in 16 solves; the tensor set at
`w - 1` certifies 12.326178 at `nx` 395785 in 15 solves and twice the time; the tensor set at `w` has
`nx` 762889. On the KYP slack faces at `w - 1` give 0.171722 at `nx` 200425, below the tensor set.
The N-D default of `get_lift_degs` is now the plain term at `w` with the `2N` faces at `w - 1`
whenever every weight is at least 2, and at `w` otherwise. The dual MOSEK solves
has the same Gram-sized slack; the proof's pointwise dual is a degree diagnostic
(`claude_tests/mincap_dual`, `test_mincap_dual`: cap 1 for `2 - max(s,t)` at `D = 0`, unbounded at
`D = 0` and below 65700 at `D = 2` for the rank-one target, the Poincare operator unbounded at
`D = 1`), under a second per level, and it explains the UNKNOWN exits as the fixed-degree
alternative of the duality note.

**The coercive form in 2-D (`hinf_tailor_2d(...,'light','P')`, `P` at the stock light degrees
plus the margin).** The reader on `P` (`like`) gives `(2,1)` per direction with four faces: `nx`
1.1e5, `m` 2788, solved to optimality in 5 s, `gamma = 0.2665` (56% above the closed form). The
product pair fails at `(2,2)`, `(3,2)`, `(2,3)` (UNKNOWN) and succeeds at `(3,3)`: `gamma =
0.171777`, 0.40% above the closed form, `nx` 1.43e6, `m` 13659, 491 s. Four faces do what the product pair cannot: `(3,1)` 0.1974 (+15%) at `nx` 3.4e5 in 25 s;
`(2,2)` `gamma = 0.171718`, 0.36% above the closed form, `nx` 5.6e5, `m` 5420, 57 s; with cap 6,
0.171719 at `nx` 4.9e5 in 36 s. The stock light
slack (`P`'s degrees + 3 in every slot, `nx` 7.6e5, `m` 13978) was stopped by the budget after
459 s at iteration 14 (objective 0.34 and falling); the baseline suite (09/24) ran the stock
executive on this plant to `gamma = 1.272` with status UNKNOWN in 283 s. The 2-D extrapolation
rests on the coercive form: `like` gives a clean bound at a seventh of the stock size, and `like` with `dw = 2` (`(2,2)` with
four faces) is within 0.36% at three quarters of it. In both dimensions the coercive form wants
`like` with `dw = 2`: `(2,3)` in 1-D at the floor, `(2,2)` faces in 2-D.

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
