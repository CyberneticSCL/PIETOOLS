# Example-library baseline, 10/28/2026

Every runnable example of the three example libraries
(`examples_PDE_library_PIETOOLS`, `examples_DDE_library_PIETOOLS`,
`examples_NDSDDF_library_PIETOOLS`), run on each settings file through the
executive the library calls. The stock results were recorded before the opvar
path is removed, so that the sopvar path can be compared row by row.

## Files

| File | Content |
|---|---|
| `stock_results.csv` | Stock path: stock executives, `lpisettings(type,0,'','mosek')`. One line per (example, settings file). |
| `sop_results.csv` | sopvar path: the `_sop` executive and form `PIETOOLS_auto_execute_sop` uses for each stock executive, `lpisettings_sop(type,0,'','mosek')`. Same columns. |
| `stock_split.csv`, `sop_split.csv` | Every row of the two results files with its time split into construction and solve, and its equality count (columns below). |
| `comparison.csv` | One line per (example, settings file) run on both paths: both results side by side, the split of each, and two verdicts (Criterion). |
| `catalog.csv` | The 142 catalog rows: executive, parameter, documented value, bisection bracket, spatial dimension, PIE size, run flag. |
| `runs.log` | One line per harness call: time, path, settings files, row count, git state, MATLAB and MOSEK versions, threads. |
| `hinf2d_results.txt` | The 2-D diagnoses (failures 2, 4, 5): the parameter cells of T and A of PDE 45, one `##FP` line per build-only SDP (PDE 45, then PDE 42 and 43 under the `hd_l1v` headers), one `##SOL` line per solve (PDE 45–47). |
| `private/` | The harness (`lb_run.m` and helpers, `catalog.mat`, and `blocks/`, the DDE and NDS/DDF library blocks, one per file and uncommented), the log parsers (`lb_mosekrows.awk`, `lb_split.awk`), the comparison (`lb_compare.awk`, `lb_augment.awk`), the diagnosis reruns (`lb_diag.m`), and the 2-D H∞ diagnosis (`hd_build.m` with `shadow/Sedumi2Mosek.m`, which records the SDP and stops before MOSEK; `hd_solve.m`; `hd_l1v.m`, the one-variable-space check). The `hd_*` scripts ran next to the harness folder, which they expect as `../library_baseline`. `private/` keeps it off the MATLAB path. |

## Columns

`library, id, settings, status, result, documented, n_solves, seconds, decision_variables, note`

- `status`. Single solve: `certified`, `inaccurate`, `infeasible`, `unknown`.
  Bisection: `converged`, `max_solves`, `lo_fail` (lower bracket not
  certified), `hi_pass` (upper bracket certified). Also `error` (message and
  file:line in `note`), `skipped` (2-D, projected over 600 s per solve) and
  `crash` (MATLAB died twice in the row).
- `result`. Stability: 1 certified, 0.5 inaccurate, 0 otherwise. Gain or norm:
  the bound the executive returns. Bisection: the largest certified parameter
  value.
- `seconds`: wall time of the row, setup and every solve.
- `decision_variables`: `numel(prog.decvartable)` of the last solve.
- `note`: the executive called; the solver exit (`fr` feasratio, `pinf`,
  `dinf`, `numerr`, `iter`, `cpu`); the SDP size (`m` rows, `nx` cone
  entries); and, for a bisection, the bracket and every step as
  `value:C/I/F/U/E(seconds)`.

Split files: `library, id, settings, path, status, n_solves, n_sdp, seconds,
solve_seconds, construct_seconds, equalities, decision_variables`.

- `n_sdp`: SDPs solved in the row, every bisection step and every sopvar
  slack-raise pass counted (`n_solves` counts executive calls).
- `solve_seconds`: MOSEK optimizer time ("Optimizer - time" in its log),
  summed over the row's SDPs. It excludes SOSTOOLS's conversion of the
  program to solver input.
- `construct_seconds`: `seconds − solve_seconds`: building the LPI,
  converting it to solver input, and reading the solution back.
- `equalities`: MOSEK's constraint count for the row's last SDP, i.e. the rows
  of the equality system. `NaN` when no SDP was solved (error, skipped).

These were parsed from the MATLAB logs of the runs (row markers `#####` and
MOSEK's solution summaries), which stay in the scratch folder.

`comparison.csv` adds to the columns of `lb_compare.awk`, per path (`stock_`,
`sop_`): `sdps, solve_seconds, construct_seconds, equalities` and `mosek`
(solution status of the row's last SDP), then `verdict_mosek`.

## Criterion

The one in the Conclusion block of `executives/2D/*.m`, applied to every
executive:
- **certified:** |feasratio − 1| ≤ 0.3, with no numerr, pinf or dinf;
- **inaccurate:** the same, with numerr;
- **infeasible:** pinf or dinf, or |feasratio + 1| ≤ 0.1;
- **unknown:** anything else.

A bisection starts from the documented bracket. It halves until
hi − lo ≤ 1e-3·|documented|, or for at most 14 solves. A step passes only if
it is certified.

For MOSEK, `sossolve.m` fills feasratio from `MSK_DINF_INTPNT_OPT_STATUS`.
MOSEK documents that measure as converging to 1 only when the problem is
well posed. It rejects 5–10 MOSEK-optimal rows per settings file on each
path. The sopvar executives classify by MOSEK's status instead
(`lpi_classify_sop`). `verdict_mosek` therefore certifies a single-call row
when its last SDP's solution status is OPTIMAL. Bisection rows keep the
feasratio verdict, since every step was classified by it. Neither rule is an
operator-level certificate: the stability LPIs have a small b, and passing
needs the residual-and-PSD gate of the acceptance tests.

## Stock coverage

138 runnable rows on 6 settings files.

- **extreme, stripped, light, heavy:** 138 rows each.
- **veryheavy:** 137 rows. DDE16 (H∞ control) was not run; it was projected at
  2.0 h per solve.
- **custom:** 107 of 137 rows. DDE16 was not run, and DDE22 through the end of
  the NDS/DDF library is paused (10/28 16:17). To resume, delete `PAUSE` in
  the scratch folder and rerun `lb_full('reduced')`.
- **2-D rows:** skipped where the projected time per solve exceeds 600 s. The
  growth floor is 1 per step into stripped or light, and 10 per step into
  custom, heavy or veryheavy.

## sopvar coverage

138 runnable rows on 4 settings files: stripped, light, heavy and veryheavy.
`custom` is `light` in `settings_sopvar` and was not run.

- **Order:** light and heavy first; then veryheavy and stripped. The eight
  2-D rows that took over 60 s per solve on stock extreme (PDE 36, 38, 39, 45,
  46, 47, 50 and the PDE 50 dual) ran last, on every settings file.
- **2-D projection:** the same rule as stock, with max 600 s per solve. Of
  those eight rows, heavy solved PDE 45 and 46 and veryheavy solved none.
  PDE 50 and its dual errored on all four settings files before any solve.

## Results

Verdict counts per settings file, feasratio rule / MOSEK rule:

| Verdict | stripped | light | heavy | veryheavy |
|---|---|---|---|---|
| same | 19 / 23 | 58 / 63 | 73 / 79 | 91 / 95 |
| sop tighter | 2 / 4 | 34 / 39 | 23 / 29 | 9 / 11 |
| sop only | 86 / 85 | 16 / 13 | 8 / 5 | 2 / 0 |
| neither | 22 / 17 | 13 / 9 | 12 / 8 | 11 / 8 |
| stock only | 1 / 1 | 6 / 3 | 5 / 0 | 1 / 0 |
| sop looser | 0 / 0 | 1 / 1 | 3 / 3 | 5 / 5 |
| sop error | 6 / 6 | 6 / 6 | 6 / 6 | 6 / 6 |
| stock error | 1 / 1 | 1 / 1 | 1 / 1 | 1 / 1 |
| skipped | 1 / 1 | 3 / 3 | 7 / 7 | 11 / 11 |

Time split over the rows both paths ran:

| Settings, dimension | Rows | Stock solve / construct (s) | sopvar solve / construct (s) | SDPs, stock → sopvar | Equalities, sop/stock | Decision variables, sop/stock |
|---|---|---|---|---|---|---|
| stripped 1-D | 120 | 5 / 118 | 107 / 153 | 121 → 473 | 1.54 | 4.17 |
| stripped 2-D | 10 | 2402 / 191 | 453 / 52 | 10 → 13 | 0.47 | 0.28 |
| light 1-D | 120 | 43 / 489 | 129 / 129 | 307 → 589 | 1.10 | 1.65 |
| light 2-D | 8 | 846 / 52 | 173 / 36 | 8 → 12 | 0.87 | 1.96 |
| heavy 1-D | 120 | 58 / 532 | 304 / 142 | 328 → 572 | 1.23 | 2.36 |
| heavy 2-D | 4 | 1975 / 419 | 149 / 16 | 4 → 4 | 0.44 | 0.38 |
| veryheavy 1-D | 119 | 947 / 8569 | 2965 / 243 | 327 → 693 | 0.98 | 0.82 |

Ratios are geometric means over rows.

- **1-D:** MOSEK takes longer on sopvar. The SDPs are larger, and there are more
  of them, because `lpi_solve_loop_sop` re-solves at a raised slack degree
  after every 'unknown' or 'infeasible' pass. An infeasible bisection step
  therefore costs `max_raise_feas + 1` SDPs.
- **Stock 1-D construction:** stock spends 89–92% of its 1-D time outside
  MOSEK, concentrated in a few delay rows. DDE16 heavy spends 236 s outside
  MOSEK on stock and 1.3 s on sopvar. PDE 29 veryheavy spends 1530 s, for an
  SDP MOSEK solves in 0.1 s.
- **2-D:** MOSEK time itself falls, because the sopvar SDPs have fewer
  equalities. PDE 35 heavy: 12,468 → 5,297 equalities, MOSEK 623 → 30 s.

## sopvar failures (diagnosed 10/28; reruns in `private/lb_diag.m`)

1. **DDE19, DDF03: `opvar2copvar:empty` in `PIETOOLS_Hinf_control_sop`.**
   - The PIEs have no regulated output: `C1.dim` is [0 3; 0 1] and
     [0 4; 0 2].
   - The DDE19 block declares z with the legacy field names `C0`, `C{1}`,
     `D1`, `D2`. `initialize_PIETOOLS_DDE` reads nz only from `C1`, `D11`,
     `D12` and their delay variants. The DDF03 block declares no output.
   - Stock returns the gain into an empty output (0.0011).
   - Renamed to `C1`, `C1i{1}`, `D11`, `D12`, DDE19 certifies 1.9639855 on
     both paths at light (documented 1.964).
2. **PDE 42, 43, 50 (and the PDE 50 dual): "Space {s1}: only R^n and L2[x;y]".**
   - `pos2d` and `eqopts2d` in `poslpivar_settings_2d_sop.m` build bases only
     for R^n and L2[s1,s2].
   - These PIEs have one-variable states: `T.dim` is [1 1; 1 1; 1 1; 1 1],
     [0 0; 1 1; 0 0; 1 1] and [0 0; 1 1; 0 0; 2 2].
   - The stock 2-D settings give such states degrees in `LF_deg.dx` and
     `LF_deg.dy`, which nothing on the sopvar path translates.
   - The gap was known since 10/08.
   - Fixed in the working tree on the evening of 10/28, uncommitted, by the
     "Posopvar structure" session. `pos2d` declares L2[x] and L2[y] bases
     from `dx`, `dy` through the 1-D translation of `settings2possopvar`;
     `eqopts2d` adds the exclude(2:7), sep(1:2) rules of `get_eq_opts_2D`.
   - Verified here (`hd_l1v.m`). The reference is stock `PIETOOLS_PDEstability`,
     whose `_2D` LPI is the P form that `PIETOOLS_stability_sop` poses in 2-D
     (see item 5). Builds with the stock slack match stock in m, cone, PSD
     blocks and coefficient sums (10 digits), with the stripped settings and
     with unequal caps in the variable and in θ. The stock settings files
     all have equal caps, so the second case is what tests the order:

     | PDE, settings | m | Cone | PSD blocks |
     |---|---|---|---|
     | 42, stripped | 4,198 | 124,053 | 24, 497 |
     | 42, unequal `dx`, `dy` caps | 4,437 | 132,247 | 28, 513 |
     | 43, stripped | 3,656 | 113,697 | 18, 476 |
     | 43, unequal `dx`, `dy` caps | 3,769 | 117,580 | 20, 484 |

   - Five new `error`/`warning` lines in the change lack the column-77
     marker.
3. **Light/heavy stock-only rows of PDE 3, PDE 40, PDE 27, DDE26, DDE11 (and
   dual), DDE25.**
   - Each is one SDP that MOSEK solved to OPTIMAL (PRIMAL_AND_DUAL_FEASIBLE).
     Equality violations are 8e-9 to 2e-6 and the Gram blocks are exactly PSD:
     the same order as the stock rows accepted.
   - Only the feasratio window rejects them (Criterion). None remains stock
     only under `verdict_mosek`.
4. **PDE 45, 46, 47 (2-D reaction-diffusion, H∞ gain, exact 0.111, 0.171,
   0.207): the LPI form, not the construction.** Diagnosed 10/28 with
   `private/hd_build.m` and `private/hd_solve.m`; output in
   `hinf2d_results.txt`.
   - Baseline: sopvar returns 2,594, 963 and 22,602 at light, and 54.8 and
     8.59 at heavy (PDE 47 heavy skipped). MOSEK ends UNKNOWN on every pass,
     with dual norm ~3e10.
   - Cause: `PIETOOLS_Hinf_gain_sop` poses a 2-D PIE in the Q form (R = T'Q
     >= 0, Q free). Stock `PIETOOLS_Hinf_gain` sends a 2-D PIE to
     `PIETOOLS_Hinf_gain_2D`, the P form (P >= eppos I, K33 = (PA)'T +
     T'(PA)). The routing dates from 25d6d82d (10/09, then to
     `PIETOOLS_Hinf_gain_2D_non_coercive_sop`); the 10/28 consolidation kept
     it. `PIETOOLS_stability_sop` poses 'Q' as 'P' in 2-D;
     `PIETOOLS_Hinf_gain_sop` does so only for 'Qd'.
   - Build-only programs of PDE 45 (the shadow `Sedumi2Mosek` stops before
     MOSEK):

     | Program | Equalities | Scalarized cone | PSD blocks | Free |
     |---|---|---|---|---|
     | stock P (`_2D`), stripped | 13,021 | 452,768 | 13, 951 | 1 |
     | sopvar P, stock slack, stock stripped degrees | 12,962 | 452,769 | 13, 951, 1 | 1 |
     | stock Q (`_2D_non_coercive`), stripped | 16,228 | 902,737 | 13, 1343 | 150 |
     | sopvar Q, stock slack, stock stripped degrees | 16,869 | 903,030 | 13, 1343, 1 | 442 |
     | sopvar P, own slack, light | 4,094 | 100,683 | 8, 334, 4×149, 1 | 1 |

   - The earlier "1,804,260 against 904,571 variables" compared the Q form
     with the P form. T of these PIEs is nonzero only in the four
     integral-integral cells, so (PA)'T + T'(PA) has no multiplier cell and
     the stock slack drops five of the nine cells. In the Q form, K33 = A'Q +
     Q'A has all nine, and the slack Gram is 1343 against 951.
   - At matched form and degrees the cones agree. The extra 1×1 block is gam
     >= 0, which sopvar declares and stock does not. The 59-row difference in
     equalities is not diagnosed.
   - Solves of the sopvar P form, every one OPTIMAL in one pass (gam, MOSEK
     seconds):

     | Program | PDE 45 | PDE 46 | PDE 47 |
     |---|---|---|---|
     | sopvar P, own slack, stripped | 0.11690, 8.1 | 0.26671, 5.7 | 0.20855, 9.0 |
     | sopvar P, own slack, light | 0.11198, 21.4 | 0.17172, 17.4 | 0.20853, 11.5 |
     | sopvar P, stock slack, stock stripped degrees | – | – | 0.20713, 392 |
     | stock P, light | 0.11197, 301 | 0.17169, 299 | UNKNOWN 0.20851, 242 |
     | stock P, stripped | UNKNOWN 0.252, 265 | UNKNOWN 1.57, 320 | 0.20707, 277 |

   - Not measured: a solve of the stock Q form. Its program has the cone of
     the sopvar Q program, whose MOSEK solve took 1354 s. That stock's Q form
     fails the same way is inferred from the cone match only.
   - Fix (not applied): pose 'Q' as 'P' in 2-D in `PIETOOLS_Hinf_gain_sop`,
     as `PIETOOLS_stability_sop` does.
   - Correction (maintainer): the P (PDE2PDE) form is the Q (PIE2PDE) form
     with Q = PT, R = T'PT. The failure is the sizing of the sopvar Q form:
     R has the settings degrees (Gram 13 at stripped, that of P in the P
     form) and Q the caps of R, so the Q program does not contain the P
     program.
   - Build-only probe of a Q form sized to contain the P form
     (`private/hd_qform.m`):
     - Q gets the caps of P*T: int [1 3] against P's [1 1], 289 free
       coefficients.
     - R is sized from the support of T'Q and `T'*Q == R` is imposed in
       every coefficient.
     - Sized 'like' R, the KYP slack has Gram 2705 at stripped and 4226 at
       light (Σ N² 16.6 M, 47.1 M), so it was not built.
     - Sized 'like' P (dw_P), as the P form sizes it:

     | PDE 45 | m | Cone | PSD blocks | nnz(A) |
     |---|---|---|---|---|
     | P form, stripped | 2,786 | 55,913 | 8, 149×5, 1 | 0.86 M |
     | Q form, stripped | 25,854 | 824,056 | 900, 400×4, 197×5, 1 | 22.3 M |
     | P form, light | 4,094 | 100,683 | 8, 334, 149×4, 1 | 1.38 M |
     | Q form, light | 27,417 | 902,456 | 900, 400×4, 442, 197×4, 1 | 23.0 M |

   - Containment holds in degree by construction. Whether R's Gram terms
     represent T'PT exactly is not verified: a solve must not return a bound
     above the P form's (0.11690 stripped, 0.11198 light).
   - Priced solve (not run), scaled in m with the MOSEK exponent 2.35–3.3
     from the sopvar Q stock-slack program (m 16,869, 1354 s) and the P
     stock-slack program (m 12,962, 392 s): 35–90 min at stripped and 40–110
     min at light, per SDP. The possible gain is bounded by the P form's
     distance to the exact value: 5.3 % at stripped, 0.87 % at light.
5. **The 2-D comparisons of all three LPI families compare different LPIs.**
   Stock dispatches a 2-D PIE to a `_2D` executive; the forms differ as follows.

   | Harness row | Stock 2-D LPI | sopvar 2-D LPI |
   |---|---|---|
   | stability (`PIETOOLS_PIE2PDEstability`) | `_2D`: Q form (free Q, T'Q = R) | 'P' (P form; no 2-D Q form) |
   | H∞ gain | `_2D`: P form | 'Q' (Q form) |
   | H2 'o', 'c' | `_2D_o`, `_2D_c`: coercive (W only) | 'o', 'c' (Q form; 'oco', 'cco' unused) |

   - The stability P form matches stock `PIETOOLS_PDEstability_2D` (item 2).
     The baseline's 2-D stability verdicts (PDE 39 among them) compare it with
     the stock Q form.
   - PDE 50 on the default sopvar path (built by the "Posopvar structure"
     session, m about 30,640) is the Q form. Stock's 2-D H2 is coercive, so
     a like-for-like run is 'oco'/'cco'.
   - The coercive forms, built by that session at light (stripped builds the
     same program; not reproduced here):
     - 'oco': m 18,700, PSD [1 1 22 540 540 540 540 1080], about 2.33 M
       variables;
     - 'cco': m 18,708, the same cone.
     Stock stripped `_2D_o`: m 16,528, PSD [1 31 1366], 1.87 M variables,
     982 s of MOSEK. Estimated 22–25 min per form (MOSEK exponent 2.35–3.3
     in m), over the 20-min cap, so not run.
   - The four stock 2-D H2 executives do not call `get_eq_opts_2D`; every
     other stock 2-D executive does. The sopvar 'stock' slack always applies
     the exclusions, so it reproduces stock's H2 programs only up to the
     slack (PDE 50: Gram 1128 against 1366).

## Environment

- MATLAB R2025b, MOSEK 11.0, 8 threads.
- git 1cf111e5, plus the uncommitted PBC fix in `convert_PIETOOLS_PDE.m`.
- The PIEs of all 99 distinct systems are bit-identical at b515277f (checked
  10/28 16:14), so the stock results hold for b515277f.
- One MATLAB crash, an R2025b access violation (PDE 50, extreme, 06:44). The
  retry ran.
- **sopvar runs, 10/28 16:18–19:34:**
  - git b515277f with the same uncommitted PBC fix, then 8354ecb0, the commit
    holding that fix.
  - From 16:57, `converters/@pie_struct/pie_struct.m` carried an uncommitted
    edit by someone else, to header comment text only.
  - One MATLAB crash, "Illegal instruction" (PDE 51, light, 16:22). The retry
    certified.
- No other MATLAB run overlapped a recorded row.
