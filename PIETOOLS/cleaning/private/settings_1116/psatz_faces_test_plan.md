# Plan: are the single faces needed, or useful, once the degrees are accounted per direction

Written 2026-11-15 on the maintainer's instruction: "It is possible the negative results we have are due to a poorly structured approach to settings which might be corrected by properly accounting for direction specific degree caps. Based on your list of cases where single faces seem to be required, create a plan to test your de novo settings architecture and see whether single faces are truly needed or useful, given the complexity they introduce in 2D and 3D." Companion to `psatz_terms_evidence.md` (the cases, labelled E1-E11 there) and `sopvar_settings_proposal.tex` (the row (d_1, d_2, d_3, patterns, F, o), Secs. 4-5). Nothing in this plan has been run.

## 1. Hypothesis and what a test can show

H0 (the maintainer's hypothesis). The negative results for the box-factor sets come from the degree structure of those runs, not from the term set:
(a) in E1-E3 every term shared one Gram basis (block 264), so a box factor g_k overshot the plain term by two degrees in its direction and the product g_1 g_2 by two in each; the failures recorded there are MOSEK statuses (numerr 2) or residual thresholds, never infeasibility certificates;
(b) in E6 the face offset lowered the Gram in every direction, not only the acting one, and the "tensor set" lowered only the acting direction, so the two were not at the same accounting;
(c) in E8 the subset caps of the stock 2-D storage bound, and the faces-versus-corners verdict moved when the storage changed.
Under per-direction lengths, the Putinar set {1, g_k} (N + 1 terms) or the product set (2^N terms) may certify what the faces (2N + 1 terms) certify, at comparable size.

What the theory fixes before any run (evidence note, T6): the faces at Gram length e lie in the box set at lengths (e + 1 plain, e for g_k), and the box set at e lies in the faces at e + 1. So no single-length comparison can show more than a one-length saving either way, and "needed" can only mean needed within a size or time budget. The test is therefore a ladder in the length with the certificate size recorded at every rung, and the question is answered by cost at equal reach, not by reach at one length.

Three sub-hypotheses, each with a measurement that decides it:
- H1: with g_k at length e - 1 in direction k (and e elsewhere), the box set certifies at the same length as the faces at e - 1 in the acting direction. Decided by the ladder of Sec. 4 on C1-C3.
- H2: the per-direction offset for the faces (acting direction only) beats the scalar offset (every direction) at equal reach. Decided by running both face variants on every case.
- H3: on the searched-P executive problems, the term set does not matter once the storage is the N-D default (E8). Decided by C4 with the storage held fixed.

## 2. Instrument (Phase 0; no solver time)

The architecture under test is the row of the proposal, read into `poscopvar_direct` arguments. The minimal code, in the scratch shadow first (`scratchpad/psatz_mix/shadow/`), then in the repo on the maintainer's word:
1. `poscopvar_direct`: accept `psatz_rows_offset` as an R x nv matrix (per row, per direction). The reduction routine already takes a per-direction row; the change is option parsing only (evidence note, Sec. 4 of the proposal).
2. A reader `row2posopts(row, spaces)` from (d_1, d_2, d_3, patterns, F, o) to `include`, the per-basis-operator `deg` (lengths to degrees, `mult` = d - 1 for the lift, `int` = d_3 - 1 for the weight), `psatz_rows` (F) and `psatz_rows_offset` (o). Named sets: `box` (N rows with one entry 3), `product` (all nonzero rows over {0, 3}^N), `faces` (2N rows with one entry in {1, 2}), `corners` (faces plus the 4 C(N, 2) two-face rows), `none`. Default o_{r,k} = 1 where r_k is nonzero; the scalar variant (o_{r,k} = 1 in every k) reproduces today's 'auto'.
3. `heatNd_lpi` (scratch copy under `psatz_mix/heat/`): a `rows` option that passes F and o through, next to its 'linear', 'product', 'tensor' and the scratch 'corner'.
4. A pinned-P harness on the container path (new, about 100 lines, from `pm_run.m` and `exec/mx_sys.m`): build T and A as copvars, Q = A'T + T'A with P = I, impose Q + Pos = 0 by `lpi_ineq_sop` on -Q with the row, bisect lam.
5. `hinf_tailor_2d` restored from `cleaning/private/definitely/.../claude_tests/slack_sizing/` into scratch for the KYP io2 case.

Verification before any experiment (all by comparison of built programs, no solves beyond seconds):
- 1-D: the row `faces` at o = 0 and the row `box` at o = 1 give the programs that `lpi_ineq_sop` gives today for 'faces' and the Markov-Lukacs pair, compared on `prog.expr.At` after sorting the decision names (bit-identical expected).
- Every term's monomial support, read from `info.terms(k).E`, equals the intended per-direction degree box (plain: 2(e - 1) per direction; g_k at o = 1: 2(e - 1) in k; face at o = 1: 2(e - 1) - 1 in k; face at o = 0: 2(e - 1) + 1 in k).
- The estimator (`options.estimate`) gives the Gram entry count of every configuration of Secs. 3-4 before a solve, so the budget of Sec. 6 is checked, not assumed.

## 3. Cases, sets and lengths

Cases, taken from the list where the faces seemed required (evidence note E1-E6) plus the controls that decide whether the instrument is right:
- C0, one variable, control. DD heat (`heatNd_lpi`, N = 1) and io1 KYP (`hinf_tailor_1d`). By T1 the box set at (e, e - 1) is complete at degree 2(e - 1); the faces must not beat it at matched degree. If they do, the instrument is wrong (memory `test-1d-before-2d`).
- C1, P = I pinned, 2-D reaction-diffusion, Dirichlet, lam* = 2 pi^2, on the container path (replaces the retired opvar2d harness of E1). Add the three plants where the product certified nothing and the faces did: L1 (reaction lam (1 - 2 s_1) x), L2 (diffusion 0.05 + s_1 (1 - s_1)) and the ODE-coupled plant of E2, each by its own margin grid.
- C2, 2-D heat DD x DN, kappa* = 12.337006, d = 0 and d = 1 (E5, E6), with the drop-one-face ablation of E5 repeated under the new accounting (the box set has no single face to drop; the ablation becomes: box set with g_2 only, and the faces with the Neumann face removed).
- C3, 2-D coercive KYP io2 at (2,2) and (3,3) (E6): gamma against 0.1711 closed form.
- C4, searched P on the executive path: S3 and S8 at 0.9 r* with the N-D default storage (route patch `psatz_mix/edits/route_edit.patch`, applied in scratch only) and, as the control for H3, with the stock settings_2d storage; the 2-D estimator PDE 48 (0.0061794 at 196k today, best 0.0024154 at 2.41M); the 2-D library rows 33-35, 37, 39, 40, 45-47 at light as a regression that nothing certified today is lost.
- C5, 3-D heat (Cor. 35 benchmark, kappa* = 1.5 pi^2) at bench degrees, on Sol: this is where the term count matters most (faces 7 blocks, box 4, product 8, corners 19).

Sets F and offsets o, every one of them built through the row:
| set | rows | o (per-direction rule) | replica of the old runs |
|---|---|---|---|
| none | - | - | psatz off |
| box, N + 1 terms | one entry 3 per row | 1 in the acting direction | o = 0: the same-Gram rows of E1 |
| product, 2^N terms | all nonzero rows over {0, 3}^N | 1 in each acting direction | the 10/08 tensor set |
| faces, 2N + 1 terms | one entry in {1, 2} per row | 1 in the acting direction (variant Fa1) and 0 (variant Fa0) | scalar offset 1 in every direction = today's 'auto' (variant FaS) |
| corners | faces plus two-face rows | 1 in each acting direction | 'corners' of 8659c287 |

Lengths: e = d_3 in {2, 3, 4} for every block (w = 1, 2, 3), lift lengths from the reader (`opts.like` rule) held fixed across sets so that only the weight side moves. On C4 the executive raise loop is replaced by two fixed runs, dw = 0 and dw = 1.

The matched-degree rule, stated once: per direction, the total degree of every weighted term is at most that of the plain term. Box and product at o = 1 meet it exactly; Fa1 is one below and Fa0 one above in the acting direction, so both are run and reported.

## 4. Measurements and verdict rules

Per (case, set, e, o): verdict, reach or gap, and size.
- Verdict: certified only with the eta gate (residual <= 1e-7, every Gram block PSD at the acceptance tolerance of `lpi-feasibility-acceptance-tests`); infeasible only with a dual (Farkas) certificate checked on the dump; everything else UNKNOWN, recorded with the MOSEK status and iterate norms. Every UNKNOWN is re-solved with SeDuMi (E3: the two solvers disagree on the psatz-off programs).
- Reach by bisection with the two-bracket scheme (memory `pietools-bisection-default`), resolution quoted with the result; gap at the fixed margins 0.9 and 0.99 of the exact value where the exact value is known (C1, C2, C3, C5).
- Size: distinct decision variables, the Gram entry count (sum of N^2 per block), rows m, nnz(At), build time, solve time.
- The ladder: for each case and set, the smallest e (and o) at which the target margin is certified, and the size there.

Rules, fixed before the runs:
- "Needed within budget": a set fails to certify the margin at every e within the budget of Sec. 6 where the faces certify it. T6 predicts this can happen only when the faces certify at the top rung e = 4.
- "Useful": at the target margin, the faces give the smallest Gram entry count (and nnz, and time, reported separately) among the certifying configurations. The faces are called useful in N dimensions if this holds on at least two of the three fixed-target cases of that dimension (C1, C2, C3 in 2-D; C5 alone in 3-D, so the 3-D verdict is one case and is said to be).
- "Artefact" (H0 confirmed): the box or product set at o = 1 certifies at the same e as the faces at Fa1 with a Gram entry count within 1.2x. Then the E1-E6 negatives were the degree structure.
- H2: Fa1 against FaS at equal verdict, by Gram entries; H3: on C4, if the term set changes no verdict with the N-D default storage, the set is not a 2-D storage knob and the tier files need not carry it per case.

## 5. Order and schedule

Order follows `test-1d-before-2d`: C0, then C1, then C2 and C3, then C4, then C5. Phase 0 precedes everything; no case runs until its programs pass the support check of Sec. 2.

| phase | content | unit cost (measured) | count | total |
|---|---|---|---|---|
| 0 | instrument and checks | seconds per build | - | half a day of editing |
| 1 | C0, C1 (5 sets x 3 lengths x 2 face variants x 6 margins, 4 plants) | 30-40 s per 2-D solve at m 4-6k (E1); < 1 s in 1-D | about 400 solves | 3-4 h, in 20-min batches while the maintainer is online, or one night of the night runner |
| 2 | C2, C3 (bisections of 16 solves) | 16-60 s per solve (E6); (3,3) KYP 491 s | 5 sets x 3 lengths x 2 cases x 16 | 8-12 h, overnight |
| 3 | C4 (S3, S8, PDE 48, library rows) | 15-100 s per S3/S8 solve; PDE 48 minutes per solve at 196k-2.4M variables | 5 sets x 2 dw x 2 storages x 3 systems + 9 library rows | one night |
| 4 | C5 on Sol | 7-10 min per solve at bench (memory `heatNd-benchmark`); heavy unsolved to date | 4 sets x 2 lengths x 8 bisection steps | 10-12 h as a Sol job array |
| 5 | the ladder to e = 5 on C1 and C2 for box and faces only | as above, growing with e | 2 cases x 2 sets x 2 lengths x 8 | one night |

Every phase writes its table before the next starts; a phase whose unit cost is more than 1.5x the estimate stops and is re-estimated (memory `estimate-runtime-before-running`).

## 6. Budget and the complexity question

Per case the budget is the largest certificate the present tiers build for it: 2-D 1M variables and 600 s per solve on the workstation; 3-D 3M variables and 15 min per solve on Sol. The complexity of the term sets, counted before any solve by the estimator: in 2-D the faces add 4 blocks of the weight Gram size against 2 (box) or 3 (product); in 3-D 6 against 3 or 7; in 4-D 8 against 4 or 15. So from N = 3 on the faces are fewer blocks than the product set and never more than twice the box set; the sets whose count grows faster than N are the corners (4 C(N, 2)) and the full preordering (4^N), neither of which this plan proposes as a default. The plan reports, for each case, the Gram entry count of the smallest certifying configuration of each set, which is the number the tier files need.

## 7. Deliverables

1. One table per case (set x e x o: verdict, gap or reach, Gram entries, m, nnz, times) under `testing_codes/private/.../psatz_faces/` with the scripts and logs, and the verdicts on H0, H1, H2, H3 by the rules of Sec. 4.
2. The recommended F and o per dimension for the tier files: in 1-D the pair (T1); in 2-D and 3-D the set with the smallest certificate at the target margin, with the length at which it certifies.
3. Updates to `psatz_terms_evidence.md` (the new rows) and to `sopvar_settings_proposal.tex` Sec. 8 (the default F, the offset rule, the caps question).

## 8. Risks and how each is handled

- UNKNOWN exits with iterate norms 1e4-1e5 (S3, S17) are solver limits, not cone evidence: every such exit is paired with a SeDuMi solve and the `mincap_dual` degree diagnostic before it enters a table.
- The stock 2-D storage caps (settings_2d, 7/16 and 9/16 pairs) confounded E8: C4 holds the storage fixed and runs both storages, so the term set is the only variable.
- Fa0 exceeds the plain degree and is not nested in the box set at the same e: it is reported as a separate variant, never pooled with Fa1.
- Bisection resolution limits reach claims to the step size; reach is reported with its bracket, as in E8.
- The per-direction offset is a repo change to `poscopvar_direct`: it stays in the scratch shadow until the maintainer approves it, and the shadow is asserted first on the path in every script (memory `pietools-matlab-path-shadowing`).
- Sol jobs end with the session that launched them unless detached by the job array (memory `sol-gpu-job-gotchas`, `detached-runs-die-with-session`).
