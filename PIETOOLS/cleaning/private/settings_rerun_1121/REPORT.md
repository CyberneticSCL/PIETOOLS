# Settings rerun of 11/21/2026 – report

Cases are named by example file and library index. The run log is RUN_LOG.md; the plan is PLAN.md.

**Applied 11/22, committed as 19103941 (not pushed):** examples_1d.patch, examples_2d.patch and dde_library.patch are in the main checkout, re-applied onto the files as comment pass 2 left them (ff7a0ee1). Their settings lines are identical to the verified worktree files, and they pass libcheck there with the same values (libcheck_applied.csv; the DDE09 block through the harness, logs/CHK_DDE09_APPLIED.log). form_39_P, optional_15_primal_Q and optional_41_form_P are not applied.

## 1. What ran

Each row of PLAN.md was searched again with ex_search from the settings its example file (or library block) sets, with three changes from the search of 11/19:

1. **The gate of 11/20.** A solve counts as certified when the harness reads it certified (measure within 0.3 of 1, no numerical error, no infeasibility flag), the executive reports optimal, and eta <= 1e-7.
2. **The fixed-storage check** (post_verify.m). A candidate refused for numerical reasons only goes to a second program. That program keeps the storage operators of the solve fixed (getsol_lpivar_sop on info.aux) and seeks only the slacks. Stage 2 must pass the same gate.
3. **The 3x rule.** The final setting of a row is the smallest certificate accepted directly by the gate, provided it costs at most 3 times the variables of the smallest certificate overall and is not near the memory limit (60% of the size gate: 1.2e6 variables or 27,000 rows). Otherwise the final setting is the smallest certificate.

Alternative starts are tried in order, cheaper and stronger first; the first that certifies is searched. The finite-energy forms 'FQ' and 'FP' come last and run only where every other start fails.

| Queue | Rows | Status |
|---|---|---|
| RR4_P1 | 13 examples in two and three variables, 21 starts | done (§4); interrupted 06:29-16:52 by a planned Windows update |
| RR4_P1D | 17 dual rows of those examples | done, informative (§4) |
| RR4_P2 | 37 one-variable rows | done, 26 min |
| RR4_P3 | alternatives for the 4 one-variable rows with no certified start | done, 5 min |
| RR4_P5 | DDE09 and DDE10 in the form their library blocks run ('P', 'Pd') | done |
| RR4_P6 | DDF03 from the settings of its library block | done |

## 2. The fixed-storage check

**Method.** The storage P and the derivative operator D are read from the solve with getsol, through two operators the stability executive now returns in info.aux. The slack of -D in the solve is rebuilt exactly with poscopvar_direct, from the degree specification, the terms, the kept basis operators and the margin that lpi_ineq_sop recorded. In the direct forms the slack of the margin of P is rebuilt the same way (info.aux.lip, Npm). Only then is -D <= 0, and the margin of P, sought in a program of slack variables alone. If the slack of the solve fails, a slack sized from the support of -D one degree higher is tried. The forms with a storage equality ('Q', 'Qd', 'FQ', H2 'c') are checked only when the solve has eta <= 1e-7.

**Residue.** Coefficients of D that no slack variable reaches were driven to zero by equality rows in the storage variables alone, so they are only as small as the residual of the solve. Fixed, they become rows "0 = residue" that make the slack program inconsistent. Zeroing coefficients by magnitude changed verdicts in both directions:

| Coefficients zeroed below | 36 'Q' | 62 'P' | 41 'Q' |
|---|---|---|---|
| 1e-9 of the largest | passes (one degree up) | passes | infeasible |
| 1e-7 of the largest | passes | infeasible | infeasible |
| 10 x the residual of the solve | infeasible | eta 7.9e-3 | infeasible |

A small coefficient that the slack can reach is a Gram diagonal entry, which supports off-diagonal entries of order sqrt(diagonal x largest). The check therefore zeroes exactly the rows without a slack variable and keeps every other coefficient. It fails if a zeroed value exceeds the equality residual of the solve, which it equals on those rows, or exceeds 1e-7 of the largest coefficient of D. The largest value zeroed in any accepted row is 1e-12.

**Validation** (val_post.m, logs VALPOST10 and SHIFT_*):

| Case | Stage 1 | Stage 2 | PIE shifted by 10 |
|---|---|---|---|
| 36 PIETOOLS_PDE_Ex_2D_Telegraph_Eq, 'Q' | refused (measure 0.37) | passes, the slacks of the solve | infeasible |
| 62 PIETOOLS_PDE_Ex_Road_Field_Berestycki, 'P' | certified | passes, the slack of the solve | infeasible |
| 40 PIETOOLS_PDE_Ex_2D_Heat_Eq_with_ODE, 'Q' | certified | passes | infeasible |
| 50 PIETOOLS_PDDE_Ex, dual H2 'c' | eta 1.7e-6 | not checked | – |
| 41 PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input, 'P' (final, 324) | unknown | passes | infeasible |
| 36, 'Q' (final, 42,308) | refused (measure 4.87) | passes | infeasible |

**Limit.** In the storage forms ('P', 'Pd', 'FP') the check takes P >= eppos I from the Gram construction of the solve, that is, from the returned Gram matrix lying in the cone (which an interior-point iterate does); it does not verify it again.

## 3. One variable

29 of the 37 rows keep the settings their example sets: each converged at its start, certified directly. The full table is table_1d_raw.md. The rows below changed, or need a decision.

**Patched** (patches/examples_1d.patch and patches/dde_library.patch; applied to the main checkout 11/22). Each example file holds one st for all its tests, so each patch passed libcheck with every test it flags (the DDE09 block through the harness, with the identical settings):

| Example | Index | Tests | Variables, before -> after | Values, before -> after |
|---|---|---|---|---|
| PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_in_BC | 9 | stability 'P', dual 'Pd' | 227 -> 187 | – |
| PIETOOLS_PDE_Ex_Wave_Eq_Tip_Damped | 21 | gain 'P', dual 'Pd' | 1,401 -> 1,257 | 0.84190244 / 0.84208439 -> 0.8419069 / 0.84222715 |
| PIETOOLS_PDE_Ex_Heat_Eq_with_Distributed_Disturbance | 22 | gain 'P' | 313 -> 301 | 0.2084336 -> 0.20843359 |
| PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance2 | 49 | H2 'o' -> 'oco', dual 'c' -> 'cco' | 11,847 -> 1,001 | 0.10157583 / 0.10157588 -> 0.10157878 / 0.10158193 |
| examples_DDE_library_PIETOOLS, H-infinity gain Example 1 | DDE09 | gain 'P', dual 'Pd' | 651 -> 623 | 25.922824 / 25.916902 -> 25.92842 / 25.925296 |

- 9: the primal alone reaches 155 (s-length 0 for the upper integral), but the dual needs s-length 1.
- 21: the dual value rises by 1.7e-4 of itself, inside the 1e-3 the search accepts.
- 22: the search also returned the slack lift degree to the stripped value (dD = 1); a first candidate that kept the file's dD = 0 gave 0.50000002 at 105 variables.
- 49: before, the dual at 11,847 variables failed the gate (eta 4.9e-7). In 'o'/'c' no setting tried certifies both tests: the dual 'c' fails or has eta 4.1e-7. The coercive pair 'oco'/'cco' certifies both at 1,001; a scan of the theta-lengths found nothing smaller, and 1,022 variables matches the recorded values to 5e-6. This changes the H2 form of the example, so it needs approval (§6).
- DDE09: the primal final certifies both tests; the dual final (596) does not certify the primal.

**Unchanged after a check:**
- DDE10 (Example A.1): the block st certifies both tests at 1,669 (0.23644325, 0.23653875). The primal alone reaches 869 and the dual 1,285, but neither final certifies the other test.
- DDE11 (Example A.2), DDE25 (H2): only the block st certifies both tests.
- DDF03 (NDS/DDF library, showering people): the block st certifies at 4,673 (0.42054329, the record 0.420543); the search saves 4 variables.
- NDS11: 68. A superseded run reached 57 through an intermediate that only the earlier magnitude cleaning accepted.
- 48 PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DNDN_EdgeObserve: 768,159, unchanged within the build budget.

**Decisions** (optional patches, not applied):
- 15 PIETOOLS_PDE_Ex_Timoshenko_Beam_2. The example runs the dual 'Qd', which fails the gate at its settings (eta 6.0e-7 at 74,593 variables; on 11/18-19 infeasible at light, heavy and extreme). The primal 'Q' with deg_Q{1} = [0 1 1] certifies at 3,229 variables, directly and through libcheck. patches/optional_15_primal_Q.patch moves the example to the primal test.
- 41 PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input. 'Q' fails at light: its margin of P is infeasible with exactly the slack of the solve, whose target matches the solve to 7e-10. 'P' at 324 variables passes only the fixed-storage check, so the library path prints unknown. patches/optional_41_form_P.patch.

## 4. Two and three variables

RR4_P1 ran 06:12 – 19:15, with a planned Windows update restart from 06:29 to 16:52 (resumed on the maintainer's word). Two fixes were made in the course of the run:
- the starts measured above the size gate were dropped: for 54, light 7.0M and heavy 29.2M variables; for 55, light 28.0M, and heavy, which is larger;
- ex_search now keeps the starts a row has tried, since a resumed row retried every start, so a start that filled a session repeated forever.

| Example | Index | Result | Path |
|---|---|---|---|
| PIETOOLS_PDE_Ex_2D_Telegraph_Eq | 36 | 'Q' 52,880 -> 42,308; no direct certificate on the path; shift check passed | stage 2 |
| PIETOOLS_PDE_Ex_2D_NS_Antonelli | 39 | 'P' with the 2-D storage of the 'light' level and dD = 0: 334,656. The file runs 'Q' at light, which fails, as on 11/20; since 'P' is the case Q = P*T of 'Q', the 'Q' settings tried exclude that point (§4a) | direct |
| PIETOOLS_PDDE_Ex | 50, H2 and dual | no start certifies | – |
| PIETOOLS_PDE_Ex_3D_Reaction_Diffusion_Eq | 53 | 'P': the search from the light tier reached 584,704; the st of the file (stripped) certifies at 580,612 and stays | direct |
| PIETOOLS_PDE_Ex_3D_Advection_Diffusion_Eq | 54 | 'P' 1,835,012 (the example st, near the memory limit; a 19-minute solve); the search reached 1,835,009 | direct |
| PIETOOLS_PDE_Ex_3D_Coupled_Heat_Eq | 55 | the example st builds 7,340,048 variables, above the size gate | – |
| PIETOOLS_PDE_Ex_2D_Delayed_Heat_Eq_Fridman | 56 | 'Q' light dD 0 is infeasible with a certificate; the example st builds 2,524,300 variables, above the gate | – |
| PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Fridman | 57 | 'Q' 23,121 -> 19,840 | direct |
| PIETOOLS_PDE_Ex_Road_Field_Berestycki | 62 | 'P' 83,269, unchanged (three stage-2 passes on the path; the direct start kept by the 3x rule) | direct |
| PIETOOLS_PDE_Ex_Heat_Edge_Cascade_Variable_Coeff | 63 | 'Q' 149,054, unchanged | direct |
| PIETOOLS_PDE_Ex_2D_Linear_Turing_Pair | 65 | 'P' 1,635,656 -> 250,280 | direct |
| PIETOOLS_PDE_Ex_Transport_Edge_Diffusion | 72 | 'Q' 256,948, unchanged | direct |

Stage 2 in RR4_P1: 7 passes, 1 failure, 0 errors.

**Patched** (libcheck_rr_C.csv; 57 and 65 applied to the main checkout 11/22):
- 57 PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Fridman: Q of degree 1 in every direction, slack weights in the directions a space does not own, the slack integral bases capped at Dmax-1. Library path: optimal at 19,840 (23,121). patches/examples_2d.patch.
- 65 PIETOOLS_PDE_Ex_2D_Linear_Turing_Pair: no slack weight raise of the feasibility passes on the spaces in one and two variables (dw_feas = [0 0 1]). Library path: optimal at 250,280 in 80 s (1,635,656 in 467 s). patches/examples_2d.patch.
- 39 PIETOOLS_PDE_Ex_2D_NS_Antonelli: the form 'Q' -> 'P', with the 2-D storage of the 'light' level and dD = 0. Library path: optimal at 334,656. This changes the stability form, so it is in its own file for approval: patches/form_39_P.patch. Not applied: see §4a.

**4a. 'Q' contains 'P'** (maintainer, 11/22: "this shouldn't happen, since P is a special case of Q"). For Q = P*T, sym(T'*Q) = T'*P*T >= eppos2*T'*T whenever P >= eppos2*I, and A'*Q + Q'*A = T'*P*A + A'*P*T. Every 'P' certificate is therefore a 'Q' certificate whenever the 'Q' program can represent Q = P*T and its two slacks. A 'Q' failure where 'P' certifies is a property of the 'Q' settings or of the 'Q' construction, not of the example. 'P' is not the only certificate of 39; it is the only one the settings tried reach. q_contains_p.m locates the exclusion. It fixes the storage of the 'P' solve and checks Q = P*T against (a) the span of the decision Q, (b) the slack of the margin of sym(T'*Q), and (c) the slack of -D, each as the 'Q' executive sizes them. Logs: logs/QCP_*.log, SPAN_39.log, Q39_DEG3.log.

| Check | 39 PIETOOLS_PDE_Ex_2D_NS_Antonelli (light) | 41 PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input (light) |
|---|---|---|
| (a) P*T in the span of Q | no at deg_Q{2} = 1 (coefficients up to 5.65e-7 of max 1.79 outside); also no at 2, [2 1 1], [2 2 1], [2 2 2] (1.1e-7); yes from 3 | yes at deg_Q [4 2 2] |
| (b) margin, slack of the 'Q' solve | certified (fr 0.99, eta 6.6e-9) | measure 2.04, eta 5.7e-9: refused |
| (c) -D, slack of the 'Q' solve | certified (fr 1.13, eta 8.2e-9) | measure 1.55, eta 2.5e-8: refused |

**39: the degree of Q excludes the 'P' point.** At the light tier Q has degree 1 in two variables, and P*T needs degree 3. The search of 11/19-11/20 tried deg_Q{2} only up to [2 2 1]. The two slacks of the degree-1 'Q' program (117,072 and 207,936 variables) hold the 'P' point with an interior.

At deg_Q{2} = 3 the 'Q' program builds 31,853,664 variables (21,945,744 with dD = 0), far above the size gate of 2e6. The executive sizes both slacks from the supports of the degree-3 decision operators. A 'Q' program that contains the 'P' point at a size the workstation solves would therefore need Q of degree 3 with slacks sized as at degree 1. The executive has no option for that.

**41: the slacks hold the 'P' point only on the boundary of their cones.** The solver reaches the 'P' point to a residual near 1e-8, but not with a measure the gate accepts. Adding 1e-4 of the target's scale times T'*T makes both (b) and (c) certify (fr 1.08, 1.13). Two controls: T'*T alone certifies on the margin slack, and the slack of the 'P' solve (172 variables, sized like P) certifies -D.
- For -D, a slack sized like the generalized storage sym(T'*Q) contains the 'P' point with an interior: certified at weight raise 0 (fr 0.97, 1,280 variables) and 2 (fr 1.16), not at 1 (fr 1.58). The executive sizes the -D slack of 'Q' from the support of D, dropping opts.like (PIETOOLS_stability_sop.m:366), so no setting reaches that sizing.
- For the margin, every slack tried holds T'PT - eppos2 T'T only on the boundary: the slack of the solve, like sym(T'*Q) at raises 0-2, and a margin of eppos2/10 or 0. The room eppos2*T'*T is about 1e-6 of the target's scale, and the room test needs about 1e-4. Scaling P does not change that ratio.

So for 41 no 'Q' setting tried contains the 'P' point with an interior. The exclusion is in the construction of the 'Q' slacks, which are sized from the supports of the generic decision operators.

**Duals (RR4_P1D, 17 rows).** Every two- and three-variable example has its dual test commented out, so these rows inform and patch nothing. The dual form is 'Qd' (direct, posed in two or more variables since 905b0d26) or 'Pd' where the example runs 'P'.

| Example | Index | Form | Variables, start -> final | Path |
|---|---|---|---|---|
| PIETOOLS_PDE_Ex_2D_Transport_Eq | 33 | Qd | 89, unchanged | direct |
| PIETOOLS_PDE_Ex_2D_KISS_Model | 34 | Qd | 72,081, unchanged | direct |
| PIETOOLS_PDE_Ex_2D_Reaction_Diffusion_Eq | 35 | Pd | no start certified | – |
| PIETOOLS_PDE_Ex_2D_Parabolic_Eq | 37 | Qd | 11,040, unchanged | direct |
| PIETOOLS_PDE_Ex_2D_Wave_Eq | 38 | Qd | 380,804 -> 42,308 (the direct start costs 9.0x) | stage 2 |
| PIETOOLS_PDE_Ex_2D_Heat_Eq_with_ODE | 40 | Qd | 8,518, unchanged | direct |
| PIETOOLS_PDE_Ex_Coupled_ODE_Heat_Eq_w_Delay | 42 | Qd | 69,915 -> 43,065 | direct |
| PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay | 43 | Qd | 1,008,600 -> 196,425 (direct start 5.1x) | stage 2 |
| PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Fridman | 57 | Qd | no start certified | – |
| PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Riemann_Fridman | 58 | Qd | 586,772 -> 150,400 (direct start 3.9x) | stage 2 |
| PIETOOLS_PDE_Ex_Wave_Delayed_Boundary_Damping_Xu | 59 | Qd | 136, unchanged | direct |
| PIETOOLS_PDE_Ex_Heat_Delayed_Robin_Boundary_Nicaise | 60 | Qd | 552, unchanged | direct |
| PIETOOLS_PDE_Ex_Heat_Laplace_Beltrami_Edge | 61 | Qd | no start certified | – |
| PIETOOLS_PDE_Ex_Road_Field_Berestycki | 62 | Pd | 427,249 -> 101,428; no direct certificate on the path | stage 2 |
| PIETOOLS_PDE_Ex_Heat_Edge_Cascade_Variable_Coeff | 63 | Qd | 652,407 -> 32,268; no direct certificate on the path | stage 2 |
| PIETOOLS_PDE_Ex_Two_Layer_Heat_Rectangles | 64 | Qd | 288,400 -> 219,044 | direct |
| PIETOOLS_PDE_Ex_Age_Structured_Diffusion_Kang | 66 | Qd | 210,361 -> 210,281 | direct |

Stage 2 in RR4_P1D: 14 passes, 3 failures, 0 errors. Shift check of the five rows accepted through stage 2: each passes at 0 and fails at 10 (logs/SHIFT_P1D_*.log).

In 38 (the wave equation, conservative, epneg 0) the derivative is zero to the residual of the solve (max|D| 8e-9), so the margin of P carries the certificate. Stage 1 of 58 reads 'infeasible' on the harness measure (-0.97), as 63 does below.

Two notes on 63:
- Its stage 1 reads 'infeasible' on the harness measure (-0.9999). The solver raised no infeasibility flag and the executive reported unknown, which the plan counts as a numerical refusal.
- Its storage is badly scaled (max|D| 4.6e5). Its stage-2 eta of 2.2e-8 is therefore an absolute residual near 1e-2: the limit of every certificate accepted on eta, largest here.

## 5. Defects found in the setup of the rerun

Each was found by a check of the run itself and fixed before the queues whose results are reported (RR4):

- rr_run.sh called ex_search without lb_setup(), so the first queues ran on the saved path of the main checkout (RR_, RR2_ invalid).
- example_st took the second output of files that return [PDE_t, PDE_b, st].
- The stage-2 slack was first built unpruned (7.7M variables on case 62); it is now the slack of the solve, rebuilt exactly.
- Residue cleaning by magnitude (§2); the descent size gate refused stage-2 programs (caps lifted during the check); the slack of the margin of P was sized from the numeric P (case 41).
- make_starts.pl took the Phase A form for DDE rows; the blocks of DDE09 and DDE10 run 'P' (RR4_P5). DDF03 started from the light tier instead of its block (RR4_P6).
- st_patch.m first diffed against the tier, which omits a field the search returned to its tier value (case 22); patches are now the diff from the file's own st.
- ex_search retried every start of a row with no certified start in each resumed session, so a start that filled a session repeated forever (case 55, 28M variables, 13.6 min to build); it now keeps the starts tried (st.altdone).

## 6. Decisions and approvals

1. 15: primal 'Q' in place of the dual test (optional patch).
2. 41: form 'P', certified only by the fixed-storage check (optional patch).
3. 49: H2 form 'oco' (dual 'cco'), part of examples_1d.patch.
4. 39: stability form 'P' in place of 'Q' (patches/form_39_P.patch). Not applied: 'Q' contains every 'P' certificate (Q = P*T), so the question is which part of the 'Q' program excludes this one (§4a).
5. 36: its final setting (42,308 variables) is certified only by the fixed-storage check, as is 41; no patch is written for it.
6. The fixed-storage check lives in the harness. Certificates accepted only through it (36, 41) print as unknown on the library path. Showing them as certified there needs the check in the toolbox: an executive option, or a routine that takes prog and info.
7. Committed 11/22: the applied example patches (9, 21, 22, 49, 57, 65, the DDE library) as 19103941, not pushed. Not committed: the instrument info.aux in PIETOOLS_stability_sop (Pm, Dm, Nm, li, form, Npm, lip) and PIETOOLS_H2_norm_sop; the harness files post_verify.m, post_shift.m, st_patch.m, val_post.m, the ex_search and lb_exec changes; this folder. 905b0d26 is on origin (pushed with 349b247f).
