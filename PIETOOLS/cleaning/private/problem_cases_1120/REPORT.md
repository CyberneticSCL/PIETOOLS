# Problem cases – report of the run of 11/20/2026 (05:10–15:22; final)

Maintainer, 11/20/2026: "create and run a plan for trying different settings variations and trades to certify these
problem cases. Run autonomously for 11 hours"; "also try cross-terms if not already"; "also try higher degree psatz
terms on 69 and 70". Plan: PLAN.md; chronology: RUN_LOG.md. Runs on a clean worktree at aad9606c; every certified
problem case rerun on 41c9210c and on f474e7de (the release comment pass, HEAD of the main checkout since 05:33) with
identical status, measure, size and executive verdict, through the harness and the library path.

Certified means both gates of the harness: MOSEK optimality measure (MSK_DINF_INTPNT_OPT_STATUS, the feasibility
ratio sossolve reports) within 0.3 of 1 without a numerical-error flag, and the executive status optimal; eta is the
normalized residual of the executive.

## 1. Outcome by case

| Case | Outcome | Setting or finding |
|---|---|---|
| 36, 2-D telegraph | **certified** | light 'Q', dD 0, wcross 1, epneg 1e-3: 52,880 variables, measure 1.098, eta 1.35e-8; library path optimal. At epneg 0: deg_Q{2} [0 1 1], 380,804, measure 1.189, eta 2.37e-8 |
| 41, heat, delayed Dirichlet boundary input | **certified in the file's own explicit-delay form** | expand_delays = true (a branch of the file) and light 'Q' with eppos 1e-2: 2,496 variables, measure 1.185, eta 1.96e-10; library path optimal. As written (the parser expands the delay): weakly infeasible at every setting tried |
| 44, wave, partly delayed boundary damping | **not certified as written at any setting tried; certified in the minimal position form ('Pd') and in energy states** | as written (ODE state, u, u_t, delayed u and u_t): infeasible in 'Q', 'P', 'Pd' and 'Qd'. Minimal form (u, u_t, delayed u_t – the same dynamics without the auxiliary states): infeasible in 'Q' and 'P'; certified in 'Pd' at epneg 0 (2,520 variables, measure 1.000, eta 9.3e-10; identical on 41c9210c), at epneg 1e-2 with light (2,520, measure 1.274, eta 6.0e-8) and at epneg 1e-3 with heavy (14,373, measure 0.844, eta 8.9e-9); at epneg 1e-1 with light 2,520, measure 0.987; library path optimal (worktree copy of the example in the minimal form, stability_dual, epneg 1e-1). Riemann states (w_s, w_t, delayed w_t): certified in 'P' (2,520, measure 1.001; the system of example 59), 'Pd' and 'Qd' |
| X3, the system of 44 on two variables | as 44 | as written infeasible (one row unknown, measure −0.50). Minimal form with the delay on s2: certified in 'Pd' at epneg 0 – 9,000 variables, measure 0.9995 (mu 0.375) and 0.9941 (mu 0.4); at epneg 1e-3 inside the measure gate with eta 2–3e-7, above 1e-7; at epneg 1e-2 certified (9,000, measure 1.001–1.002, eta below 1e-9); library path optimal (worktree copy, stability_dual, epneg 1e-2). Riemann states certify: 'Q' 5,214, 'P' 9,000 (mu 0.375 and 0.4) |
| 39, 2-D Navier–Stokes | **certified in 'P'** | light 'P' with the stock 2-D storage and no slack degree raise, at either stock level: st.settings_2d = settings_2d_stock('light'); st.lpi.ineq.dD = 0 – 334,656 variables, measure 1.12, eta 4.9e-9, 64 s solve; or settings_2d_stock('heavy') with dD 0 – 335,601 variables, measure 1.015, eta 4.39e-8, 52 s. Each identical on 41c9210c and optimal through the library path (form 'P' and the two statements in the worktree copy of the example). With dD 1 (heavy), 1,692,801 (measure 0.9998, eta 5.5e-11, 790 s). Every 'Q' variation: optimal, eta about 1.3e-8, measure 2.4–7.3 (refused by the measure gate) |
| 70, 2-D wave, acoustic boundary | **not certified as written at any setting tried; certified in a mixed position–strain form** | as written: infeasible with certificate at every setting tried ('Q', 'P', 'Pd'; to 1.8M variables). Mixed form – ph with its s1-derivatives only (order (2,0)), the strain e2 = ph_s2 a state (order (0,1)), the acoustic coupling a Dirichlet datum e2(s1,0) = −et – certified in 'Q' light with st.lpi.ineq.dD = 0 at epneg 0: 187,116 variables, measure 1.034, eta 7.1e-10, pass 1 (with deg_Q{2} [0 1 1]: 186,809, measure 0.9996); identical on 41c9210c; library path optimal (worktree copy of the example in the mixed form, pc_lib70_edit.pl) |
| 69, 2-D wave–heat | **not certified as written at any setting tried; certified in a mixed position–strain form** | as written: infeasible with certificate at every setting tried ('Q', 'P', 'Pd'; to 1.7M variables); without the heat interface certified (928,900). Mixed form – u with its s2-derivatives only (order (0,2)), the strain e1 = u_s1 a state (order (1,2), e1 = 0 at s2 = 0, 1), the heat-flux coupling a Dirichlet datum e1(0,s2) = −w_s1(0,s2) – certified in 'Q' light with st.lpi.ineq.dD = 0 at epneg 0: 267,060 variables, measure 0.722 (inside the gate, 0.022 from its edge), eta 4.3e-8, pass 1; identical on 41c9210c; library path optimal (worktree copy of the example in the mixed form, pc_lib69_edit.pl) |
| 56, 2-D heat with a 3-D delay state | **not certified** | infeasible at every local setting tried (266,148–1,587,084 variables; 'Q', 'P', 'Pd'; the delay state of order 0 in space); the degree that certifies its 1-D twin 43 (deg_Q{3} [1 2 2]) gives 7,857,280 variables, m 135,522; the default degree fails at every delay tried (0.25–0.9 tau*) |
| 50 dual, H2 | not certified | 'cco' infeasible with certificate at every start (60,645 to 1.66M); 'c' unknown at every stripped start (66,923–887,267; measure −0.38 to −0.28); 'c' light/heavy 5.2–6.0M (above the gate) |
| 55, 3-D | not certified: no solve finished | stripped 'P' (7,340,048 variables, m 105,717) on Sol at 64 threads: cancelled at the 6 h limit inside the solve, 191 GiB peak memory (job 2067511); locally dD 0 infeasible at slack weight 1–3 (114,704–684,304), and the caps give 1,808 (infeasible) or the full program |

## 2. Findings that hold beyond the settings

1. **The margins are a lever the tiers never move.** Every tier sets eppos = 1e-4 and epneg = 0. A positive epneg
   (exponential decay, a stronger claim) or a larger eppos brought 36 and 41 inside the measure gate where every degree,
   term and slack variation left them inconclusive.
2. **The measure gate decides 36 and 39.** For MOSEK, sossolve reports MSK_DINF_INTPNT_OPT_STATUS as the feasibility
   ratio; lb_exec certifies only within 0.3 of 1. On 36 (epneg 0) and 39 the solver reports optimal with eta near
   1e-8 while the measure lies at 0.375 (36) and 2.4–7.3 (39); neighbouring settings of 36 move the measure between 0.36
   and 7.0 at equal eta. The gate policy is D93, the maintainer's.
3. **Formulation, not settings, blocks 41, 44 and X3** (and 69, 70: item 7).
   - 41: the delayed input enters the heat equation as a Dirichlet value. The same system with the input as a flux,
     or decoupled, certifies in 'Q' and 'P' at 2,496–4,000 variables; with the delay as an explicit transport state
     (the file's expand_delays = true) the Dirichlet form certifies too (Section 1). The parser's own delay
     expansion gives a larger program (19,776 in 'Q') that no setting resolves.
   - 44 and X3: as written, every form fails ('Q', 'P', 'Pd', 'Qd'). Without the auxiliary ODE state and delayed
     position (the minimal form u, u_t, delayed u_t) the primal forms still fail and the dual form 'Pd' certifies;
     the energy states certify in every form.
   - Reading for 41 (inferred, not proven): a Lyapunov functional for a heat equation driven at a Dirichlet boundary
     needs the boundary value of the driving state, which a PI operator on the L2 state cannot express; the solver
     then sees a weakly infeasible program (primal residual to 1e-10, dual objective settling at 0.145 > 0).
4. **The parser's delay expansion against an explicit transport state** (41): same system, 19,776 against 2,496
   variables in 'Q', inconclusive against certified. A question for the parser (Declan), not changed here.
5. **In one variable, the dual form certifies the position state** (DIAG69). The 1-D twins of 69 (u, v, w along s1)
   and 70 (ph, ps and the edge states as ODEs along s2) are infeasible with certificate in 'Q' and 'P' at light and
   stripped and certified in 'Pd' (2,520 and 1,356 variables, measure 1.022 and 1.003); their energy-state forms (68;
   e, v, de, et) certify in every form. PDE 18 (1-D damped wave in u) and the minimal forms of 44 and X3 (DIAG44D,
   DIAG44E) show the same pattern. The dual forms had never run on a problem case before this run. In two variables
   the pattern reverses: 'Pd' is infeasible on the uncoupled wave of PDE 38, which 'Q' certifies with the velocity
   of order (2,2), (1,2) or (2,1) alike (DIAG38), on 69, 70 and 56 (W2VS), and on the PDE 50 system, which 'P'
   certifies (W50S). On self-adjoint PIEs (the 2-D heat equation) 'Pd' and 'P' build the same programs with the same
   results (DIAGPD): these failures are a property of the dual form on non-self-adjoint systems, not of its
   implementation. Measured (DIAG2C, item 7): in two variables 69 and 70 fail on their couplings – the heat interface
   of 69, the edge states of 70 – not on the position state or the order of the velocity (70 with the velocity of
   order (2,2) also fails; DIAG70).
6. **Two settings fields do not act on the stability executive.** PIETOOLS_stability_sop passes its own dw to
   lpi_ineq_sop (base 0 in form 'Q'), so st.lpi.ineq.dw and st.lpi.ineq.wmin are not read there; its weight knobs are
   st.sop.dw_feas, dw_P, dw_Q and the raises max_raise_feas, max_raise_Q. Every harness row reports the last pass of
   the executive (for 'Q' light in two variables, pass 2: the slack weight raised by one).
   By design since the floor was added (11/10): a dw in the call switches off the tier dw and wmin (lpi_ineq_sop,
   the fill of dD, dw, wmin; its header on wmin), and every _sop executive passes dw, so the two fields reach only
   programs that call lpi_ineq_sop without dw. The tier file (lpisettings_sop) describes lpi.ineq as the slack
   options of an operator inequality and does not state this exception; neither does the executive header.
7. **A boundary coupling of a wave in the position state defeats 'Q' and 'P'; a strain state removes it.** Measured:
   in one variable, the wave u_tt = u_ss with u(0) = 0, u_s(1) = −u_t(1) in the states u, u_t is infeasible in 'Q'
   (to deg_Q{1} [8 6 6], 155,457 variables) and 'P', and certified in 'Pd' at 1,120 (DIAG44F). In two variables, 70
   without its edge states certifies in 'Q' (Neumann at s2 = 0; 380,880), 70 with the damping datum ph_s2 = ps at
   s2 = 0 does not, 69 without its heat interface certifies (928,900), and 70 with the strain e2 = ph_s2 as a state
   – the coupling a Dirichlet datum of e2 – certifies at 187,116 (DIAG2C, DIAG2M). By hand for the 1-D wave (corrected
   11/21): its energy E = (|u_s|^2 + |v|^2)/2 carries the boundary term −u(1)v(1), so no linear Q into L2 gives
   <Tx, Qx> = E (restricted to states with v = 0, the functional u ↦ u(1) is unbounded in L2). The boundary value is
   itself a PI operator on the fundamental state, u(1) = subs(T_u,s,1) x, and E = <Tx, Qx> + u(1) (Q_b x) exactly,
   with Q x = (−x_u, T_v x)/2 and Q_b = −subs(T_v,s,1)/2. The obstruction is the parametrization T'Q, not the PI
   family: a storage with the added term subs(T_u,s,1)' Q_b contains the energy (not built or tested). The earlier
   wording, "outside the PI family", was wrong (maintainer, 11/21).
8. **The rewritten forms pass soundness checks** (late checks, f474e7de). Each rewrite is a new PIE, so each was
   tested where it must fail:
   - 70 (mixed form): at epneg 1e-3 and 1e-2 infeasible (1,288,544) – 70 is not uniformly stable; with the edge
     damping reversed, et_t = +c*et − de − ps(s1,0), infeasible at c = 1 and 0.1 (DIAG70S, SND70).
   - 69 (mixed form): at epneg 1e-3 and 1e-2 infeasible at the certified program (267,060) – 69 is not exponentially
     stable (SND69); with the heat part made unstable, w_t = w_s1s1 + w_s2s2 + 30w, infeasible at the certified
     program (267,060, pass 1; pass 2 above the m gate, 1,926,185, m 45,866) while the example (r = 0) certifies
     (DIAG69S).
   - 44 (minimal form, 'Pd'): for k = 1 and tau = 1, with z = exp(lam) and u = sinh(lam s), the characteristic
     equation reduces to (2 − mu) z^3 + mu z^2 + mu z − mu = 0, Re lam = ln|z| (derived by hand; it reproduces the
     Riemann witness, −0.289 at mu = 0.4 and −0.0163 at mu = 0.65). The complex pair reaches |z| = 1 only at
     mu = 2/3: there the real root equals the product of the roots, mu/(2 − mu), and 2w^2 + w − 1 = 0 with
     w = mu/(2 − mu) gives w = 1/2; the real root never reaches |z| = 1 (z = 1 and z = −1 give 2 and −2). So at
     tau = 1 the spectral abscissa ln max|z| is negative for mu < 2/3 and positive for mu > 2/3 – the bound
     mu < 1/2 of the example holds for every delay. The dual certificate holds for mu = 0.45, 0.48, 0.49, fails
     from 0.55 (conservative to 0.65) and at 0.7, 0.8, 0.9 (unstable): no certificate past the boundary (DIAG44G,
     DIAG44H).
   Heavier tiers do not extend the reach: at mu = 0.55, 0.6, 0.65 the minimal form in 'Pd' and the Riemann form in
   'Pd' and 'P' are infeasible at heavy and veryheavy (DIAG44I) – the certificates cover the delay-independent
   bound mu < 1/2, none of the delay-dependent range 1/2 < mu < 2/3.
   - 70 (mixed form) certifies at every edge-oscillator setting tried, (damping, stiffness) = (0.1, 1), (10, 1),
     (1, 10), (1, 1): 187,116 variables, measure 1.002–1.138 (DIAG70P).
   - A second solver, SeDuMi 1.3 at its defaults, on the same programs: 44 and X3 ('Pd') and 70 (mixed form)
     reach feasibility ratios 0.98–1.00 with eta 5.6e-10 to 5.3e-7 and status inaccurate (DIAGSED); 36 and 69
     (two-variable 'Q') diverge (eta 4e-2 to 6e-2 in pass 1; SED36, SED69). No infeasibility certificate on any of
     them; SeDuMi alone passes the harness gate on none.
9. **What a certificate of a rewritten form says about the example as written** (added 11/20 evening). The
   rewrites change the state variables, not the PDE.
   - 69, 70: let J append to a state of the example the strain of the coupled direction (70: e2 = ph_s2; 69:
     e1 = u_s1). For every classical solution z of the example, J z solves the mixed form: each equation and boundary
     condition of the mixed form is one of the example, or follows from one by differentiation (e2_t = ps_s2;
     e2(s1,0) = ph_s2(s1,0) = −et), and the mixed form omits the Dirichlet condition of the position in that
     direction, so it admits more solutions, not fewer. A 'Q' certificate of the mixed form (V = <x', T'Q x'>
     non-increasing on every solution, V >= eppos2 ||T x'||^2) therefore gives, for every solution of the example,
     ||z(t)|| <= ||J z(t)|| <= C ||x(0)||, with x the fundamental state of the example (the lower derivatives in
     J are PI operators of x): the notion the 'Q' test of the example itself certifies.
   - 44, X3: every classical solution of the example restricts to a solution of the minimal form (drop the ODE
     state x1 and the delayed position x4; the remaining equations and boundary conditions appear in the example
     verbatim). The certificate bounds (u, u_t, delayed u_t). x1 and x4 are driven by u(1,t) and feed nothing
     back; their decay needs a bound on u_s (u(1) = int_0^1 u_s), which depends on the norm of the duality result
     behind 'Pd' – not checked.
10. **The finite-energy forms, and the dual H2 test without its margin, certify none of the open cases** (11/20
   evening; pc_fe.m on 17d68747, results/fe_17d68747/, logs_fe/). 'FP' (P >= 0, T'PA + A'PT + 2 epneg T'T <= 0)
   and 'FQ' (P >= 0, T'Q = P, A'Q + Q'A + 2 epneg T'T <= 0) certify int_0^inf ||Tx||^2 dt <= C ||Tx(0)||^2. Both are homogeneous in (P, epneg), so
   feasibility does not depend on the value of epneg > 0; the runs use epneg = 1, since at the default 1e-5 the
   operator P = 0 violates the inequality by less than the solver tolerance. For 50, 'cco' at eppos 0 (W >= 0)
   admits the controllability Gramian, which W >= eppos I excludes for a single input (the Gramian is compact);
   the bound then rests on the 'P' certificate of the system (98,080 variables, W50S).

   'FQ' is built for one-variable PIEs only: for a PIE in two or more variables, PIETOOLS_stability_sop poses 'FP'
   in its place (line 146), with no warning, and prints the label 'FQ' (set at line 134, before the substitution).
   Every open case is such a PIE – 43 and 50 have two variables (space and delay), 55 and 56 three – so 'FQ' was not
   tested on any of them: the 43 'FQ' rows built the 'FP' program (the same variable and row counts). 'FP' is the
   case Q = PT of 'FQ', so 'FQ' is feasible wherever 'FP' is (at matching degrees of Q).

   | Case | Test | Variables | Result |
   |---|---|---|---|
   | 55 | 'FP' stripped, dD 0 | 114,704 | infeasible with certificate |
   | 43 (one space variable and the delay) | 'FP' light, dD 1 and dD 0 ('FQ' requested, posed as 'FP') | 428,185; 91,145 | infeasible with certificate |
   | 43 | 'FP' heavy ('FQ' requested, posed as 'FP') | 3,497,665 | above the gate, not solved |
   | 56 | 'FP' stripped, dD 0, at 0.9 and 0.5 tau* | 258,409 | infeasible with certificate |
   | 50 dual | 'cco' at eppos 0: stripped with the cap, stripped, light (pass 1, pass 2) | 16,337; 96,417; 503,633; 1,660,065 | infeasible with certificate |

   The margin of W is therefore not what blocks the dual H2 test of 50.
11. **A posteriori check of the gate-rejected solves: the storage fixed, only the slack sought** (11/21; maintainer:
   "fix the storage output of the previous result and just find the slack", and of the projection check, "you can
   try it"). pc_post.m on 17d68747 in a worktree (C:/wt1121v, removed; the two instruments are in
   results/post_17d68747/: info.aux returned by the executives, and a sossolve copy that keeps the SDP data);
   results in results/post_17d68747/, logs in logs_post/.
   - Stage 1: the row as the harness ran it. Stage 2: the storage fixed – P in 'P' and 'Pd'; P and Q in 'Qd'; Q in
     'Q'; R, Q, W and gam in the H2 form 'c' – and only the slacks sought, with the degree specification and terms
     of the stage-1 slack, unpruned. Pruning reads the nonzero cells of the numeric target and kept 1,352 of the
     1,568 stage-1 entries on a control; unpruned, the stage-2 slack contains the stage-1 slack. The margin of the
     storage is checked by a slack of its own: P − eppos I >= 0 in 'P' and 'Pd', P − eppos2 T T' in 'Qd',
     P − eppos2 T'T in 'Q'.
   - Before it is fixed, the storage is projected onto the stage-1 rows that involve it alone: coefficients of D
     that no slack monomial reaches, and T'Q = Q'T or T Q = P. Its Gram blocks alternate with the PSD projection.
     Fixed off those rows, D keeps components outside the span of the slack, and sossolve declares such a row
     infeasible before the solve (the first control run).
   - Stage 2 passes when it passes the harness gate: measure within 0.3 of 1, numerr 0, executive status optimal,
     eta <= 1e-7.
   - Controls: two certified rows pass stage 2 – the explicit-delay form of case 41 in 'Q' with eppos 1e-2 (measure
     1.118, eta 1.6e-8) and 44min in 'Pd' at epneg 0 (0.993, 5.8e-9). The 44min storage against the PIE at
     mu = 0.8 (unstable; the bound is 2/3), with a slack sized from the support of its own derivative (44,649
     entries), is infeasible.

   Systems: case 41 is PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input with its delay as an explicit transport
   state (pc_diag41b.m); case 36 is PIETOOLS_PDE_Ex_2D_Telegraph_Eq; 69mix2 is
   PIETOOLS_PDE_Ex_Wave_Heat_Transmission_2D_Batty (case 69) in the mixed position–strain form (pc_diag69m.m);
   44min is the minimal position form of PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay (case 44) at mu = 0.4, and
   x3min the same with the delay on s2, the layout of PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D (pc_diag44e.m);
   case 50 is PIETOOLS_PDDE_Ex, the dual H2 test. Cases are indices of
   PIETOOLS_examples/examples_PDE_library_PIETOOLS.m.

   | Row | Stage 1: status, measure, eta, variables | Storage repair | Stage 2: status, measure, eta | Verdict |
   |---|---|---|---|---|
   | case 41, 'Q' light | unknown, 0.340, 2.0e-10, 2,496 | none (42 rows already exact) | certified, 1.073, 1.5e-8 | passes |
   | case 36, 'Q' light dD 0 wcross 1, epneg 0 | unknown, 0.375, 1.8e-8, 52,880 | Q onto 1,607 rows | certified, 1.121, 1.0e-8 | passes |
   | case 36, as above with epneg 1e-3, deg_Q{2} [0 1 1] | unknown, 1.524, 7.0e-7, 52,804 | Q onto 1,571 rows | certified, 0.788, 5.1e-9 | passes |
   | 69mix2, 'Q' light dD 0, deg_Q{2} [0 1 1] | unknown, 10.71, 7.5e-8, 266,756 | Q onto 6,182 rows | certified, 1.078, 1.3e-8 | passes |
   | 44min, 'Qd' light | unknown, 0.560, 9.5e-9, 3,287 | P, Q onto 628 rows | certified, 0.934, 1.2e-8 | passes |
   | 44min, 'Pd' stripped, epneg 1e-3 | unknown, 1.316, 2.1e-7, 2,520 | P moved: its Gram block λ_min −4.4e-3 after the repair | certified, 1.061, 9.3e-9 | passes, on the repaired P |
   | 44min, 'Pd' light, epneg 1e-4 | unknown, 1.552, 1.3e-8, 2,520 | 74 rows; after the repair the Gram block of P has λ_min −0.31 | infeasible (measure −0.91) | fails |
   | x3min mu 0.375, 'Pd' heavy, epneg 1e-3 | unknown, 0.027, 2.0e-5, 63,621 | Gram λ_min −0.69 after the repair | infeasible | fails |
   | x3min mu 0.4, 'Pd' heavy, epneg 1e-3 | inaccurate, 1.042, 4.7e-5, 63,621 | Gram λ_min −0.13 after the repair | infeasible | fails |
   | case 50, 'c' stripped dD 0, cap_inc 'auto' | unknown, −0.383, 1.7e-6, 66,923 | R, Q, W, gam onto 19,863 rows | infeasible (−1) | fails |
   | case 50, 'c' stripped dD 0 | unknown, −0.284, 4.7e-7, 236,155 | onto 18,214 rows | infeasible | fails |
   | case 50, 'c' light dD 0, cap_inc 'auto' | unknown, −0.193, 3.1e-6, 242,485 | onto 19,923 rows | infeasible | fails |
   | SeDuMi in both stages: 44min 'Pd' light at epneg 0 and 1e-1, x3min mu 0.4 'Pd' light at epneg 0 and 1e-2 (the rows of pc_diagsed.m) | inaccurate (numerr 1), 0.983–1.004, 5.6e-10 to 2.5e-9 | rows met to 8e-17 | certified, 0.977–1.009, 1.8e-11 to 5.5e-10 | all four pass |

   - The projection check x̂ = x − A'(AA')⁻¹(Ax − b), then λ_min of every Gram block, certifies nothing. λ_min after
     the projection is negative on every row of both stages, the controls included; on the rows that pass stage 2
     it is −9e-12 to −7.6e-9, or −2e-11 to −1.4e-6 relative to the largest eigenvalue of the block. Before the
     projection their smallest Gram eigenvalue is 1e-12 to 9e-10, the size of the projection step: the solutions
     lie on the boundary of the cone. An exact check needs a strictly feasible slack first (facial reduction, or a
     margin the solve can keep).
   - "Passes" means that the slack-only program passes the gate of a harness certificate with the storage held at
     the repaired stage-1 value: as reliable as a solve the gate accepts, not exact. Measured: the projection of
     the whole stage-1 point onto its equality rows moves it by at most 1e-6 of its norm for 36, 41, 69mix2, 44min
     'Qd' and the SeDuMi rows, by 5e-3 for the 44min 'Pd' row at epneg 1e-3, and by 8e-2 at epneg 1e-4. Inferred
     from that: in the first group the measure gate misread certificates; the 44min 'Pd' row at epneg 1e-3 is
     certified by the repaired storage, not by the stage-1 one.
   - New against §1: case 36 at epneg 0 passes at 52,880 stage-1 variables, where the certificate of §1 needed
     380,804; case 41 passes at the default eppos; SeDuMi alone certifies the minimal forms of 44 and X3. Case 50
     (dual H2) stays open: the fixed storage admits no slack at any of the three starts.
## 3. The trade pass on the 80 one-variable rows (TRADE3)

TRADE1 (11/19) stopped at the start of 80 of 115 rows (PHASE_C.md §13). TRADE3, with the statement split fixed
(ex_search apply_mods, worktree copy), converged on all 80 and cut 31 (103,755 -> 62,584 variables). Every test that
shares the example st was then verified (EXVT3, the Phase C harness verification): 22 example settings pass every
test, 87,285 -> 53,357 variables in sum (trade3_verified.txt); 9 do not (the cut of one test breaks its dual or a
second executive: duals of PDE 7, 9, NDS01; gain duals of DDE10, 11, 12, PDE 21; the DDE09 primal; DDE17 stability).
Library path on 41c9210c: the 13 PDE files (16 tests) and the 8 DDE/NDS blocks (11 tests, settings received as
written) pass; EXVT3, the library path and the blocks give identical results on f474e7de (the release comment
pass); Saba 4.4 (per-case st) is verified by the harness only. Largest cuts: PDE 32 29,227 -> 13,971; DDE16
23,449 -> 13,486; DDE13 11,759 -> 5,791; DDE23 3,057 -> 1,812; PDE 52 2,471 -> 1,755. Candidates for the example
files, not written (the release comment pass is editing them).

## 4. What was tried, by case

Jobs of this run (tags of pc_jobs.csv and the pc_diag* scripts); the 11/19 diagnostics (A0, AmT, A − cT;
PHASE_C.md §9) are cited where they decide a reading. Variable counts in parentheses.

| Case | Jobs | Varied | Result |
|---|---|---|---|
| 36 | W36 15, W36F 11, W36V 7, W2S 1 | eppos 1e-3 to 1e-1, epneg 1e-4 to 1e-1, both; deg_Q{2} [1 1 1], [0 1 1], [2 2 2]; psatz{2} 'faces1', 'none'; the cap; stripped, heavy; SDPT3, SeDuMi; the stock 2-D storage in 'P' | certified at epneg 1e-3 (52,880; also stripped, with the cap, with deg_Q{2} [1 1 1]), epneg 1e-1 (380,880), deg_Q{2} [0 1 1] at epneg 0 (380,804); the rest unknown, measure 0.36–9.1; SDPT3 no solution information; the stock storage: MATLAB access violation in pass 2 (751,936), not rerun |
| 41 | W41 42, DIAG41 8, DIAG41B 20, W1V 4 | 'Q' and 'P' at light and stripped, seven margin settings each; deg_Q{1} [2 2 2] to [8 8 8] (to 286,624); len_sos, len_com, dw_feas; heavy, veryheavy; SDPT3, SeDuMi; the input as a flux, decoupled; the explicit delay state with margins; the dual forms 'Pd', 'Qd' at light and stripped | as written unknown at every primal setting (measure −0.89 to 0.53) and infeasible in 'Pd', 'Qd' (4,000 / 7,888); flux and decoupled variants certify; the explicit delay state certifies |
| 44 | W44 17, DIAG44 13, W1V 4, DIAG44D 9, DIAG44E 5, DIAG44F 3, DIAG44G 10, DIAG44H 7 | deg_Q{1} to [10 6 6] (1,251,116; [8 8 8] above the gate); margins; 'P' with len_com, len_sos [4;4;4], heavy, veryheavy; SDPT3; 13 formulations; 'Pd', 'Qd' on the example, minimal and Riemann forms; the minimal form in 'Pd' at epneg 1e-4 to 1e-1, stripped, heavy; in 'Q' at deg_Q{1} [8 6 6] | as written infeasible in every form ('Q', 'P', 'Pd', 'Qd'); minimal form infeasible in 'Q' (to 354,007) and 'P', certified in 'Pd' (2,520) at epneg 0 to 1e-1 (light) and 1e-3 (heavy, 14,373), unknown at epneg 1e-4 (measure 1.55) and at stripped epneg 1e-3 (1.32); in the mu sweep certified to 0.49, not from 0.55 (stable to 2/3 at tau = 1; finding 8); Riemann states certify in every form |
| X3 | WX3 9, WXC 6, DIAGX3 8, W1V 2, DIAG44E 5, DIAG44F 4 | deg_Q{1}, deg_Q{2}; margins; the stock storage (light); heavy 'P'; cross terms in both forms; Riemann states; 'Pd' on the example; the minimal form with the delay on s2 in 'Pd' at mu 0.375 and 0.4, epneg 0 to 1e-2, light, stripped, heavy | as written infeasible (one unknown, measure −0.50), also in 'Pd'; minimal form certified in 'Pd' at epneg 0 and 1e-2 (9,000; measure 0.994–1.002), at epneg 1e-3 inside the measure gate with eta 2–3e-7, at heavy epneg 1e-3 unknown (63,621, numerr 2); Riemann states certify |
| 39 | W39 8, WXC 2, W2S 1, W2X 2, W39S 5, W2SH 1, W39C 1, W39CL 1, LIB_39, LIB_39L | 'Q': deg_Q{2} [1 1 1], [2 1 1], [2 2 1], margins, psatz{2} 'faces1', cross terms; 'P': cross terms, the stock storage at light and heavy, slack cuts | 'Q' unknown everywhere (measure 2.4–7.3, optimal, eta about 1e-8); 'P' with cross terms infeasible; 'P' with the stock storage at light or heavy and dD 0 certified (334,656; 335,601), identical on 41c9210c, library path optimal; heavy at dD 1 certified (1,692,801) |
| 69 | W2Q 5, W2S 1, W2X 2, W2SH 3, W2P 8, W2E 3, W2VS 2, W2VL 1, DIAG69 (1-D twin), DIAG2C 2, DIAG2M 1, DIAG69M 2, DIAG69N 5, DIAG69P 3, DIAG69T 2, LIB_69, SND69 2, DIAG69S 2 | 'Q': deg_Q{2} [2 1 1], [0 1 1], [0 2 2] and 2 (above the gate), psatz{2} 'faces1', 'corners', box quadratics, faces with the box product, cross terms, eppos 1e-2, the cap; 'P': stripped, the stock storage at light and heavy, cross terms; 'Pd' at stripped and light, dD 0 and 1; the 1-D twin in the position state in every form; the interface removed; the mixed position–strain form | as written infeasible with certificate on every solved job (336,681–1,696,905; box quadratics and faces with the box product above the gate, 3.4–3.9M); the 1-D twin certifies in 'Pd' only; without the interface certified (928,900); mixed form certified in 'Q' light or stripped dD 0 (267,060, measure 0.722, the same program under eppos, cap and wcross changes), identical on 41c9210c, library path optimal |
| 70 | W2Q 5, W2S 1, W2X 2, W2SH 3, W2P 8, W2E 3, W2VS 2, W2VL 1, DIAG69 (1-D twin), DIAG70 3, DIAG2C 4, DIAG2M 3, DIAG70N 5, DIAG70T 4, LIB_70, SND70 2, DIAG70S 3 | as 69; the velocity of order (2,2); the edge states removed (Neumann at s2 = 0) or replaced by the damping datum ph_s2 = ps; the mixed form at light, stripped, with the cap, without the Psatz terms; W2D (deg_Q{2} 2, 2.0M) stopped before its solve | as written infeasible with certificate on every solved job (362,596–1,802,596); with the velocity of order (2,2) infeasible; Neumann edge certified (380,880); damping datum infeasible (numerr 2); mixed form certified in 'Q' light or stripped dD 0 (187,116, measure 1.034), unknown at dD 1 (758,364, measure 0.55) and without the Psatz terms (707,680), identical on 41c9210c, library path optimal |
| 56 | W2X 3, W56 4, W2VS 1, W56D 2, DIAG56 2 | 'Q' with cross terms and the cap, epneg 1e-3, eppos 1e-2; 'P' stripped, cross terms; 'Q' stripped with the cap; 'Pd' light dD 0; deg_Q{3} [1 2 2] and 2 (the setting that certifies the 1-D twin, PDE 43); the delay state of order 0 in space | infeasible with certificate on every solved job (266,148–1,587,084); deg_Q{3} [1 2 2] and 2 above the gate (7,857,280 and 7,859,628; m 135,522 and 138,686); 'Q' stripped with the cap above the gate (2,524,300); the delay state of order 0 in space larger (12,781,675 at [1 2 2]) |
| 50 dual | EXS_50d 9 starts, EXS_50d2 3 starts, EXS_50d3 3 starts, W50S 4 | 'cco' searched from stripped to heavy starts; 'c' from stripped starts (dD 0 with the cap, dD 0, as is) and light and heavy dD 0 starts; the stability test of the system in 'Pd', 'Qd' and 'P' | 'cco' infeasible with certificate at every start (60,645 to 1.66M); 'c' unknown at every start within the gate (stripped 66,923–887,267, light dD 0 242,485 and 1,585,440; measure −0.38 to −0.19; gam 305–629); 'c' at light (dD 1) and heavy 5.2–6.0M, above the gate; the system certifies in 'P' (98,080) and fails in 'Pd' (96,416–501,968) |
| 55 | Sol job 2067511; W55 2; W55C 3 | stripped 'P' at dD 1 on Sol (64 threads); 'P' stripped dD 0 with the slack weight in three variables raised (st.sop.dw_feas [0 1 2], [0 1 3]); numeric caps at dD 1 | Sol: cancelled at the 6 h limit inside the solve (7,340,048 variables, m 105,717; 191 GiB; about 327 core-hours); dD 0 infeasible at weight 1, 2, 3 (114,704; 284,944; 684,304); caps [2 2 2], [3 3 3], [4 4 4] give the uncapped program (7,340,048), the 'auto' cap 1,808 (infeasible) |

## 5. Not tried, and why

1. **Programs above the local gate** (2e6 variables or m 45,000): 69 'Q' with deg_Q{2} 2 or [0 2 2] (4.46M, m
   70,475), 70 at deg_Q{2} 2 (2,015,785; W2D stopped before its solve once the mixed form certified), 44 with deg_Q{1}
   [8 8 8] (2.4M), 56 at deg_Q{3} [1 2 2] (7.86M, m 135,522), 55 at dD 1 (7.34M; the Sol job ran out its 6 h
   limit), the 50 dual in 'c' at light (dD 1) and heavy (5.2–6.0M). Sol jobs above 10 CHE wait for the maintainer.
2. **A full energy-state form of 69 and 70 in two variables.** The Riemann construction does not map onto the PIE
   format in two variables (e2 = u_s2 has no datum at s2 = 0 or 1 when u = 0 there). Built instead: the mixed
   position–strain forms, the strain a state only in the coupled direction (Section 1); both certify.
3. **The stock 2-D storage on 56.** exec_tools_sop reads settings_2d only for programs in two variables (nv == 2);
   the delay state of 56 adds a third.
4. **SeDuMi on 2-D cases.** The 39 SeDuMi job raised to 1,828,161 variables (Cholesky factor nnz 2.3e8, about 20 min
   per iteration, 24 GB) and was stopped; SeDuMi left the 2-D grid. SDPT3 returned no solution information on 36,
   41, 44 and X3.
5. **The parser delay expansion of 41** (19,776 variables against 2,496 with an explicit transport state): a
   question for the parser (Declan); not changed.
6. **The gate policy (D93):** not changed; certified means both gates throughout.

## 6. Open cases, candidate statements, defects

**Open – no certificate within the local gate:**

1. **56** (0.9 tau*): infeasible at every local setting tried (to 1,587,084 variables). The degree that certifies its
   1-D twin (43: deg_Q{2} [1 2 2]) gives 7,857,280 variables, m 135,522 at deg_Q{3} [1 2 2] – above the PDE 55
   program that ran out the 6 h Sol limit. The default degree fails at every delay tried (DIAG56T: 43 at 0.5,
   0.75, 0.9 tau*; 56 at 0.25, 0.5, 0.75 tau*): both systems decay at rate 0.1 without the delay against a delayed
   feedback of gain 1, so the degree raise is needed at any delay, not near the threshold only; 43 in 'P' does not
   certify at any setting tried (W43P), so 'P' offers 56 no cheaper route.
2. **55** (3-D coupled heat): only the stripped 'P' program at dD 1 (7,340,048 variables, m 105,717) remains; dD 0
   is infeasible at slack weight 1, 2 and 3, and the caps give either 1,808 variables (infeasible) or the full
   program. The 6 h Sol job (64 threads) reported no iterate; a longer job is the maintainer's decision. Drafts, not
   submitted: sol_drafts/ (55: 48 h, cap 6,144 CHE; 56 at deg_Q{3} [1 2 2]: 48 h, 384 GB, cap 7,680 CHE). The
   maintainer deferred both (11/20 16:34, "Not now").
3. **50 dual (H2)**: 'cco' infeasible with certificate at every start, where the primal 'oco' certifies; 'c' unknown
   (measure −0.38 to −0.19) at stripped and at light dD 0 (to 1,585,440), 5.2–6.0M at light (dD 1) and heavy. A 'cco' certificate W (W >= eppos I,
   (A*W)*T' + T*(W*A') + B1*B1' <= 0) is a 'Pd' stability certificate (P = W), and the stability test of the
   PDE 50 system is infeasible in 'Pd' at every size within the gate (96,416–501,968) while 'P' certifies it at
   98,080 (W50S). Inferred, not measured: the dual H2 test fails for the reason the dual stability test fails in
   two variables (as for PDE 38); the two programs parameterize W and P with their own degrees. The 'c' form
   certifies the dual H2 test of PDE 49, a one-variable system (0.10157588 at veryheavy, 11,847 variables; 11/19).
   The 'c' points that ended unknown do not verify with their storage fixed: the slack program is infeasible at all
   three starts tried (§2 item 11, 11/21).

**Certified only in another formulation (as written, no form tried certifies):** 41 (the file's explicit delay
state), 44 and X3 (the minimal position form in 'Pd'; the Riemann form, example 59), 69 and 70 (the mixed
position–strain forms). Changing the formulation of an example is the maintainer's call.

**Candidate statements for the example files** (not written – the release comment pass is editing them; each
verified on aad9606c, 41c9210c and f474e7de, and through the library path on 41c9210c and f474e7de):

| Case | Form | Statements | Variables | Library path |
|---|---|---|---|---|
| 39 | 'P' | st.settings_2d = settings_2d_stock('light') (or 'heavy'); st.lpi.ineq.dD = 0 | 334,656 (335,601) | optimal |
| 36 | 'Q' | st.lpi.ineq.dD = 0; st.lpi.ineq.wcross = 1; st.epneg = 1e-3 | 52,880 | optimal |
| 41 | 'Q' | expand_delays = true (the file's branch); st.eppos = 1e-2 | 2,496 | optimal |
| 70 | 'Q' | the mixed form (pc_lib70_edit.pl); st.lpi.ineq.dD = 0 | 187,116 | optimal (worktree copy) |
| 69 | 'Q' | the mixed form (pc_lib69_edit.pl); st.lpi.ineq.dD = 0 | 267,060 (measure 0.722) | optimal (worktree copy) |
| 44 | 'Pd' | the minimal form (pc_lib44x3_edit.pl); stability_dual = 1 with stability_form 'P'; st.epneg = 1e-1 | 2,520 | optimal (worktree copy) |
| X3 | 'Pd' | the minimal form, delay on s2 (pc_lib44x3_edit.pl); stability_dual = 1 with stability_form 'P'; st.epneg = 1e-2 | 9,000 | optimal (worktree copy) |

TRADE3 (Section 3): 22 verified cuts, trade3_verified.txt.

**Defects found (reported, not changed):**

1. pie_witness_sop, mode 'rate' (the "numerical spectrum max Re" of PIETOOLS_stability_sop): positive values on
   every 1-D position-state PIE built here, growing with the Legendre degree N on two of them (dw1 0.91, 2.36, 4.99,
   10.12 at N = 12, 24, 48, 96, though dw1 has no eigenvalues) – spurious eigenvalues of the Galerkin pencil (A,T).
   The certificates are unaffected; the printed value misleads.
2. st.lpi.ineq.dw and st.lpi.ineq.wmin, set in every tier, do not act on PIETOOLS_stability_sop (it passes dw;
   base 0 in 'Q'). A row of 11/19 (55, wmin [0 0 1]) and six rows of this run tested nothing for that reason.
3. The forms 'Q' and 'P' do not certify a wave in the position state whose Neumann datum is the trace of another
   state (finding 7), in one and two variables; in one variable 'Pd' does, in two variables nothing tried does
   except the mixed forms.
4. The parser's delay expansion of 41 (19,776 variables, inconclusive) against the explicit transport state (2,496,
   certified): for the parser (Declan).
