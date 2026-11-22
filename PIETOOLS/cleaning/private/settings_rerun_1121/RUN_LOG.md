# Settings rerun – run log (11/21/2026)

Maintainer, 11/21/2026: "Yes to all three, run the plan. But preference is for certificates which directly
certify if the complexity cost is relatively low -- i.e. less than 3 times, not counting cases near the memory
threshold of 32GB RAM." Then: "Only use FQ and FQd when other tests fail. I anticipate that there may be no cases
where these tests succeed and others fail."

Decisions (PLAN.md §6): (1) a stage-2 pass counts as a certificate for the example settings; (2) the instruments
go into the toolbox; (3) dual and finite-energy tests for two- and three-variable examples, finite-energy forms
only where every other test fails. No dual finite-energy form exists in PIETOOLS_stability_sop ('FQd' is not a
form); the finite-energy rows use 'FQ' and 'FP'.

Rule implemented in ex_search (opts.post): a row ends on its smallest direct certificate (harness gate of 11/20)
when that costs at most 3 times the decision variables of its smallest certificate, direct or a posteriori, and
lies below 60% of the size gate (1.2e6 variables, 27,000 rows; the size gate of 2e6 / 45,000 is the workstation
limit of 32 GB); else on the smallest certificate.

- 02:10 Instruments (uncommitted, main checkout): info.aux in PIETOOLS_stability_sop and PIETOOLS_H2_norm_sop;
  the opt-in SOS_KEEP_SDP flag in SOSTOOLS400/sossolve.m (worded in SOSTOOLS terms). Harness
  (testing_codes/private/baselines/2026-11-18_studies/private): post_verify.m (from pc_post.m, extended to 'FP',
  'FQ' and the direct 'Qd'/'FQ'; soundness by the shifted PIE A + delta*T), lb_exec.m (LB_KEEP_LAST), ex_search.m
  (opts.post, the rule, start keyword examplest, 'FQ' in the 'Q' moves), val_post.m.
- 02:36 Validation on worktree C:/wt1121r (905b0d26 + the changes; logs/VALPOST.log): case 36
  (PIETOOLS_PDE_Ex_2D_Telegraph_Eq) 'Q' light dD 0 wcross 1 – stage 1 unknown (measure 0.375, eta 1.8e-8),
  post_verify passes (measure 1.121, eta 1.0e-8, as pc_post); shifted by 10 infeasible. Case 40
  (PIETOOLS_PDE_Ex_2D_Heat_Eq_with_ODE) certified directly – post_verify passes; shifted by 10 infeasible. Case 50
  (PIETOOLS_PDDE_Ex) dual H2 'c' stripped – fails (as pc_post). ex_search opts.post on case 40 from its example st:
  start certified directly, final path direct.
- 02:37 Starts (make_starts.pl): starts_P1.csv (13 keys: 36, 39, 50, 50/H2 dual, 53, 54, 55, 56, 57, 62, 63, 65,
  72; extra alternatives for 36, 39, 53–56), starts_P1d.csv (the dual test of 17 certified two- and three-variable
  stability rows: 'Q>Qd' or 'Q>Pd'; example st, light, heavy), starts_P2.csv (37 one-variable rows). Forms: PDE
  rows the form the example names (results_SWAPSURVEY.csv), DDE/NDS/DDF rows the Phase A form; starts: PDE rows
  the example st, other rows the tier and changes of PHASE_C_TABLE.md. Not searched: 69 and 70 as written (no form
  certifies them; their mixed forms are separate systems, problem_cases_1120 §2 item 7).
- 02:37 Queues started (rr_run.sh): RR_P1 then RR_P1D (one MATLAB process), RR_P2 (the other); sessions of 17
  minutes of builds; logs/RR_RR_P1.log, RR_RR_P2.log; results exs_RR_*.csv and exs_RR_*_final.csv in the worktree harness folder.
- 02:48 Audit done (stability_form_swap_1121/AUDIT.md): five conclusions rest on the 'P' programs of the 11/17 HEAD reference (D91, D93, q_p_differentiation_proposal, EXAMPLE_SETTINGS.md:354, PHASE_C.md:527); correction lines appended there; they bear on 36, 39, 53, all in starts_P1.csv.
- 03:11 Maintainer: "this info is already in the prog structure". The SOS_KEEP_SDP flag is dropped: post_verify rebuilds the SDP from prog (prog.expr.At/b, prog.var, prog.extravar, prog.solinfo.x – the assembly of sossolve's processvars; after simplify sossolve restores the full matrices and maps x back). SOSTOOLS400/sossolve.m restored to HEAD in the main checkout; lb_exec no longer sets the flag. The running queues keep the flagged copy on wt1121r; the new post_verify (compatible: a captured SDP is used when present) copied there. val_sdp.m compares the rebuilt SDP with the captured one on cases 36, 40, 50 when a MATLAB slot frees.
- 03:18 val_sdp (logs/VALSDP.log, results_VALSDP.csv): the SDP rebuilt from prog equals the one sossolve received on
  cases 36, 40 and 50 (K, RR, At identical; b within 4e-17 relative; x identical; simplify off), and post_verify
  gives the same verdicts on both (36 passes, measure 1.121; 40 passes; 50 'c' fails). The flag is not needed.
- 03:11 RR_P2 done (37 one-variable rows; exs_RR_P2_final.csv on the worktree): 21 rows converged on their start,
  direct; 7 smaller – 49 H2 'o' 11,847 -> 1,331 and its dual 'c' 11,847 -> 2,507, DDE25 dual 'c' 865 -> 658, DDE11
  732 -> 596, NDS11 68 -> 57, DDF03 5,669 -> 5,661 (all direct), and 41 (PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_
  Input) 19,776 -> 2,200 through stage 2 only (no direct certificate on its path). DDE09's start fails the gate of
  11/20 (unknown, measure -0.31, eta 1.1e-8; stage 2 does not rescue it). Eight PDE rows and 13.2 errored: their files
  return [PDE_t,PDE_b,st] and example_st took the second output; fixed (the last output), rerun as RR_P2B.
- 03:25 INVALIDATED: RR_P1 (partial), RR_P2, RR_P2B, DBG15. rr_run.sh called ex_search without lb_setup(), so
  each session ran on MATLAB's saved path – the main checkout's toolbox and example files, which hold other
  sessions' uncommitted edits – with only the harness files from the worktree (the current folder). After the
  main-checkout sossolve was restored at ~03:05, lb_exec found no recorded SDP and never kept the solve, so the
  stage-2 path was silently off (RR_P2B: 0 post notes; RR_P1 rows after 03:05). Found by an instrumented copy of
  ex_search (LB_LAST absent after lb_exec). Their results are not used. Fix: rr_run.sh calls lb_setup() and logs
  which() of the executive, sossolve and post_verify; the worktree now holds the final design (sossolve at HEAD,
  lb_exec without the flag, post_verify rebuilding the SDP from prog). Smoke test SMOKE2 on the worktree: case 15
  'Qd' – stage 2 runs and fails (infeasible), case 36 – stage 2 passes (measure 1.12, eta 1e-8). The RR_P1 processes
  (bash loops and MATLAB) were stopped; 13.2 now starts in 'Q'/'Qd', the forms its file sets for index 13.2.
- 03:26 Relaunched as RR2_P1 then RR2_P1D (one MATLAB) and RR2_P2 (37 rows, the other); logs/RR2_*.log.
- 03:40 Maintainer: the rebuild "is fragile to changes in preprocessing. But its not really needed for the
  fixed-storage verification. All we need is the storage operator itself, which is obtained from getsol". RR2
  queues stopped (mine). post_verify rewritten without SDP data: the operators of the solve from
  getsol_lpivar_sop(prog, info.aux.*); -D sought first with the slack of the solve rebuilt exactly at the operator
  level (poscopvar_direct with the degree specification, terms, kept basis operators li.include and margin that
  lpi_ineq_sop records, on the spaces of the solve's Dm), then, if that fails the gate, with a slack sized from the
  support of -D one degree up (inside try: the size gate may refuse it). The margin of P is checked only in the
  direct forms (P = sym(T'Q)); in the storage forms P is the Gram parametrization, positive by construction (a
  slack for it built 7.7M variables on case 62). Coefficients of D below 1e-9 of its largest are set to zero:
  the residue of equality rows the slack of the solve cannot reach ("0 = 1e-12" rows made case 62 inconsistent).
  Forms with a storage equality ('Q', 'Qd', 'FQ', H2 'c') are checked only when the solve has eta <= 1e-7.
  val_sdp.m removed (nothing is rebuilt now).
- 04:15 Validation (val_post 'a', logs/VALPOST6.log): case 36 passes (the slack of the solve did not pass this time,
  the one-degree-up slack did: measure 1.12, eta 5.9e-9); case 62 (PIETOOLS_PDE_Ex_Road_Field_Berestycki) 'P'
  passes with exactly the slack of the solve (82,980; dropped entries <= 9.8e-12 of 0.015); case 40 passes; every
  row shifted by 10 is infeasible; case 50 'c' not checked (eta 1.7e-6).
- 04:17 Relaunched as RR3_P1 then RR3_P1D, and RR3_P2; logs/RR3_*.log.
- 04:45 RR3_P2 done (one-variable rows, on the worktree with lb_setup): 9 rows smaller, every one direct, and the
  3x rule chose the direct certificate each time – 49 H2 11,847 -> 1,331; DDE25 dual 865 -> 658; DDE11 732 -> 596;
  9: 227 -> 155; 9 dual 227 -> 187; 21 gain dual 1,401 -> 1,257; 22 313 -> 301; NDS11 68 -> 57; DDF03 5,669 ->
  5,661. Not certified at their start: DDE09 (fails the gate of 11/20, measure -0.31); 49 dual and 15 (eta 4.9e-7
  and 6.0e-7 > 1e-7: Q-type forms are not checked a posteriori); 41 (stage 2 infeasible). DDE25 also showed a post
  error LB_SIZE: the size gate of the descent refused the stage-2 program. Fix in ex_search try_post: the caps are
  lifted to EXS_NDVCAP/EXS_MCAP during post_verify and restored after.
- 04:48 RR3_P1 stopped (mine) after part of PDE:36, run with the 1e-9 cleaning; its rows are rerun in RR4.
- 04:50-05:01 The cleaning threshold. Zeroing D by magnitude changed verdicts both ways:
  | zeroed below (absolute)    | 36 'Q' (rmax 1.6e-10) | 62 'P' (rmax 3.7e-11) | 41 'Q' (rmax 5e-9) |
  | 1e-9 of the largest        | fallback passes        | passes                 | infeasible          |
  | 1e-7 of the largest        | passes                 | infeasible             | infeasible          |
  | 10 x residual of the solve | infeasible             | unknown, eta 7.9e-3    | infeasible          |
  (logs/VALPOST6-8.log; rmax = largest entry of R = -D + margin*eD - N, N the slack of the solve.) Every row fits
  one mechanism: entries the slack cannot reach are residue no larger than rmax and must go; entries it can reach
  must stay however small, since a small Gram diagonal entry supports off-diagonal entries of order
  sqrt(diag*max) (62: zeroing entries of 3.5e-10 moved eta from 1e-10 to 7.9e-3). Replaced by structural
  cleaning: slack_like zeroes the constant of each row of P - Pop == 0 that carries no slack variable (params.A
  entries whose params.B column is zero; a fixed block entirely), keeps everything else, and the check fails if a
  zeroed value exceeds rmax + 1e-12 max|D| (on those rows the residue equals R exactly). H2 'c' returns no slack
  operator, so any nonzero residue fails it. Validation VALPOST9 (case 41 smoke + val_post 'a').
- 05:10 VALPOST9 (structural cleaning): 36 passes with the slack of the solve (fr 1.05, eta 7.2e-9, residue 0);
  62 passes with the slack of the solve (residue 8.8e-13 <= rmax 3.7e-11); 40 passes; every row shifted by 10
  fails; 50 'c' not checked (eta 1.7e-6). 41 still infeasible with residue 0. (The slack counts in the notes
  differ from aux.li.ndec by about half: lpi_ineq_sop counts sum(N.^2) per Gram block, the rebuild counted the
  unique names, N(N+1)/2 – the same slack.)
- 05:14 Case 41 (PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input, 'Q' light, example st; stage 1 unknown, fr
  0.078, eta 5.1e-10) split by target (scratchpad dbg41): -D alone passes with exactly the slack of the solve
  (certified, fr 1.005, eta 8.5e-9); the margin of P alone, P - eppos2 T'T >= 0, is infeasible. post_verify sized
  that slack from the numeric P with info.dw (no dwbase, a numeric cap and pruning), not as the executive did.
  Fix: PIETOOLS_stability_sop returns the slack of the margin of P and its info (info.aux.Npm, lip; direct forms
  only, two fields of the uncommitted aux instrument); post_verify rebuilds it exactly with slack_like, residue
  rows zeroed under the same guard. Validation VALPOST10.
- 05:20 VALPOST10 (exact rebuild of both slacks): 36 passes (residue [0 0]), 62 passes (residue 8.8e-13), 40
  passes; every row shifted by 10 fails; 50 'c' not checked. Case 41 (DBG41c): -D passes, the margin of P fails
  with exactly the slack of the solve (3,256 names; residue [0 0]); the target matches the solve's (P-margin
  residual of the solve 7e-10, eppos2 1e-4). Verdict: the start of 41 (example st, 'Q' light; stage 1 unknown, fr
  0.078) is not certified – its storage P = sym(T'Q) does not admit P >= 1e-4 T'T with the slack of the solve.
  (The projection check of the morning passed 41 from another stage-1 point, measure 0.34, after repairing the
  storage: a different verification.)
- 05:20 Relaunched on C:/wt1121r (which() of the executive, sossolve and post_verify: the worktree): RR4_P1 then
  RR4_P1D in one MATLAB (rr_run.sh, at most 20 and 10 sessions of 17 minutes), RR4_P2 in the other (6 sessions);
  logs/RR4_*.log. The interactive MATLAB desktop open since 11/18 (DEMO3) is the maintainer's and is left alone.
- 05:46 RR4_P2 done (37 one-variable rows, 26 min). Smaller, all direct: 49 (PIETOOLS_PDE_Ex_Transport_Eq_with_
  Disturbance2) H2 'o' 11,847 -> 1,331; 9 (PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_in_BC) 'P' 227 -> 155 and its dual 'Pd'
  227 -> 187; 22 (PIETOOLS_PDE_Ex_Heat_Eq_with_Distributed_Disturbance) 313 -> 301; DDE11 gain 'P' 732 -> 596;
  DDE25 H2 dual 'c' 865 -> 820; DDF03 5,669 -> 5,661; 21 (PIETOOLS_PDE_Ex_Wave_Eq_Tip_Damped) gain dual 'Pd' 1,401 ->
  1,257. 29 rows keep their start. Against the superseded RR3_P2: NDS11 now 68 (RR3 57) and DDE25 dual 820 (RR3
  658) – RR3 reached them through intermediates kept by the 1e-9 cleaning (NDS11: 57 after a capinc step). No
  start certified: DDE09 gain 'Q' (fr -0.31), 49 dual 'c' (eta 4.9e-7), 15 (PIETOOLS_PDE_Ex_Timoshenko_Beam_2)
  'Qd' (eta 6.0e-7), 41 'Q' (the margin of P fails). 17 stage-2 passes, 27 stage-2 failures, 0 post errors; the
  largest residue zeroed is 1e-12.
- 05:48 Guard ceiling: zeroed residue must also stay below 1e-7 of the largest coefficient of the target (the
  gate tolerance). An NDS11 intermediate (61 variables) passed stage 2 after a solve with max|R| 1.4e-4 of max|D|;
  nothing was zeroed there (its attempt B passed), but the guard alone would have admitted residue of that size.
  In effect for RR4_P1 calls after 05:48:39 (MATLAB reloads the changed file); earlier stage-2 notes print the
  residue and the largest D, so they can be audited.
- 05:50 catalog_keys.csv (key, library, index, example file, executive) from catalog.mat; report_table.pl names
  every row by example file and library index from it.
- 05:52 RR4_P3 (alternatives for the four uncertified rows; starts_P3.csv): the cheaper and stronger starts
  first, 'FQ'/'FP' last for the two stability rows (maintainer: "Only use FQ and FQd when other tests fail"; no
  'FQd' form exists). The value rows keep their recorded best (gain 25.907746; H2 dual 0.10157588).
- 05:57 RR4_P3 done (5 min). Every row certifies with an earlier alternative; 'FQ'/'FP' never ran:
  - DDE09 (examples_DDE_library_PIETOOLS, Hinf_gain 'Q'): the light tier without Phase A's reductions certifies
    (1,040 variables, gam 25.912755 against the recorded 25.907746); descent to 948 (cap_inc 'auto'), direct.
  - 49 dual (PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance2, H2_norm_c): 'c' at light passes stage 2 but bounds
    0.10586, 4% above the recorded 0.10157588, so the value test refuses it; 'cco' with the example st certifies
    0.10157511 at 7,397 variables; descent to 245, direct.
  - 15 (PIETOOLS_PDE_Ex_Timoshenko_Beam_2, PIE2PDEstability_dual): the primal 'Q' with the example st certifies
    (272,977 variables, eta 6.0e-9); descent to 3,229, direct. The dual 'Qd' (the example's test) stays
    uncertified under the gate of 11/20 (eta 6.0e-7 at the example st; 11/18-19: infeasible at light, heavy,
    extreme). Moving the example to the primal test changes its test: maintainer's decision.
  - 41 (PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input): 'Q' fails at light (the margin of P); 'P' with the
    example st passes stage 2 (stage 1 unknown, fr -0.097, eta 6e-15; 4,000 variables); descent to 324 through
    stage 2 (no direct certificate on its path).
- 05:58 post_shift RR4_P3 (logs/SHIFT_P3.log, results_SHIFT_RR4_P3.csv): 41 'P' at 324 – delta 0 passes (fr 0.84,
  eta 3.4e-9), delta 10 infeasible. Note: in the storage forms the check takes P >= eppos I from the Gram
  construction, i.e. from the returned Gram being in the cone, which an interior-point iterate is; it does not
  re-verify it.
- 05:59-06:14 Example patches for the one-variable rows (worktree only; patches/examples_1d.patch, libcheck
  libcheck_rr_A.csv and libcheck_rr_B.csv, logs/LIBCHECK_*.log). One st serves every test of a file, so each
  candidate went through libcheck, which runs all flagged tests:
  - 9 PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_in_BC: len_com{1} = 1 (the dual's final) – 'P' and 'Pd' optimal at 187 (227).
  - 21 PIETOOLS_PDE_Ex_Wave_Eq_Tip_Damped: len_sos{1} = [2;1;2] – gain 0.8419069, dual 0.84222715 at 1,257 (0.84190244,
    0.84208439 at 1,401).
  - 22 PIETOOLS_PDE_Ex_Heat_Eq_with_Distributed_Disturbance: len_inc{1} = 0 and the file's dD = 0 removed – 0.20843359
    at 301 (0.2084336 at 313). The first candidate kept dD = 0 and gave 0.50000002 at 105: the search had also
    raised dD to the stripped value 1, which a diff from the tier omits. st_patch.m now prints the diff from the
    start (the file's st) as well.
  - 49 PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance2: no common st in 'o'/'c' ('c' fails or eta 4.1e-7 at every st
    tried). H2_norm_form 'oco' (dual 'cco') with the 'o' final st and len_sos{1} = [0;2;2]: 0.10157878 and 0.10158193
    at 1,001 (11,847, dual eta 4.9e-7); a scan of len_sos{1} (scratchpad dbg49b) found no smaller common st; [1;2;2]
    at 1,022 matches the recorded values to 5e-6.
  - Not patched, maintainer's decision: 15 PIETOOLS_PDE_Ex_Timoshenko_Beam_2 – the primal 'Q' with deg_Q{1} =
    [0 1 1] certifies at 3,229 (libcheck optimal), the example's dual 'Qd' does not pass the gate; 41
    PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input – 'P' at 324 certifies only a posteriori, so the library path
    prints unknown (libcheck: unknown).
- 06:14 Form mismatch in the DDE starts: make_starts.pl took the Phase A form for DDE rows, but the library blocks
  of DDE09 (gain Example 1, answer 25.907 at tau 1.5) and DDE10 (Example A.1, answer .2364) run Hinf_gain_form 'P'
  (Phase A chose 'Q', noting 'P' reaches the bound at 1,292 and 4,789 variables). RR4 tested 'Q'/'Qd' for both,
  so the library's own tests were not run. The other DDE/NDS rows match their blocks: DDE11 'P'/'Pd', NDS11
  'P'/'Pd', DDE25/26 H2; the controller and estimator rows (DDE13-24, DDF03) take no form. RR4_P5
  (starts_P5.csv): DDE09 and DDE10, primal and dual, in 'P'/'Pd' with the blocks' settings.
- 06:16 RR4_P5 (the blocks' 'P'/'Pd'): every block st certifies under the gate of 11/20 – DDE09 gain 651 (gam
  25.922824), dual 651 (25.916902); DDE10 1,669 (0.23644325), dual 1,669 (0.23653875) – and 'P' is cheaper than
  the 'Q' RR4 tested: descents to DDE09 623 / dual 596, DDE10 869 / dual 1,285 (all direct). The 'Q' rows of DDE09
  and DDE10 in RR4_P2/P3 are superseded.
- 06:18 Starts against the current library blocks (block_st.pl locates each block by its system lines, since
  lb_meta's line ranges predate the 11/19-11/20 edits): all DDE/NDS starts match their blocks except DDE09 and
  DDE10 (form) and DDF03 (light tier against the block's heavy, len_sos{1} [4;2;2], cap_inc 'auto', var.deg{1}
  [1 2 2]). RR4_P6: DDF03 from the block st certifies at 4,673 (gam 0.42054329, the record 0.420543); descent to
  4,669 – no patch. The light-tier DDF03 of RR4_P2 is superseded.
- 06:20 One st per block for primal and dual (scratchpad dbgdde): DDE09 – the primal final (623) passes both
  (25.92842, 25.925296): patch; DDE10 – neither final passes both: block kept (1,669 certifies both); DDE11 and
  DDE25 – only the block st passes both: no patch. patches/dde_library.patch (DDE09 block, worktree; the main
  checkout's DDE library carries another session's uncommitted edits, so the patch needs rebasing there).
- 06:22 Optional patches, maintainer's decision: patches/optional_15_primal_Q.patch (Timoshenko_Beam_2 runs the
  primal 'Q', deg_Q{1} [0 1 1]) and patches/optional_41_form_P.patch (Heat_Eq_w_Delayed_Boundary_Input in 'P',
  certified only a posteriori). Not applied in the worktree.
- 06:27 RR4_P1 stopped by a planned Windows update restart (System log 1074 at 06:29 and 06:31: MoUsoCoreWorker,
  TrustedInstaller, "Operating System: Upgrade (Planned)"), during PDE:54 (PIETOOLS_PDE_Ex_3D_Advection_Diffusion_Eq):
  its starts 'P' light (7,004,225 variables) and 'P' heavy (29,223,041) were built and refused by the size gate;
  the third, the example st, was building. Rows done before: 36, 39, 50 and its dual, 53.
- 16:52 Resumed (maintainer: "windows rebooted, resume processes"): RR4_P1 then RR4_P1D in one MATLAB. The two
  measured-oversize starts of PDE:54 were dropped from starts_P1.csv (the copy before is starts_P1.pre1650.csv),
  since a resumed row retries every start from the first. post_shift on PDE:36 (accepted through stage 2) in the
  other process; post_shift.m takes an optional key list and appends to its results file.
- 16:53 post_shift PDE:36 (PIETOOLS_PDE_Ex_2D_Telegraph_Eq, 'Q' at 42,308; logs/SHIFT_P1a.log, results_SHIFT_RR4_P1.csv):
  stage 1 refused by the gate (measure 4.87, executive optimal, eta 6.6e-8); delta 0 passes with both slacks of the
  solve (fr 1.01, eta 8.1e-9, residue [0 0]); delta 10 infeasible.
- 17:36-17:54 PDE:54 done (PIETOOLS_PDE_Ex_3D_Advection_Diffusion_Eq): the example st certifies directly at
  1,835,012 variables (fr 0.957, eta 4.6e-10; 1,134 s, near the memory limit); descent to 1,835,009 – unchanged in
  effect. PDE:55 (PIETOOLS_PDE_Ex_3D_Coupled_Heat_Eq): 'P' light built 28,016,900 variables in 13.6 min and the size
  gate refused it; the next session started the same build again. A row with no certified start retried every
  start from the first in each session, so a start that fills a session repeats forever. My RR4_P1 processes were
  stopped (bash chain, rr_run loop, MATLAB at 14.5 GB; the maintainer's MATLAB desktop left alone). Fix in
  ex_search: the starts tried are kept in st.altdone and skipped on resume. The two oversize starts of PDE:55 were
  dropped from starts_P1.csv as for PDE:54. Resumed at 17:54: RR4_P1 then RR4_P1D.
- 17:55 PDE:55 (PIETOOLS_PDE_Ex_3D_Coupled_Heat_Eq): the example st builds 7,340,048 variables, 105,717 rows, above
  the size gate (2e6, 45,000); no start certified on the workstation – open, as on 11/20.
- 17:59 PDE:56 (PIETOOLS_PDE_Ex_2D_Delayed_Heat_Eq_Fridman): 'Q' light dD 0 is infeasible with a certificate
  (840,156 variables; no stage 2, which only takes numerical refusals); the example st builds 2,524,300 variables,
  51,558 rows, above the size gate. Open, as on 11/20.
- PDE:57 (PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Fridman): 'Q' 23,121 -> 19,840, direct (dD 0, wcross 1,
  cap_inc 'auto', deg_Q{2} [1 1 1]); patch candidate.
- PDE:62 (PIETOOLS_PDE_Ex_Road_Field_Berestycki): 'P' 83,269, unchanged, direct.
- PDE:63 (PIETOOLS_PDE_Ex_Heat_Edge_Cascade_Variable_Coeff): 'Q' 149,054, unchanged, direct (the example st).
- PDE:65 (PIETOOLS_PDE_Ex_2D_Linear_Turing_Pair): 'P' 1,635,656 -> 250,280, direct (dw_feas [0 0 1] among the
  changes; the 2-D storage struct differs from the light tier, so the patch is read from the start); patch candidate.
- 19:15 RR4_P1 done. Also: PDE:62 (PIETOOLS_PDE_Ex_Road_Field_Berestycki) 'P' 83,269 unchanged (three stage-2 passes
  on its path, the direct start kept by the 3x rule); PDE:63 (PIETOOLS_PDE_Ex_Heat_Edge_Cascade_Variable_Coeff) 'Q'
  149,054 unchanged; PDE:65 (PIETOOLS_PDE_Ex_2D_Linear_Turing_Pair) 'P' 1,635,656 -> 250,280 direct; PDE:72
  (PIETOOLS_PDE_Ex_Transport_Edge_Diffusion) 'Q' 256,948 unchanged. Stage 2 in RR4_P1: 7 passes, 1 failure, 0 errors.
  RR4_P1D started (dual rows; every two- and three-variable example has its dual test commented out, so these
  rows inform but patch nothing).
- 19:20 Patches of the two-variable examples (worktree), each the diff from the st its file returns:
  57 (st.sop.deg_Q{2} [1 1 1], wcross 1, cap_inc 'auto'), 65 (dw_feas [0 0 1]), 39 (form 'P', settings_2d_stock
  ('light'), dD 0 – a form change: 'Q' certifies at no setting tried). 53: its file returns the stripped tier, which
  RR4 never solved (its 'P' light start certified first); measuring the file st before deciding (scratchpad dbg53).
  report_table.pl now shows the start that certified (from the per-build log).
- 19:21 PDE:53 (PIETOOLS_PDE_Ex_3D_Reaction_Diffusion_Eq): the st its file returns (stripped tier, 'P') certifies
  directly at 580,612 variables (fr 0.986, eta 9.7e-10; 177 s), fewer than the 584,704 the search reached from the
  light tier – the file keeps its st, no patch. libcheck of the patches of 57, 65, 39 (libcheck_rr_C.csv).
- 19:24 libcheck_rr_C.csv: 57 optimal at 19,840 (23,121 on 11/19); 65 optimal at 250,280 in 80 s (1,635,656 in 467 s);
  39 optimal at 334,656 in 'P'. Exported patches/examples_2d.patch (57, 65) and patches/form_39_P.patch (39: the
  form 'Q' -> 'P', for approval).
- 19:44 RR4_P1D so far (duals; informative): 33 89 unchanged; 34 72,081 unchanged; 35 no start certified; 37 11,040
  unchanged; 42 69,915 -> 43,065 direct; 43 'Qd' 1,008,600 -> 196,425 through stage 2 – post_shift (logs/SHIFT_P1D_43.log):
  delta 0 passes (fr 1.19, eta 3.5e-8; stage 1 inaccurate, eta 1.6e-8), delta 10 infeasible.
- 19:59 RR4_P1D: 57 dual no start certified; 62 dual (PIETOOLS_PDE_Ex_Road_Field_Berestycki, 'Pd') 427,249 -> 101,428
  through stage 2 – post_shift (logs/SHIFT_P1D_62.log): delta 0 passes (fr 0.87, eta 7.8e-10, residue 7.4e-11 <=
  the solve's 4.7e-8; stage 1 unknown, fr 0.06), delta 10 infeasible (the slack of the solve; the one-degree-up
  fallback is above the size gate).
- 20:00 RR4_P1D: 63 dual (PIETOOLS_PDE_Ex_Heat_Edge_Cascade_Variable_Coeff, 'Qd') 652,407 -> 32,268 through stage 2 –
  post_shift (logs/SHIFT_P1D_63.log): delta 0 passes (fr 1.03, eta 2.2e-8), delta 10 infeasible. Two notes: stage 1
  reads 'infeasible' on the harness measure (fr -0.9999) with no infeasibility flag (pinf = dinf = 0) and the
  executive 'unknown', which the plan counts as a numerical refusal; and the storage is badly scaled (max|D| 4.6e5,
  residual of the solve 0.042), so eta 2.2e-8 is an absolute residual near 1e-2 – the limit of every certificate
  accepted on eta, direct or not, largest here.
- 20:18 RR4_P1D done (17 dual rows; table in REPORT.md §4): 64 288,400 -> 219,044 and 42 69,915 -> 43,065 direct;
  stage 2 for 38 (380,804 -> 42,308), 43, 58 (586,772 -> 150,400), 62, 63; no start certified for 35, 57, 61; the
  rest unchanged. Stage 2: 14 passes, 3 failures, 0 errors. post_shift on 38 and 58 running.
- 20:20 post_shift on 38 and 58 duals (logs/SHIFT_P1D_38_58.log): both pass at 0, fail at 10. All seven rows
  accepted through stage 2 (36, 41; duals 38, 43, 58, 62, 63) pass the shift check. 38: max|D| 8e-9 (the wave
  equation is conservative; epneg 0). 58: stage 1 'infeasible' on the measure (-0.97), executive 'unknown'.
  REPORT.md completed.
- 11/22 00:20-00:26 Applied to the main checkout (maintainer, 11/22: "Yes, apply them"; not committed):
  examples_1d (9, 21, 22, 49), examples_2d (57, 65), dde_library (DDE09 block). HEAD had moved to ff7a0ee1 (comment
  pass 2 rewrote the header notes and markers of 9, 21, 22, 49 and the DDE library), so the edits were re-applied
  to the current versions (scratchpad apply_rr.pl): header entry after the "Changes since PIETOOLS 2025" note
  (after the Initial coding line in 57, 65, created after the release), markers 11/22/2026 at column 77, the
  continuation indent of the new notes; st and form lines identical to the verified worktree files (checked line
  by line). Verified in the main checkout: libcheck_applied.csv – 9 'P'/'Pd' optimal 187; 21 0.8419069/0.84222715
  at 1,257; 22 0.20843359 at 301; 49 0.10157878/0.10158193 at 1,001; 57 optimal 19,840; 65 optimal 250,280; the DDE09
  block through the harness (logs/CHK_DDE09_APPLIED.log): 25.92842 / 25.925296 at 623, both pass. Not applied:
  form_39_P, optional_15_primal_Q, optional_41_form_P.
- 11/22 00:27-00:34 'Q' contains 'P' (maintainer: "this shouldn't happen, since P is a special case of Q"),
  q_contains_p.m (this folder; logs/QCP_41*.log, q_contains_p_41*.csv). Case 41, P from its final 'P' solve (324
  variables), Q = P*T against the 'Q' program at the light st its file runs:
  (a) P*T lies in the span of the decision Q (deg_Q [4 2 2]): no coefficient outside (max |P*T| 87.7).
  (b) T'PT - eppos2 T'T >= 0 on the margin slack of the 'Q' solve (3,256): executive optimal, eta 5.7e-9, measure
      2.04 – refused; the same at eppos2/10 and 0. Control: T'T alone certifies on that slack (fr 1.03).
  (c) -D on the -D slack of the 'Q' solve (6,723, sized from the support of D): optimal, eta 2.5e-8, measure 1.55 –
      refused; on the slack of the 'P' solve (172, sized like P): certified (fr 0.85, eta 3.4e-9); on slacks sized
      from the targets' own supports (288 each): infeasible.
  Room: adding 1e-4 of the target's scale times T'T certifies both (b) (fr 1.08) and (c) (fr 1.13), and 1e-2, 1 too.
  So the 'P' point is in the 'Q' program, on the boundary of the cones of its slacks: the solver reaches it to 1e-8
  without a measure the gate accepts. The slack of -D sized like the generalized storage sym(T'Q) contains it
  with an interior at raise 0 (certified, fr 0.97, 1,280 variables) and 2 (fr 1.16), not at 1 (fr 1.58).
- 11/22 00:28-00:38 Case 39 (logs/QCP_39.log, SPAN_39.log): P from its 'P' certificate (334,656), the 'Q' program at
  light (648,900: optimal, eta 2.4e-9, measure 7.14). (a) P*T is outside the span of Q at deg_Q{2} = 1 (5.65e-7 of
  max 1.79), 2, [2 1 1], [2 2 1], [2 2 2] (1.1e-7); inside from 3 (Q 11,664 variables at 3, 900 at 1). (b) the
  margin and (c) -D on the slacks of the 'Q' solve both certify (fr 0.99, 1.13). So in 39 the degree of Q excludes
  the 'P' point; the search of 11/19-11/20 stopped at deg_Q{2} [2 2 1]. Running 'Q' light with deg_Q{2} = 3 (and dD 0)
  (logs/Q39_DEG3.log). Case 41 more (logs/QCP_41v5.log): -D on a slack sized like sym(T'Q) certifies at raise 0 and
  2, not 1; the margin on such slacks at raises 0-2 stays degenerate (measure 1.69-2.04); its room from eppos2 is
  about 1e-6 of the target scale against the 1e-4 that certifies, and scaling P does not change the ratio.
- 11/22 00:41 Case 39, 'Q' light with deg_Q{2} = 3 (logs/Q39_DEG3.log): 31,853,664 variables, 237,369 rows (dD 0:
  21,945,744, 164,247) – above the size gate; the slacks are sized from the supports of the degree-3 decision
  operators. The slacks of the degree-1 solve (117,072 and 207,936 variables) certify the 'P' point, so a 'Q'
  program with Q of degree 3 and slacks sized as at degree 1 would contain it at about the size of the light 'Q'
  program; no executive option sizes the slacks apart from deg_Q. REPORT.md §4a updated.
- 11/22 00:45 Committed the seven applied files as 19103941 (maintainer: "Yes, commit them"); not pushed. 905b0d26 is
  already on origin/ndopvar (it went out with 349b247f); REPORT.md corrected where it said otherwise.
- 11/22 01:24-01:36 Committed 63579184 (functional-inequality library: PIETOOLS_FI.m, DEMO3/DEMO14 help, README,
  testing_codes/private/PIETOOLS_examples/test_examples_FI_library.m). Merged origin/ndopvar (Declan b1c1625c,
  41a8787c, 2d9bc49f; clean merge tree fdf00cd4, tested from its export: FI test 17/17, libcheck 9, 21, 22, 49,
  57, 65 at the verified values, 3 and 7 optimal, DEMO3 and DEMO14 run) as 1c41b872 by the shared-checkout
  procedure; the one dirty file it changes, PIETOOLS_H2_norm_sop.m (info.aux, uncommitted), was merged 3-way
  (both 11/21 log entries kept) and written back; status set and other sessions' diffs byte-equal. Pushed
  2d9bc49f..1c41b872 (maintainer: "Yes, all of it").
- 11/22 01:40 results/ holds what lived only in the worktree C:/wt1121r: the 42 result files of the rerun's runs
  (exs_*, results_SHIFT_*, results_VALPOST_a, results_VALSDP, valpost_starts) and the row states of the RR4 runs
  (exs_state_rr4p1, p1d, p2, p3, p5, p6). The worktree's tracked edits are in the main checkout (its st_patch.m is
  the older copy).
