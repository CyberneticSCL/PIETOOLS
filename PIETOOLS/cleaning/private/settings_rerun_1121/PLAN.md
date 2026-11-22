# Rerun of the per-example settings search – plan (11/21/2026)

Maintainer, 11/21/2026: "I suppose we now have to rerun the customized settings search for the 2D-3D problems";
"then create a plan for rerunning the optimized settings search on the impacted examples. You can also integrate
the a posteriori LPI verification technique from the Fixed-storage subagent."

Cases are named by library index and example file (PIETOOLS_examples/examples_PDE_library_PIETOOLS.m and the DDE
library). Nothing in this plan has been run.

## 1. What the form fix changes, and what it does not

Two substitutions existed in PIETOOLS_stability_sop for PIEs in two or more variables (a delay adds a variable):

| Substitution | Active | Ended |
|---|---|---|
| 'Q' -> 'P' | 10/28 04:12 (2e3fa94e) | 11/17 06:12 in the main working tree (direct 'Q', DECISIONS D68; committed bb367fe1, 11/18 03:11) |
| 'Qd' -> 'Pd', 'FQ' -> 'FP'; gain 'Qd' -> 'Pd' | 10/28 (gain: 10/30 in N-D) | 11/21 02:05 (905b0d26) |

- 905b0d26 changes no program of an admitted form. Of the 140 catalog tests, none requests 'Qd', 'FQ' or gain
  'Qd' on a PIE in two or more variables, a dual form on a PIE with K, or gain 'Q' with Tw nonzero
  (stability_form_swap_1121/results_SWAPSURVEY.csv). No setting of 11/19 builds a different program after the fix.
- The 11/19 search did not run on substituted forms: Phase A took its two- and three-variable choices from the
  11/18 studies, Phase B ran on 11/19, and the library-path checks of Phase C on aad9606c, 41c9210c and f474e7de,
  all after 11/17 06:12. The records and decisions that do rest on substituted results are in
  stability_form_swap_1121/AUDIT.md.

So the rerun is not forced by the fix. It is worth doing for three reasons: (1) decisions that rested on
substituted results (AUDIT.md, class (b)); (2) the a posteriori check, which certifies candidates the 11/19 gate
refused for numerical reasons; (3) 'Qd', 'FQ' and gain 'Qd' now run in two or more variables.

## 2. The a posteriori check (cleaning/private/problem_cases_1120/REPORT.md §2 item 11; pc_post.m)

Stage 1 is the candidate as the search solves it. Stage 2 fixes the storage at its stage-1 value and solves for
the slacks only, with the stage-1 degree specification and terms:
- 'P', 'Pd': the Gram blocks of P clipped to the PSD cone (P >= eppos I exactly), then the slack of -D;
- 'Q': Q projected onto the rows that involve Q alone (T'*Q = Q'*T, which the derivative formula needs exactly),
  then the slacks of P - (eppos2/2) T'*T and of -D;
- 'Qd': P and Q onto their rows (T*Q = P);
- H2 'c': R, Q, W and gam fixed.
A candidate passes when the slack-only program passes the harness gate (measure within 0.3 of 1, numerr 0,
executive optimal, eta <= 1e-7). Measured 11/21: 7 of 15 gate-refused rows pass, among them case 36
(PIETOOLS_PDE_Ex_2D_Telegraph_Eq) at epneg 0 with 52,880 stage-1 variables where its certified setting needed
380,804, and case 41 (PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input) at the default eppos. Controls: two
certified rows pass; a storage tested against an unstable parameter fails. The projection check certified nothing
(the solutions lie on the boundary of the cone) and is reported only, not used to accept.

"Passes" means as reliable as a solve the gate accepts, not exact.

## 3. Rows to rerun

Counts from the 11/19 search logs (testing_codes/private/baselines/2026-11-18_studies/private/exs_*.csv): 979
builds on two- and three-variable rows (2.7 solver-hours) with 16 refusals by the measure gate alone (executive
optimal); 3,386 builds on one-variable rows (1.0 solver-hour) with 218 such refusals.

Priority 1 – two and three variables:
- 1a. Rows with no certified setting or with disagreeing gates in Phase A: 36 PIETOOLS_PDE_Ex_2D_Telegraph_Eq,
  39 PIETOOLS_PDE_Ex_2D_NS_Antonelli, 63 PIETOOLS_PDE_Ex_Heat_Edge_Cascade_Variable_Coeff, 72
  PIETOOLS_PDE_Ex_Transport_Edge_Diffusion (gates disagree); 50 PIETOOLS_PDDE_Ex (H2 and its dual), 53
  PIETOOLS_PDE_Ex_3D_Reaction_Diffusion_Eq, 54 PIETOOLS_PDE_Ex_3D_Advection_Diffusion_Eq, 55
  PIETOOLS_PDE_Ex_3D_Coupled_Heat_Eq, 56 PIETOOLS_PDE_Ex_2D_Delayed_Heat_Eq_Fridman, 65
  PIETOOLS_PDE_Ex_2D_Linear_Turing_Pair, 69 PIETOOLS_PDE_Ex_Wave_Heat_Transmission_2D_Batty, 70
  PIETOOLS_PDE_Ex_2D_Wave_Acoustic_Boundary_Beale (nothing certified). Since 11/20: 36, 39 certified with
  candidate st (problem_cases_1120 §6), 69 and 70 in mixed position–strain forms.
- 1b. Rows whose descent stopped at a measure-gate refusal: 62 PIETOOLS_PDE_Ex_Road_Field_Berestycki (5), 57
  PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Fridman (4), 39 (4), 36 (2), 54 (1).
- 1c. The decisions of AUDIT.md class (b) that bear on a setting (stability_form_swap_1121/AUDIT.md, 11/21):
  D91, D93 and EXAMPLE_SETTINGS.md line 354 read the 'P' programs of the 11/17 HEAD reference as the old 'Q'. They
  bear on 39 (no st set because "the storage form of HEAD cannot be selected"; it was 'P'), 53 ('P' light
  certified on HEAD) and 36; all three are in starts_P1.csv, 39 and 53 with 'P' alternatives.
- 1d. Optional: 'Qd' or 'FQ' tests in two or more variables, only for examples that should carry a dual or
  finite-energy test (none does now: every 2-D and 3-D example has its dual test commented out).

Priority 2 – one variable:
- 2a. Rows whose descent stopped at measure-gate refusals: DDE18 (30), PDE 31 (26), DDE20 (19), PDE 32 (11),
  PDE 48 (10), DDE14 (8), DDE24 (7), DDE11/gain dual (7), PDE 49 (6), NDS11 (6), and 30 more rows with 1–6 each.
- 2b. The seven rows Phase A took from the 11/17 tree (DDE13, DDE18, PDE 13.2/dual, 31, 32, 49, 49/dual), which
  Phase A marked for a rerun on the current tree; the Phase B and C records do not say they were rerun.

## 4. Method

- Search: ex_search (the Phase B tool) from each row's current example st (the 11/19 final settings), with the
  same knobs, moves and budget (40 builds or 20 minutes per row), and these changes:
  1. The gate of 11/20 (harness 'certified', executive optimal, eta <= 1e-7) in place of 11/19's (no eta).
  2. A candidate refused for numerical reasons only – executive optimal, inaccurate or unknown; no infeasibility
     certificate (pinf = dinf = 0); measure outside 0.7–1.3 or eta above 1e-7 – goes to the a posteriori check,
     and is accepted when stage 2 passes the gate. Both stages are logged.
  3. Each row accepted through stage 2 is rechecked with its storage against an unstable parameter where the
     catalog has a setter; the slack program must be infeasible. Rows without a setter are flagged.
- Forms: each example keeps the form it names (Phase C); 1d adds 'Qd' or 'FQ' rows only on request.
- Tooling (each item needs approval; none is in place):
  - the internals of the executive's last pass returned as info.aux (one line in PIETOOLS_stability_sop and in
    PIETOOLS_H2_norm_sop; the 11/21 instrument);
  - the SDP data of a solve: the harness copy of sossolve used by pc_post (post/sossolve.m, debug only), or a
    four-line option in sossolve, worded in SOSTOOLS terms;
  - post_verify.m in the harness, from pc_post stage 2 (stability forms and H2 'c'; gain and estimator forms keep
    the plain gate);
  - the ex_search change above.
- Verification: every changed st through the library path (libcheck.m); the example files themselves are edited
  only by patch (the release comment pass edits them).

## 5. Order and resources

1. post_verify built and checked on the 11/21 controls and on the 7 rows that passed on 11/21 (reproduce them).
2. Priority 1a and 1b, then 1c.
3. Priority 2.
4. Library path and soundness rows; report: per row the old and new st, variables, verdict, and the path that
   accepted it (gate or stage 2); example-file patches.

Workstation only, at most two MATLAB processes, 20-minute sessions while the maintainer is online, a clean
worktree off the current HEAD (905b0d26 or later). Rows above the local gate (2e6 variables or m 45,000; 55 and 56
at heavy) only on Sol with the maintainer's word. Estimate: Priority 1 about 3 to 5 hours (the 11/19 N-D search
took 2.7 solver-hours; stage 2 costs about one stage-1 solve per refused candidate), Priority 2 about 1.5 hours.

## 6. Decisions for the maintainer

1. Accept a stage-2 pass as a certificate for the example settings ("as reliable as a solve the gate accepts")?
2. The instruments in the toolbox (info.aux; a sossolve option) or harness-only copies?
3. Dual or finite-energy tests ('Qd', 'FQ') for any two- or three-variable example (1d)?
