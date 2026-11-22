# Run log – problem cases (11/20/2026)

Plan: PLAN.md. Worktree C:/wt1120 (aad9606c). Runner pc_qrun.sh; queues in the session scratchpad (qPC1-3.txt);
logs in logs/<TAG>_<n>.log; results in the worktree harness folder (results_<TAG>.csv, exs_<TAG>*.csv).

- 05:10 worktree created; ex_search apply_mods fixed (worktree copy); pc_run.m, pc_sanity.m, pc_jobs.csv (110 jobs),
  pc_starts_50d.csv; ST_STOPMIN 60. Sanity 05:11: 13 functions resolve to the worktree; split keeps [2;1;2];
  settings_2d_stock light and heavy exist; PDE 41 'P' stripped eppos 1e-3: unknown, fr 0.101, numerr 2, 2,272 vars,
  2.6 s. SeDuMi is off the path after lb_setup: pc_run adds it for jobs that name it.
- 05:20 queues started: 1 (W41, W44, WX3, EXS_50d, TRADE3), 2 (-s: W36, W39), 3 (-x after queue 2: W2S, W2Q, W56).
- 05:15 W41 done (42 jobs, 2.5 min): no certificate. MOSEK unknown on 38 (feasibility ratio -0.89 to 0.53, numerr 2
  on most), infeasible on 1 (deg_Q{1} [8 4 4], dD 0); SDPT3 returned no solution (termcode -7, -5); SeDuMi
  unknown (0.42, 0.32). Margins (eppos 1e-3 to 1e-1, epneg 1e-4 to 1e-2) do not move it.
- 05:19 PDE 36 CERTIFIED in W36: light 'Q', dD 0, wcross 1, epneg = 1e-3: fr 1.098, numerr 0, 52,880 variables,
  m 5,643, 6.3 s. Same program at epneg 0 / 1e-4: ratio 0.375 / 7.02 (unknown); eppos 1e-3 to 1e-1: 0.375-0.379.
  Queued W36F (epneg 2e-3 to 1e-1, cheaper neighbours) and EXS_36 (descent from it) ahead of W39.
- 05:30 maintainer: "also try cross-terms if not already". Checked: on 41 and 44 wcross and len_cross{2} leave the
  program unchanged (DIAGA9: same sizes, same verdicts) – inert, no new jobs; X3, 69, 70 had cross terms only
  without the new levers; the two 'P' cross-term starts of 39 never ran (EXU stopped the row after two starts);
  56 with wcross was refused at the EXU m cap (36,021) and never solved; the 50 dual had none. Added WXC (queue 1,
  after WX3: X3 x6, 39 x2), W2X (queue 3, exclusive, after W2Q: 39 'P' x2, 69/70 'P' and 'Q' x4, 56 x3) and three
  cross-term starts of the 50 dual (pc_starts_50d.csv; maxbuild 16).
- 05:35 maintainer: "also try higher degree psatz terms on 69 and 70"; then "run autonomously for 11 hours" (window
  now to ~16:40). Added W2P (queue 3, exclusive, after W2X): on 69 and 70, 'Q' light with dD 0 and the triangular
  cap, with the term set raised to faces plus box quadratics (codes [0 3 4 5 6 7 8]), 'corners' (faces plus
  products of two faces of different directions), box product plus faces ([0 1 3 4 5 6]), and the weight degree
  of every term raised (wmin [3 3 0], [4 4 0]; dw 2); 'P' stripped with dD 0 and the box quadratics, or wmin [4 4 0].
  Not tried on 69/70 before: every earlier start used the tier term sets ('auto', 'faces1') at the tier weights.
- 05:41 W44 done (17 jobs, 26 min): no certificate. Infeasible 11 (Q degree [8 4 4], [2 6 6], [6 6 6] with dD 2;
  margins in 'Q' and 'P'; 'P' storage raises), unknown 3 ([10 6 6] at 1.25M, ratio -0.77; 'P' epneg 1e-3, 0.075;
  'P' veryheavy len_sos [4;4;4], -0.53), SDPT3 no solution, [8 8 8] above the gate (2.41M). The Riemann form (59)
  stays the only certified form of this system.
- 05:4x W36 (14 of 15): a second certificate, at epneg 0: light 'Q', dD 0, wcross 1, deg_Q{2} [0 1 1]: fr 1.189,
  380,804 variables, m 10,791, 44.5 s. deg_Q{2} [2 2 2]: 1.78M, unknown 2.81; psatz{2} 'none': 9.14; SDPT3 no
  solution after 717 s; SeDuMi unknown 0.67 (380,880).
- 05:44 WX3 done (9 jobs, 3 min): infeasible on 8 (Q degree raises – deg_Q{1} does not act on X3, deg_Q{2} 2 gives
  85,383 –, margins, stock 2-D storage in 'P', 'P' heavy), SDPT3 no solution. W36 last job (heavy 'Q', dD 0,
  triangular cap: 1,776,564 variables, m 29,895, about 12 s per MOSEK iteration) running at 05:45.
- 05:48 WXC done (8 jobs, 4 min): X3 infeasible on 5, unknown on 1 ('P' cross terms with epneg 1e-3, -0.50); PDE 39
  'Q' with cross terms, dD 0, default deg_Q{2}: MOSEK optimal, optimality measure 6.67 (epneg 1e-3) and 2.40
  (epneg 1e-2), refused by the gate. W36 last job (heavy 'Q', triangular cap, 1.78M): optimal, measure 2.81.
- 05:50 W36F (3 of 11): epneg 2e-3, 5e-3, 1e-2 give measures 1.90, 1.85, 0.36 – the epneg 1e-3 certificate (1.098)
  is an isolated point of the gate. The gate: sossolve sets feasratio = MSK_DINF_INTPNT_OPT_STATUS for MOSEK; lb_exec
  certifies only |fr - 1| <= 0.3 without numerr. The 11/19 logs show the executive status 'optimal' with eta
  1.8e-8 (36, 52,880, measure 0.375) and 1.3e-8 / 2.4e-9 (39, measures 2.64 / 7.14): by the executive flags and the
  eta gate (1e-7) both rows were accepted on 11/19; only the harness gate refused them (D93: the policy is the
  maintainer's). diag_amt (worktree copy) now records exst and eta per job (header log, (was) lines); W36V reruns
  the 36 settings to record them; queued after EXS_36.
- 06:03 PDE 36 settled under both gates (W36F, W36V; exst and eta recorded): light 'Q', dD 0, wcross 1 with
  epneg 1e-3: measure 1.098, optimal, eta 1.35e-8, 52,880 variables, m 5,643 (stripped gives the same program;
  deg_Q{2} [1 1 1] and the triangular cap leave it unchanged); with epneg 0 and deg_Q{2} [0 1 1]: 1.189, optimal,
  eta 2.37e-8, 380,804; with epneg 1e-1: 0.965, 380,880 (epneg >= 3e-2 builds the larger program). Neighbours fail
  only the measure gate (epneg 2e-3: 1.90, 5e-3: 1.85, 1e-2: 0.36, 3e-2: 1.65; all optimal, eta about 1.7e-8).
  EXS_36 (38 builds) found no smaller program than 52,880. Library path next (LIB_36, the st in the worktree copy
  of the example: dD 0, wcross 1, epneg 1e-3).
- 06:03 EXS_50d (7 of 9 starts): 'cco' infeasible with certificate at every start (stripped 404,021; with cross
  terms 423,137; eppos 1e-2 + dD 0: 96,417; light + triangular cap: 60,645); 'c' light with cross terms (6.03M)
  and 'c' heavy dD 0 (5.21M) above the gate.
- 06:13 LIB_36 (library path, the worktree example with dD 0, wcross 1, epneg 1e-3): stability optimal, 52,880
  variables, 9.6 s, st from the file. EXS_50d done: 9 starts, none certified ('cco' light dD 0 eppos 1e-2: infeasible
  at 1.66M; 'c' light eppos 1e-2: 5.90M, above the gate). TRADE3 started 06:13.
- 06:10 W39 (6 of 9): every job optimal with eta 1.26e-8 to 1.34e-8 at 316,737 variables, measure 2.62-3.30 (margins
  eppos 1e-3 to 1e-1, epneg 1e-4, 1e-3; psatz{2} 'faces1'): refused by the measure gate only.
- 06:26 TRADE3 done (12.5 min): the 80 one-variable rows TRADE1 missed converged; 31 cut, 103,755 -> 62,584 variables
  in sum (DDE13 11,759 -> 5,791; DDE16 23,449 -> 13,486; DDE17 Hinf_control 7,901 -> 4,781; DDE10 1,329 -> 869;
  DDE12 1,509 -> 869; DDE23 3,057 -> 1,812; PDE 21 1,197 -> 681; PDE 20 633 -> 361; NDS01 680 -> 440; NDS12 595 -> 371).
  Optimizations of certified rows, not certifications. EXVT3 (queue 4, beside the exclusive queue) verifies them as
  EXV1B verified Phase C: plan pc_plan_trade3.csv (pc_plan_t3.pl), every test of each plan row under the new st
  (primal and dual share the example st); st_run in the worktree gains the tag (fix_st_run.pl).
- 06:29 EXVT3 done (1.3 min): 46 tests under the 31 cut settings, 37 verified (verify_check.pl), 9 not: the cut of one
  test breaks a test sharing the example st – the duals of PDE 7, 9, NDS01 (infeasible / unknown), the gain duals of
  DDE10, 11, 12 and PDE 21 (+0.6% to +7.2%), the DDE09 primal (+0.12%), the DDE17 stability bisection. 22 example
  settings pass every test: 87,285 -> 53,357 variables in sum (trade3_verified.txt): PDE 32 29,227 -> 13,971; DDE16
  23,449 -> 13,486; DDE13 11,759 -> 5,791; DDE23 3,057 -> 1,812; PDE 52 2,471 -> 1,755; DDE02 2,380 -> 1,964; ...
  Candidates for the maintainer (example files not edited during the comment pass); library path not yet run.
- 06:31 W39 SeDuMi job stopped: the executive raised it to 1,828,161 variables, m 30,203 (nnz(L) 2.3e8, one SeDuMi
  iteration in about 20 min, 24 GB, 11.8 GB free) and it held the exclusive wave behind queue 2. Dropped from the grid;
  the two unrun W39 jobs (deg_Q{2} [1 1 1]; [2 2 1] with the triangular cap) requeued in queue 5. Queue 3 (exclusive
  2-D wave) released at 06:31.
- 06:31 W2S crashed (MATLAB access violation 0xc0000005, the R2025b JIT/evalc crash of memory matlab-evalc-handle-crash)
  in the second pass of PDE 36 'P' with the stock 2-D storage: pass 1 (dw 0) infeasible with certificate at 103,936
  variables (14.9 s), pass 2 built 751,936, crash. 36 moved to tag W2S36 (36 certifies in 'Q'); W2S (69, 70, 39)
  requeued in queue 3 after W2Q.
- 06:43 W39 done (8 jobs; SeDuMi dropped): no job inside the measure gate. 'Q' light, dD 0, wcross 1: deg_Q{2} [2 1 1]
  (316,737) optimal, eta 1.26e-8 to 1.34e-8, measure 2.62-3.30 under every margin and psatz{2} 'faces1'; [1 1 1]
  (316,548) optimal, eta 1.8e-8, 3.51; [2 2 1] with the triangular cap (1,313,577, m 40,377, 643 s) optimal, eta
  1.78e-7, 7.26. With epneg 1e-3 / 1e-2 and the default deg_Q{2} (WXC): 6.67 / 2.40, optimal.
- 06:43 DIAG41 (pc_diag41.m, 25 s; light, 'Q' and 'P'): PDE 41 with the delay as an explicit transport state and the
  Dirichlet coupling x1(1) = x2(1): 'Q' optimal, eta 1.95e-10, measure 0.34 at 2,496 variables (refused by the gate),
  'P' unknown; with the delayed input as a flux x1_s(1) = x2(1): certified in 'Q' (0.996, 2,496) and 'P' (0.979, 4,000);
  decoupled (x1(1) = 0): certified in both (1.005, 1.013); the example file (delay expanded by the parser): unknown in
  both (19,776 / 4,000). The Dirichlet boundary input is what defeats the LPI; the explicit delay state (the file's own
  branch expand_delays = true) brings 'Q' to MOSEK optimal at 2,496. DIAG41B: settings on that form (20 jobs).
- 06:45 DIAG41B (pc_diag41b.m, 36 s): PDE 41 with the explicit delay state CERTIFIED under both gates in 'Q' light at
  2,496 variables, m 199: epneg 1e-2 (1.009, optimal, eta 2.25e-10), epneg 1e-3 (1.034, 8.8e-11), eppos 1e-2 (1.185,
  1.96e-10; epneg 0, the default notion), and dD 2 (0.933, 3,900 variables). Unknown: light/stripped as is (0.34,
  optimal), heavy/veryheavy (numerr 2), deg_Q{1} raises (unchanged program), dD 0; 'P' unknown at every setting
  (numerr 2). The example as written (expand_delays = false, the parser expands the delay: 19,776 variables) is
  inconclusive at every setting (W41, 42 jobs). LIB_41 next: the worktree example with expand_delays = true, eppos 1e-2.
- 06:45 W2Q (69): infeasible with certificate at about 1.0M (deg_Q{2} [2 1 1] with the cap; psatz{2} 'faces1' with
  wcross; wcross with eppos 1e-2); deg_Q{2} 2 with the cap 4.46M, above the gate.
- 06:46 LIB_41 (library path; worktree example with expand_delays = true and st.eppos = 1e-2): stability optimal,
  2,496 variables, m 199, 4.8 s. DIAG44 next (pc_diag44.m): PDE 44 as written, with the state orders swapped (x2
  order 2, x3 order 1), minimal (u, u_t, delayed u_t), and in Riemann states; 'Q', 'P', 'Q' epneg 1e-3, eppos 1e-2.
- 06:48 DIAG44 (pc_diag44.m, 1.6 min): PDE 44 as written infeasible in 'Q', 'P', 'Q' epneg 1e-3, 'Q' eppos 1e-2;
  the state orders swapped (x2 order 2, x3 order 1): the converter refuses it (boundary inputs with more conditions
  than the order); minimal (u, u_t, delayed u_t; no ODE state, no delayed position): infeasible in every form
  (42,323 / 9,000); Riemann (w_s, w_t, delayed w_t): CERTIFIED in 'Q' (0.995, 12,564), 'P' (1.001, eta 4.8e-11,
  2,520), 'Q' epneg 1e-3, 'Q' eppos 1e-2. The position state u defeats the LPI; the energy states certify (= example
  59, the same system). Not a settings matter.
- 06:50 69 and 70 carry the position of the wave (69: u, v = u_t, heat w; 70: phi, phi_t, delta, delta_t); the 1-D
  twin 68 uses energy states and certifies at 324. A 2-D energy-state form does not map onto the PIE format directly
  (e2 = u_s2 has Neumann data only at s2 = 0, 1): not tested; the position-state reading stays an inference for them.
  W2Q (70): deg_Q{2} 2 with the cap 2.02M (above the gate); [2 1 1] with the cap, faces1 with wcross, wcross with eppos
  1e-2: infeasible with certificate (457,585-494,572). DIAGX3 next: X3 in Riemann states.
- 06:51 DIAGX3 (pc_diagx3.m): X3 in Riemann states (w_s, w_t on s1; delayed velocity on s2) CERTIFIED at mu 0.375 and
  0.4: 'Q' 5,214 variables (1.018 / 1.013), 'P' 9,000 (0.999 / 0.997); as written infeasible in both forms. W2Q done:
  69 and 70 infeasible with certificate on every solved job, 'P' stripped dD 0 eppos 1e-2 included.
- 06:55 recheck on 41c9210c: worktree C:/wt1120b (detached 41c9210c), the run files and the two worktree-test example
  edits (36: dD 0, wcross 1, epneg 1e-3; 41: expand_delays = true, eppos 1e-2) copied; queue qPCB (SANITY, W36V,
  LIB_36, LIB_41, DIAG41B, DIAG44, DIAGX3, EXVT3) by pc_qrun_b.sh, logs in logs_41c9210c/.
- 07:01 recheck on 41c9210c (queue qPCB, 9.5 min; lpisettings_sop resolved to C:\wt1120b): every result identical to
  aad9606c – W36V 7/7, LIB_36, LIB_41, DIAG41B 20/20, DIAG44 13/13, DIAGX3 8/8 (status, measure to 3 digits, sizes,
  executive status) and EXVT3 (the same 37 of 46 tests verified). Declan's dpvar subsref and compress change none.
- 07:05 LIB_T3: the 13 PDE example files among the 22 verified TRADE3 cuts (Saba 4.4 left out: per-case st) carry the
  cut after their st block in the 41c9210c worktree copy (pc_t3_files.pl); libcheck runs them at the file defaults.
- 07:03 LIB_T3 (41c9210c, 36 s): all 16 tests of the 13 files optimal at the cut sizes (11: 1,876; 12: 1,629; 18 dual:
  336; 20: 361, gains 0.50046845 / 0.50021558; 22: 301, 0.20843359; 29: 520, 1.0000333; 3: 612; 32: 13,971,
  0.5180508; 5: 469; 52: 1,755, 1.0909044; 6: 401; 67: 825; 68: 297). With EXVT3, 22 example settings are verified
  by the harness, 13 of them also on the library path; the 8 DDE/NDS blocks and Saba 4.4 by the harness only.
- 07:15 BLK_T3: the 8 DDE/NDS blocks among the verified TRADE3 cuts carry the cut after their 11/19 settings lines in
  the 41c9210c worktree library copies (pc_t3_blocks.pl, marked 11/20 "worktree test"); the worktree blockcheck takes
  lines of either marker into the block (one line, annotated). blockcheck runs them through PIETOOLS_DDE.m / _DDF.m.
- 07:06 BLK_T3 (41c9210c, 18 s): 11 tests of the 8 blocks certified with the cut settings received (settings_match
  true) at the cut sizes: DDE02 1,964 and DDE08 2,656 (primal and dual), DDE13 controller 1.3421124 at 5,791, DDE16
  0.9630435 at 13,486, DDE21 estimator 2.328165 at 1,431, DDE22 1.3506131 at 2,141, DDE23 1.8095232 at 1,812, NDS12
  371 (primal and dual). 21 of the 22 verified cuts pass the library path on 41c9210c; Saba 4.4 by the harness only.
- 07:07 W2S: PDE 69 'P' light with the stock 2-D storage infeasible with certificate in both passes (233,856;
  1,691,856, m 26,723, 933 s); PDE 70: pass 1 infeasible (332,704), pass 2 (1,708,624, m 18,226) solving. Queued at
  the end of queue 3: EXS_50d2, the 50 dual in 'c' at stripped (dD 0 with the cap; dD 0; as is), never sized.
- 07:28 W2S done: PDE 39 CERTIFIED in 'P' light with the stock 2-D storage at heavy (settings_2d_stock('heavy')):
  pass 1 (dw 0), measure 0.9998, optimal, eta 5.5e-11, residual 4.6e-6, 1,692,801 variables, m 26,851, 790 s solve.
  Both gates; no gate-policy question. PDE 70 'P' light stock light: infeasible with certificate in both passes
  (332,704; 1,708,624, 607 s). W39S queued after W2X (exclusive): slack cuts on the 39 start (dD 0 with the cap; dD 0;
  stripped; the cap; psatz{2} 'none').
- 08:30 W2X (cross terms: len_cross{2} 2 in 'P', wcross 1 in the slack) done, no case certified. 39 'P' light (the
  container storage) with cross terms: measure -0.9998, executive optimal with eta 0.999 – no solution (238,905); with
  epneg 1e-3 infeasible with certificate (1,696,905, 1,163 s). 69 and 70 'P' light, cross terms, eppos 1e-2: infeasible
  with certificate (1,696,905, 998 s; 1,802,596, 693 s). 69 and 70 'Q' dD 0, deg_Q{2} 2, the cap, wcross 1: above the
  gate (4,456,629, m 70,475; 2,129,860, m 35,909). 56 'Q' light dD 0 wcross 1 with the cap, or with epneg 1e-3:
  infeasible with certificate (1,492,444); 56 'P' light dD 0 with cross terms: infeasible with certificate (471,977).
- 08:32 W39S (slack cuts on the 39 certificate), first two rows: dD 0 with the cap infeasible with certificate (59,841);
  dD 0 alone CERTIFIED – 'P' light, settings_2d_stock('heavy'), dD 0: measure 1.015, optimal, eta 4.39e-8, 335,601
  variables (1,692,801 with dD 1), m 10,572, 52 s solve. Rows 3-5 (stripped; the cap; psatz{2} 'none') running.
- 08:40 queues rewired so the 39 checks run before the long tags: queue 3 ends after W39S; then qPCB2 on 41c9210c
  (pc_qrun_b.sh -x): W39C (the dD 0 row rerun) and LIB_39 (the 41c9210c worktree copy of the example with form 'P',
  st.settings_2d = settings_2d_stock('heavy'), st.lpi.ineq.dD = 0; pc_lib39_edit.pl); then qPC7 on aad9606c (runner
  4, -x): W2P (higher-degree Psatz terms on 69 and 70: box quadratics, corner products, faces with the box product,
  wmin [3 3 0] and [4 4 0], dw 2), W56, EXS_50d2.
- 08:41 W2SH added first in qPC7: the stock 2-D storage at heavy (the storage of the 39 certificate) on 69 and 70,
  never tried (W2S ran them at the stock light level) – with dD 0 and as is; 39 at the stock light level with dD 0
  (whether heavy is needed). pc_jobs_gen.pl regenerated pc_jobs.csv: only the 5 W2SH rows differ.
- 08:58 W39S done (all 'P' light with settings_2d_stock('heavy') unless noted): dD 0 CERTIFIED (335,601, measure
  1.015, eta 4.39e-8) – the cheapest; stripped tier CERTIFIED with the program of light at dD 1 (1,692,801, 0.9998,
  eta 5.5e-11, 889 s); the cap, with or without dD 0, infeasible with certificate (59,841, the same program); psatz{2}
  'none' above the gate (2,987,505, m 46,075: without the box term the slack degrees rise). LIB_39 uses dD 0 as set.
- 08:59 PDE 55 on Sol (job 2067511, 64 threads): stripped 'P' at 7,340,048 variables, m 105,717 (Gram blocks 7 x
  1,024); still solving at 5 h 56 min, 4 min before the 6 h limit; MOSEK prints no iterations into the job log.
- 09:00 W2SH extended: 69 and 70 with the stock storage at veryheavy and dD 0 (7 rows). W39C began on 41c9210c 08:59.
- 09:01 qPCB2 on 41c9210c (3 min): W39C – the dD 0 row identical to aad9606c (certified, measure 1.015, optimal, eta
  4.39e-8, 335,601, m 10,572); LIB_39 – the 41c9210c worktree copy of PIETOOLS_PDE_Ex_2D_NS_Antonelli.m with
  stability_form 'P', st.settings_2d = settings_2d_stock('heavy') and st.lpi.ineq.dD = 0 after its st line: optimal
  through the library path at the same program (335,601, m 10,572, 62 s). PDE 39 certified on both trees and on the
  library path; the three statements are the candidate for the example file (not written: the release comment pass).
- 09:06 ported to the main-checkout harness (uncommitted; both files there were identical to aad9606c, so the worktree
  files were copied): ex_search.m (apply_mods splits at a ';' outside brackets; header entry, (was) line) and
  diag_amt.m (exst and eta columns; header entry, (was) lines). The st_run EXVT3 tag stays in the worktree (it reads
  pc_plan_trade3.csv of this run).
- 09:05 PDE 55 Sol job 2067511: TIMEOUT – cancelled at the 6 h limit (09:03:05) inside the solve of the stripped 'P'
  program (7,340,048 variables, m 105,717), no iterate reported; TotalCPU 13 d 15 h (about 327 core-hours at 64
  threads), MaxRSS 200,704,644 KiB (191 GiB). A longer or larger job needs the maintainer (far above 10 CHE).
- 11:13 qPC7 done (runner 4; W2SH 09:02-10:09, W2P 10:09-10:58, W56 10:58-11:08, EXS_50d2 11:08-11:13):
  - W2SH ('P' light): 39 with settings_2d_stock('light') and dD 0 CERTIFIED (334,656 variables, m 10,129, measure
    1.12, optimal, eta 4.9e-9, 64 s) – heavy is not needed; 69 and 70 with the stock storage at heavy and dD 0
    infeasible with certificate (1,337,841; 686,608); at heavy with dD 1 and at veryheavy with dD 0 above the gate
    (69: 6,766,641 and 6,773,220, m 46,040-46,725; 70: 3,312,720 and 5,598,880).
  - W2P (higher-degree Psatz terms, 'Q' light with dD 0 and the cap unless noted): every solved row infeasible with
    certificate – 69: corners 1,448,852; wmin [3 3 0], wmin [4 4 0], dw 2 all 1,005,460 – the same program: the stability
    executive passes its own dw to lpi_ineq_sop (PIETOOLS_stability_sop.m, the struct o; base 0 in form 'Q'), so
    st.lpi.ineq.dw and st.lpi.ineq.wmin are never read and these rows tested nothing (corrected 11:23). Every
    reported 'Q' program is pass 2 of the executive, the slack weight raised by one (max_raise_Q [2 1 0]; 69:
    pass 1 139,460, pass 2 1,005,460); box quadratics [0 3 4 5 6 7 8] and faces with the box product [0 1 3 4 5 6]
    above the gate (3,929,076, m 40,776; 3,367,908, m 42,907). 70: [0 3 4 5 6 7 8] 1,730,391; corners 667,073;
    [0 1 3 4 5 6] 1,483,278; wmin and dw rows 457,393 (unchanged program). 'P' stripped with dD 0: [0..8] 1,308,969
    (69), 1,214,308 (70); wmin [4 4 0] 336,681 (69), 362,596 (70), the programs without it.
  - W56: 'P' stripped infeasible (828,713); 'Q' stripped with the cap above the gate (2,524,300, m 51,558); 'Q' light
    dD 0 with eppos 1e-2 or epneg 1e-3 infeasible (840,156).
  - EXS_50d2 (the 50 H2 dual in 'c' from stripped starts): every start unknown (66,923; 236,155; 887,267; measure
    -0.38, -0.28, -0.31; eta 1.7e-6, 4.7e-7, 6.2e-6; gam 458, 629, 398); no start certified, no search.
- 11:20 PDE 38 (2-D undamped wave in the position state u, 'Q' light dD 0, deg_Q{2} = [0 1 1]: certified at 380,804;
  11/19 EXS2) shows that the position state alone does not defeat the 2-D LPI. deg_Q{2} = [0 1 1] was never run on 69
  or 70. Queued qPC8 (runner 8, -x, aad9606c): W2E – 69 and 70 'Q' light dD 0 with deg_Q{2} [0 1 1] (with and without
  the cap) and [0 2 2]; then W2D – 70 at deg_Q{2} = 2 with the cap (2,015,785, m 35,390; ndvcap 2.1e6). qPCB3 (runner
  B3, -x, 41c9210c): W39CL – the stock-light 39 row; LIB_39L – the worktree copy of the 39 example with
  settings_2d_stock('light') in place of 'heavy' (pc_lib_39_results.csv keeps the heavy check).
- 11:25 qPCB3 on 41c9210c: W39CL – the stock-light row identical to aad9606c (certified, measure 1.12, optimal, eta
  4.9e-9, 334,656, m 10,129); LIB_39L – the worktree copy of the 39 example with settings_2d_stock('light') and dD 0:
  optimal through the library path (334,656, m 10,129; pass 1 at 46,656 infeasible, pass 2 optimal). PDE 39 has two
  verified statements, stock light or heavy with dD 0 (334,656 or 335,601 variables; measure 1.12 or 1.015).
- 11:23 W2E row 1: 70 'Q' light dD 0, deg_Q{2} [0 1 1], the cap – infeasible with certificate in both passes
  (68,856; 457,225).
- 11:26 the dual forms 'Pd' and 'Qd' were never run on a problem case (no row in any results file of 11/17-11/20).
  PIETOOLS_stability_sop: 'Pd' P >= eppos I, T*P*A' + A*P*T' + 2*epneg*T*P*T' <= 0, every dimension; 'Qd' 1-D only
  (posed as 'Pd' in two or more variables). PDE 18 (1-D damped wave in the position state) certifies only in 'Pd'.
  Queued: W1V on runner 9 (qPC9, -x, aad9606c) – 44 and 41 in 'Pd' and 'Qd' at light and stripped, X3 in 'Pd' at light
  and stripped; W2V on runner 8 ahead of W2D and W2R (qPC8) – 70, 69 'Pd' stripped dD 0, 70, 69, 56 'Pd' light dD 0,
  70, 69 'Pd' light.
- 11:28 W1V done (runner 9, 1.5 min): every dual row infeasible with certificate at every pass. 44: 'Pd' light
  (4,707 / 9,411 / 16,419; numerical spectrum max Re -0.289), 'Qd' light (6,227 / 12,167 / 21,023), 'Pd' and 'Qd'
  stripped (9,411; 12,167). 41: 'Pd' and 'Qd' at light and stripped (4,000 / 7,888; eta 2.0e-2 / 1.3e-2). X3: 'Pd'
  light and stripped (16,614 / 32,406). The dual forms do not resolve 44, 41 or X3.
- 11:30 DIAG69 written (pc_diag69.m; pc_run branch 'DIAG69'): 69 and 70 in one variable, with the position state and
  with the energy states – '68' (the example file, e, v, w), '69pos' (u, v, w along s1 of 69), '70riem' (e = ph_s,
  v = ps, ODE states de, et), '70pos' (ph, ps, de, et along s2 of 70) – each in 'Q', 'P', 'Pd', 'Qd' at light and
  'P', 'Pd' at stripped. Queued on runner 9 after W1V.
- 11:29 DIAG69 done (runner 9, 73 s; results_DIAG69.csv). In one variable the position state defeats the primal
  forms and the dual form 'Pd' certifies it:
  | variant | 'Q' light | 'P' light | 'Pd' light | 'Qd' light | 'P' stripped | 'Pd' stripped |
  |---|---|---|---|---|---|---|
  | 68 (e, v, w) | cert. 12,564, 1.253 | cert. 2,520, 1.037 | cert. 2,520, 1.037 | cert. 3,492, 0.919 | cert. | cert. |
  | 69pos (u, v, w) | infeasible 42,323 | infeasible 9,000 | **cert. 2,520, 1.022** | unknown 3,287, 1.512 (optimal, eta 1.2e-8) | infeasible 5,112 | **cert. 2,520** |
  | 70riem (e, v, de, et) | cert. 6,196, 1.041 | cert. 1,356, 1.167 | cert. 1,356, 0.969 | cert. 1,852, 0.875 | cert. | cert. |
  | 70pos (ph, ps, de, et) | infeasible 23,685 | infeasible 4,428 | **cert. 1,356, 1.003** | cert. 1,705, 0.974 | infeasible 2,604 | **cert. 1,356** |
  The pattern of PDE 18 (1-D damped wave in u: 'Pd' only). In two variables 'Pd' had never run on 69 or 70.
- 11:31 W2V split by size (memory: runner 8 holds 10.9 GB on the 69 [0 1 1] row): W2VS – 70, 69 'Pd' stripped dD 0,
  70, 69, 56 'Pd' light dD 0 – on runner 9 now; W2VL – 70, 69 'Pd' light (dD 1) – first on runner 8 after W2E, then
  W2D, W2R.
- 11:34 W2E rows 1-3: deg_Q{2} [0 1 1] gives 70 the same program with and without the cap (68,856 / 457,225) and
  69 1,005,289 with the cap (139,289 in pass 1); all infeasible with certificate. W2VS rows 1-2: 70 and 69 'Pd'
  stripped dD 0 infeasible with certificate (75,172 / 362,596; 48,681 / 336,681).
- 11:35 PDE 38 declares its velocity x2 of order (2,2) like the position ("x2 must be 2nd order differentiable wrt
  both s1 and s2"); 69 declares v of order (1,2), 70 ps of order (2,1). With equal orders the energy is <Tx,Qx> with
  Q = blkdiag(-(I(x)G2 + G1(x)I), T_v), kernels of degree 1 (inferred, not measured). 70 admits equal orders: the
  second condition of ps in s2 is the time derivative of ph_s2(s1,0) = -et, ps_s2(s1,0) = et + de + ps(s1,0). 69 does
  not (v of order 2 in s1 needs v_s1(0), a third derivative of w). Written: pc_diag70.m (tag DIAG70; variant '70m', in
  'Q' light dD 0 with deg_Q{2} [0 1 1], 'Q' light dD 0, 'Pd' light dD 0); pc_diag44 gains inputs J and tag, and the
  dual mapping (tag DIAG44D: the 44 variants in 'Pd', 'Qd', 'Pd' epneg 1e-3); W2V38 (38 in 'Pd' stripped and light dD
  0, a calibration). Runner 9 queue: DIAG44D, DIAG70, W2V38. Runner 8: W2E, then W2VL, W56D (56 at deg_Q{3} [1 2 2]
  and 2, the analogue of the setting that certifies PDE 43), W2D, W2R.
- 11:36 W2E done (runner 8): deg_Q{2} [0 1 1] infeasible with certificate on 70 (457,225) and 69 (1,005,289), with
  and without the cap (the same programs); [0 2 2] above the gate before any solve – 70 2,015,105 (m 34,675; within
  0.1% of the deg_Q{2} = 2 program that W2D runs at ndvcap 2.1e6, which contains [0 2 2]), 69 4,455,909 (m 70,343).
  W2VS row 3: 70 'Pd' light dD 0 infeasible (77,156 / 364,580).
- 11:37 W2VS done (runner 9): every 2-D dual row with dD 0 infeasible with certificate in both passes – 70 'Pd'
  stripped 362,596, light 364,580; 69 'Pd' stripped 336,681, light 339,705; 56 'Pd' light 266,148 (m 19,916). The 1-D
  result of DIAG69 does not carry to two variables at these sizes; W2VL (dD 1) running on runner 8.
- 11:40 pc_diag38.m (tag DIAG38) in place of W2V38: the wave of PDE 38 without couplings, the velocity of order
  (2,2) ('38', the control), (1,2) ('38v12', as 69) or (2,1) ('38v21', as 70), each in 'Q' light dD 0 with deg_Q{2}
  [0 1 1] and 'Pd' light dD 0. Separates the order of the velocity from the couplings of 69 and 70. Runner 9 queue:
  DIAG44D (running), DIAG70, DIAG38.
- 11:38 DIAG44D done (runner 9, 80 s; results_DIAG44D.csv; the 44 variants of pc_diag44 in the dual forms):
  'minimal' (u order 2, v = u_t order 1, z = u_t(1,t-tau s); no ODE state, no delayed position) CERTIFIED in 'Pd' at
  epneg 0 – 2,520 variables, measure 1.000, optimal, eta 9.3e-10; 'Qd' unknown (measure 0.560, optimal, eta 9.5e-9);
  'Pd' with epneg 1e-3 unknown (measure 1.316, optimal, eta 2.05e-7 – outside both gates). 'example' infeasible in
  'Pd' and 'Qd' (16,419; 21,023), 'Pd' epneg 1e-3 unknown (measure -1.63, numerr 2). 'riemann' certified in all three
  (2,520-3,492). 'orders' does not convert (the boundary condition on the ODE state x1 has more conditions than
  orders). The minimal position form of 44 – the dynamics of the example without its auxiliary states – certifies in
  the dual, where the auxiliary states are the obstruction; in 'Q' and 'P' the minimal form is infeasible as well
  (DIAG44).
- 11:43 DIAG70 done (4.6 min): 70m (ps of order (2,2) with ps_s2(s1,0) = et + de + ps(s1,0)) infeasible with
  certificate in 'Q' dD 0 deg_Q{2} [0 1 1] (457,225, m 13,821), 'Q' dD 0 (457,393), 'Pd' dD 0 (364,580): the order of
  the velocity does not rescue 70.
- 11:45 pc_diag44e.m (tag DIAG44E) queued after DIAG38: 44min in 'Pd' at epneg 1e-4 and 1e-2 (light), stripped (epneg 0
  and 1e-3), heavy (1e-3); x3min (z on s2, the X3 layout) at mu 0.375 and 0.4 in 'Pd' light (epneg 0, 1e-3), stripped.
- 11:45 DIAG38 rows 1-2: the control '38' certified in 'Q' dD 0 [0 1 1] (380,804, measure 1, eta 3.3e-9; pass 1
  52,804 unknown, divergent) and infeasible with certificate in 'Pd' dD 0 (22,980 / 150,980). In two variables the
  dual form is weaker than 'Q' on a wave in the position state; the 1-D result of DIAG69 does not carry over.
- 11:47 DIAG38 done: '38v12' and '38v21' certified in 'Q' dD 0 [0 1 1] like '38' (380,804; measure 1.005, 0.990;
  eta 2.6e-10, 5.5e-10); all three infeasible in 'Pd' (150,980). The order of the velocity does not decide the 2-D
  'Q' test of an uncoupled wave: 69 and 70 fail on their couplings (heat interface at s1 = 0; edge states at s2 = 0).
- 11:48 DIAG44E done (44 s): 44min 'Pd' – light epneg 1e-2 CERTIFIED (2,520, measure 1.274, eta 6.0e-8), heavy epneg
  1e-3 CERTIFIED (14,373, measure 0.844, eta 8.9e-9), stripped epneg 0 certified (2,520, 1.000), light epneg 1e-4
  unknown (1.552, optimal, eta 1.3e-8), stripped epneg 1e-3 unknown (1.316, eta 2.05e-7). x3min (z on s2) 'Pd' at
  epneg 0 CERTIFIED – mu 0.375: 9,000, measure 0.9995, eta 2.4e-10 (light and stripped); mu 0.4: 9,000, 0.9941, eta
  9.5e-10; at epneg 1e-3 inside the measure gate (0.849; 1.101) with eta 1.95e-7 and 3.11e-7, above 1e-7: not
  counted.
- 11:50 DIAG44D rerun on 41c9210c (runner B3): identical to aad9606c in every column but time (minimal 'Pd' certified,
  2,520, measure 1, eta 9.32e-10). Queued: DIAG44E on B3; DIAG44F on runner 9 (x3min 'Pd' light epneg 1e-2 and heavy
  epneg 1e-3 at mu 0.375, 0.4; 44min 'Pd' epneg 3e-2, 1e-1; 'dw1' – u_tt = u_ss, u(0) = 0, u_s(1) = -u_t(1) in the
  position state – in 'Q' light, 'Q' deg_Q{1} [8 6 6], 'P', 'Pd'; 44min 'Q' deg_Q{1} [8 6 6]). pc_diag44f.m.
- 11:53 pc_diag2c.m (tag DIAG2C) queued after DIAG44F on runner 9: the coupling explanation for 69 and 70 measured –
  '70neu' (70 without edge states, ph_s2(s1,0) = 0), '70damp' (ph_s2(s1,0) = ps(s1,0), boundary damping in place of
  the edge states), '69dec' (69 with u_s1(0,s2) = 0 and w(0,s2) = 0), each in 'Q' light dD 0 at [0 1 1] and at the
  default deg_Q. REPORT §1 rows 44, X3 and §2 items 3, 5 updated with DIAG38, DIAG44D, DIAG44E.
- 11:52 DIAG44F done (96 s): x3min 'Pd' light epneg 1e-2 CERTIFIED at mu 0.375 (9,000, measure 1.001, eta 3.2e-10) and
  0.4 (9,000, 1.002, eta 8.2e-10); heavy epneg 1e-3 unknown / inaccurate (63,621, numerr 2). 44min 'Pd' light epneg
  3e-2 and 1e-1 CERTIFIED (2,520; 0.992, 0.987; eta 9.6e-10, 1.3e-9). 'dw1' (u_tt = u_ss, u(0) = 0, u_s(1) = -u_t(1),
  position state): 'Q' light infeasible (18,341; measure -0.99, numerr 2, exst unknown, eta 4.1e-10), 'Q' deg_Q{1}
  [8 6 6] infeasible (155,457; -0.99, numerr 2), 'P' infeasible (4,000), 'Pd' CERTIFIED (1,120, 0.996, eta 7.6e-10).
  44min 'Q' deg_Q{1} [8 6 6] infeasible (354,007; -0.998, numerr 2). The direct forms 'Q' and 'P' do not certify a
  boundary-damped 1-D wave in the position state at any tested degree; the dual 'Pd' does, at 1,120 variables.
- 11:55 The numerical witness of the executive (pie_witness_sop 'rate': largest real part of the finite eigenvalues
  of the Legendre-Galerkin pencil (A,T), N = 24) reports max Re > 0 on every position form built here – dw1 2.365,
  44min 2.149, 69pos 0.844 – and negative values on the energy forms (68 -0.110; 44 and its Riemann form -0.289).
  dw1 has no eigenvalues (u = sinh(lam s): cosh(lam) + sinh(lam) = exp(lam) = 0 has no root), and 69pos, 44min are
  dissipative (dE/dt = -|w_s|^2; boundary terms <= 0 for mu < 1/2), so these values are spurious eigenvalues of the
  discretized pencil, not instability (dw1 exact; the other two by the energy identity, not by a computation).
  Measured next: the witness of dw1 at N = 12, 24, 48, 96 (pc_witness.m).
- 11:57 DIAG44E rerun on 41c9210c (runner B3): identical to aad9606c in every column but time. Queued: WITNESS
  (pc_witness.m) after DIAG2C on runner 9; DIAG44F on B3 (the exponential rows of x3min and 44min, and dw1).
- 12:07 DIAG2C done (14 min): '69dec' (69 without the heat interface) CERTIFIED in 'Q' dD 0 – [0 1 1] 928,729,
  measure 0.888; default deg_Q 928,900, 0.908 (eta 3e-10, 1e-10); '70neu' (70 without edge states, Neumann at s2 = 0)
  CERTIFIED at the default deg_Q (380,880, measure 1.226, eta 2.9e-10), unknown at [0 1 1] (380,804, measure 1.753,
  optimal); '70damp' (ph_s2(s1,0) = ps(s1,0) in place of the edge states) infeasible at both (429,524; 429,600;
  measure -0.9999, numerr 2, exst unknown). Measured: 69 fails on its heat interface; 70 fails on a dissipative
  boundary condition that feeds the position the trace of another state – in 1-D the same structure (dw1) fails in
  'Q' and 'P' and certifies in 'Pd'; in 2-D 'Pd' fails on every wave in the position state.
- 12:07 WITNESS (17 s): max Re of the Legendre-Galerkin pencil at N = 12, 24, 48, 96 – dw1 0.911, 2.365, 4.987,
  10.115 and 44min 0.743, 2.149, 4.767, 9.901 (growing about linearly in N: spurious eigenvalues; dw1 has no
  eigenvalues); 69pos 1.422, 0.844, 1.082, 0.907 (no trend; dissipative by the energy identity); 68 (control, energy
  states) -0.159, -0.110, -0.076, -0.053 (toward 0 from below, as a system that is not exponentially stable). The
  executive prints a positive "numerical spectrum max Re" for these position forms although they are stable: a defect
  of the witness (pie_witness_sop 'rate'), not of the certificates.
- 12:09 DIAG44F rerun on 41c9210c (B3): identical to aad9606c in every column but time.
- 12:12 W2VL done: 70 'Pd' light (dD 1) infeasible (336,820 / 1,712,740; 724 s), 69 'Pd' light infeasible
  (238,905 / 1,696,905; 1,413 s). W56D began (runner 8).
- 12:15 By hand, for dw1 (exact energy E = (|v|^2 + |u_s|^2)/2 in the 'Q' form, T from u(0) = 0, u_s(1) = -v(1),
  v(0) = 0): Q_uu = -I/2 matches |u_s|^2, but the cross term (v(1)) * int theta x_u needs a point evaluation at s = 1
  (a delta) in Q_vu or Q_uv – not in the PI family. A hand computation for one system, not a proof that no other
  storage exists. Where the coupling is a Dirichlet datum of a strain state (the energy forms), no such term arises.
  Built: pc_diag2m.m (tag DIAG2M, runner 9 after WITNESS) – '70mix' (ph (2,0), e2 = ph_s2 (0,1), ps (2,1), de, et;
  e2(s1,0) = -et) and '69mix' (u (0,2), e1 = u_s1 (1,0), v (1,2), w (2,2); e1(0,s2) = -w_s1(0,s2)), each in 'Q' dD 0
  [0 1 1], 'Q' dD 0, 'Pd' dD 0. Both extended systems keep e2 - ph_s2 (e1 - u_s1) constant, so only epneg = 0 applies.
- 12:16 DIAG2M row 1: '70mix' CERTIFIED in 'Q' light dD 0 deg_Q{2} [0 1 1] – 186,809 variables, m 14,376, measure
  0.9996, optimal, eta 8.0e-11, pass 1 (dw 0), 26 s solve. PDE 70 certified (epneg 0) in the mixed position-strain
  form. Rows 2-6 running; pc_diag2m.m copied to wt1120b (41c9210c) with its pc_run branch for the recheck.
- 12:19 W2D stopped at 12:19 (rc 127; its MATLAB process, started by runner 8, ended by hand) before its solve:
  70 is certified in the mixed form at 186,809 variables, and the 2.0M program of the example form held 12 GB and
  growing while 69mix waits for memory (15.8 GB free). W2R removed from the queue. Runner 8 idle.
- 12:20 Library path for 70: pc_lib70_edit.pl edited the wt1120b (41c9210c) copy of
  PIETOOLS_PDE_Ex_2D_Wave_Acoustic_Boundary_Beale.m to the '70mix' states and PDE_t (old lines commented as (was),
  CRLF kept) and added st.lpi.ineq.dD = 0 and st.sop.deg_Q{2} = [0 1 1] after its st line (form 'Q' as in the
  file); pc_lib_70.csv; queued LIB_70 on B3 after the DIAG2M recheck (began 12:19). Runner 9: DIAG2M rows 2-6, then
  DIAG56 (pc_diag56.m: 43 and 56 with the delay state of order 0 in space – '43lib' control, '43thin', '56thin' –
  in 'Q' light dD 0 cap at deg_Q{nv} [1 2 2] and default; W56D: 56 at deg_Q{3} [1 2 2] or 2 is 7.86M variables,
  m 135,522-138,686, not solved).
- 12:23 DIAG2M done (7 min): '70mix' also CERTIFIED in 'Q' light dD 0 at the default deg_Q (187,116 variables,
  m 14,833, measure 1.034, optimal, eta 7.1e-10) – the statement st.lpi.ineq.dD = 0 alone suffices; 'Pd' dD 0
  infeasible (768,065; eta 0.27). '69mix' does not convert: "The boundary conditions along s1 are order
  incompatible" (e1 of order 0 in s2 against w_s1 of order 2 in s2). pc_diag69m.m ('69mix2': e1 of order (1,2) with
  e1 = 0 at s2 = 0, 1, which u = 0 there implies; 'Q' light dD 0, and with deg_Q{2} [0 1 1]) queued on runner 8 (idle).
- 12:24 The 70 worktree copy (wt1120b) keeps only st.lpi.ineq.dD = 0 after its st line (deg_Q{2} [0 1 1] removed
  before LIB_70 starts; pc_lib70_edit.pl and its cleaning copy updated to match). LIB_70 checks the statement to be
  recommended: form 'Q', light, dD 0, the mixed states.
- 12:26 DIAG2M rerun on 41c9210c (B3): identical to aad9606c in every column but time ('70mix' certified at 186,809
  and 187,116; 'Pd' infeasible; '69mix' does not convert). LIB_70 began on B3.
- 12:26 DIAG56 rows 1-5: '43lib' (control) certified at deg_Q{2} [1 2 2] (132,609, measure 0.998, eta 2.9e-9),
  infeasible at the default deg_Q (50,991); '43thin' certified at [1 2 2] (174,209, 0.997, eta 1.2e-10), infeasible at
  the default (75,470); '56thin' at deg_Q{3} [1 2 2] above the gate before any solve (12,781,675, m 101,520). The delay
  state of order 0 in space does not shrink the program (1-D: 174,209 against 132,609; 3-D: 12.8M against 7.86M).
  PDE 56 at the degree that certifies its 1-D twin needs about 7.9M variables (m 135,522): above the program of the
  PDE 55 Sol job that ran out its 6 h limit (7.34M, m 105,717).
- 12:27 LIB_70 (B3, 41c9210c): the worktree copy of the 70 example in the mixed form (pc_lib70_edit.pl) with form
  'Q', st = lpisettings_sop('light') and st.lpi.ineq.dD = 0: optimal through the library path – 187,116 variables,
  m 14,833, eta 7.1e-10, pass 1, 33 s solve. PDE 70 certified (epneg 0) in the mixed form on both trees and on the
  library path; the candidate change to the example is its state declarations and PDE_t (pc_lib70_edit.pl) and the
  one statement st.lpi.ineq.dD = 0 (not written: the release comment pass; a change of formulation, the
  maintainer's call).
- 12:32 DIAG69M done (runner 8, 4 min): '69mix2' (u of order (0,2), e1 = u_s1 of order (1,2) with e1 = 0 at s2 = 0, 1,
  v (1,2), w (2,2); e1(0,s2) = -w_s1(0,s2)) CERTIFIED in 'Q' light dD 0 at epneg 0 – 267,060 variables, m 23,714,
  measure 0.722 (inside the gate, 0.022 from its edge), optimal, eta 4.3e-8, pass 1, 107 s solve; with deg_Q{2}
  [0 1 1] optimal (eta 7.5e-8) but measure 10.71 (unknown). PDE 69 certified in the mixed form. Queued on B3
  (41c9210c): DIAG69M recheck, then LIB_69 – the wt1120b copy of PIETOOLS_PDE_Ex_Wave_Heat_Transmission_2D_Batty.m
  edited by pc_lib69_edit.pl to '69mix2' with st.lpi.ineq.dD = 0 (form 'Q' as in the file); pc_lib_69.csv.
- 12:36 DIAG69N queued on runner 8 (pc_diag69m gains inputs J, tag, T): '69mix2' in 'Q' dD 0 with eppos 1e-3, eppos
  1e-2, at stripped, with the cap, with wcross 1 – a setting with the measure nearer 1 than 0.722.
- 12:33 DIAG56 row 6: '56thin' at the default deg_Q with the cap infeasible with certificate (1,587,084, m 32,521,
  422 s) – larger than 56 as declared (840,156) and infeasible. 56 stays open.
- 12:37 DIAG70N queued on runner 9 (pc_diag2m gains inputs J with a tier per row, tag, V): '70mix' in 'Q' at
  stripped dD 0, light dD 0 with the cap, light dD 0 with psatz{2} 'none', stripped dD 0 with the cap, light as is –
  the cheapest certified setting for the example statement.
- 12:37 DIAG69M rerun on 41c9210c (B3): identical to aad9606c in every column but time ('69mix2' certified at light
  dD 0, 267,060, measure 0.7223, eta 4.27e-8; [0 1 1] measure 10.71).
- 12:39 LIB_69 (B3): the worktree copy of the 69 example in the mixed form (pc_lib69_edit.pl) with form 'Q',
  st = lpisettings_sop('light') and st.lpi.ineq.dD = 0: optimal through the library path – 267,060 variables,
  m 23,714, eta 4.3e-8, pass 1, 95 s solve. PDE 69 certified (epneg 0) in the mixed form on both trees and on the
  library path; candidate change as for 70 (state declarations and PDE_t of pc_lib69_edit.pl; st.lpi.ineq.dD = 0).
- 12:41 W55 queued on runner 8 after DIAG69N: PDE 55 'P' stripped dD 0 with st.sop.dw_feas [0 1 2] and [0 1 3] – the
  slack weight in three variables, never raised (the 11/19 row with st.lpi.ineq.wmin [0 0 1] changed nothing: wmin is
  not read by the executive). P = I is an energy certificate for 55 at lam = 1 (C2*3*pi^2 = 8.9 > lam), so a
  low-degree storage exists; the question is the slack.
- 12:49 W55 done (6 min): 55 'P' stripped dD 0 with dw_feas [0 1 2] infeasible with certificate (284,944; Gram blocks
  [432 128 x6]), with [0 1 3] infeasible (684,304; [432 288 x6]). The slack weight does not resolve 55 at dD 0.
- 12:54 DIAG69N done (runner 8): '69mix2' with eppos 1e-3, eppos 1e-2, at stripped, with the cap, with wcross 1 – the
  same program each time (267,060, m 23,714), measure 0.7223-0.7226, eta 4.3e-8: certified, reproducible, and no
  knob tried moves the measure off 0.722.
- 12:54 DIAG70N done (runner 9): '70mix' at stripped dD 0, light dD 0 with the cap, stripped dD 0 with the cap – the
  program of light dD 0 (187,116, measure 1.034); psatz{2} 'none' unknown (707,680, measure 0.316, numerr 2); light
  as is (dD 1) unknown (758,364, measure 0.552, optimal, eta 1.1e-9). The statement st.lpi.ineq.dD = 0 is needed and
  sufficient at light or stripped.
- 12:55 W55C queued on runner 8: 55 'P' stripped at dD 1 with cap_inc [2 2 2], [3 3 3], [4 4 4] (one cap per
  direction; 'auto' gave 1,808 infeasible, no cap 7.34M).
- 12:57 DIAG69P queued on runner 9: '69mix2' in 'Q' light as is (dD 1), light dD 0 with psatz{2} 'faces1', heavy
  dD 0 – other programs, for a measure nearer 1 than 0.722.
- 12:59 W55C done (3 min): 55 'P' stripped at dD 1 with cap_inc [2 2 2], [3 3 3], [4 4 4] – each the uncapped program
  (7,340,048, m 105,717; above the gate): no program lies between the 'auto' cap (1,808) and no cap. DIAG69P:
  '69mix2' with psatz{2} 'faces1' the dD 0 program again (267,060, measure 0.7223); light at dD 1 above the m gate
  (1,150,980, m 54,558); heavy dD 0 above the gate (8,340,061, m 123,647). Every local program of 69mix2 gives 0.722.
- 13:01 pc_diag56t.m (tag DIAG56T) queued on runner 8: the delay sweep tau = f*tau* – 43 at f = 0.5, 0.75, 0.9 ('Q'
  light dD 0 with the cap, default deg_Q) and 56 at f = 0.25, 0.5, 0.75 ('Q' light dD 0, default deg_Q) – how far
  below the threshold the default degree certifies (the library delays: 43 at 0.966 tau*, 56 at 0.9 tau*).
- 13:13 Library path for 44 and X3: pc_lib44x3_edit.pl edited the wt1120b copies of
  PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay.m (44) and ..._2D.m (X3) to the minimal form (u, v = u_t, delayed v; z on
  s2 for X3; old lines as (was), CRLF kept), the dual test (stability_dual = 1 and stability_form 'P', i.e. 'Pd', in
  place of stability = 1 and 'Q'; PIETOOLS_auto_execute_sop runs form stability_form followed by 'd' for that
  flag) and st.epneg = 1e-1 (44) or 1e-2 (X3) after the st line; pc_lib_44.csv, pc_lib_X3.csv; LIB_44 and LIB_X3
  queued on B3. REPORT §4 (all finished jobs), §5 items 1-2 and §6 (open cases, candidate statements, defects) written.
- 13:04 LIB_44, LIB_X3 (B3, 41c9210c): the worktree copies of the 44 and X3 examples in the minimal form with the dual
  test (stability_dual, form 'P') and st.epneg 1e-1 (44) or 1e-2 (X3): stability_dual optimal through the library
  path – 44: 2,520 variables, m 315; X3: 9,000, m 291. 44 and X3 certified (exponential, dual) in the minimal form
  on both trees and on the library path.
- 13:11 DIAG56T done (10 min): at the default deg_Q, 43 ('Q' light dD 0 with the cap) infeasible at tau = 0.5, 0.75,
  0.9 tau* (50,991) and 56 ('Q' light dD 0) infeasible at 0.25, 0.5, 0.75 tau* (840,156). The default degree fails
  at every delay tried, not near the threshold only: both systems decay at rate 0.1 without the delay (first mode
  r - 1 - 2 = -0.1 for 56; a0 - 1 - 1 = -0.1 for 43) against a delayed feedback of gain 1. 56 needs the degree raise
  of 43 (deg_Q{3} [1 2 2], 7.86M variables).
- 13:13 W50S queued on runner 9: the stability test of the PDE 50 system in 'Pd' light, 'Pd' stripped dD 0, 'Qd'
  light and 'P' light dD 0. Reason: the H2 dual form 'cco' (W >= eppos I, (A*W)*T' + T*(W*A') + B1*B1' <= 0) is
  infeasible at every start, and any such W is a 'Pd' certificate (P = W: T*P*A' + A*P*T' <= -B1*B1' <= 0), so 'cco'
  can be feasible only where 'Pd' is (a two-line argument, not yet measured on 50).
- 13:18 W50S done (5 min): the PDE 50 system certifies in 'P' light dD 0 (98,080 variables, measure 1.003, optimal,
  eta 2.6e-10, pass 1); 'Pd' stripped dD 0 infeasible with certificate (96,416 / 501,968); 'Pd' light and 'Qd'
  light infeasible in pass 1 (405,684; 'Qd' posed as 'Pd': two variables) and above the gate in pass 2 (2,158,560).
  The H2 dual form 'cco' needs a coercive W that satisfies the 'Pd' inequality, and 'Pd' is infeasible here: inferred
  (the two programs parameterize W and P with their own degrees), the dual H2 failure is that of the dual stability
  test, as for PDE 38 (two variables, 'Pd' infeasible where the
  primal certifies).
- 13:24 pc_diagpd.m (tag DIAGPD) queued on runner 9: 'P' against 'Pd' (light dD 0) on self-adjoint PIEs, where the
  two conditions coincide (T, A self-adjoint for u_t = c*Lap(u) + r*u with Dirichlet data): 'heat1d' (control),
  'heat2d', 'rd2d' (r = 15, threshold 2*pi^2). A 'Pd' failure where 'P' certifies would point to the two-variable
  'Pd' path, not to the systems.
- 13:21 DIAGPD done (18 s): on self-adjoint PIEs 'P' and 'Pd' build the same programs with the same results –
  heat1d certified (180, measure 1.027), heat2d certified (5,745, 1.096, eta 5.4e-10), rd2d (r = 15, 76% of the
  threshold 2*pi^2) infeasible in both (5,745 / 37,745). The two-variable path of 'Pd' agrees with 'P' where the
  conditions coincide: its failures on 38, 50, 56, 69, 70 are a property of the dual form on those (non-self-adjoint)
  systems, not an implementation defect. (Side result: the 2-D reaction-diffusion equation at r = 15 is not certified
  at light dD 0 in either form.)
- 13:26 Storage-degree cuts for the example statements: DIAG70T on runner 9 ('70mix' 'Q' light dD 0 with deg_Q{1}
  [2 2 2], [2 1 1] – the 1-D edge-state blocks, [4 2 2] at light –, deg_Q{2} 0, and deg_Q{1} [2 2 2] with deg_Q{2}
  [0 1 1]); DIAG69T on runner 8 ('69mix2' 'Q' light dD 0 with deg_Q{2} 0 and [1 0 0]).
- 13:42 DIAG70T: '70mix' with deg_Q{1} [2 2 2] or [2 1 1] the program of light dD 0 (187,116; the 1-D storage
  degree does not change it), with deg_Q{2} 0 infeasible (108,705), with deg_Q{1} [2 2 2] and deg_Q{2} [0 1 1] the
  [0 1 1] program (186,809). DIAG69T: '69mix2' with deg_Q{2} 0 infeasible (843,604) and [1 0 0] infeasible
  (1,924,841) – lower storage degrees give larger, infeasible programs (the slack is sized from the storage). The
  example statements stay st.lpi.ineq.dD = 0 for both.
- 13:43 W43P queued on runner 9: PDE 43 in 'P' (light dD 0, light, heavy dD 0, stripped) – a 'P' route for its 3-D
  twin 56 would show here first.
- 13:48 W43P done (6 min): PDE 43 in 'P' infeasible at light dD 0 (91,145), light (428,185), stripped (427,689),
  unknown at heavy dD 0 (1,508,273; measure -1.20, numerr 2). No 'P' route for 43, so none for 56; 56 stays the
  'Q' program at deg_Q{3} [1 2 2] (7.86M).
- 13:50 DIAGRD queued (pc_diagrd.m): the side result of DIAGPD – the 2-D reaction-diffusion equation at r = 15 (76%
  of 2*pi^2) not certified at light dD 0 – at 'P' light (dD 1), 'P' heavy, 'Q' light dD 0, 'Q' light, and at r = 10
  and 5 in 'P' light dD 0.
- 13:51 DIAGRD done (side result, not a problem case): the 2-D reaction-diffusion equation at r = 15 (76% of 2*pi^2)
  certifies in 'P' light at dD 1 (188,545, measure 0.968), 'P' heavy (752,976), 'Q' light dD 0 (16,740, measure
  0.985) and 'Q' light (72,100); 'P' light dD 0 certifies r = 5 (37,745) and not r = 10. dD 0 in 'P' is the
  restriction there; 'Q' dD 0 certifies at a ninth of the 'P' program.
- 13:54 Final check on the current HEAD. The main checkout is at f474e7de (the release comment pass, 7 commits after
  41c9210c, 05:31-05:33, 415 files under the code folders); every verification so far ran on 41c9210c. Worktree
  C:/wt1120c (detached f474e7de) with the harness of this run (pc_*.m, the patched diag_amt.m and ex_search.m,
  pc_jobs.csv with two VER rows, pc_lib_*.csv; no results files carried over) and the example edits applied by the
  same scripts (pc_lib39_edit.pl with the stock-light statements, pc_lib70/69_edit.pl, pc_lib44x3_edit.pl; 36: dD 0,
  wcross 1, epneg 1e-3; 41: expand_delays = true, eppos 1e-2). Runner C (pc_qrun_c.sh, logs in logs_f474e7de),
  queue qPCC: VER (39, 36), DIAG2M, DIAG69M, DIAG44E, DIAG44F, LIB_39L, LIB_36, LIB_41, LIB_70, LIB_69, LIB_44,
  LIB_X3. candidate_patches/: the diffs of the 41c9210c worktree copies (39, 44, X3, 69, 70, and every example edit).
- 14:07 Final check done (runner C, 13:53-14:07): on f474e7de every candidate is identical to the earlier trees in
  every column but time – VER (39 stock light dD 0: certified, 334,656, measure 1.12; 36 dD 0 wcross 1 epneg 1e-3:
  certified, 52,880, 1.098), DIAG2M, DIAG69M, DIAG44E, DIAG44F (to aad9606c), and the library path optimal for 39,
  36, 41, 70, 69, 44, X3 at the same sizes (LIB_39L, LIB_69, LIB_70, LIB_44, LIB_X3 to 41c9210c; LIB_36, LIB_41 to
  aad9606c). The candidates hold on the release code.
- 14:09 TRADE3 on f474e7de: wt1120c receives pc_plan_trade3.csv, verify_EXVT3.csv, pc_lib_T3.csv and the run copies of
  st_run.m (EXVT3) and blockcheck.m (either marker); pc_t3_files.pl set the cut in the 13 PDE example files and
  pc_t3_blocks.pl in the 8 DDE/NDS blocks (both scripts applied cleanly to the commented release files); pc_run gains
  the BLK_T3 branch. Queue qPCC2 (runner C): EXVT3, LIB_T3, BLK_T3, STOP.
- 14:11 TRADE3 on f474e7de done (2 min): EXVT3 identical to 41c9210c (46 tests: library, id, settings, status, result,
  documented, variables), LIB_T3 identical (16 tests of the 13 files, id..m), BLK_T3 identical (11 tests of the 8
  blocks, feasibility ratios included). Every candidate of this run – the 7 problem-case statements and the 22
  TRADE3 cuts – holds on the release code.
- 14:16 Wrap-up: runners 8, 9, B3 stopped (STOP; no MATLAB -batch process left). Copied into results/<commit>/
  (aad9606c 93 files, 41c9210c 54, f474e7de 51): the run files of each worktree harness (untracked and modified –
  pc_*.m, pc_jobs.csv, pc_lib_*.csv, results_*.csv, exs_*.csv, the patched diag_amt.m, ex_search.m, st_run.m,
  blockcheck.m) and examples_worktree.patch (every example edit of that worktree). pc_jobs_gen.pl, pc_qrun*.sh and
  pc_bycase.pl copied beside them. Worktrees C:/wt1120, C:/wt1120b, C:/wt1120c removed (git worktree remove --force;
  pruned); C:/wt1117 and the master-wt of another session left as they are.
- 14:17 Late checks (the window runs to 16:40; feedback 11/19: keep the runners busy): worktree C:/wt1120d at f474e7de
  with the harness of results/f474e7de (no results files), runner D (pc_qrun_d.sh, logs_f474e7de_d), queue qPCD:
  SND70, SND69 – the mixed-form storages at epneg 1e-3 and 1e-2 must fail (69 not exponentially, 70 not uniformly
  stable; a certificate there would make the rewrite suspect); DIAG44G (pc_diag44g.m) – 44min and the Riemann form
  in 'Pd' at mu = 0.45, 0.48, 0.49, 0.55, 0.6 (stable whenever mu < 1/2; the Riemann witness as reference), x3min
  at 0.45, 0.55; EXS_50d3 – the 50 H2 dual in 'c' at light (dD 0 with and without the cap) and heavy dD 0, never
  run (pc_starts_50d3.csv).
- 14:35 Late checks done (runner D, 18 min). SND70: '70mix' at epneg 1e-3 and 1e-2 infeasible with certificate
  (1,288,544): no exponential certificate for 70 (not uniformly stable) – consistent. SND69: '69mix2' at epneg 1e-3 and
  1e-2 infeasible at the certified program (267,060, pass 1); pass 2 above the m gate (1,926,185, m 45,866) –
  consistent. DIAG44G: 44min and the Riemann form certified in 'Pd' at mu = 0.45, 0.48, 0.49 (2,520), infeasible at
  0.55 and 0.6 (9,000); the Riemann witness -0.228, -0.194, -0.182, -0.117, -0.066 – at tau = 1 the system is still
  stable past mu = 1/2 (the bound holds for every delay), so the failures at 0.55, 0.6 are conservatism, and the
  soundness test is not decided there; x3min certified at 0.45, infeasible at 0.55. EXS_50d3: the 50 H2 dual in 'c'
  unknown at light dD 0 with the cap (242,485; measure -0.19; gam 508) and without (1,585,440; -0.30; gam 305.5),
  above the gate at heavy dD 0 (5,209,966).
- 14:37 DIAG44H queued (runner D, pc_diag44h.m): mu = 0.65, 0.7, 0.8, 0.9, 0.95 – the Riemann witness at N = 24, 48
  and 44min, riemann in 'Pd': a decisive soundness test where the witness turns positive.
- 14:37 DIAG44H (runner D): MATLAB segmentation fault (rc 139) after 7 rows, at mu = 0.9; the rows in hand decide the
  test. By hand (k = 1, tau = 1; u = sinh(lam s), z = exp(lam)): the characteristic equation of 44 reduces to
  (2 - mu) z^3 + mu z^2 + mu z - mu = 0, Re lam = ln|z|; it gives -0.289 at mu = 0.4 and -0.0163 at mu = 0.65, the
  Riemann witness values (DIAG44D, DIAG44H: -0.0162 at N = 24 and 48), and its complex pair reaches |z| = 1 exactly
  at mu = 2/3 (real root mu/(2 - mu); 2w^2 + w - 1 = 0 with w = mu/(2 - mu)). So at tau = 1 the system is
  exponentially stable iff mu < 2/3. The witness at mu >= 0.7 grows with N (38 / 144 at N = 24 / 48): a neutral
  system, eigenvalues of equal real part at every frequency, resolved by larger N. 44min and riemann in 'Pd':
  infeasible at 0.65 (stable, conservatism), 0.7, 0.8, 0.9 (unstable). With DIAG44G: the dual certificate holds for
  mu <= 0.49, fails from 0.55, and never certifies an unstable case.
- 14:46 DIAG70S (runner D, 8.5 min): '70mix' with the edge damping reversed, et_t = +c*et - de - ps(s1,0) (dE/dt =
  c*|et|^2), infeasible with certificate at c = 1 and c = 0.1 (pass 2, 1,288,544); the control c = -1 (the example)
  certified (187,116, measure 1.034). The mixed form does not certify an energy-increasing variant.
- 14:49 DIAG69S (runner D, 2.5 min): '69mix2' with w_t = w_s1s1 + w_s2s2 + 30w (r = 30 > 2*pi^2) infeasible with
  certificate at the certified program (267,060, pass 1); pass 2 above the m gate (1,926,185, m 45,866); the control
  r = 0 certified (267,060, measure 0.722). Every soundness check of the rewritten forms failed where it must.
- 14:52 REPORT read through: finding 8 (soundness), the §4 rows of 44, 69, 70, 50 dual and 55, §5 item 1, §6 brought
  up to date; the 44 boundary stated as the sign of the spectral abscissa. Queued on runner D (qPCD5): DIAG44I
  (pc_diag44i.m: 44min 'Pd', riemann 'Pd' and 'P' at heavy and veryheavy, mu = 0.55, 0.6, 0.65 – the range between
  the certified 0.49 and the boundary 2/3) and DIAG70P (pc_diag70p.m: '70mix' at edge damping and stiffness (d, kk) =
  (0.1, 1), (10, 1), (1, 10), (1, 1) – robustness of the certificate).
- 14:56 DIAG44I (2 min): at mu = 0.55, 0.6, 0.65 every row infeasible – 44min 'Pd', riemann 'Pd' and riemann 'P' at
  heavy (20,853) and veryheavy (113,364; numerr 2 on most). The certificates of 44 reach mu = 0.49 at every tier and
  stop before 0.55: they cover the delay-independent bound mu < 1/2 and none of the delay-dependent range
  1/2 < mu < 2/3 at tau = 1.
- 14:56 DIAG70P (2 min): '70mix' certified at every edge oscillator setting – (d, kk) = (0.1, 1) measure 1.038,
  (10, 1) 1.002, (1, 10) 1.138, (1, 1) 1.034 – each at 187,116 variables, light dD 0, eta 4e-11 to 1.1e-8.
- 14:58 PDE 49, a one-variable system, certified its dual H2 test in 'c' at veryheavy on 11/19 (0.10157588, 11,847):
  the 'c' form works in one variable; on 50 (two variables) it fails as the dual stability test does (REPORT §6).
- 14:59 Late-check files copied to results/f474e7de_late/ (43 files: the run harness, results of SND70, SND69,
  DIAG44G, DIAG44H, DIAG44I, DIAG69S, DIAG70S, DIAG70P, EXS_50d3); worktree C:/wt1120d removed; pc_qrun_d.sh copied.
  No MATLAB batch process left. End of the run: nothing committed, nothing pushed; no example or code file of the
  main checkout changed by this run except the harness port of 09:06 (ex_search.m, diag_amt.m, uncommitted) and the
  records under cleaning/private/.
- 14:58 sol_drafts/ (NOT submitted; above 10 CHE, the maintainer's word needed): sol55b.slurm (55 'P' stripped,
  general partition, 64 cores, 256 GB, 48 h; cap 6,144 CHE) and sol56.slurm, sol56.m, sol56_starts.csv (56 'Q' light
  dD 0 with the cap at deg_Q{3} [1 2 2]; 64 cores, 384 GB, 48 h; cap 7,680 CHE); README with the basis (job 2067511:
  no iterate in 6 h; one dense Cholesky of the Schur complement about 3.9e14 flops for 55, 8.3e14 for 56) and the
  note to add a MOSEK log file before resubmitting.
- 14:59 SeDuMi confirmation (feedback 11/19: confirmation runs need no decision): worktree C:/wt1120e at f474e7de with
  the harness of results/f474e7de_late; pc_diagsed.m (SeDuMi 1.3 on the path after lb_setup; solver 'sedumi'): 44min
  'Pd' light at epneg 0 and 1e-1, x3min 'Pd' light at epneg 0 and 1e-2 (mu = 0.4), then 70mix 'Q' light dD 0
  (187,116 with MOSEK; SeDuMi time unknown, the 39 SeDuMi job of this morning took about 20 min per iteration at
  1.8M). Runner E (pc_qrun_e.sh, logs_f474e7de_sedumi).
- 15:00 DIAGSED rows 1-4 (SeDuMi 1.3): 44min 'Pd' at epneg 0 and 1e-1 and x3min 'Pd' at epneg 0 and 1e-2 – the
  MOSEK programs (2,520; 9,000), feasibility ratio 0.983 to 1.000, eta 5.6e-10 to 2.5e-9, status inaccurate (SeDuMi
  numerr 1: accuracy flag), no infeasibility certificate. SeDuMi does not contradict the MOSEK certificates; by the
  harness gate (numerr 0) it certifies none of them. Row 5 (70mix) running.
- 15:08 DIAGSED row 5: '70mix' with SeDuMi – the MOSEK program (187,116), feasibility ratio 0.984, status inaccurate
  (numerr 1), eta 5.3e-7 (residual 6.1e-4), 500 s solve. On all five programs SeDuMi reaches feasibility ratios 0.98
  to 1.00 and no infeasibility certificate; at its default accuracy none passes the harness gate (numerr 0; eta
  1e-7 for 70mix).
- 15:21 SED36, SED69 (SeDuMi): 36 ('Q' light dD 0 wcross 1 epneg 1e-3) divergent in pass 1 (52,880; eta 4.1e-2) and
  pass 2 (380,880; eta 3.0e-5; measure 0.765); 69mix2 divergent in pass 1 (267,060; eta 5.9e-2; 446 s), pass 2 above
  the m gate. SeDuMi gives no information on the two-variable 'Q' programs of 36 and 69 (as for 39 this morning);
  on 44, X3 and 70 it reaches the MOSEK programs at ratio 0.98-1.00. No program gave an infeasibility certificate.
- 15:22 End of the experiments: results/f474e7de_sedumi/ (34 files); worktree C:/wt1120e removed; no MATLAB batch
  process left. sossolve passes no MOSEK parameters (mosekopt('minimize info', prob)), so a tolerance test of the 69
  certificate (measure 0.722) needs a code change, outside this run. Everything left needs the maintainer: the Sol
  jobs (sol_drafts/), the example rewrites (candidate_patches/), the dual H2 form of 50.
- 16:34 Maintainer answered "Not now" to the 48-hour Sol jobs for 55 and 56: the drafts stay in sol_drafts/,
  not submitted. No local run started: the window closes at 16:40.
- 16:39 Maintainer: "Commit the harness fix". Committed as add6e3ef on ndopvar (ex_search.m and diag_amt.m only),
  on top of 8382d71a (another session merged origin at 16:38); not pushed.
- 17:42 Maintainer: "Push it". add6e3ef pushed to origin/ndopvar (fast-forward from 8382d71a).
- 19:20 Maintainer, on "only after I rewrote their equations": "Is this justified?" REPORT corrected: "not certifiable
  as written" and "at every setting" now read "not certified as written at any setting tried" and "at every setting
  tried" (8 places, §1 and §6); §2 item 9 added: certificates of the mixed forms transfer to 69 and 70 as written
  (every solution of the example, with the strain appended, solves the mixed form); for 44 and X3 the dropped states
  x1, x4 are not covered (their decay needs a bound on u_s). Item 6: dw/wmin unread by design since 11/10, undocumented.
- 19:22 Finite-energy tests (maintainer: "maybe test if there is a certificate using finite energy stability test?"),
  worktree C:/wt1120f (17d68747), pc_fe.m, logs_fe/: 'FP' (P >= 0, T'PA + A'PT + 2 epneg T'T <= 0) and 'FQ' at
  epneg 1, and the dual H2 'cco' of 50 at eppos 0 (W >= 0). By 19:27: 55 'FP' stripped dD 0 infeasible (114,704);
  43 (1-D twin of 56) 'FQ' and 'FP' light, dD 1 and 0, infeasible (428,185; 91,145); 50 'cco' eppos 0 stripped
  infeasible (16,337 with the cap; 96,417). Running: 43 heavy, 56 'FP' stripped at 0.9 and 0.5 tau*, 50 'cco' light.
- 19:32 Finite-energy queue A done: 43 'FQ' and 'FP' at heavy above the size gate (3,497,665 variables, m 10,706;
  not solved); 56 'FP' stripped dD 0 infeasible with certificate at 0.9 tau* and at 0.5 tau* (258,409 each, 31-36 s).
  Every finite-energy row solved so far is infeasible with certificate. Still running: 50 'cco' eppos 0 at light.
- 19:34 Finite-energy queue C done: 50 'cco' at eppos 0 at light infeasible with certificate in both passes (503,633;
  1,660,065). No finite-energy row certifies (REPORT §2 item 10). Files in results/fe_17d68747/ (pc_fe.m, three
  CSVs); worktree C:/wt1120f removed. Nothing committed.
- 11/21 00:51 Correction: for a PIE in two or more variables PIETOOLS_stability_sop poses 'FP' when 'FQ' is requested
  (line 146, no warning; the printed label still reads 'FQ'). 43 has two variables (space, delay), so its 'FQ' rows
  repeat the 'FP' program (428,185 / 91,145 variables, 4,626 / 2,019 rows in both). 'FQ' was not tested on any open
  case. REPORT §2 item 10 corrected.
- 11/21 01:05 Correction to 12:15 (11/20): "not in the PI family" was wrong (maintainer: "the subs operation acting on
  T should provide this operator"). u(1) = subs(T_u,s,1) x is a PI operator on the fundamental state, and
  E = <Tx,Qx> + u(1)(Q_b x) exactly with Qx = (-x_u, T_v x)/2, Q_b = -subs(T_v,s,1)/2. What holds: no linear Q into L2
  gives <Tx,Qx> = E (with v = 0, u -> u(1) is unbounded in L2) – a limit of the parametrization T'Q. REPORT §2 item 7
  corrected. Then: fixed-storage slack check (maintainer: "fix the storage output of the previous result and just
  find the slack"), and the projection check.
- 11/21 01:16-01:36 A posteriori check (side session; maintainer: "fix the storage output of the previous result and
  just find the slack"; projection check: "you can try it"). pc_post.m on 17d68747, worktree C:/wt1121v (removed after
  the runs; instruments and results in results/post_17d68747/, logs in logs_post/). First control runs failed for two
  mechanical reasons, both fixed: the stage-2 slack was pruned on the numeric target (1,352 of 1,568 entries), and a
  storage fixed off the rows that involve it alone left D with components no slack reaches (sossolve: infeasible
  before the solve). Final: controls pass (case 41 explicit delay 'Q' eppos 1e-2; 44min 'Pd' epneg 0), the 44min
  storage against mu = 0.8 is infeasible. Pass the gate with the storage fixed: case 41 'Q' light (stage-1 measure
  0.34), case 36 'Q' at epneg 0 (0.375; 52,880) and at epneg 1e-3 (1.52), 69mix2 (10.7), 44min 'Qd' (0.56), 44min 'Pd'
  epneg 1e-3 (on a repaired storage), the four SeDuMi rows in both stages. Fail: 44min 'Pd' epneg 1e-4, x3min heavy
  epneg 1e-3 (both mu), case 50 dual H2 'c' at three starts (open case stays open). Projection: negative lambda_min
  on every row, controls included (Grams at 1e-12 to 1e-9 before the step). REPORT §2 item 11. Case 50 is
  PIETOOLS_PDDE_Ex (not PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion). Nothing committed.
