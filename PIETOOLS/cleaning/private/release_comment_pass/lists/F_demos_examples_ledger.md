# Ledger – group F_demos_examples (11/19/2026)

57 files read, 57 edited, comment-only. codeident.sh on all 57 at the end: **57 OK, 0 DIFF, 0 NEW**.
Line endings kept (CRLF files stay CRLF; heatNd_apply.m stays LF). No 12/31/2021 import dates in this
group; the GPL block is copied from HEAD of executives_sopvar/PIETOOLS_stability_sop.m.

## Edited files

Common to every file below: "For support" line added where missing; new = 0 files get the log line
`% MMP, 11/19/2026: Release comment pass (comment-only).`; new = 1 files lose all dated log entries, per-line
markers and (was) lines, keeping the Initial coding line.

DEMOs (PIETOOLS_demos/)
- DEMO1_stability_and_simulation.m (0): manual pointer removed; file name in line 1 corrected; header model
  corrected to the code (no input u in the ODE, -b*X2, BCs X2(0) = 0 and X1(1) = x, output r); "design
  controllers" dropped (the demo tests stability only); paragraph on steps and workspace; step comment.
- DEMO3_poincare_inequality.m (0): manual pointer removed; steps/workspace paragraph; step comments.
- DEMO4_stability_parameter_bisection.m (0): manual pointer removed; bisection paragraph; step comments.
- DEMO5_Hinf_optimal_estimator.m (0): manual pointer and pointer to the retired PIETOOLS_Hinf_estimator
  replaced by the commented lpiscript_sop call; steps/workspace paragraph; step comments.
- DEMO6_Hinf_optimal_control.m (0): as DEMO5; PDE inputs s*w, s*u as in the code; prime added in the LPI
  block T*(C1*P+D12*Z)'; comment "don't plot final solution" corrected to "plot the solution".
- DEMO7_observer_based_control.m (0): manual pointer removed; inputs s*w, s*u and output y added to the
  header; sign -L*D21 in the closed-loop formula; plot comment corrected; steps/workspace paragraph.
- DEMO8_H2norm.m (0): manual pointer removed; file name corrected; description corrected (the demo bounds the
  H2 norm only); LPI formula and outputs added; step comments.
- DEMO9_L2_gain_analysis.m (0): manual pointer removed; LPI block matched to the code (it was the code's LPI
  up to the congruence diag(-I,I,I), with the typo P*B'); steps/workspace paragraph; step comments.
- DEMO10_PBC.m (1): 5 log entries, markers, (was) lines and a dead commented block (undefined F) removed;
  steps/workspace paragraph; step comments.
- DEMO11_poincare_inequality_2D.m (1): license and authors added (Initial coding DB - 11/18/2026); two
  references to DEMO3 restated; paragraph; step comments.
- DEMO12_2DControl.m (1): name line added; log (incl. AT 10/28 entry) and header log line removed; "see
  pie2copvar" restated; markers and (was) lines removed; paragraph; step comments.
- DEMO13_Hinf_OFcontrol.m (1): 6 log entries and ~150 markers/(was) lines removed; pointers to
  PIETOOLS_Hinf_OFcontrol_sop (not called), README.md, EXEC_TOOLS_SOP, lpi_classify_sop, tl.toop/opvar and
  "lines 138 and 141" restated as what the code does; header describes the inlined synthesis and the
  commented hinfsyn block; step comments.
- DEMO14_trace_inequality.m (1): DJ log entry and marker removed; "modeled on DEMO3" dropped; H1op -> H0op
  in the header formula; paragraph; step comments.
- sopvar_demos/heatNd_apply.m (1): header pointers to sopvar.m, the sopvar document and CLAUDE.md removed,
  and the test rationale; license added; log entry, markers, (was) lines removed; step comments; header
  (what/INPUT/OUTPUT) for the local function monvec.

Manual snippets (PIETOOLS_demos/snippets_from_manual/, all new = 0)
- GettingStarted.m: pointer to the .pdf removed; license and authors added (Initial coding DJ - 12/15/2022).
- PIETOOLS_Code_Illustrations_Ch11_Demos.m: manual pointers replaced by a description of each block.
- ..._Ch2_Introduction.m: manual pointer removed; header BCs and -b*phi2 matched to the code; code comments
  phi2_{ss}, phi1_{ss} -> first derivatives (the code uses diff(phi,s,1)); "t = 15" -> "t = 12"; paragraph.
- ..._Ch3_PI_Operators.m: manual pointer removed; "opvar object" -> "copvar object" (3x, the code declares
  copvar); domain of C corrected to q=(x1,x2,v) in R^2 x L2[0,3]; paragraph.
- ..._Ch4_Declaring_PDEs.m: manual pointer removed; each block now states the system it declares.
- ..._Ch5_PIEs.m: manual pointer removed; paragraph; comments on the piess calls (PIE3 as the code builds it).
- ..._Ch6_PIESIM.m: manual pointer removed; each block's system and exact solution stated; step comments.
- ..._Ch7_LPI_Programming.m: manual pointer removed; descriptions of the retired poslpivar, lpivar,
  lpi_ineq, lpi_eq, getObserver, getController renamed to the functions the code calls (poscopvar,
  lpivar_cdopvar, lpi_ineq_sop, lpi_eq_sop, getObserver_sop, getController_sop, also in the commented
  getController line); paragraph.
- ..._Ch8_Declaring_PDEs.m: manual pointer removed; a comment per block on what it builds and prints.

Examples (PIETOOLS_examples/)
- Examples_Library/Stability_Examples/PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m (0): what-it-does line
  (incl. the base-workspace variables it sets); license; authors (DJ - 08/20/2022); step comment.
- command_line_input_Examples/*.m (33 scripts; new = 0: example_2D, l2gain_1-4, stability_1-14 except 13_2,
  13_3; new = 1: l2gain_23, l2gain_26, stability_13_2, 13_3, 16_1, 16_3, 16_4, 17, 18, 19, 36, 39, 42):
  header (system as declared, defaults, what is displayed, workspace PDE/PIE, the optional workspace
  variables params and GUI); license; authors from git (SS 2022, DJ 01/07/2025, DB 05/31 or 06/19/2026).
  Also: stability_13_2/13_3 – stale Timoshenko 4th-order comment block and commented-out code removed,
  inline comments matched to the code (x, not x2; -c*x1; header damping term +c*bet*w_t);
  stability_16_1 – four wrong inline BC comments matched to eq_BC; stability_16_4 – BC comment without gam
  as in the code.
- examples_FI_library_PIETOOLS.m (1): license added; pointer to functional_inequalities/
  heisenberg_ehrling_sop.m, two DEMO2 references, "based on DEMO3, DEMO13" and "the reference figure"
  restated; no-input prompt documented; headers for solve_min, degrees, process_inputs; step comments.

## Skipped
None (all 57 files clean in git status at edit time).

## Kept by decision (for the maintainer)
- "See also" line pointing to the manual: example_2D.m ("See also Chapter 4 and Chapter 8 of the manual").
- Literature citations: DEMO3 (arXiv 1910.01338), DEMO4 (Ahmadi 2015, Valmorbida et al. bibtex), DEMO10
  (Jagt et al., CDC 2025), DEMO13 (formulation of Jungbae Chun), DEMO14 (Barraza Martinez et al.,
  arXiv:2108.13926), Das example (bibtex), FI library [1]-[3]; command-line examples cite "Lamare 2016 [1]"
  ... "Tang 2011 [13]" with no reference list in the files.

## Candidates for the maintainer (not fixed)
1. stability_11.m: `R = 2,9;` sets R = 2 (and prints it), not 2.9.
2. l2gain_2.m: the code declares x_t = -x_ss + s*w (backward heat equation); the comments say x_ss.
3. stability_2-14, 13_2, 16_1, 16_3, 17, 19, 36, 39 run `clc; clear;` before
   `if exist('params','var')`, so a caller's params is deleted and the override never applies (only
   stability_1, 13_3, 18 have the clear commented out).
4. stability_7.m: a params override of lam does not update c (computed on the same line, before).
5. stability_16_1.m: last BC entry has no `==0`; eq_BC prints (no semicolon); the header block's BC list
   does not match eq_BC. stability_13.m (eq_PDE) and stability_17.m (eq_BC) also print.
6. stability_10.m: header comment x(0) = xo, code x(0) = -xo (inline comment matches the code).
7. stability_8.m and PIETOOLS_PDE_Example_..._Das.m use +Bxr*x_s(a); the uncommitted 11/19 log of
   PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m finds this unstable at the default lam = pi^2-1 and switches to -Bxr.
8. PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m: unreferenced older duplicate of the _Ex_ file; loads a
   nonexistent ..._Example_..._GUI.mat; its term format couples x(0) (loc = 0, no D; x(0) = 0 by the BC)
   where its batch format couples x_s(a).
9. DEMO3: `st = settings_sopvar_heavy;` is never used; commented-out duplicate of the BVP declaration.
10. DEMO5: `opts.method='qp'` has no semicolon; the same opts is then passed to PIESIM.
11. DEMO12: the printed "||Z - KP||" is params{1} of the first block of err*err', not a norm.
12. DEMO13: prints the label 'Hinf_OFcontrol_demo12' (file is DEMO13 since the 11/18 renumbering).
13. Ch5 snippet: PIE3 passes B1 = s(2-s) as Bu and 0 as Bw, while the 5.4 PDE has s(2-s) multiplying w.
14. Ch7 snippet: the manual chapter presumably still names the retired functions; the license block asks
    code and manual to match.
15. Ch2 snippet: comment "L2-gain bound ... from 5.2 to 0.78" not re-checked with the _sop executives.
16. DEMO9: license title line reads "PIETOOLS - DEMO8" (license kept as is).
17. new = 0 command-line examples keep blocks of commented-out code in the retired state()/sys() format.

## Unresolved
1. DEMO12 header: "v = d^2/(ds1 ds2) x is the fundamental state". For x twice differentiable in each
   variable the Ch7 snippet uses v = x_{s1s1s2s2}; not verified, left as written.
2. DEMO13 comment "U := (Q - S)^{-1}, such that V = T U" against the imposed equality T'*T*U = V; left.
3. DEMO11 uses c for both the Poincare constant and the lower bound of s2; DEMO14 header alternates xa
   and w for the same input. Left (other authors, no rule violated).
