# dC ledger – report for Declan, 11/21/2026

Agent dC. Items: C37, C41, C42, C43, C45, C47 (converters); C53 part (a), C58 (PIESIM); C101
(copvar/cdopvar messages). HEAD 905b0d26 (fa0deda4 + 905b0d26, which touches only three executives;
none of the files below changed after 78641839 except PIESIM_input_check.m, where the C53(a) line moved
from 150 to 152). No repository file edited.

Scratch (scripts and outputs, rerunnable unchanged):
`C:\Users\mpeet\AppData\Local\Temp\claude\363c499c6332a503\874c85cb-c18d-4876-bea3-39cdd053e21f\scratchpad\dC\`
(below `<S>`): C37_repro.m, C41_repro.m, C43_repro.m (with the helper dC_hit.m), C53a_repro.m, C58_repro.m,
C101_repro.m; each writes `<item>_repro.txt`; the console of each run is in `<item>_console.log`. Every
function under test resolved in the working tree (`which` lines at the top of each .txt).

| Item | Status | Class |
|---|---|---|
| C37 | CONFIRMED (run) | error on a reachable path |
| C41 | CONFIRMED (run) for the DJ files; NOT DJ for the DB and SS lines | error on a rare path |
| C42 | CONFIRMED (reading) | error / wrong selection / misleading display, unreachable at HEAD |
| C43 | CONFIRMED (run, direct call of the helper) | error on a rare or unreachable path |
| C45 | CONFIRMED (reading) | no effect |
| C47 | CONFIRMED (reading) | misleading text |
| C53(a) | CONFIRMED (run) | error on a reachable path |
| C58 | CONFIRMED (run); one more DJ defect in the same file (run) | error on a rare path; misleading legend |
| C101 | CONFIRMED (run); "warning at each access" NOT REPRODUCED | misleading message |

## C37 – CONFIRMED (run)

- HEAD: `converters/convert_PIETOOLS_DDF.m:101, 119, 131`; `converters/legacy/convert_PIETOOLS_DDE2PIE_legacy.m:108,
  145-159` (cited 145-148). Documentation of the kernels as functions of pvar s: `PIETOOLS_DDE.m:159-178`,
  `PIETOOLS_DDF.m:144, 421-434`, and the manual (PIETOOLS_2025_User_Manual.pdf, Tables 4.1, 9.1, 9.2, 9.3:
  "Distributed Delay Terms: May be functions of pvar s"). The help of convert_PIETOOLS_DDF (comment pass
  11/19) says kernels in s1.
- Author: DJ, 8e178278, 2024-12-07 ("Update system default variables"): `pvar s theta` -> `pvar s1 s1_dum`,
  and subs/int in s -> s1, in both files (`git blame -w -C -C -M`, `git show 8e178278`). The 2025 release
  (origin/PIETOOLS-2025) ships the s1 code; its manual says s.
- Reproduction `<S>/C37_repro.m` (measured, C37_repro.txt): x' = -x + int_{-1}^0 k(s) x(t+s) ds, tau = 1.
  k = 0.5 or k = s1: convert_PIETOOLS_DDE and convert_PIETOOLS_DDE2PIE_legacy return a pie_struct.
  k = s (pvar s): convert_PIETOOLS_DDE errors "Input contains independent variables."
  [double:24 < convert_PIETOOLS_DDF:119 < convert_PIETOOLS_DDE:92]; the legacy route errors with the same
  text [double:24 < pi4_copvar_sop:99 < convert_PIETOOLS_DDE2PIE_legacy:181] (Ab0 of line 148 keeps s,
  and the R->R block is converted with double). NDS with Adi{1} = s: same error
  [convert_PIETOOLS_DDF:119 < convert_PIETOOLS_NDS:115]; DDF input with Cvd{1} = s: same error at
  convert_PIETOOLS_DDF:119; with s1 both run.
- Consequence: error on a reachable path – any DDE, NDS or DDF whose distributed-delay kernel depends on
  pvar s, as the manual and the two scripts document; constant kernels and kernels in s1 convert.
- Report entry:

**C37** `converters/convert_PIETOOLS_DDF.m:101, 119, 131`, `converters/legacy/convert_PIETOOLS_DDE2PIE_legacy.m:108, 145-159` (DJ, 8e178278, 2024-12-07) – the kernels are substituted and integrated in s1, while the manual (Tables 4.1, 9.1-9.3), `PIETOOLS_DDE.m:159-178` and `PIETOOLS_DDF.m:144, 421-434` declare them as functions of pvar s. A kernel in s is not integrated: convert_PIETOOLS_DDE stops with "Input contains independent variables." at convert_PIETOOLS_DDF:119 (double of a polynomial in s), and the legacy route of PIETOOLS_DDE.m (DDE_minimal_rep = 0) stops with the same message in pi4_copvar_sop:99, since Ab0 (line 148) keeps s. Error on a reachable path: every DDE, NDS or DDF with a distributed-delay kernel in s, as documented; constant kernels and kernels in s1 convert. Fix: rename the kernel variable s to s1 before the integration, or document s1 in the manual and the two scripts.
Reproduction:
```
pvar s s1
DDE = struct('A0',-1,'tau',1);   DDE.Adi{1} = s;
PIE = convert_PIETOOLS_DDE(DDE,'pie');        % error at convert_PIETOOLS_DDF:119
PIE = convert_PIETOOLS_DDE2PIE_legacy(DDE);   % error at pi4_copvar_sop:99
DDE.Adi{1} = s1;   PIE = convert_PIETOOLS_DDE(DDE,'pie');   % runs
```

## C41 – CONFIRMED (run) for the DJ files; NOT DJ for the DB and SS lines

- HEAD (lines as cited, unchanged since 78641839): `converters/legacy/convert_PIETOOLS_PDE_2D.m:583-587, 612-616`
  (opvar/opvar2d), `:927` (opvar); `convert_PIETOOLS_PDE_nonlinear_2D.m:438-442, 467-471, 844`;
  `legacy/2D/Compute_Tmap_PDE_2D_Separable.m:545, 1230`; `legacy/2D/PIETOOLS_FTC_Expansion_2D.m:229-233, 246-247`;
  `legacy/2D/compute_PIETOOLS_PDE_periodic_map_2D.m:149, 199, 200, 211, 214`;
  `converters/internal/compute_PIETOOLS_PDE_periodic_map.m:190, 201, 208`. Also called and absent:
  clean_opvar (convert_PIETOOLS_PDE_2D.m:572 and others), mat2opvar (727, 757) – no such file at HEAD.
- Authors (`git blame -w -C -C -M`): DJ – 675e5b8e (2022-08-10, "Update 2D PDE to PIE converter"),
  d8918672 (2024-09-07), 0214a3fd (2024-10-16, "Combine 1D and 2D PIE conversion"), 77a4d086
  (2026-02-21); the nonlinear file holds a copy of these lines by DJ (f61e020e, 2026-02-24). DJ moved the old
  converter to legacy: b4256de6 (2026-11-02, "Create ND converter", internal/legacy) and a222bbc7
  (2026-11-16, converters/legacy). NOT DJ: convert_PIETOOLS_PDE_2D.m:927 (DB; marker "DB, 10/27/2026" in
  the file until the comment pass 98790f2e removed it), compute_PIETOOLS_PDE_periodic_map_2D.m (DB,
  e6cd8114, 2026-10-29), compute_PIETOOLS_PDE_periodic_map.m (DB, 82172a74, 2026-10-28); the parser/legacy
  lines are SS (@state/state2terms.m:39 9b700eff 2022-07-24, @state/int.m:89 ab4d569f 2022-10-30,
  @terms/terms.m:38 cd6a8315 2022-02-10), under the C15/C16 decision.
  The classes were retired on 11/08/2026.
- Reachability (grep over the tree without cleaning/ and testing_codes/, and the GUI document.xml): no
  function calls convert_PIETOOLS_PDE_2D or convert_PIETOOLS_PDE_nonlinear_2D except the legacy files
  themselves; the only other caller is testing_codes/private/PIESOS/Testfolder/Test_convert_PDE_nonlinear_varnames.m
  (off the path).
- Reproduction `<S>/C41_repro.m` (measured, C41_repro.txt): `which opvar`, `which opvar2d` empty (exist 0).
  convert_PIETOOLS_PDE_2D(PDE,'silent'): (a) 1-D heat: "Undefined function 'opvar' for input arguments of
  type 'char'." [PIETOOLS_FTC_Expansion_2D:229 < convert_PIETOOLS_PDE_2D:416]; (b) 2-D heat: "Undefined
  function 'opvar2d' ..." [PIETOOLS_FTC_Expansion_2D:246]; (c) 2-D heat with the boundary input w(s1):
  "The expansion of the PDE state in terms of the PIE state fails; perhaps the PDE is ill-posed, or not
  currently supported." [Compute_Tmap_PDE_2D_Separable>construct_1D_Tops:524 < :214 <
  convert_PIETOOLS_PDE_2D:380] – the catch at lines 522-524 replaces the missing-class error; (d) 1-D
  nonlinear (x_t = x_ss - x.*x): 'opvar' [PIETOOLS_FTC_Expansion_2D:229 < convert_PIETOOLS_PDE_nonlinear_2D:353
  < convert_PIETOOLS_PDE_2D:232]; convert_PIETOOLS_PDE_nonlinear_2D directly: the same. convert_PIETOOLS_PDE
  converts (a), (b), (c) (copvar T).
- Consequence: error on a rare path – a direct call of either legacy converter with a pde_struct; for a
  2-D PDE with spatially varying boundary inputs the message blames the PDE.
- Report entry:

**C41** `converters/legacy/convert_PIETOOLS_PDE_2D.m:583-616`, `convert_PIETOOLS_PDE_nonlinear_2D.m:438-471, 844`, `legacy/2D/Compute_Tmap_PDE_2D_Separable.m:545, 1230`, `legacy/2D/PIETOOLS_FTC_Expansion_2D.m:229-247` (DJ, 675e5b8e 2022-08-10, d8918672 2024-09-07, 0214a3fd 2024-10-16, 77a4d086 2026-02-21) – the files construct objects of the classes opvar and opvar2d, retired on 11/08/2026, and call clean_opvar and mat2opvar, which no longer exist. A pde_struct passed to convert_PIETOOLS_PDE_2D or convert_PIETOOLS_PDE_nonlinear_2D reaches PIETOOLS_FTC_Expansion_2D or Compute_Tmap_PDE_2D_Separable (lines 377-416) and stops: measured, "Undefined function 'opvar'" ('opvar2d' in 2-D) at PIETOOLS_FTC_Expansion_2D:229 (246) for 1-D and 2-D heat equations and a nonlinear 1-D PDE; a 2-D PDE with a spatially varying boundary input stops instead with "The expansion of the PDE state in terms of the PIE state fails; perhaps the PDE is ill-posed, or not currently supported." (Compute_Tmap_PDE_2D_Separable:522-524 catches the missing-class error). Error on a rare path: no function calls these converters except each other (moved to legacy in b4256de6, 11/02/2026, and a222bbc7, 11/16/2026); convert_PIETOOLS_PDE converts the same PDEs. The opvar lines of DB (convert_PIETOOLS_PDE_2D.m:927, compute_PIETOOLS_PDE_periodic_map_2D.m, internal/compute_PIETOOLS_PDE_periodic_map.m) belong to the same path.
Reproduction:
```
pvar t s
x = pde_var('state',1,s,[0,1]);
PDE = initialize([diff(x,t)==diff(x,s,2); subs(x,s,0)==0; subs(x,s,1)==0]);
PIE = convert_PIETOOLS_PDE_2D(PDE,'silent');   % Undefined function 'opvar' [PIETOOLS_FTC_Expansion_2D:229]
PIE = convert_PIETOOLS_PDE(PDE);               % runs
```

## C42 – CONFIRMED (reading)

- HEAD `converters/legacy/2D/Compute_Tmap_PDE_2D_Separable.m`: local compute_input_diffs lines 888-893
  (Dop_x, Dop_y: one entry per state component, then repelem to one entry per scalar state variable),
  960-961 and 1015-1016 (`xdiffs_k = unique(Dop_x(tmp_log))`, tmp_log 4 x 1: one flag per row block R,
  L2[x], L2[y], L2[x,y], lines 914-956); local print_reorder_diff_summary lines 1409
  (`varnames_RHS_t = varname_LHS_t;` – no line of the function assigns varname_LHS_t), 1444 (`Ldiff = [...]`
  assigns inside the loop over the variables, so a derivative in s1 and s2 keeps only the last), 1453
  (`ceil(log(dval(jj)/log(10)))`: the division is inside the log).
- Authors: DJ – d8918672 (2024-09-07, "Add support for boundary forcing") for 892-893, 961, 1016, 1444,
  1453 (`git blame -w -C -C -M`, `git log -S` for 1444 and 1453); d7fd1c0e (2026-02-15, "Update PDE to
  PIE converter") introduced line 1409 (`git log -S'varname_LHS_t'`). The comments at 910-911 are MMP
  (comment pass).
- Reachability (measured in C41 case (c)): Compute_Tmap_PDE_2D_Separable stops at line 214
  (construct_1D_Tops, the opvar error of C41) before compute_input_diffs (lines 262, 271) and
  print_reorder_diff_summary (309-355). Unreachable at HEAD. Inferred from convert_PIETOOLS_PDE_2D:377-380:
  before the retirement of the classes (11/08/2026) 2-D PDEs with spatially varying boundary inputs
  reached them (line 1409 only since d7fd1c0e, 2026-02-15, after the 2025 release).
- Semantics: 1409 – "Unrecognized function or variable 'varname_LHS_t'" for any component without spatial
  variables that is not a boundary value; 961/1016 – a logical index shorter than the array selects
  entries 1-4 of Dop_x (state variables 1-4), not the states w(k) maps to, and errors if a true flag lies
  beyond numel(Dop_x); 1444, 1453 – the printed summary only.
- Report entry:

**C42** `converters/legacy/2D/Compute_Tmap_PDE_2D_Separable.m:961, 1016, 1409, 1444, 1453` (DJ, d8918672 2024-09-07; line 1409 d7fd1c0e 2026-02-15) – compute_input_diffs indexes Dop_x, Dop_y (one entry per scalar state variable after the repelem of lines 892-893) with tmp_log (4 x 1, one flag per row block R, L2[x], L2[y], L2[x,y]), so it reads the derivative orders of state variables 1-4, not of the states that w(k) enters; print_reorder_diff_summary reads varname_LHS_t (line 1409; varnames_LHS_t elsewhere), keeps only the last variable in Ldiff (line 1444 assigns instead of appending), and computes ceil(log(dval(jj)/log(10))) (line 1453). Unreachable at HEAD: the file stops at line 214 first (C41); before the retirement of the classes, 2-D PDEs with spatially varying boundary inputs reached these lines (inferred). Fix (1409, 1444, 1453): varnames_LHS_t; `Ldiff = [Ldiff,'(',...]`; log(dval(jj))/log(10).
Reproduction: by reading – line 1409 reads a variable that no line of the function assigns; a logical index (961, 1016) selects by position, so a 4 x 1 mask on a vector of sum(size_x) entries selects entries 1-4.

## C43 – CONFIRMED (run, direct call of the helper)

- HEAD `converters/legacy/2D/expand_PIETOOLS_PDE_BCs.m:545`
  `Rindx = find(PDE.([Robj,'_tab'])(:,1)==BC_new_ii{kk}.term{jj}.x);` – Robj is set only in the term loop
  of lines 206-220 ('x', 'w' or 'u' of the last term of BC ii); line 547 then reads PDE.x_tab(Rindx(1),:).
- Author: DJ, 675e5b8e (2022-08-10, "Update 2D PDE to PIE converter") for 220 and 545 (`git blame -w -C -C -M`).
- Reproduction `<S>/C43_repro.m` (measured, C43_repro.txt): the callers stop earlier (C41), so
  expand_PIETOOLS_PDE_BCs is called directly; a conditional breakpoint whose condition logs and returns
  false (`<S>/dC_hit.m`) records each pass through line 545, with line 220 as control (44-85 passes per
  case). States x1, x2 on [0,1]^2 (x2 index 2), one input w (index 1). BCs x2(s1,0) + w = 0 (x2 of order
  [2;2] or [1;2]): line 545 never reached. BC int_0^s1 x2(th,0) dth + w = 0, x2 of order [0;2], input
  term last: line 545 reached with Robj = 'w', term x = 2; "Index exceeds array bounds."
  [expand_PIETOOLS_PDE_BCs:547]. Same BC with the input term first: line 545 with Robj = 'x', runs
  (11 BCs). Inferred: with the state and the input both of index 1 the lookup returns the right row from
  the wrong table.
- Consequence: error on a rare or unreachable path – only a direct call of this helper reaches it at
  HEAD; the outcome depends on the order of the terms of a BC.
- Report entry:

**C43** `converters/legacy/2D/expand_PIETOOLS_PDE_BCs.m:545` (DJ, 675e5b8e, 2022-08-10) – the integral-removal loop looks up a state term in PDE.([Robj,'_tab']), Robj being left by the term loop of lines 206-220 at the object of the last term of the BC. Measured with a direct call: for the BC int_0^s1 x2(th,0) dth + w = 0 (x2 the second state, w the only input) line 545 runs with Robj = 'w', finds no row, and line 547 stops with "Index exceeds array bounds."; the same BC written w + int(...) = 0 expands. Error on a rare or unreachable path: the two callers (convert_PIETOOLS_PDE_2D:432, convert_PIETOOLS_PDE_nonlinear_2D:368) stop earlier at HEAD (C41). Fix: PDE.x_tab at line 545.
Reproduction:
```
pvar t s1 s2
x1 = pde_var('state',1,[s1;s2],[0,1;0,1]);  x2 = pde_var('state',1,[s1;s2],[0,1;0,1],[0;2]);  w = pde_var('input',1);
PDE = initialize([diff(x1,t)==diff(x1,s1,2)+diff(x1,s2,2); subs(x1,s1,0)==0; subs(x1,s1,1)==0; subs(x1,s2,0)==0; subs(x1,s2,1)==0;
                  diff(x2,t)==diff(x2,s2,2); int(subs(x2,s2,0),s1,0,s1)+w==0; subs(x2,s2,1)==0]);
PDEx = expand_PIETOOLS_PDE_BCs(PDE);   % Index exceeds array bounds. [expand_PIETOOLS_PDE_BCs:547]
```

## C45 – CONFIRMED (reading)

- HEAD: local `print_convert_summary(PIE,obj,tdiff_list)` at `converters/legacy/convert_PIETOOLS_PDE_2D.m:1942`
  (called at 905) and `convert_PIETOOLS_PDE_nonlinear_2D.m:1649` (called at 676); the file
  `converters/print_convert_summary.m` defines print_convert_summary(PIE,obj,aux_idcs,F_aux), called by
  convert_PIETOOLS_PDE.m:625, 629.
- Authors: the local functions DJ (a7ad0239, 2025-06-18, in convert_PIETOOLS_PDE.m; carried into the
  legacy files with the copies of that converter, b4256de6 for convert_PIETOOLS_PDE_2D). The same-named
  file was created by DB on 11/18/2026 (7764389f, "Move print_convert_summary to its own file"); its
  signature comes from convert_PIETOOLS_PDE.m, changed by DB 87236133 (11/07/2026, 6 arguments) and DJ
  dbf458a2 (11/13/2026, 4 arguments).
- Semantics: a call inside a file resolves to a local function before any function on the path, so 905
  and 676 call the local 3-argument function; both calls lie after the opvar errors of C41.
- Report entry:

**C45** `converters/legacy/convert_PIETOOLS_PDE_2D.m:1942`, `converters/legacy/convert_PIETOOLS_PDE_nonlinear_2D.m:1649` (DJ, a7ad0239, 2025-06-18) – each file keeps a local print_convert_summary(PIE,obj,tdiff_list), while `converters/print_convert_summary.m` (moved out of convert_PIETOOLS_PDE.m by DB on 11/18/2026) defines the same name with the arguments (PIE,obj,aux_idcs,F_aux). No effect: a local function takes precedence for calls in its file (lines 905, 676), and both calls are unreachable at HEAD (C41).
Reproduction: by reading – MATLAB resolves a call inside a file to its local function before the path.

## C47 – CONFIRMED (reading)

- HEAD `converters/operator_converters/dpvars2sdvars.m:64` "% PIETOOLS - quad_expand" (the licence title of
  PIESOS/internal/quad_expand.m). `help dpvars2sdvars` prints lines 2-61 only (blank line 62), so the text is
  seen in the source.
- Author: DJ, 6d39aee9, 2026-10-05 ("Create dopvar to sdopvar converters"; then in sopvar/Testfolder/converters,
  moved by MMP in 43c97ba7).
- Report entry:

**C47** `converters/operator_converters/dpvars2sdvars.m:64` (DJ, 6d39aee9, 2026-10-05) – the licence title line reads "PIETOOLS - quad_expand", copied from PIESOS/internal/quad_expand.m. Misleading text in the source (help prints lines 2-61). Fix: "PIETOOLS - dpvars2sdvars".
Reproduction: by reading – line 64.

## C53 (a) – CONFIRMED (run)

- HEAD `PIESIM/PIESIM_routines/PIESIM_Setup/PIESIM_input_check.m:152` (150 at 78641839)
  `PDE = initialize_PIETOOLS_PDE(psize,true);` – psize is `struct()` with the field dim (lines 135-136);
  the block (lines 140-155) adds a temporary 1-D state to PDE and then initializes psize instead.
- Author: DJ, d6d2d88e, 2022-10-10 ("Update PIESIM new terms format.", the block copied from the 2-D
  converter, where the call is initialize_PIETOOLS_PDE(PDE,true)); `git log -S'initialize_PIETOOLS_PDE(psize,true)' --all`.
  Plain blame names YP (web re-uploads, ed2fa4fc). Correction to D_piesim_ledger.md: it cites d6e76870,
  which is the merge of PR #34 (2024-12-16); the line comes from d6d2d88e.
- Reproduction `<S>/C53a_repro.m` (measured, C53a_repro.txt): x' = -x, y' = -2y as a pde_struct (dim 0),
  PIESIM with N 8, tf 0.5, dt 1e-3, BDF, plot 'no': "The input PDE is not appropriately specified. Please
  define your PDE as a "pde_struct" class object, ..." [initialize_PIETOOLS_PDE:75 < PIESIM_input_check:152
  < PIESIM:194]. With PDE in place of psize (scratch copy `<S>/c53_shadow/PIESIM_input_check.m`, first on
  the path): the input check and the conversion pass, then "Unrecognized function or variable
  'acheb_glob'." [PIESIM_discretize_icf:104 < PIESIM_discretize_all:66 < PIESIM:253] – acheb_glob is
  assigned only in the loop over the nx PDE states (lines 95-99), nx = 0; that code is of the 2021b
  import (YP), not DJ.
- Consequence: error on a reachable path – PIESIM on a pde_struct without spatial variables, a case the
  block of lines 140-155 exists to handle.
- Report entry:

**C53(a)** `PIESIM/PIESIM_routines/PIESIM_Setup/PIESIM_input_check.m:152` (DJ, d6d2d88e, 2022-10-10) – for a pde_struct without spatial variables, the block of lines 140-155 adds a temporary 1-D state to PDE and then calls initialize_PIETOOLS_PDE(psize,true), psize being struct() with the field dim; PIESIM stops with "The input PDE is not appropriately specified. ..." (initialize_PIETOOLS_PDE:75). With PDE passed instead (scratch copy), the run stops later at PIESIM_discretize_icf:104 ("Unrecognized function or variable 'acheb_glob'": assigned only in the loop over PDE states; 2021b code), so ODE-only input needs more than this change. Error on a reachable path: PIESIM on an ODE-only pde_struct. Fix: PDE for psize at line 152 (necessary, not sufficient).
Reproduction:
```
X = pde_var(2);   PDE = initialize(diff(X,'t')==[-1 0;0 -2]*X);
opts = struct('N',8,'tf',0.5,'dt',1e-3,'intScheme',1,'plot','no');
sol = PIESIM(PDE,opts,struct('ic',[1 1]));   % error [initialize_PIETOOLS_PDE:75 < PIESIM_input_check:152]
```

## C58 – CONFIRMED (run)

- HEAD `PIESIM/PIESIM_routines/Postprocessing/PIESIM_plot_solution.m:229` (colour index
  `colors_y(psize.no0-ii+1,:)`, 7 rows, and label `$z_`), `:148-150` (return before `figs = {}` at 166);
  `PIESIM/PIESIM_routines/2D/Postprocessing/PIESIM_plot_solution_2D.m:220` (same index and label).
- Authors (`git log -S` over all refs; plain blame names YP, web re-uploads): DJ d00ccc83 (2024-12-26,
  "Change PIESIM plot_solution": the observed-output plot with `colors_y(psize.noo-ii,:)` and `$z_`, and
  `function figs = ...`), 770742f5 (2024-12-28, "Bugfix": index noo-ii+1), fe5c062a (2025-01-01, "Update
  PIESIM plot functions": the early return, `figs = {}`, the 2-D file); YP a8d9e617 renamed noo -> no0.
  D_piesim_ledger.md names fe5c062a for the colour index; it predates it (d00ccc83, 770742f5).
- Reproduction `<S>/C58_repro.m` (measured, C58_repro.txt): (a) `figs = PIESIM_plot_solution(...)` with
  plot 'no': 'Output argument "figs" (and possibly others) not assigned ...'; (b) direct calls with
  no0 = 8, BDF, plot 'yes': "Index in position 1 exceeds array bounds. Index must not exceed 7."
  [PIESIM_plot_solution:229], [PIESIM_plot_solution_2D:220]; (c) no0 = 2: runs, figure "Observed Output
  Evolution", legend labels $z_1(t)$ $z_2(t)$ (1-D and 2-D); (d) PIESIM on the heat equation with 8 sensed
  outputs y = (1:8)'*int x, plot 'yes': the error of (b) [PIESIM_plot_solution:229 < PIESIM:285].
  Found in passing (DJ, d00ccc83): (e) PIESIM on 8 decoupled heat equations, plot 'yes': "Index in position
  1 exceeds array bounds. Index must not exceed 7." [PIESIM_plot_solution:314 < PIESIM:285] –
  `colors_PDE(n+psize.n0,:)` (lines 314, 316) has 7 rows and no guard (the 2-D file guards the same index,
  lines 319, 364); plot 'no' runs. By reading (DJ, fe5c062a): line 359 `num2str(ii)+psize.n0` adds n0 to a
  character code (wrong subplot titles for ii+n0 >= 10), masked by the line 314 error.
- Consequence: error on a rare path (more than 7 finite-dimensional observed outputs, or more than 7 ODE
  plus PDE states, with opts.plot = 'yes'; the default is 'no'); misleading legend for any number of
  observed outputs; figs unassigned only for a direct call with an output (PIESIM.m:285 requests none).
- Report entry:

**C58** `PIESIM/PIESIM_routines/Postprocessing/PIESIM_plot_solution.m:229, 148-150`, `PIESIM/PIESIM_routines/2D/Postprocessing/PIESIM_plot_solution_2D.m:220` (DJ, d00ccc83 2024-12-26, 770742f5 2024-12-28, fe5c062a 2025-01-01) – the finite-dimensional observed outputs 1-7 are coloured colors_y(psize.no0-ii+1,:) from 7 rows: with 8 outputs and opts.plot = 'yes', PIESIM stops with "Index in position 1 exceeds array bounds. Index must not exceed 7." (line 229 via PIESIM.m:285; 2-D line 220); with fewer, the legend labels them $z_1(t)$, $z_2(t)$, the regulated-output symbol, in the figure "Observed Output Evolution". The 1-D file returns at line 149, before figs = {} (line 166), when opts.plot is not 'yes', so a direct call with an output errors (PIESIM requests none). Also: colors_PDE(n+psize.n0,:) (lines 314, 316, d00ccc83) has 7 rows and no guard – 8 PDE states with opts.plot = 'yes' stop at line 314 (the 2-D file guards this index); line 359 adds psize.n0 to the character code of num2str(ii) (wrong titles from ii+n0 = 10; masked by line 314). Error on a rare path (opts.plot = 'yes'; default 'no'); misleading legend. Fix: colors_y(ii,:) or a modulo index; $y_; figs = {} before the return; the guard of the 2-D file at 314, 316; num2str(ii+psize.n0) at 359.
Reproduction:
```
pvar t s;   syms sx
x = pde_var('state',1,s,[0,1]);   y = pde_var('sense',8);
PDE = initialize([diff(x,t)==diff(x,s,2); y==(1:8)'*int(x,s,[0,1]); subs(x,s,0)==0; subs(x,s,1)==0]);
opts = struct('N',6,'tf',0.01,'dt',1e-3,'intScheme',1,'plot','yes');
sol = PIESIM(PDE,opts,struct('ic',sin(pi*sx)));   % error [PIESIM_plot_solution:229 < PIESIM:285]
```

## C101 – CONFIRMED (run); sub-claim "warning at each access" NOT REPRODUCED

- HEAD: copvar `sopvar/@copvar/subsref.m:104-111, 117-118`, `subsasgn.m:113-120, 151-152`,
  `getfield_opvar1d.m:79, 132`, `getfield_opvar2d.m:89, 124`, `setfield_opvar1d.m:72, 98, 104`,
  `setfield_opvar2d.m:96, 177`; cdopvar `sopvar/@cdopvar/subsref.m:104, 107, 110, 117`,
  `subsasgn.m:115, 118, 121, 153`, `getfield_opvar1d.m:85, 138`, `getfield_opvar2d.m:94, 129`,
  `setfield_opvar1d.m:72, 98, 104`, `setfield_opvar2d.m:96, 181`.
- Authors (verified): every one of these lines DJ – 3b92e68b (2026-11-15, "Support old opvar fields in
  copvar subsref"; subsref, getfield_*, and the parts of setfield_* copied from getfield_*) and 2f2af6a1
  (2026-11-16, "Update copvar/cdopvar subsasgn"); `git blame -w -C -C -M` on each line, and
  `git log -S'copvar:opvarField' --all` finds only these two commits. The file headers of subsref.m and
  subsasgn.m name MMP (Initial coding 09/29/2026), hence the MMP credit in CANDIDATES.md; the A_sopvar
  evidence holds.
- Reproduction `<S>/C101_repro.m` (measured, C101_repro.txt; T = PIE.T of an ODE coupled to a heat
  equation, T2 = T([2 1],[2 1]) with the L2 block first, Td = cdopvar(T)): each message names 'opvar',
  'opvar2d' (copvar) or 'dopvar', 'dopvar2d' (cdopvar), and "opvar-like"/"opvar2d-like", e.g.
  copvar:notOpvar1d "'P' is a field of an 'opvar', but this copvar is not opvar-like: ...",
  copvar:badValue "An 'opvar' has one spatial variable: 'I' must have one row.", cdopvar:notOpvar2d
  "'R00' is a field of a 'dopvar2d', ...". The warning copvar:opvarField is given at the first read
  (T.P) and at the first assignment (T.P = 2) of the session only: T.Q1, T.R.R0 and T.Q2 = ... give none
  (persistent `warned` in subsref and subsasgn since 3b92e68b, 2f2af6a1); cdopvar likewise.
  Reachability: the subsref/subsasgn messages and setfield_opvar1d:98 (T.I = zeros(0,2)) by indexing;
  the structure and field-name messages of getfield_*/setfield_* only by direct calls of these methods
  (subsref and subsasgn test the structure and the name first).
- Proposed wording of A_sopvar_ledger.md, checked: names no retired class (is_opvar1d, is_opvar2d are
  functions at HEAD); its conditions match the help of is_opvar1d (at most one spatial variable, R blocks
  before L2 blocks) and is_opvar2d (at most two, blocks ordered R, L2[x], L2[y], L2[x,y]); one %s each,
  identifiers unchanged. Gaps: it covers subsref and subsasgn only – the getfield_*/setfield_* messages
  listed above need the same change ("does not have the 1-D form R^m x L2^n[a,b] -> R^p x L2^q[a,b]",
  "is not a field of the 1-D form", "The 1-D form has one spatial variable: ..."; 2-D and cdopvar
  likewise); and it writes the 1-D form as the domain R^m x L2^n[a,b], where getfield_opvar1d:80 and the
  help of is_opvar1d write the map R^m x L2^n[a,b] -> R^p x L2^q[a,b].
- Report entry:

**C101** `sopvar/@copvar/subsref.m:104-118`, `subsasgn.m:113-120, 151-152`, `getfield_opvar1d.m:79, 132`, `getfield_opvar2d.m:89, 124`, `setfield_opvar1d.m:72, 98, 104`, `setfield_opvar2d.m:96, 177`, and the @cdopvar twins (DJ, 3b92e68b 2026-11-15, 2f2af6a1 2026-11-16) – the errors and warnings of the field interface describe the operator through the retired classes: "'P' is a field of an 'opvar', but this copvar is not opvar-like: ...", "The operator does not have the structure of an 'opvar2d', ...", "An 'opvar' has one spatial variable: 'I' must have one row.", the warning copvar:opvarField "... as for an 'opvar' or 'opvar2d'.", and 'dopvar', 'dopvar2d' for a cdopvar. The warning is given once per session (first read, first assignment), not at each access. Misleading message: the classes no longer exist; the user meets the subsref/subsasgn texts by indexing (P.Q1, P.R00, P.I = ...), the getfield/setfield texts by direct calls. Fix: name the form instead, e.g. "'%s' is a field of the 1-D form R^m x L2^n[a,b] -> R^p x L2^q[a,b] (fields P, Q1, Q2, R), but this copvar does not have that form: ...", identifiers unchanged.
Reproduction:
```
pvar t s;   x = pde_var(1);   z = pde_var('state',1,s,[0,1]);
PIE = convert_PIETOOLS_PDE(initialize([diff(x,t)==-x; diff(z,t)==diff(z,s,2); subs(z,s,0)==0; subs(z,s,1)==x]));
T = PIE.T;   v = T.P;   v = T.Q1;     % copvar:opvarField "... as for an 'opvar' or 'opvar2d'", at T.P only
T2 = T([2 1],[2 1]);   v = T2.P;      % copvar:notOpvar1d "'P' is a field of an 'opvar', ..."
Td2 = cdopvar(T2);   v = Td2.R00;     % cdopvar:notOpvar2d "'R00' is a field of a 'dopvar2d', ..."
```
