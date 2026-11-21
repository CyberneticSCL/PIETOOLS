# Release comment pass – candidates for the maintainer

Source: the "Candidates for the maintainer" sections of the 13 ledgers of the release comment pass of 11/19–11/20/2026 (commits ece48c96..f474e7de; `lists/*_ledger.md`, `lists/resume/*_ledger.md`) and the list of user-visible error strings in `lists/B_lpi_exec_settings_ledger.md`. Status: every item was found by reading and re-read against HEAD 78641839 (after lrmultiply moved into `sopvar/@sdopvar/mtimes.m`, weightPIE into `PIESIM/PIESIM_routines/PIESIM.m`, clencurt to `cleaning/private/dead_code`); none was run or reproduced; none is fixed. Paths are relative to `PIETOOLS/`, line numbers those of HEAD. Evidence "semantics": the defect follows from MATLAB language rules at that line, whatever the data; "read": it follows from the control or data flow, not run. Author: the Initial coding line or log of the file, or `git blame` of the cited lines with the pass commits ignored; "2021b" marks code imported with PIETOOLS 2021b (12/31/2021), credited to the authors named in the file header.

## Outcome of the bug-fix pass, 11/20/2026

Maintainer, 11/20/2026: "fix the unfixed bugs". Records: `cleaning/private/bugfix_1120/` (BRIEF.md and the ledgers A_sopvar, B_conv, D_piesim, E_examples). Each fix was reproduced before the edit and verified after it against a reference independent of the repaired code: closed forms, exact solutions, a second conversion route, collocation outside PIETOOLS. Valid inputs are unchanged: 46 delay-system conversions and the 36 GUI files of the example library, 41 1-D and 32 2-D PIESIM examples and 5 DDE or PIE-input cases are bit-identical; test_executives_sop passes 20 of 20 checks (2-D included); the parse check passes 75 of 75 edited files. Annotation: new files (`release_new_files_1120/new_files_classified.tsv`) carry no log line or marker; existing files carry a log line "MMP, 11/20/2026", column-77 markers and (was) lines; help text is corrected in place.

| Outcome | Items |
|---|---|
| Fixed: 52, and 3 new | C1, C3–C6, C8–C12, C14, C17, C31, C32, C35, C36, C38–C40, C46 (text), C48–C52, C55–C57, C59, C61, C74–C77, C79, C84 (beam file), C85–C87, C89, C90, C92 (header of stability_8), C93–C95, C102–C105, C107–C109; new: C112, C113, C116 |
| Partly fixed | C53: (b) fixed; (a) held, a line of Declan (d6e76870) |
| Held: code of Declan | the 33 items credited DJ below; C101 (`git blame -w -C`: every message line is from 3b92e68b, 11/15/2026, or 2f2af6a1, 11/16/2026, both Declan; the file headers name MMP); C58 (fe5c062a); the Saba file of C84; line 62 of stability_8 and the Das files of C92; new: C115 |
| Held: PIESOS (Declan; @quadPoly was ndopvar code) | C63, C64, C66–C71, C73, C111 |
| Held: 2-D settings mapping | C7, C106 |
| Needs decision | C15, C16, with C114: parser/legacy needs the retired opvar class – port or retire the folder; C21–C23, C44 and convert_PIETOOLS_PDE_terms_internal (C46): no caller – move to cleaning/private/dead_code; C36: ODE.A overrides an explicit nx; C60: `T.I = I` in transl_mod was removed by YP on 5/2/2026 ("Fixed rescaling") – restore it (no result changed in 46 runs) or correct the help; C8 (fixed: the 2-D retry uses 'inv', the alternative method): or delete the retry, as the other synthesis executives do |
| No change | C2, C91, C110 (decided 11/20/2026) |
| Not reproduced | C54: subs accepts the 1 x 1 cell (Symbolic Math Toolbox) |

Limit left: the symbolic integrator (opts.intScheme = 2, C52) still accepts no sampled inputs.

New items found in the pass:

**C112** `converters/internal/poly_decomp_1D_econ.m:75` (MMP; 2021b) – `monomials(Qt(:,i))` with one argument, while `SOSTOOLS400/monomials.m` needs two, so every 'polynomial' distributed-delay kernel errored on the minimal-DDF routes (minimize_PIETOOLS_DDE2DDF, minimize_PIETOOLS_DDF); a polynomial without variables (A + 0*s1) errored at Q.Varname(1) – FIXED: the monomials of each column from its degmat, a constant as 'double'; measured: every route matches the closed form to 1.8e-15.

**C113** `converters/initialize_PIETOOLS_DDE.m:296` (MMP; 2021b) – tests exist('Bdi'), not B1di, as in C39 – FIXED: a DDE whose nw is set by B1di alone gets nw = 1 (before: a dimension error).

**C114** `parser/legacy/@sys/sys.m` (SS) – declares the method convert with one output, so the second output documented in convert.m cannot be returned – with C15, C16.

**C115** `EL/Stability_Examples/PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m` (DJ) – with -Bxr (11/19) the coupled system loses stability at lam = 9.7946 (eigenvalues by Chebyshev collocation and by sine-Galerkin, outside PIETOOLS), so "stable for lam < pi^2" holds for neither sign; with +Bxr it is stable only for lam < -5.7799 – report to Declan.

**C116** `CL/l2gain_2.m:40` (SS) – the comment stated gamma = 0.3333; after C75 the system is the heat equation of library case 22, whose H-infinity norm is 5/24 = 0.2083, attained at w = 0 (Chebyshev collocation, independent of PIETOOLS) – FIXED (comment).

## Report for Declan, 11/21/2026

Maintainer, 11/21/2026: "draft the report for Declan". Draft: `cleaning/private/bugfix_1120/DECLAN_REPORT.md` – 35 items in code of Declan and 10 PIESOS items of other authors, each located on fa0deda4, its author traced with `git blame -w -C -C -M` and `git log -S`, and reproduced unless marked "by reading"; ledgers and brief in `cleaning/private/bugfix_1120/declan/`, every reproduction rerun by the coordinator with the same output.

Out of the hold after the author tracing (fixable here): C13 (SS, adf9f3f9; Declan added the header only), C25 (2021b install script), C24 for `PIETOOLS_DDF.m:272` (2021b text), C81 (DB, 9efb8fff; the header "Initial coding DJ" is C120), C98 (DB, 7e551368, both files), C115 (the claim is 2021b text; -Bxr is 99f7dff5, MMP, 11/19/2026; the sign question went into the report as C92).

Corrections to the items: C80 is an error on a reachable path, not a silent wrong result; C68 is a silent wrong result; C81 fails at load for GUI = 1; C53(a) comes from d6d2d88e, not from the merge d6e76870; C58 from d00ccc83, 770742f5 and fe5c062a; C101 warns once per session, not at each access.

**C117** `PIESOS/polyopvar/@tensopvar/initialize.m:62-64` (DJ, 6fd03c15) – converts factors of the retired classes opvar and opvar2d; the branch never runs – in the report.

**C118** `sopvar/@sopvar/inv_1D.m:497` (SS, 0992dfd0), local zeroKernelLike – calls quadPoly() with no input, which errors ("quadPoly must have 1 or 6 inputs, nothing else.", measured), on the branch of zero kernels (not run).

**C119** `PIESOS/polyopvar/@quadPoly/subsasgn_1.m:2-10` (MMP, 8e335783, release comment pass) – the help says that assignment throws an error; the overload is inactive and assignment succeeds (C69).

**C120** `EL/PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion.m`, header line "Initial coding DJ - 08/18/2022" (MMP, 978570e9, release comment pass) – the file is by DB (9efb8fff, 2025-01-06).

**C121** `EL/Stability_Examples/PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das_GUI.mat` (2021b) – stores +Bxr, while the code uses -Bxr since 99f7dff5 (MMP, 11/19/2026): GUI = 1 loads the system that is unstable at the default lam – with C92.

Fixed 11/21/2026 (maintainer: "Fix them"), each reproduced before the edit and checked after it:
- C13 `pietools_path_update.m`: adds the folder of the file. Called by name from a user folder it now adds the 89 PIETOOLS entries that it adds from the PIETOOLS folder or with run, and not the user folder (before: the user folder and its subfolder, no PIETOOLS folder).
- C24, the line of `PIETOOLS_DDF.m`: sosineq_on removed (kept as a (was) line).
- C25 `pietools_update.m`: startup.m closed in both branches (stubbed copy: no handle open afterwards; before: one).
- C81, C120 `EL/PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion.m`: GUI = 1 raises "No GUI file associated with this example", as the _Control example (before: load error); help corrected; Initial coding DB - 01/06/2025.
- C98: the Parabolic example prints no warning (SeDuMi certifies gam = 23.7312 at its st; MOSEK 23.7285); the reaction-diffusion example names the tier that certifies with SeDuMi (light 8.2133, heavy 8.2140, veryheavy 8.2094 certified; its st, tuned for MOSEK: 6.238, not certified; MOSEK at st 8.2156).
- C115, C121: the GUI file stores -Bxr (data.eqns{5-8}, nothing else changed); the PIE of the GUI route equals the PIE of the code to 3.4e-6 (lam rounded to 8.8696 in the file) and differs from the +Bxr PIE by 2.67.
- C118 `sopvar/@sopvar/inv_1D.m:497`: quadPoly(0); kernels with an empty basis reach the branch and give M^{-1} = -1 and zero kernels, as zero kernels give on the general path (before: "quadPoly must have 1 or 6 inputs").
- C119: the help of subsasgn_1.m states that the file is not a method subsasgn and that assignment is unchecked.

The report entry C92 records the decision of 11/19/2026 (the sign of the paper is a typo; -Bxr) instead of asking for one.

## Summary

150 ledger entries (130 candidates, 20 string entries): 8 merged as duplicates, 33 not confirmed on re-reading, 111 items below (3 entries split into several items, 2 grouped with others).

| Area | Silent wrong result | Error, reachable path | Error, rare or unreachable path | Dead code or no effect | Misleading text or message | Items |
|---|---|---|---|---|---|---|
| sopvar | 0 | 0 | 1 | 3 | 1 | 5 |
| LPI programming, executives, settings | 0 | 2 | 0 | 3 | 1 | 6 |
| parser, GUI, root scripts | 2 | 3 | 4 | 10 | 4 | 23 |
| converters | 2 | 1 | 6 | 2 | 2 | 13 |
| PIESIM | 3 | 4 | 5 | 2 | 0 | 14 |
| PIESOS | 1 | 2 | 4 | 3 | 2 | 12 |
| demos and examples | 10 | 0 | 1 | 7 | 9 | 27 |
| strings naming retired functions or files | 0 | 0 | 0 | 0 | 11 | 11 |
| **Total** | **18** | **12** | **21** | **30** | **30** | **111** |

## sopvar

**C1** `sopvar/@sopvar/end.m:43-48` (AT), `sopvar/@sdopvar/end.m:41-46` (MMP) – no branch for n > 2, so `out` stays unassigned – error on a rare path: "Output argument not assigned" for an index of three or more subscripts that uses `end` – semantics – fix: an else branch with an error, as in `sopvar/@copvar/end.m`.

**C2** `sopvar/@copvar/copvar.m:115-122` (MMP) – the protected properties leftCommonBasis, leftCommonC, rightCommonBasis, rightCommonC are written by no method on the path (only by `orphaned_file/private/@copvar/canonicalize.m`, off the path); the help (lines 29-31) says they stay empty – dead code – read.

**C3** `sopvar/@sdopvar/mtimes.m:1517` (SS; line 1522 AT), local function lrmultiply – `Bout = sparse(...)` is overwritten by `Bout = B*LRkron` (line 1522) before any read – dead code – semantics – fix: delete line 1517.

**C4** `sopvar/@sopvar/ctranspose.m:51-57, 61-62` (MMP, SS) – the field assignments to At are discarded by the constructor call at line 115 (the comment at lines 49-50 says so); S3Var and S3Vardummy are never read – dead code – semantics.

**C5** `sopvar/@sdopvar/sdopvar.m:64` (SS) – the help says "Zd is a vector of decision variables, sorted"; the constructor stores Zd as given (line 192) and merges keep first occurrences (`unique(...,'stable')`, `sopvar/@sdopvar/private/sync_basis.m:108`) – misleading text – read.

## LPI programming, executives, settings

**C6** `executives_sopvar/PIETOOLS_auto_execute_sop.m:151-153, 191, 196` (MMP) – with no executive flag in the workspace the prompt offers 'stability' and 'stability_dual', but the dispatch tests 'PIE2PDEstability' and 'PIE2PDEstability_dual' – error on a reachable path: these two answers end in "Unknown executive" (line 241) – read – fix: accept the prompt names in the two tests.

**C7** `settings_sopvar/settings_2d_stock.m:143-148` (MMP) – use_sosineq is overwritten from the base workspace (sosineq_on, else use_sosineq, by evalin), so a workspace variable changes the returned struct; the levels light, stripped and custom set use_sosineq = 1 themselves (line 77), and 1 makes the slack role error (`executives_sopvar/utility_functions/poslpivar_settings_2d_sop.m:126-127`) – error on a reachable path (2-D PIE, settings.sop.slack = 'stock'); the result depends on the base workspace, where C24 leaves sosineq_on = 0 – read.

**C8** `executives_sopvar/PIETOOLS_Hinf_estimator_sop.m:161-166` (MMP) – the 2-D catch branch repeats getObserver_sop with struct('method','qp'), the default method of getObserver_sop (`executives_sopvar/utility_functions/getObserver_sop.m:92`), so the retry repeats the computation that failed – no effect – read.

**C9** `lpi_programming_sopvar/misc/lpi_classify_sop.m:84-85` (MMP) – both branches of the if/else are `cl.obj = double(lpigetsol_sop(prog,obj))` – dead code – semantics.

**C10** `executives_sopvar/utility_functions/pie_disc_sop.m:231` (MMP), local ENTRY – the 'double' branch and the fall-through (line 232) both return X(i,k) – dead code – semantics.

**C11** `executives_sopvar/PIETOOLS_H2_control_sop.m:142`, `PIETOOLS_H2_estimator_sop.m:149`, `PIETOOLS_Hinf_control_sop.m:153`, `PIETOOLS_Hinf_estimator_sop.m:179` (MMP) – info.inv holds the third output of getController_sop or getObserver_sop, which under the default method 'qp' is a struct with the field qp (diagnostics of the quadratic program), not data of an inverse of P – misleading field name – read.

## parser, GUI, root scripts

**C12** `GUI/parsePDEParams.m:699-706` (MMP, SS, DJ; 2021b), local eval_b – any evaluation error other than one exact MATLAB message text returns 0, so a coefficient that fails to evaluate becomes 0 without notice; line 693 reads exp(end) before the isempty test (line 696), an index error for an empty coefficient – silent wrong result (GUI input) – read (line 693: semantics) – reachability: PIETOOLS_PDE_GUI.mlapp calls parsePDEParams.

**C13** `pietools_path_update.m:39-40` (DJ) – adds pwd and all subfolders, not the folder of the file (help, lines 6-7: pwd must be the PIETOOLS folder) – silent wrong result: called by name from another folder, it adds that folder tree – semantics – fix: fileparts(mfilename('fullpath')).

**C14** `PIETOOLS_PDE.m:246, 252` (YP) – `solution.timedep.pde{1,2,3}=solution.timedep.primary{2,3,4}` addresses one element with three subscripts on each side; primary is 1 x 4 in 2-D (PIESIM.m help), so the right side is out of range; line 252 then tests length(pde)==2; X (line 239) is unused – error on a reachable path: Step 4b, 2-D PDE, opts.plot not 'yes' (the script sets 'yes', line 176) – semantics (indexing), read (size of primary).

**C15** `parser/legacy/@state/int.m:97` (SS) – `terms(T,objC)` passes the whole state vector where line 87 passes objC(i); line 89 (`opvar T;`) needs the retired 'opvar' class and errors first; commented merge-conflict markers (<<<<<<<, =======, >>>>>>>) and dead code remain at lines 56-71 – error on a reachable path (legacy state interface, a state that depends on the integration variable) – read.

**C16** `parser/legacy/@sys/addequation.m:40-41` (SS) – a 'state' equation is converted by state2terms but never appended; state2terms errors first (`parser/legacy/@state/state2terms.m:39`, `opvar T;`) – error on a reachable path (legacy sys interface); without that error the equation would be dropped without notice – read.

**C17** `parser/legacy/@sys/convert.m:56-58, 65-78` (SS) – type 'pie' with convertTo 'pie' returns before assigning the first output params; convertTo 'ddf' with a type other than 'nds' or 'dde' assigns no output and prints "Conversion to ddf was successful" – error on a rare path ("output argument not assigned" when outputs are requested), false message – read.

**C18** `pietools_update.m:170` (DJ) – the mkdir-failure error formats with `c`, assigned only at line 179 (`d` holds the folder) – error on a rare path: an undefined-variable error, or a base-workspace c (the file is a script), replaces the message – read – fix: d for c.

**C19** `pietools_update.m:126, 234-235` (DJ) – default_dir is assigned only in the fresh-install branch; the savepath fallback reads it in the update branch as well – error on a rare path (update, savepath fails, folder dialog cancelled) – read.

**C20** `parser/internal/display_PDE.m:553-643` (DJ), local construct_term – the output use_Cij is assigned only inside `if isfield(PDE_term,'C') && ~isempty(PDE_term.C)`; the caller (line 275) requests it; the input `name` of display_PDE is unused (help, line 8) – error on a rare or unreachable path (a term without C; initialize assigns C) – read.

**C21** `GUI/isvalidTerm.m:60-83` (MMP, SS, DJ; 2021b) – no caller (none in a .m file or in PIETOOLS_PDE_GUI.mlapp); coeff is made a string (line 60) before `['(',coeff,')*',...]`, which yields a 1 x 7 string array, so `if contains(coeff,'theta')` (line 66) is never true and the `&&`, `||` at lines 71 and 83 error on non-scalar operands; inputs derOrder, highestDer, a, b unused – dead code – semantics.

**C22** `parser/legacy/@terms/identifyEqnType.m:40-51` (SS) – no caller; the first loop has empty branches; statevec.state.type and statevec.diff are rejected by the state subsref (`parser/legacy/@state/subsref.m:65, 82`), so a call errors – dead code – read.

**C23** `parser/legacy/@terms/serializeStateVec.m:37` (SS) – no caller; `zeros(length(statevec):4200)` passes the vector n:4200 as a size list (an array of 4201-n dimensions), an out-of-memory error if called – dead code – semantics.

**C24** `PIETOOLS_DDE.m:282` (DJ), `PIETOOLS_DDF.m:272` (PIETOOLS Team; 2021b) – sosineq_on = 0 is set in the base workspace and read by no executive (`executives_sopvar/utility_functions/exec_sop_settings.m:134` rejects a settings field of that name); settings_2d_stock reads the workspace variable (C7) – no effect on the DDE or DDF run; it changes later 2-D settings of the session – read.

**C25** `pietools_update.m:196-207` (DJ) – startup.m is opened in append mode and fid is closed only when the restorepath line is absent – no effect on results (a file handle left open) – read.

**C26** `parser/@pde_struct/expand_delays.m:547` (DJ), local add_delay_state – `PDE.BC_tab(2+nvars+tau_var_idx) = 1` is a linear index into the whole table (the last row was intended); expand_delays rebuilds BC_tab through initialize (line 168; `parser/@pde_struct/initialize.m:530-533, 641-642`) – no effect – semantics (index), read (rebuild).

**C27** `parser/pde_var.m:293` (DJ) – local getPDEvarID has no caller (the call at line 258 is commented out) – dead code – read.

**C28** `parser/@pde_struct/times.m:56-59`, `parser/@pde_struct/mtimes.m:65-67, 84-185` (DJ) – branches on isa(...,'opvar'), a class absent from the repository – dead code – read.

**C29** `parser/@pde_struct/times_pdes.m:180` (DJ), local times_term – obj2_size reads PDE.(obj1){obj2_num}; obj1 = obj2 = 'x' by lines 163-178 – no effect – read.

**C30** `parser/@pde_struct/combine_vars.m:194` (DJ) – `fprintf('\n','No variables have been merged.\n')`: the format has no conversion, so the text argument is ignored and only a newline prints – no effect (message lost) – semantics – fix: one format string.

**C31** `PIETOOLS_PDE.m:11-12` (SS) – the description says the script determines stability and the H-infinity norm and designs an observer; Step 3 runs no executive by default (all calls commented out; OUTPUT, line 21, says none) – misleading text – read.

**C32** `parser/legacy/@state/state.m:56, 116, 122-123` (SS) – the declared method signature state2terms(obj,operator,var,val) differs from state2terms.m (obj only); the constructor accepts 4 inputs (line 116), yet the error for 5 or more says "only takes 3 inputs" – misleading message – read.

**C33** `parser/@pde_struct/display_PDE_crude.m:177, 179, 199` (DJ) – the first term, and every factor after the first, do not request use_Cjj (only line 197 does), so the hint "Call PDE.C{i,j}" (line 216) prints only if a later term uses C_ij – misleading (missing) message – read.

**C34** `parser/@pde_struct/mpower.m:46-47` (DJ) – accepts deg = 0 (message: "nonnegative integer"), while power rejects it (`parser/@pde_struct/power.m:50-51`, "Degree-0 powers ... not supported") – misleading message (power still raises the error) – read.

## converters

**C35** `converters/internal/convert_PIETOOLS_NDS2DDF.m:112-114` (MMP) – the distributed-delay coefficients (Adi, B1di, ..., C2dei) go to a local Cvd{i}, never to DDF.Cvd{i}; DDF_out (line 129) is unused; the header (lines 14-16) states the omission – silent wrong result: the DDF and the PIE of an NDS lose every distributed-delay term – read – reachability: `converters/convert_PIETOOLS_NDS.m:97`, `PIETOOLS_DDF.m:306` – fix: DDF.Cvd{i}.

**C36** `converters/legacy/initialize_PIETOOLS_PDE_terms_legacy.m:158-181, 289-302` (MMP, SS; 2021b) – nx comes from ODE.A only if one of nx, Bxu, Bxw, Bxr, Cy, Cz, Cv is defined: an ODE given by A alone gets nx = 0, and A is replaced by zeros(0) with the notice "A has incorrect dimension" (lines 368-370); A also overrides an explicit nx; nr is read as size(PDE.PDE.Crp,1) or size(PDE.PDE.Drb,1), which count cells of terms – silent wrong result (an ODE dropped, a wrong nr; console notice only) – read – reachability: the struct format with fields n, ODE, PDE (`converters/initialize_PIETOOLS_PDE.m:69-73`).

**C37** `converters/convert_PIETOOLS_DDF.m:101, 119, 131`, `converters/legacy/convert_PIETOOLS_DDE2PIE_legacy.m:108, 145-148` (DJ, 12/07/2024) – the kernels are substituted and integrated in s1, while `PIETOOLS_DDE.m:163-178` and `PIETOOLS_DDF.m:144, 424-429` document them as functions of pvar s; a kernel in s is not integrated, and line 119 applies double to a polynomial in s – error on a reachable path (distributed-delay kernels that depend on s, as documented; constant kernels unaffected) – read.

**C38** `converters/internal/minimize_PIETOOLS_DDE2DDF.m:103-104`, `converters/internal/convert_PIETOOLS_NDS2DDF.m:99-100` (MMP) – with no delay both call convert_PIETOOLS_DDE, which calls minimize_PIETOOLS_DDE2DDF again (`converters/convert_PIETOOLS_DDE.m:86`): unbounded recursion; DDF_out is never assigned in that branch – error on a rare path (a DDE or NDS without delays; recursion limit) – read.

**C39** `converters/initialize_PIETOOLS_NDS.m:898` (MMP) – A0 is never defaulted or checked, so `NDS_out.A0=A0` fails without A0; the Ei/Edi (lines 222-225), C1ei/C1dei (442-445) and C2ei/C2dei (479-482) branches of the dimension searches are unreachable (the outer conditions, lines 200, 420, 457, omit these fields), and the last two assign nz, not ny; the numel errors for C1ei, C2ei, Edi, C1dei, C2dei (lines 167-184) name "Ei"; line 357 tests 'Bdi', not a field name, where B1di is meant – error on a rare path (A0 omitted), misleading messages – read.

**C40** `converters/initialize_PIETOOLS_DDF.m:118, 350-469` (MMP) – `nx=size(B1,1,'var')` passes the char 'var' as a dimension argument of size (reached with B1 and without A0); the checks `if size(X)~=[r,c]` (Cr, Br1, Br2, Drv, Cv, Cvd, C1, D11, D12, D1v, C2, D21, D22, D2v) are true only if both dimensions differ – error or a wrong nx on a rare path; a block with one wrong dimension passes the check – semantics.

**C41** `converters/legacy/convert_PIETOOLS_PDE_2D.m:583-616, 927`, `converters/legacy/convert_PIETOOLS_PDE_nonlinear_2D.m:438-471, 844`, `converters/legacy/2D/Compute_Tmap_PDE_2D_Separable.m:545, 1230`, `converters/legacy/2D/PIETOOLS_FTC_Expansion_2D.m:229-247` (DJ), `converters/legacy/2D/compute_PIETOOLS_PDE_periodic_map_2D.m:149-214`, `converters/internal/compute_PIETOOLS_PDE_periodic_map.m:190-208` (DB) – construct opvar and opvar2d objects; no 'opvar' or 'opvar2d' class exists in the repository (likewise `parser/legacy/@state/state2terms.m:39`, `parser/legacy/@state/int.m:89`, and the property of class opvar at `parser/legacy/@terms/terms.m:38`) – error on a rare path: these converters are called only by each other (entry: a direct call of convert_PIETOOLS_PDE_2D) – read.

**C42** `converters/legacy/2D/Compute_Tmap_PDE_2D_Separable.m:891-892, 960-961, 1409, 1444, 1453` (DJ) – compute_input_diffs indexes Dop_x, Dop_y (one entry per state component, lines 891-892) with tmp_log (one flag per row block, 4 x 1); print_reorder_diff_summary reads varname_LHS_t (undefined; varnames_LHS_t elsewhere), overwrites Ldiff for a second variable, and computes ceil(log(dval(jj)/log(10))) – unreachable at HEAD (the file errors first, C41) – read.

**C43** `converters/legacy/2D/expand_PIETOOLS_PDE_BCs.m:545` (DJ) – the integral-removal loop looks up Rindx in PDE.([Robj,'_tab']) with Robj left from the last term of the loop at lines 211-218 ('w' or 'u' possible), for a state term – unreachable at HEAD (called only by the converters of C41) – read – fix: PDE.x_tab.

**C44** `converters/legacy/convert_PIETOOLS_PDE_batch2terms_legacy.m:4` (SS) – the function line names convert_PIETOOLS_PDE_batch2terms; MATLAB uses the file name; no caller – no effect – semantics.

**C45** `converters/legacy/convert_PIETOOLS_PDE_2D.m:1942`, `converters/legacy/convert_PIETOOLS_PDE_nonlinear_2D.m:1649` (DJ) – a local print_convert_summary(PIE,obj,tdiff_list) duplicates, with another signature, `converters/print_convert_summary.m` (PIE,obj,aux_idcs,F_aux); the local function takes precedence in each file – no effect (two definitions) – read.

**C46** `converters/legacy/initialize_PIETOOLS_PDE_terms_legacy.m:75` (MMP, SS; 2021b), `converters/legacy/convert_PIETOOLS_PDE_terms_internal.m:39-44` (SS) – header notes say "The workspace is cleaned at the end of the file" (both are functions); the second keeps a developer note asking whether a PIE structure should be created; convert_PIETOOLS_PDE_terms_internal has no caller – misleading text – read.

**C47** `converters/operator_converters/dpvars2sdvars.m:64` (DJ) – the license title line reads "PIETOOLS - quad_expand" – misleading text – read.

## PIESIM

**C48** `PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_PI_discretize.m:61` (YP) – a one-variable kernel counts as a function of s when the variable name has at most 2 characters, else of theta; pie2piesim keeps the PIE variable names (dummy = name + '_dum', `converters/operator_converters/pie2piesim.m:364-368`), so a kernel in a primary variable named with 3 or more characters (e.g. eta) is sampled as a function of theta; two-variable kernels (line 56) take the orientation from the alphabetical order of R.varname – silent wrong result (transposed kernel samples) – read – reachability: 1-D integral kernels (PIESIM_3PI_opint_discretize), variable names of 3 or more characters; s1, s1_dum are classified correctly.

**C49** `PIESIM/PIESIM_routines/PIESIM_Setup/PIESIM_input_check.m:608-609` (YP), PIE input, 1-D – the control-input count tests ~isempty(PIE.w_tab) before reading PIE.u_tab – silent wrong result: a PIE with distributed control inputs and no disturbance gets nux = 0, the inputs counted as finite-dimensional – read – fix: test PIE.u_tab.

**C50** `PIESIM/PIESIM_routines/PIESIM_Setup/PIESIM_initial_setup.m:187, 252` (YP), 2-D corner values – `if (length(uinput.w{i}==1))` (u: line 252) takes the length of a logical array, nonzero for any nonempty input, so the branch for 2 x m samples is never taken and the whole sample array is scaled by the corner value; in the sample format [t; w] (`PIESIM/examples_pde_library_PIESIM_1D.m:217`) both branches scale the time row – silent wrong result (2-D, separable corner input given as samples) – semantics (condition), read (format).

**C51** `PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_discretize_icf_2D.m:186, 197, 202, 207` (YP) – nu is used (`coeff.u=zeros(Nforce,nu)`) but never assigned (only nw, line 65); the u branch passes uinput.uspace(k), a 1 x 1 cell, where the w branch passes uinput.wspace{k}, and lacks the non-separable fallback of the w branch (lines 157-158) – error on a reachable path: every 2-D simulation with control inputs (PIESIM_disturbance_check sets uspace) – read – fix: nu = psize.nu; uspace{k}.

**C52** `PIESIM/PIESIM_routines/Time_Integrators/PIESIM_int_symbolic.m:71, 99-135` (YP) – reads psize.no, which no file assigns (psize.no0 exists: `PIESIM_input_check.m:264, 296, 428`); lines 99-135 index uinput.w(i), wdot(i), u(i), udot(i) with parentheses, while PIESIM_initial_setup stores cells (line 143) – error on a reachable path: opts.intScheme = 2 with a diagonalizable Atotal (the other case: C55) – read.

**C53** `PIESIM/PIESIM_routines/PIESIM_Setup/PIESIM_input_check.m:150, 317` (YP) – (a) PDE.dim == 0: initialize_PIETOOLS_PDE(psize,true) receives psize (struct() with the field dim, lines 133-134) instead of the augmented PDE and errors (`converters/initialize_PIETOOLS_PDE.m:75`); (b) 2-D PDE without uinput.ic: line 317 uses n, never assigned (ns, line 233, holds the count) – error on a reachable path (an ODE-only pde_struct; a 2-D simulation without initial conditions) – read – fix: (a) PDE for psize; (b) n0+ns.

**C54** `PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_discretize_icf.m:148, 166` (YP) – passes uinput.wspace(k), uspace(k), 1 x 1 cells (`PIESIM_disturbance_check.m:261-302` assigns wspace{k}), to PIESIM_NonPoly2Mat_cheb, which applies subs to Rop(m,i) (`PIESIM_NonPoly2Mat_cheb.m:59`) – error on a reachable path unless subs accepts a cell first argument (no subs.m outside class folders in the repository; toolbox dispatch not checked) – read – reachability: 1-D, an input that depends on s (psize.nwx or nux > 0) – fix: wspace{k}, uspace{k}.

**C55** `PIESIM/PIESIM_routines/Time_Integrators/PIESIM_time_integrate.m:107` (YP) – the BDF2 fallback (Atotal not diagonalizable) calls PIESIM_int_bdf with 5 arguments; it takes 6 (gridall before Dop), so Dop is unset inside (`PIESIM_int_bdf.m:116`) – error on a rare path (opts.intScheme = 2, Atotal not diagonalizable) – semantics – fix: pass gridall.

**C56** `PIESIM/PIESIM_routines/PIESIM_Setup/PIESIM_input_check.m:530, 1085` (YP) – (c) the DDE branch compares the number of exact solutions with psize.n0+ns, ns unassigned in that branch (error with opts.ifexact and uinput.exact); (e) the legacy terms check (local PIESIM_input_check_PDE_terms_legacy) compares with psize.n0+psize.n, a vector, so the elseif fires only if the count differs from every entry (ns, line 1078, unused there) – error on a rare path; a count mismatch passes – read (e: semantics).

**C57** `PIESIM/PIESIM_routines/Postprocessing/PIESIM_transform_to_solution.m:197, 328` (YP) – the NaN branch (diverged solution) tests ns>0, ns unassigned in the function; line 328 writes solution.timepde.pde (timedep intended), a stray field; the input PIE is unused (help, line 19; likewise the 2-D file) – error on a rare path (a diverged simulation) – read.

**C58** `PIESIM/PIESIM_routines/Postprocessing/PIESIM_plot_solution.m:229`, `PIESIM/PIESIM_routines/2D/Postprocessing/PIESIM_plot_solution_2D.m:220` (YP) – finite-dimensional observed outputs 1-7 are colored colors_y(psize.no0-ii+1,:) (7 rows: out of range for no0 > 7) and labelled $z_i$, the regulated-output symbol; the 1-D file returns (line 149) before figs is assigned when opts.plot is not 'yes' (PIESIM.m:285 requests no output) – error on a rare path (more than 7 observed outputs, BDF, plotting on), misleading legend – read.

**C59** `PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_1Dto2D2Mat_cheb_2D.m:93-108`, `PIESIM_2Dto1D2Mat_cheb_2D.m:92-107`, `PIESIM_9PI2Mat_cheb_2D.m:108-145` (YP) – A_opmult and A_opint_block_* are assigned only for nonempty kernels but summed unconditionally (undefined, or left from the previous block) – error or a stale sum on a rare or unreachable path (pie2piesim fills zero kernels with zero matrices, `converters/operator_converters/pie2piesim.m:241-249`) – read.

**C60** `PIESIM/PIESIM_routines/PIESIM_Setup/transl_mod.m:77-101` (YP; MMP, SS) – translates the kernels to the interval I but leaves T.I at the original interval – no effect in PIESIM (no routine reads .I after rescalePIE); the struct carries a stale interval – read.

**C61** `PIESIM/PIESIM_routines/PIESIM_Setup/PIESIM_disturbance_check.m:57` (YP) – `nargin=length(varargin)` shadows the function nargin with the same value – no effect – read.

## PIESOS

**C62** `PIESOS/polyopvar/@intop/isvalid.m:291` (DJ), local is_consistent – `safe_vars = [safe_vars,find(safe_idcs)]` appends positions in the reduced var_list, not the variable numbers rem_vars(safe_idcs), from the second pass on – silent wrong result: for chains of three or more dependent integration limits a variable can be marked resolved wrongly, and isvalid returns a wrong verdict – read – fix: rem_vars(safe_idcs), taken before line 292.

**C63** `PIESOS/polyopvar/@quadPoly/horzcat.m:51`, `PIESOS/polyopvar/@quadPoly/vertcat.m:49` (SS) – `elseif isnumeric(X)`: X is never assigned (varargin{k} intended) – error on a reachable path: concatenation of a quadPoly with a numeric array – read – fix: isnumeric(varargin{k}).

**C64** `PIESOS/custom/filter_Pcell.m:115` (CR) – the second output reads entry.coefficient, which errors for a numeric Pcell{1,1}, a block type the first loop filters (line 89) and the help admits – error on a reachable path (numeric Pcell{1,1}) – read.

**C65** `PIESOS/internal/otimes_params.m:656` (DJ), local otimes_quadPolyCell – the check of B tests `~isa(A,'struct')` (B intended) – error on a rare path: an invalid B with a struct A passes the check and fails later with a less informative error – read – fix: isa(B,'struct').

**C66** `PIESOS/polyopvar/@quadPoly/subs.m:73` (AT) – `if size(vars) ~= size(vals)` is true only if both dimensions differ, so 1 x 2 vars with 1 x 3 vals pass – error on a rare path: the count mismatch is not caught here – semantics.

**C67** `PIESOS/polyopvar/@quadPoly/reshape.m:1, 57` (AT) – only RESHAPE(A,M,N) is accepted (header corrected in the pass); RESHAPE(A,[M,N]), a form MATLAB arrays accept, errors at line 57 (N not supplied) – error on a rare path – read.

**C68** `PIESOS/polyopvar/@quadPoly/plus.m:51-88` (SS) – no dimension check (the comment at line 51 asks for one); the result takes F.dim, and liftIndex lifts both operands with F.dim (lines 72, 79, 85-86); for operands of different dimensions the outcome (a MATLAB size error or a result of the size of F) is not determined by reading – error or a wrong size on a rare path – read.

**C69** `PIESOS/polyopvar/@quadPoly/subsasgn_1.m:1` (SS) – defines `function P = subsasgn(P,S)` in a file named subsasgn_1.m, so MATLAB registers the method subsasgn_1, never the class subsasgn; the signature also lacks the assigned value – dead code – semantics.

**C70** `PIESOS/polyopvar/@quadPoly/quadPoly.m:98-99` (SS) – `if ~presorted && 0` never runs, so the input presorted has no effect; the branch would call clean(obj), and quadPoly has no clean method – dead code – semantics.

**C71** `PIESOS/polyopvar/@quadPoly/var_swap.m:125-126, 232` (SS) – sPosMap and tPosMap are never read; rebuildBasisFromDegrees returns posMap as empty cells – dead code – read.

**C72** `PIESOS/polyopvar/@tensopvar/isvalid.m:127-128` (DJ) – the test accepts 'sopvar', 'tensopvar' and 'intop' factors; the message says factors must be 'sopvar' – misleading message – read.

**C73** `PIESOS/custom/kernel_degree.m:28` (CR) – the license block says GPL "version 3"; every other PIETOOLS file says version 2 – misleading text (license) – read.

## demos and examples

Below, `EL/` stands for `PIETOOLS_examples/Examples_Library/` and `CL/` for `PIETOOLS_examples/Examples_Library/command_line_input_Examples/`.

**C74** `CL/stability_11.m:48` (SS) – `R = 2,9;` assigns R = 2 and displays it; the header (line 9) gives R = 2.9 – silent wrong result (the example runs at R = 2) – semantics – fix: `R = 2.9;`.

**C75** `CL/l2gain_2.m:44` (SS) – the code declares x_t = -x_ss + s*w, a backward heat equation; the inline comment says x_ss – silent wrong result or wrong comment (the intent is not stated) – read.

**C76** `CL/stability_2.m` to `CL/stability_15.m` (SS), `CL/stability_13_2.m`, `16_1`, `16_3`, `17`, `19`, `36`, `39` (DB) – `clc; clear;` runs before `if exist('params','var')` (e.g. stability_11.m:43, 49), so the params of the caller are deleted and the override documented in each header never applies; stability_15.m:65 also sets k, aa, II, g, which no equation uses – silent wrong result (user parameters ignored) – semantics – stability_1, 13_3 and 18 keep params.

**C77** `CL/stability_7.m:48` (SS) – c is computed from lam on the line before the params loop, so a params override of lam leaves c unchanged – silent wrong result once C76 is fixed (masked at HEAD) – read.

**C78** `EL/Stability_Examples/PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat.m:87-96` (DJ) – c is computed from lam (line 90) before the params loop (line 96); `EL/Stability_Examples/PIETOOLS_PDE_Ex_Diagne.m:81-92, 130-136` (DJ) – the first params loop precedes the defaults (lines 91-92) and is overwritten, the second follows the derived coefficients (Lm, al1-al3, a21, a32, a13, K), so overrides of Hs, Vs, Bs, g, Cf do not reach the PDE – silent wrong result (params override ignored) – read.

**C79** `EL/2D/PIETOOLS_PDDE_Ex.m:82-93` (DB) – tau is set (line 93) after the params loop, so params cannot change it; `pvar s t` (line 91) is unused – silent wrong result (tau override ignored) – read.

**C80** `EL/2D/PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D.m:126-127` (DJ) – a, b are parameters (line 100), but these BCs evaluate x2, x3 at s1 = 1, not at b – silent wrong result for b ≠ 1 – read – fix: b for 1.

**C81** `EL/PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion.m:87` (DJ) – GUI = 1 loads PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq_GUI.mat, the file of another example – silent wrong result (the GUI shows another PDE) – read.

**C82** `PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch5_PIEs.m:147` (DJ; operators DB) – `PIE3 = piess({T,0,Tu},A,{0,B1},...)` puts B1 = s(2-s) on u and 0 on w, while T, Tu, A, C1, C2, D12, D22 are those of the Sec. 5.4 PDE (lines 83-95), in which s(2-s) multiplies w – silent wrong result if PIE3 is meant to be that PDE (PIE3 and the closed loops of lines 154-174 are then another system) – read – fix: {B1,0}.

**C83** `EL/Stability_Examples/PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m:119-122, 150` (DJ) – unreferenced older copy of PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m; the term format couples x(0) (loc = 0, no D; x(0) = 0 by the BC) where the comment and the batch format couple x_s(0); GUI = 1 loads PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das_GUI.mat, which does not exist – silent wrong result (term format: the coupling vanishes), error (GUI) – read.

**C84** `EL/Stability_Examples/PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq.m:191` (DB), `EL/Stability_Examples/PIETOOLS_PDE_Ex_Saba.m:110` (DJ) – the error branch for an unknown INDEX calls num2str(indx); the input is index – error on a rare path: an undefined-variable error replaces the message – read – fix: index.

**C85** `PIETOOLS_demos/DEMO3_poincare_inequality.m:100` (MMP) – `st = settings_sopvar_heavy;` is never used: lpiprogram (line 103) gets no settings, so the program runs at the default tier 'light'; lines 79-82 keep a commented-out copy of the BVP declaration – dead code – read.

**C86** `PIETOOLS_demos/DEMO5_Hinf_optimal_estimator.m:160` (MMP) – `opts.method='qp'` lacks a semicolon (opts prints); the same opts goes to PIESIM (line 194), which ignores the field – no effect (console output) – semantics.

**C87** `EL/PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion.m:104`, `EL/PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion_Control.m:99` (MMP) – set st.sop.deg_Q, read only by PIETOOLS_stability_sop (line 173), in H2 synthesis examples; neither file has a params loop (the help of the _Control file states PARAMS is unused) – no effect – read.

**C88** `EL/2D/Hinf_Gain_Examples/PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DDDD.m:155-161`, `_DDDN.m:155-160`, `_DNDN.m:156-161` (DJ) – r_lim and L2_gain are computed and neither returned nor displayed – dead code – read.

**C89** `CL/stability_16_1.m:86` (DB) – the last BC lacks `==0`; pde_struct vertcat turns a term concatenated with equations into term == 0 (`parser/@pde_struct/vertcat.m:105-111`), so the BC holds; the missing semicolon prints eq_BC, as in stability_13.m:59 (SS) and stability_17.m:63 (DB) – no effect (console output) – read.

**C90** `EL/PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance2.m:68` (DB) – `pvar s1 s2` is unused (s is used) – dead code – read.

**C91** `CL/l2gain_1.m` to `CL/l2gain_4.m`, `CL/stability_1.m` to `CL/stability_15.m` (SS) – blocks of commented-out code in the retired state()/sys() format (e.g. stability_13.m:79-83) – dead code (comments) – read.

**C92** `CL/stability_8.m:1, 62` (SS; line 62 DJ), `EL/Stability_Examples/PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m:25, 119` (DJ) – use +Bxr*x_s(a), with the header claim "stable for lam < pi^2" at lam = pi^2-1; `EL/Stability_Examples/PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m` switched to -Bxr on 11/19, with the log statement that +Bxr is stable only for lam < -5.7799 (lines 38-39, 68-72) – misleading text (stability claim; two signs among the copies) – read (stability not re-derived).

**C93** `CL/stability_10.m:11` (SS) – the header block gives x(s=0) = xo; the code (line 58) and line 3 give x(0) = -xo – misleading text – read.

**C94** `PIETOOLS_demos/DEMO12_2DControl.m:110-112` (DB) – prints "||Z - KP|| = %.3e" with full(norm_err.C{1}.params{1}), the first parameter of the first block of err*err', not a norm (lines 108-109 say so) – misleading message – read.

**C95** `PIETOOLS_demos/DEMO13_Hinf_OFcontrol.m:220` (DB) – prints the label 'Hinf_OFcontrol_demo12' (the file is DEMO13 since the renumbering of 11/18) – misleading message – read.

**C96** `PIETOOLS_demos/DEMO9_L2_gain_analysis.m:32` (DJ) – the license title line reads "PIETOOLS - DEMO8" – misleading text – read.

**C97** `EL/Hinf_Optimal_Controller_Examples/PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_with_Control.m:96` (DJ) – z = int (1-s) v2 ds has no control term (D12 = 0), a singular H-infinity synthesis problem; the header before 11/20 gave z = [u; int 0.5(1-s)^2 v2 ds] – misleading example (formulation) – read.

**C98** `EL/Hinf_Gain_Examples/PIETOOLS_PDE_Ex_Parabolic_Eq_with_Disturbance.m:149`, `EL/Hinf_Gain_Examples/PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Disturbance.m:149` (DJ; st MMP) – print "With sedumi solver, this example requires at least settings veryheavy" while st (line 153) is 'stripped' or 'light' – misleading message – read.

**C99** `EL/2D/Hinf_Optimal_Observer_Examples/PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DNDN_EdgeObserve.m:112-113` (DJ) – the file name and the commented term-based block observe the edges; the code senses integrals over s2 and s1 (header corrected on 11/20) – misleading name – read.

**C100** `PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch6_PIESIM.m:272` (DJ) – the x3 plot title shows s_{3}; x3 is a function of s2 – misleading text – read.

## User-visible strings naming retired functions, files, dates or reports

All: misleading text or message; evidence read.

**C101** `sopvar/@copvar/subsref.m:104-117`, `sopvar/@copvar/subsasgn.m:113-151` (MMP); `sopvar/@copvar/getfield_opvar1d.m:79, 132`, `getfield_opvar2d.m:89, 124`, `setfield_opvar1d.m:72, 98, 104`, `setfield_opvar2d.m:96, 177` (DJ); the @cdopvar twins (`sopvar/@cdopvar/subsref.m:104-117`, `subsasgn.m:115-153`, `getfield_opvar1d.m:85, 138`, `getfield_opvar2d.m:94, 129`, setfield_opvar1d, setfield_opvar2d) – the errors and warnings of the opvar-field interface name the retired classes 'opvar', 'opvar2d', 'dopvar', 'dopvar2d'; the warning copvar:opvarField (cdopvar:opvarField) is issued at each opvar-style field access.

**C102** `executives_sopvar/utility_functions/exec_sop_settings.m:131` – "ST should be a settings struct (lpisettings)"; `executives_sopvar/utility_functions/settings_2d_sop.m:74` – "an lpisettings struct" (MMP) – lpisettings is retired (no file).

**C103** `executives_sopvar/utility_functions/exec_sop_settings.m:139-140, 144-145, 158-160` (MMP) – "settings.%s was retired on 11/12/2026" (and the settings_nd message with the same date), "(settings report 11/16/2026, Sec. 5)": dates and an unpublished report.

**C104** `lpi_programming_sopvar/poscopvar.m:258`, `lpi_programming_sopvar/misc/copquadvar.m:287` (MMP) – the 'psatz' error names poslpivar_2d (retired) for psatz = 2 ("... ball, not the box").

**C105** `executives_sopvar/utility_functions/poslpivar_sop.m:129` – "poslpivar degrees describe one spatial variable per space" (poslpivar retired); `executives_sopvar/utility_functions/poslpivar_settings_2d_sop.m:127` – "use_sosineq selects lpi_ineq_2d" (lpi_ineq_2d retired) (MMP).

**C106** `executives_sopvar/utility_functions/settings2possopvar.m:188-191, 249-252` (MMP) – lossy-report texts, issued as warnings by poslpivar_settings_2d_sop (lines 194-198, 269-272), describe the behaviour of possopvar (in orphaned_file/private, off the path); lines 235-238 name poslpivar_2d (report.lines, not displayed).

**C107** `sopvar/misc/dpvar_op_copvar.m:149, 165` (MMP) – "(legacy: nnz(dim(:,1)) > 1)", "(legacy: nnz(dim(:,2)) > 1)": the dim convention of the retired opvar.

**C108** `sopvar/misc/mat2copvar_grid.m:164` (MMP) – "refused here as mat2opvar refuses Q1": mat2opvar is retired.

**C109** `sopvar/misc/cat_copvar_grid.m:107, 119` (MMP) – "convert others first, e.g. with opvar2copvar": opvar2copvar exists but converts objects of the retired 'opvar' class.

**C110** See also lines of the help: `sopvar/misc/gk_grid.m:38`, `sopvar/misc/chebfit_mono.m:28` – GETCONTROLLER_DIRECT_SOP (cleaning/private/dead_code); `converters/operator_converters/opvar2copvar.m:45` – COPVAR2OPVAR, OPVAR2D2COPVAR (dead code), RAND_COPVAR (testing_codes/private); `converters/operator_converters/pie2copvar.m:34` – COPVAR2OPVAR, COPVAR2OPVAR2D (dead code) (MMP).

**C111** `PIESOS/custom/DMB.m:110` (CR) – "Def. 4 product has an inconsistent output component count.": a definition of an external paper.

## By author

- **MMP**: C1, C2, C4, C6–C12, C21, C35, C36, C38–C40, C46, C60, C85–C87, C98, C102–C110, C112, C113 (C101: the file headers only; the message lines are from commits of Declan, see the outcome section).
- **DJ** (reported to Declan, not patched): C12, C13, C18–C21, C24–C30, C33, C34, C37, C41–C43, C45, C47, C62, C65, C72, C78, C80–C84, C88, C92, C96–C101, C115.
- **SS**: C3–C5, C12, C15–C17, C21–C23, C31, C32, C36, C44, C46, C60, C63, C68–C71, C74–C77, C89, C91–C93, C114, C116.
- **YP**: C14, C48–C61.
- **AT**: C1, C3, C66, C67.
- **DB**: C41, C76, C79, C82, C84, C89, C90, C94, C95.
- **CR**: C64, C73, C111.
- **PIETOOLS Team** (2021b import, no author named): C24 (PIETOOLS_DDF.m).

Items credited to "MMP, SS, DJ; 2021b" (C12, C21) follow the header of the imported file; the original author of the cited lines is not established.

## Not confirmed on re-reading

- degbalance_core and the possopvar documentation of the 'int'/'mult' slots (A, 11/19): possopvar.m lies in orphaned_file/private, off the path; no user-visible text is affected.
- settings_sopvar_custom deg_Q{1} = [2 2 2] against light [4 2 2] (B, 11/19): the help states the difference (`settings_sopvar/settings_sopvar_custom.m:4-5, 60`), and the lpisettings_sop tables (lines 19, 29) list only columns where custom equals light; a choice, not a defect.
- poscopvar fallback for a program without settings (B, 11/19): documented (`lpi_programming_sopvar/poscopvar.m:208-209, 220`); lpiprogram always stores settings (`lpi_programming_sopvar/lpiprogram.m:107`), so only programs built without lpiprogram reach it.
- PIETOOLS_stability_sop help listing dd1, dd12, options1, options12, override1 (B, 11/19): corrected in the pass; no such text at HEAD.
- lpisettings_sop psatz 'auto' (B, 11/19): corrected; line 66 gives {'auto','faces1','auto'}.
- getController_sop, getObserver_sop help leading with K = Z*P^{-1} (B, 11/19): corrected; the method 'qp' leads (lines 6-8).
- op2copvar_sop passing a 'dopvar' to opvar2copvar (D, 11/19): needs a 'dopvar' object, and no 'dopvar' or 'opvar' class exists in the repository; `converters/operator_converters/op2copvar_sop.m:50-51` is unreachable.
- author line "Claude" in compute_PIETOOLS_PDE_periodic_map_2D (D, 11/19): now `% Initial coding DB - 10/29/2026` (line 80), decided 11/20.
- third-party files without the GPL block (D, 11/19): decided 11/20 (no block); clencurt moved to cleaning/private/dead_code in 78641839.
- contact lines of the PIESIM files (D, 11/19): decided 11/20 (two-line support line).
- marker "% Claude, 10/28/2026" in Compute_Tmap_PDE_2D_Separable (D1, 11/20): absent at HEAD; the markers read "DJ (CC)".
- removed TODO of piesos_LF (piesos_polyvar against piesos_quadvar) (E, 11/19): a design question, no code defect.
- degree-infeasibility notes kept in PIESOS_stability_local and PIESOS_sublevel (E, 11/19): comment content, no code defect.
- Initial coding lines that name other files in PIESOS (E, 11/19): attribution, no code defect.
- no cleaning after quadPoly diff and int (E, 11/20): the coefficients are correct; zero terms stay stored; a request, not a defect.
- Cmat of size 1 x nZ^2 in the local quad_expand_quadPoly of `PIESOS/polyopvar/@intop/mtimes_tensopvar.m` (E, 11/20): the output monomials are the right basis Zt of ZL (`quadPoly(C,Zs,Zt,dim,ns,nt)`, lines 949-950), which requires 1 x nZ^2 coefficients for dim [1,nZ]; the size is correct.
- Ch7 snippet against the manual (F, 11/19): concerns the manual text, not code.
- Ch2 snippet gain values "from 5.2 to 0.78" (F, 11/19): numbers; reading cannot confirm or refute them.
- PDE_b listed in OUTPUT of 13 2-D examples (F2, 11/20): the comments were corrected in the pass; no code defect.
- "we" and a possessive in comments of DB and DJ (A2, 11/20): wording, not functional.
- 13 string entries of `lists/B_lpi_exec_settings_ledger.md` that name functions on the path, none retired (lpiprogram, lpisolve, lpidecvar, lpi_eq_sop, lpisettings_sop, LPISETTINGS_STOCK_SOP, sossolve, copquadvar, poscopvar, int_semisep, unpack_sheets, settings_2d_stock, cdopvar/verify, the stock settings): `lpi_programming_sopvar/getsol_lpivar_sop.m:87, 137-139`, `lpigetsol_sop.m:67`, `lpidecvar.m:73`, `lpiscript_sop.m:99`, `lpisetsettings.m:60`, `lpivar_cdopvar.m:150`, `misc/lpigetdual_sop.m:88, 91`, `misc/lpi_eq_cdopvar.m:108`, `misc/private/collect_eq_rows.m:103, 188`, `misc/poscopvar_lift.m:329, 333`, `misc/private/unpack_sheets.m:65`; `executives_sopvar/utility_functions/exec_tools_sop.m:157, 169`, `poslpivar_sop.m:148`, `poslpivar_settings_sop.m:96`, `poslpivar_settings_2d_sop.m:213, 276`.
