# Ledger – C_parser_gui_root (release comment pass, 11/19/2026)

127 files read, 127 edited, 0 skipped. codeident.sh on all 127 at the end: 127 OK / 0 DIFF / 0 NEW.

Abbreviations: rel = log line `% MMP, 11/19/2026: Release comment pass (comment-only).` (new = 0 files only);
sup = `% For support, contact M. Peet, Arizona State University at mpeet@asu.edu` added after the license;
ini = `% Initial coding <initials> - <date>` added from the git first add; lic = PIETOOLS GPL block added.

## Edited files

- GUI/getVarList.m – title line named getPDEstructure.m (corrected); type options 'regulated', 'observed', 'disturb', 'control' and the fallback added to INPUT; sup, ini, rel; 2 step comments.
- GUI/isvalidPolynomial.m – INPUT/OUTPUT; sup, ini, rel; 2 step comments.
- GUI/isvalidTerm.m – INPUT/OUTPUT (derOrder, highestDer, a, b marked not used); sup, ini, rel.
- GUI/parsePDEParams.m – "exclusively used by the PDE_GUI" replaced by what it does; input classes and OUTPUT; header for local eval_b; sup, ini, rel; 5 step comments.
- PIETOOLS_PDE.m – INPUT/OUTPUT; lic (between the description and the MMP log); sup, ini, rel; 4 User Manual / PIESIM file pointers removed; example count 52 -> 72 (n_examples at HEAD); stale "comment line 28 / uncomment line 65" replaced; step comments in Steps 4 and 4b.
- pietools_path_update.m – INPUT/OUTPUT; lic, sup, ini, rel; step comment.
- pietools_update.m – INPUT/OUTPUT; lic, sup, ini, rel; 2 step comments; section comment "Install latest version of SOSTOOLS" corrected (the section installs tbxmanager, which installs SeDuMi).
- parser/@npie_struct/display.m – pointers to 'get_display' and 'print_terms' removed; sup, rel.
- parser/@npie_struct/eq.m – sup, rel.
- parser/@npie_struct/expand_tderivatives.m – two pointers to `initialize` removed; sup, rel; step comment.
- parser/@npie_struct/horzcat.m – pointer to `initialize` removed; sup, rel.
- parser/@npie_struct/initialize.m – pointer to 'npie_struct' removed; sup, rel.
- parser/@npie_struct/npie_struct.m (new) – 2 change-log entries and 3 per-line markers removed; pointers to pie_struct, pie2npie_struct, convert_PIETOOLS_PDE removed; f_pde added to the class description; INPUT/OUTPUT for the constructor, get.x, the alias / vars-dom-dim / table get-set groups, empty_pie2pde, empty_pde2pie; sup.
- parser/@npie_struct/private/pie_compat_tab.m (new) – sup.
- parser/@npie_struct/private/pie_compat_vars.m (new) – sup.
- parser/@npie_struct/private/pie_space_canon.m (new) – sup; 2 step comments.
- parser/@npie_struct/subsasgn.m – sup, rel.
- parser/@npie_struct/vertcat.m – pointer to `initialize` removed; sup, rel.
- parser/@pde_struct/combine_vars.m – pointer to initialize_PIETOOLS_PDE removed; INPUT/OUTPUT for local replace_vars_terms; rel.
- parser/@pde_struct/convert.m – pointers to @pde_struct/initialize, the 'sys' convert, convert_PIETOOLS_PDE, convert_PIETOOLS_BVP (header, local BVP2pie_struct, 3 body comments) replaced by what the code does; flag 'Top' described; sup, rel.
- parser/@pde_struct/diff.m – INPUT/OUTPUT for local compute_d_obj and addvars; sup, rel.
- parser/@pde_struct/display.m – pointers to initialize_PIETOOLS_PDE removed; "phi_i" corrected to gamma_i; crude display of a PDE that is not initialized stated; rel; 2 step comments.
- parser/@pde_struct/display_PDE_crude.m – pointer to display_PDE removed; OUTPUT added; sup, rel.
- parser/@pde_struct/display_term.m – sup, rel.
- parser/@pde_struct/eq.m – sup, rel.
- parser/@pde_struct/expand_delays.m – INPUT/OUTPUT for local add_delay_state; rel.
- parser/@pde_struct/expand_tderivatives.m – rel.
- parser/@pde_struct/expand_tderivatives_noBCs.m – rel.
- parser/@pde_struct/horzcat.m – sup, rel.
- parser/@pde_struct/initialize.m – INPUT/OUTPUT for local extract_static_constraints; converter-usage sentence and two "see <local function>" pointers removed; rel.
- parser/@pde_struct/int.m – sup, rel.
- parser/@pde_struct/is_pde_term.m – sup, rel.
- parser/@pde_struct/is_pde_var.m – sup, rel.
- parser/@pde_struct/minus.m – sup, rel.
- parser/@pde_struct/mpower.m – sup, rel.
- parser/@pde_struct/mtimes.m – INPUT/OUTPUT for local mat2cell_poly; sup, rel.
- parser/@pde_struct/pde_common_basis.m – sup, rel.
- parser/@pde_struct/pde_struct.m – pointers to initialize_PIETOOLS_PDE (class header, Hidden-properties comment, constraint_ID comment) removed; INPUT/OUTPUT for the constructor and get.dim; sup, rel.
- parser/@pde_struct/plus.m – sup, rel.
- parser/@pde_struct/power.m – sup, rel.
- parser/@pde_struct/private/is_LHS_term.m – sup, rel.
- parser/@pde_struct/private/var2term.m – sup, rel.
- parser/@pde_struct/reorder_comps.m – pointer to @pde_struct/initialize removed; sup, rel.
- parser/@pde_struct/setControl.m – sup, rel.
- parser/@pde_struct/setObserve.m – sup, rel.
- parser/@pde_struct/size.m – sup, rel.
- parser/@pde_struct/static_constraints.m (new, DB) – lic, sup; pointers to 'initialize', the converter and the display routines replaced by what the function does; 4 step comments.
- parser/@pde_struct/subs.m – sup, rel.
- parser/@pde_struct/subsasgn.m – sup, rel.
- parser/@pde_struct/subsref.m – sup, rel.
- parser/@pde_struct/times.m – sup, rel.
- parser/@pde_struct/times_pdes.m (new) – change-log entry removed; pointer to expand_delays removed; sup.
- parser/@pde_struct/uminus.m – sup, rel.
- parser/@pde_struct/uplus.m – sup, rel.
- parser/@pde_struct/vertcat.m – sup, rel.
- parser/@pie_struct/ctranspose.m – sup, rel.
- parser/@pie_struct/display.m – pointers to "pde_struct/convert", 'get_display', 'print_terms' removed; sup, rel.
- parser/@pie_struct/eq.m – sup, ini (DJ - 12/31/2021), rel.
- parser/@pie_struct/expand_tderivatives.m – sup, rel; step comment.
- parser/@pie_struct/horzcat.m – sup, rel.
- parser/@pie_struct/initialize.m – sup, rel.
- parser/@pie_struct/pie_struct.m – pointers to npie_struct, npie2pie_struct, convert_PIETOOLS_PDE removed; INPUT/OUTPUT for the constructor, the alias / vars-dom-dim / table get-set groups, empty_pie2pde, empty_pde2pie; sup, rel.
- parser/@pie_struct/private/pie_compat_tab.m (new) – sup.
- parser/@pie_struct/private/pie_compat_vars.m (new) – sup.
- parser/@pie_struct/private/pie_rename_dummy.m (new) – INPUT/OUTPUT for local rename_in; sup; step comment.
- parser/@pie_struct/private/pie_space_canon.m (new) – sup; 2 step comments.
- parser/@pie_struct/subsasgn.m – sup, rel.
- parser/@pie_struct/subsref.m – sup, rel.
- parser/@pie_struct/vertcat.m – sup, rel.
- parser/internal/display_PDE.m – two pointers to initialize_PIETOOLS_PDE replaced; INPUT/OUTPUT for local construct_term; rel.
- parser/internal/print_terms.m (new) – change-log entry removed; NOTES naming the 'display' / 'get_display' callers removed; sup; 2 step comments.
- parser/internal/state2pde_struct.m – sup, rel; step comment.
- parser/internal/stateNameGenerator.m – what/INPUT/OUTPUT, lic, sup, ini, rel; step comment.
- parser/internal/terms2pde_struct.m – sup, rel.
- parser/internal/use_unicode.m (new) – pointer to 'varname2display' removed; sup; 3 step comments.
- parser/internal/varname2display.m (new) – 2 change-log entries removed; caller classes and design rationale removed from the description; sup; 2 step comments.
- parser/legacy/@state/combine.m – ini, rel; 2 step comments.
- parser/legacy/@state/diff.m – ini, rel; 2 step comments.
- parser/legacy/@state/eq.m – ini, rel; step comment.
- parser/legacy/@state/horzcat.m – what/INPUT/OUTPUT, lic, sup, ini, rel.
- parser/legacy/@state/int.m – INPUT named obj (argument objC) and said "differentiation" (integration), corrected; ini, rel; 3 step comments.
- parser/legacy/@state/isdot.m – ini, rel; step comment.
- parser/legacy/@state/isequal.m – ini, rel; step comment.
- parser/legacy/@state/ismember.m – ini, rel; step comment.
- parser/legacy/@state/isout.m – ini, rel; step comment.
- parser/legacy/@state/minus.m – ini, rel; 2 step comments.
- parser/legacy/@state/mtimes.m – ini, rel; 2 step comments.
- parser/legacy/@state/ne.m – ini, rel.
- parser/legacy/@state/numArgumentsFromSubscript.m – INPUT/OUTPUT, lic, sup, ini, rel.
- parser/legacy/@state/plus.m – ini, rel; 2 step comments.
- parser/legacy/@state/state.m – INPUT/OUTPUT for the constructor, header for set.dom / set.maxdiff; ini, rel; 2 step comments.
- parser/legacy/@state/state2terms.m – INPUT named objA (argument obj), corrected; ini, rel; step comment.
- parser/legacy/@state/subs.m – ini, rel; 2 step comments.
- parser/legacy/@state/subsref.m – ini, rel; step comment.
- parser/legacy/@state/uminus.m – description "converts objA" corrected to -objA; OUTPUT named out (sumTerms); ini, rel; step comment.
- parser/legacy/@state/uplus.m – OUTPUT named out (sumTerms); ini, rel; step comment.
- parser/legacy/@sys/addequation.m – ini, rel; step comment.
- parser/legacy/@sys/convert.m – OUTPUT listed obj only (outputs params, out); ini, rel; step comment.
- parser/legacy/@sys/disp.m – ini, rel; step comment.
- parser/legacy/@sys/getPDEparams.m – INPUT named obj (argument pdeObj); OUTPUT said 'sys' (returns a 'pde_struct'); rel.
- parser/legacy/@sys/getParams.m – ini, rel; step comment.
- parser/legacy/@sys/getStatesFromEquations.m – header for local combine_statename; ini, rel; step comment.
- parser/legacy/@sys/removeControl.m – ini, rel; step comment.
- parser/legacy/@sys/removeObserve.m – ini, rel; step comment.
- parser/legacy/@sys/removeequation.m – ini, rel.
- parser/legacy/@sys/setControl.m – ini, rel; step comment.
- parser/legacy/@sys/setObserve.m – ini, rel; step comment.
- parser/legacy/@sys/summarize_sys.m – varargin and OUTPUT added; ini, rel; 2 step comments.
- parser/legacy/@sys/sys.m – headers for the constructor, the get/set methods and local validateEquation; ini, rel.
- parser/legacy/@terms/eq.m – OUTPUT named terms_obj (obj); ini, rel; step comment.
- parser/legacy/@terms/identifyEqnType.m – ini, rel; step comment.
- parser/legacy/@terms/int.m – ini, rel; 4 step comments.
- parser/legacy/@terms/isequal.m – ini, rel; step comment.
- parser/legacy/@terms/length.m – ini, rel; step comment.
- parser/legacy/@terms/minus.m – ini, rel; step comment.
- parser/legacy/@terms/mtimes.m – ini, rel; 2 step comments.
- parser/legacy/@terms/plus.m – rel; 3 step comments.
- parser/legacy/@terms/serializeStateVec.m – what/INPUT/OUTPUT, lic, sup, ini, rel; usage note "not currently used" removed.
- parser/legacy/@terms/subsref.m – ini, rel; step comment.
- parser/legacy/@terms/terms.m – INPUT/OUTPUT for the constructor; ini, rel.
- parser/legacy/@terms/uminus.m – INPUT/OUTPUT; ini, rel; step comment.
- parser/legacy/@terms/uplus.m – INPUT/OUTPUT; ini, rel.
- parser/legacy/@terms/vertcat.m – rel; 2 step comments.
- parser/legacy/poly2double.m – INPUT/OUTPUT, lic, sup, ini, rel.
- parser/pde_var.m – call forms without type or size, and D_order as fifth argument, added to the header; sup, rel.
- parser/pielft.m – pointer "(see also 'initialize')" removed; INPUT/OUTPUT for local same_space and extract_ops; sup, rel.
- parser/piess.m – pointers to 'convert', 'piess_opvar', INITIALIZE replaced by what the code does; sup, rel.

## Skipped

None: no listed file had uncommitted edits when edited (git status checked before each folder).

## Candidates for the maintainer (functional issues, not fixed)

1. pietools_update.m: the error after a failed mkdir formats with `c`, undefined there (`d` intended).
2. pietools_update.m: `default_dir` exists only in the fresh-install branch; the savepath fallback errors in the update branch.
3. pietools_update.m: startup.m is opened in append mode; `fid` is closed only when the restorepath line is absent.
4. pietools_path_update.m: adds pwd, not the folder of the file; run from another folder, it adds that folder.
5. PIETOOLS_PDE.m, Step 4b (runs only when opts.plot is not 'yes'): `solution.timedep.pde{1,2,3}=solution.timedep.primary{2,3,4}` assigns one cell element, not three; the 2-D branch then tests length(pde)==2; `X` unused.
6. PIETOOLS_PDE.m: Step 3 runs no executive by default (all calls commented out), while the description says the script determines stability, H-infinity norm and observers.
7. GUI/isvalidTerm.m: coeff is made a string before `['(',coeff,')*',...]`, which then yields a string array; `if contains(coeff,...)` on it is true only if every element matches, so the theta/s checks after the concatenation do not fire (inferred from MATLAB concatenation rules, not run). Inputs derOrder, highestDer, a, b unused.
8. GUI/parsePDEParams.m, eval_b: `exp(end)` read before the isempty test (error on an empty coefficient); any other evaluation error returns 0 silently.
9. parser/legacy/@state/int.m: `out = [out;terms(T,objC)]` passes the whole vector objC, not objC(i). Commented merge-conflict markers (<<<<<<< / ======= / >>>>>>>) with dead code remain in its comments.
10. parser/legacy/@sys/addequation.m: a 'state' equation is converted to terms but never appended.
11. parser/legacy/@sys/convert.m: for a 'pie' system and convertTo 'pie' it returns before assigning the first output params; for convertTo 'ddf' with another type no output is assigned.
12. parser/legacy/@terms/identifyEqnType.m: no caller; reads statevec.state and statevec.diff, which 'state' subsref rejects; first loop has empty branches.
13. parser/legacy/@terms/serializeStateVec.m: no caller; `zeros(length(statevec):4200)` passes the vector n:4200 as a size.
14. parser/legacy/@state/state.m: declared private signature state2terms(obj,operator,var,val) differs from the file (obj only); 5+ inputs error with "only takes 3 inputs".
15. parser/internal/display_PDE.m, construct_term: output use_Cij assigned only when the term has a nonempty C; input `name` of display_PDE unused.
16. parser/pde_var.m: local getPDEvarID has no caller (its call is commented out).
17. parser/@pde_struct/expand_delays.m, add_delay_state: `PDE.BC_tab(2+nvars+tau_var_idx) = 1` is a linear index into the whole table; `PDE.BC_tab(end,2+nvars+tau_var_idx)` appears intended (inferred, not run).
18. parser/@pde_struct/combine_vars.m: `fprintf('\n','No variables have been merged.\n')` prints only a newline.
19. parser/@pde_struct/display_PDE_crude.m: the C_ij hint is printed only if a term after the first uses C_ij (use_Cij of the first term is not requested).
20. parser/@pde_struct/times_pdes.m, times_term: obj2_size reads PDE.(obj1){obj2_num} (harmless while both factors must be states).

## Unresolved

1. PIETOOLS_PDE.m, Step 4b, 2-D branch: the step comment states the apparent intent (states in s1, then states in (s1,s2) at the middle grid point in s2); the code does not do this (candidate 5).
2. GUI/isvalidTerm.m: the meaning of the unused inputs derOrder, highestDer, a, b is not stated (marked "not used").
3. parser/@pde_struct/expand_delays.m, add_delay_state: the encoding of loc_Robj (0 interior, 1 lower, 2 upper boundary) is inferred from its indexing.
4. parser/legacy/@terms/identifyEqnType.m: the step comment states the intended classification; the code errors on 'state' subsref (candidate 12).

## Notes

- Licence block copied from executives_sopvar/PIETOOLS_stability_sop.m, from `% Copyright (C) 2026 PIETOOLS Team` through `% authorship, and a brief description of modifications`; the end text in BRIEF.md ("authors, and a brief summary of the modification") occurs in no file.
- Initial coding lines follow the git first add: GUI/*.m and PIETOOLS_PDE.m read "DJ - 12/31/2021" (the repository import), the legacy files "SS - <date>"; the actual authors or dates may differ.
- The editor dropped the trailing blank of 9 code lines; restored, so the code lines match HEAD byte for byte apart from the marker removals in npie_struct.m and the trailing comment of constraint_ID in pde_struct.m. Four blank lines lost their indentation.
- Kept under rule 4 (new = 0): change-log entries and markers that name other functions or retired classes (e.g. "see also update to pde_var", "moved to piess_opvar", the 11/05 markers in PIETOOLS_PDE.m).
- Literature citation kept: [Jagt, 2025] in the change log of parser/@pie_struct/pie_struct.m.
