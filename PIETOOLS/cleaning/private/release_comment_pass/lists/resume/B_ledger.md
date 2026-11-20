# B – release comment pass, resume ledger (11/20/2026)

Package B: the 21 files of `B.tsv` and the 27 files of `B_runs.txt` – 44 distinct files (`@npie_struct/display.m`,
`expand_tderivatives.m`, `initialize.m` and `PIETOOLS_auto_execute_sop.m` are on both lists). All edits comment-only;
every file codeident OK after its edits. Before the first edit: the 41 files already passed were codeident OK;
`PIETOOLS_auto_execute_sop.m`, `PIETOOLS_DDE.m`, `PIETOOLS_DDF.m` were clean in git.

## Edited files – B.tsv

- executives_sopvar/PIETOOLS_H2_control_sop.m – log (6 entries) and BEGIN/END body markers removed; help: the retired settings fields dd1, dd12, options1, options12, override1, ddZ (rejected by the settings check) replaced by lpi.pos, lpi.var.deg{1}, eppos (eppos2); opvar PIE input and EXEC_PIE_SOP / EXEC_SOP_SETTINGS / EXEC_INFO_SOP / GETCONTROLLER_SOP / lpisettings_sop pointers replaced by what they hold; K described as the bounded-degree gain minimizing ||Z - K*P||^2, an approximation of Z*P^{-1} (help said K = Z*P^{-1}; the default method is 'qp'); info.inv described as diagnostics; body: "opvar or container PIE", possessive, "(EXEC_TOOLS_SOP)", "container path" comments restated; local header without LPI_INEQ_SOP; 'stock' quoted.
- executives_sopvar/PIETOOLS_H2_estimator_sop.m – as H2_control (L the bounded-degree gain minimizing ||P*L - Z||^2, an approximation of P^{-1}*Z).
- executives_sopvar/PIETOOLS_H2_norm_sop.m – log (9 entries) removed, Initial coding line cut to its date (consolidation history); help as H2_control, plus opvar2d, POSCOPVAR, SETTINGS_2D_STOCK pointers removed and settings_2d described; body: "(as the stock 1-D files)" removed, "poslpivar's default degrees" restated as the degrees {1,[1 1 1],[1 1 1]}, "(EXEC_TOOLS_SOP, eye)" removed; 3 step comments in the form switch.
- executives_sopvar/PIETOOLS_Hinf_OFcontrol_sop.m (DB; rules only) – log (DB 11/01, MMP 11/05–11/12) removed, Initial coding line kept; help: retired fields replaced by lpi.pos; opvar/opvar2d input and EXEC_PIE_SOP / EXEC_SOP_SETTINGS / EXEC_INFO_SOP / lpisettings_sop pointers removed; the NOTES pointer to CLOSEDLOOPPIE / PIETOOLS_Hinf_gain_sop turned into a See also line, the limitation kept; body: the tl.toop/opvar and OP2COPVAR_SOP history comments restated; field form added to the local header; 1 step comment.
- executives_sopvar/PIETOOLS_Hinf_control_sop.m – as H2_control; info.Psol, info.Zsol described; 2 step comments.
- executives_sopvar/PIETOOLS_Hinf_estimator_sop.m – log (8 entries) and BEGIN/END removed; help as H2_control (Z degrees lpi.var.deg{min(nv,3)}, settings_2d); L described as the 'qp' gain in one or two variables, G "(1-D)" corrected to one or two variables; body: the 2-D comment corrected (it described the container inverse of a separable P first; the code calls the default method 'qp' and repeats 'qp' after an error); 1 step comment.
- executives_sopvar/PIETOOLS_Hinf_gain_sop.m – log (9 entries) and BEGIN/END removed; help: the 'Q' form named the generalized storage T'*Q, the 'Pd' default margin stated (1e-4 on R^n, 1e-6 on the space in every variable), retired fields, opvar, POSCOPVAR, SETTINGS_2D_STOCK, EXEC_* pointers removed; 4 step comments.
- executives_sopvar/PIETOOLS_well_posedness_sop.m – log (6 entries) removed; help: retired fields replaced, eppos (also the margin of the surjectivity condition), dw_feas, eppos_stab described, opvar/EXEC_* pointers removed; body: "as the stock executive", a possessive and pointers restated; 1 step comment.
- executives_sopvar/PIETOOLS_auto_execute_sop.m (new; full pass) – log (5 entries), 18 per-line markers and 6 (was) lines removed; Initial coding line without the "stock script" history (attribution kept); header: the flags and the inputs stability_form, Hinf_gain_form, H2_norm_form described by what they do (H2 dual: 'c' followed by the letters of H2_norm_form after the first), settings (passed to every executive; tier prompt), the GUI stop and the cleared variables stated, "as its header states them" replaced; 12 step comments (runs of 43 and 20 lines).
- PIETOOLS_DDE.m (2025 file) – header: truncated sentence completed (the executives the flags select, with the forms), "one of 5 settings scripts" (retired settings_PIETOOLS_*) corrected to the tier lpisettings_sop('light') and the example settings st, the PIESIM simulation step added, "minimal Delay-Differential Equation" corrected to the minimal DDF representation, INPUT/OUTPUT added; PIETOOLS GPL block added; support line moved from the NOTES of the user-input block to the end of the help; "Initial coding PIETOOLS Team (PIETOOLS 2021b; in the repository since 12/31/2021)" added; release log line appended; body: the PIESIM.m / User's manual pointer replaced by the fields of solution, the settings-script list rewritten as the tiers, the legacy-conversion comment corrected (DDE directly to PIE, not to a DDF), library file name case corrected. Log, markers and (was) lines kept.
- PIETOOLS_DDF.m (2025 file) – as PIETOOLS_DDE.m; examples_DDF_library_PIETOOLS (no such file) corrected to examples_NDSDDF_library_PIETOOLS (the file called), NDS and DDE inputs described, "For most people, you should use PIETOOLS_DDE" removed, "the default in the DDE executable" restated as the default DDF_min = 1.
- parser/@npie_struct/display.m (DJ; newfix) – log (4 entries) and release line removed, Initial coding entry cut to "DJ, 10/30/2026: Initial coding;"; the time-derivative symbol rule of a removed entry stated in NOTES; 6 step comments.
- parser/@npie_struct/eq.m (newfix) – release line removed; 4 step comments.
- parser/@npie_struct/expand_tderivatives.m (newfix) – Initial coding line without "based on the pie_struct version"; release line removed; "default names assigned by initialize" restated; local header: PIE 'pie_struct' corrected to 'npie_struct'; 13 step comments.
- parser/@npie_struct/horzcat.m (newfix) – Initial coding line cut; release line removed; three pointers to `initialize` restated.
- parser/@npie_struct/initialize.m (newfix) – Initial coding line cut; DJ 11/08 entry (fact already in NOTES) and release line removed; 4 step comments.
- parser/@npie_struct/subsasgn.m (newfix) – release line removed.
- parser/@npie_struct/vertcat.m (newfix) – as horzcat; 1 step comment.
- parser/@pde_struct/mpower.m (newfix) – release line removed; title POWER corrected to MPOWER.
- parser/@pde_struct/power.m (newfix) – DJ 11/09 entry ("Support terms with delay") restated in INPUT; release line removed; pointer to expand_delays in the local function restated.
- parser/@pde_struct/times.m (newfix) – release line removed.

## Edited files – B_runs.txt (step comments only, unless stated)

- executives_sopvar/utility_functions/qp_copvar.m (AT) – 1.
- executives_sopvar/utility_functions/settings2possopvar.m – 1.
- lpi_programming_sopvar/misc/lift_copvar.m – 2; duplicated header line "% linear indices." of local ALPHA_SELECT removed.
- lpi_programming_sopvar/misc/poscopvar_direct.m – 7.
- GUI/getVarList.m – 3.
- GUI/parsePDEParams.m – 32 (one label per term loop, one per derivative-order parse; inserted by script).
- parser/@pde_struct/combine_vars.m – 5; convert.m – 2; display_PDE_crude.m – 3; display_term.m – 1; initialize.m – 5; pde_struct.m – 3; subsref.m – 1.
- parser/@pie_struct/display.m – 4; expand_tderivatives.m – 5.
- parser/legacy/@state/isequal.m – 1; state.m – 2; subsref.m – 1.
- parser/legacy/@sys/convert.m – 2; getPDEparams.m – 1; sys.m – 2.
- parser/pde_var.m – 5; parser/piess.m – 1.
- (@npie_struct/display.m, expand_tderivatives.m, initialize.m, PIETOOLS_auto_execute_sop.m: under B.tsv.)

## Skipped

None.

## Candidates for the maintainer (functional, not fixed)

1. PIETOOLS_auto_execute_sop.m: with no flag set, the prompt offers 'stability' and 'stability_dual', but the branches test 'PIE2PDEstability' and 'PIE2PDEstability_dual'; these two answers end in "Unknown executive" (inferred from the code, not run).
2. PIETOOLS_Hinf_estimator_sop.m, 2-D: getObserver_sop(Psol,Zsol) already uses its default method 'qp'; the catch branch repeats the same computation with struct('method','qp') – a redundant fallback.
3. H2_control, H2_estimator, Hinf_control, Hinf_estimator: info.inv holds the diagnostics of the 'qp' method (field qp), not of an inversion of P; the name suggests the 'inv' method.
4. parser/@pde_struct/mpower.m: deg = 0 passes its check (message "nonnegative integer"), but power rejects it ("Degree-0 powers ... not supported").
5. parser/@pde_struct/times.m: branches on isa(fctr1,'opvar'), a retired class, remain in the code.
6. PIETOOLS_DDE.m, PIETOOLS_DDF.m: sosineq_on = 0 is set but read by no executive (the settings check rejects a field sosineq_on; only settings_2d_stock reads the workspace variable).

## Unresolved

1. PIETOOLS_Hinf_OFcontrol_sop.m (DB): the NOTES keep, without the function pointers, the statement that the L2 gain of the closed loop is tested in the gain form 'Q' because TK carries factors of T; whether 'P' fails or is only less suited is not verified.
2. PIETOOLS_well_posedness_sop.m: the help says the LPI tests well-posedness, but not which notion a feasible LPI certifies (not stated with confidence).
3. PIETOOLS_DDE.m, PIETOOLS_DDF.m: the qualitative remarks on the levels ("works well for most problems", "often infeasible"), written for the retired settings scripts, now label the tiers of lpisettings_sop; not verified for the tiers.
4. PIETOOLS_auto_execute_sop.m: "Initial coding MMP, 10/08/2026, from a script of DJ, 01/12/2025" keeps the attribution of the removed log text; git first add: MMP 10/09/2026.

## Literature or author attributions kept

- PIETOOLS_Hinf_OFcontrol_sop.m: the LPI of Jungbae Chun and "Theorem 11(a) of Shivakumar et al." (help); "based on formulation and code of Jungbae Chun" (Initial coding line).
- PIETOOLS_DDE.m: "our Automatica Paper" in the comment of the legacy conversion.
- PIETOOLS_auto_execute_sop.m: the script of DJ it derives from (Initial coding line).

## Notes

- Line endings kept per file (CRLF and LF); the Edit tool rewrote qp_copvar.m and poscopvar_direct.m, which held 34 and 62 LF-only lines added on 11/19, as CRLF throughout; git (text=auto) normalizes on commit. Final-newline state of every file as at HEAD.
- The Edit tool drops trailing blanks of the lines it writes: one code line (PIETOOLS_DDE.m `close all; clear; %path(pathdef); `) lost its trailing blank and was restored; a check over the 44 files finds no other HEAD code line without its trailing blanks.
- No comment line follows a continuation (`...`) that HEAD did not have (checked over the 44 files). Comments placed before an `elseif` or between `switch` cases are valid MATLAB. MATLAB not run, so no syntax check.
- parser/@pde_struct/power.m: the local function complete_term carries the PIETOOLS GPL block and an Initial coding line (existing; kept).
- After the scratchpad collision warning, every final check (codeident, runs, trailing blanks over 80 HEAD code lines, final newline, line endings, comments after continuations) was re-run from live files with temp files in scratchpad/B/: results as below.

## codeident.sh (all 44 files of both lists, end of session)

OK 44 / DIFF 0 / NEW 0

## runs.awk, N = 20 (all 44 files of both lists)

0 runs (49 at the start); 54 runs of 15–19 lines left by judgment.
