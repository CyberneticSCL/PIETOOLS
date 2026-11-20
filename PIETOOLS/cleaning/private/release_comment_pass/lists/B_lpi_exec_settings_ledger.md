# B_lpi_exec_settings – release comment pass ledger (11/19/2026)

Status: PAUSED on the maintainer's instruction. 62 of 70 files fully passed; the 8 executives listed under "Not yet passed" have only the mechanical first step (per-line markers and (was) lines removed). All 70 files: comment-only edits, codeident.sh OK.

Common to every new file (new = 1) unless stated: the change log after the Initial coding line, per-line markers, (was) lines and orphaned commented-out continuation lines, and BEGIN/END change markers removed; text after the date on the Initial coding line cut. The license blocks, "For support" lines and See also lines kept.

## Edited files
- lpi_programming_sopvar/poscopvar.m: log (66 lines) and per-line markers/(was) lines removed; header: spec section, LPIPROGRAM / LPI_INEQ_SOP pointers, PIE usage sentence and settings-report references removed, psatz_rows described, session cache noted; body: report sections, test and function references removed, step comments added; FIT_LEN header completed.
- lpi_programming_sopvar/getsol_lpivar_sop.m: log and markers removed (Initial coding line cut to the date); references to LPISOLVE, sosgetsol, getequation, rand_sdopvar, sopvar2sdopvar, verify, a spec section and a former local removed; possessives reworded.
- lpi_programming_sopvar/lpi_eq_sop.m: log removed; help and step comments without the names of the routines called, LPIPROGRAM / LPIGETDUAL_SOP pointers and the retired-class history.
- lpi_programming_sopvar/misc/get_lift_degs.m: log and markers removed; header: lift-notes section, POSCOPVAR_DIRECT and heatNd_lpi pointers, the unpublished-note theorem behind Dmax/cap_inc removed (options described by what they set), output fields cap_inc, Dmax_rule added; body: BEGIN/END history block, test/measurement notes (DEMO3, S18, withdrawn corner rows, session times) removed; corner products named face products.
- lpi_programming_sopvar/lpi_ineq_sop.m: log, markers, BEGIN/END history blocks removed; help: routine pointers (GET_LIFT_DEGS, POSCOPVAR_DIRECT, LPI_EQ_SOP, SOSINEQ, LPIPROGRAM, LPISETTINGS_SOP, LPIGETDUAL_SOP, LIFT_COPVAR) replaced by what they compute, the lift/weight rule stated, wcross/dw_space/Dmax/cap_inc/psatz_rows described, info.degrees fields completed; body: lift-notes section, heatNd_lpi, retired LPI_INEQ history removed, step comments added.
- lpi_programming_sopvar/lpidecvar.m (2025 file, DJ): "shell for sosdecvar/sospolymatrixvar" restated as what is declared, "(see also 'lpiprogram')" removed, support line added; release log line appended.
- lpi_programming_sopvar/lpigetsol_sop.m: log removed; help without GETSOL_LPIVAR_SOP / SOSGETSOL / LPISOLVE pointers; body: LPIGETSOL history, orphaned commented-out elseif and the trailing note on deleted locals removed.
- lpi_programming_sopvar/lpiprogram.m (2025 file, DJ; MMP additions): SETTINGS input without LPISETTINGS_SOP / LPI_INEQ_SOP / LPISOLVE / LPISETSETTINGS pointers; domain row "[1,b]" corrected to "[a,b]"; support line added; release log line appended.
- lpi_programming_sopvar/lpiscript_sop.m: log removed; help without the LPISCRIPT comparison and lpisettings_sop pointers (dispatch table of executives kept as the description of what is run); body: LPISCRIPT/lpisettings_sop references removed, a step comment added.
- lpi_programming_sopvar/lpisetsettings.m: help without LPI_INEQ_SOP / LPISOLVE / LPIPROGRAM / LPISETTINGS_SOP / SOSSOLVE pointers, support line added; Initial coding line cut to the date (rationale and maintainer quote removed); step comments added.
- lpi_programming_sopvar/lpisetobj.m (2025 file, DJ): "shell for 'sossetobj'" and "(see also 'lpiprogram')" replaced by what is done; support line added; release log line appended.
- lpi_programming_sopvar/lpisolve.m (2025 file, DJ; MMP additions): help without the sossolve / LPISETTINGS_SOP / LPIPROGRAM / LPISETSETTINGS / SOSSOLVE pointers, the empty-solver rule stated; two marked comments (MMP 11/12, 11/18) without the EXEC_SOP_SETTINGS / DEMO3 / date references (markers kept); a step comment added; release log line appended.
- lpi_programming_sopvar/lpivar_cdopvar.m: log, markers and the trailing note on deleted locals removed; help: state-feedback usage example and LPIPROGRAM pointers removed; body: poscopvar/verify/tensor_index/canonicalize_multiplier/get_lpivar_degs references and former-local notes removed, possessives reworded.
- lpi_programming_sopvar/misc/lpi_classify_sop.m: log removed; lpisolve/sossolve/SOSSOLVE pointers removed from help and body.
- lpi_programming_sopvar/misc/lpi_shape_sop.m: log and the Initial coding remark (test helper) removed; sossolve pointers in help and body restated.
- lpi_programming_sopvar/misc/lpigetdual_sop.m: log removed; LPI_EQ_CDOPVAR / LPI_EQ_SOP pointers removed from INPUT.
- lpi_programming_sopvar/misc/lpi_soseq.m: log removed (Initial coding line cut to the date); getequation/soseq-internals comparisons and the sosprogram/lpiprogram pointer removed; the soseq fallback, which the code calls, still named.
- lpi_programming_sopvar/misc/lpi_eq_sdopvar.m: log (initial entry rewritten as "Initial coding MMP, 08/28/2026") and the BEGIN/END move note removed; lpiprogram pointer removed; one step comment.
- lpi_programming_sopvar/misc/lpi_eq_cdopvar.m: log removed; lpiprogram, collect_eq_rows, impose_eq_rows, lpi_soseq, lpigetdual_sop, sossolve and class-documentation pointers removed from help and body comments; possessives reworded.
- lpi_programming_sopvar/misc/eq_opts_sopvar.m: log, former-local notes and the trailing note on deleted locals removed; POSSOPVAR/possopvar pointers replaced by what the options mean.
- lpi_programming_sopvar/misc/lift_copvar.m: log removed (Initial coding line cut); copquadvar pointers removed; LIFT_BLOCK comment corrected to the kernel the code builds (I_mk kron V^alpha, component outer, as the header states; it said V^alpha kron I_mk).
- lpi_programming_sopvar/misc/poscopvar_lift.m: log removed (Initial coding line cut); lift_copvar/copquadvar pointers removed from comments.
- lpi_programming_sopvar/misc/posmult_cdopvar.m: log removed; copquadvar and @cdopvar/mtimes pointers removed from the pruning comment.
- lpi_programming_sopvar/misc/private/expand_full.m: log removed; Initial coding line added from the date the log gave (MMP, 09/12/2026).
- lpi_programming_sopvar/misc/private/dvar_rows.m: log removed (its initial entry rewritten as the Initial coding line); "as before" and a possessive removed.
- lpi_programming_sopvar/misc/private/merge_monomial_product.m: log removed; initial entry "MP, 08/22/2026" written as "Initial coding MMP, 08/22/2026".
- lpi_programming_sopvar/misc/private/build_exponent_grid.m: log removed; initial entry written as "Initial coding MMP, 08/22/2026".
- lpi_programming_sopvar/misc/private/process_degrees_one.m: log removed; Initial coding line added from the date the log gave (MMP, 08/22/2026); former-private, build_exponent_grid, poslpivar_2d and settings2possopvar pointers removed.
- lpi_programming_sopvar/misc/private/impose_eq_rows.m: log and BEGIN/END markers removed; collect_eq_rows/lpi_soseq pointers restated.
- lpi_programming_sopvar/misc/private/unpack_sheets.m: log and BEGIN/END markers removed; Initial coding line added from the date the log gave (MMP, 08/22/2026); former-private and copquadvar history removed.
- lpi_programming_sopvar/misc/private/collect_eq_rows.m: log and BEGIN/END markers removed (Initial coding line cut); getequation, soseq, lpi_soseq, impose_eq_rows, LPIGETDUAL_SOP and canonicalize_multiplier pointers and the one-soseq-per-parameter history removed; possessives reworded.
- lpi_programming_sopvar/misc/copquadvar.m: log (72 lines) and five BEGIN/END change blocks removed; spec sections (sopvar document Sec. 6.1, 8.4, 9.x), sopquadvar/poslpivar(_2d)/settings2possopvar/sosquadvar/int_semisep/unpack_sheets/plus_batch comparisons, former-local/private notes, the deleted alpha_grid note and benchmark names removed; possessives reworded.
- lpi_programming_sopvar/misc/poscopvar_direct.m: log (80 lines), BEGIN/END blocks and a log note inside TERM_WEIGHT removed; LPIPROGRAM/SOSQUADVAR/SOSGETSOL/POSCOPVAR_LIFT/PROCESS_DEGREES_ONE/copquadvar/lift_copvar pointers, the settings report and proposal references, the unpublished-note citation behind cap_inc and a dated test note removed; "Products of faces" named face products; possessives reworded; step comments added to the DIR_LIST table.
- settings_sopvar/settings_sopvar_{light,heavy,veryheavy,stripped,extreme,custom}.m: logs (76-87 lines each), BEGIN/END blocks and the orphaned "deg_Q line read before" note removed; help rewritten without the settings report, LPISETTINGS_SOP/POSCOPVAR/LPI_INEQ_SOP/LPIVAR_CDOPVAR/EXEC_SOP_SETTINGS pointers and the "stock storage" history, field list completed (Dmax, cap_inc, deg_Q, max_raise_Q, raise_Q); body: section comments without routine names, report notation (d1, d2, d3, Pi, c, F, o) replaced by what each field is (s- and theta-lengths of the letters a, b, c; caps; term sets).
- settings_sopvar/lpisettings_sop.m: log removed (Initial coding line cut: session and README references); help: settings report, POSCOPVAR/SETTINGS_SOPVAR_<TYPE>/EXEC_SOP_SETTINGS/LPI_INEQ_SOP/LPIVAR_CDOPVAR/LPISETTINGS_STOCK_SOP/CUADMM_SEDUMI pointers and the tier history removed; the table notation d1, d2, d3, Pi defined; stale "psatz 'auto'" corrected to the tier value {'auto','faces1','auto'}; Dmax, cap_inc, deg_Q, max_raise_Q, raise_Q added; body comment restated.
- settings_sopvar/lpisettings_stock_sop.m: log removed (Initial coding line cut); help: baseline/stock-executive usage and the "held until 11/12/2026" history removed, the slack N = N1 + N2 stated, LPISETTINGS_SOP/EXEC_TOOLS_SOP/POSLPIVAR_SOP pointers removed, the inline pointer to SETTINGS_2D_STOCK turned into a See also line (added); attribution of the degree levels (Peet, Shivakumar, Jagt 2021-2022) kept; step comments added.
- settings_sopvar/settings_2d_stock.m: log removed, its initial entry written as "Initial coding MMP, 11/17/2026"; help: the "carried until 11/17/2026" history, POSCOPVAR/EXEC_TOOLS_SOP/SETTINGS_2D_SOP/POSLPIVAR_SETTINGS_2D_SOP pointers removed, attribution of the 2-D levels (Peet, Shivakumar, Jagt 2021-2022) kept; body: "stock ... level" labels and routine references restated, one step comment added.
- executives_sopvar/utility_functions/exec_sop_settings.m: log, BEGIN/END block and orphaned (was) lines removed; help: lpisettings_sop/lpisolve/LPI_INEQ_SOP/LPI_SOLVE_LOOP_SOP/PIE_WITNESS_SOP/LPIGETDUAL_SOP/GET_LIFT_DEGS/LPISETTINGS_STOCK_SOP/POSCOPVAR pointers, the settings report and the tier dates removed; the fields deg_Q, max_raise_Q, raise_Q (read but undocumented) added; body comments restated, one step comment added.
- executives_sopvar/utility_functions/exec_info_sop.m: log removed; LPI_SOLVE_LOOP_SOP/LPI_CLASSIFY_SOP/LPI_INEQ_SOP/EXEC_SOP_SETTINGS/LPI_SHAPE_SOP pointers and the "every container executive" usage replaced by what the inputs hold; possessives reworded.
- executives_sopvar/utility_functions/exec_kmult_sop.m: Initial coding line cut to the date (list of executives using it and convert_PIETOOLS_PDE pointer removed).
- executives_sopvar/utility_functions/exec_pie_sop.m: log removed (Initial coding line cut); help without the executives_sopvar/PIE2COPVAR/INITIALIZE pointers; body comments without "as every executive requires".
- executives_sopvar/utility_functions/exec_tools_sop.m: log (103 lines, MMP and DJ entries), BEGIN/END block and a DJ note with commented-out header lines removed; help rewritten without EXEC_PIE_SOP/EXEC_SOP_SETTINGS/POSLPIVAR_SETTINGS(_2D)_SOP/SETTINGS_2D_SOP/POSCOPVAR/LPIVAR_CDOPVAR/LPI_INEQ_SOP/GET_LIFT_DEGS pointers, the settings report and the "stock light 2-D storage" history; max_raise_Q added to the reduced sop fields (the code reduces it); local headers without MAT2COPVAR_SOP/EYE_COPVAR_SOP/SETTINGS2POSSOPVAR/opvar2d references; possessives reworded.
- executives_sopvar/utility_functions/closedLoopPIE.m (2025 file; SS, DJ, DB, MMP): the PIETOOLS GPL block (missing) added from PIETOOLS_stability_sop.m between the help and the DEVELOPER LOGS lines; support line added; the blank line splitting the help (INPUTS/OUTPUTS) made a comment line so help prints the whole block; "(as returned by PIETOOLS_Hinf_OFcontrol_sop)" removed; release log line appended. Log, markers and (was) lines kept.
- executives_sopvar/utility_functions/getController_sop.m (MMP, AT): log (MMP and AT entries) and markers removed, incl. the special marker "AT, 11/17/2026 (was 'inv' default)"; help restated so the default method 'qp' leads (it described K = Z*P^{-1}, the 'inv' method, as what is returned); COPVAR/INV and QP_COPVAR pointers removed (See also kept).
- executives_sopvar/utility_functions/getObserver_sop.m (MMP, AT): as getController_sop: log and markers removed, help restated with the default method 'qp' first, COPVAR/INV and QP_COPVAR pointers removed.
- executives_sopvar/utility_functions/get_lpivar_degs_sop.m: log removed (Initial coding line cut: test-folder origins); LPIVAR_CDOPVAR and retired GET_LPIVAR_DEGS references removed from help, body and the QDEG_ND header.
- executives_sopvar/utility_functions/lpi_solve_loop_sop.m: log and orphaned (was) lines removed; lpisolve/EXEC_SOP_SETTINGS/LPI_CLASSIFY_SOP pointers and a former-local note removed.
- executives_sopvar/utility_functions/pie_disc_sop.m: log removed (Initial coding remark on PIESIM removed); PIE2PIESIM pointers removed.
- executives_sopvar/utility_functions/pie_witness_sop.m: log (MMP and DJ entries) removed; PIE_DISC_SOP pointer, former-local notes and the opvar-test history removed; possessives reworded.
- executives_sopvar/utility_functions/poslpivar_sop.m: log removed (Initial coding line cut: test-folder origins); retired POSLPIVAR history and the POSCOPVAR/poslpivar.m pointers removed.
- executives_sopvar/utility_functions/poslpivar_settings_sop.m: log removed (Initial coding line cut); help: the operator built stated in place of the POSLPIVAR_SOP call, LPISETTINGS_SOP/LPISETTINGS_STOCK_SOP pointers and the retired POSLPIVAR history removed.
- executives_sopvar/utility_functions/poslpivar_settings_2d_sop.m: log removed; help and local headers without POSCOPVAR/SETTINGS2POSSOPVAR/POSLPIVAR_2D/get_eq_opts_2D/settings_PIETOOLS_*_2D/LPISETTINGS references (the operations described as translations); possessives reworded.
- executives_sopvar/utility_functions/qp_copvar.m (AT; MMP, DJ log entries): log removed; support line added and a blank line set between help and license so help stops at the help block; apply_sopvar pointers and "as before"/"old product"/former-local notes removed; seven step comments added (option checks, assembly, solve).
- executives_sopvar/utility_functions/settings2possopvar.m: log removed (Initial coding line cut); help and comments describe the translation without the possopvar/poslpivar/poslpivar_2d/build_monoms/sopquadvar references, the retired settings_PIETOOLS_* history and the test reference; See also kept.
- executives_sopvar/utility_functions/settings_2d_sop.m: log removed; lpisettings_sop pointer, the stock extract_psatz_deg/extract_psatz_opts reference and a note on the dual gain executive removed.
- executives_sopvar/utility_functions/trace_rn_sop.m: log removed (Initial coding line cut: test-folder origin).
- executives_sopvar/PIETOOLS_stability_sop.m: log (118 lines), markers, (was) lines and BEGIN/END blocks removed; help: the opvar/opvar2d PIE input and the retired stock-format settings fields (dd1, dd12, options1, options12, override1, which EXEC_SOP_SETTINGS now rejects) removed, lpisettings_sop/POSCOPVAR/SETTINGS_2D_STOCK/EXEC_PIE_SOP/EXEC_SOP_SETTINGS/EXEC_INFO_SOP pointers replaced by what they hold; body: report sections, test cases (S3, S17, S18) and counts, "header 11/17/2026" pointers and "No storage" wording removed (the Q form named the generalized storage T'*Q).

## Error and warning strings that reference files, functions, dates or reports (code, not changed)
- error string: getsol_lpivar_sop:badBlockZd message "(see cdopvar/verify)"; notSolved message names "lpisolve".
- error string: "use 'lpiprogram' to initialize your program".
- error string: "either 'lpisolve' has not been run".
- error string: "a settings structure similar to the output of lpisettings_sop()".
- error string: "Settings should be a 'struct' (LPISETTINGS_SOP) or the name of a tier."
- error string: "'prog' should be an LPI program structure; use 'lpiprogram'."
- error strings: lpigetdual_sop:tag "was the equality written by lpi_eq_sop?"; lpigetdual_sop:y "(solve it first; sossolve keeps y)".
- error string: "First input must be an LPI program structure; see 'lpiprogram'."
- error strings: "A 'subset' cap grades the basis of 'copquadvar'... Use 'poscopvar'." and the 'joint' message (copquadvar, poscopvar).
- error string: "unpack_sheets takes the triplets of int_semisep(...,'triplets') ...".
- error strings: "see 'lpiprogram'"; "declare it with 'lpidecvar' before imposing constraints".
- error string: "'psatz' should be 0, 1, or a face code ... 2 is poslpivar_2d's ball, not the box."
- error strings: exec_sop_settings:st "(lpisettings)"; :retired "settings.%s was retired on 11/12/2026" with "use sop.stock (LPISETTINGS_STOCK_SOP)"; the settings_nd message "retired on 11/12/2026"; :row "(settings report 11/16/2026, Sec. 5)".
- error strings: exec_tools_sop:slack messages name settings_2d_stock(level) and "the stock 1-D and 2-D settings".
- error strings: poslpivar_sop:dim "poslpivar degrees describe one spatial variable per space"; poslpivar_sop:excludeRn "poscopvar needs one basis operator per space".
- error string: poslpivar_settings_sop:stock "take the settings from LPISETTINGS_STOCK_SOP".
- error strings: :sosineq "use_sosineq selects lpi_ineq_2d, which has no container counterpart"; :exclude1 and :exclude1d "poscopvar needs one basis operator per space".
- report strings (sprintf into out.report.lines/lossy, code): name possopvar and poslpivar_2d.
- error string: settings_2d_sop:st "an lpisettings struct (field settings_2d)".
- lpi_programming_sopvar/poscopvar.m, misc/copquadvar.m: error "'psatz' should be 0, 1, or a face code ... 2 is poslpivar_2d's ball, not the box."

## Not yet passed (PAUSE)
Markers and (was) lines already removed (codeident OK); still to do: change log, BEGIN/END blocks, references (report, executives, routines, tests), header and step comments.
- executives_sopvar/PIETOOLS_H2_control_sop.m
- executives_sopvar/PIETOOLS_H2_estimator_sop.m
- executives_sopvar/PIETOOLS_H2_norm_sop.m
- executives_sopvar/PIETOOLS_Hinf_OFcontrol_sop.m
- executives_sopvar/PIETOOLS_Hinf_control_sop.m
- executives_sopvar/PIETOOLS_Hinf_estimator_sop.m
- executives_sopvar/PIETOOLS_Hinf_gain_sop.m
- executives_sopvar/PIETOOLS_well_posedness_sop.m

## Skipped
- none: no file of the list had uncommitted edits of another session when its first edit was made.

## Candidates for the maintainer (functional, not fixed)
- lpi_classify_sop.m: the if/else that reads back OBJ has two identical branches (cl.obj = double(lpigetsol_sop(prog,obj)) for 'dpvar'/'polynomial' and for any other OBJ).
- settings_sopvar_custom.m: help says "as shipped, the light tier" and lpisettings_sop says "custom as light", but deg_Q{1} is [2 2 2] in custom and [4 2 2] in light (light changed 11/18). Comments now state the difference; the values were not touched.
- settings_2d_stock.m: use_sosineq is overridden from the base workspace (sosineq_on / use_sosineq via evalin); with 1 the struct selects an SOS inequality that the executives reject, so a stray base-workspace variable changes the returned settings.
- pie_disc_sop.m, local ENTRY: the 'double' branch and the fall-through both return X(i,k) (redundant branch).
- poscopvar.m: a program without settings falls back to the mult/int/joint rule (struct ps), which the old comment called 'light'; the 'light' tier is now a storage row (len_com).
- PIETOOLS_stability_sop.m (help, now corrected): it listed the stock-format settings fields dd1, dd12, options1, options12, override1, which EXEC_SOP_SETTINGS rejects as retired.
- lpisettings_sop.m (help, now corrected): it gave psatz 'auto' for the tiers; every tier sets {'auto','faces1','auto'}.
- getController_sop.m, getObserver_sop.m (help, now corrected): the description led with the 'inv' formula K = Z*P^{-1} though the default method is 'qp'.

## Unresolved
- lpi_ineq_sop.m: the weight rule in the help is stated for wR = 0; an undocumented opts.wR still reaches GET_LIFT_DEGS.
- get_lift_degs.m, poscopvar_direct.m: the justification of Dmax / cap_inc (an unpublished note, Shivakumar 11/17/2026, Thm 1.4, Lemma 1.2) was removed with the reference; the options are described by what they set. Restore as a citation if the note is published.

## Added text of note
- closedLoopPIE.m: PIETOOLS GPL block added (copied from PIETOOLS_stability_sop.m) and the help made contiguous; release log line appended (2025 file).
- lpisettings_stock_sop.m: the inline pointer to SETTINGS_2D_STOCK turned into an added See also line.
- Release log line "% MMP, 11/19/2026: Release comment pass (comment-only)." appended in the 2025 files lpidecvar.m, lpiprogram.m, lpisetobj.m, lpisolve.m, closedLoopPIE.m; support lines added where missing (lpidecvar, lpiprogram, lpisetobj, lpisetsettings, closedLoopPIE, qp_copvar).
- Initial coding lines written from dates the logs gave (no git-import dates used): expand_full (MMP, 09/12/2026), unpack_sheets, process_degrees_one, merge_monomial_product, build_exponent_grid (MMP, 08/22/2026), dvar_rows (09/22/2026), lpi_eq_sdopvar (08/28/2026), settings_2d_stock (11/17/2026).

## Literature or author attributions kept
- lpisettings_stock_sop.m, settings_2d_stock.m: degree levels of M. Peet, S. Shivakumar and D. Jagt (2021-2022).

## codeident.sh (all 70 files of the list, end of session)
OK 70 / DIFF 0 / NEW 0
