# Ledger – group D_converters_piesim (converters/, PIESIM/), 11/19/2026

Status: PAUSED on maintainer instruction. All 115 files edited; codeident.sh on all 115: **115 OK, 0 DIFF, 0 NEW**.
Line endings and end-of-file state match HEAD (CRLF files kept CRLF; op2copvar_sop.m and opvar2copvar.m were LF and stay LF).

Conventions applied: license block copied from PIETOOLS_stability_sop.m (Copyright ... through "authorship, and a brief
description of modifications"), with a banner and `% PIETOOLS - NAME` title line; where a file already had the
"If you modify" lines, only the part up to the closing banner was inserted. Added support line
`% For support, contact M. Peet, ...` placed after the license, directly before the Initial coding line (rule 3).
Old files (new = 0): `% MMP, 11/19/2026: Release comment pass (comment-only).` appended at the end of the log.
No "Initial coding DJ - 12/31/2021" line was written (coordinator rule); the one 2021b-import file without an
Initial coding line got `% Initial coding MMP, SS, DJ (PIETOOLS 2021b; in the repository since 12/31/2021)`.

## Edited files – new = 1 (dev notes removed: change logs, per-line markers, (was) lines, Start/End and BEGIN/END markers, commented-out old code)

- PIESIM_1Dto1D2Mat_cheb_2D: formula; inputs prow/pcol (header named px/py), N as 1x2; "Called by" removed; license; support; log removed; step comments.
- PIESIM_1Dto2D2Mat_cheb_2D: formula; inputs prow/pcol (header named p/p1), Rop as cell {R0,R1,R2}; "Called by" removed; license; support; log removed; step comments.
- PIESIM_2Dto0D2Mat_cheb_2D: formula; p as 2 x ns; "Called by" removed; license; support; log removed; step comments.
- PIESIM_2Dto1D2Mat_cheb_2D: formula; "1D to 2D" corrected to 2D to 1D; prow/pcol; license; support; log removed; step comments.
- PIESIM_3PI_op_discretize_2D: formula; nonexistent output A_nonsquare removed from header; license; support; log removed; step comments.
- PIESIM_3PI_op_line_int_discretize_2D: wrong file name in title fixed; formula; license; support; log removed; step comments.
- PIESIM_9PI_opint_area_discretize_2D: formula; A_nonsquare and size "(N+1)x(N+1)" corrected (prod(N+1) square, not truncated); license; support; log removed; step comments.
- PIESIM_9PI_opint_line_discretize_2D: formula; size corrected; license; support; log removed; step comments.
- PIESIM_9PI_opmult_discretize_2D: formula; A_nonsquare and size corrected; "Called by" removed; license; support; log removed; step comments.
- PIESIM_Mat_Truncate: index rule stated; N, p shapes; blank line in help block closed; license; support; step comments.
- PIESIM_3PI_opint_discretize: formula; nonexistent input p removed; license; support; log removed; step comment.
- PIESIM_3PI_opmult_discretize: formula and one-variable restriction; blank line in help closed; license; support; log removed; step comments.
- PIESIM_PI2Mat_opint_cheb: formula; "see Peet & Peet 2021 paper" replaced by the definition of p; license; support; log removed; step comments.
- PIESIM_PI_discretize: description corrected (returns Chebyshev coefficients, not grid values); name-length rule stated; "Used by" removed; license; support; step comments.
- PIESIM_disturbance_check: psize fields and the output fields wspace/uspace/wsep/usep/w_order/u_order documented; license; support; MMP log removed; step comments.
- weightPIE: input structure (struct of kernels, no ODE states) stated; usage sentence removed; license; support; MMP log removed.
- clencurt (third party): formula and INPUT/OUTPUT added; step comments; no license (see candidates).
- convert_PIETOOLS_BVP: references to compose_Tmaps, split_BVPs, decouple_BVPs, compute_Tmap_1D, lift/eliminate_BVP_sources, fund_inputs, convert_PDE_term replaced by statements; supported-BC conditions stated; solve_BVP_group header rewritten without references; support; log removed.
- internal/BVP/compose_Tmaps: inputs maps/source described without references; support; log removed.
- internal/BVP/compute_Tmap_1D: references removed; support; log removed.
- internal/BVP/compute_canonical_BCs: support; log removed.
- internal/BVP/decouple_BVPs: "as initialize requires" reworded; support; log removed.
- internal/BVP/eliminate_BVP_sources: header and local-function inputs described without references; support; log removed.
- internal/BVP/fund_inputs: inputs_info fields stated (was a pointer to compose_Tmaps); support; log removed.
- internal/BVP/lift_BVP_sources: references removed; support; log removed.
- internal/BVP/split_BVPs: references removed; limitation stated; support; log removed.
- internal/compute_PIETOOLS_PDE_periodic_map: signature line, INPUT/OUTPUT added; license; support; usage note replaced by the rank requirement; step comments.
- internal/convert_PDE_term: Tops fields described without references; support; log removed.
- legacy/2D/compute_PIETOOLS_PDE_periodic_map_2D: signature line; info fields stated; references removed; local-function INPUT/OUTPUT; support; log removed.
- legacy/convert_PIETOOLS_PDE_2D: outputs (initialized pie_struct, copvar operators) stated; references and inline DB/Claude markers removed; 1-D periodic rank criterion and no-equation error stated; "formulas below" corrected to above; commented-out code removed; local-function headers; support; log removed.
- legacy/convert_PIETOOLS_PDE_nonlinear_2D: references removed; history wording removed; Start/End markers and commented-out code removed; support; log removed.
- legacy/pie_pde_maps_2D: references removed; pde2pie fields stated; support; log removed.
- operator_converters/dopvar2sdopvar: block mapping and accepted classes stated; commented-out code removed; support.
- operator_converters/dpvar2sdvar: function name in help fixed; P.n formula fixed; "row vector" comment corrected to column; license; support; Initial coding DJ - 08/17/2026 (git); log removed.
- operator_converters/dpvars2sdvars: help signature fixed (named quad_expand); input name vars; NOTES corrected (A, B as outputs); "row vector" corrected (3x); support.
- operator_converters/mat2sopvar: multiplier role and variable restriction stated; support; log removed.
- operator_converters/npie2pie_struct: references removed; support; step comment.
- operator_converters/op2copvar_sop: header rewritten with INPUT/OUTPUT (no legacy framing, no references); license; support; local-function header; step comments; log removed.
- operator_converters/opvar2copvar: references and history removed from NOTES; BEGIN/END markers and orphan commented-out lines removed; local-function I/O; log removed.
- operator_converters/opvar2d2sopvar: dummy-name and opvar handling stated; test/function references removed; support; log removed.
- operator_converters/opvar2sopvar: P and dummy-name restrictions stated; history comments rewritten; support; log removed.
- operator_converters/pie2copvar: header rewritten to what the code does (calls initialize); references removed; body comment without history; DB/DJ log removed.
- operator_converters/pie2npie_struct: references removed; support; log removed.
- operator_converters/pie2piesim: references removed; error conditions stated; support; log removed.
- operator_converters/sdvar2dpvar: license; support; error condition stated; MMP/DJ log (with CLAUDE.md and test references) removed.
- print_convert_summary: OUTPUT stated; reference removed; DB marker removed; support; Initial coding DB - 11/18/2026 (git).

## Edited files – new = 0 (logs and markers kept; release log line appended)

Full pass (header, references, license/support as needed, step comments where reviewed):
- PIESIM_integral_projection1D: formula; "Called by" removed; license; support; step comments.
- PIESIM_discretize_icf_2D: license; support; step comments on the forcing blocks.
- PIESIM_9PI2Mat_cheb_2D: inputs prow and N; license; support.
- PIESIM_NonPoly2Mat_cheb_2D: title name fixed; license; support.
- PIESIM_plot_solution_2D, PIESIM_transform_to_solution_2D, PIESIM_plot_solution, PIESIM_transform_to_solution: blank lines inside the help block closed; title names fixed (2D files); "Requires" references removed; output figs documented (1D plot); license; support.
- PIESIM.m: references to solver_PIESIM.m and the manual removed; license.
- PIESIM_4PI2Mat_cheb, PIESIM_NonPoly2Mat_cheb: paper reference / "Requires" removed; license; support.
- vec_order2elem_order: title name fixed; license; support.
- transl_1D, transl_2D, transl_3PI, transl_mod, rescalePIE, rescalePIE_2D: OUTPUT added, wrong input names (T for PIE) and signatures fixed, help blank line closed (transl_1D).
- solver_PIESIM: references to other files removed; release line placed at the true end of the log; license.
- fcgltran4d (written by YP): license; support.
- convert_PIETOOLS_PDE: references removed in header and body; support.
- initialize_PIETOOLS_PDE: INPUT/OUTPUT added; manual reference removed; title fixed; support.
- legacy/initialize_PIETOOLS_PDE_batch, legacy/initialize_PIETOOLS_PDE_terms_legacy: INPUT/OUTPUT; license; support.
- legacy/convert_PIETOOLS_PDE_terms: INPUT/OUTPUT; manual reference removed; Initial coding (2021b form) added.
- legacy/convert_PIETOOLS_PDE_terms_internal: title fixed; INPUT/OUTPUT; reference removed; help blank lines closed; license; support; Initial coding SS - 12/05/2022 (git).
- legacy/var_swap: help function name fixed (said subs); support.
- internal/poly_decomp_1D_econ: license.

Mechanical pass only (license where missing, support line, release log line):
PIESIM_Poly2Mat_cheb_2D, PIESIM_discretize_all_2D, PIESIM_discretize_domain_2D, PIESIM_discretize_ops_2D,
PIESIM_fullPI2Mat_cheb_2D, PIESIM_3PI2Mat_cheb, PIESIM_Poly2Mat_cheb, PIESIM_discretize_all, PIESIM_discretize_domain,
PIESIM_discretize_icf, PIESIM_discretize_ops, PIESIM_initial_setup, PIESIM_input_check, PIESIM_options_check,
PIESIM_int_bdf, PIESIM_int_symbolic, PIESIM_stability_check, PIESIM_time_integrate, examples_pde_library_PIESIM_1D,
examples_pde_library_PIESIM_2D, fcgltran2d, fcgltran2dperm, fcht, ifcht, zwgll (release line only),
convert_PIETOOLS_DDE, convert_PIETOOLS_DDF, convert_PIETOOLS_NDS, initialize_PIETOOLS_DDE, initialize_PIETOOLS_DDF,
initialize_PIETOOLS_NDS, internal/convert_PIETOOLS_NDS2DDF, internal/minimize_PIETOOLS_DDE2DDF, minimize_PIETOOLS_DDF,
legacy/2D/Compute_Tmap_PDE_2D_Separable, legacy/2D/PIETOOLS_FTC_Expansion_2D, legacy/2D/expand_PIETOOLS_PDE_BCs,
legacy/convert_PIETOOLS_DDE2PIE_legacy, legacy/convert_PIETOOLS_PDE_batch, legacy/convert_PIETOOLS_PDE_batch2terms_legacy
(license added), legacy/convert_PIETOOLS_PDE_terms_legacy (license added).

## Not yet passed (to resume)

1. The 42 "mechanical pass only" files above: rules 1 (header completeness), 5 (step comments) and 6 (references)
   not reviewed. Known reference hits still open: convert_PIETOOLS_DDF (2), convert_PIETOOLS_PDE_batch (12),
   convert_PIETOOLS_PDE_terms_legacy (3), convert_PIETOOLS_DDE2PIE_legacy (2), Compute_Tmap_PDE_2D_Separable (2),
   initialize_PIETOOLS_DDE/DDF, minimize_PIETOOLS_DDF, convert_PIETOOLS_DDE/NDS, convert_PIETOOLS_NDS2DDF,
   minimize_PIETOOLS_DDE2DDF, PIESIM_input_check, examples_pde_library_PIESIM_1D/2D.
2. Rule 5 (step comments) for the remaining old PIESIM files of the full-pass list (PIESIM.m, the postprocessing
   files, transl_2D, rescalePIE/_2D, PIESIM_9PI2Mat_cheb_2D); longest uncommented stretches: PIESIM_input_check
   (49 lines), transl_2D (50), PIESIM_plot_solution_2D (40), PIESIM_transform_to_solution_2D (40), PIESIM_int_bdf (39),
   PIESIM_stability_check (36), PIESIM_plot_solution (35).
3. Large new DJ files (convert_PIETOOLS_BVP, compose_Tmaps, fund_inputs, convert_PIETOOLS_PDE_2D,
   convert_PIETOOLS_PDE_nonlinear_2D): headers, logs and references done; step comments not added where the
   existing comments were judged sufficient (no line-by-line review).

## Skipped

None (no file of the list had uncommitted edits of other sessions).

## Candidates for the maintainer (functional issues seen, not fixed; observed by reading, nothing run)

1. PIESIM_discretize_icf_2D, `uspace` branch: `nu` is used (`coeff.u=zeros(Nforce,nu)`) but never defined in the
   function (only `nw=psize.nw`): a 2-D simulation with control inputs would stop there; the branch also passes
   `uinput.uspace(k)` (a 1x1 cell) where the w branch passes `uinput.wspace{k}`, and has no non-separable fallback.
2. PIESIM_1Dto2D2Mat_cheb_2D, PIESIM_2Dto1D2Mat_cheb_2D: `A_opmult`, `A_opint_block_1/2` are assigned only if R0,
   R1, R2 are nonempty but summed unconditionally (undefined, or stale from the previous block).
3. PIESIM_PI_discretize: a one-variable kernel counts as a function of s when the variable name has at most 2
   characters (else of theta); two-variable kernels rely on the alphabetical order of R.varname for orientation.
4. convert_PIETOOLS_PDE_2D and convert_PIETOOLS_PDE_nonlinear_2D contain a local print_convert_summary
   (PIE,obj,tdiff_list) that duplicates, with another signature, the file converters/print_convert_summary.m.
5. op2copvar_sop passes a 'dopvar' to opvar2copvar, which throws 'opvar2copvar:badInput' unless the object is
   also an 'opvar' (class hierarchy not checked).
6. compute_PIETOOLS_PDE_periodic_map_2D: the author line reads "Initial coding: Claude 10/28/2026" (git first add:
   DanBraghini 10/29/2026); kept as is.
7. Third-party files clencurt, fcgltran2d, fcgltran2dperm, fcht, ifcht, zwgll: the PIETOOLS GPL block was NOT added
   (it would state "Copyright (C) 2026 PIETOOLS Team" on code by Trefethen, von Winckel and unknown authors); no
   support line added. Decision needed. fcgltran4d (written by YP) got the block.
8. PIESIM files now name M. Peet as contact (brief wording); PIESIM.m and solver_PIESIM.m name Y. Peet.
9. dpvars2sdvars: the existing license title line reads "% PIETOOLS - quad_expand" (kept: existing block).
10. initialize_PIETOOLS_PDE_terms_legacy and convert_PIETOOLS_PDE_terms_internal: header notes say "The workspace is
    cleaned at the end of the file" (functions); the latter also keeps a developer note ("Perhaps we should create a
    PIE structure ... I don't know ..."), left as other-author text in an old file.
11. PIESIM_disturbance_check: `nargin=length(varargin)` shadows nargin (harmless).

## Unresolved (stated with limited confidence)

1. mat2sopvar: header says Pmat is the multiplier term in the variables shared by vars.in and vars.out; the
   treatment of input-only variables is not stated.
2. convert_PIETOOLS_BVP, info.fund.p: wording copied from compose_Tmaps ("xf_c stacks the states for each auxiliary
   value if p is the total size of the states"); meaning not verified.
3. PIESIM_9PI_opint_area_discretize_2D: axis order of the 4-D sample array not stated.
4. convert_PIETOOLS_PDE_2D OUTPUTS: "initialized pie_struct with copvar operators" taken from the code comment at
   the initialize call; the field list still names the opvar/opvar2d objects they are built from.

## Literature citations kept

compute_PIETOOLS_PDE_periodic_map (Jagt, Chernyshenko and Peet, arXiv:2503.22896v5, Def. 1, Lemma 5, Prop. 6,
Lemma 7); compute_PIETOOLS_PDE_periodic_map_2D and convert_PIETOOLS_PDE_2D ([Jagt, 2025]); clencurt (Trefethen,
source URL); examples_pde_library_PIESIM_1D examples 35-36 ("Amritam's paper").
