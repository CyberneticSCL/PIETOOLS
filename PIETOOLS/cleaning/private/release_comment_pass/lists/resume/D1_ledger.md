# Ledger – D1 (converters/), resume 11/20/2026

Scope: the 53 files of `lists/resume/D1.tsv`. Before the first edit: codeident 53 OK (no file changed by another
session). Edited this session: 39 files (16 mechanical-only, 23 passed). Not edited: 14 passed files without runs or
open items (compute_canonical_BCs, poly_decomp_1D_econ, legacy/convert_PIETOOLS_PDE_terms, var_swap, dopvar2sdopvar,
dpvar2sdvar, dpvars2sdvars, npie2pie_struct, op2copvar_sop, opvar2copvar, opvar2sopvar, pie2copvar, pie2npie_struct,
sdvar2dpvar). Interrupted once by the usage limit (04:00): codeident on all 53 then OK, nothing to repair.

Final: **codeident 53 OK / 0 DIFF / 0 NEW**; **runs.awk (N = 20): 0** (100 at the start); runs of 15-19 lines: 95
(by judgment; the ones met while working were broken). Line endings and final-newline state as before (CRLF files
CRLF; op2copvar_sop.m, opvar2copvar.m LF), checked against the saved state and HEAD. All new text ASCII.

## Edited – mechanical-only (16; new = 0: change logs, markers, (was) lines and commented-out code kept)

- convert_PIETOOLS_DDE: references to minimize_PIETOOLS_DDE2DDF / convert_PIETOOLS_DDF replaced by the construction;
  PIE formula corrected (w and u swapped, stray "+ +"); "type" -> out_type; INPUT (DDE formula, fields) / OUTPUT.
- convert_PIETOOLS_DDF: references to solver_PIETOOLS_DDF / initialize_PIETOOLS_DDF removed; PIE formula corrected;
  PIE state [x; phi_s] stated; INPUT (DDF formula; Cvd in s1) / OUTPUT; comments "number of states" (nK) and
  "number of inputs" (nv) corrected, "DDE" -> DDF; step comments with the formulas (Chv, DI, CI, Tb0, Itau, Tba, Tbb).
- convert_PIETOOLS_NDS: references removed; PIE formula corrected; INPUT (NDS formula) / OUTPUT (DDF has no Cvd).
- initialize_PIETOOLS_DDE: header with the DDE formula and the dimension rule, INPUT/OUTPUT ("as in PIETOOLS_DDE"
  removed); "less than nK" -> "at most nK"; 49 step comments; stray comment after `end` ("Now for the hard part -
  setting up the PIE") removed.
- initialize_PIETOOLS_DDF: header with the DDF formula; references to convert_PIETOOLS_DDF removed; "length 1"
  corrected to 0; INPUT/OUTPUT; "at most nK"; an aspiration comment ("make sure no states are used") replaced by the
  fact; 28 step comments (size checks stated as erroring only if both dimensions differ); stray trailing comment removed.
- initialize_PIETOOLS_NDS: header ("as described at the end of this file" was false) with the NDS formula, A0
  required, INPUT/OUTPUT; unreachable branches stated; 54 step comments; stray trailing comment removed.
- internal/convert_PIETOOLS_NDS2DDF: header rewritten (irrelevant PIE text and the convert_PIETOOLS_DDE comparison
  removed): channel construction, Cvd not stored, INPUT/OUTPUT; D21 -> D22 typo; TODO and references replaced.
- internal/minimize_PIETOOLS_DDE2DDF: title name fixed; construction (rank factorization) replaces the stale
  comparison; INPUT/OUTPUT; "number of states" (nK) corrected; D22 typo; TODO and "This makes no sense..." replaced
  by the fact (no delays: recursion); step comments.
- minimize_PIETOOLS_DDF: header likewise (r_i includes v); INPUT/OUTPUT; initialize_PIETOOLS_DDE reference removed;
  nv comment, D22 typo, TODO replaced.
- legacy/2D/Compute_Tmap_PDE_2D_Separable: references (@pde_struct/initialize, reorder_comps) replaced; formula typo
  Twop*uf -> Tuop*uf (2x); INPUT/OUTPUT for construct_1D_Tops and check_state_diffs; compute_input_diffs INPUTS
  corrected (Dvals_v, size_x; Dvals_vy); get_opdim note (nvars = 2); print_reorder_diff_summary signature line fixed;
  comments in the two Twop branch chains.
- legacy/2D/PIETOOLS_FTC_Expansion_2D: reference removed (input must be initialized); opvar output for at most one
  variable stated; step comments.
- legacy/2D/expand_PIETOOLS_PDE_BCs: header was a copy of PIETOOLS_FTC_Expansion_2D (outputs Pop_f2x, Pop_b2x,
  bc_state_tab, which the function does not return); rewritten to what it does (FTC expansion of BC terms, split into
  plane/edge/corner BCs, sorting) with INPUT/OUTPUT; NOTES kept.
- legacy/convert_PIETOOLS_DDE2PIE_legacy: title fixed; references removed (incl. "as convert_PIETOOLS_DDF does" on
  a marked line, marker kept at column 77); PIE formula corrected; INPUT/OUTPUT; step comments with the formula.
- legacy/convert_PIETOOLS_PDE_batch: references removed (solver_PIETOOLS_PDE; "as PI4_COPVAR_SOP names it"; "as the
  pde_struct converter" on a marked line; the pi4 handle described instead of "= PI4_COPVAR_SOP"); PIE formula
  corrected, 11 -> 12 operators; INPUT/OUTPUT; release log line moved from the body to the end of the change log;
  "undefined operators to empty opbjects" -> parameters set to zero; step comments (T, Q, K, core values).
- legacy/convert_PIETOOLS_PDE_batch2terms_legacy: title fixed; reference to the two initializers replaced by the
  conversion summary and INPUT/OUTPUT (added above the first code line; the field list stays below it, code not moved).
- legacy/convert_PIETOOLS_PDE_terms_legacy: title fixed; "the new version" corrected; blank line in the help closed;
  solver_PIETOOLS_PDE reference removed; 11 -> 12; two false notes removed (workspace "cleaned at the end": the clear
  lines are commented out; "Perhaps we should create a PIE structure ... I don't know": it returns a pie_struct);
  INPUT/OUTPUT incl. the 'char' outputs; references in the pi4 header and a marked line removed; step comments.

## Edited – passed (23)

Runs broken (step comments, inserted with a placement check): convert_PIETOOLS_BVP (9 runs), convert_PIETOOLS_PDE (2),
initialize_PIETOOLS_PDE, BVP/compose_Tmaps (3), BVP/compute_Tmap_1D (2), BVP/decouple_BVPs, BVP/eliminate_BVP_sources
(3), BVP/fund_inputs, BVP/lift_BVP_sources, BVP/split_BVPs (2), compute_PIETOOLS_PDE_periodic_map, convert_PDE_term
(2), legacy/2D/compute_PIETOOLS_PDE_periodic_map_2D (3), legacy/convert_PIETOOLS_PDE_2D (7),
legacy/convert_PIETOOLS_PDE_nonlinear_2D (4), legacy/convert_PIETOOLS_PDE_terms_internal (2),
legacy/initialize_PIETOOLS_PDE_batch (6), legacy/initialize_PIETOOLS_PDE_terms_legacy (8), legacy/pie_pde_maps_2D,
opvar2d2sopvar, pie2piesim, print_convert_summary. Also:
- compute_PIETOOLS_PDE_periodic_map_2D: `% Initial coding DB - 10/29/2026` (maintainer decision 11/20).
- convert_PIETOOLS_PDE_2D: OUTPUTS rewritten (Unresolved 3); NOTES formula corrected (T*xf, {T, Tu, Tw}, time
  derivatives, field names); 'Top' flag output stated; blank line in the help block closed.
- convert_PIETOOLS_BVP: info.fund.p reworded (Unresolved 2). mat2sopvar: input-only variables stated (Unresolved 1).
- Pointers to functions removed in passed files: eliminate_BVP_sources (2x "see source_input_map"), fund_inputs ("see
  choose_pieces"; "see fund_inputs" in two local headers, replaced by the fields), convert_PDE_term ("see
  kernel_operator", "as in the input Tops of the main function").
- initialize_PIETOOLS_PDE_batch, initialize_PIETOOLS_PDE_terms_legacy: INPUT/OUTPUT for the local isdefined.

## Skipped

None.

## Candidates for the maintainer (functional issues seen by reading; nothing run, nothing fixed)

1. internal/convert_PIETOOLS_NDS2DDF: the distributed-delay coefficients (Adi, B1di, B2di, Edi, C1di, ..., C2dei) go
   to a local `Cvd{i}`, never to `DDF.Cvd{i}`; the DDF (and PIE) of an NDS drop every distributed-delay term
   (initialize_PIETOOLS_DDF then sets Cvd to zero). `DDF_out=DDF;` is unused.
2. minimize_PIETOOLS_DDE2DDF and convert_PIETOOLS_NDS2DDF with no delays call convert_PIETOOLS_DDE, which calls
   minimize_PIETOOLS_DDE2DDF again: unbounded recursion (DDF_out never assigned).
3. initialize_PIETOOLS_NDS: A0 is never defaulted or checked (NDS_out.A0 = A0 fails without A0); the Ei/Edi,
   C1ei/C1dei, C2ei/C2dei branches of the dimension searches are unreachable, and the C2ei/C2dei ones assign nz, not
   ny; the numel errors for C1ei..C2dei name "Ei".
4. initialize_PIETOOLS_DDF: `nx=size(B1,1,'var')` errors (A0 undefined, B1 defined); the checks
   `if size(X)~=[r,c]` (Cr, Br1, Br2, Drv, Cv, Cvd, C1, ..., D2v) error only if both dimensions differ.
5. convert_PIETOOLS_DDF and convert_PIETOOLS_DDE2PIE_legacy integrate the kernels in s1 (since 12/07/2024), while
   PIETOOLS_DDE.m and PIETOOLS_DDF.m document kernels in pvar s; a kernel in s is not integrated.
6. Compute_Tmap_PDE_2D_Separable: compute_input_diffs indexes Dop_x/Dop_y (one entry per state variable) with
   tmp_log (one flag per row block of Twop, 4 x 1); print_reorder_diff_summary uses the undefined `varname_LHS_t`
   (component without variables), overwrites Ldiff for a second variable, and has `ceil(log(dval/log(10)))`.
7. expand_PIETOOLS_PDE_BCs: in the integral-removal loop, Rindx is looked up in PDE.([Robj,'_tab']) with Robj left
   from the last term processed ('w' or 'u' possible) for a state term.
8. initialize_PIETOOLS_PDE_terms_legacy: nr read as size(PDE.PDE.Crp,1) / size(Drb,1) (cells of terms: number of
   cells); nx from ODE.A only when one of nx, Bxu, ..., Cv is defined (then it overrides an explicit nx).
9. The legacy 2-D converters (convert_PIETOOLS_PDE_2D, _nonlinear_2D, Compute_Tmap_PDE_2D_Separable,
   PIETOOLS_FTC_Expansion_2D, compute_PIETOOLS_PDE_periodic_map_2D) build opvar/opvar2d objects; no @opvar or
   @opvar2d class is in the repository.
10. convert_PIETOOLS_PDE_batch2terms_legacy: the function line names convert_PIETOOLS_PDE_batch2terms (file name
    differs; MATLAB uses the file name).
11. Compute_Tmap_PDE_2D_Separable (new = 0) keeps a per-line marker `% Claude, 10/28/2026` – same authorship question
    as the periodic_map_2D line decided on 11/20; not changed.

## Unresolved (11/19 items settled; remaining confidence notes)

1. mat2sopvar – settled from the sopvar definition: input-only variables S1 are integrated over their domain with a
   kernel constant in them (ZR of degree 0): (Pop*x)(S2,S3) = Pmat(S2,S3)*int_{S1} x(S3,S1) dS1.
2. convert_PIETOOLS_BVP info.fund.p – settled from compose_Tmaps: p = ns (each block stacks the states) or the number
   of auxiliary values of the non-scalar rank-deficient direction; the block count (product over the scalar
   rank-deficient directions outside the subset) is stated in compressed form.
3. convert_PIETOOLS_PDE_2D OUTPUTS – settled: the initialized pie_struct fields T, Tw, Tu, A, B1, B2, C1, D11, D12,
   C2, D21, D22 as 'copvar' operators, space and the dependent dim/dom/vars (old field names Tx, Bw, Cz, ... removed).
4. Left as written: a duplicated license fragment in the body of minimize_PIETOOLS_DDE2DDF; "(B_T in paper)" in
   convert_PIETOOLS_PDE_terms_legacy (unattributed).

## Decisions applied

- New = 0 files: change-log entries and BEGIN/END change notes (with their marked explanation lines) kept verbatim
  under rule 4 although they name functions (PI4_COPVAR_SOP, OP2COPVAR_SOP, INITIALIZE) or describe opvar; plain
  step comments and headers cleaned under rule 6. Commented-out code kept.
- Literature citations kept: arXiv:1910.03881 (DDE/DDF/NDS converters), [Jagt, 2025] (convert_PIETOOLS_PDE_2D).
- In-file pointers ("see NOTES", "listed at the end of this file") kept.
