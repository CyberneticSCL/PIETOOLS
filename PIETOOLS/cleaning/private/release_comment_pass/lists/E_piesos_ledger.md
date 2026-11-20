# Ledger – group E_piesos (PIESOS/, 236 files)

Status: complete for headers, license, authors, developer notes and pointers; step comments partly (see Unresolved).

## Notes on the file list

- PIESOS/ is absent from origin/PIETOOLS-2025 (no path under PIESOS, polyopvar, quadPoly, intop or
  tensop* in that tree). The new = 0 flag of 105 files is name-based (blkdiag, ctranspose, ... exist in
  other classes in 2025), so these files are new by path. Their change logs and markers were kept, as the
  flag directs; a maintainer decision whether to strip them too.

## Edited files

Common to every edited file: contact line added where missing (before the `Initial coding` line); new = 1:
change log (except the `Initial coding` entry), per-line markers, `(was)` lines removed; new = 0: release
log line appended. Below: what else changed.

### PIESOS/ (root, 18 files, all new = 1, DJ)
- Liediff.m – pointers to piesos_quadvar/innerprod/convert_PIETOOLS_PDE/piesos_LF/subs removed; log facts
  restated (adjoint rejected; repeated PDE names get suffix _i).
- dmonomials.m – step comment on the input check.
- innerprod.m – pointers (lift_output/common_outputs, ctranspose/mtimes, int, set.C) removed; step comments.
- monopmat.m – comment-line markers removed; INPUT/OUTPUT for local process_degs.
- piesos_LF.m – pointers removed (piesos_program, convert_PIETOOLS_PDE, piesos_quadvar, piesos_LF_alt,
  V3_Liediff/Examples); 'degs' option described; repeated PDE names (_i) restated; TODO question removed.
- piesos_eq.m, piesos_lpivar.m – piesos_program pointer removed; step comments.
- piesos_getsol.m, piesos_setobj.m, piesos_solve.m, piesos_program.m – "see also" pointers replaced
  (SOSTOOLS program format).
- piesos_ineq.m – pointers removed (piesos_program, innerprod, piesos_poslpivar, piesos_quadvar,
  constant2polyopvar, monopmat); psatz option described; INPUT/OUTPUT for 4 local functions.
- piesos_intopvar.m – pointers removed; "former local SELECT_DEGS" history wording removed; INPUT/OUTPUT
  for 6 local functions; step comment in build_intmat.
- piesos_polyvar.m – intmat/degs options described instead of pointing to piesos_intopvar.
- piesos_poslpivar.m, piesos_sosvar.m – monopmat/piesos_quadopvar/piesos_poslpivar pointers removed.
- piesos_quadopvar.m – exclude/sep/otype described instead of "see monopmat"; pointers (int, sopvar ZL,
  monopmat) removed; comment-line markers removed; INPUT/OUTPUT for split_LR, declare_Zs.
- piesos_quadvar.m – pointers to piesos_quadopvar/monopmat/dmonomials/set.C/dpvar subsref removed;
  L/R options described; INPUT/OUTPUT for 3 local functions.

### PIESOS/custom (6 files, new = 1, CR and DJ)
- DMB.m – incomplete license block completed (warranty, address and "If you modify" lines were missing);
  "Def. 4 in Paper1" references removed; commented-out varmat line removed; contact line; step comment.
- SOS_DP.m – "Def. 9 CDC paper" and piesos_sosvar/piesos_quadvar pointers removed.
- filter_Pcell.m – V3_Liediff pointer and piesos_getsol usage example removed; step comments.
- kernel_degree.m – header states the formula instead of "opdeg in SOS_DP"; SOS_DP pointers and
  "former local MAX_KERNEL_DEGREE" removed.
- polyopvar_equal.m – piesos_getsol usage note, clean_poly pointer and a malformed marker ("% DJ), ...")
  removed; step comments.
- weighted_sobolev_norm.m – 'diff' pointer removed; step comments.

### PIESOS/demos_and_examples (4 files, new = 1, DJ)
- PIESOS_DEMO1/2/3 – INPUT (none) and OUTPUT (workspace) blocks added; DEMO3 header no longer claims
  that the script computes the exact minimum from the BVP (it prints a fixed value for r = 12).
- examples_PDE_library_PIESOS.m – license block added; "Inputs:" renamed INPUT, OUTPUT block added.
  Literature citations [1]-[5] (Hale & Verduyn Lunel 1993, Chafee & Infante 1974, Liu & Krstic 2001,
  Skellam 1951, van den Berg & Jaquette 2018) attribute examples to their sources and were kept.

### PIESOS/executives (4 files, new = 1, DJ)
- PIESOS_settings.m – pointers (piesos_LF, piesos_sosvar, piesos_solve) removed; step comment.
- PIESOS_stability_global.m, PIESOS_stability_local.m – "as output by convert", "see PIESOS_settings"
  and weighted_sobolev_norm pointers replaced by the settings fields used; a dated "(11/2026)" remark
  removed.
- PIESOS_sublevel.m – usage pointers to PIESOS_stability_local, convert, PIESOS_settings removed.

### PIESOS/internal (8 files, new = 1, DJ)
- clean_dpvar.m – step comments.
- coeff_error.m – piesos_eq/piesos_getsol pointers removed; INPUT/OUTPUT for 2 local functions.
- is_eye.m – 'ispvar' pointer removed.
- otimes_params.m – INPUT/OUTPUT for 7 local functions.
- common_vars_quadPolyCell.m, find_common_idcs.m, quad_expand.m, set_vars_quadPolyCell.m – contact
  line only.

### PIESOS/polyopvar/*.m (9 files; int_simple.m new = 0, others new = 1)
- int_simple.m (new = 0, DJ/MP/SS) – first help line named `int(E,var,L,U)`: now the actual call
  INT_SIMPLE(E,IVAR,X1,X2); release log line appended (contact line already present).
- copvar2tensopmat.m, tensopmat2copvar.m – step comment on the input check ("See also" kept).
- tensopvar2polynomial.m – apply_op, 'tensopmat' and otimes_native pointers removed; INPUT/OUTPUT
  for 3 local functions; step comments in convert_op.
- constant2polyopvar.m, mat2tensopmat.m, mat2tensopvar.m, pdevar2polyopvar.m, sopvar2tensopvar.m –
  contact line and log removal only.

### PIESOS/polyopvar/@intop (44 files; 26 new = 1, 18 new = 0)
- All: contact line; new = 1 logs/markers removed (mtimes_tensopvar, otimes, apply_op, combine_terms,
  intop, kernels); new = 0 release log line appended.
- private/Zdims, Zmap, Zpermmap, enumerate_valid_orders, get_int_order, kernels2params, lmult_params,
  prune_basis, remap_cols, remap_dvars, sync_basis – license block added.
- Pointers removed: apply_op ('polyopvar', get_tensorvars), combine_terms (canonicalize_intmat,
  int_order, unique_intmat), ctranspose/transpose ('polyopvar/ctranspose'), display/get_display
  (get_display, display_term, 'intop', is_id_op), intop (get_tensorvars, canonicalize_intmat), isvalid,
  kernels and kernels2params ("documentation of the class", get_tensorvars), mtimes
  (mtimes_tensopvar), mtimes_tensopvar ('sopvar', otimes_params.m, "former local" wording),
  enumerate_valid_orders ("Used by switch_ints", switch_ints-relative input descriptions rewritten,
  "the original" history clauses removed).
- INPUT/OUTPUT added: intop get.* methods (8), apply_op (3), canonicalize_intmat (2), combine_terms,
  isvalid, unique_intmat, mtimes_tensopvar (6), enumerate_valid_orders/backtrack.

### PIESOS/polyopvar/@polyopvar (38 files; 16 new = 1, 20 new = 0 – incl. private/factor_key)
- All: contact line; logs/markers removed (new = 1) or release log line appended (new = 0).
- Pointers removed from descriptions and step comments ("see 'polyopvar'", "see set.C", "see
  'ctranspose'", "see @tensopmat/...", "see num_factors/factor_key/subs_const/clean_poly/horzcat",
  "as in plus/constant2polyopvar", "see there") in about 30 files; error-message strings unchanged.
- polyopvar.m – header added to the constructor (calling forms, INPUT, OUTPUT); INPUT/OUTPUT for
  get.matdim, get.vartab, get.degmat.
- INPUT/OUTPUT added to local functions: combine_factors (3), common_basis, initialize, otimes, subs (5,
  and OUTPUT for subs_const).

### PIESOS/polyopvar/@tensopmat (38 files; 15 new = 1) and @tensopvar (37 files; 13 new = 1)
- All: contact line; logs/markers removed (new = 1) or release log line appended (new = 0).
- Pointers removed (display/get_display, initialize, is_id_op, lift_output, mtimes, merge_out_spaces,
  tensopmat constructor, ctranspose/transpose, num_factors, opower, subsasgn, diff, otimes_op note on
  intop/otimes, constant2polyopvar, parse_diff_inputs).
- INPUT/OUTPUT added to local functions: tensopmat common_vars, initialize, mtimes (2), otimes (3),
  get.matdim; tensopvar apply_op/split_inputs, common_vars, diff (2), mtimes_nested, get.matdim,
  get.deps.

### PIESOS/polyopvar/@quadPoly (40 listed files; 13 new = 1; SS, AT, DJ)
- License block, contact line and `Initial coding` line added to the 28 files without them (authors
  and dates from git: SS 01/16–02/15/2026, AT 02/18/2026; numArgumentsFromSubscript.m SS 06/14/2022,
  the git --follow origin of its content); release log line appended to the 18 new = 0 ones.
- disp.m (new = 0) – its MMP log entry moved from the help block to the end of the header log.
- quadPoly.m – DJ log entry removed from the class header; constructor header with INPUT/OUTPUT
  (states the one-input conversion of a matrix to a constant quadPoly); header for static `create`.
- unionBasis, quadPoly2polynomial, randquadPoly (new = 1) – help-block MMP log entries, markers,
  (was) lines and "former local" wording removed.
- INPUT/OUTPUT added to main functions without any: numel, end, isempty, uminus, uplus, minus,
  transpose, ctranspose, combine, disp, subsref, var_swap, tensExps, liftIndex, quadPoly2polynomial,
  numArgumentsFromSubscript; reshape.m – wrong Input/Output lines (copied from trace) replaced.
- subsasgn_1.m – header added (had none).
- Headers (what/INPUT/OUTPUT) added to the 12 local functions without any: combine (lin2subv,
  subv2lin, basisSize), disp (format, monomStr, termStr, joinTerms), basisSize in
  polynomial2quadPoly, quadPoly2polynomial, var_swap, and subs/isa_char_or_string.
- Pointers fixed: Sym2quadPoly usage line named `quadPoly_fromSym`; blkdiag "as in vertcat/horzcat";
  quadPoly2polynomial "consistent with quadPoly_eval" (no such function).
- hessian.m (AT, new = 1) – already compliant; not edited.

## Skipped

- None of the listed files had uncommitted edits when its folder was started.
- Not in the list (created after the census), not edited: PIESOS/internal/get_tensorvars.m,
  PIESOS/polyopvar/@quadPoly/times.m. (get_tensorvars.m received the contact line by a folder-wide
  run; the line was removed again and git status shows the file clean.)

## Candidates for the maintainer

- PIESOS/custom/filter_Pcell.m: the second output reads `entry.coefficient`, which errors when
  Pcell{1,1} is numeric, although numeric blocks are documented as supported.
- PIESOS/custom/kernel_degree.m: license block says GPL "version 3"; every other file says version 2
  (block kept as is).
- PIESOS/custom/DMB.m: error message text still cites "Def. 4" (code string, not changed).
- PIESOS/piesos_LF.m: removed TODO comment asked whether piesos_polyvar would be better than
  piesos_quadvar for the non-symmetric (PIE-state) case.
- PIESOS/internal/otimes_params.m, otimes_quadPolyCell: the input check for B tests
  `~isa(B,'quadPoly') && ~isa(A,'struct')` (A instead of B in the second test).
- PIESOS/polyopvar/@quadPoly/horzcat.m and vertcat.m: `elseif isnumeric(X)` uses an undefined X
  (should test varargin{k}); a numeric input errors instead of being converted.
- PIESOS/polyopvar/@quadPoly/reshape.m: the header promises RESHAPE(A,[M,N]); the code only takes
  scalar M, N. Commented-out test code at the end of the file kept (new = 0).
- PIESOS/polyopvar/@quadPoly/subsasgn_1.m: function named `subsasgn` in a file subsasgn_1.m, so it
  is never the class subsasgn (dead file; field assignment of quadPoly stays possible).
- PIESOS/polyopvar/@quadPoly/quadPoly.m: `if ~presorted && 0` never runs (dead branch; the
  presorted input has no effect).
- PIESOS/polyopvar/@quadPoly/plus.m: comment "should add error checking (e.g. matching matrix
  dimensions)" – no dimension check (new = 0, kept).
- PIESOS/executives/PIESOS_stability_local.m and PIESOS_sublevel.m: notes on degrees observed to be
  infeasible kept as limitation notes (test-derived; maintainer may want them removed).
- Initial coding lines that name other files were kept as attribution: piesos_LF.m ("based on
  V3_Liediff by CR"), PIESOS_stability_local.m ("based on LocalStability_v2 by CR"),
  @intop/private/enumerate_valid_orders.m ("moved from 'mtimes_tensopvar'").

## Unresolved

- Step comments (rule 5) not completed: 52 runs of 15 or more consecutive code lines without a comment
  remain, mostly display/validation code (@intop combine_terms 4, isvalid 3, display_term 3;
  @tensopmat get_display 3, initialize 2; @tensopvar get_display 2, display 2; @polyopvar
  get_display 2; mtimes_tensopvar 2; tensopmat2copvar 2; PIESOS_stability_local 2; one each in 25
  other files incl. @quadPoly diff, jacobian, hessian, subs, disp, unionBasis).
- Local functions with a descriptive header but no INPUT/OUTPUT block (17, all @quadPoly):
  combine_one_side, tensorExpAt, leftNumericTimes, rightNumericTimes, productBasis1D,
  pairToProdIndex, kronStrides, mergeSortedCellstr, mergeSortedUnique, tensorIndexMap,
  fullDegreeTable, randExpCell, tensorDecomp, subsToDegrees, rebuildBasisFromDegrees, degreesToSubs,
  subv2lin_rows.
- Nine SS/AT @quadPoly headers describe inputs and outputs under "Inputs:/Outputs:" or "Input/Output"
  (Sym2quadPoly, decomposeLR, horzcat, plus, polynomial2quadPoly, unionBasis, randquadPoly, trace,
  vertcat); left as written (other authors, no restyle).
- Error-message strings that cite other functions or "Def. 4" (code) were not changed.

## codeident counts

All 235 edited files: OK 235, DIFF 0, NEW 0 (hessian.m unchanged). Line endings: CRLF kept, no
final-newline changes. After the pass, all 236 files carry the license block, the contact line and an
`Initial coding` line.
