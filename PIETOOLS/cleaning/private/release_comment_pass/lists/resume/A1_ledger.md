# Ledger – A1 (sopvar, the 67 files unfinished on 11/19) – 11/20/2026

Status: all 67 files of `A1.tsv` passed (23 "partial", 44 "untouched"; all new = 1). **codeident.sh over
all 67: OK 67, DIFF 0, NEW 0.** **runs.awk (N = 20) over all 67: 0 runs** (16 runs of 15–19 lines remain,
left by judgment). No MATLAB run, no git state changed. Every untouched file was clean in git before its
first edit; every partial file was codeident OK before its first edit (re-checked after the 04:00 resume:
all 67 OK, the last file edited, gk_grid.m, intact). Line endings and final-newline state of every file
unchanged (original state read off the index-recorded file size: 40 CRLF, 27 LF). Temp files in
`scratchpad/A1/`.

## Method
1. `strip_markers.pl` (the 11/19 script) on the untouched files that carried per-line markers: mtimes.m
   (151 markers, 12 `(was)` lines), subsref.m (23, 1), repmat.m (1); plus.m by hand (13, 1).
2. By hand, every file: change logs out (every dated entry; `Initial coding` lines kept, history text in
   them cut to author and date), BEGIN/END and Start/End markers out, references out (sopvar spec sections,
   sopvar.pdf, sopvar_implementation_notes.pdf, other files/functions, usage), history wording and
   possessives out of MMP-session text; support line added where missing (35 files); license block added
   where missing (2); INPUT/OUTPUT for main and local functions; step comments.
   Facts a user needs that only a removed log entry stated were restated as description (listed below).

## Edited files
@sopvar
- sopvar.m: logs out; dims doc corrected (see Factual); multi-index labels gamma 1/2/3 = alpha 0/1/-1
  stated; `kernels` property documented; usage sentence and pointers (canonicalize_multiplier, eq,
  lpi_eq_sdopvar, display, get.kernels) out; support line; get.kernels INPUT/OUTPUT; step comment.
- int.m, is_id_op.m, kron.m, lift_output.m, unkron.m, subs.m (DJ): pointers (int2sopvar/mtimes, "see sopvar",
  "see also kron", "as in 'sopvar'") out; support; step comments (subs: 1 run).
- inv_1D.m: "(see @sopvar/inv)", decomposeLR pointer, 11/06 log (with caller list), "lambda removed",
  "replaces build_Phi" out; support; step comments (2 runs). Gohberg-Krein named as the method.
- minus.m, uminus.m, uplus.m (11/06 log out), transpose.m ("by calling 'ctranspose'" out): support, step.
- mtimes.m: markers/(was) stripped; log out; Sec. 6 / Sec. 5.4 refs, dpvar_op_copvar pointer, "legacy",
  3 "former local X", "important correction", int_semisep_AT pointer, possessives out; dpvar-factor
  semantics stated in INPUT (from the 09/29 log); step comments (5 runs).
- numArgumentsFromSubscript.m, subsasgn.m: history ("as before", "pre-existing", "unchanged", "still"),
  subsref/@sdopvar pointers, mixin.Scalar rationale, possessives out; Initial coding cut; step comments.
- op_slice.m (DJ): "moved from @sopvar/subsref" cut from its Initial coding line; support; step comment.
- plus.m: markers/(was) out; log out (AT 11/17 scalar summand stated in INPUT; zero-block fact kept as
  step comment); dpvar_op_copvar/"legacy"/sopvar.m pointers, possessives out.
- repelem.m (factual fix), repmat.m: support.
- size.m (AT): 11/06 log out; step comment. sopvar2sdopvar.m: Sec. 8 and sync_basis pointer out (rule
  stated); logs out; step comment.
- subsref.m (AT): markers/(was) out; log out; op_slice pointer, history out; support; step comments.
- private/decode_params.m, display_term.m, params2kernel.m (DJ): usage/pointers (params2kernel,
  get_display, varname2display, decode_params, gamma_of_cell, display) out (the integrals per multi-index
  entry stated instead); display_term 11/06 log out; support; step comments (display_term: 2 runs).
- private/leftShiftMonomials_SS.m (SS): logs, Sec. 5.1, @sdopvar-copy pointer, "zeroed" history out;
  classes in INPUT/OUTPUT; support and `Initial coding SS - 04/15/2026` (git; header had none); step
  comments (1 run).
- private/UnionBasisMonomials.m: support; step comments.

constructors
- eye_copvar_sop.m, mat2copvar_sop.m, zeros_copvar_sop.m: mat2opvar/opvar comparisons, POSCOPVAR /
  lpivar_cdopvar / H-infinity executives usage, Sec. 4, "as for 'mat2copvar_sop'" (inputs written out),
  "passing 'verify'" out; MAT2OPVAR dropped from the See also line (no such file outside cleaning/); logs
  out; step comments.
- int2copvar.m (DJ): "as for 'int2sopvar'" replaced by the four inputs; "wrapper" sentence out; support.
- int2sopvar.m (DJ): "as for 'copvar'" (4), "sopvar document, Sec. 4 (see 'gamma_of_cell')", "(sopvar.m)"
  out; logs and Start/End of changes markers out (s_dum rule already in INPUT); support.
- pi4_copvar_sop.m: opvar2copvar NOTES and Initial-coding history out; local header I/O; step comments.

misc
- conventions/ (10 files): spec-PDF references out (conventions stated), Initial-coding history cut,
  cross-function pointers/usage lists, possessives out; parse_copvar_spaces: "Differences from the parsers
  this replaces" block out, DJ 11/02 empty-dims fact and struct-array refusal stated in INPUT; step
  comments (parse_copvar_spaces: 3 runs).
- canonical_form/ (6 files): logs, BEGIN/END (canonicalize_multiplier; measured facts kept), trailing
  "Deleted: ..." history blocks (3 files), "former local", conventions-folder and spec pointers,
  UnionBasisMonomials/eq/lpi_eq_sdopvar/plus pointers, sossolve/SeDuMi wording out; local headers with
  INPUT/OUTPUT (9 local functions); support; step comments (8 runs).
- apply_sopvar.m (factual fix; commented-out `%x = subs(...)` line removed), compose_sopvar.m (mtimes
  rationale out; 3 runs), split_copvar.m (1 run), parse_diff_inputs.m, parse_int_inputs.m (ordering fact
  from the Initial-coding text stated in OUTPUT; 1 run), pvar2varname.m (factual fix) – DJ: usage out,
  support, step comments.
- cat_copvar_grid.m, derive_copvar_meta.m, eq_copvar.m (factual fix), merge_copvar_registry.m,
  verify_copvar_meta.m, mat2copvar_grid.m, dpvar_op_copvar.m: usage, legacy-opvar comparisons, Sec. refs,
  pointers (merge_dvar_lists, sync_basis, subsref, poscopvar, mat2opvar, getsol), BEGIN/END markers,
  possessives, logs out; MAT2OPVAR dropped from mat2copvar_grid See also; local headers I/O; step comments.
- chebfit_mono.m, gk_grid.m: Initial-coding history, "former local/nested", usage (@sopvar/inv,
  getController_direct_sop, gain construction), "map Sec. 8.2" out; local headers I/O; step comments.
- int_semisep.m (AT, MMP): Sec. 6.1/6.2 and notes-PDF refs, copquadvar/unpack_sheets/possopvar/heatNd
  usage, "four near-identical copies" history, BEGIN/END markers (restated as step comments), 13 log
  entries out; Initial coding cut to `AT, 2026: Initial coding`; partial index lists stated in INPUT;
  local headers in call form; support; step comments (2 runs).
- on_registry_sop.m, spIkron.m: **license block added** (copied exactly from PIETOOLS_stability_sop.m;
  verified by diff); support line; logs and Initial-coding history out; step comments.

## Factual errors corrected (comments only, verified from the code)
1. @sopvar/sopvar.m: P.dims documented as "2x1 ... [p,q]"; the code uses 1 x 2 [q,p] (size.m, the
   constructor INPUT, coefficient size dims(1)*nL x dims(2)*nR).
2. @sopvar/repelem.m: step comment said "repmat of the parameters"; the code repeats each operator
   row/column index (repelem).
3. misc/apply_sopvar.m: INPUT said x in L2^{m}; x is n x q, so L2^{n}.
4. misc/eq_copvar.m (local `cut`): said slicing errors on a [] / scalar-0 parameter cell; @sopvar/op_slice
   and @sdopvar/op_slice both keep such cells as zero cells.
5. misc/pvar2varname.m: license title line read "PIETOOLS - nopvar"; now "PIETOOLS - pvar2varname".

## Skipped
None.

## Candidates for the maintainer (not fixed)
1. misc/dpvar_op_copvar.m: two error messages (code strings) say "(legacy: nnz(dim(:,1)) > 1)" /
   "(legacy: nnz(dim(:,2)) > 1)" – history wording in user-visible text.
2. misc/mat2copvar_grid.m: error 'mat2copvar_grid:integral' says "refused here as mat2opvar refuses Q1";
   mat2opvar.m no longer exists outside cleaning/.
3. misc/cat_copvar_grid.m: the BadOperand messages advise "convert others first, e.g. with opvar2copvar"
   (the opvar route; converters/operator_converters/opvar2copvar.m still exists).

## Unresolved
1. misc/int_semisep.m: `AT, 2026: Initial coding` has no month/day; git first commit of this name is
   Talitsky 08/17/2026, after a rename from int_semisep_AT, so the original date is not established.

## Other notes for the maintainer
- Seven files show a first git commit of 12/31/2021 (rename-follow): minus, size, uminus, uplus,
  int2copvar, check_reserved_names, kron_split. Each header carries its own 2026 Initial coding line,
  kept; no 12/31/2021 date written.
- MAT2OPVAR removed from two `See also` lines (retired function, no file) – a deviation from "See also
  lines kept"; reversible.
- The Edit tool dropped one trailing blank on a comment line of leftShiftMonomials_SS.m (comment only).
- Other authors (DJ, AT, SS, MP): rule violations only; their wording ("we", one "constructor's" in DJ's
  op_slice.m) kept.
- Citations kept: gk_grid.m (Shivakumar, Das, Peet, arXiv 2208.13104, Lemma 16 and Cor. 17); inv_1D.m
  (Gohberg-Krein method).

## Final counts (all 67 files of A1.tsv)
codeident.sh: OK 67, DIFF 0, NEW 0. runs.awk N = 20: 0 (N = 15: 16).
