# Ledger – group A_sopvar (sopvar/, 178 files) – PAUSED 11/19/2026

Status at the pause: 111 files fully passed, 23 files touched only by the marker strip (step 1 below), 44
files not touched. codeident.sh over every modified file: **134 OK, 0 DIFF, 0 NEW**. No MATLAB run, no
git state changed. git status of sopvar/ was clean before the first edit; every modified file is from
this pass.

## Method
1. Script pass over the 94 new files (`strip_markers.pl`, scratchpad): 931 per-line markers cut from line
   ends, 145 commented-out `(was)` lines deleted (4 of them the first lines of a continued old
   statement). codeident: 94 OK.
2. Per file, by hand: logs (new files), BEGIN/END markers, references, history wording, missing support
   and Initial coding lines, INPUT/OUTPUT gaps, step comments; log line
   `% MMP, 11/19/2026: Release comment pass (comment-only).` on every edited old file. Marker-bearing
   comment lines of old files reworded with the marker kept at column 77 (`retext.pl`).
3. License end line copied as it stands in the source files ("authorship, and a brief description of
   modifications"); no license block was added yet (the two files without one are not yet passed).

## Edited files (fully passed)
- @cdopvar (30 files, all): Sec. 8/8.1/8.2/8.3.1 and Sec. 4 references, pointers to cat_copvar_grid,
  merge_dvar_lists/pair, merge_copvar_registry, dpvar_op_copvar, eq_copvar, verify_copvar_meta,
  op_slice, getfield/setfield/is_opvar*, plus_batch, @sopvar/@sdopvar methods out; dpvar and scalar
  operands stated in INPUTS (taken from the logs); cdopvar.m: logs and BEGIN/END block out, row/column
  consistency and dpvar dispatch stated, constructor step comments; merge_dvar_pair.m: PDF Sec. 8.1.4
  reference replaced by the Combine-Basis formula; plus_batch.m: heatNd bench reference replaced by
  "sums of 2N+1 terms in N variables"; DJ files (display, end, get_display, getfield/setfield,
  is_opvar1d/2d, op_slice): rule violations only, support lines added; log lines on old files.
- @copvar (33 files, all): same treatment as the @cdopvar twins; copvar.m: logs (kept MP initial coding),
  BEGIN/END, opvar2d sentence, file references out, INPUTS/OUTPUTS added (were missing), protected
  properties described as never filled; inv.m: @opvar/inv and inv_opvar comparisons, @sopvar/inv
  pointers (options and info fields now listed), opvar2d/GETOBSERVER_SOP note out; subsref.m: demos
  reference out; verify.m: logs and BEGIN/END out.
- @sdopvar (34 files, all): mtimes.m: sopvar.pdf sec. 6 and @sopvar/mtimes refs out, two step comments;
  ctranspose, eq, mtimes, plus: support line and "Initial coding DB - 08/10/2026 or 08/17/2026" (git log;
  headers had none); sdopvar.m: logs out, "Initial coding SS - 08/01/2026" (git), kernels property and
  dpvar dispatch stated, canonicalize/eq/lpi_eq pointers out; ChangeDecVar.m: logs out (LOC limitation,
  row sum of repeated names, Pi formula, memory stated), spec sec. 8.1.4 out, "Initial coding DB -
  08/10/2026" (git); sync_basis.m, setdvars.m, plus_batch.m, apply_basis_map.m, leftShiftMonomials_SD/SS,
  lrmultiply.m: logs and history ("replaces", "former local", "public entry point") out;
  numArgumentsFromSubscript.m, subsasgn.m: "as before", "pre-existing", "today" out; DJ files
  (display, display_term, decode_params, params2kernel, get_display, op_slice): violations only.
- @sopvar (13 of 40): blkdiag, cat, clean_op, ctranspose (sopvar document Sec. 4 out), degbalance
  (@opvar/@dopvar2d counterpart and degbalance_core pointers out, rule stated), diff, display, end, eq,
  get_display (usage pointers, see-refs out, log entry out, step comments), horzcat, inv (sopvar document
  Sec. 4, GK_GRID/CHEBFIT_MONO pointers out), vertcat.
- misc/degree_sizing/degbalance_core.m: logs out; @opvar/@dopvar2d/possopvar/eq_opts_sopvar and
  verification-script references out; rule restated (ceil of the per-subset maximal degree), opts and
  output fields described; a clause contradicting the code on the slot order removed (code: int <-
  target input degree, mult <- target output degree); "Deleted: degree_table" note out.

## Not yet passed (to resume)
Marker strip only (logs, references, headers still to do; code identical to HEAD), 23 files:
@sopvar/private/leftShiftMonomials_SS.m, @sopvar/repelem.m, @sopvar/sopvar.m, @sopvar/sopvar2sdopvar.m,
constructors/eye_copvar_sop.m, int2sopvar.m, mat2copvar_sop.m, zeros_copvar_sop.m,
misc/canonical_form/canonical_adjoint_map.m, canonical_var_order.m, canonicalize_multiplier.m,
is_canonical_multiplier.m, monomial_gather.m, misc/cat_copvar_grid.m, chebfit_mono.m,
conventions/parse_copvar_spaces.m, dpvar_op_copvar.m, eq_copvar.m, gk_grid.m, int_semisep.m,
mat2copvar_grid.m, merge_copvar_registry.m, spIkron.m.

Not touched, 44 files:
@sopvar/int.m, inv_1D.m, is_id_op.m, kron.m, lift_output.m, minus.m, mtimes.m,
numArgumentsFromSubscript.m, op_slice.m, plus.m, private/UnionBasisMonomials.m, private/decode_params.m,
private/display_term.m, private/params2kernel.m, repmat.m, size.m, subs.m, subsasgn.m, subsref.m,
transpose.m, uminus.m, unkron.m, uplus.m; constructors/int2copvar.m, pi4_copvar_sop.m;
misc/apply_sopvar.m, canonical_form/prune_zero_monomials.m, compose_sopvar.m, conventions/cell_of_gamma.m,
check_reserved_names.m, copvar_space_list.m, gamma_of_cell.m, kron_split.m, kron_strides.m,
monomial_position.m, monomial_position_map.m, multiindex_grid.m, derive_copvar_meta.m,
on_registry_sop.m (no license block), parse_diff_inputs.m, parse_int_inputs.m, pvar2varname.m,
split_copvar.m, verify_copvar_meta.m.
Also pending: spIkron.m has no license block. @sopvar/blkdiag, horzcat, vertcat carry their own initial
coding line (AT 01/21/2026), so the 12/31/2021 import commit was not used.

## Skipped
None (no file in the list had uncommitted edits of another session).

## Candidates for the maintainer (not fixed)
1. @copvar/copvar.m: the protected properties leftCommonBasis, leftCommonC, rightCommonBasis, rightCommonC
   are written by no method on the path (only orphaned_file/private/@copvar/canonicalize.m): dead.
2. @sdopvar/private/lrmultiply.m: `Bout = sparse(...)` before the product is dead, overwritten by
   `Bout = B*LRkron` (noted on 11/06, still present).
3. @sopvar/ctranspose.m: the assignments to At before the constructor call and S3Var, S3Vardummy are
   never read (noted on 11/06, still present).
4. @sdopvar/sdopvar.m help (SS): "Zd is a vector of decision variables, sorted" – merged lists are stable
   unions (unique 'stable'), not sorted.
5. misc/degree_sizing/degbalance_core.m: the removed comment said the possopvar documentation of the
   'int'/'mult' slot correspondence is wrong; possopvar is outside group A, not checked.

## Unresolved
None among the passed files. Citations kept (literature attribution): Shivakumar, Das, Peet, arXiv
2208.13104, Lemma 18 (@copvar/inv.m) and Sec. VII, Lemma 16, Cor. 17 (@sopvar/inv.m).

## codeident (all modified files, at the pause)
OK 134, DIFF 0, NEW 0.
