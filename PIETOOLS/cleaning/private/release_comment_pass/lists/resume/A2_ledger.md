# Ledger – package A2 (sopvar newfix + step comments) – 11/20/2026

82 files: 68 in `A2_newfix.txt`, 22 in `A2_runs.txt` (8 in both). Before the first edit every file was
modified by the 11/19 pass and codeident printed OK for all 82, so none was skipped. No MATLAB run, no git
state changed, no file outside the two lists edited. Line endings kept (68 CRLF, 14 LF; the four files
without a final newline still have none). Work files: scratchpad `A2/` (backup of the 82 files as they
stood before this package: `A2/backup/before_strip.tar`).

## Method
1. Script (`strip_markers.pl` of 11/19, unchanged) over the 68 newfix files: 487 per-line markers cut,
   82 commented-out `(was)` lines deleted (5 of them first lines of a continued old statement). Then the
   line `% MMP, 11/19/2026: Release comment pass (comment-only).` deleted from all 68 (perl, exact match).
2. By hand: change logs (every dated entry; `Initial coding` lines kept, history appended to them cut),
   BEGIN/END markers, `Deleted MMP ...` notes, one leftover continuation line of a deleted old statement
   (@sdopvar/mtimes.m), "former local" and "as before" wording. User facts from removed entries restated in
   the present tense (listed per file). Then rules 1–7, local-function headers included.
3. Step comments wherever `runs.awk -v N=20` reported a run, including one run the marker strip created
   (@copvar/inv.m, kmat).

## Edited files – newfix (68)
All 68: markers, `(was)` lines and the 11/19 log line out (step 1). Further changes:
- @cdopvar/blkdiag: prose dash -> parentheses.
- @cdopvar/ctranspose: log (rename, size) and "Split out of @copvar" out; prose dash.
- @cdopvar/display (DJ): log entry out; 2 step comments (run 124-144).
- @cdopvar/end (DJ), op_slice (DJ), uminus: step 1 only.
- @cdopvar/eq: log out (the 4 (was) lines of the old OUTPUTS entry went in step 1).
- @cdopvar/horzcat, vertcat, size: logs out (dpvar entry already stated in INPUTS).
- @cdopvar/minus: INPUTS now state the numeric/'dpvar' scalar and 'dpvar' matrix operands that plus accepts.
- @cdopvar/mtimes: log out; "Deleted MMP" note replaced by the fact it held (promoted fixed factor has an
  empty list, the decision factor keeps its list, O(blocks), no sort); prose dashes; rewrap.
- @cdopvar/plus: log out; "Deleted MMP" note out; prose dash; rewrap.
- @cdopvar/numArgumentsFromSubscript: log and "with subsref/subsasgn" out; NOTES in present tense (was
  "silently kept"); possessive; classes quoted.
- @cdopvar/subsasgn: log/history out; first line names the 'dopvar'/'dopvar2d' field exception;
  possessives; "(measured)" out; rewrap; 2 step comments (run 129-148).
- @cdopvar/subsref: log (incl. matlab.mixin.Scalar rationale) out; OUTPUT covers slice, field, builtin;
  possessives; rewrap.
- @cdopvar/transpose, @copvar/transpose: past-tense "Why" note restated in the present tense.
- @copvar/blkdiag: possessives. @copvar/ctranspose: log out; dash; rewrap.
- @copvar/diff, end, int, op_slice, subs (DJ), uminus: step 1 only.
- @copvar/display (DJ): log entry out; 2 step comments (run 123-144).
- @copvar/eq: log out; caller variable in "'~(Twop==0)'" replaced by ~(P==0).
- @copvar/horzcat, vertcat: logs out; possessive. @copvar/size: log out.
- @copvar/inv: log, BEGIN/END markers and three "the former local ..." notes out; See also COPVAR2OPVAR
  dropped (only in cleaning/private/dead_code); inv_sep INPUT/OUTPUT written out instead of "(see INV)";
  INPUT/OUTPUT added to 10 local and nested functions; possessive; step comments (runs 143-165, 422-450,
  kmat run from the strip, and three 16-19 line runs).
- @copvar/minus: INPUTS gain the scalar/'dpvar' operands. @copvar/plus: log out; dashes; rewrap.
- @copvar/mtimes: log out; restated "a numeric matrix M is an error; dpvar(M) multiplies by it"; dash; rewrap.
- @copvar/numArgumentsFromSubscript, subsasgn, subsref: as the @cdopvar twins (step comments, run 127-146).
- @sdopvar/blkdiag: log out; classes quoted; possessive.
- @sdopvar/ctranspose (DB): log, "(see the help text)", "the former local SIZE_DECVARS" out.
- @sdopvar/degbalance, @sopvar/degbalance: "thin wrapper ... matching '@opvar' ..." out; See also POSSOPVAR
  dropped (orphaned_file/private, off the path).
- @sdopvar/display (DJ): log out; step comment (run 107-127).
- @sdopvar/end: OUTPUT states the linear-index case (n = 1); pointer to 'size' out; classes quoted.
- @sdopvar/eq (DB): log out; "see the note in eq" replaced by the fact.
- @sdopvar/horzcat, vertcat: logs out; classes quoted; vertcat "Cross-check" restated as the identity
  [A;B] = ([A',B'])'.
- @sdopvar/minus: INPUTS gain scalar/'dpvar' operands; rewrap.
- @sdopvar/mtimes (DB): log out; facts restated (entries that cancel to rounding level are stored, not
  pruned; canonical order vars.out = [S2,S3], vars.in = [S3,S1], checked against canonical_var_order);
  leftover continuation line of a deleted old statement removed; an MMP comment saying kron(MR,I) "is built"
  contradicted DJ's code (kron not formed) and was removed, DJ's comment kept; two "former local" notes
  out; possessives (B's, A's, adjoint's, int_2b's, MATLAB's); step comments for runs 388-423, 461-507,
  631-660, 809-847, 1261-1283, 1471-1497 and four 16-19 line runs.
- @sdopvar/numArgumentsFromSubscript: history appended to Initial coding out; possessive; class quoted.
- @sdopvar/op_slice (DJ): "(slicing moved from @sdopvar/subsref)" out of the Initial coding line.
- @sdopvar/plus (DB): log and BEGIN/END out; the BEGIN note restated as a present-tense step comment.
- @sdopvar/private/display_term (DJ): log entry out.
- @sdopvar/repmat: classes quoted. @sdopvar/size: pointer to 'end' out; classes quoted; step comment.
- @sdopvar/subsasgn: history/rationale out; possessives; "(measured)" out; rewrap.
- @sdopvar/subsref: log and "as before" out; classes quoted; rewrap; For support line added (was missing).
- @sdopvar/transpose, uminus: step 1 only.
- @sopvar/blkdiag, horzcat, vertcat (AT): logs out ("Initial AT 01/21/2026" kept); past-tense guard comment
  restated; possessive; horzcat/vertcat OUTPUT: 'cdopvar' with a decision 'dpvar' entry, else 'copvar'.
- @sopvar/ctranspose (MMP, SS): log out; "(see the help text)" out.
- @sopvar/diff (DJ), eq (AT): logs out. @sopvar/display (DJ): log out; step comment (run 103-123).
- @sopvar/end (AT): OUTPUT states the linear-index case.
- @sopvar/inv: log ("Container counterpart of @opvar/inv", file moves, test) out; See also
  GETCONTROLLER_DIRECT_SOP (cleaning/private/dead_code), INV_OPVAR_2, INV_OPVAR (orphaned_file/private)
  dropped; "see header" replaced by the reason.

## Edited files – step comments only (14)
@cdopvar/get_display (6 comments; runs 126-145, 218-250, 346-365), @cdopvar/getfield_opvar1d (2;
118-140), @cdopvar/setfield_opvar1d (3; 82-106), @cdopvar/setfield_opvar2d (7; 81-101, 213-239, 377-404),
@copvar/get_display (6; 124-143, 216-248, 344-363), @copvar/getfield_opvar1d (4; 112-134, 161-184),
@copvar/getfield_opvar2d (2; 179-202), @copvar/setfield_opvar1d (3; 82-106), @copvar/setfield_opvar2d (3;
81-101, 209-238), @sdopvar/get_display (2; 167-192), @sdopvar/private/leftShiftMonomials_SD (2; 65-90),
@sdopvar/private/leftShiftMonomials_SS (5; 99-122, 127-159, 203-230), @sdopvar/private/sync_basis (2;
89-108; INPUT/OUTPUT added to the local basis_positions), @sopvar/get_display (2; 163-188).

## Skipped
None.

## Candidates for the maintainer (not fixed)
1. The opvar-field interface of 'copvar'/'cdopvar' (subsref, subsasgn, numArgumentsFromSubscript,
   getfield_opvar1d/2d, setfield_opvar1d/2d, is_opvar1d/2d) names the retired 'opvar', 'opvar2d', 'dopvar',
   'dopvar2d' in help text and in error/warning strings (code). Comments left consistent with the strings.
   Decide whether these names stay in the release.
2. @sopvar/ctranspose.m: the assignments to At before the constructor and S3Var, S3Vardummy are still dead
   (A_sopvar ledger, candidate 3).
3. Outside A2: sopvar/misc/gk_grid.m and chebfit_mono.m (A1) name getController_direct_sop, which is in
   cleaning/private/dead_code; converters/operator_converters/opvar2copvar.m and pie2copvar.m name
   copvar2opvar, also dead code.
4. @sopvar/end.m, @sdopvar/end.m: n > 2 leaves out unassigned (MATLAB "output not assigned" error).
5. @sdopvar/mtimes.m (DB) uses "we" in the int_2b derivation comments; @sdopvar/op_slice.m (DJ) has
   "constructor's". Other authors' wording, left.

## Unresolved
None. Literature citations kept: Shivakumar, Das, Peet, arXiv 2208.13104, Lemma 18 (@copvar/inv.m);
Sec. VII, Lemma 16, Cor. 17, Gohberg-Krein (@sopvar/inv.m).

## Decisions taken (for review)
- See also lines kept, but targets not on the release path dropped: COPVAR2OPVAR, POSSOPVAR (2 files),
  GETCONTROLLER_DIRECT_SOP, INV_OPVAR_2, INV_OPVAR.
- Prose dashes replaced by colons, parentheses or commas rather than en dashes: every sopvar .m file is
  ASCII (display_term.m declares its characters by code "to keep the file ASCII").
- DJ-format initial lines (`% DJ, 11/06/2026: Initial coding;`) and `% Initial AT 01/21/2026` kept as written.

## Final counts (all 82 files)
codeident: OK 82, DIFF 0, NEW 0. runs.awk N = 20: 0 runs (38 at the start); N = 15: 38 runs of 15-19
lines left by judgment (about half in DJ files, 7 in @sdopvar/mtimes.m).
