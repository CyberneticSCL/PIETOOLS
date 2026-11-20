# Ledger – E (PIESOS/), resume 11/20/2026

Lists: `E.tsv` (104 newfix, 2 new files), `E_runs.txt` (83 runs of 20+ uncommented code lines in 52
files), the 17 @quadPoly local functions without INPUT/OUTPUT. 141 files edited, all in PIESOS/.
Every file was codeident OK (or, for `internal/get_tensorvars.m`, `@quadPoly/times.m` and
`@quadPoly/hessian.m`, clean in git) right before its first edit.

## Edited files

Common to the 104 newfix files: the line `% MMP, 11/19/2026: Release comment pass (comment-only).`
removed (script, exact match); every other dated log entry (the `Initial coding` line kept), every
per-line marker (trailing blanks trimmed) and every `(was)` line removed. No Start/End markers
existed. Afterwards no newfix file has a dated comment other than its `Initial coding` line. Files
not listed below received only the release-line removal and were checked against rules 1–7
(headers, local headers, license, contact and Initial coding lines, pointers).

### Newfix files with more than the release-line removal (paths below PIESOS/polyopvar/)
- @intop/display.m: log entry removed; 'get_display' pointer removed from OUTPUTS; 3 step comments.
- @intop/isvalid.m: log entry removed; NOTES: the intervals in each intmat{i} must be ordered
  (restated from the log); step comments for 5 runs, incl. the cycle check in is_consistent.
- @intop/private/display_term.m: 3 log entries removed; pointers removed from the header ("for use
  in 'display'", "using get_int_order") and step comments (get_int_order, varname2display); step
  comments for 3 runs.
- @intop/times.m: inlining-pass log entry, a (was) line, "former local TIMES_PARAM" wording and 12
  markers removed; step comment on the coefficient layout (rows m*(nd+1), columns n*nZ).
- @intop/subsref.m: 4 step comments (index checks and selection).
- @polyopvar/diff.m, vertcat.m: log entries removed; NOTES: an adjoint input gives an error
  (restated); vertcat: constant2polyopvar pointer replaced, step comment.
- @polyopvar/display.m: 3 log entries removed (adjoint display, call_name already in the header);
  'get_display' pointer removed; step comments.
- @polyopvar/horzcat.m: log entry removed; constant2polyopvar pointer replaced.
- @polyopvar/initialize.m: log entry removed; 2 step comments.
- @polyopvar/mtimes.m: 3 log entries, 1 marker removed; "computed using times" -> "the elementwise
  product A.*B".
- @polyopvar/plus.m: 2 log entries, 1 marker removed; 'lift_output' pointer removed.
- @polyopvar/repmat.m: 3 step comments.
- @polyopvar/subs.m: 4 log entries removed; NOTES: an adjoint f gives an error (restated);
  "(clean_poly)" replaced by what the cleaning removes; subs_const header: apply_op, subs_factor,
  clean_poly pointers removed; 2 step comments.
- @polyopvar/subsref.m: 4 log entries removed; NOTES: '()' and '{}' on an adjoint give an error
  (restated); dangling " ," left by the 11/19 removal of "(see 'size')" repaired.
- @polyopvar/times.m: log entry removed; step comment.
- @quadPoly/combine.m (SS): 6 commented-out index formulas of a superseded ordering removed (they
  contradicted the code); combine_one_side INPUT/OUTPUT block; step comment.
- @quadPoly/diff.m, int.m (AT): commented-out code (2 lines each) and the TODO line "we need to call
  clean quadpoly after diff" removed (candidate 3); INPUT/OUTPUT completed ('double' input, unknown
  names, limits: one limit is the upper limit with lower limit 0); int.m: comment "Update
  left-monomials" over the right-monomial loop corrected; step comments (diff: run 71-138).
- @quadPoly/disp.m (SS): inlining-pass log entry, 2 (was) lines, "former local(s)" wording and 15
  markers removed (unicode string lines verified unchanged); first help line named FORMAT (a local
  function), now DISP; "(unchanged, now fed ... from tensorExpAt)" removed; tensorExpAt
  INPUT/OUTPUT; step comments.
- @quadPoly/eq.m (AT): stale "Date: 6/13/19" removed (contradicts Initial coding AT 01/21/2026);
  8 lines of commented-out ndopvar code removed; INPUT: scalar 0 allowed, tol optional (1e-14), the
  test stated as every coefficient of P1-P2 below tol (was "max(P1-P2)<tol"); step comments.
- @quadPoly/isempty.m (SS): first help line "almost trivial function. may have its uses." replaced
  by what the function does.
- @quadPoly/jacobian.m (AT): INPUT named P1 for the argument F – corrected; vars optional (default
  all variables), stacking of H stated; step comments.
- @quadPoly/mtimes.m (DJ, SS): INPUT/OUTPUT for leftNumericTimes, rightNumericTimes;
  productBasis1D "Inputs:/Outputs:" -> INPUT/OUTPUT; pairToProdIndex heading typo INPTUS; find_idcs:
  csubs was described as "k x m ... in ZL" – corrected to k x n in ZR.
- @quadPoly/repmat.m (AT): equalled HEAD after the release-line removal; its header was copied from
  blkdiag ("varargs as the blocks on the diagonal", INPUT varargin) – replaced by what the code does
  (tiles P1 N x M; N, M optional); step comments.
- @quadPoly/reshape.m (AT): 19 lines of commented-out test code removed; header no longer claims
  RESHAPE(A,[M,N]) (candidate 4); trailing comments "rows: n blocks of size ds", "cols: m blocks of
  size dt" corrected to M blocks of size n_Zs, N blocks of size n_Zt; step comments.
- @quadPoly/subs.m (AT, DJ): 2 log entries and 2 commented-out lines removed; "per quadPoly.m's ...
  convention" pointer removed; step comments (left/right substitution factors, hatC = L'*C*R).
- @quadPoly/var_swap.m (SS): 6 commented-out index formulas removed; comment citing the unused
  sPosMap/tPosMap corrected; INPUT/OUTPUT for tensorDecomp, subsToDegrees,
  rebuildBasisFromDegrees (posMap stated as returned empty; the old line claimed a map),
  degreesToSubs, subv2lin_rows.
- @tensopmat/display.m: log entry removed; 'get_display' pointer removed; step comments.
- @tensopmat/initialize.m: 6 step comments (output space per row; block conversion; lookup_dom).
- @tensopmat/isvalid.m: 1 step comment.
- @tensopmat/mtimes.m: log entry removed; "(per @tensopmat/initialize.m's row check)" removed; step
  comments.
- @tensopmat/plus.m: 2 step comments.
- @tensopmat/subsasgn.m: 2 log entries removed; pointers to 'subsref' and 'initialize' removed.
- @tensopmat/subsref.m: log entry removed; "(like subsasgn for C{I,J}=RHS)" removed.
- @tensopvar/display.m: log entry removed; 'get_display' pointer removed; step comments.
- @tensopvar/eq.m: 2 step comments.
- @tensopvar/int.m: 2 log entries and 2 markers (one on the function line) removed; "as named by
  `apply_op`" removed.
- @tensopvar/isvalid.m: truncated comment "Check that the different factors" completed; step
  comment.
- @tensopvar/mtimes.m: 2 log entries, 1 marker removed; NOTES: for different numbers of factors the
  factors acting on R are merged into their neighbors, else an error (restated from the log).
- @tensopvar/subsref.m: log entry removed.

### New files (full pass, other-author rules)
- internal/get_tensorvars.m (DJ): NOTES listing where the naming convention is used (intop,
  apply_op, piesos_quadopvar, piesos_quadvar, apply_poly, subs, sopvar) reduced to the convention;
  "as Kop.pvars for 'intop'" removed; contact line added; 1 step comment.
- @quadPoly/times.m (DJ): "the same restrictions as in 'mtimes'" replaced by the restriction;
  contact line added (Initial coding line kept as written, 01/21/2026; first git commit 11/19/2026);
  step comment.

### Local-function headers (17 @quadPoly functions; INPUT/OUTPUT added)
combine_one_side (combine.m), tensorExpAt (disp.m), leftNumericTimes, rightNumericTimes,
productBasis1D, pairToProdIndex (mtimes.m), kronStrides (polynomial2quadPoly.m, which also lost 6
lines of an alternative stride loop and 6 lines of commented-out missing-variable checks the code no
longer performs), mergeSortedCellstr, mergeSortedUnique, tensorIndexMap (private/unionBasis.m, which
also lost the question "maybe we do not need PA and PB, ever?"), fullDegreeTable
(quadPoly2polynomial.m, 2 commented-out lines removed), randExpCell (randquadPoly.m), tensorDecomp,
subsToDegrees, rebuildBasisFromDegrees, degreesToSubs, subv2lin_rows (var_swap.m).

### Runs (step comments only, unless noted)
custom/DMB.m (CR), custom/polyopvar_equal.m, demos_and_examples/examples_PDE_library_PIESOS.m,
executives/PIESOS_stability_local.m (DJ, CR), executives/PIESOS_sublevel.m,
internal/otimes_params.m, monopmat.m, piesos_ineq.m, piesos_intopvar.m (also 2 lines of trailing
blanks left by 11/19 marker removal trimmed), piesos_poslpivar.m, piesos_quadopvar.m,
piesos_quadvar.m, @intop/canonicalize_intmat.m, @intop/combine_terms.m, @intop/get_display.m,
@intop/mtimes_tensopvar.m, @polyopvar/combine_factors.m, @polyopvar/get_display.m,
@polyopvar/opower.m (also "as in 'otimes'", "followed by combine_terms" removed),
@polyopvar/repelem.m (also "rejected in repmat.m" replaced by the reason), @quadPoly/hessian.m (AT:
also INPUT P1 -> F, vars optional, H(i,j) = d^2F/(d vars{i} d vars{j})), @quadPoly/quadPoly.m (one
comment states the sorting branch never runs), @tensopmat/get_display.m, @tensopmat/opower.m and
@tensopvar/opower.m (also "as in 'otimes'" replaced by the meaning of type), @tensopmat/otimes.m,
@tensopmat/tensopmat.m, @tensopvar/get_display.m, @tensopvar/otimes.m,
polyopvar/tensopmat2copvar.m, polyopvar/tensopvar2polynomial.m.

## Skipped
- None.
- PIESOS/piesos_quadvar.m: HEAD moved to 41c9210c (upstream removed the local vec_dpvar); the
  coordinator merged these comment edits into the new version; re-checked: codeident OK, no run.

## Candidates for the maintainer (not fixed; "inferred" = read, not run)
1. @intop/isvalid.m, is_consistent: `safe_vars = [safe_vars,find(safe_idcs)]` appends positions in
   the reduced var_list, not the variable numbers rem_vars(safe_idcs), from the second pass on.
   Inferred: for chains of three or more dependent variables a wrong variable can be marked
   resolved; hand traces with cycles of length <= 3 still return false.
2. @quadPoly/subs.m: `if size(vars) ~= size(vals)` errors only if both dimensions differ (if on a
   logical vector), so 1x2 vars with 1x3 vals pass. Inferred.
3. @quadPoly/diff.m, int.m: the removed TODO asked for a cleaning step after differentiation;
   none is applied, and no quadPoly clean method exists.
4. @quadPoly/reshape.m: RESHAPE(A,[M,N]) is not supported (header corrected, code unchanged).
5. @quadPoly/var_swap.m: rebuildBasisFromDegrees returns posMap as empty cells; sPosMap, tPosMap
   are unused.
6. @quadPoly/quadPoly.m: `if ~presorted && 0` never runs; it would call clean(obj), which does not
   exist for quadPoly.
7. @intop/mtimes_tensopvar.m (local function building ZL, ZR): with output variables,
   `Cmat = sparse(ridcs(:),cidcs(:),1)` is 1 x nZ^2, while a quadPoly of dim [1,nZ] with a left
   basis of nZ monomials needs nZ x nZ coefficients. Inferred; not reached for functionals (no
   output variables).
8. @tensopvar/isvalid.m: the error text says factors must be 'sopvar', but 'tensopvar' and 'intop'
   factors are accepted (string unchanged).
9. The 11/19 candidates (E_piesos_ledger.md) stand; not re-examined.

## Unresolved
- The nine SS/AT "Inputs:/Outputs:" headers stay as written (decision); horzcat.m among them keeps
  "Uses unionBasis convention ... but we lift via maps" (pointer, "we").
- @quadPoly/numArgumentsFromSubscript.m (SS): rationale with "we" and the MATLAB overloads numel,
  subsref, subsasgn – left.
- Other authors' possessives left ("MATLAB's" in @polyopvar/subsref.m and @quadPoly/times.m,
  "C_new's" in @quadPoly/subs.m, "B's" in @tensopmat/mtimes.m).
- @tensopmat/initialize.m and @tensopvar/initialize.m (DJ) name 'opvar'/'opvar2d' in NOTES: they
  describe a conversion the code performs; kept.
- Not on my lists, seen with pointers, not edited: @polyopvar/clean_poly.m (header cites
  combine_factors.m, combine_terms.m, @tensopmat/clean_op.m).
- 88 runs of 15–19 lines remain in 46 of the 141 files (judgment; each lies under a block comment).

## Counts
- codeident.sh over the 141 edited files: OK 141, DIFF 0, NEW 0; over all 238 PIESOS .m files:
  OK 238 (HEAD 41c9210c).
- runs.awk (N = 20) over every file on my lists (141): 0 runs (start: 83 in 52 files); over all of
  PIESOS: 0.
- CRLF kept (238 w/crlf); the final-newline state of every PIESOS file equals HEAD.
- Only PIESOS/polyopvar/int_simple.m (new = 0) carries the release line; no per-line marker or
  (was) line remains in PIESOS.
