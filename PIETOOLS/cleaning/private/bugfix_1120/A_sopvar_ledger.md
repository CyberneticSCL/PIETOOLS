# Bug-fix pass 11/20/2026 – agent A_sopvar – ledger

Items: C1, C3, C4, C5, C6, C8, C9, C10, C11, C101 (copvar/cdopvar subsref, subsasgn only), C102, C103,
C104, C105, C107, C108, C109. HEAD 17d68747. Scratch folder (scripts, outputs):
`C:\Users\mpeet\AppData\Local\Temp\claude\363c499c6332a503\874c85cb-c18d-4876-bea3-39cdd053e21f\scratchpad\A_sopvar\`
(below: `<S>`).

Annotation regime: every file edited in this pass is listed in column 4 of
`cleaning/private/release_new_files_1120/new_files_classified.tsv` (classes New or Replacement), so the
new-file regime applies – no header log line, no dated markers, no (was) lines. `git status --porcelain`
was clean for every file before the first edit; line-ending states before editing are in
`<S>/eol_before.txt`, original copies in `<S>/orig/`.

Test scripts (all in `<S>`): `t1_class.m` (messages, `end`, Zd order, fixtures of ctranspose and
products), `t2_tests.m` (Testfolder tests), `t3_exec.m` (executives, C6 dispatch, C11), `t4_est2d.m`
(C8, 2-D estimator); outputs `t1_class_<tag>.txt`, `t2_<tag>_<k>.txt`, `t3_exec_<tag>.txt`,
`t4_est2d_<tag>.txt`, tag `before` / `after`. Every function under test resolved inside the working
tree (`which` lines at the top of each output).

## C1 – FIXED
- Files: `sopvar/@sopvar/end.m` (AT; lines 43-48 by `git blame`: AT), `sopvar/@sdopvar/end.m` (MMP).
  Each: an `else` branch with an error after `elseif n==2`, as in `sopvar/@copvar/end.m`:
  `error('sopvar:endIndex','A sopvar is indexed as P(I,J), or P(I) for one row or column.')`,
  `error('sdopvar:endIndex','An sdopvar is indexed as P(I,J), or P(I) for one row or column.')`.
  +2 lines each.
- Reproduction (measured, `t1_class_before.txt`): `Ps(1,1,end)` on a 2 x 2 sopvar and `Pd(1,1,end)` on
  a 2 x 2 sdopvar both raise `MATLAB:unassignedOutputs` "Output argument "out" (and possibly others)
  not assigned a value in the execution with "sopvar/end" function." (sdopvar/end likewise).
- Verification: see "after" below.

## C3 – FIXED
- File: `sopvar/@sdopvar/mtimes.m`, local `lrmultiply`: deleted `Bout = sparse(size(B,1),size(L,1)*size(R,2));`
  and the blank line after it (HEAD lines 1517-1518; author SS by `git blame -C`, line 1522 AT).
- Semantics (no run needed for the defect): `Bout` is assigned again by `Bout = B*LRkron;` with no read in
  between, so the first value is never used. Cost of the deleted line: an empty sparse with q rows
  (column pointers only), so no measurable change in time or memory is expected.
- Verification: see "after" below (C3/C4 share the chain runs).

## C4 – FIXED
- File: `sopvar/@sopvar/ctranspose.m`: deleted the comment and the assignments `At = A;
  At.vars_S2 = ...; At.vars_S1 = ...; At.dom_2 = ...; At.dom_1 = ...; At.dims = A.dims';` (HEAD lines
  49-58) and `S3Var = A.vars_S3; S3Vardummy = strrep(S3Var,'s','t');` (lines 61-62); authors MMP, SS.
  `nvars3 = numel(A.vars_S3);` kept (read below). -12 lines.
- Semantics: `At` is reassigned whole by `At = sopvar(Ct,At_vars,leftZ,rightZ,At_dom,At_dim);` at the
  end, with no read of `At` in between; inside a method of the class the dot assignments use the builtin
  subsasgn, and the properties have no set methods or validators (`sopvar/@sopvar/sopvar.m:116-134`),
  so the deleted lines have no side effect. `S3Var`, `S3Vardummy` are never read.
- Verification: see "after" below.

## C5 – FIXED (text)
- File: `sopvar/@sdopvar/sdopvar.m:64` (SS), class help: "Zd is a vector of decision variables, sorted
  (all monomial degrees 1)" -> "Zd is a vector of decision variables, in the order given, not sorted /
  (all monomial degrees 1)" (1 line -> 2).
- Measured (`t1_class_before.txt`): `sdopvar(...,{'d2';'d1'},...)` stores Zd = d2,d1; the sum of that
  operator and one on {'d3';'d1'} has Zd = d2,d1,d3 (issorted 0), and so does `[Q1,Q2]` – the stable
  union of `sopvar/@sdopvar/private/sync_basis.m:108`. No code reads Zd as sorted (grep for issorted /
  sort on Zd in sopvar, lpi_programming_sopvar, executives_sopvar: none).
- Verification: text only; the measured behaviour above is the reference.

## C6 – FIXED
- File: `executives_sopvar/PIETOOLS_auto_execute_sop.m:191, 196` (MMP by `git blame -C`): the two tests
  accept the prompt answers too: `if strcmpi(exec{j},'PIE2PDEstability') || strcmpi(exec{j},'stability')`,
  `elseif strcmpi(exec{j},'PIE2PDEstability_dual') || strcmpi(exec{j},'stability_dual')`. 2 lines changed.
- Reproduction (measured, `t3_exec_before.txt`; no prompt run): the dispatch part of the file (from
  "% The form of the stability LPI" to the end, read from the file under test and written to
  `<S>/dispatch_<tag>.m`) run in the base workspace with `exec` set to the cell that
  `split(input(...))` returns, PIE = 1-D heat equation, settings = 'light':
  `{'stability'}` -> "Unknown executive 'stability'."; `{'stability_dual'}` -> "Unknown executive
  'stability_dual'."; the flag names `{'PIE2PDEstability'}`, `{'PIE2PDEstability_dual'}` -> run, status
  optimal, forms Q and Qd.
- Verification: see "after" below.

## C8 – FIXED
- File: `executives_sopvar/PIETOOLS_Hinf_estimator_sop.m:159-165` (MMP): after an error of the default
  method, the 2-D branch calls `getObserver_sop(Psol,Zsol,struct('method','inv'))` instead of repeating
  `'qp'`; the verbose message says "inverse (inv)" instead of "least squares (qp)"; the comment says
  so; the help (line 18) gains one line: "(in two variables, should that computation fail, L = P^{-1}*Z
  for a separable P)". 3 lines changed, 1 added.
- Intent from history (`git log -L`): the branch was written on 11/04/2026 (e3d9c332, MMP) as "the
  container inverse of a separable P (COPVAR/INV); for a non-separable P the least-squares gain of the
  method 'qp'", when the default of getObserver_sop was 'inv'; on 11/17/2026 (33d30877, AT, "Change
  default inverse option for getController, getObserver") the default became 'qp', so the first call
  became 'qp' and the retry repeated it; the 11/20 comment pass described that state. The method list
  of getObserver_sop is {'inv','qp'} (line 94), so the alternative to the default is 'inv'. The other
  three synthesis executives call the default method once (no retry).
- Reproduction (measured): (1) `t3_exec_before.txt`: on the solved 1-D estimator,
  `getObserver_sop(P,Z)` and `getObserver_sop(P,Z,struct('method','qp'))` return equal gains
  (isequaln 1) – one computation. (2) `t4_est2d_before.txt`: 2-D heat equation with w, z, y,
  'stripped', qp_copvar replaced by a stub that errors (`<S>/stub_qp/`, first on the path): the
  executive prints "gain extraction: qp_copvar disabled by the test stub; least squares (qp)." and then,
  from the outer catch, the same error again; L = [] and info.inv = []. Without the stub: L copvar,
  info.inv {qp}; `getObserver_sop(P,Z,struct('method','inv'))` on that solved 2-D P, Z returns a copvar
  gain (info fields method, d, relrms, condIKH, Kdiff, nF, note).
- Alternative for the maintainer: delete the retry instead (the outer catch already reports the error),
  as in the other synthesis executives.
- Verification: see "after" below.

## C9 – FIXED
- File: `lpi_programming_sopvar/misc/lpi_classify_sop.m:84-86` (MMP): the if/else whose two branches are
  both `cl.obj = double(lpigetsol_sop(prog,obj));` replaced by that one statement, with the comment
  "a 'dpvar', 'polynomial' or name". -2 lines.
- Semantics: `isa` has no side effect and both branches evaluate the same expression; lpigetsol_sop
  takes 'dpvar', 'polynomial' and char alike (`lpigetsol_sop.m:82-84`). Duplicated since the file was
  created (25d6d82d, 10/09/2026).
- Verification: exercised by every executive run of `t3_exec_after` (lpi_solve_loop_sop classifies each
  solve); see "after".

## C10 – FIXED
- File: `executives_sopvar/utility_functions/pie_disc_sop.m:231` (MMP), local `entry`: deleted
  `if isa(X,'double'),     e = X(i,k);     return,     end`; the next line `e = X(i,k);` returns the same.
- Semantics: for a double X both lines evaluate `X(i,k)`; duplicated since creation (25d6d82d).
- Verification: pie_disc_sop runs in the 1-D stability witness (`settings.sop.witness`, maxre); see
  "after".

## C11 – FIXED (help text; field name kept)
- Files (MMP): `executives_sopvar/PIETOOLS_H2_control_sop.m:54-55`, `PIETOOLS_H2_estimator_sop.m:61-62`,
  `PIETOOLS_Hinf_control_sop.m:58-59`, `PIETOOLS_Hinf_estimator_sop.m:68-70` (HEAD lines). The description
  of info.inv now says what it holds: "(output INFO of GETCONTROLLER_SOP [GETOBSERVER_SOP], method 'qp':
  the field qp; no inverse of P is formed)"; in the H-infinity estimator also "for the 2-D fallback
  'inv' the data of the inverse" (C8). +2 lines in each of three files, +2 in the fourth. info.inv not
  renamed.
- Measured (`t3_exec_before.txt`): 1-D H-infinity estimator (the PDE of DEMO5), 'light': info.inv is a
  'struct' with the single field qp; `t4_est2d_before.txt`: 2-D, no stub, likewise {qp}. With the method
  'inv' (1-D) the third output has fields N, h, r1, r2, n, rcondU22, condR0, d, relrms, note.
- Verification: text only; the measured fields are the reference.

## C101 – HELD (lines of Declan Jagt)
- Files: `sopvar/@copvar/subsref.m:104-118`, `sopvar/@copvar/subsasgn.m:113-120, 151-152`,
  `sopvar/@cdopvar/subsref.m:104-118`, `sopvar/@cdopvar/subsasgn.m:115-122, 153-154`. Not edited.
- Reason: `git blame -w -C` attributes every message line of the field interface in these four files to
  Declan (3b92e68b, 11/15/2026, "Support old opvar fields in copvar subsref"; 2f2af6a1, 11/16/2026,
  "Update copvar/cdopvar subsasgn"); the brief holds any line whose author is DJ by `git blame`. The
  file headers name MMP (Initial coding 09/29/2026), hence the MMP credit in CANDIDATES.md.
- Measured (`t1_class_before.txt`, unchanged in `_after`): `T.P` on a 1-D copvar warns
  `copvar:opvarField` "A copvar has no field 'P'; it is computed from the blocks as for an 'opvar' or
  'opvar2d'. ..."; `P.P` on a cdopvar warns `cdopvar:opvarField` (" ... as for a 'dopvar' or
  'dopvar2d'."); `T.Q1` on a 2-D copvar raises `copvar:notOpvar1d` "'Q1' is a field of an 'opvar', but
  this copvar is not opvar-like: ...".
- Proposed wording for Declan (identifiers and behaviour unchanged): notOpvar1d "'%s' is a field of the
  1-D form R^m x L2^n[a,b] (fields P, Q1, Q2, R), but this copvar does not have that form: it needs at
  most one spatial variable, and R blocks before L2 blocks (see is_opvar1d)."; notOpvar2d "'%s' is a
  field of the 2-D form R x L2[x] x L2[y] x L2[x,y] (fields R00 to R22), but this copvar does not have
  that form: ... (see is_opvar2d)."; notOpvar "'%s' is a field of the 1-D or the 2-D form, but this
  copvar has neither: ..."; opvarField (read) "A copvar has no field '%s'; it is computed from the
  blocks of the 1-D (P, Q1, Q2, R) or 2-D (R00 to R22) form. ..."; opvarField (set) "... setting it
  changes the blocks of that form. ..."; cdopvar likewise.

## C102 – FIXED (text)
- Files (MMP): `executives_sopvar/utility_functions/exec_sop_settings.m:131` "ST should be a settings
  struct (lpisettings)." -> "(lpisettings_sop).";
  `executives_sopvar/utility_functions/settings_2d_sop.m:74` "ST should be an lpisettings struct (field
  settings_2d) ..." -> "ST should be a struct of LPI settings (field settings_2d) or a 2-D settings
  struct.". Identifiers unchanged.
- Triggered before/after (`t1_class_*.txt`): `exec_sop_settings(5)`, `settings_2d_sop(struct('x',1))`;
  after: the new texts, identifiers exec_sop_settings:st, settings_2d_sop:st.

## C103 – FIXED (text)
- File (MMP): `executives_sopvar/utility_functions/exec_sop_settings.m:139, 145, 160`: " on 11/12/2026"
  removed from both 'retired' messages; "(settings report 11/16/2026, Sec. 5)." -> "(see
  lpisettings_sop)." (its help documents the storage row, `settings_sopvar/lpisettings_sop.m:53-56`).
  The field named and the replacement stay in each message; format arguments unchanged.
- Triggered after (`t1_class_after.txt`): "settings.dd1 was retired; use lpi.pos."; "settings.sosineq_on
  was retired; remove it: ..."; "settings_nd.LF_int, LF_mult and Zop_deg were retired; use entry 3 of
  lpi.pos and lpi.var."; "lpi.pos is a storage row (len_com): its fields mult, int, joint, sep are not
  read; set len_inc, len_com, len_sos and cap (see lpisettings_sop)." – %s filled, no specifier left.
- This file has mixed line endings (100 CRLF of 204 lines): the Edit tool rewrote it as LF, so the file
  was rebuilt from the original copy by a byte-exact literal replacement (`<S>/lit_replace.pl`, pairs in
  `<S>/work/`), checked with `diff` against the original (4 lines differ, terminators kept) and copied
  back; CR count 100 and `w/mixed` restored. The same procedure was used for the other two mixed files
  (C104 copquadvar.m, C105 poslpivar_settings_2d_sop.m), without the Edit tool.

## C104 – FIXED (text)
- Files (MMP): `lpi_programming_sopvar/poscopvar.m:258`, `lpi_programming_sopvar/misc/copquadvar.m:287`:
  the last part of the 'psatz' error, "2 is poslpivar_2d's ball, not the box." -> "2 is a ball weight,
  not supported (1 is the box weight)." (no retired name, no possessive).
- Triggered after: `poscopvar(prog,T,1,struct('psatz',2))` and
  `copquadvar(prog,1,{{'s'}},[0 1],1,struct('psatz',2))` give "'psatz' should be 0, 1, or a face code
  2d+1 / 2d+2 with d in 1..nv, here nv = 1; 2 is a ball weight, not supported (1 is the box weight)."

## C105 – FIXED (text)
- Files (MMP): `executives_sopvar/utility_functions/poslpivar_sop.m:129` "poslpivar degrees describe one
  spatial variable per space; ..." -> "The degrees D = {d1,d2,d3} describe one spatial variable per
  space; a space over %d variables has no such form." (%d argument unchanged);
  `executives_sopvar/utility_functions/poslpivar_settings_2d_sop.m:127` "use_sosineq selects
  lpi_ineq_2d, which has no container counterpart." -> "use_sosineq = 1 (an SOS inequality in place of
  the slack) is not supported; set use_sosineq = 0." Text only; the psatz mapping of that file is
  untouched.
- Triggered after: `poslpivar_sop(prog,T2,1,struct())` with a 2-D T2 -> "... a space over 2 variables
  has no such form."; `poslpivar_settings_2d_sop(prog,T2,st,'slack','out',PIE2)` with
  settings_2d.use_sosineq = 1 -> the new text.

## C107 – FIXED (text)
- File (MMP): `sopvar/misc/dpvar_op_copvar.m:149, 165`: "(legacy: nnz(dim(:,1)) > 1)" and
  "(legacy: nnz(dim(:,2)) > 1)" removed.
- Triggered after: `dpvar({'d1','d2'})*X` with X: L2[s] -> R x L2[s], and `X'*dpvar({'d1';'d2'})` ->
  "A matrix times an operator whose output spaces differ is ambiguous; build the multiplier with
  mat2copvar_sop instead." / "An operator whose input spaces differ times a matrix is ambiguous; ...".

## C108 – FIXED (text)
- File (MMP): `sopvar/misc/mat2copvar_grid.m:164`: "refused here as mat2opvar refuses Q1." -> "only
  multipliers are accepted here." (callers with mult_only true: mat2copvar_sop with options.mult_only,
  the sum of a matrix and an operator in dpvar_op_copvar, eye_copvar_sop).
- Triggered after: `mat2copvar_sop(1,struct('out',1,'in',1),struct('out',{{{}}},'in',{{{'s'}}}),[0 1],
  struct('mult_only',true))` -> "Block (1,1) maps out of a space with variable(s) {s} the output space
  lacks, so a constant block there is an integral, not a multiplier; only multipliers are accepted
  here." (%d, %d, %s filled).

## C109 – FIXED (text)
- File (MMP): `sopvar/misc/cat_copvar_grid.m:107, 119`: "e.g. with opvar2copvar." -> "e.g. with
  mat2copvar_sop." (the constructor of a container from a matrix, the operand this error meets in
  practice).
- Triggered after: `[T1, ones(1,1)]` and `[Pc1, ones(1,1)]` -> "Cannot horzcat a 'copvar' ('cdopvar')
  with a 'double'. Operands must be copvar/cdopvar containers or sopvar/sdopvar blocks; convert others
  first, e.g. with mat2copvar_sop."

## After the edits – verification (measured)
- C1 (`t1_class_after.txt`): `Ps(1,1,end)` -> `sopvar:endIndex` "A sopvar is indexed as P(I,J), or P(I)
  for one row or column."; `Pd(1,1,end)` -> `sdopvar:endIndex` "An sdopvar is indexed as ...".
  Unchanged branches: `P(end,1)`, `P(1,end)` give 1 x 1 slices; on a one-column operator `P(end)` equals
  `P(2,1)` (isequaln 1), sopvar and sdopvar.
- C3, C4 (reference: the 'before' run, same seeds; `t_cmp.txt`): 40 random sopvar adjoints (1 to 3
  variables, every S1/S2/S3 role), 30 sdopvar*sopvar and 30 sopvar*sdopvar products: all bitwise equal
  (isequaln), max abs difference 0. Testfolder tests (`t2_after_<k>.txt`): Test_copvar 31 passed
  (adjoint property <y,Ax> = <A*y,x>, err 2.27e-13; products against composed actions), Test_cdopvar 23
  passed (cdopvar*copvar, copvar'), test_sdopvar all passed, Test_mtimes_sopvar composition err 0
  (needs `testing_codes/private/misc` for rand_poly); outputs identical to 'before' line by line.
  Chain (sopvar -> LPI -> executive): PIETOOLS_stability_sop on the 1-D heat equation x_t = x_ss + 5x,
  Dirichlet, 'light', forms 'P' and 'Q' (0.9 s, 0.5 s; maxre -4.8696 = -(pi^2-5)), and on the 2-D heat
  equation x_t = x_s1s1 + x_s2s2, Dirichlet, 'stripped' with the settings of the 2-D reaction-diffusion
  example, form 'P' (0.5 s): status optimal before and after; solution vector RRx, SDP shape and every
  block of the solved P bitwise equal, max abs difference 0. The 1-D H-infinity estimator: gam
  1.002308616 before and after, RRx bitwise equal.
- C6 (`t3_exec_after.txt`): `{'stability'}` -> runs, status optimal, form Q; `{'stability_dual'}` ->
  runs, status optimal, form Qd; the flag names still run (Q, Qd).
- C8 (`t4_est2d_after.txt`): with the qp_copvar stub the executive prints "gain extraction: qp_copvar
  disabled by the test stub; inverse (inv)." and returns L, a copvar, with info.inv fields method, d,
  relrms, condIKH, Kdiff, nF, note; that L equals `getObserver_sop(P,Z,struct('method','inv'))` called
  directly on the P, Z of the run without the stub (P, Z of the two runs equal; isequaln 1). Without the
  stub: unchanged (L copvar, info.inv {qp}, gam 0.00040697058 in both runs and both states).
- C9, C10: the executive runs above pass through both edited lines (lpi_solve_loop_sop calls
  lpi_classify_sop, which reads gam of the estimator; pie_witness_sop calls pie_disc_sop, which gives
  maxre); results bitwise equal to 'before'.
- Line endings: `git ls-files --eol` gives the same `w/` state for all 25 candidate files before and
  after (`<S>/eol_before.txt`, `<S>/eol_after.txt`); CR counts change only by the CRLF lines added or
  deleted in CRLF files (`<S>/crlf_counts_before.txt`, `<S>/crlf_counts_after.txt`).

## Parse check (measured)
- `syntax_check` of the release pass over the 21 edited files (list `<S>/edited_files.txt`, output
  `<S>/syntax_check_out.txt`): "21 files checked, 0 with parse problems".
- checkcode before (original copies) and after (`<S>/checkcode_diff.txt`): no new message in any file;
  removed: IFBDUP (identical if/else branches) in lpi_classify_sop.m (C9), NASGU (value never used,
  S3Vardummy) in sopvar/@sopvar/ctranspose.m (C4), one message fewer in sopvar/@sdopvar/mtimes.m (C3).

## Summary
| Item | Status | Files |
|---|---|---|
| C1 | FIXED | sopvar/@sopvar/end.m, sopvar/@sdopvar/end.m |
| C3 | FIXED | sopvar/@sdopvar/mtimes.m |
| C4 | FIXED | sopvar/@sopvar/ctranspose.m |
| C5 | FIXED (text) | sopvar/@sdopvar/sdopvar.m |
| C6 | FIXED | executives_sopvar/PIETOOLS_auto_execute_sop.m |
| C8 | FIXED (fallback 'inv'; alternative: delete the retry) | executives_sopvar/PIETOOLS_Hinf_estimator_sop.m |
| C9 | FIXED | lpi_programming_sopvar/misc/lpi_classify_sop.m |
| C10 | FIXED | executives_sopvar/utility_functions/pie_disc_sop.m |
| C11 | FIXED (help text) | executives_sopvar/PIETOOLS_H2_control_sop.m, PIETOOLS_H2_estimator_sop.m, PIETOOLS_Hinf_control_sop.m, PIETOOLS_Hinf_estimator_sop.m |
| C101 | HELD (message lines by DJ, `git blame`) | none |
| C102 | FIXED (text) | executives_sopvar/utility_functions/exec_sop_settings.m, settings_2d_sop.m |
| C103 | FIXED (text) | executives_sopvar/utility_functions/exec_sop_settings.m |
| C104 | FIXED (text) | lpi_programming_sopvar/poscopvar.m, lpi_programming_sopvar/misc/copquadvar.m |
| C105 | FIXED (text) | executives_sopvar/utility_functions/poslpivar_sop.m, poslpivar_settings_2d_sop.m |
| C107 | FIXED (text) | sopvar/misc/dpvar_op_copvar.m |
| C108 | FIXED (text) | sopvar/misc/mat2copvar_grid.m |
| C109 | FIXED (text) | sopvar/misc/cat_copvar_grid.m |
