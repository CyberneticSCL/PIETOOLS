# G2_conv_root – ledger of the comment pass of 11/21/2026

Agent G2_conv_root. List `G2_conv_root.tsv` (38 files). Procedure and tools: `BRIEF.md`.

Conventions of this group:
- TAG authors: everyone whose change since the 2025 release persists in the file (git blame of the
  lines not in the 2025 version, `git log --follow` for deletions), first change first. A commit that
  only renamed or moved the file (similarity 100 %, e.g. a222bbc7 of 11/16/2026) does not count as a
  change: prep proposes DJ for such files; the TAG here is MMP alone (open question Q1).
- Files with mixed CRLF/LF line endings (`w/mixed`: the four lpi_programming_sopvar files and
  closedLoopPIE.m): the Edit tool rewrites every line ending to CRLF (tested on a scratch copy), which
  changes the `w/` state. Their header edits use `G2/blockrep.pl` (scratch folder): replaces one exact
  block of lines, every other line keeps its bytes and line ending; fails unless the block occurs once.
- Prep run from a private copy `G2/cpass_g2.pl` that writes its diff files to `G2/` instead of the shared
  folder (the other agents run prep at the same time); otherwise identical to `cpass.pl`.

## Files

### Comment-only files (code identical to 2025): 13 files, class condense

Each: TAG `MMP, 11/21/2026`; note `% MMP, 11/21/2026: Changes since PIETOOLS 2025: help and comments
revised.` in place of the new header entries listed; apply: 0 (was) dropped, 0 kept, 0 retagged; new
comment lines read (help text and step comments of the release comment pass): none removed; check PASS,
codeident OK, EOL unchanged, numstat = numstat -w.

| File | New header entries replaced |
|---|---|
| GUI/getVarList.m | 11/19 release pass |
| GUI/isvalidPolynomial.m | 11/19 release pass |
| GUI/isvalidTerm.m | 11/19 release pass |
| converters/convert_PIETOOLS_DDE.m | 11/19 release pass |
| converters/convert_PIETOOLS_NDS.m | 11/19 release pass; 11/20 help entry (field Cvd of the DDF output; the code change is in convert_PIETOOLS_NDS2DDF) |
| converters/initialize_PIETOOLS_PDE.m | 11/19 release pass |
| converters/legacy/2D/PIETOOLS_FTC_Expansion_2D.m | 11/19 release pass; prep proposed `DJ, MMP` (DJ: rename only, a222bbc7) – Q1 |
| converters/legacy/convert_PIETOOLS_PDE_terms_internal.m | 11/19 release pass; DJ rename only – Q1 |
| converters/legacy/initialize_PIETOOLS_PDE_batch.m | 11/19 release pass; DJ rename only – Q1 |
| converters/legacy/var_swap.m | 11/19 release pass; DJ rename only – Q1 |
| converters/minimize_PIETOOLS_DDF.m | 11/19 release pass |
| lpi_programming_sopvar/lpidecvar.m | 11/04 entry (relocation, code unchanged); 11/19 release pass; w/mixed: blockrep |
| lpi_programming_sopvar/lpisetobj.m | 11/04 entry `MMP - 11/04/2026` (relocation; prep did not list it as an entry, its format is not matched); 11/19 release pass; w/mixed: blockrep |

### pietools_path_update.m – condense
- TAG `MMP, 11/21/2026`. Note (replaces the 11/19 release-pass entry and the 11/21 fix entry):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: adds the folder of this file`
  `%   and its subfolders to the path, not the current folder; help and comments`
  `%   revised.`
- apply: 0 (was) dropped, 1 kept (`current_folder = pwd;`, 2025 code), 2 markers retagged.
- Internal notes removed: none. Help line reworded (check flagged "session"): "the MATLAB path of the
  current session is extended" -> "the current MATLAB path is extended".
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### pietools_update.m – condense
- TAG `MMP, 11/21/2026`. Note (replaces the 11/19 release-pass entry and the 11/21 fix entry):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: startup.m is closed also when`
  `%   it already holds the line "tbxmanager restorepath"; help and comments`
  `%   revised.`
- apply: 0 (was) dropped, 1 kept (`fclose(fid);` inside the if, 2025 code), 1 marker retagged.
- Internal notes removed: none. Help line reworded (check flagged "sessions"): "the path with PIETOOLS
  is saved for future sessions" -> "the MATLAB path with PIETOOLS is saved".
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### converters/internal/poly_decomp_1D_econ.m – condense
- TAG `MMP, 11/21/2026` (prep proposed `DJ, MMP`; DJ: rename only, a222bbc7 – Q1). Note (replaces the
  11/19 release-pass entry and the 11/20 fix entry):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: a polynomial without variables`
  `%   is treated as a constant (H of class 'double'); the monomial basis of each`
  `%   column of Q is formed from its degree matrix, rows sorted.`
- apply: 0 (was) dropped, 3 kept (2025 code), 7 markers retagged.
- Internal notes rewritten: the 2-line note on the failure at `Q.Varname(1)` -> what the branch does
  (constant Q: Zth identity, H = Q as 'double'); the 3-line note naming the monomials function of
  multipoly/SOSTOOLS and its error message -> one line (basis of column i from its degmat, by degree).
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### lpi_programming_sopvar/lpisolve.m – condense (w/mixed: blockrep)
- TAG `MMP, 11/21/2026`. Note (replaces the entries of 11/04 relocation, 11/08 comment pass, 11/12
  settings, 11/18 empty solver – with its session narrative "Found by DEMO3 ..." – and 11/19 release pass):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: OPTS may be a settings`
  `%   structure (its field sos_opts is used); without OPTS, the options are`
  `%   prog.settings.sos_opts if present, else the solver defaults; an empty`
  `%   solver field is set to 'mosek' when mosekopt is on the path, else to`
  `%   'sedumi'; help and comments revised.`
- apply: 0 (was) dropped, 5 kept (the 2025 if/else around sossolve), 20 markers retagged.
- Internal notes removed: those of the header entries only; body step comments kept.
- check PASS, codeident OK, EOL same (mixed), numstat = numstat -w.

### converters/internal/convert_PIETOOLS_NDS2DDF.m – condense
- TAG `MMP, 11/21/2026`. Note (replaces the 11/19 release-pass entry and the two 11/20 fix entries):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: the distributed delay kernels`
  `%   are stored in the field Cvd of the output; an NDS without delays (tau`
  `%   empty) throws an error; help and comments revised.`
- apply: 0 (was) dropped, 2 kept (2025 code), 4 markers retagged.
- Internal notes rewritten: "No delays: not supported; an error (the DDE conversion recursed without
  end)." (2 lines, bug history) -> "No delays: not supported, an error".
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### converters/internal/minimize_PIETOOLS_DDE2DDF.m – condense
- TAG `MMP, 11/21/2026`. Note (replaces the 11/19 release-pass entry and the 11/20 fix entry; also
  covers the unlogged error check of a8637201, 08/14/2026, `if rT==0`):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: a DDE without delays (tau`
  `%   empty) or with a delay without coefficients throws an error; help and`
  `%   comments revised.`
- apply: 0 (was) dropped, 1 kept (2025 code), 2 markers retagged.
- Internal notes rewritten: "... an error, since convert_PIETOOLS_DDE calls this function again
  (unbounded recursion)." (2 lines) -> "No delays: not supported, an error".
- Open: the three lines of the `rT==0` check carry no marker (unlogged change of 08/14/2026); left so.
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### converters/legacy/2D/expand_PIETOOLS_PDE_BCs.m – condense
- TAG `DJ, MMP, 11/21/2026` (DJ: new error message, 77a4d086 of 02/21/2026, by blame; MMP: comments).
  Note (replaces the 11/19 release-pass entry):
  `% DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: a BC with an input that`
  `%   varies in space throws the error "Conversion to PIE is not supported for`
  `%   the proposed PDE."; help and comments revised.`
- apply: 0 / 0 / 0 (no markers). The two commented-out lines of the old error call (2025 code, no
  marker, DJ) kept as 2025 code.
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### converters/initialize_PIETOOLS_DDE.m – condense
- TAG `MMP, 11/21/2026`. Note (replaces the 11/19 release-pass entry and the 11/20 fix entry):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: the number of disturbances nw`
  `%   is also read off B1di; help and comments revised.`
- apply: 0 (was) dropped, 1 kept (the 2025 condition with 'Bdi'), 2 markers retagged.
- Internal notes rewritten: "B1di, the field name (Bdi is none)" -> "Any disturbance coefficient
  defined, B1di included".
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### converters/initialize_PIETOOLS_DDF.m – condense
- TAG `MMP, 11/21/2026`. Note (replaces the 11/19 release-pass entry and the 11/20 fix entry):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: nx can be read off B1; a`
  `%   coefficient with either dimension wrong throws the size error; help and`
  `%   comments revised.`
- apply: 0 (was) dropped, 15 kept (2025 code: `size(B1,1,'var')` and the 14 `size(X)~=[r,c]` tests),
  19 markers retagged (prep: 15 KEEP-WAS, 19 RETAG before the hand edits; 3 two-line comments merged).
- Internal notes rewritten: inline "'var' is no dimension argument (an error)" removed; "by isequal: if
  size(X)~=[r,c] is true only when both dimensions differ", "by isequal, as in Group 2" (3x) -> "in both
  dimensions".
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### converters/initialize_PIETOOLS_NDS.m – condense
- TAG `MMP, 11/21/2026`. Note (replaces the 11/19 release-pass entry and the 11/20 fix entry):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: A0 is set to zero when`
  `%   undefined and its size is checked; nx is also read off Ei, Edi, nw off`
  `%   B1di, nz off C1ei, C1dei, and ny off C2ei, C2dei; the errors for an`
  `%   element without a matching delay name the field (C1ei, C2ei, Edi, C1dei,`
  `%   C2dei); help and comments revised.`
- apply: 0 (was) dropped, 11 kept (2025 code), 28 markers retagged (prep: 29 RETAG; two comment lines
  merged into one).
- Internal notes rewritten: "X, Y included, else their branches below are never reached" (3x) and "B1di,
  the field name (Bdi is none)" -> "Any ... coefficient defined, X and Y included"; "A0: ..., as in
  initialize_PIETOOLS_DDE; without it NDS_out.A0=A0 failed on a missing A0" -> "A0: zero if undefined,
  else size checked (nx x nx)".
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### lpi_programming_sopvar/lpiprogram.m – condense (w/mixed: blockrep)
- TAG `MMP, DJ, 11/21/2026` (MMP from 10/01/2026, DJ 11/06/2026 cellstr names). Note (replaces the
  entries of 10/01, 11/04, DJ 11/06, 11/12 and the 11/19 release pass; they held a measured error
  message, a cleaning/ path, a maintainer quote and the history of a retired adapter):
  `% MMP, DJ, 11/21/2026: Changes since PIETOOLS 2025: optional last input`
  `%   SETTINGS (a settings structure or a tier name, default 'light'), stored`
  `%   in prog.settings, the other inputs declared by the new local function`
  `%   declare_program; any number of spatial variables, also none; variable`
  `%   names as cellstr accepted; each dummy variable is checked to be a`
  `%   polynomial variable; help and comments revised.`
- apply: 0 (was) dropped, 9 kept (2025 code), 20 markers retagged (DJ marker of 11/06 included).
- Internal notes removed: none in the body (step comments only).
- check PASS, codeident OK, EOL same (mixed), numstat = numstat -w.

### PIETOOLS_DDE.m – condense
- TAG `YP, DB, MMP, 11/21/2026` (YP 02/18/2026: `uinput.ic`; DB 05/14/2026: legacy converter name, second
  initialization commented out; MMP from 11/04/2026). Note (replaces the entries 11/04 (2), 11/18 (2),
  11/19 (2); they held measurements, D-numbers, DECISIONS, cleaning/ paths, maintainer quotes):
  `% YP, DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the LPIs are solved`
  `%   by the container executives, with the settings of the 'light' tier or`
  `%   the settings st of a library example that sets them; stability_form =`
  `%   'P' and Hinf_gain_form = 'Q' set the forms of the stability and`
  `%   H-infinity gain tests (the dual appends 'd'); the DDE is initialized`
  `%   once; the direct conversion calls convert_PIETOOLS_DDE2PIE_legacy; the`
  `%   PIESIM initial condition is uinput.ic; help and comments revised.`
- apply: 5 (was) dropped (`% DDE_minimal_rep=0;` and the four intermediate calls with 'Q'/'Qd'), 11 kept
  (2025 executive calls), 18 markers retagged (prep: 22 RETAG; 7 marked comment lines rewritten as 3).
- By hand: the marker removed from `DDE_minimal_rep=1;`, which equals the 2025 line (2025 had 1; the
  11/04 entry described a development state 0), so the line is the 2025 line again.
- Internal notes rewritten: "(they start from a tier; the example file names why)" -> one line; the
  5-line comment on the container executives (stock executives stop on a copvar, PIETOOLS_H2_controller
  "called before, does not exist") -> 2 lines (what the calls take and return).
- DB commented-out `initialize_PIETOOLS_DDE` call kept (2025 code, no marker).
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### PIETOOLS_DDF.m – condense
- TAG `MMP, 11/21/2026`. Note (replaces the entries 11/04 (2), 11/12, 11/18, 11/19 (3), 11/21; they held
  measurements, D-numbers, DECISIONS, cleaning/ paths, maintainer quotes):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: the LPIs are solved by the`
  `%   container executives, with the settings of the 'light' tier or the`
  `%   settings st of a library example that sets them; stability_form = 'P'`
  `%   and Hinf_gain_form = 'Q' set the forms of the stability and H-infinity`
  `%   gain tests (the dual appends 'd'); controller and estimator: separable`
  `%   storage, settings.lpi.pos.sep = 1, unless lpi.pos has a field len_com; a`
  `%   DDE (DDE2DDF = 1) is converted by convert_PIETOOLS_DDE2PIE_legacy or,`
  `%   with DDF_min = 1, through its minimal DDF; sosineq_on removed; help and`
  `%   comments revised.`
- apply: 6 (was) dropped (intermediate 'Q'/'Qd' calls, the intermediate `settings.lpi.pos.sep = 1;` x2),
  14 kept (2025 code), 22 markers retagged (prep: 27 RETAG; 10 marked comment lines rewritten as 5).
- Internal notes rewritten: "settings_2d_stock would" dropped; "(they start from a tier; the block names
  why)"; "the two conversions of PIETOOLS_DDE.m. Both called convert_PIETOOLS_DDE2PIE, which does not
  exist."; the 5-line container-executive comment ("called before, do not exist") -> 2 lines;
  "exec_sop_settings rejects it" (2x) -> "P separable for the synthesis; a storage row (len_com) takes no
  sep".
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### PIETOOLS_PDE.m – condense
- TAG `YP, MMP, DB, 11/21/2026` (YP 02/18-22/2026: one-vector `uinput.ic`, PIESIM outputs `primary`/`grid`,
  no `ndiff`, example index; MMP from 11/05/2026; DB 11/17/2026: example 17). Note (replaces the entries
  11/04, 11/05 (2), 11/19, 11/20; they held bug narratives and PIESIM test cases):
  `% YP, MMP, DB, 11/21/2026: Changes since PIETOOLS 2025: Step 1 loads library`
  `%   example 17; Step 3 takes the 'heavy' settings of the container`
  `%   executives and lists these executives; Step 4 passes the initial`
  `%   conditions as one vector uinput.ic, sized from PIE.space (for PIESIM of`
  `%   a pde_struct, in the order of PDE.x_tab), with a zero initial estimate`
  `%   per kind of space for an estimator, and calls PIESIM without ndiff;`
  `%   Step 4b plots the PIESIM output fields primary and grid, in 2-D the`
  `%   states in s1 and in (s1,s2) at the middle node in s2; help and comments`
  `%   revised.`
- apply: 22 (was) dropped, 3 kept (`lpisettings('heavy')`, the two `uinput.ic.ODE/PDE` estimator lines:
  2025 code), 52 markers retagged (prep: 60 RETAG; see below). The 22 dropped: 2 (was) lines of the
  2025 description (help text); 12 (was) lines of the 2025 commented template calls of Step 3
  (`% PIETOOLS_auto_execute`, `% [prog, P] = PIETOOLS_PIE2PDEstability(...)`, ...; comment lines in 2025,
  so not 2025 code for apply; they name executives retired from the release); 8 (was) lines of the
  YP versions of Step 4/4b (post-2025 code).
- By hand: markers removed from the 3 new help lines (the description of the script).
- Internal notes rewritten: Step 3 executive comment (stock executives "stop on it") 5 -> 4 lines; the
  4-line PIE.space comment (opvar fields, legacy format history) -> 2 lines; the 6-line x_tab comment
  (bug narrative) -> 4 lines.
- First apply attempt died at `open(... '>:raw')` (file briefly locked, likely OneDrive); nothing was
  written; the retry succeeded.
- check PASS, codeident OK, EOL same (lf), numstat = numstat -w.

### GUI/parsePDEParams.m – condense
- TAG `MMP, 11/21/2026`. Note (replaces the entries 11/05, 11/08, 11/19, 11/20; they held GUI test counts
  "33 of the 47", "46 of the 47", a commit hash):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: the loop over the orders of`
  `%   differentiability runs to max(Xvals(:,1)), so two or more PDE states`
  `%   parse; each PDE state gets the variables and domain of the PDE`
  `%   (PDE_struct.vars, PDE_struct.dom); an empty coefficient text gives 1,`
  `%   and a coefficient text that fails to evaluate, also with a leading '(',`
  `%   throws an error naming it; help and comments revised.`
- apply: 0 (was) dropped, 7 kept (2025 code), 9 markers retagged.
- BEGIN/END of eval_b (a 6-line BEGIN comment with the bug history) -> `% MMP, 11/21/2026: BEGIN retry
  with a leading '(', else an error` and the matching END line.
- Internal notes rewritten: inline "'' errored at exp(end)" -> "empty text: no exp(end)".
- check PASS, codeident OK, EOL same (lf), numstat = numstat -w.

### converters/convert_PIETOOLS_DDF.m – condense
- TAG `DJ, MMP, 11/21/2026` (DJ 10/29/2026: PIE.space; MMP 11/04/2026: containers). Note (replaces DJ
  10/29, MMP 11/04, the 11/19 release pass):
  `% DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: the twelve operators of`
  `%   the PIE are 'copvar' containers of 'sopvar' blocks, not 'opvar' objects;`
  `%   the spatial variable and domain are set in PIE.space (vars, dom); help`
  `%   and comments revised.`
- apply: 28 (was) dropped (the 2025 section comments of the opvar block, commented once more – not code
  – and the intermediate DJ line `PIE.space.vars = ...(Top.var1)...`), 97 kept (the 2025 opvar code), 19 markers
  retagged (prep: 24 RETAG; 8 marked comment lines rewritten as 3).
- BEGIN/END -> `% DJ, MMP, 11/21/2026: BEGIN the operators as 'copvar' containers` and the matching END.
- Internal notes rewritten: the 8-line description (former opvar, OP2COPVAR_SOP in INITIALIZE, "Former
  code ... (MP 10_01_2020; variables DJ 12/07/2024), dims (SS 9/28)") -> 3 lines (operator form and
  arguments); inline "the former opvar field I" and "s1 was Top.var1" removed. "% Former code:" kept as
  the head of the (was) block.
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### converters/legacy/convert_PIETOOLS_DDE2PIE_legacy.m – condense
- TAG `MMP, 11/21/2026` (prep proposed `DJ, MMP`; DJ: rename only – Q1). Note (replaces 11/04 and the
  11/19 release pass; the 11/04 entry held a library count "12 of the 26"):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: the twelve operators of the`
  `%   PIE are 'copvar' containers of 'sopvar' blocks, not 'opvar' objects;`
  `%   D12, D22 and Tu act on u, with nu columns, also when nu ~= nw; the`
  `%   spatial variable and domain are set in PIE.space and the PIE is`
  `%   initialized; help and comments revised.`
- apply: 29 (was) dropped (2025 section comments of the opvar block, not code), 97 kept (2025 opvar
  code), 22 markers retagged (prep: 27 RETAG; 9 marked comment lines rewritten as 4).
- BEGIN/END -> `% MMP, 11/21/2026: BEGIN/END the operators as 'copvar' containers`.
- Internal notes rewritten: the 9-line description (former opvar, OP2COPVAR_SOP in INITIALIZE, "the
  former dims gave nw ...", former code authors) -> 4 lines; inline "the former opvar field I" removed.
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### converters/legacy/convert_PIETOOLS_PDE_batch.m – condense
- TAG `MMP, 11/21/2026` (prep proposed `DJ, MMP`; DJ: rename only, and the DJ initials prep saw are in a
  2025 comment and in the 11/04 entry text – Q1). Note (replaces 11/04, which held a maintainer quote and
  the derivation history, and the 11/19 release pass; it states that the DJ 09/29/2021 change is no longer
  in the code, rule 4 of for_developers.txt):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: the twelve operators of the`
  `%   PIE are 'copvar' containers of 'sopvar' blocks, not 'opvar' objects: each`
  `%   is the product of a row slice of Apop and a column slice of Phfop (new`
  `%   local function mult_or_zero: the zero operator for an empty channel), or`
  `%   a column slice of Tbigop; the dummy variable is the spatial variable`
  `%   with '_dum'; the PIE is initialized; the opvar code, the zero filling of`
  `%   DJ, 09/29/2021 included, is commented out; help and comments revised.`
- apply: 0 (was) dropped, 1 kept (2025 code), 3 markers retagged.
- The 2025 code commented out without markers (Apop; Phfop, Ptop, Tbigop, op_slice, DJ zero filling)
  kept. Its 4 BEGIN/END pairs (2-4 line BEGIN comments) -> one line each: `% MMP, 11/21/2026: BEGIN/END
  the opvar Apop, commented out`, `... the opvars Phfop, Ptop, Tbigop and their slices, commented out`,
  `... the operators as 'copvar' containers`, `... new local function mult_or_zero`.
- By hand: the commented 2025 comment `% % DJ 09/29/2021: fill in empty parameters ...` read as a marker
  by apply/check (it would have become `% % MMP, 11/21/2026`); restored to its 2025 text `% DJ 09/29/2021:
  ...`, which the tools skip as a 2025 line.
- Internal notes removed: inline "Top has no var1, var2".
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### Rerun with the updated cpass.pl (coordinator fix of 06:14: "DJ (CC)" initials, markers without a
blank, per-process diff files)
- The 31 files finished before the fix: fresh prep code diffs identical to the first prep (the first prep
  ran from the private copy, free of the shared-file race); apply changed no file; check PASS, codeident
  OK, EOL unchanged, numstat = numstat -w for all 31. No note needed correction.
- converters/legacy/2D/Compute_Tmap_PDE_2D_Separable.m reclassified review by the coordinator (3 of 4 new
  entries DJ); it had not been edited.

### converters/legacy/convert_PIETOOLS_PDE_batch2terms_legacy.m – condense
- TAG `MMP, 11/21/2026` (prep proposed `DJ, MMP`; DJ: rename only – Q1). Note (replaces 11/05, which held a
  test count "51 batch PDEs, before: 11 differ" and the bug history, and the 11/19 release pass):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: the batch coefficients are`
  `%   rewritten in s1, s1_dum; the columns of E0, C10 and C20 are reordered`
  `%   from the batch order of the boundary values, [x1(a) x1(b) x2(a) x2(b)`
  `%   x2s(a) x2s(b)], to [x1(a) x2(a) x2s(a) x1(b) x2(b) x2s(b)]; help and`
  `%   comments revised.`
- apply: 0 (was) dropped, 3 kept (2025 code), 3 markers retagged.
- BEGIN/END: two pairs with 2- and 4-line BEGIN comments -> `% MMP, 11/21/2026: BEGIN/END the
  coefficients, polynomials in PDE.vars, in s1, s1_dum` and `... columns of E0, C10, C20: batch order ->
  order in r`.
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### converters/legacy/convert_PIETOOLS_PDE_terms.m – condense
- TAG `MMP, 11/21/2026` (prep proposed `DJ, MMP`; DJ: rename only – Q1). Note (replaces 11/05, which held
  two GUI test cases with numbers, and the 11/19 release pass):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: the coefficients of the terms`
  `%   are rewritten from PDE.vars in s1, s1_dum (new local function`
  `%   terms_vars); the coefficients of terms of one equation on the same block`
  `%   (an ODE state, an input, a boundary value or derivative of a PDE state)`
  `%   add up; help and comments revised.`
- apply: 0 (was) dropped, 13 kept (2025 assignments), 18 markers retagged.
- BEGIN/END of the new local function -> `% MMP, 11/21/2026: BEGIN/END new local function terms_vars`.
- Internal notes rewritten: "(INITIALIZE_PIETOOLS_PDE_terms_legacy)" dropped from a step comment.
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w (no newline at end of file, as in HEAD).

### converters/legacy/convert_PIETOOLS_PDE_terms_legacy.m – condense
- TAG `MMP, 11/21/2026` (prep proposed `DJ, MMP`; DJ: rename only, DJ initials only in a 2025 comment and
  the 11/05 entry text – Q1). Note (replaces 11/05, which held a usage count and the derivation history,
  and the 11/19 release pass; it states that the DJ 09/29/2021 change is no longer in the code):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: the twelve operators of the`
  `%   PIE are 'copvar' containers of 'sopvar' blocks, not 'opvar' objects: each`
  `%   is P0 + L*M*V, with M the PDE subsystem, P0 a block of the ODE matrices,`
  `%   L = [Bxr, I], Dzr or Dyr, and V = [Cv, I], Dvw or Dvu (new local function`
  `%   pie_op: the zero operator for an empty channel); the PIE is initialized;`
  `%   the opvar code, the zero filling of DJ, 09/29/2021 included, is`
  `%   commented out; help and comments revised.`
- apply: 0 (was) dropped, 1 kept (2025 code), 3 markers retagged.
- The 2025 code commented out without markers kept; its 5 BEGIN/END pairs (2-line BEGIN comments) -> one
  line each (`... the empty opvars of the PIE, commented out`, `... the opvar construction of the
  operators, commented out`, `... the operators as 'copvar' containers`, `... the zero filling of DJ,
  09/29/2021, commented out`, `... new local function pie_op`).
- By hand: `% % DJ 09/29/2021: fill in ...` (read as a marker by the tools) restored to its 2025 text.
- Internal notes removed: inline "Top has no var1, var2".
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### converters/legacy/initialize_PIETOOLS_PDE_terms_legacy.m – condense
- TAG `MMP, 11/21/2026` (prep proposed `DJ, MMP`; DJ: rename only – Q1). Note (replaces the 11/19 release
  pass and the two 11/20 entries):
  `% MMP, 11/21/2026: Changes since PIETOOLS 2025: nx is read off ODE.A also when`
  `%   A is the only ODE field; nr is read off the rows of a Crp or Drb`
  `%   coefficient, not of the cell of terms; help and comments revised.`
- apply: 0 (was) dropped, 3 kept (2025 code), 9 markers retagged.
- Internal notes rewritten: "(it got nx = 0 and was replaced by zeros(0))" and "(the rows of the cell
  array were read)" removed from two marked step comments.
- check PASS, codeident OK, EOL same (crlf), numstat = numstat -w.

### executives_sopvar/utility_functions/closedLoopPIE.m – asnew (18 % of the 2025 code left)
- Removed: all log entries after `% Initial coding SS - 5/30/2021` (29 lines: the SS 1/18/2024 and DJ
  12/07/2024, 01/04/2025 entries of 2025, DJ 10/29/2026, DB 11/03-11/04/2026 (4), MMP 11/08, 11/09, 11/19);
  the 9 (was) lines (opvar/op2copvar_sop branches); the dated note line `% 11/04/2026: xhat row of PIEcl.Bu
  ...`; 60 trailing dated markers (DB, MMP; in the help, the INPUT/OUTPUT text and the code), the inline
  DB notes on the observer rows included. Kept: license block with its "If you modify this code" lines,
  help, "For support" line, Initial coding line.
- Facts of log entries restated in the help (release rule 4): the dynamic closed loop requires Tw = 0 and
  Tu = 0 besides D22 = 0 (the code errors otherwise); "PIEs with an algebraic constraint (K, Kw, Ku) are
  not supported." (DJ 10/29/2026; the code errors). The observer sign and u-path conventions of the DB
  entries are in the help formula already.
- w/mixed: edited by the script `G2/asnew_closedLoopPIE.pl` (per-line expected text checked, each line
  keeps its line ending; dry run diffed before installing).
- codeident OK, EOL same (mixed), numstat = numstat -w, no dated marker, (was) line or log entry left
  except the Initial coding line.

### converters/convert_PIETOOLS_PDE.m – review (not edited)
- New notes: DJ 8 (DJ 11/02, 11/06 x2, 11/08, 11/11, DJ (CC) 11/13, 11/15, 11/20), DB 2 (11/18 entry;
  the body note `% DB 11/07/2026: outputs flagged as static constraints ...` at line 139), MMP 1 (11/19
  release pass). 12 % of the 2025 code left: by the 25 % rule the file would be asnew; the review class
  (notes mainly DJ) takes precedence. Decision for the maintainer: condense (draft below) or asnew.
- Proposed TAG `DJ, DB, MMP, 11/21/2026` (blame of the new lines: DJ, DB, MMP; YP commit 98d89c7c carries
  the same change as DJ d7fd1c0e).
- Draft note (in place of the 10 header entries, lines 90-104):
  `% DJ, DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: one converter for any`
  `%   number of spatial variables; linear PDEs give a 'pie_struct' with`
  `%   'copvar' operators, nonlinear PDEs (or the flag 'nonlinear') an`
  `%   'npie_struct' with 'tensopmat' operators; constraints K*xf + Kw*wf +`
  `%   Ku*uf = 0 (static constraints PDE.constraint_ID, conditions that are no`
  `%   BC along one direction, states without dynamics); delays and higher time`
  `%   derivatives as additional states; the PIE describes the PDE by space,`
  `%   space_pde, pde2pie and pie2pde; comp_order is no longer supported; the`
  `%   summary is printed by print_convert_summary; help and comments revised.`
- Prep verdicts (updated tool): 11 RETAG (markers DJ 11/02 (CC), 11/08, 11/11, DJ (CC) 11/13, 11/15,
  11/20 x5), 0 KEEP-WAS, 0 DROP-WAS, 0 REVIEW-BEGIN/END; the body note of line 139 (DB 11/07/2026) is
  read as a header entry: as a step comment it would lose its initials and date. 180 new unmarked comment
  lines (help and step comments; a scan for internal wording found none).
- If condensed: apply with the TAG retags the 11 markers; codeident, EOL (crlf) and numstat checks as
  for the other files.

### converters/legacy/2D/Compute_Tmap_PDE_2D_Separable.m – review (reclassified by the coordinator; not edited)
- New entries: DJ 02/15/2026, DJ (CC) 10/29/2026, DJ (CC) 11/02/2026, MMP 11/19/2026 (release pass);
  87 % of the 2025 code left.
- Proposed TAG `DJ, DB, MMP, 11/21/2026` (blame of the new lines: DJ from the 02/2026 merges, DB 10/29/2026
  one code and one comment line, MMP comments; prep proposed `DJ, MMP, DB` from the commit order without
  --follow).
- Draft note (in place of the 4 entries, lines 151-157):
  `% DJ, DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: outputs Top_f and`
  `%   Dvals (structs, the maps and derivative orders of x, w and u) instead of`
  `%   Top, Twop, Tuop, Dvals_w, Dvals_u; optional input suppress_summary; the`
  `%   boundary inputs are differentiated in every direction they need, with`
  `%   boundary values of the inputs added (new local function expand_input);`
  `%   the 1-D sub-PDEs are converted by convert_PIETOOLS_PDE_2D; periodic BCs`
  `%   throw an error; help and comments revised.`
- Prep verdicts (updated tool): 8 RETAG (DJ 02/15/2026 x3, DJ 02/01/2026 x1 – a date that differs from
  the header entry, DJ (CC) 11/02/2026, DB 10/29/2026, DJ (CC) 10/29/2026 x2), 1 KEEP-WAS (line 521, the
  2025 call of convert_PIETOOLS_PDE; its marker also reads "earlier: DJ, 12/16/2..."), 1 DROP-WAS (line
  529, `if ~isempty(PIE_1D.F)`, not 2025 code), 0 REVIEW-BEGIN/END.

## Summary

- condense: 35 edited (13 comment-only, 22 with code changes since 2025); asnew: 1 edited
  (closedLoopPIE.m); review: 2 not edited, drafts above (convert_PIETOOLS_PDE.m; Compute_Tmap_PDE_2D_
  Separable.m, reclassified by the coordinator before any edit); deferred/skipped: 0 (all 38 files were
  clean in `git status` before their first edit; at the end exactly the 36 edited files are modified).
- Final checks on the 36 edited files: check PASS (35 condense), codeident OK (36), `w/` EOL state as at
  the start (36), numstat = numstat -w (36), one note of at most 9 lines per condense file, no
  "Release comment pass" entry left; MATLAB parse check (lock runner, syntax_check): 36 files, 0 parse
  problems.

## Open questions

- Q1. TAG of files whose only DJ commit is the rename a222bbc7 (similarity 100 %): TAG `MMP, 11/21/2026`
  in 11 files where prep proposed `DJ, MMP` (PIETOOLS_FTC_Expansion_2D, convert_PIETOOLS_PDE_terms_
  internal, initialize_PIETOOLS_PDE_batch, var_swap, poly_decomp_1D_econ, convert_PIETOOLS_DDE2PIE_legacy,
  convert_PIETOOLS_PDE_batch, convert_PIETOOLS_PDE_batch2terms_legacy, convert_PIETOOLS_PDE_terms,
  convert_PIETOOLS_PDE_terms_legacy, initialize_PIETOOLS_PDE_terms_legacy). Adding DJ is a retag of the
  note and the markers (apply with the new TAG).
- Q2. PIETOOLS_PDE.m: apply dropped 12 (was) lines that kept the 2025 commented template calls of Step 3
  (`% PIETOOLS_auto_execute`, `% [prog, P] = PIETOOLS_PIE2PDEstability(...)`, ...): comment lines in 2025,
  hence not 2025 code for the tool; they name executives retired from the release. Restorable from HEAD.
- Q3. convert_PIETOOLS_PDE.m (review) has 12 % of its 2025 code: asnew treatment or the condensed note
  drafted above.
- Q4. Mixed CRLF/LF files (4 lpi_programming_sopvar files, closedLoopPIE.m) were edited by scripts that
  keep each line ending (G2/blockrep.pl, G2/asnew_closedLoopPIE.pl), not by the Edit tool, which rewrites
  every ending to CRLF and changes the `w/` state.
- Notes: minimize_PIETOOLS_DDE2DDF.m carries an unlogged, unmarked change of 08/14/2026 (rT==0 error),
  now in the note, no markers added. In convert_PIETOOLS_PDE_batch.m and convert_PIETOOLS_PDE_terms_
  legacy.m the tools read the commented 2025 comment `% % DJ 09/29/2021: ...` as a marker; it was restored
  to its 2025 text. Two help lines (pietools_path_update.m, pietools_update.m) were reworded only because
  check flags the word "session".
