# B_conv ledger – bug-fix pass 11/20/2026

Agent B_conv. Items: C12, C15, C16, C17, C32 (legacy parser), C21, C22, C23, C44 (no caller), C35, C36,
C38, C39, C40, C46; follow-up: C112, C113 and two cleanups. HEAD 17d68747. Scratch: `/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/B_conv/`.

Annotation convention used for header documentation lines (help text and the description block
before the licence): edited in place without per-line markers or (was) lines, the change stated in
the header log line – the header rule of 10/28/2026 (history only in the log). Code and body comments
carry the col-77 markers. (Measured: in files that start with a %%%% banner, `help` prints the banner
only, e.g. initialize_PIETOOLS_PDE_terms_legacy, `B_conv/help_legacy.txt`.)

DJ hold (coordinator clarification): every changed line that plain `git blame` credits to Declan traces,
by `git blame -w -C -C -M` on the block, to the bulk commits e15f2543 ("PIETOOLS 2021b Release.",
content of the 2021b import) or 4cefdf8f ("Move old initializers/converters to legacy"); none to a DJ
feature or bug-fix commit. No item held; nothing restored. Lines of DJ commits next to the edits (the nds
branch of @sys/convert, c88f2787; `if ~suppress`, 3a997fb3; the 12/07/2024 default-variable lines) are
untouched.

## C35 – FIXED
- Files: `converters/internal/convert_PIETOOLS_NDS2DDF.m` (help lines 14-16 and 37-39 corrected;
  log lines 77-78; body: comment line 116, (was) line 117, `DDF.Cvd{i}=[...` line 118);
  `converters/convert_PIETOOLS_NDS.m` (help lines 45-47: the DDF now has Cvd; log lines 77-78; help
  lines by MMP, comment pass 11/20, per git blame – the file is credited DJ, these lines are not).
  Cited lines 112-114 are 2021b import (e15f2543), not a DJ commit after 12/31/2021.
- Reproduction (measured): `B_conv/c35_repro.m`, output `B_conv/c35_before.txt`. Five systems: S1 (1
  delay, constant Adi), S2c (2 delays, constant kernels in every block x,w,u -> x,z,y), S2p (S2c with
  kernels polynomial in s1), S3c/S3p (NDS only: distributed neutral Edi, C1dei). Routes: NDSmin =
  convert_PIETOOLS_NDS(.,'pie'); NDSraw = convert_PIETOOLS_DDF(convert_PIETOOLS_NDS2DDF(.)) (the
  PIETOOLS_DDF.m path with DDF_min = 0); DDE = convert_PIETOOLS_DDE(.,'pie'). Each PIE compared on its
  transfer matrix [z;y]<-[w;u] at 4 complex points and its rightmost spectrum (Galerkin, pie2piesim +
  pie_disc_sop, Legendre N = 16 and 28) against the closed-form characteristic matrix
  Delta(l) = l(I - sum Ei e^{-l tau} - int Edi e^{ls}) - A0 - sum Ai e^{-l tau} - int Adi e^{ls} and
  transfer matrix of the delay system (Gauss-Legendre, 40 nodes).
  Before: the DDF from the NDS has no field Cvd; NDS routes differ from the closed form by
  8.4e-2 (S1), 4.3e-1 (S2c), 7.3e-1 (S2p), 2.1e-1 (S3c), 1.4e-1 (S3p); their eigenvalues are not roots
  of det Delta (sigma_min/sigma_max of Delta 0.03-0.6); the DDE route matches the closed form (1e-15),
  so max|G_NDS - G_DDE| = 8.4e-2 (S1), 4.3e-1 (S2c).
- Verification (measured, `B_conv/c35_after.txt`): field Cvd present; S1, S2c: NDSmin, NDSraw, DDE agree
  pairwise to <= 1.5e-15 and with the closed form to <= 1.2e-15, eigenvalues roots of det Delta
  (sigma ratio <= 6e-14), identical at N = 16 and 28; S2p, S3p (polynomial kernels): NDSraw matches the
  closed form to 1.8e-15; S3c: NDSmin and NDSraw match to 1.4e-15. Reference: closed form (independent
  of every converter) and the DDE route (does not use NDS2DDF).
- Finding on the polynomial-kernel routes (measured): `converters/internal/poly_decomp_1D_econ.m:75`
  called `monomials(Qt(:,i))` with one argument, so every 'polynomial' kernel errored in
  minimize_PIETOOLS_DDE2DDF and minimize_PIETOOLS_DDF – assigned as C112 and FIXED (below); with it the
  NDS route converts polynomial distributed kernels correctly as well.
- Not changed: `DDF_out=DDF;` (last line, unused, no effect).

## C38 – FIXED (reproduced on a narrower path than stated)
- Files: `converters/internal/minimize_PIETOOLS_DDE2DDF.m` (log line 68; comment lines 102-103, (was)
  line 105, error line 106); `converters/internal/convert_PIETOOLS_NDS2DDF.m` (log line 79; comment
  lines 101-102, (was) line 104, error line 105). Cited lines are 2021b import (the nK==0 branch with
  "%%% ADD case logic for when there are no delays", e15f2543); the recursion arose when DJ renamed the
  2021b DDE-to-PIE converter (5da2a25f, 09/29/2022, comment "This makes no sense..."), which did not
  write these lines.
- What convert_PIETOOLS_DDE does with a delay-free system (measured, `B_conv/c38_repro.m`,
  `c38_stack.m`): it errors in initialize_PIETOOLS_DDE – with tau missing, the explicit "you forgot to
  define the delays" error; with tau = [], "Unrecognized function or variable 'Ai'"
  (initialize_PIETOOLS_DDE line 772, initialize_PIETOOLS_NDS line 908: the K-loops create no cells for
  K = 0). Same for convert_PIETOOLS_NDS and convert_PIETOOLS_NDS2DDF. Hence the fix: an error in place
  of the recursive call, in both branches.
- Reproduction (measured): the nK==0 branches are reached only when tau = [], every cell field is given
  (as {}) and the dimensions nx, nw, nu, nz, ny are given (`B_conv/c38_allcells.m`: cells without
  dimensions fail with "Array indices must be positive integers" in the dimension scan;
  `c38_allcellsdims.m`): then convert_PIETOOLS_DDE, minimize_PIETOOLS_DDE2DDF, convert_PIETOOLS_NDS,
  convert_PIETOOLS_NDS2DDF all end in MATLAB:lang:StackOverflow "Out of memory. The likely cause is an
  infinite recursion", stack depth 58247 (minimize < convert_PIETOOLS_DDE < minimize ...).
- Verification (measured, `c38_after.txt`, `c38_allcellsdims_after.txt`): the same four calls stop at
  once with "The DDE has no delays (tau is empty). Define at least one delay." (NDS: "The NDS has no
  delays ..."), stack depth 6-7; a DDE and an NDS with one delay still convert (pie_struct, copvar T).
- Remaining (measured, not changed): tau = [] without cell fields still errors in the initializers with
  "Unrecognized function or variable 'Ai'"; a clear message there would need a K = 0 check in
  initialize_PIETOOLS_DDE/NDS (outside this item).

## C39 – FIXED
- File: `converters/initialize_PIETOOLS_NDS.m` (code lines from the 2021b import, e15f2543). Help lines
  16-25: A0 no longer described as required/unchecked; log lines 57-58; numel messages: (was)/new lines
  168-190; nx scan: comment 206, (was) 207, condition 208 (+ Ei, Edi), comment 216-217; A0 block 243-250
  (default zeros(nx) with a message, size check), copied from initialize_PIETOOLS_DDE (message names
  "states" only); nw scan: comment 374, (was) 375, condition 376 ('B1di' for 'Bdi'); nz scan: comment
  439, (was) 440, condition 441 (+ C1ei, C1dei), comment 449-450; ny scan: comment 478, (was) 479,
  condition 480 (+ C2ei, C2dei), comment 488-489, (was)/new lines 503-504 and 506-507 (ny, not nz).
- Reproduction (measured, `B_conv/c39_repro.m`, `c39_before.txt`): no A0 -> "Unrecognized function or
  variable 'A0'"; A0 2x3 and A0 2x2 with nx = 3 accepted silently; states from Ei only or Edi only,
  outputs from C1ei/C1dei/C2ei/C2dei only, w from B1di only -> each errors "The dimensions of X{1} are
  inconsistent ..." (dimension read as 0, since the outer conditions omit these fields); numel errors
  for C1ei, C2ei, Edi, C1dei, C2dei all say "element of Ei".
- Verification (measured, `c39_after.txt`): no A0 -> A0 = zeros(2), passes; A0 2x3 and A0 2x2 with
  nx = 3 -> "The dimensions of A0 are inconsistent with the detected number of states"; Ei-, Edi-,
  C1ei-, C1dei-, C2ei-, C2dei-, B1di-only systems pass with nx = 2, nz = 1, ny = 1 (C2ei of 2 rows:
  ny = 2), nw = 1; each numel error names its field. Valid inputs (reference: HEAD): `B_conv/lib_regress.m`
  with HEAD copies of the four edited converter files first on the path (`lr_base.mat`) and with the
  working tree (`lr_after.mat`): 46 conversions (DDE01-26, NDS01-12 less NDS06 (fragment), DDF01-03 by
  the minimized and the direct route, the active blocks of examples_DDE_library_PIETOOLS and
  examples_NDSDDF_library_PIETOOLS) give bit-identical Galerkin matrices of all 12 operators
  (`lr_compare.txt`); DDE_RAND and NDS06 fail in the block itself, identically in both runs.
- Same 'Bdi' typo in `converters/initialize_PIETOOLS_DDE.m:296`: assigned as C113 and FIXED (below).

## C40 – FIXED
- File: `converters/initialize_PIETOOLS_DDF.m` (code lines from the 2021b import). Log lines 56-57;
  (was) 120 / new 121 (`nx=size(B1,1);`); the 14 checks `elseif size(X)~=[r,c]` -> `elseif
  ~isequal(size(X),[r,c])`, each with its (was) line (Cr, Br1, Br2, Drv, Cv, Cvd, C1, D11, D12, D1v, C2,
  D21, D22, D2v); the comments of groups 2-5 that stated "error only if both dimensions differ"
  rewritten (markers).
- Semantics: `if v` with a vector v is true only if every element is nonzero, so `size(X)~=[r,c]` errs
  only when both dimensions differ; `size(B1,1,'var')` passes a char as a dimension argument.
- Reproduction (measured, `B_conv/c40_repro.m`, `c40_before.txt`): a DDF with B1 and without A0 ->
  "Dimension argument must be a positive integer scalar within indexing range."; each of the 14 blocks
  given with one extra column passed initialize_PIETOOLS_DDF silently.
- Verification (measured, `c40_after.txt`): the DDF without A0 passes (nx = 2, A0 = zeros(2)); each of
  the 14 one-dimension-wrong blocks stops with "error in size of <name><i>"; the valid DDF passes. Valid
  inputs: the library regression of C39 (DDF01-03 by both routes and the active NDS/DDF library block)
  is bit-identical to HEAD.
- Not changed: the messages "error in size of Br11" (name and index run together), outside the item.

## C36 – FIXED (nx, nr); explicit-nx precedence not changed
- File: `converters/legacy/initialize_PIETOOLS_PDE_terms_legacy.m`. Log lines 112-113; comment lines
  158-161; (was) 164 / new 165 (`if ~isempty(logVal) || isfield(PDE.ODE,'A')`); comment 302-303;
  (was)/new 305-306 and 308-309 (`size(PDE.PDE.Crp{end}.coeff,1)`, `Drb{end}`), as nBC is read
  (`PDE.BC.Ebb{end}.coeff`). Blame: plain blame names DJ 4cefdf8f ("Move old initializers/converters
  to legacy", 2022-08-18) for the changed lines; `git blame -w -C -C -M` on the blocks traces them to
  the 2021b import e15f2543 (PIETOOLS_2021b/converters/initialize_PIETOOLS_PDE_terms.m, Initial coding
  MMP 5/29/2021, SS 6/1/2021), where `size(PDE.PDE.Crp,1)` and `if ~isempty(logVal)` already stand
  (lines 252, 131). The DJ lines 178-180 (`if ~suppress`) and the DJ 12/07/2024 lines are untouched.
  'A' was not added to the isdefined list: the list also holds the field names of PDE.PDE, whose term
  structure PDE.PDE.A every PDE has.
- Reproduction (measured, `B_conv/c36_repro.m`, `c36_before.txt`): heat equation on [0,1] (Dirichlet)
  plus ODE.A = [0.5 1; 0 -2] alone: nx = 0, "A has incorrect dimension. Defaulting to zero", the PIE has
  no finite states (T.dim [0 0;1 1]) and its rightmost eigenvalue is -9.8696 (stable), while the system
  has the eigenvalue 0.5 – a stability test would certify an unstable system. With PDE.n.nx = 2 the same
  system converts correctly (0.5, -2, -pi^2, -4 pi^2). r given only by Crp (2 rows) or only by Drb:
  nr = 1 (rows of the cell array) and the terms replaced by zeros with a notice; the PIE is unchanged,
  since nothing uses r when Bxr, Dzr, Dyr are absent (they define nr first).
- Verification (measured, `c36_after.txt`): A alone gives nx = 2, T.dim [2 2;1 1], eigenvalues 0.5, -2,
  -9.8696, -39.4784 = eig(ODE.A) and -(k pi)^2 (closed form), identical to the explicit-nx route; Crp-
  or Drb-only gives nr = 2 without notices; PIE spectrum unchanged.
- NEEDS DECISION (not changed): ODE.A also takes precedence over an explicit PDE.n.nx (documented in
  the comment); with an inconsistent nx this initializer zeroes the mismatched matrices with a notice in
  either order. An error on inconsistency would change the convention of the legacy format.

## C12 – FIXED
- File: `GUI/parsePDEParams.m` (w/lf, kept). eval_b code lines from the 2021b import (e15f2543 by
  `git blame -w -C -C -M`; no DJ commit after 12/31/2021 wrote them; `git log -L` lists DJ fbb76aee,
  2022-08-18, only through whole-file hunks). Log lines 59-60; eval_b help 681-685 corrected; (was) 695
  / new 696 (`if ~isempty(exp) && strcmp(exp(end),'*')`); BEGIN/END block 705-721: (was) lines of the
  if/else, new try/catch retry and error.
- What the exact message text catches: the retry prepends '(' – a text with an unmatched ')'; the text
  names mismatched delimiters (inferred: an older MATLAB wording of that parse error). Measured
  (`B_conv/c12_census.m`, `c12_census.txt`): in R2025b eval('2*s)') fails with
  MATLAB:m_improper_grouping "Error: Unmatched ')'. Check for missing '('.", not that text, so the
  retry never runs and such a text also became 0. In the 36 GUI files of the example library, 318
  eval_b calls: 162 evaluate, 156 are '*'/'+' (1), none retries, none fails, none is empty.
- Reproduction (measured, HEAD copy first on the path, `c12_base.txt`): transport GUI file with the
  coefficient text replaced: 'k*' (k undefined) -> parsed without notice, the X1 equation left with 0
  terms (coefficient 0, term dropped: X1_t = 0); '' (term without '*') -> "Array indices must be
  positive integers or logical values."; '2)*' -> 0 terms (the retry not taken in R2025b).
- Verification (measured, `c12_after.txt`, `c12_compare.txt`): 'k*' -> "The coefficient 'k' of a term
  could not be evaluated: Unrecognized function or variable 'k'."; '' -> 1 term, coefficient 1 (as
  the saved '*' form, the reference); '2)*' -> coefficient 2; the 36 library GUI files give isequaln
  pde_struct objects with the HEAD copy and the working tree. Every failure is now retried with a
  leading '(', which can succeed only for an unmatched ')' (a text that parses but fails at run time
  becomes a parse error with the extra '(').

Whitespace note: the Edit tool dropped trailing blanks on five untouched neighbouring lines (three
blank lines in initialize_PIETOOLS_DDF.m group 2, `    elseif logVal==5   ` in the legacy initializer,
`    try ` in eval_b); restored byte-exact in the follow-up (cleanup 1, below).

## parser/legacy reachability (C15, C16, C17, C32) – measured
`B_conv/legacy_reach.m`, `legacy_reach_before.txt`: `which opvar` is empty (exist 0). state('pde'),
state('ode'), diff(x,t) and sys() run. Every operation that builds terms errors on the retired class:
int(x,s,[0 1]) at @state/int.m:89, int of a state independent of s at @state/mtimes.m:64, x+x at
plus.m:58, 2*x at mtimes.m:64, -x at uminus.m:43, diff(x,t)==x at minus.m:61, terms() ("'opvar' is not
recognized as a class": property operator of @terms/terms.m:38), addequation(sys(),x) at
state2terms.m:39. So no terms object (no equation) can be built. sys objects without equations reach
@sys/convert; the sys properties reject type 'nds' (mustBeMember) and a DDE/NDS/DDF struct as params
(mustBeA pde_struct/pie_struct), so the dde, nds and ddf branches of convert cannot receive their data.

## C15 – NEEDS DECISION
`parser/legacy/@state/int.m:97` is reached only after `opvar T;` (line 89), which errors (measured
above); the merge-conflict comments (56-71) are dead text. "parser/legacy depends on the retired opvar
class: port or retire the folder." Not edited.

## C16 – NEEDS DECISION
`parser/legacy/@sys/addequation.m:40-41`: the state branch errors in state2terms (opvar, line 39); the
terms branch cannot be reached because no terms object can be constructed (measured above). Same
decision as C15. Not edited.

## C17 – FIXED
- File: `parser/legacy/@sys/convert.m` (SS lines; the nds branch, line 52, is DJ c88f2787 and is not
  touched). Log 39-40; new line 60 (`params = obj.params;` in the 'pie' branch); lines 79-82 (`else` +
  error for convertTo 'ddf' from any other type).
- Reproduction (measured, HEAD copies of @sys and @state first on the path, `B_conv/c17_repro.m`,
  `c17_base.txt`): a 'pie' sys holding a PIE: `o = convert(P,'pie')` -> "Output argument "params" ...
  not assigned"; `o = convert(sys(),'ddf')` -> the same; with nargout 0 it prints "Conversion to ddf was
  successful" although nothing was converted.
- Verification (measured, `c17_after.txt`): the 'pie' sys returns its PIE (isequal(o,PIE) = 1, the PIE
  stored, the reference); `convert(sys(),'ddf')` -> "A sys object of type 'pde' cannot be converted to a
  DDF." with no success message; a 'pde' sys holding an initialized pde_struct still returns a
  pie_struct (unchanged).
- New findings (measured, not changed, outside the item): sys.m declares `obj = convert(obj,convertto)`
  (one output) while convert.m defines `[params,out]`, so `[p,o] = convert(...)` errors "Too many output
  arguments" and the documented second output cannot be obtained; `P.params` of a 'pie' sys without a
  PIE errors in getParams (no 'pie' branch).

## C32 – FIXED (message); declaration left
- File: `parser/legacy/@state/state.m` (SS). Log 41; (was) 124 / new 125: "State class definition takes
  at most 3 inputs (a 4th, the state name, is internal)".
- Reproduction (measured, `c17_base.txt`): `state('pde',1,[t,s],'nm',5)` -> "State class definition
  only takes 3 inputs", although the 4-input form (line 116) is accepted (measured: OK).
- Verification (measured, `c17_after.txt`): the new message for 5 inputs; 3 and 4 inputs construct.
- Not changed: the declaration `objterms = state2terms(obj,operator,var,val)` (line 56) against
  state2terms.m (obj only) – no runtime effect, and state2terms errors on opvar (C15/C16 decision).

## Files reported to have no caller (C21, C22, C23, C44, C46 part) – measured
`grep -rIl --exclude-dir=.git <name> .` over the whole PIETOOLS tree, and `unzip -p
PIETOOLS_PDE_GUI.mlapp matlab/document.xml | grep -c <name>` (the app code; parts listed with unzip -l:
document.xml, appModel.mat, metadata; no .mlx file exists). For isvalidTerm, identifyEqnType,
serializeStateVec, convert_PIETOOLS_PDE_batch2terms(_legacy) and convert_PIETOOLS_PDE_terms_internal the
only hits are the file itself and lists or reports under cleaning/ (release_comment_pass lists,
example_settings_1119/sol55/package_files.txt) and testing_codes/private/baselines/2026-11-02_opvar_phase0/
opvar_removal_plan.md; the mlapp code has 0 hits each. identifyEqnType and serializeStateVec are
@terms methods, and no terms object can be constructed at HEAD (measured above).

## C21 – NEEDS DECISION
`GUI/isvalidTerm.m` (2021b import, credited MMP, SS, DJ): no caller. "Dead file of another developer:
move to cleaning/private/dead_code per the dead-code rule." Not edited, not moved.

## C22 – NEEDS DECISION
`parser/legacy/@terms/identifyEqnType.m` (SS, 07/18/2022): no caller; unreachable (no terms object).
Same decision. Not edited, not moved.

## C23 – NEEDS DECISION
`parser/legacy/@terms/serializeStateVec.m` (SS, 02/14/2022): no caller; unreachable (no terms object).
Same decision. Not edited, not moved.

## C44 – NEEDS DECISION
`converters/legacy/convert_PIETOOLS_PDE_batch2terms_legacy.m` (SS 5/30/2021, 2021b import): no caller
(neither name). Same decision. Not edited, not moved.

## C46 – FIXED (text, live file); dead file NEEDS DECISION
- `converters/legacy/initialize_PIETOOLS_PDE_terms_legacy.m` (live: initialize_PIETOOLS_PDE): the note
  "The workspace is cleaned at the end of the file ..." (2021b import, e15f2543; header documentation)
  replaced by lines 75-76: the cleared variables are local to the function. Log lines 114-115. The
  col-77 markers first put on lines 75-76 were removed in the follow-up (cleanup 2): header text.
- `converters/legacy/convert_PIETOOLS_PDE_terms_internal.m` (SS, 12/05/2022): no caller (measured
  above) – treated as C21: NEEDS DECISION (dead file: move to cleaning/private/dead_code); its note and
  the developer question left as they are.

## Regression of the GUI route (C12, C36) – measured
The GUI converts with `converters/legacy/convert_PIETOOLS_PDE_terms.m` (mlapp document.xml line 1391),
which calls initialize_PIETOOLS_PDE_terms_legacy. `B_conv/gui_regress.m`: the 36 GUI files of the example
library parsed (parsePDEParams) and converted along that route, with HEAD copies of parsePDEParams and
initialize_PIETOOLS_PDE_terms_legacy first on the path (`gr_base.mat`) and with the working tree
(`gr_after.mat`): 36 of 36 give identical Galerkin matrices of all 12 operators (`gr_compare.txt`), 0
errors in either run.

# Follow-up (coordinator, 11/20/2026)

## C112 – FIXED
- File: `converters/internal/poly_decomp_1D_econ.m` (MMP 6/26/20; code lines from the 2021b import
  e15f2543 by `git blame -w -C -C -M`). Log 52-53; comment 55-56, (was) 57 / new 58
  (`if ~isa(Q,'polynomial') || isempty(Q.varname)`), (was) 61 / new 62 (`H = double(Q);`); (was) 81,
  comment 82-84, new 85-86: the monomials of column i built from its degmat, sorted by degree,
  `polynomial(speye(Nt),sortrows(Qi.degmat),Qi.varname,[Nt 1])` – what the one-argument monomials(p)
  of the multipoly toolbox did (PIETOOLS_2020a/multipoly/monomials.m, removed with 2021b).
- Probe (measured, `B_conv/c112_probe.m`): for Q = [s1 1 0; 2 s1^2+3 0; 0 4s1 0] the rebuilt bases
  are [1; s1], [1; s1; s1^2], [1] (zero column: one constant term), all monomial vectors (ismonom), and
  poly2basis reproduces each column exactly (e = 0). multipoly drops the variable of A+0*s1: such a
  kernel is a polynomial without variables, and `Q.Varname(1)` (now line 68, unused afterwards) failed
  "Index exceeds array bounds"; it now takes the constant branch, with H as 'double' (the svd of the
  callers needs a double).
- Reproduction (measured, `B_conv/c35_repro.m` with S2q added: S2c with every kernel a polynomial
  constant; `c35_c112before.txt`): S2p, DDE and NDSmin routes -> "Not enough input arguments" (in
  monomials); S3p, NDSmin -> the same; S2q, DDE and NDSmin -> "Index exceeds array bounds" (subsref).
- Verification (measured, `c35_c112after.txt`): S2p: NDSmin, NDSraw, DDE match the closed form to
  <= 1.8e-15 and each other to <= 3.1e-15, eigenvalues roots of det Delta (sigma ratio <= 1.4e-14);
  S2q: all three routes match the closed form to <= 1.2e-15; S3p: NDSmin matches to <= 5.6e-16; S1, S2c,
  S3c unchanged. Constant (double) kernels take the unchanged branch (H = double(Q) returns a double Q
  as is): the 46 library conversions are bit-identical with HEAD copies of all six edited converter
  files first on the path (`lr2_base.mat`, `lr2_after.mat`, `lr2_compare.txt`).

## C113 – FIXED
- File: `converters/initialize_PIETOOLS_DDE.m` (2021b import line). Log 54; comment 297, (was) 298,
  condition 299 ('B1di' for 'Bdi').
- Reproduction (measured, `B_conv/c113_repro.m`, `c113_before.txt`): a DDE in which B1di alone
  defines nw (constant kernel, and polynomial in s1) -> "The dimensions of B1di{1} are inconsistent with
  the detected number of inputs or states" (nw read as 0).
- Verification (measured, `c113_after.txt`): both pass with nw = 1 (B1 2x1, D11 1x1); the PIE transfer
  function w -> z matches the closed form C1 (lI - A0)^{-1} int_{-1}^0 B1di(s) e^{ls} ds to 2.5e-16
  (constant) and 4.0e-16 (polynomial) at three points; library DDEs bit-identical (`lr2_compare.txt`).

## Cleanup 1 – trailing whitespace restored byte-exact
Lines whose only difference from HEAD was trailing whitespace (found by pairing removed and added
lines of `git diff -U0 HEAD` modulo trailing blanks, over all eleven edited files): initialize_PIETOOLS_DDF.m
358, 366, 375 (`    `), initialize_PIETOOLS_PDE_terms_legacy.m 304 (`    elseif logVal==5   `),
parsePDEParams.m 702 (`    try `). Restored with `B_conv/wsfix/restore_ws.pl` (binary-safe perl; each
target line must read the stripped text, its CRLF or LF kept) on copies in `B_conv/wsfix/`; `diff` of
copy and result shows exactly these lines (sizes +12, +3, +1 bytes); copied back after `cmp` showed the
files unchanged since the copy. After: for each edited file `git diff HEAD` and
`git diff HEAD --ignore-space-at-eol` count the same changed lines (no whitespace-only change); CR/LF
counts 522/522, 886/886, 0/724 (LF file), as before.

## Cleanup 2 – markers removed from the header note
`converters/legacy/initialize_PIETOOLS_PDE_terms_legacy.m` lines 75-76 (the '~ NOTE' text of C46): the
col-77 markers removed; the log line (114-115) records the change.

## Files edited and final checks
converters/internal/convert_PIETOOLS_NDS2DDF.m, converters/convert_PIETOOLS_NDS.m,
converters/internal/minimize_PIETOOLS_DDE2DDF.m, converters/initialize_PIETOOLS_NDS.m,
converters/initialize_PIETOOLS_DDF.m, converters/legacy/initialize_PIETOOLS_PDE_terms_legacy.m,
GUI/parsePDEParams.m, parser/legacy/@sys/convert.m, parser/legacy/@state/state.m,
converters/internal/poly_decomp_1D_econ.m, converters/initialize_PIETOOLS_DDE.m.
`git ls-files --eol`: w/crlf for all, w/lf for GUI/parsePDEParams.m – as before the edits. Marker
columns checked (col 77, or one space after code reaching col 75). syntax_check over the eleven files
(`B_conv/edited_files.txt`): "11 files checked, 0 with parse problems" (`B_conv/syntax_check_out.txt`).
No commit, stash, reset or move.
