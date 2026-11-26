# Phase 0 review – code that names the 2025 classes

Maintainer, 11/26/2026: "each case needs independent review" (D2), and `poslpivar_sop` "needs
review" (D3). One section per case: what the code does, what reaches it, its behaviour today and
after a plain rename, the options, the proposal. (m) marks a measured result; the rest is read from
the code. State: HEAD f75d2f6c = origin/ndopvar, no local edits to tracked files.

- Probe: `inventory/probe_legacy.m` – MATLAB R2025b, one batch, under 2 minutes. The post-rename
  behaviour is simulated with the current name `copvar`, since the rename changes names only.
- Rescan: `inventory/h1scan.pl`, quote-aware, string literals kept; output `inventory/h1_rescan.txt`.
  73 lines: the 57 of `old_code.tsv` plus 16 message lines in `@copvar`/`@cdopvar`. The inventory
  parser took the `%s` inside those strings for a comment. Toolbox H1 is 49 lines, not 33.

## C1 – PIESIM `transl_mod.m:70–71` (YP; a 2025 file)

- Code: `if ~isa(T,'opvar') && ~isstruct(T)`, error "First argument must be a struct of kernels or
  an opvar object" (YP, MMP 11/21).
- Reached only by `PIESIM.m:226` → `rescalePIE` (17 calls) → `transl_mod`. Before that, PIESIM
  converts any `pie_struct` to structs of kernels (`pie2piesim`, `PIESIM.m:195, 205, 213`);
  `PIESIM_input_check.m:545` admits a `pie_struct` or a struct with field `dim`.
- Today: the `opvar` clause is dead; the 2025 class is off the path.
- After a plain rename: a container passes the check, and the body reads and writes `T.I`,
  `T.var1`, `T.Q1`, `T.R.R0`, … through the 1-D field emulation. Reachable only from a plain-struct
  PIE whose operators are containers; not measured.
- Proposal: test `~isstruct(T)` only; message "First argument must be a struct of kernels"; drop
  "'opvar' object, or" from the help of T. Full protocol: header entry (the 11/21 entry, "as well
  as an 'opvar'", no longer holds), `(was)` lines, markers.

## C2 – low-rank helpers, `PIETOOLS_demos/lowrank_2d_stability/private/` (CC; new files)

**C2a `pielr_evalop.m:34`** – `isa(X,'opvar') || isa(X,'opvar2d') || isa(X,'double') ||
isa(X,'polynomial')` → return.
- Containers return earlier, at lines 27–32 (CC 11/08), so line 34 never sees one, today or after
  the rename: no behaviour change. Lines 24–26 ("the opvar classes below left the path") become
  false after the rename.
- Proposal: drop the two 2025-class tests and that sentence.

**C2b `pielr_maxop.m:38–66`** – enumeration of the 2025 properties, with an assert on the count.
- The container branch (line 24) tests `isa(X,'copvar')` only. Today a `cdopvar` reaches line 38
  and stops with 'pielr_maxop:input'. After a plain rename a `dopvar` container matches line 38 and
  fails the assert of line 64, "enumeration changed" – a misleading error. `pielr_lpi.m:212` passes
  substituted operators; `pielr_opcheck.m:160, 296` pass residual-check operators, not traced.
- Proposal: delete lines 38–66 and the local `mxc`; another class raises 'pielr_maxop:input' naming
  the container class; revise help lines 5–17.

**C2c `pielr_norm_pie.m:23–41`** – repackages a bare struct whose T and A are 2025 classes.
- `pde_struct`, `pie_struct` and the sys wrapper return earlier. Today every bare struct stops with
  'pielr_solve:input'. Callers: `pielr_solve.m:162` (user input); `cuadmm/heat2d_compare.m:181` and
  `pielr_bench.m:107` pass a `pie_struct`.
- After a plain rename: a bare struct of 1-D containers passes with `S.dim = 1`, and `pielr_solve`
  runs on it – a path never exercised.
- Proposal: delete the branch; other inputs raise 'pielr_solve:input' naming the three accepted
  forms.

## C3 – `converters/legacy/convert_PIETOOLS_PDE_terms_internal.m:227–301`

SS 12/05/2022; in 2025 under `converters/internal/legacy/`; moved by Declan in a222bbc7 (11/16).
- No caller in the live tree; on the path, so callable directly. Never ported to the containers.
- Today: stops at line 227, `opvar Thatop Tvop;`, with "Undefined function 'opvar'" (m, the same
  construct in `@state`).
- After a plain rename: `opvar A B;` builds empty containers (command syntax, `copvar.m:141`);
  `.dim` assignments fail with 'copvar:setDim' (m); `.R.R0`, `.P`, `.I` succeed with a warning (m).
- Options: a. move to `orphaned_file/private`; b. report to Declan, Phase 1 waits; c. leave.
- Proposal: a, named in the announcement (D7).

## C4 – `parser/legacy/@state` (7 lines), `@terms` (4 lines): the legacy state/terms format

Moved by Declan; in 2025 under `parser/`.
- No live code constructs `state`, `terms` or `sys` objects (no call outside their class folders).
  `terms2pde_struct` raises "Objects of type 'terms' (a legacy input format) are no longer
  supported. Declare the PDE using 'pde_var' instead." (DJ 11/25).
- Today (m): `state('pde')` works; the first operation fails – `x+x` and `-x` give "Undefined
  function 'opvar' for input arguments of type 'char'".
- After a plain rename (m): `opvar T; T.R.R0 = eye(2)` succeeds with a warning, so `plus`,
  `uminus`, `uplus` and `state2terms` would run; `@terms/int` fails at `T.var2 = var`
  ('copvar:dummyName'); `@terms/vertcat` fails at `.dim` ('copvar:setDim'). A legacy script gets
  further than today, and through `@sys/convert` → `convert_PIETOOLS_PDE_terms_legacy` may even
  finish; not measured end to end.
- The rename script leaves these lines alone: they already read `opvar`.
- Options: a. the `state` constructor raises the `terms2pde_struct` message (one line; the 11 lines
  become unreachable); b. move `@state` and `@terms` to `orphaned_file/`; c. report to Declan,
  Phase 1 waits; d. leave.
- Proposal: a, reported to Declan in the announcement.

## C5 – messages of `@copvar`/`@cdopvar` that name the 2025 classes (new files)

- 34 message lines in 12 files (`getfield_opvar1d`, `getfield_opvar2d`, `setfield_opvar1d`,
  `setfield_opvar2d`, `subsasgn`, `subsref`, in each class). 22 name `opvar`/`dopvar`; 12 name only
  `opvar2d`/`dopvar2d`. (First stated as 36 and 14; recounted line by line when editing.)
- After a plain rename the 22 contradict themselves: "'%s' is a field of an 'opvar', but this opvar
  is not opvar-like"; "The operator does not have the structure of an 'opvar'"; "An opvar has no
  field '%s'; … as for an 'opvar' or 'opvar2d'" (today a warning, m).
- Proposal: name the layout, not the class – "the 1-D field layout (P, Q1, Q2, R)" for
  'opvar'/'dopvar', "the 2-D field layout (R00 … R22)" for 'opvar2d'/'dopvar2d', "of that layout"
  for "opvar-like" – in all 34, so each family of messages reads alike. Identifiers
  ('copvar:notOpvar1d', …) unchanged; text only. Alternative: the 22 only.

## C6 – test code (24 lines)

**C6a `sopvar_demos/test_heatNd_pie.m:253`** – `if isa(P,'opvar') || (isa(P,'copvar') &&
is_opvar1d(P))`.
- After a plain rename the first test admits every container: 2-D and 3-D cases take the 1-D branch
  and fail at `P.R.R0`.
- Proposal: `if isa(P,'copvar') && is_opvar1d(P)`.

**C6b `cuadmm/bl_suites.m:69`** – text: "a change to dpvar, opvar, poslpivar, lpi_eq or lpivar that
alters any assembled program changes one of these".
- No effect on any computation; after the rename the names denote the container successors, which
  build these programs.
- Proposal: leave.

**C6c phase-0 harness, `baselines/2026-11-02_opvar_phase0/private`.**
- `p0_f2.m:108, 123`: `legacy = isa(Xin,'opvar') || isa(Xin,'opvar2d')` sends 2025 operators
  through `op2copvar_sop` and a 2025-field reference evaluation. Dead today. After the rename,
  containers take that branch on the rename export and not on the parent export, so F2 differs for
  a harness reason (H10); `op2copvar_sop` is also retired by D4.
- `p0_classes.m:39` (a class listing of 11/02) and `p0_compare.m:53` (a report label) change no
  output.
- Proposal: delete the 2025 branches of `p0_fixop` in `p0_f2.m`; leave the other two.

**C6d class-ops baseline, `baselines/2026-11-05_class_ops/private` (18 lines).**
- `bench_class_ops.m` times the 2025 `opvar`/`dopvar` against `sopvar`/`sdopvar` (`D = dopvar();`,
  `P = opvar();` at 209, 225); `analyze_class_ops.m` and `plot_scaling.m` label the recorded
  results. No container is named. Today the 2025 part cannot run (classes retired 11/08); after the
  rename it would time containers under the 'opvar'/'dopvar' labels.
- Consequence for PLAN Phase 3: this baseline cannot be the timing spot check, since neither class
  it times is renamed.
- Proposal: leave the scripts as the 11/05 record; one README line states both facts; Phase 3 times
  container class operations on both exports instead.

## C7 – `poslpivar_sop` (D3)

- Not a second construction. It translates the 1-D degree format of the stock settings – D =
  {d1,d2,d3} with the 2025 fill rules, options psatz, exclude, sep – into the `deg`/`include` form of
  `poscopvar`, and calls `poscopvar` once (`poslpivar_sop.m:160`). It adds the SIDE choice with a
  squareness check and allows one variable per space.
- Callers: `poslpivar_settings_sop` (the 'lf' role and the stock slack, reading the stock fields
  dd1/options1 … dd3/options3; called by `exec_tools_sop`) and `PIETOOLS_H2_norm_sop.m:300, 304`
  (W, default degrees). Tests: `test_translators_sop` (11 calls).
- `poscopvar` cannot take this format directly: a 1 x 3 cell is its per-space degree form for three
  spaces.
- Options: a. keep in this pass, named in the `_sop` pass; b. fold it into `poslpivar_settings_sop`
  and route W of the H2 executive through `poslpivar`, retiring `poslpivar_sop` (identity check on
  the H2 and settings tests); c. merge into `poslpivar` – rejected, the degree forms collide.
- Proposal: a.

## Decisions (maintainer, 11/26/2026)

| Case | Decision |
|---|---|
| C1 `transl_mod` | Leave it if compatible with the container format, else modify it. Measured incompatible (below): modify. |
| C2a `pielr_evalop` | Drop the two 2025-class tests and the sentence, as proposed. |
| C2b `pielr_maxop` | Delete the 2025 branch and `mxc`, as proposed. |
| C2c `pielr_norm_pie` | Delete the bare-struct branch, as proposed. |
| C3 terms converter | Move to `orphaned_file/private/converters/legacy/`, as proposed. |
| C4 `@state`, `@terms` | Move both folders to `orphaned_file/private/parser/legacy/` – not the proposal; announced after the rename (D7). |
| C5 messages | Reword all of them (34; asked as 36) to name the layout, as proposed. |
| C6a `test_heatNd_pie` | Drop the `opvar` test, as proposed. |
| C6b `bl_suites` | Leave, as proposed. |
| C6c phase-0 harness | Delete the 2025 branches of `p0_fixop`; leave `p0_classes`, `p0_compare`; as proposed. |
| C6d class-ops baseline | README line; Phase 3 times container class operations instead; as proposed. |
| C7 `poslpivar_sop` | Rename to `poslpivar_legacy` (the name the maintainer gave; free in the repository and the 2025 release), in Phase 1 with the rename table. |

**C1, measured** (`inventory/probe_c1*.m`, `transl_mod_c*.m`; MATLAB batches under 2 minutes each).
- `transl_mod` with the check widened to the container (as `isa(T,'opvar')` behaves after the
  rename) fails on a 1-D container with every kernel nonzero, on [0,2]: 'copvar:badKernelVar',
  "'Q1' depends on 's_bar'". PIE.T of a heat equation on [0,2] fails at 'R.R1' the same way. PIE.A
  (constant R0 only) runs. Cause: each kernel passes through the field in the auxiliary variables
  `s_bar`, `theta_bar` before the back-substitution, and the container setter rejects that kernel.
- Change: translate each kernel in a local `K` and assign the field once (12 lines). The check line
  stays; after the rename it admits the containers. `T.I` stays unset, as since 8d717205 (YP
  05/02/2026, "Fixed rescaling").
- With the change: on the same three cases the container results equal the struct path to 0
  (every kernel), and the struct path is bit-identical to the original (`isequal`), two of the
  three being the `pie2piesim` structs that PIESIM uses.
- Implementation checks for C4: what in `parser/legacy/@sys` and elsewhere still calls `@state` or
  `@terms` methods.

## Implementation, 11/26/2026 (maintainer: "Yes, implement")

Committed 11/26 on f75d2f6c (maintainer: "commit them"): a67b402f (C1, C2, C5, C6), c341990e (C3,
C4, D4: the moves, their call sites, the two READMEs). Not pushed. Annotation: the full protocol in files of the 2025 release (`transl_mod`, DEMO7) and in
non-release files that carry developer notes (low-rank helpers, test harnesses); none in the new
release files (`@copvar`/`@cdopvar` messages, DEMO12).

| Case | Files | Change |
|---|---|---|
| C1 | `PIESIM/.../transl_mod.m` | each kernel through a local `K`, assigned once (5 pairs); `(was)` lines, markers, log entry |
| C2a | `pielr_evalop.m` | the `opvar`/`opvar2d` pass-through tests and their sentence deleted; help note |
| C2b | `pielr_maxop.m` | the 2025 enumeration (lines 38–66) and `mxc` deleted, replaced by the input error; help revised |
| C2c | `pielr_norm_pie.m` | the bare-struct branch deleted, replaced by an error naming the accepted inputs; help revised |
| C3 | `converters/legacy/convert_PIETOOLS_PDE_terms_internal.m` | moved to `orphaned_file/private/converters/legacy/` |
| C4 | `parser/legacy/@state` (20 files), `@terms` (14) | moved to `orphaned_file/private/parser/legacy/` |
| C5 | 12 files of `@copvar`/`@cdopvar` | 34 messages name the layout |
| C6a | `test_heatNd_pie.m` | `if isa(P,'copvar') && is_opvar1d(P)`; log entry, `(was)` line |
| C6c | `p0_f2.m` | the 2025 branches of `p0_fixop` and `p0_apply_opvar` deleted; `legacy` always false |
| C6d | class-ops `README.md` | dated note |
| D4 | `op2copvar_sop.m`, `pie2copvar.m` | moved to `orphaned_file/private/converters/operator_converters/` |
| D4 | DEMO7, DEMO12 | the calls deleted (DEMO12 reads `PIE.T`, `PIE.A`, `PIE.Bu`; help revised) |
| D4 | `hd_qform.m` (3 baselines) | `Tm = PIE.T;` … in place of `op2copvar_sop(PIE.T)` … |
| D4 | `lb_setup.m` (4 baselines), `p0_path.m` | `pie2copvar` out of the path assertions; `op2copvar_sop` into the retired names |
| – | `orphaned_file/private/README.md` | section "Moved 11/26/2026" |
| – | `testing_codes/README.md` | `pie2copvar` (retired now) and `opvar2d2copvar` (retired earlier, still listed) out of the list of functions only tests call |

Line endings: the working tree mixes CRLF and LF; `core.autocrlf` is true and `.gitattributes` sets
`text=auto`, so the repository receives LF and `git diff` shows content changes only. Leftover
references: `op2copvar_sop` in the retired-name list of `p0_path.m` (intended); `state2terms` in a
branch of `@sys/addequation.m` that `isa(eqn,'state')` no longer enters.

### Tests, 11/26/2026 (MATLAB R2025b, MOSEK 11; 03:09–03:19)

Edited tree against an export of f75d2f6c; harness `p0t/` in the session scratchpad (`p0t_setup`
asserts one PIETOOLS tree on the path).
- Parse check (mtree) of the 28 edited .m files: 0 errors.
- Moved names: `op2copvar_sop`, `pie2copvar`, `convert_PIETOOLS_PDE_terms_internal`, `state`,
  `state2terms` resolve nowhere; `terms` resolves only to `Terms.m` methods of the Model-Based
  Calibration toolbox.
- Messages (C5): all eight families and the warning triggered on 1-D, 2-D and 3-D containers; each
  prints the layout wording, none names a 2025 class.
- `transl_mod` (C1): struct path edited = HEAD (`isequal`, PIE.T and PIE.A through `pie2piesim`);
  containers through the widened check = struct path to 0; a `copvar` is still refused before the
  rename.
- `pielr_*` units (C2): double and copvar pass through `pielr_evalop`; `pielr_maxop` on a copvar
  gives m = 1, npar = 3, on a cdopvar the input error; `pielr_norm_pie` initializes a pie_struct and
  refuses a bare struct with the new message.
- DEMO1, DEMO7, DEMO12, edited against HEAD, every saved variable compared: DEMO1 identical; DEMO7
  identical apart from solver CPU times and figure handles; DEMO12 identical apart from the CPU time,
  with `PIEc` absent by design.
- `test_heatNd_pie`: all 70 checks pass, on the edited tree and on HEAD. Its `heatNd_path` asserts
  that `cx_resid` resolves, which needs `.../claude_tests/cx_exec` on the path (pre-existing).
- `test_pde_var_components`: completes.
- Low-rank `pielr_bench` case `rd1d-lam0.5` (tier 1): not certified, rel NaN, while the IPM
  reference reaches 1.188e-8 – the same on HEAD, so not caused by these edits.
- Old-class code after Phase 0 (`inventory/h1scan.pl`): 23 lines, all decided to stay –
  `transl_mod.m:70–71` (C1), `bl_suites.m:69` (C6b), `p0_classes.m`, `p0_compare.m` (C6c), the
  class-ops record (C6d, 18).
