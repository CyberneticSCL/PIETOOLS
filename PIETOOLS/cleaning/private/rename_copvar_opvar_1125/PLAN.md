# Rename `copvar` → `opvar`, `cdopvar` → `dopvar` – plan

Maintainer, 11/25/2026: "create a plan for renaming the copdvar and cdopvar objects to opvar and
dopvar. Leave sopvar/sdopvar as is. this should be done throughout, including comments. also compile
a list to functions with the cop/copvar suffix with a goal to remove these suffixes to try and match
the 2025 release naming conventions where possible"; then "Also applies to subroutines".

Status: plan; decisions D1–D8 taken on 11/26 (Section 4). Phase 0 steps 1 (C1–C6), 2 and 3
implemented, tested and committed on 11/26 (a67b402f, c341990e, records 8198c35d; `REVIEW_PHASE0.md`).
The H2 comment pass of step 1 (done 11/26) committed 11/27 (1d9737e1); step 4 done on 11/26; step 5:
origin unchanged at f75d2f6c. Case C8 (`tests_1d`) moved to `orphaned_file/private` on 11/26,
committed 11/27 (3315e9f9). Phase 1 dry run on an export, 11/26, and applied to the checkout,
11/27 (`phase1/README.md`): every check of step 3 passes on both; committed 11/27 (ae2b95c1).
Phase 1b committed 11/27 (abfb0e83, `phase1b/README.md`). None pushed. Every count below was measured at f75d2f6c (local = origin/ndopvar). Scripts and raw outputs: `inventory/`. Function and
subroutine table: `rename_table.tsv` (Section 6).

## 1. Scope

In scope:
- the classes `copvar` → `opvar` and `cdopvar` → `dopvar`: class folders, constructors,
  `InferiorClasses` lists, `isa`/`class` strings, error identifiers, help text, comments;
- every function and subroutine (local function) whose name carries `copvar`, `cdopvar` or `cop`:
  the 2025 release name where one exists, else the same name with the class token replaced
  (Section 6);
- every live file: the toolbox (all of `PIETOOLS/` except `testing_codes/`, `cleaning/`,
  `orphaned_file/`), test code under `testing_codes/`, `testing_codes/README.md`,
  `.claude/agents/pietools-applications.md`, and the class table of the untracked `CLAUDE.md`.

Unchanged:
- `sopvar`, `sdopvar` and every other class (`nopvar`, `ndopvar`, `polyopvar`, `tensopvar`,
  `tensopmat`, `intop`, `pde_struct`, `pie_struct`);
- SOSTOOLS400 – no file names any of these classes (measured);
- records, as in the 09/25 rename `mopvar` → `copvar` (ffc3a939): `cleaning/` (496 files, 4,972
  lines name a container), `orphaned_file/` (17 files, 214 lines), the dated records under
  `testing_codes/` (30 csv, md, log, sh and txt files, 737 lines), and every header change-log
  entry. A translation note goes into `cleaning/README.md` and `testing_codes/README.md`.

## 2. Inventory at f75d2f6c

Lines that name `copvar` or `cdopvar` (case-insensitive). A line counts as code when the name
appears before the first `%` outside a quoted string.

| Area | Files | Lines | Code lines | Comment or text lines |
|---|---:|---:|---:|---:|
| sopvar | 94 | 732 | 306 | 426 |
| executives_sopvar | 23 | 289 | 111 | 178 |
| lpi_programming_sopvar | 17 | 176 | 67 | 109 |
| converters | 17 | 124 | 61 | 63 |
| PIETOOLS_demos | 24 | 121 | 55 | 66 |
| parser | 16 | 117 | 56 | 61 |
| PIESOS | 2 | 24 | 10 | 14 |
| PIETOOLS_examples | 3 | 13 | 2 | 11 |
| PIETOOLS_DDE.m, PIETOOLS_DDF.m, PIETOOLS_PDE.m | 3 | 3 | 0 | 3 |
| .claude/agents | 1 | 2 | 0 | 2 |
| **Toolbox** | **200** | **1,601** | **668** | **933** |
| testing_codes, .m files | 83 | 698 | 356 | 342 |

Further measurements:
- 1,882 occurrences in the toolbox. 152 occurrences of "a copvar" or "a 'copvar'" become "an
  opvar" (117 toolbox, 24 test code, 11 records).
- 98 live files carry a container name in their path: 63 in `sopvar/@copvar` and
  `sopvar/@cdopvar`, 29 toolbox function files, 6 test files.
- 67 distinct error identifiers `'copvar:…'` and `'cdopvar:…'`.
- `InferiorClasses` lists naming a container: `parser/@pde_struct/pde_struct.m:1` (`?copvar`,
  DJ 11/25) and `sopvar/@cdopvar/cdopvar.m:1` (`?copvar`).
- 35 of the 200 toolbox files existed in the 2025 release (key: `@class/method` inside a class
  folder, else basename, against origin/PIETOOLS-2025): parser 13, demos and manual snippets 11,
  converters 7, root scripts 3, `closedLoopPIE` 1. The other 165 are new in this release.
- Stored data: none of the 42 tracked .mat files holds a container or a 2025 `opvar` (byte scan
  of the inflated MAT elements, `inventory/matscan.pl`). `PIETOOLS_PDE_GUI.mlapp` names neither
  class.
- Name resolution, MATLAB R2025b, on the saved path and again after `addpath(genpath(PIETOOLS))`:
  `opvar`, `dopvar` and all 39 proposed names resolve nowhere (`inventory/namecheck.m`). Retired
  copies of 8 of them (`opvar`, `dopvar`, `poslpivar`, `lpivar`, `mat2opvar`, `rand_opvar`,
  `pie2opvar`, `apply_opvar`) lie in private folders of `cleaning/` or `orphaned_file/`, off the
  path.
- Branches: only origin/ndopvar moves (5 Declan commits on 11/25); origin/master last moved 08/19.
- Rerun 11/26 on the tree after Phase 0 (8198c35d plus the H2 pass; `inventory/areas_1126.tsv`,
  `callers_1126.txt`): toolbox 202 files, 1,594 lines, 661 code lines. converters 15 files
  (`op2copvar_sop`, `pie2copvar` moved); PIESIM 1 file and PIETOOLS_demos 3 more files from header-log
  lines of 11/26, which Phase 1 leaves as written. Rerun again on the base of Phase 1.

## 3. Hazards

A plain substring swap breaks in ten places.

**H1 – code that names the retired classes.** After the rename, `isa(x,'opvar')`,
`isa(x,'dopvar')` and `opvar X;` refer to the containers. 33 code lines in the toolbox:

| File(s) | Lines | Code | Today | After a plain rename |
|---|---:|---|---|---|
| `PIESIM/PIESIM_routines/PIESIM_Setup/transl_mod.m:70–71` | 2 | `if ~isa(T,'opvar') && ~isstruct(T)` | a container raises the error | a container passes; the body reads `T.I`, `T.var1`, `T.Q1` through the 1-D field emulation – unmeasured |
| `PIETOOLS_demos/lowrank_2d_stability/private/pielr_evalop.m`, `pielr_maxop.m`, `pielr_norm_pie.m` | 6 | `isa(X,'opvar')` dispatch | dead branches | taken by containers |
| `converters/legacy/convert_PIETOOLS_PDE_terms_internal.m:227–301` | 8 | `opvar Aop_PDE Bop_PDE;`, then `.dim`, `.R.R0` | "Unrecognized function 'opvar'" if reached | builds containers, then assigns 2025 fields |
| `parser/legacy/@state` (7 files), `parser/legacy/@terms` (4 files) | 11 | `opvar T; T.R.R0 = …` | as above | as above |
| `sopvar/@copvar`, `@cdopvar`: `getfield_opvar1d.m`, `setfield_opvar1d.m` | 6 | error text "not the structure of an 'opvar'" | correct | self-contradictory |

Test code adds 24 lines: class-ops baseline scripts 18, phase-0 harness 4, `bl_suites` 1,
`test_heatNd_pie` 1 (`inventory/old_code.tsv`). A quote-aware rescan on 11/26
(`inventory/h1scan.pl`) found 16 more message lines in `@copvar`/`@cdopvar`: the toolbox holds 49
such lines, not 33. Cases and decisions: `REVIEW_PHASE0.md`.

**H2 – comments that name the retired classes become ambiguous.** In toolbox .m files, 201 lines
in 48 files name `opvar` or `dopvar`: 44 header-log lines (kept as written), 108 comment lines, 49
lines of code or with a trailing comment (`inventory/old_release_lines.tsv`). 22 toolbox lines name
both `copvar` and `opvar`, 13 both `cdopvar` and `dopvar`, 3 both `poscopvar` and `poslpivar`;
after a swap they read "converts an opvar to an opvar". Upper-case names of retired converters
(`OPVAR2COPVAR`, `COPVAR2OPVAR`, `COPVAR2OPVAR2D`, See-also entries kept by the maintainer decision
of 11/20) sit in `pie2copvar.m`, `@copvar/inv.m` and two test files.

**H3 – four new names are 2025 names with other signatures.**

| New name | 2025 signature | Current signature |
|---|---|---|
| `poslpivar` | `[prog,Pop,Qmat,Zop,gs] = poslpivar(prog,n,d,options)` | `[prog,Pop,Qcell,basis_list] = poscopvar(prog,dims,spaces,dom,deg,options)`, or `poscopvar(prog,X,deg,options)` |
| `lpivar` | `[prog,Zop] = lpivar(prog,n,d,opts)` | `[prog,Pop] = lpivar_cdopvar(prog,dims,spaces,dom,deg,options)`, or `(prog,X,deg,options)` |
| `mat2opvar` | `Pop = mat2opvar(Pmat,Pdim,vars,dom)` | `Pop = mat2copvar_sop(Mat,dims,spaces,dom,options)` |
| `rand_opvar` (test utility) | `L = rand_opvar(nmat,pvar_deg,var1,var2,I)` | `Pop = rand_copvar(spaces,dims,dom,deg,dnsty,occ)` |

A 2025 script that calls these, or that writes `opvar T; T.P = …; T.R.R0 = …` (the containers
emulate the 2025 1-D fields through `getfield_opvar1d` and `setfield_opvar1d`), reaches the new
code. Whether each such call stops with a clear error or runs with another meaning is unmeasured;
D6 makes the three functions stop with a migration error. `rand_opvar` lies in
`testing_codes/private`, off the path, so 2025 scripts cannot reach it.
A PIETOOLS 2025 installation on the MATLAB path shadows the new classes, or is shadowed by them,
by path order and without an error.

**H4 – class loading.** MATLAB refuses to load a class whose `InferiorClasses` names a missing
class. The folder renames, both constructors and the two lists of Section 2 land in one commit.

**H5 – error identifiers.** The 67 identifiers change with the class; tests that match them change
in the same commit.

**H6 – saved objects.** An object whose class no longer exists loads as a 0x0 double without a
warning (measured 11/04). No tracked .mat holds a container. Untracked result files under
`cleaning/` and ignored baseline fixtures that hold one load empty after the rename and stay
readable from a pre-rename checkout only.

**H7 – mechanics.** A folder rename fails with "Permission denied" while a MATLAB session watches
the folder (11/08), and the checkout lies under OneDrive. Close every MATLAB session and pause
OneDrive sync for the rename; `git mv` file by file.

**H8 – concurrent edits.** The parser (16 touched files, 13 of them in the 2025 release) and
PIESOS (2 files) are Declan code, and origin/ndopvar receives his commits daily. Git follows a
committed rename in a merge; an uncommitted edit to a renamed file has to be carried over by hand.

**H9 – documentation.** `PIETOOLS_2025_User_Manual.pdf` (two copies) names the 2025
`opvar`/`dopvar` 316 times; after the rename those pages describe another class under the same
name. The manual chapters drafted for the containers follow the rename.
`sopvar/sopvar_lift_notes.tex` (20 lines) and `sopvar/sopvar_duality_pruning_notes.tex` (5 lines)
change, and their PDFs are rebuilt (MiKTeX pdflatex is installed).

**H10 – baseline harnesses.** Phase 3 compares runs before and after the rename on exports of both
commits. Harness code that names the classes (`p0_classes`, `p0_compare`, `p0_f2`, the class-ops
scripts, `post_verify`) has to accept both names during that comparison.

## 4. Decisions (maintainer, 11/26/2026)

**D1 – annotation: the release rule.** The 165 new files get no log line and no markers; the 35
existing files get one log line, `MMP, <date>: copvar, cdopvar renamed opvar, dopvar, with the
functions of the rename table; names only.` No per-line markers and no `(was)` lines. Phase 0 and
Phase 1b edits change behaviour: in the 2025 files they follow the full protocol, in new files the
release rule.
- The 35 files used the 2025 names in 2025 but were rewritten for the containers since: 230 of
  their lines name `copvar` or `cdopvar`, and after the map only 26 of them equal a line of the
  2025 file (18 of 108 code lines, 8 of 122 comment lines; per-line markers removed before the
  comparison, `inventory/d1_norm.pl`). 16 of the 35 carry earlier log entries that name a
  container; the new line keeps those entries consistent with the code (CLAUDE.md §1).
- The classification is taken before the rename. Afterwards some class-folder files and the files
  of `poslpivar`, `lpivar` and `mat2opvar` carry keys of 2025 files (`@opvar/mtimes`,
  `poslpivar`); their content is new in this release, and they stay new.

**D2 – the H1 code: each case reviewed on its own.** No batch fix. Phase 0 reviews each file of the
H1 table, and the 24 test lines, separately: what the line does, what reaches it, the proposed
change. Each change is approved on its own, the Declan files included. The 11/25 proposals
(`isstruct(T)` in `transl_mod`, deleting the `pielr_*` branches, a report on the legacy files,
rewording the six `getfield_opvar1d`/`setfield_opvar1d` messages) are starting points, not
decisions. Reviewed 11/26: twelve cases, C1–C7, each decided by the maintainer
(`REVIEW_PHASE0.md`).

**D3 – `poscopvar` → `poslpivar`** (first answered `posopvar`, then reverted); `lpivar` and
`mat2opvar` kept. Review item: `poslpivar_sop` (callers `PIETOOLS_H2_norm_sop`,
`poslpivar_settings_sop`), thought likely an oversight. Reviewed 11/26 (C7 of `REVIEW_PHASE0.md`):
it translates the 1-D degree format of the stock settings and calls `poscopvar` once
(`poslpivar_sop.m:160`). The maintainer chose a descriptive name, `poslpivar_legacy`, given in
Phase 1 with the rename table; the later `_sop` pass then meets no clash.

**D4 – retire `op2copvar_sop` and `pie2copvar`.** Phase 0 moves both to `orphaned_file/private` and
deletes their calls: DEMO7:140, DEMO12:81 (Dan Braghini), 5 and 4 test files.
- `op2copvar_sop` has returned its input unchanged since 11/22.
- `pie2copvar` returns a container PIE unchanged and otherwise calls `initialize`, which accepts
  only containers since b1211608.
- Renamed, they would read `op2opvar` and `pie2opvar`; `pie2opvar` named the opposite conversion
  (containers to the 2025 classes) until 6757b6d3, 11/08.

**D5 – the lift chain is renamed with the rest:** `poscopvar_lift` → `poslpivar_lift`,
`posmult_cdopvar` → `posmult_dopvar`, `copquadvar` → `opquadvar`. No release code calls them
(`poscopvar_lift` has no caller, `posmult_cdopvar` is called only by it, `copquadvar` appears in
its error text and in 2 test files); retiring them is a separate question.

**D6 – a 2025 call stops with a clear error.** `poslpivar`, `lpivar` and `mat2opvar` recognize the
2025 argument forms of H3 and stop with a message that gives the new call – Phase 1b, a separate
commit, since Phase 1 must pass code identity. `rand_opvar` needs nothing: it lies in
`testing_codes/private`, off the path. The 2025 field assignments (`opvar T; T.R.R0 = …`), which
the containers emulate on purpose, are outside D6; Phase 3 records what they do.

**D7 – no freeze.** The Declan files (16 parser files with `mtimes_copvar` and `pde_struct.m:1`, the
two PIESOS converters) change in the rename commit, and the announcement follows the push. An
uncommitted edit of his to a renamed file is carried over by hand when he pulls (H8). Phase 1
still waits for Phase 0.

**D8 – locals kept.** `is_cop` (4 uses) and `cop` (3 uses) in `parser/piess.m` keep their names; the
class string inside, `isa(P,'copvar')`, changes. Correction of 11/25: `ncop` (`piesos_quadopvar.m`,
`Test_piesos_tensor_gram.m`) counts variable copies, not containers.

## 5. Phases

**Phase 0 – preconditions (separate commits).**
1. The changes approved in `REVIEW_PHASE0.md` (C1–C6, reviewed and decided 11/26). Then the H2
   comment lines: name a retired class "the 2025 `opvar` class" where the history matters, else
   drop the mention. Header logs stay as written. Done 11/26, committed 11/27 (1d9737e1): 143 edits in 34 files,
   comment-only (`REVIEW_PHASE0.md`, H2).
2. The `poslpivar_sop` review of D3: done 11/26 (C7); its rename joins Phase 1.
3. The D4 retirements.
4. Check that no toolbox line names `opvar` or `dopvar` once comments are stripped, then rerun the
   inventory scripts. Done 11/26: `transl_mod.m:75–76` only (C1); inventory in Section 2.
5. No freeze (D7): merge origin, close MATLAB sessions, pause OneDrive. 11/26: origin/ndopvar
   unchanged at f75d2f6c, nothing to merge; MATLAB and OneDrive at the start of Phase 1.

**Phase 1 – the rename (one commit, because of H4 and H5).**
1. `git mv` the 98 files: `sopvar/@copvar/` → `sopvar/@opvar/` with `copvar.m` → `opvar.m`, the
   same for `cdopvar`, and the function files per `rename_table.tsv`. At the Phase 1 base: 97
   (D4 retired two, C7 added `poslpivar_sop`).
2. One script edits the live text files, in this order:
   - a. function and subroutine names, longest first (`poscopvar_direct` before `poscopvar`), in
     lower and upper case (`POSCOPVAR_DIRECT`);
   - b. remaining identifiers that contain `copvar` or `cdopvar` (test locals such as
     `cg_copvar_pd`), by substring;
   - c. the class names in lower, upper and capitalized form (`copvar`, `COPVAR`, `Copvar`);
   - d. the article: "a copvar" → "an opvar";
   - e. the D1 log lines, on the 35 files classified existing before the rename.

   Left as written, as records (H2 pass, 11/26): header-log entries (Section 1); `(was)` lines,
   which quote the code of their date (27 live lines name a container; the map would turn
   `test_heatNd_pie.m:258` into `isa(P,'opvar') || (isa(P,'opvar') && …)`); and the retired names
   that hold a container token, in any case: `op2copvar_sop`, `pie2copvar`, `opvar2copvar`,
   `copvar2opvar`, `opvar2d2copvar`, `cg_copvar_pd`, `Test_cg_copvar_nd`, `test_copvar_kyp`,
   `test_poscopvar_vs_poslpivar_2d` (`inventory/comment_only_tokens_1126.tsv`). Checks 3a and 3c
   allow these.
3. Checks before the commit:
   - a. no `copvar` or `cdopvar` remains outside the records of Section 1;
   - b. code identity: in each .m file, the code without comments equals the map applied to the
     old code without comments (the `codeident` tool of the release comment pass). Any other
     difference is a defect;
   - c. comment identity: the comments equal the map plus the article rule applied to the old
     comments; the only added lines are the D1 log lines;
   - d. `which -all`: `opvar`, `dopvar` and each new name resolve exactly once; each old name
     resolves nowhere;
   - e. `checkcode` and the parse check pass on every edited file.

**Phase 1b – 2025 call forms (D6, separate commit).** In `poslpivar`, `lpivar` and `mat2opvar`,
detect the 2025 argument forms of H3 and raise an error that gives the new call. The detection
tests the number and types of the arguments once per call. The files are new in this release, so
no log line and no markers (D1); a comment states what is detected and why. Done 11/27 (abfb0e83):
31 calls in the 2025 forms stop, 23 current-form outputs unchanged (`phase1b/README.md`).

**Phase 2 – documentation.**
- `testing_codes/README.md`, `.claude/agents/pietools-applications.md`, and the class table in
  `CLAUDE.md`.
- The two sopvar notes: tex sources and rebuilt PDFs.
- The manual drafts.
- The translation notes in `cleaning/README.md` and `testing_codes/README.md`.
- A memory entry with the map.

**Phase 3 – verification.** A rename changes no arithmetic, so every number must match.
- The chain (CLAUDE.md §3), 1-D before 2-D and N-D:
  1. SOSTOOLS400, `sopvar` and `sdopvar` are untouched;
  2. `opvar` and `dopvar`: `Test_opvar`, `Test_dopvar`, the converter and indexing tests;
  3. the `lpi_programming_sopvar` tests;
  4. `test_executives_sop`;
  5. demos and examples;
  6. PIESIM.
- Identity before and after, on exports of the parent commit and the rename commit:
  - runs: the 24-row runner (7 min on 11/24); the container LPI and executive tests (1.6 min on
    11/07); F1–F4 of the phase-0 baseline (22 min on 11/08, so after 04:30);
  - what must match: At, b, c, K, x and y bit for bit. Solver times are excepted.
- Cost: no computation changes, so time and memory in q and in the number of spatial variables stay
  the same by construction. Identical SDP data plus one timing spot check of container class
  operations on both exports confirm it. The class-ops baseline of 11/05 cannot serve: it times
  the 2025 classes and `sopvar`/`sdopvar`, neither renamed (C6d).
- D6 check: each 2025 call form of `poslpivar`, `lpivar` and `mat2opvar`, from their 2025 headers,
  raises the migration error; the current forms run unchanged (covered by the identity runs).
  Record what the 2025 1-D field assignments of the manual (Chapter 3) do.

**Phase 4 – integration.** Commit and push, each on approval. Then the announcement (D7), drafted
for the maintainer to send to Declan and the other developers: the rename table, and the
instruction to commit or push local work before pulling.

## 6. Functions and subroutines with `cop` in the name

Full table with locations: `rename_table.tsv`. Callers are files whose code (comments stripped)
names the function, its own file excluded, as toolbox / tests. New names follow the 2025 patterns:
`poslpivar`, `lpivar`, `lpi_<op>` for LPI declarations; `<op>_opvar` (2025: `inv_opvar`,
`diff_opvar`, `rand_opvar`) and `X2opvar` (2025: `mat2opvar`, `dopvar2opvar`) for operator
routines.

New name equal to a 2025 name:

| Current | Kind | Callers | New | 2025 file |
|---|---|---:|---|---|
| `copvar` | class | 86 / 36 | `opvar` | `opvar/@opvar` |
| `cdopvar` | class | 47 / 11 | `dopvar` | `opvar/@dopvar` |
| `poscopvar` | function | 17 / 10 | `poslpivar` | `lpi_programming/poslpivar.m` |
| `lpivar_cdopvar` | function | 19 / 17 | `lpivar` | `lpi_programming/lpivar.m` |
| `mat2copvar_sop` | function | 10 / 1 | `mat2opvar` | `opvar/mat2opvar.m` |
| `rand_copvar` | test utility | 0 / 6 | `rand_opvar` | `opvar/rand_opvar.m` |

No 2025 counterpart:

| Current | Kind | Callers | New |
|---|---|---:|---|
| `eye_copvar_sop` | constructor | 8 / 6 | `eye_opvar` |
| `zeros_copvar_sop` | constructor | 3 / 1 | `zeros_opvar` |
| `pi4_copvar_sop` | constructor | 4 / 0 | `pi4_opvar` |
| `int2copvar` | constructor | 0 / 1 | `int2opvar` |
| `qp_copvar` | utility (Talitsky) | 4 / 1 | `qp_opvar` |
| `op2copvar_sop` | converter | 1 / 5 | retire (D4) |
| `pie2copvar` | converter | 1 / 4 | retire (D4) |
| `lpi_eq_cdopvar` | internal | 1 / 9 | `lpi_eq_dopvar` |
| `poscopvar_direct` | internal | 2 / 6 | `poslpivar_direct` |
| `poscopvar_lift` | internal | 0 / 0 | `poslpivar_lift` (D5) |
| `posmult_cdopvar` | internal | 1 / 0 | `posmult_dopvar` (D5) |
| `copquadvar` | internal | 1 / 2 | `opquadvar` (D5) |
| `lift_copvar` | internal | 3 / 0 | `lift_opvar` |
| `cat_copvar_grid` | internal | 6 / 0 | `cat_opvar_grid` |
| `derive_copvar_meta` | internal | 2 / 0 | `derive_opvar_meta` |
| `dpvar_op_copvar` | internal | 16 / 0 | `dpvar_op_opvar` |
| `eq_copvar` | internal | 2 / 0 | `eq_opvar` |
| `mat2copvar_grid` | internal | 6 / 1 | `mat2opvar_grid` |
| `merge_copvar_registry` | internal | 6 / 0 | `merge_opvar_registry` |
| `split_copvar` | internal | 1 / 0 | `split_opvar` |
| `verify_copvar_meta` | internal | 2 / 0 | `verify_opvar_meta` |
| `copvar_space_list` | internal | 19 / 10 | `opvar_space_list` |
| `parse_copvar_spaces` | internal | 10 / 0 | `parse_opvar_spaces` |
| `copvar2tensopmat` | PIESOS (Declan) | 1 / 0 | `opvar2tensopmat` |
| `tensopmat2copvar` | PIESOS (Declan) | 1 / 0 | `tensopmat2opvar` |
| `mtimes_copvar` | parser private (Declan) | 1 / 0 | `mtimes_opvar` |
| `poslpivar_sop` | executive internal | 2 / 1 | `poslpivar_legacy` (C7) |

Subroutines (local functions):

| Current | Defined in | Uses | New |
|---|---|---:|---|
| `is_id_cdopvar` | `sopvar/@cdopvar/get_display.m:381` | 2 | `is_id_dopvar` |
| `check_copvar` | `Test_int_sopvar.m:164` | 10 | `check_opvar` |
| `maxabs_copvar` | `baselines/2026-11-18_studies/private/post_verify.m:160` | 6 | `maxabs_opvar` |
| `apply_cop` | three converter tests (`pie_test_utils.m:62`, `test_mixed_domain_pde.m:954`, `test_parser_edge_cases.m:367`) | 29 | `apply_opvar` – a local `apply_op` would shadow the `@intop` method of that name |

Tests: `Test_copvar` → `Test_opvar`, `Test_cdopvar` → `Test_dopvar`, `rand_cdopvar` →
`rand_dopvar`, `Test_qp_copvar_nd` → `Test_qp_opvar_nd`, `test_copquadvar_faces` →
`test_opquadvar_faces`.

No method or private function inside `@copvar` or `@cdopvar` carries `cop`. The methods named after
the 2025 layout (`getfield_opvar1d`, `setfield_opvar1d`, `is_opvar1d`, `is_opvar2d`,
`getfield_opvar2d`, `setfield_opvar2d`) keep their names.

## 7. Effort and run time

- Phase 1 edits about 2,300 live lines (1,601 toolbox, 698 test code) plus 35 log lines. The 09/25
  rename edited 518 lines plus 214 log lines.
- Phase 3 runs about 9 minutes inside the 20-minute window (runner and tests); F1–F4 run after
  04:30.
- Rollback: the rename is one commit, and `git revert` undoes it until other commits build on it.
  The freeze keeps that window short.

## 8. Out of scope, listed for a later pass

- The `_sop` suffix (41 toolbox functions). The 2025 names of the user-level ones are free on the
  path: `lpi_eq`, `lpi_ineq`, `getsol_lpivar`, `lpigetsol`, `lpiscript`, `get_lpivar_degs`,
  `getController`, `getObserver`, `PIETOOLS_stability`, `PIETOOLS_Hinf_gain`, and others. The
  Phase 1 machinery serves a second pass. `poslpivar_sop` becomes `poslpivar_legacy` in this pass
  (C7), so that pass meets no clash with `poslpivar`.
- Folder names `lpi_programming_sopvar` and `executives_sopvar` (2025: `lpi_programming`,
  `executives`). The class folders stay in `sopvar/` (2025: `opvar/@opvar`).
