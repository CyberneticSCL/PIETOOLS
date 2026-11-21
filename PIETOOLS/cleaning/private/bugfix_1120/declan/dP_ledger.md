# Ledger dP – parser, GUI, root scripts (report for Declan, 11/21/2026)

Agent dP. Items: C13, C18, C19, C20, C24, C25, C26, C27, C28, C29, C30, C33, C34 of
`cleaning/private/release_comment_pass/CANDIDATES.md` (line numbers there: 78641839). HEAD at the
start: 905b0d26 (fa0deda4 + "Executives: never replace a requested form", which touches none of
these files). Since 78641839, among the files read here, only three changed: display_PDE.m and
initialize.m (0a1ed7f0, DJ: a skip for states without an equation at display_PDE.m:177-180; the
interior-variation error of BCs removed at initialize.m:624-635) and exec_sop_settings.m (4f44a021:
message text). No cited line changed; the line numbers that CANDIDATES.md gives for display_PDE.m
and initialize.m are those of HEAD, not of 78641839 (e.g. display_PDE.m:553 was 548). Every other
cited line number holds at HEAD. pietools_update.m: `git diff -w 7ced3647 HEAD` has no code line.

Scripts and outputs: scratch folder
`C:\Users\mpeet\AppData\Local\Temp\claude\363c499c6332a503\874c85cb-c18d-4876-bea3-39cdd053e21f\scratchpad\dP\`
(`<item>_repro.m` writes `<item>_repro.txt`). Git paths below are relative to the PIETOOLS folder.

Author evidence: `git blame -w -C -C -M`, `git log -S`, and the diff of the commit that introduced
the defect where the cited line itself is older (C18, C19, C24).

---

## C13 – NOT DJ (code by SS); defect CONFIRMED (run)

- HEAD: `pietools_path_update.m:39-40` (`current_folder = pwd;` / `addpath(genpath(current_folder));`), unchanged.
- Author: blame credits DJ, 7c435804 (2025-01-06 17:12, "Create pietools_path_update.m"). The two code lines are
  byte-identical (trailing space after `pwd;`, no final newline) to adf9f3f9, Sachin Shivakumar, 2025-01-06 13:53,
  "Create pietools_path_update.m", branch NonCoerciveLyap (not an ancestor of HEAD, hence plain blame); 7c435804
  added only the three comment lines above them. `git log --all -S'current_folder = pwd'` finds these two commits only.
- Reproduction: `C13_repro.m` (session-local path changes only). Measured:
  - A: default path + the PIETOOLS folder, pwd = a user folder, `pietools_path_update` by name: the user folder and its
    subfolder are added, no PIETOOLS subfolder (parser 0, sopvar 0), `which('pde_var')` = "".
  - B: same, pwd = the PIETOOLS folder: parser 1, sopvar 1, user folder 0.
  - C: `run(<full name>)` from the user folder: as B (run changes to the folder of the script).
- Consequence: silent wrong result (wrong path, no message) when called from any folder other than the PIETOOLS folder.
- Entry (for the maintainer, not for Declan):
  **C13** `pietools_path_update.m:39-40` (SS, adf9f3f9, 2025-01-06; entered the main line with DJ, 7c435804, the same
  day) – adds pwd and its subfolders, not the folder of the file. Called by name from a user folder with only the
  PIETOOLS folder on the path, it added the user folder and its subfolder and no PIETOOLS subfolder
  (`which('pde_var')` empty); called from the PIETOOLS folder, or with run and the full file name, it added the
  PIETOOLS tree. Silent wrong result (wrong path, no message), reached from any folder other than the PIETOOLS
  folder. Fix: `addpath(genpath(fileparts(mfilename('fullpath'))))`.
  Reproduction:
  ```
  restoredefaultpath; addpath(root);            % root: the PIETOOLS folder
  mkdir work; cd work; pietools_path_update;
  contains(path,[pwd,';']), which('pde_var')    % 1, '' (no PIETOOLS subfolder added)
  ```

## C18 – CONFIRMED (run)

- HEAD: `pietools_update.m:169-170` (error with `c`), `:179` (`[~, c] = urlwrite(...)`), unchanged since 7ced3647
  (`git diff -w 7ced3647 HEAD` has no code line).
- Author: line 170 is from the 2021b install script (blame -C -C: 1d9b0634, the 2021b import; in the 2020a script,
  c587fd94, `c = uigetdir(pwd)` at line 22 defined c). The error became possible in 7ced3647 (DJ, 2022-12-15,
  "Update install file."), which replaced `c = uigetdir(pwd)` by `download_dir = uigetdir(pwd)` in the fresh-install
  branch and kept the use of c. Defect: DJ.
- Reproduction: `C18_repro.m` – scratch copy of pietools_update.m (`make_update_stub.m`: the file text unchanged,
  stubs appended for websave, unzip, delete, movefile, rmdir, mkdir, urlwrite, tbxmanager, addpath, savepath, input,
  uigetdir, cd, and exist/which/fopen restricted; a probe checks that the local stubs shadow before the copy runs).
  Update branch, `exist('tbxmanager')` -> 0, mkdir -> false. Measured:
  - (a) no c: stopped at line 170: "Unrecognized function or variable 'c'." (no assignment of c before line 170;
    `which('c')` = "").
  - (b) c = 'C:\an\unrelated\folder' in the workspace: "An error appeared when trying to create the folder
    "C:\an\unrelated\folder". ..." – the folder that failed is d = <run>\PIETOOLS\tbxmanager.
- Consequence: error on a rare path (mkdir of PIETOOLS/tbxmanager fails, tbxmanager not installed); the intended
  message is replaced, or names an unrelated base-workspace value.
- Entry:
  **C18** `pietools_update.m:169-170` (DJ, 7ced3647, 2022-12-15; the line is from the 2021b install script) – the
  mkdir-failure error formats with c, which the script assigns only at line 179; 7ced3647 replaced the earlier
  definition `c = uigetdir(pwd)` by download_dir. On a scratch copy with stubbed downloads and a failing mkdir:
  "Unrecognized function or variable 'c'." at line 170; with a variable c in the workspace (the file is a script),
  the message names that value instead of the folder d. Error on a rare path (mkdir of PIETOOLS/tbxmanager fails
  while tbxmanager is not installed); the intended message is lost. Fix: d for c.
  Reproduction: C18_repro.m (scratch copy of pietools_update.m with every download, file and path call stubbed;
  `exist('tbxmanager')` -> 0, `mkdir` -> false).

## C19 – CONFIRMED (run)

- HEAD: `pietools_update.m:126` (`default_dir = pwd;`, inside the fresh-install branch, lines 115-150), `:234-235`
  (savepath fallback, after both branches); update branch lines 57-114 assign no default_dir.
- Author: lines 234-235 are 2021b text (blame -C -C: 1d9b0634); before 7ced3647, `default_dir = pwd;` was line 21 of
  pietools_install.m, executed unconditionally. 7ced3647 (DJ, 2022-12-15) created the update branch and moved
  `default_dir = pwd;` into the fresh-install branch (blame: line 126 = 7ced3647). Defect: DJ.
- Reproduction: `C19_repro.m` (stubbed copy as in C18; tbxmanager present; savepath() -> 1; uigetdir -> 0). Measured:
  - (U) update branch: stopped at line 234: "Unrecognized function or variable 'default_dir'." (no second savepath
    call in the stub log);
  - (F) fresh-install branch: no error; savepath(<start folder>\pathdef.m) called.
  Inferred (reading): in a real run the error also skips line 246, `cd(current_dir)`, so the current folder stays in
  the new PIETOOLS folder (cd at line 154).
- Consequence: error on a rare path: an update (PIETOOLS on the path) by a user who cannot write the default
  pathdef.m and cancels the folder dialog.
- Entry:
  **C19** `pietools_update.m:126, 234-235` (DJ, 7ced3647, 2022-12-15) – default_dir is assigned only in the
  fresh-install branch (line 126; 7ced3647 moved it there from the top of the script); the savepath fallback reads
  it in both branches. On the stubbed copy, savepath failing and the folder dialog cancelled: the update branch
  stops at line 234 with "Unrecognized function or variable 'default_dir'.", so the path is not saved (by reading,
  the current folder also stays in the new PIETOOLS folder: line 246 is skipped); the fresh-install branch saves to
  <start folder>\pathdef.m. Error on a rare path: an update by a user who cannot write the default pathdef.m and
  cancels the dialog. Fix: current_dir (line 54, the start folder in both branches) for default_dir at lines 234-235.
  Reproduction: C19_repro.m (stubbed copy of pietools_update.m; savepath() -> 1, uigetdir -> 0).

## C20 – CONFIRMED (run)

- HEAD: `parser/internal/display_PDE.m:553-643` (C block of local construct_term; use_Cij at 562, 576, 593), caller
  line 275 (two outputs), help line 8 (`name` "not currently used...").
- Author: DJ – the if at 553: d295c40f (2022-07-08, "Upload initialization for new terms format."); the use_Cij
  assignments: 74790c13 (2022-11-29, "Update display_PDE."); caller line 275 and the signature at 322: 423134af
  (2026-02-17, "Update PDE display.").
- Reproduction: `C20_repro.m`. Measured:
  - (a) parser route: after initialize, 0 of 4 factors (x_t = x_ss + s*x with Dirichlet BCs) and 0 of 7 factors
    (x_t = x_ss + 5x - x^3) lack C; display_PDE runs without error on both.
  - (b) C removed from PDE.x{1}.term{1}, is_initialized set back to true (the edit had reset it to 0): "Output
    argument "use_Cij" (and possibly others) not assigned a value in the execution with "display_PDE>construct_term"
    function." at display_PDE:275.
- Consequence: error on a rare or unreachable path: initialize sets C for every term (initialize.m:1534-1541,
  2068-2082 default to the identity), and any edit of a pde_struct resets is_initialized (subsasgn.m:164-166), so the
  path needs a PDE edited by hand and marked initialized. The unused input `name` is documented (help line 8).
- Entry:
  **C20** `parser/internal/display_PDE.m:553-643`, local construct_term (DJ, d295c40f, 2022-07-08; use_Cij
  74790c13, 2022-11-29; caller line 275, 423134af, 2026-02-17) – use_Cij is assigned only inside
  `if isfield(PDE_term,'C') && ~isempty(PDE_term.C)`, and line 275 requests it. A term without C in a PDE marked
  initialized ends in "Output argument "use_Cij" (and possibly others) not assigned a value in the execution with
  "display_PDE>construct_term" function." Error on a rare or unreachable path: initialize gives every term a C
  (0 of 4 and 0 of 7 factors without C for a linear and a nonlinear PDE, both displayed without error), so only a
  PDE edited by hand and marked initialized reaches it. The input name is unused (help, line 8, says so). Fix:
  `use_Cij = false;` before the if, as in display_term.m:306.
  Reproduction:
  ```
  pvar s t; x = pde_var(s,[0,1]);
  PDE = initialize([diff(x,t)==diff(x,s,2)+s*x; subs(x,s,0)==0; subs(x,s,1)==0],true);
  PDE.x{1}.term{1} = rmfield(PDE.x{1}.term{1},'C');   PDE.is_initialized = true;
  display_PDE(PDE)      % Output argument "use_Cij" ... not assigned
  ```

## C24 – CONFIRMED (run); DJ for PIETOOLS_DDE.m:282, NOT DJ for PIETOOLS_DDF.m:272

- HEAD: `PIETOOLS_DDE.m:282`, `PIETOOLS_DDF.m:272` (`sosineq_on=0; % binary variable indicating whether to use
  ineqaulity or equality constraint`), the only uncommented lines of each script that name sosineq_on; settings from
  `lpisettings_sop('light')` at DDE:284, DDF:274.
- Author: DDE:282 – 9abf156e (DJ, 2025-01-06, "Add H2 norm DDE example."), which replaced
  `settings.sosineq_on=0;` (after `settings = lpisettings('light');`) by the variable `sosineq_on=0;` before that
  line. DDF:272 – 2021b text (blame: 1d9b0634, import; in 2021b the settings scripts read the workspace variable);
  it became unread when DanBraghini replaced `settings_PIETOOLS_light;` by `settings = lpisettings('light');`
  (0d2deec4, 2022-12-03). NOT DJ for the DDF line.
- Reproduction: `C24_repro.m` (the scripts are not run: both call an example library that prompts with input()).
  Measured: (a) `isfield(lpisettings_sop('light'),'sosineq_on')` = 0; (b) `exec_sop_settings` with
  settings.sosineq_on = 0: "settings.sosineq_on was retired; remove it: the container executives pose no sosineq.";
  (c) the statement of each line evaluated in the base workspace: `settings_2d_stock('light').use_sosineq` 1 -> 0,
  ineq_opts replaced by eq_deg (settings_2d_stock.m:144-145 read the base workspace).
- Consequence: dead code or no effect on the DDE (DDF) run; the variable left in the base workspace changes a 2-D
  settings struct built later in the session by settings_2d_stock (read only with sop.slack = 'stock').
- Entry:
  **C24** `PIETOOLS_DDE.m:282` (DJ, 9abf156e, 2025-01-06) – sosineq_on = 0 is a script variable: 9abf156e turned
  `settings.sosineq_on=0;`, placed after `settings = lpisettings('light')`, into this variable before that line. No
  executive reads it: the settings come from lpisettings_sop (line 284), which has no such field, and the container
  executives reject the field ("settings.sosineq_on was retired; remove it: the container executives pose no
  sosineq.", exec_sop_settings.m:133-140). The variable stays in the base workspace, where settings_2d_stock reads it
  (settings_2d_stock.m:144-145): once the line has run, `settings_2d_stock('light').use_sosineq` is 0 instead of 1.
  No effect on the DDE run; it changes a 2-D settings struct built later in the session. The same line in
  PIETOOLS_DDF.m:272 is 2021b text, unread since 0d2deec4 (DB, 2022-12-03). Fix: delete the line.
  Reproduction:
  ```
  s0 = settings_2d_stock('light');   evalin('base','sosineq_on=0;');    % the statement of line 282
  s1 = settings_2d_stock('light');   [s0.use_sosineq, s1.use_sosineq]   % [1 0]
  st = lpisettings_sop('light');   st.sosineq_on = 0;   exec_sop_settings(st)   % error: retired
  ```

## C25 – NOT DJ (2021b install script); defect CONFIRMED (run)

- HEAD: `pietools_update.m:196-207` – `fid = fopen(p,'a')` (196), `fclose(fid)` only inside `if ~pathexist` (203-207).
- Author: blame -w -C -C: lines 191-207 from 1d9b0634 (the 2021b import). The oldest version on this branch with the
  block, 90ce5fa5 (Sachin, 2021-09-29, "Replaced files with PIETOOLS 2021b", itself an import of the 2021b files),
  already has fclose inside the if; the 2020a script (c587fd94, Sachin, 2020-07-06, branches PIETOOLS-2025 and
  others) had `fclose(fid)` after the if. 7ced3647 (DJ) kept the block unchanged. NOT DJ; original author of the
  move not established (2021b release).
- Reproduction: `C25_repro.m` (stubbed copy; which('startup.m') -> a startup.m in the scratch folder; update branch,
  tbxmanager present, savepath succeeds). Measured: (a) startup.m holding "tbxmanager restorepath": the script runs
  to the end, handle 4 (the scratch startup.m) still open (openedFiles), file unchanged; (b) without the line: the
  line is appended, no handle left open.
- Consequence: dead code or no effect on results (a file handle left open until fclose('all') or exit; every rerun
  of the update on an installation whose startup.m already holds the line).
- Entry (for the maintainer, not for Declan):
  **C25** `pietools_update.m:196-207` (2021b install script; kept by DJ in 7ced3647) – startup.m is opened in append
  mode and closed only when the line "tbxmanager restorepath" is absent. On the stubbed copy, with the line present,
  the script ends with the handle of startup.m open; without it, the line is appended and the handle closed. No
  effect on results (a file handle left open), on every rerun after the first installation. Fix: `fclose(fid)` after
  the if, as in the 2020a script.
  Reproduction: C25_repro.m (stubbed copy of pietools_update.m; a startup.m holding the line).

## C26 – CONFIRMED (reading)

- HEAD: `parser/@pde_struct/expand_delays.m:547`, local add_delay_state; row appended at 545-546; re-initialization at
  line 168; `parser/@pde_struct/initialize.m:530-533` (BC_tab rebuilt from zeros), `641-642` (BC_diff_tab appended).
- Author: DJ, 00a638da (2022-10-16, "Update pde_struct files."); `git log -S` finds no earlier commit with the line.
- Reading: a single subscript indexes a matrix in column-major order. With m rows and k = 2+nvars+tau_var_idx
  (tau_var_idx <= nvars, the table has 2+2*nvars columns, so k never exceeds numel), the statement sets row
  mod(k-1,m)+1, column floor((k-1)/m)+1: for m >= k the ID (column 1) of BC row k, for m = 1 the intended cell
  (column k of the appended row, as the comment says). No effect: between line 547 and line 168, BC_tab entries are
  only moved or appended (369, 545-546); initialize sorts BC_tab rows by column 1 (283-284) and discards the table at
  line 530 (PDE.BC is not reordered).
- Consequence: dead code or no effect (an overwritten cell of a table that is rebuilt).
- Entry:
  **C26** `parser/@pde_struct/expand_delays.m:547`, local add_delay_state (DJ, 00a638da, 2022-10-16) –
  `PDE.BC_tab(2+nvars+tau_var_idx) = 1` is a single (column-major) index, so it sets a cell of row mod(k-1,m)+1
  (for m >= k the ID in column 1 of row k), not column k of the row appended at lines 545-546, except for a one-row
  table. No effect: expand_delays re-initializes the PDE (line 168), and initialize rebuilds BC_tab from zeros
  (initialize.m:530-533, 641-642); in between the entries are only moved, appended, or sorted in a table that line
  530 discards. Fix: `PDE.BC_tab(end,2+nvars+tau_var_idx) = 1;`, or delete the line.
  Reproduction: by reading – a single subscript indexes in column-major order; BC_tab is rebuilt by initialize.

## C27 – CONFIRMED (reading)

- HEAD: `parser/pde_var.m:293-346`, local getPDEvarID; the only call, line 258, is commented (`%idx = getPDEvarID(obj);`),
  followed by `idx = stateNameGenerator;` (259).
- Author: DJ, 02f99523 (2024-06-23, "Add declaration of PDE vars."), which added the function and the call already
  commented out (diff lines `+%idx = getPDEvarID(obj);`, `+idx = stateNameGenerator;`).
- Reading: a local function is visible only within its file; no other call and no handle `@getPDEvarID` exists
  (grep, all .m files outside cleaning/).
- Consequence: dead code since its creation.
- Entry:
  **C27** `parser/pde_var.m:293-346`, local getPDEvarID (DJ, 02f99523, 2024-06-23) – no caller: the only call (line
  258) has been commented out since the commit that added the function, with stateNameGenerator in its place; a
  local function is visible only within its file, and no handle to it is created. Dead code. Fix: delete
  getPDEvarID.
  Reproduction: by reading – the call at line 258 is a comment; local functions are file-scoped.

## C28 – CONFIRMED (run)

- HEAD: `parser/@pde_struct/times.m:56, 58` and `mtimes.m:65, 67, 84` (the lines that name 'opvar'); branch
  `mtimes.m:84-185`.
- Author: DJ – times.m:56-59 and mtimes.m:65-67: a85b45da (2026-02-17, "Update pde_var multiplication.");
  mtimes.m:84-185: 4e29b50f (2024-06-29, "Allow multiplication with opvar.").
- Reproduction: `C28_repro.m`. Measured after pietools_path_update: `which('opvar','-all')` empty,
  `exist('opvar','class')` = 0, `meta.class.fromName('opvar')` empty; copvar, cdopvar, sopvar, sdopvar have no
  superclass. The only opvar classdef in the checkout is `cleaning/private/opvar/opvar_class/opvar.m`, under a private
  folder that genpath skips.
- Consequence: dead code: isa(.,'opvar') is false for every input, so the branch mtimes.m:84-185 is never taken, and
  times.m:58 and mtimes.m:67 act as else.
- Entry:
  **C28** `parser/@pde_struct/times.m:56-59`, `parser/@pde_struct/mtimes.m:65-67` (DJ, a85b45da, 2026-02-17) and
  `mtimes.m:84-185` (DJ, 4e29b50f, 2024-06-29) – the tests isa(...,'opvar') name a class that is not on the path
  (`which('opvar','-all')` empty, `exist('opvar','class')` = 0 after pietools_path_update), so they are false for
  every input: the branch mtimes.m:84-185 is never taken, and times.m:58 and mtimes.m:67 act as else. Dead code.
  Fix: delete the opvar tests and the branch mtimes.m:84-185.
  Reproduction: `which('opvar','-all')`, `exist('opvar','class')` – empty, 0.

## C29 – CONFIRMED (reading)

- HEAD: `parser/@pde_struct/times_pdes.m:180`, local times_term: `obj2_size = PDE.(obj1){obj2_num}.size;`; obj1, obj2
  set at 164-168, 174-178.
- Author: DJ, b8fad485 (2026-02-13, "Create mtimes_pdes.m"); `git log -S` finds no earlier commit with the line.
- Reading: obj1 = 'x' and obj2 = 'x' are the only values (lines 167, 177 raise an error otherwise), so
  PDE.(obj1){obj2_num} is PDE.(obj2){obj2_num}; `if ~strcmp(obj1,obj2)` (185) is never true.
- Consequence: dead code or no effect.
- Entry:
  **C29** `parser/@pde_struct/times_pdes.m:180`, local times_term (DJ, b8fad485, 2026-02-13) – obj2_size reads
  PDE.(obj1){obj2_num}.size; obj1 and obj2 are both 'x' (lines 164-168 and 174-178 raise an error for any other
  object), so the value is that of PDE.(obj2){obj2_num}, and `if ~strcmp(obj1,obj2)` (line 185) is never true. No
  effect while only states can be multiplied. Fix: obj2 for obj1 at line 180.
  Reproduction: by reading – lines 164-178 admit only 'x' for both factors.

## C30 – CONFIRMED (run)

- HEAD: `parser/@pde_struct/combine_vars.m:194`: `fprintf('\n','No variables have been merged.\n')`.
- Author: DJ, ed59ca7c (2022-10-17, "Update "pde_struct" files."); `git log -S` finds no earlier commit.
- Reproduction: `C30_repro.m`. Measured: (a) the statement prints one character, "\n"; (b) combine_vars on a PDE in
  one variable prints "\n\nAll spatial variables have been rescaled to exist on the interval [-1,1].\n", without the
  message; (c) with three mergeable variables (manual, Sec. 8.2.4d) the merge message prints.
- Consequence: misleading (missing) message: every combine_vars call without a merge and without silent_merge, also
  through expand_delays (line 190, combine_vars(PDE)).
- Entry:
  **C30** `parser/@pde_struct/combine_vars.m:194` (DJ, ed59ca7c, 2022-10-17) – `fprintf('\n','No variables have been
  merged.\n')` prints one newline: the format has no conversion, so the text argument is not printed. combine_vars on
  a PDE in one variable printed only the rescaling line. Message lost in every call without a merge and without
  silent_merge, also from expand_delays (line 190). Fix: `fprintf('\nNo variables have been merged.\n')`.
  Reproduction:
  ```
  pvar s t;   x = pde_var(s,[0,1]);
  PDE = initialize([diff(x,t)==diff(x,s,2); subs(x,s,0)==0; subs(x,s,1)==0],true);
  PDE2 = combine_vars(PDE);    % no "No variables have been merged."
  ```

## C33 – CONFIRMED (run)

- HEAD: `parser/@pde_struct/display_PDE_crude.m:177, 179, 199` (display_term with two outputs), 197 (three outputs),
  208 (`use_C = use_C || use_Cjj`), 216-217 (hint).
- Author: DJ – 177-179, 197-199: 423134af (2026-02-17, "Update PDE display."); 208-216: 4324069c (2024-06-23,
  "Update PDE structure display."); 217: 8fd693bd (2026-11-06, "Update display functions.").
- Reproduction: `C33_repro.m` (uninitialized equations, displayed by display_PDE_crude via display.m:65-68; the same
  with Dirichlet BCs, initialized, displayed by display_PDE). Measured, hint printed (crude / initialized):
  - x_t = s*x_ss, C11(s) in term 1: 0 / 1;
  - x_t = x_ss + s*x, C12(s) in term 2: 1 / 1;
  - x_t = x_ss + x*(s*x), C122(s) in factor 2 of term 2: 0 / 1 (display_PDE: 'Call "PDE.C{i,j,k}"');
  - x_t = x_ss + (s*x)*x, C121(s) in factor 1 of term 2: 1 / 1, but the crude hint names PDE.C{i,j} for a coefficient
    displayed with three indices (display_PDE names PDE.C{i,j,k}).
- Consequence: misleading (missing) message when an unfinished PDE is displayed.
- Entry:
  **C33** `parser/@pde_struct/display_PDE_crude.m:177, 179, 199` (DJ, 423134af, 2026-02-17) – only line 197 requests
  use_Cjj, for the first factor of terms 2 to n; the first term and every later factor are not counted. On
  uninitialized PDEs, x_t = s*x_ss (shown with C11(s)) and x_t = x_ss + x*(s*x) (C122(s)) print no hint 'Call
  "PDE.C{i,j}"', while display_PDE prints it for both once initialized; x_t = x_ss + s*x prints it. When printed for
  a product term (x_t = x_ss + (s*x)*x, C121(s)), the hint names PDE.C{i,j}, where display_PDE names PDE.C{i,j,k}.
  Missing or misleading hint when an unfinished PDE is displayed (display.m:65-68). Fix: request the third output at
  lines 177, 179 and 199 and OR it into use_C, as at lines 197 and 208.
  Reproduction:
  ```
  pvar s t;   x = pde_var(s,[0,1]);
  E = diff(x,t)==s*diff(x,s,2)     % shows C11(s), no 'Call "PDE.C{i,j}"' line
  ```

## C34 – CONFIRMED (run)

- HEAD: `parser/@pde_struct/mpower.m:46-47`; `parser/@pde_struct/power.m:47-48` ("nonnegative integers"), `50-51`
  (degree 0 refused).
- Author: DJ, e5b96e56 (2026-02-18, "Update PDE power functions.") for mpower.m:46-47 and power.m:47-51.
- Reproduction: `C34_repro.m`. Measured: x^2 no error; x^0 and x.^0 -> "Degree-0 powers of PDE terms are currently
  not supported." (power:51); x^(-1), x^1.5 -> "Degrees must be specified as real nonnegative integer." (mpower:47);
  x.^(-1) -> "Degrees must be specified as real nonnegative integers." (power:48).
- Consequence: misleading message: both checks state "nonnegative", and degree 0 is refused.
- Entry:
  **C34** `parser/@pde_struct/mpower.m:46-47` (DJ, e5b96e56, 2026-02-18) – the check accepts deg = 0 and its message
  asks for a "real nonnegative integer"; power rejects 0 (power.m:50-51) after a first check whose message also says
  "nonnegative" (power.m:47-48). x^(-1) gives "Degrees must be specified as real nonnegative integer.", x^0 gives
  "Degree-0 powers of PDE terms are currently not supported.". Misleading message: the stated admissible degrees
  include 0, which is refused. Fix: test deg<1 in mpower and ask for a positive integer in both messages.
  Reproduction:
  ```
  pvar s;   x = pde_var(s,[0,1]);
  x^(-1)    % "... real nonnegative integer."
  x^0       % "Degree-0 powers of PDE terms are currently not supported."
  ```

---

## Summary

| Item | Status | Consequence class |
|---|---|---|
| C13 | NOT DJ (SS, adf9f3f9); defect confirmed (run) | silent wrong result (wrong path, no message) |
| C18 | CONFIRMED (run) | error on a rare path (message replaced) |
| C19 | CONFIRMED (run) | error on a rare path |
| C20 | CONFIRMED (run) | error on a rare or unreachable path |
| C24 | CONFIRMED (run); DDE line DJ, DDF line NOT DJ | dead code or no effect (base-workspace side effect) |
| C25 | NOT DJ (2021b install script); defect confirmed (run) | no effect (file handle left open) |
| C26 | CONFIRMED (reading) | no effect |
| C27 | CONFIRMED (reading) | dead code |
| C28 | CONFIRMED (run) | dead code |
| C29 | CONFIRMED (reading) | no effect |
| C30 | CONFIRMED (run) | misleading (missing) message |
| C33 | CONFIRMED (run) | misleading (missing) message |
| C34 | CONFIRMED (run) | misleading message |

Runs (MATLAB R2025b, lock runner, 11/21/2026 02:45-02:54, each process under 20 s): parser group (C20, C24, C28,
C30, C33, C34), update group (C18, C19, C25), C13 alone (it resets the path), C18 again standalone with the final
helpers. Nothing was written outside the scratch folder: `git status` shows no new file from these runs; the stubs
refuse fopen outside the scratch folder and replace which('startup.m').

Rerun: each script runs alone, e.g.
`bash .../scratchpad/mlrun.sh "run('C:\...\scratchpad\dP\C18_repro.m')"`; or a group with the driver:
`bash .../scratchpad/mlrun.sh "dP_list = {'C20','C24','C28','C30','C33','C34'}; run('C:\...\scratchpad\dP\dP_run.m')"`.
C18, C19 and C25 use `make_update_stub.m` and `run_update_stub.m` in the same folder. Run C13 in its own process.
