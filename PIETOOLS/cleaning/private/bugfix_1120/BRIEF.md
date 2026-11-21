# Bug-fix pass, 11/20/2026 – brief for every agent

Maintainer, 11/20/2026 ~19:00: "fix the unfixed bugs" – the items of
`cleaning/private/release_comment_pass/CANDIDATES.md` (read its header: evidence conventions, authors).
Paths below are relative to the PIETOOLS folder
`C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS` (git paths carry a `PIETOOLS/` prefix).
HEAD is 17d68747; CANDIDATES.md line numbers refer to 78641839 and may have shifted in the 21 files
changed by 0a1ed7f0 (Declan, algebraic constraints). Re-read every cited line at HEAD.

## Scope and holds

Fix only the items assigned to you. Never edit, even in passing:
- code of Declan Jagt (DJ): any item credited DJ, and any line whose author is DJ by `git log -L` / `git blame`;
- all of `PIESOS/` (Declan owns it), including `PIESOS/polyopvar/@quadPoly`; `orphaned_file/`;
- `SOSTOOLS400/`;
- behaviour of the 2-D settings mapping: `settings_sopvar/settings_2d_stock.m`,
  `executives_sopvar/utility_functions/settings2possopvar.m`, the psatz mapping in
  `poslpivar_settings_2d_sop.m`, the stock 2-D executives (text-only message fixes are allowed where
  your item list says so);
- files with uncommitted edits by another session. Before the first edit of a file run
  `git status --porcelain -- <file>`; if it shows a change you did not make, skip the item and record
  it. Known: `PIETOOLS_demos/DEMO2_volterra_operator_norm.m`, `PIETOOLS_demos/DEMO15_trace_inequality_2D.m`.

If a fix would need an edit in a held file, record the item as HELD with the reason; do not edit.

## Per item

1. Re-read the cited lines at HEAD.
2. Reproduce. Where the path is reachable, write a minimal MATLAB script that shows the defect (wrong
   value, error, lost message) BEFORE editing, and record the measured output. For semantics-only items
   (dead line, unreachable branch) state the MATLAB rule; no run needed. A hedged item ("not determined
   by reading", "unless") is a hypothesis: test it; if it does not reproduce, record NOT REPRODUCED and
   do not edit. Keep "measured" and "inferred" apart in the ledger.
3. Fix with the minimal change (CLAUDE.md s.1): fewest changed lines; no refactor, no new helper file,
   no rename of anything a user calls or reads (functions, options, output fields). If the only correct
   fix is not minimal, or the intent of the code cannot be established from the code, its help, the
   user manual (`PIETOOLS_documentation/`) or git history (`git log -L`, the 2025 release
   `origin/PIETOOLS-2025`), record NEEDS DECISION with the options and do not edit.
4. Verify. Rerun the reproduction and check the corrected output against an independent reference: the
   semantics of the class or the manual, a closed form, or a second route that is not the inverse of the
   first. Never test mutually inverse routines against each other.
5. Record the item in your ledger.

## Annotation – two regimes

**New file** = its path is listed in column 4 of
`cleaning/private/release_new_files_1120/new_files_classified.tsv` (any class). Release rule
(`PIETOOLS_documentation/for_developers.txt`, "Special Rules for public releases"): no developer notes
in new functions – no header log line, no dated per-line markers, no (was) lines. Change the code; where
the corrected line is not obvious, add a terse comment saying what it does.

**Existing file** = every other file. Full protocol (CLAUDE.md s.1):
- one header log line in the format the file already uses, after the last log line (usually
  `% MMP, 11/19/2026: Release comment pass (comment-only).`), e.g.
  `% MMP, 11/20/2026: <what changed and why>.` (two lines at most);
- each changed line ends with `% MMP, 11/20/2026`, the `%` at column 77 (pad with spaces; if the code
  reaches column 75 or beyond, one space then the marker); never on a line ending in `...`;
- a minor change keeps the old line commented out above the new one, ending `% MMP, 11/20/2026 (was)`;
- a larger change gets BEGIN/END marker lines.

Initials MMP throughout. Edit with the Edit tool only – never `sed -i` on the working tree (it rewrites
CRLF files as LF). After editing, `git ls-files --eol -- <file>` must show the same `w/` state as before.
Check marker columns with
`awk '{i=index($0,"% MMP, 11/20/2026"); if(i) print FILENAME":"FNR": col "i}' <file>`.

## Wording of comments and messages

Terse, what and why. No "I", no possessives ('s), spaced en dashes ` – `, never em dashes. User-visible
strings name no retired class or function (opvar, opvar2d, dopvar, dopvar2d, lpisettings, poslpivar,
poslpivar_2d, lpi_ineq_2d, mat2opvar, possopvar), no dates, no unpublished reports, no cleaning/ paths.
Do not restyle comments of other authors beyond the item.

## MATLAB

- The maintainer is online: every MATLAB run under 5 minutes; no large SDPs.
- At most ONE batch MATLAB at a time across all agents. Run every MATLAB command through the lock:
  `bash /c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/mlrun.sh "run('C:\...\yourscript.m')"`
  with Bash timeout 600000. Exit code 75 = lock busy for 240 s: retry. Batch several cases per script to
  save start-up time (~15 s), but one case group per process (R2025b occasionally crashes in long runs);
  write results to a file as they come.
- Start each script with
  `cd('C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS'); pietools_path_update;`
  The warnings "Name is nonexistent or not a directory: ...\sopvar\psatz_faces_results..." are a stale
  saved path; ignore them. Assert with `which` that each function under test resolves inside the
  working tree (not cleaning/, orphaned_file/, another checkout).
- Test scripts and outputs go in your own scratch folder
  `/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/<agent>/`,
  never in the repository. Library example files may call input(); call example files directly.
- At the end, one parse check over your edited files:
  `addpath('C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS\cleaning\private\release_comment_pass'); syntax_check('C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS', '<list file>', '<out file>')`
  (list file: one path per line relative to the PIETOOLS folder).

## Ledger and report

Ledger: `cleaning/private/bugfix_1120/<agent>_ledger.md`. Per item: status (FIXED / NOT REPRODUCED /
NEEDS DECISION / HELD), files and lines changed, reproduction (script, measured output before),
verification (measured output after, reference used), measured vs inferred. Update it after each item,
so that an interruption loses at most one item.

Do not commit, stash, reset, checkout or restore. Edit nothing outside your item list except your ledger
and scratch folder.

Final message: per item, status in one line; then the list of files edited; at most 40 lines.

## Coordinator review, 11/20/2026 (after the four agents)

Edits made in review:
- `lpi_programming_sopvar/poscopvar.m`, `lpi_programming_sopvar/misc/copquadvar.m` (C104): psatz 1 is
  named "the box product" (the vocabulary of the code), not "the box weight".
- Help blocks corrected in place, without (was) lines or markers (they print with `help` and, for the
  first line, `lookfor`): `CL/stability_7.m`, `CL/stability_8.m` (first line), `CL/stability_10.m`,
  `EL/2D/PIETOOLS_PDDE_Ex.m`, `PIESIM/.../Time_Integrators/PIESIM_int_symbolic.m` (output list).
- C116 (new): `CL/l2gain_2.m:40`, gamma = 0.3333 -> 0.2083, with the log line.

Independent checks (measured):
- C35: with the HEAD copy of convert_PIETOOLS_NDS2DDF first on the path, an NDS with the kernel
  2*int x(t+s) ds gives Cvd{1} = [0;0] ("Cvd{1} undefined, defaulting to 0"); the working tree gives
  [-0.9701;0], and the PIEs differ (T and A).
- C116: sup_w |G(i w)| = 0.208333 at w = 0 = 5/24 (Chebyshev collocation, N = 80, outside PIETOOLS).
- test_executives_sop(false): 20 of 20 checks, 2-D included, 17.6 s.
- syntax_check: 75 of 75 edited files, 0 parse problems; `git diff --numstat` equals `git diff -w
  --numstat` (no whitespace-only change); `git ls-files --eol` states as before the pass.
