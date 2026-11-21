# Comment pass over files of the 2025 release, 11/21/2026 – brief for every agent

Maintainer, 11/21/2026: "go through all comments in the repo which are not in new files. For these
modifications to pre-existing files, condense the developer notes to a single note (less than 10 lines).
Again, respect the rules on commenting. Comments should indicate changes with respect to previous release.
Internal notes pertaining to development and testing are not needed. These comments should be enough to
indicate what has changed. There are some files which have been completely changed, such as the settings
files, and these should be treated as new files."

Clarifications of the maintainer, 11/21/2026:
- leave PIESOS/ and parser/ alone;
- "leave commented code for code that existed in the prior release, but commented code created and removed
  in the internal development process can be removed. All XY date markers should match the single header
  note.";
- only the notes added after the 2025 release are condensed; the notes of the 2025 release and the Initial
  coding lines stay as they are;
- a file with less than 25% of its 2025 code left is treated as a new file;
- a file whose new notes are mainly by DJ or SS is flagged for the maintainer, not edited.

Paths are relative to the PIETOOLS folder `C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS`.
The 2025 release is `origin/PIETOOLS-2025` (fb7f1209, a `PIETOOLS/` prefix inside). Commenting rules: read
`PIETOOLS_documentation/for_developers.txt` (both lists) and the wording rules of
`cleaning/private/release_comment_pass/BRIEF.md`.

Tools (`$T` = `/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/cpass1121`,
`$RC` = `cleaning/private/release_comment_pass`):
- `perl $T/cpass.pl prep FILE` – the 2025 version, authors since 2025 and a proposed TAG, the header
  entries new since 2025, every new marked line with a verdict (RETAG, KEEP-WAS, DROP-WAS,
  REVIEW-BEGIN/END), the new comment lines without a marker, and the code diff against 2025. "New" means
  that no line of the 2025 version equals it.
- `perl $T/cpass.pl apply FILE "TAG"` – drops the new (was) lines whose code is not 2025 code (steps of
  development), keeps the (was) lines of 2025 code and retags them `% TAG (was)`, retags every other new
  marker `% TAG`, puts the markers at column 77; lines of the 2025 version are never touched.
- `perl $T/cpass.pl check FILE "TAG"` – PASS, or the failures.
- `bash $RC/codeident.sh FILE` – OK when the code with comments removed equals HEAD.

## Per file of class `condense`

1. Run prep.
2. TAG: the initials of everyone who changed the file since 2025, first change first (prep proposes them
   from the commits and the new notes; MMP covers every commit of MATTHEW M PEET), then `, 11/21/2026`, e.g.
   `YP, DJ, MMP, 11/21/2026`. The date is 11/21/2026 for every file.
3. Edit the header with the Edit tool: delete every new header entry with its continuation lines, among
   them "Release comment pass (comment-only)", and write the single note in their place, after the notes
   of the 2025 release:
   `% TAG: Changes since PIETOOLS 2025: <text>` and continuation lines `%   <text>`, 9 lines at most in all.
   The note says what changed with respect to the 2025 release, from the code diff of prep and the new
   entries: new, removed or renamed inputs, outputs, options, defaults, cases supported, behaviour, user
   messages. State the result, not the history: no dates of intermediate steps, no "moved", "replaces the
   change of", no quotes of the maintainer, no tests, studies, measurements, DECISIONS, cleaning/,
   sessions, reports; no function or file names unless the change is their addition, removal or renaming.
   A file changed in comments only: `help and comments revised`. A 2025 entry whose wording changed (its
   initials and date are those of a 2025 entry) stays; a (was) signature line at the end of the log is a
   marked line (step 4).
   BEGIN/END lines of a change in the body become `% TAG: BEGIN <what>` and `% TAG: END <what>`, one line
   each.
4. Run apply with the TAG. Compare its counts with the verdicts of prep; resolve REVIEW lines by hand.
5. Read the new comment lines of prep: remove or rewrite internal notes about development and testing
   (measurements, experiments, who asked for what, test cases, sessions, dates of trials); keep help text
   and step comments that say what the code does (release rules 2, 5, 6). Commented-out code added and
   removed since 2025 without a marker: remove; commented-out code that is in the 2025 version: keep.
6. Run check (PASS) and codeident (OK); `git ls-files --eol -- FILE` shows the same `w/` state as before;
   `git diff --numstat -- FILE` equals `git diff -w --numstat -- FILE` (no whitespace-only change).
7. Ledger entry (below).

## Class `asnew` (less than 25% of the 2025 code left)

Treat as a new file: remove every developer note from header and body – all log entries, those of the
2025 version included, all dated markers, all (was) lines and BEGIN/END lines – and keep the license, the
help header and, at the end of the header, the authors, contact and date of initial coding (release rules
1-4). Then codeident OK and no dated marker or log entry left except the Initial coding line.

## Class `review` (new notes mainly by DJ or SS)

Do not edit. Put in the ledger the TAG and note you would write, and the prep verdicts, for the maintainer.

## Class `deferred`

Another session has uncommitted edits in the file. Do not edit; record it.

## Never

Edit a file outside your list; touch PIESOS/, parser/, SOSTOOLS400/, testing_codes/, orphaned_file/, or
cleaning/ other than your ledger; change code (codeident must say OK). Before the first edit of a file run
`git status --porcelain -- FILE`: a change you did not make means skip and record. No commit, stash, reset,
checkout or restore. Hand edits with the Edit tool only, never `sed -i`; the Edit tool can drop trailing
blanks of unchanged lines in its strings – the numstat check of step 6 catches it.

## Parse check

At the end, one check over your edited files through the MATLAB lock runner (rules in
`cleaning/private/bugfix_1120/BRIEF.md`, section MATLAB):
`addpath('C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS\cleaning\private\release_comment_pass'); syntax_check('C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS', '<list file>', '<out file>')`
(list file in your scratch folder `$T/<agent>/`, one path per line). No other MATLAB run is needed.

## Wording

Terse; no "I", no possessives ('s), spaced en dashes ` – `, never em dashes. Other authors: no restyle
beyond the rules.

## Ledger

`cleaning/private/comment_pass_1121/<group>_ledger.md`. Per file: class, TAG, the note verbatim, apply
counts (dropped and kept (was) lines, retagged markers), internal notes removed, open questions. Update it
after each file.

Final message: counts per class, files skipped, open questions; at most 30 lines.
