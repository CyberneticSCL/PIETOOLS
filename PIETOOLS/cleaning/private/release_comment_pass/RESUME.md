# Release comment pass – resume (11/20/2026)

Read `BRIEF.md` first. Its rules hold, with the corrections and additions below. Paths are relative to
`P = C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS`; tools are in
`T = P\cleaning\private\release_comment_pass\`; your lists are in `T\lists\resume\`.

## Corrections to BRIEF.md

1. **License end line.** The GPL block ends with `% authorship, and a brief description of modifications`
   (the text in `executives_sopvar/PIETOOLS_stability_sop.m`; BRIEF.md misquotes it).
2. **Initial coding.** A first git commit dated 12/31/2021 is the import of PIETOOLS 2021b, not the coding
   date: never write `Initial coding DJ - 12/31/2021` or any 12/31/2021 date. For such a file keep the
   authors the header names and write `% Initial coding <initials> (PIETOOLS 2021b; in the repository
   since 12/31/2021)`; with no author named, `% Initial coding PIETOOLS Team (PIETOOLS 2021b; in the
   repository since 12/31/2021)`.
3. **new.** The flag in your list is final. It was recomputed with class-method keys (`@class/method.m`)
   against origin/PIETOOLS-2025; the first census compared basenames and marked 198 new files as
   existing. PIESOS/ is wholly new, except `PIESOS/polyopvar/int_simple.m` (moved from the 2025
   `opvar/2D/int_simple.m`; new = 0).

## Dirty-file check (replaces the last sentence of BRIEF.md, "Scope")

Files the pass already edited show as modified in git. Before the first edit of such a file run
`bash T/codeident.sh <path>`: OK – proceed; DIFF – another session changed its code: skip it and record it.
A file the pass has not edited (status "untouched" or "not yet edited" in your list) must be clean:
`git -C <repo> status --porcelain -- PIETOOLS/<path>` prints nothing; otherwise skip it and record it.

## newfix: new files that were passed as existing

They still carry developer notes and the line `% MMP, 11/19/2026: Release comment pass (comment-only).`
added on 11/19. Apply rule 4 for new = 1: remove the change log (every dated entry; keep the `Initial
coding` line), per-line markers at line ends (`% XY, mm/dd/yyyy`, with or without `(was)` or a note),
commented-out old code (`(was)` lines, including all continuation lines of an old statement),
Start/End and BEGIN/END markers, and the release log line. Keep the license block, the `For support` and
`Initial coding` lines, and `See also` lines. Trim the trailing blanks a removed marker leaves. A fact a
user needs that only a removed entry states (a limitation, a convention, a sign or ordering rule) is
restated as present-tense description in the header or a step comment. Then check the file against rules
1–7, local-function headers included. A script may strip markers and `(was)` lines; run codeident on
every file it touched.

## Step comments (rule 5): objective floor

`awk -v N=20 -f T/runs.awk <file>` prints `FILE:START-END:LENGTH` for each run of 20 or more code lines
with no comment (blank lines neither extend nor break a run; a comment line or a trailing comment breaks
it). After your edits no file on your list may have such a run. Break each run by commenting its steps –
what each block does, in terms of the header – never by restating single statements. Runs of 15–19
lines: by judgment. Starting list: `lists/resume/<agent>_runs.txt` (measured 11/20/2026 01:20).

## Factual errors in existing comments

A comment that contradicts the code is a rule 1 or rule 5 violation in any file, including other
authors' files: correct it to what the code does once verified from the code, and record it in the
ledger. Never change code to match a comment; a code defect goes to "Candidates".

## Ledger and report

Write `lists/resume/<agent>_ledger.md` (not the 11/19 group ledger): one line per edited file (what
changed); Skipped (with reason); Candidates for the maintainer (functional issues seen, not fixed);
Unresolved; final `codeident.sh` counts (OK / DIFF / NEW) and the final `runs.awk` (N = 20) count over
every file on your lists. Report back in under 300 words.

## Unchanged from BRIEF.md

Comment-only edits; codeident OK after every file; never run MATLAB; never run a git command that changes
state (add, commit, stash, checkout, reset, restore, merge); edit only files on your lists (another agent
owns every other file of the pass at the same time); other authors' files: rules only, no restyling;
terse wording; no "I", no possessives, spaced en dashes ` – `; no history framing (write as if the
retired opvar classes never existed); psatz terms named by their factor; the Q form is a generalized
storage T'Q, never "no storage"; only the 'P' form with epneg > 0 certifies exponential stability.
Keep each file's line endings (CRLF stays CRLF) and its final-newline state.
