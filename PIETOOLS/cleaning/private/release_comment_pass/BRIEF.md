# Release comment pass – rules for each group (11/19/2026)

Sources: `PIETOOLS_documentation/for_developers.txt` ("Special Rules for public releases"), the developer
annotation protocol (repository CLAUDE.md, section 1), the header protocol of 10/28/2026, the comment-style
rule (terse), and the maintainer decisions of 11/19/2026 (other authors: rules only, no restyle; SOSTOOLS400
excluded; "See also" lines kept; PIETOOLS_demos/cuadmm and lowrank_2d_stability skipped).

All paths are relative to `C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS` (call it P).
Tools in `P\cleaning\private\release_comment_pass\`: `codeident.sh` (code of a file with comments removed,
compared with HEAD; prints OK, DIFF or NEW), `codeident.awk`.

## Scope

Your file list `lists/<group>.tsv` (columns: path, new). Edit nothing else. Never touch SOSTOOLS400/,
cleaning/, orphaned_file/, testing_codes/, PIETOOLS_demos/cuadmm/, PIETOOLS_demos/lowrank_2d_stability/.
Other sessions edit this checkout: right before editing a file run `git -C <repo> status --porcelain -- PIETOOLS/<path>`;
if it prints anything, skip the file and record it in the ledger as "skipped: uncommitted edits".

## Hard rule: comment-only

Never change code: no code token, no indentation or spacing inside a code line, no line split or join, not
even an obvious bug. Functional issues go to the ledger. After editing a file run
`bash P/cleaning/private/release_comment_pass/codeident.sh <path>`: it must print OK. If it prints DIFF,
restore the code (only comments may differ from HEAD) before moving on.

## Rules per file

1. **Header** – every function file, every script, and every local function (subroutine):
   a. what it does, with enough detail and context to apply it outside its original use; when it
      performs a mathematical operation, the formula (plain-text notation, as the file already uses);
   b. INPUT: every input – its class or structure (class in quotes: `'copvar'`, `'struct'` with its fields,
      `'double'` m x n, ...), optional inputs and defaults – related to the description;
   c. OUTPUT: every output, likewise.
   Keep the file's existing header placement (after the function line in most files; above it in
   PIESIM-style files). A file without a header gets one directly after the function line, in the form
   `% [OUT1,OUT2] = NAME(IN1,IN2) <what it does>`, then `% INPUT` and `% OUTPUT` with one entry per
   argument (model: `executives_sopvar/utility_functions/exec_info_sop.m`). Local functions: a short
   what/INPUT/OUTPUT block, no license. Keep the help block contiguous so `help NAME` prints it first.
2. **License** – every file (not local functions) carries the PIETOOLS GPL block. If missing, copy it
   exactly from `executives_sopvar/PIETOOLS_stability_sop.m`, from `% Copyright (C) 2026 PIETOOLS Team`
   through `% authors, and a brief summary of the modification`. Keep an existing block as it is.
3. **Authors** – at the end of the header (after the license): names, contact, date of initial coding.
   Keep existing `% For support, contact ...` and `% Initial coding ...` lines. Where missing, add
   `% For support, contact M. Peet, Arizona State University at mpeet@asu.edu` and
   `% Initial coding <initials> - <mm/dd/yyyy>`, taken from
   `git -C <repo> log --follow --diff-filter=A --format="%an %ad" --date=format:%m/%d/%Y -- PIETOOLS/<path> | tail -1`
   (MATTHEW M PEET -> MMP, Declan or Declan Jagt -> DJ, Talitsky -> AT, DanBraghini -> DB,
   Sachin Shivakumar -> SS; any other name as given).
4. **Developer notes**
   - **new = 1** (not in the 2025 public release): remove ALL developer notes – the change log (dated
     entries such as `% DJ, 11/16/2026: ...`, `% 05/04/23 - DJ: ...`), per-line markers at line ends
     (`% MMP, 11/14/2026`, with or without `(was)`), commented-out old code (the `(was)` lines), and
     `Start/End of changes` section markers. Keep the license block (its "If you modify this code ..."
     lines included) and the `Initial coding` line. Where a log entry holds a fact a user needs (a
     limitation, a convention, a sign or ordering rule), restate that fact in the description or a step
     comment as what the code does – never as history.
   - **new = 0** (in the 2025 release): keep the change log, the markers and the `(was)` lines; append one
     log line at the end of the log: `% MMP, 11/19/2026: Release comment pass (comment-only).` No per-line
     markers for comment edits.
5. **Body** – every step of every function body has a terse comment saying what is done, in terms of the
   header description. Add missing step comments; keep existing ones.
6. **References** – descriptions and comments do not point to external documents or other files or
   functions: remove references to the manual, the Overleaf spec or its sections, README/NOTES/.md/.tex/.pdf
   files, tests, session notes, cleaning/ paths, and phrases such as "as in X.m", "see Y", "called by Z",
   "used by the executives". Replace the pointer with what the code does. Keep `See also` lines (decision).
   A literature citation that attributes an example or a method to its source may stay; record it in the
   ledger. Describe what the code does, not how it is used. Brief notes on observed performance,
   bottlenecks or limitations may stay.
7. **Wording of any text you write** – terse: state the fact and stop (one line where one suffices; no
   "Note that"; no restating the code). No "I", no possessives ('s), spaced en dashes ` – `. No history
   framing ("moved from", "formerly", "the container version of", "replaces opvar"): write as if the retired
   opvar classes never existed. Psatz terms are named by their factor (plain term, box product,
   face (s_d - a_d), box quadratic, face products), never a bare "face" or "plain". Never "no storage" (the
   Q form is a generalized storage T'Q). Only the 'P' form with epneg > 0 certifies exponential stability.
   No test numerics or cleaning/ references in example comments.
8. **Other authors** (files whose header names an author other than MMP, or whose comments carry other
   initials): apply rules 1-6 only where a rule is violated; never reword their comments for style.
9. **Never** run MATLAB, never run git add/commit/stash/checkout/reset or anything that changes git state,
   never edit outside your list.

## Ledger and report

Write `lists/<group>_ledger.md`: one line per edited file (what changed); "Skipped" (with reason);
"Candidates for the maintainer" (functional issues seen: bugs, dead code, misleading code – not fixed);
"Unresolved" (headers or formulas you could not state with confidence). Run codeident.sh on every edited
file at the end and put the counts (OK / DIFF / NEW) in the ledger. Report back in under 300 words: files
read, files edited, codeident counts, the number of candidates and unresolved items, and anything that
needs the maintainer.
