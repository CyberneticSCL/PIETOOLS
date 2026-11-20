# Release comment pass – plan (11/19/2026)

## Instructions applied

1. `PIETOOLS_documentation/for_developers.txt`, "Special Rules for public releases": developer notes out of
   new functions; a header on every function and subroutine (what it does, with formulas; every input;
   every output); the public license in every file; authors, contact and date of initial coding at the end
   of the header; every step of the body commented in terms of the header; no references to external
   documents, other functions or files; describe what the code does, not how it is used.
2. Developer rules 1-5 of the same file and the repository CLAUDE.md section 1 (header log, per-line
   markers at column 77, `(was)` lines, start/end markers, what and why) – kept for files that were in the
   2025 public release.
3. Header protocol (maintainer, 10/28/2026): help block, then license, then history; INPUT/OUTPUT as data
   objects; subroutines too. Comment style: terse, coverage complete.
4. Wording feedback (11/18-11/19/2026): no history framing (opvar retired), psatz terms by factor, Q form as
   generalized storage, exponential stability only for 'P' with epneg > 0, no test numerics or cleaning/
   references in example comments.
5. Decisions of 11/19/2026: other authors' files get the rules only, no restyling; SOSTOOLS400 excluded
   (checked only for mentions of PIETOOLS); `See also` lines kept; PIETOOLS_demos/cuadmm and
   lowrank_2d_stability skipped.

Rules as given to the agents: `BRIEF.md`.

## Scope (census of 11/19/2026, origin/ndopvar 7f7a0300 checkout)

921 .m files outside SOSTOOLS400; "new" = function name absent from origin/PIETOOLS-2025.

| group | folders | files | new |
|---|---|---|---|
| A | sopvar | 178 | 94 |
| B | lpi_programming_sopvar, executives_sopvar, settings_sopvar | 70 | 65 |
| C | parser, GUI, root scripts | 127 | 13 |
| D | converters, PIESIM | 115 | 46 |
| E | PIESOS | 236 | 131 |
| F | PIETOOLS_demos (without cuadmm, lowrank_2d_stability), PIETOOLS_examples | 57 | 20 |

Deferred: 84 files with uncommitted edits of other sessions (`lists/deferred_dirty.txt`, 79 of them in
PIETOOLS_examples) – passed once their owners commit.

## Execution

Six agents, one per group, in parallel; comment-only edits in the main checkout; each file checked with
`codeident.sh` (code with comments removed identical to HEAD); a ledger per group with candidates for the
maintainer (functional issues are reported, not fixed). Then: a MATLAB syntax check (checkcode) over every
edited file, codeident over the whole edited set, and one commit per group on approval.
