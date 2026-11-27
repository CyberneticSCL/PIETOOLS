# Phase 1 – tools and dry run (11/26/2026)

Maintainer, 11/26/2026: "proceed to phase 1". The rename was run on an export outside OneDrive
(`C:/wt1126ren/base`: 8198c35d plus the uncommitted H2 pass and the C8 move; renamed copy
`C:/wt1126ren/new`). The checkout is not renamed.

## Tools

- `ren_apply.pl ROOT paths`: the moves, `old<TAB>new`.
- `ren_apply.pl ROOT check|apply DATE existing35.txt REPORTDIR`: the map of PLAN.md, Phase 1
  step 2, on every live .m file at its path after the moves.
  - Token map: the 40 rows of `rename_table.tsv` that are renamed, in lower, upper and
    capitalized form. A token naming a container that is neither mapped nor retired stops the
    run (UNKNOWN).
  - Records kept as written: header-log lines, lines holding `(was)`, and the nine retired names
    of PLAN.md, Phase 1.
  - Article: "a copvar" becomes "an opvar", also when "a" ends the comment line before.
  - A trailing comment after 2 or more blanks keeps its column, so markers stay at column 77.
  - D1: one log line, after the header log, in each file of `existing35.txt` that changes.
- `ren_check.sh BASE NEW paths.tsv DATE OUT` with `ren_cmp.pl`: checks 0 and 3a–3c, written apart
  from `ren_apply.pl`. Their map is a sed script (`map.sed`) generated from `rename_table.tsv`.
- `ren_mcheck.m`: checks 3d and 3e. `ren_ccmp.m`: `checkcode` messages of each changed file,
  before and after.

## Dry run (`dryrun_1126/`)

- Moves: 97 files – `@copvar` 33, `@cdopvar` 30, 28 toolbox function files, 6 test files. The
  98 of Section 2 lose `op2copvar_sop` and `pie2copvar` (D4) and gain `poslpivar_sop` →
  `poslpivar_legacy` (C7). No target existed.
- Map: 265 live .m files, 2,213 lines; 2,580 tokens of 79 kinds (`tokens.tsv`); UNKNOWN 0.
- Left as written: 21 retired-name tokens; 138 record lines that name a container (107
  header-log, 31 `(was)`; `records.tsv`).
- Article: 140 insertions, 15 of them at a line end. 18 lines with an insertion are code, all
  inside error, warning or assert strings.
- Realigned: 93 trailing comments, every one at its old column.
- D1 lines: 32 files. DEMO2, `@pie_struct/horzcat` and `@pie_struct/vertcat` name a container
  only in header logs, so they do not change.
- A second run changes nothing.

Checks:
- 0: every file outside the change list is byte-identical; each moved file exists.
- 3b: code identity in 265 of 265 files.
- 3c: comment identity, 0 defects.
- 3a: no container name outside the records and the retired names.
- 3d, on the factory path plus `genpath` of each tree:
  - each of the 29 new names resolves once, inside the new tree; the 29 old names resolve
    nowhere;
  - the 7 private files exist and resolve only to themselves;
  - `opvar()` and `dopvar()` construct, and `pde_struct` loads (its `InferiorClasses` names
    `?opvar`);
  - the method lists equal those of `copvar` and `cdopvar` up to the constructor (34 and 29).
- 3e: no parse or syntax error in the 265 files, and their `checkcode` messages are identical
  before and after.

## Applying it to the checkout

After the Phase 0 commits, with every MATLAB session closed and OneDrive sync paused (H7):
1. `ren_apply.pl PIETOOLS paths` must give `dryrun_1126/paths.tsv`. Then `git mv` each pair,
   file by file.
2. `ren_apply.pl PIETOOLS apply MM/DD/YYYY existing35.txt REPORTDIR`.
3. Every live file must equal its copy in `C:/wt1126ren/new` byte for byte; only the date of the
   32 D1 lines may differ. Then `ren_mcheck.m` on the checkout, and the commit on approval.
