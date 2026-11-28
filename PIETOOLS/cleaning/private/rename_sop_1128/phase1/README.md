# Phases 0 and 1 – review, tools, dry run, checkout (11/28/2026)

Maintainer, 11/28/2026: "Lets do this now"; decisions E1–E4 in `../PLAN.md`, Section 5. Phase 1
committed as fe6d3530.

## Phase 0 – code that names a 2025 name of the table

Live code that already used a new name, which resolves nowhere today and binds to the renamed
function after Phase 1 (`oldcode.pl` of the inventory session; the variable `settings` aside):
30 lines, all in dated baselines under `testing_codes/private/baselines/` (`2026-10-28_example_library`,
`2026-11-02_opvar_phase0`, `2026-11-05_class_ops`, `2026-11-17_tiers`, `2026-11-18_studies`):
`lpisettings(...)` 21, `lpi_ineq` 4, `PIETOOLS_Hinf_gain` 3 (labelled "stock"), `lpi_eq`,
`getsol_lpivar`, `settings_PIETOOLS_heavy` 1 each. They called the stock functions, off the path
since 11/08. Kept as written, as C6b and C6d of the container rename (dated records); a note in
`testing_codes/README.md` (Phase 2) states that these names now resolve to the live functions.
In the same folders `lb_run.m` and `ab_bisect.m` test `contains(r.called,'_sop')`, also kept.

Names built at run time: `test_executives_sop.m` builds `['PIETOOLS_' nm{1} '_sop']` (lines 115,
116) – a rule of `sop_apply.pl` drops the suffix in code; the commented-out copies stay. The text
`PIETOOLS_*_sop` (5 comment lines) becomes `PIETOOLS_*`. The labels the executives print
(`'H2_control_sop'`, `['stability_sop ' form]`, ...) lose the suffix.

Not changed, pre-existing: `See also ... GETCONTROLLER_DIRECT_SOP` in `sopvar/@sopvar/inv.m`,
`sopvar/misc/chebfit_mono.m`, `sopvar/misc/gk_grid.m` names a function retired by a3fec6b0.

## Tools

- `sop_apply.pl ROOT TABLE paths` and `... check|apply DATE EXISTING REPORTDIR`: `ren_apply.pl`
  of the container rename without the article rule; the map from `../rename_table.tsv` (renamed
  rows and folders), the nine labels and two abbreviations of `bl_sop.m`; a token holding `sopvar`
  names the class and stays; any other token with `_sop` must be listed in `@KEEP`, else UNKNOWN.
- `sop_check.sh`, `sop_cmp.pl`: checks 0 and 3a–3c with a sed map generated apart (`map.sed`).
- `sop_mcheck.m`: 3d (each moved function resolves once, at its new path; no old name
  resolves; off the path a file resolves nowhere but to itself) and 3e (no mtree or syntax
  error). `ren_ccmp.m` of the container rename: `checkcode` before and after.
- `sop_cmp_tree.pl`: the checkout against an export, outside `cleaning/` and `orphaned_file/`.

## Dry run (`dryrun_1128/`, export of 91ea628e at `C:/wt1128sop/base`, renamed in `.../new`)

- Moves: 92 files – `executives_sopvar` 29, `lpi_programming_sopvar` 33, `settings_sopvar` 9,
  the test mirrors 1 and 19, `sopvar/misc/on_registry_sop.m`; no target existed; the five
  folders emptied and removed.
- Map: 266 live .m files, 1,718 lines; 1,850 tokens of 129 kinds; kept 196, class 441; realigned
  trailing comments 215; record lines left as written 211; D1 lines 72 – the 74 files of the
  2025 release that change, less two baseline copies of `lpisolve.m` (`gate/`) without a header
  log; UNKNOWN 0. A second run changes nothing.

Checks:
- 0: no unlisted file changed, none missing.
- 3b: code identity in 266 of 266 files.
- 3c: comment identity, 0 defects; 72 D1 blocks; 134 record lines kept.
- 3a: no old name outside the records (135 record lines keep one).
- 3d: 92 moves (64 on the path, 28 off it), 0 defects. 3e: 266 files, 0 errors.
- `checkcode`: identical but for IFBDUP in the four dated copies of `lb_run.m` (line 149 and
  near): its `'sop'` and stock branches now both call `lpisettings` – the Phase 0 finding above.

## Applied to the checkout, 11/28/2026

No MATLAB session running (the desktop of 11/21 had been closed; `PIETOOLS_H2_norm_sop.m` equal
to HEAD, so its buffer was not saved); OneDrive running. Before: the checkout equals
`C:/wt1128sop/base` in all 2,578 files outside `cleaning/` and `orphaned_file/` (1,364 git-ignored
files aside). `sop_apply.pl paths` gives `dryrun_1128/paths.tsv`; the 92 files moved with `mv`, the
five folders removed with `rmdir`; `apply 11/28/2026` (`apply_1128/`): every report equals its
dry-run copy; a second run changes nothing. After: the checkout equals `C:/wt1128sop/new` in all
2,578 files; `git status` lists 92 deleted, 92 new and 196 modified files, exactly the expected
sets; no new path is git-ignored. Committed fe6d3530: 288 paths (92 renames, 196 modified),
+1,862 −1,718.
