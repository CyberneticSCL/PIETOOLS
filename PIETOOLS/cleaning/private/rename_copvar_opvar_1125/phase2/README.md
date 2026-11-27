# Phase 2 – documentation, 11/27/2026

Autonomous run of 11/27/2026 (maintainer: "run autonomously for the next 12 hours to finish up
this copvar/opvar switch"; for the manual drafts, "Skip the manual").

## Changed

- `sopvar/sopvar_lift_notes.tex` (31 lines) and `sopvar/sopvar_duality_pruning_notes.tex` (5
  lines): the map of Phase 1, by `tex_map.pl` (names written with `\_` and with `_`; the article
  rule; the retired `test_poscopvar_lift` and `test_poscopvar_direct` kept). Five mentions of the
  2025 functions, ambiguous after the rename, now name the release: in the lift notes, the 1-D
  legacy `poslpivar` (introduction), the `lpi_ineq` pair (degree count in 1-D), `d{2}(3)` of
  `poslpivar` and the first column of Table `tab:slots`, `@opvar/degbalance`. Kept: mentions true
  of the container functions (the `lpivar` family of Q; the identity basis of `poslpivar`).
- Their PDFs, rebuilt with MiKTeX pdflatex in a scratch folder (two passes, no error, no undefined
  reference). The rebuild also carries the edits of 010e8bd0 (11/22; 2 and 1 lines), which the
  PDFs of 10/08 and 11/04 lacked.
- `.claude/agents/pietools-applications.md`: the user-level LPI functions are those of the root of
  `lpi_programming_sopvar`, by their new names; the retired stock functions and `lpiscript` (no
  such file; the function is `lpiscript_sop`) are dropped. Not changed: four other mentions of the
  stock path (description, steps 4, 7 and 9), retired on 11/08 and outside the rename.
- `testing_codes/README.md`, `cleaning/README.md`: a translation note; the dated records keep the
  names of their date.
- Memory: `opvar-rename-map.md`, the map.

Not in git: `CLAUDE.md` of the repository root (untracked). The class table and the dependency
chain had, before:

```
opvar / dopvar  ->  lpi_programming 1D + 2D  ->  executives  ->  examples
... exercise it through `lpi_programming` and an actual executive or example.
| `opvar`, `dopvar` | `opvar/@opvar`, `opvar/@dopvar` |
| `opvar2d`, `dopvar2d` | `opvar/2D/@opvar2d`, `opvar/2D/@dopvar2d` |
| `copvar`, `cdopvar` (containers; renamed from `mopvar`, `mdopvar` 2026-09-25) | `sopvar/@copvar`, `sopvar/@cdopvar` |
```

and now name `lpi_programming_sopvar` and `executives_sopvar`, the containers as `opvar`,
`dopvar` in `sopvar/@opvar`, `sopvar/@dopvar`, and the 2025 classes as retired to
`cleaning/private/opvar/`.

## Kept as written (records)

READMEs of the dated baselines, results and logs of dated runs, the F2 indices,
`opvar_removal_plan.md`, `MIXED_DOMAIN_REPORT.md`, `BENCH.md`, `night_now.sh`, and the tex reports
under `cleaning/private/settings_1116/`.

## Skipped

The manual drafts (maintainer). `PIETOOLS_2025_User_Manual.pdf` describes the 2025 classes (H9).
