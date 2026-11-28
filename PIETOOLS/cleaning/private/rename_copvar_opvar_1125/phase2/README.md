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

## Follow-up, 11/27/2026: the stock path in the agent file

Maintainer, on the four mentions left above: "Yes, update it". The file named code retired on
11/08 (`executives/`, `lpi_programming/`, `settings/`, now under `cleaning/private/`) in more
places than four; each now names what is on the path:

- description: "(stock or container path)" -> "of executives_sopvar"; `lpi_programming*` ->
  `lpi_programming_sopvar`.
- Step 5: the stock executives and `lpiscript` dropped; `lpiscript_sop`, whose help maps each name
  to an executive and form; `executives_sopvar/README.md` (moved to `orphaned_file/claude` by
  010e8bd0) replaced by the help of each executive; the stability dispatch stated by the forms of
  `PIETOOLS_stability_sop` ('stability' of `lpiscript_sop` = 'Q'; default form 'P'); settings by
  the tiers of `lpisettings_sop`, with the rule of 11/19 that an example names its settings and
  executive form; the dual-form note measured again (below).
- Step 7: `getController`, `getObserver` dropped; the methods 'qp' and 'inv' of the `_sop`
  versions, as their help states.
- Running MATLAB: `which -all lpiscript` -> `which -all opvar` (a PIETOOLS 2025 installation
  shadows the container class) and `which -all lpiscript_sop`; `convert` lists two files, the
  methods of `pde_struct` and of the legacy `sys`, so "resolve once" was wrong for it.
- "Size b mismatch": the comparison with the stock routine dropped.

Not changed: `claude_tests/` in "Writing scripts" (no such folder on the path).

Measurement (`dualw_check.m`, `dualw_1127.tsv`): 1-D reaction-diffusion, lam = pi^2/2,
Dirichlet, u(t) entering uniformly, z = [int x; u], tier light; w(s,t) entering the PDE directly,
or w(t) through s(1-s).

| run | w(s,t) | w(t) |
|---|---|---|
| `PIETOOLS_Hinf_gain_sop` 'Q' | optimal 0.18248168 | optimal 0.033310667 |
| 'P' | optimal 0.18248066 | optimal 0.03331102 |
| 'Qd' | optimal 0.18248329 | optimal 0.033310846 |
| 'Pd' | infeasible | optimal 0.033313368 |
| `PIETOOLS_H2_norm_sop` 'c' | optimal 0.28774368 | optimal 0.052428229 |
| 'o' | unknown 1.2e-3 (floor) | optimal 0.052428372 |
| 'cco' | infeasible | optimal 0.052396853 |
| 'oco' | optimal 7.3e-5 (floor) | optimal 0.052393641 |
| `PIETOOLS_Hinf_control_sop` | infeasible | optimal 0.032866021 |
| `PIETOOLS_H2_control_sop` | infeasible | optimal 0.052219933 |

The w(t) column equals the 1-D table of the moved `executives_sopvar/README.md` (10/2026, column
'like') to 5 or more digits, which identifies its names: `Hinf_gain_dual` = 'Qd',
`Hinf_gain_dual_coercive` = 'Pd', `H2_norm_c_coercive` = 'cco', `H2_norm_o_coercive` = 'oco'. Its
note listed `Hinf_gain_dual` among the dual forms that fail for a distributed w, with no 1-D row
for it (the evidence was the 2-D row, INFEASIBLE at `max_raise` 0); in 1-D, 'Qd' certifies. 2-D
not measured again.
