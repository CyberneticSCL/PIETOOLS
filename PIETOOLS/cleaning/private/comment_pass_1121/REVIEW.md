# Comment pass over files of the 2025 release, 11/21/2026 – outcome and review items

Maintainer request and clarifications: `BRIEF.md`. Records: `census_2025_files.tsv` (275 files of the 2025
release outside PIESOS/, SOSTOOLS400/, cleaning/, testing_codes/, orphaned_file/ and the two non-release demo
folders), group lists `G1`-`G5` (`.tsv`), one ledger per group (`<group>_ledger.md`, every note verbatim),
and the tools (`cpass.pl`, `review.sh`, `kept2025.pl`, `wasaudit.pl`, `restore_was.pl`; scratch paths).

## Outcome

| Class | Files | Result |
|---|---|---|
| condense | 169 | one note `% <initials>, 11/21/2026: Changes since PIETOOLS 2025: ...` (at most 9 lines) in place of the notes added since 2025; markers of lines changed since 2025 retagged to the same initials and date; (was) lines of 2025 code kept; commented-out code of development steps removed; internal notes on development and testing removed |
| asnew (under 25% of the 2025 code left) | 3 | treated as new files, no developer notes: `PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch3_PI_Operators.m`, `command_line_input_Examples/stability_13.m`, `executives_sopvar/utility_functions/closedLoopPIE.m` |
| review (new notes mainly by DJ) | 2 | not edited; draft note and verdicts in `G2_conv_root_ledger.md`: `converters/convert_PIETOOLS_PDE.m` (8 of 11 new entries by DJ; also only 12% of its 2025 code left, so it would be asnew), `converters/legacy/2D/Compute_Tmap_PDE_2D_Separable.m` (3 of 4) |
| deferred | 1 | `PIETOOLS_demos/DEMO2_volterra_operator_norm.m`: another session holds uncommitted edits since 11/18 |
| unchanged since 2025 | 3 | `PIESIM/` files identical to their 2025 version; nothing to do |

parser/ and PIESOS/ were left alone, as asked.

Checks over the 172 edited files: code with comments removed identical to HEAD (codeident); `cpass.pl check`
PASS (the single note, every new marker carries its tag, every (was) line is 2025 code or 2025 commented
code, markers at column 77); no line shared by HEAD and the 2025 version lost; line endings unchanged
(w/crlf, and the w/lf or w/mixed files as they were); no whitespace-only change; MATLAB parse check 172
files, 0 problems.

## Edits of the coordinator, after the agents

- `cpass.pl` first treated only live 2025 code as 2025 code, so apply dropped (was) lines whose text was
  commented-out code already in 2025. Audited every file (`wasaudit.pl`); 13 such lines put back, retagged:
  12 commented template calls of Step 3 in `PIETOOLS_PDE.m`, 1 commented `eq_ODEPDE` line in
  `Stability_Examples/PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m`. Left out (old help prose, corrected in place
  as help text is): 3 lines in `2D/PIETOOLS_PDDE_Ex.m`, `Hinf_Gain_Examples/
  PIETOOLS_PDE_Ex_Heat_Eq_with_Distributed_Disturbance.m` ("gamma = 0.3333") and the Das file.
- Comment words damaged by the rename no -> n0 after 2025 restored in 5 PIESIM files (not, now, denotes,
  polynomial, note, enough).
- Rules decided in the pass: the initials of a note are those of the authors of changes to the content since
  2025 (a move without change, or a revert of a merge, does not count); lines changed since 2025 without a
  marker get none (only existing markers are retagged); revised help text keeps no (was) lines.

## Questions for the maintainer

1. `converters/convert_PIETOOLS_PDE.m`, `converters/legacy/2D/Compute_Tmap_PDE_2D_Separable.m`: approve the
   drafted notes, or (first file) treat as new.
2. `PIETOOLS_PDE.m`: the 12 restored (was) lines call executives that left the release (stock
   PIETOOLS_*). They follow the rule (2025 commented code stays); drop them?
3. `examples_DDE_library_PIETOOLS.m`, `examples_NDSDDF_library_PIETOOLS.m`: 7 (was) lines of 2025 commented
   example code kept by the same rule.
4. `PIETOOLS_demos/DEMO5_Hinf_optimal_estimator.m`: `popts.sep = 0` since 11/18 (bb367fe1, no note); the
   2025 demo used sep = 1. The note and the step comment describe sep = 0. Intended?
5. Lines changed since 2025 without a marker (code of YP, DB, DJ in PIESIM, demos and examples; the notes name
   the changes): add the tag marker to them?
6. In 11 converter files whose only commit by DJ is the move of 11/16 (a222bbc7, content unchanged) the tag
   is `MMP, 11/21/2026` (list in the G2 ledger).

## Code problems seen in the pass (comment-only pass: not changed)

- `PIESIM/PIESIM_routines/Postprocessing/PIESIM_transform_to_solution.m:155`: the message reads "does n0t
  match" (rename damage in a string).
- `examples_PDE_library_PIETOOLS.m`: the error message of case 16 names a date ("removed ... on 11/19/2026").
- `PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch6_PIESIM.m`: `ax.FontSize = 14` without
  a semicolon prints the axes object.
- `Stability_Examples/PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi.m`: `if logvals` (2025 typo) errors for
  GUI = 1.
- `command_line_input_Examples/stability_3.m`: Lm is built with the constant 10 instead of g, so a params
  override of g does not reach Lm.
- `PIESIM/.../PIESIM_Poly2Mat_cheb_2D.m`: the padded coefficient array grows past (N_x+1) x (N_y+1) when the
  degree plus 2 exceeds N+1 (by reading, not run).
- Stale text kept: comments on rsize/csize in `PIESIM_9PI2Mat_cheb_2D.m`; the PIESIM help lists
  ic.pde/ic.pie/ic.ode while the code reads one uinput.ic.
