# G3_demos ledger – comment pass of 11/21/2026 (18 files under PIETOOLS_demos/)

Tool note: `cpass.pl` does not recognize markers and header entries of the form `% DJ (CC), 10/30/2026`
(the initials pattern stops at " (CC)"); 58 such lines are in this group. G3 ran prep/apply/check with a
private copy, `$T/G3/cpass_g3.pl`, whose initials pattern also accepts an optional ` (CC)` (and whose
diff scratch files are under `$T/G3/`); `check` was run with both the original and the private copy,
and `$T/G3/xcheck.pl` lists any remaining dated comment other than the TAG. Counted with `(CC)`, no file
of the group has more than half of its new entries by DJ or SS (DEMO3 and DEMO7 have exactly half), so
the classes of G3_demos.tsv stand.

Each edited file: check PASS (original and private copy), xcheck clean, codeident OK, `w/crlf`
unchanged, numstat equal to numstat -w.

## Summary
- condense: 16 edited; asnew: 1 edited (Ch3); deferred: 1 not touched (DEMO2); review: 0.
- Final verification (all 17): codeident 17 OK; check PASS 16/16 (both tool copies); eol unchanged;
  no whitespace-only change; no trailing blanks added. xcheck leaves only Copyright years, the
  `Initial coding` line of GettingStarted and the reworded 2025 entry `DB, 12/29/2024` of DEMO5/6.
- Parse check (MATLAB lock runner, syntax_check): 17 files checked, 0 with parse problems
  (`$T/G3/syntax_out.txt`).
- Open for the maintainer: (1) DEMO5 `popts.sep = 0` (2025: sep = 1), see the DEMO5 entry; (2) several
  code changes since 2025 carry no per-line marker (PIESIM inputs/outputs of 01–02/2026 in DEMO1, 5, 6,
  7, 9 and Ch2, Ch6, Ch11; `y==subs(X,s,1)` in Ch4; the x/phi order in DEMO1 and Ch2); the pass only
  retags existing markers, so these lines stay unmarked; (3) cpass.pl misses `DJ (CC)` markers and
  markers without a preceding blank (tool note above).
- Candidate (code, not changed): Ch6 `ax.FontSize = 14` without semicolon.

## PIETOOLS_demos/DEMO1_stability_and_simulation.m
- Class condense. TAG `YP, DJ, MMP, 11/21/2026`.
- Note:
  ```
  % YP, DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: x is declared before
  %   phi; PIESIM takes the initial values of x and phi as one vector
  %   uinput.ic and returns them in solution.timedep.primary{1} and
  %   primary{2}, the output in regulated{1}; the LPI is declared by
  %   lpiprogram(pie.vars,pie.dom), P = poscopvar(prog,T,2) +
  %   1e-4*eye_copvar_sop(P) and lpi_ineq_sop, in place of poslpivar and
  %   lpi_ineq; the predefined test is lpiscript_sop; help text revised.
  ```
- Header entries removed: 7 (YP 02/17, DJ (CC) 10/29 and 10/30, MMP 11/04, 11/12, 11/14, 11/19).
- apply: dropped 4 (was), kept 0 (was), retagged 5 (2 of them `DJ (CC)`).
- Internal notes removed: the maintainer quotes and "Same program." of the removed entries.
- Open: none.

## PIETOOLS_demos/DEMO3_poincare_inequality.m
- Class condense. TAG `DJ, MMP, AT, 11/21/2026` (AT: the typo fix "Demonstration" in the help).
- Note:
  ```
  % DJ, MMP, AT, 11/21/2026: Changes since PIETOOLS 2025: H2op is the solution
  %   operator of the BVP x_{ss} = v, x(a) = x(b) = 0, returned by convert, and
  %   H1op = diff(H2op,s), in place of the operators T and C1 of the PIE of an
  %   auxiliary heat equation; the LPI is declared on the variables and domain
  %   of H2op and imposed by lpi_ineq_sop, which allows the inequality to fail
  %   outside [a,b] by default; the solution is read by lpigetsol_sop; help
  %   text revised.
  ```
- Header entries removed: 6 (DJ (CC) 10/30, MMP 11/04, DJ 11/06, DJ 11/11, MMP 11/19, MMP 11/20).
- apply: dropped 2 (was) (the removed `st = settings_sopvar_heavy` and its comment), kept 0,
  retagged 2 (`DJ (CC)`).
- By hand: removed 4 lines of unmarked commented-out code added on 11/18/2026 (`%pvar t s s2`,
  `%x = pde_var(...)`, `%BVP = [...]` over two lines) and the blank line after them, and the second
  blank line left by the dropped (was) lines; step comment `x = Top * v` -> `x = H2op * v` (`Top`, the
  output of the retired BVP solver call, no longer exists).
- Unmarked code change: the leading blank of the ` x = pde_var(...)` / ` BVP = [...]` lines (11/18)
  stays (code).
- Open: none.

## PIETOOLS_demos/DEMO4_stability_parameter_bisection.m
- Class condense. TAG `DJ, MMP, 11/21/2026`.
- Note:
  ```
  % DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: the LPI of each bisection
  %   step is declared by lpiprogram(PIE.vars,PIE.dom), P =
  %   poscopvar(prog,T,2) + 1e-4*eye_copvar_sop(P) and lpi_ineq_sop, in place
  %   of poslpivar and lpi_ineq; the predefined executive is lpiscript_sop;
  %   help text revised.
  ```
- Header entries removed: 6 (DJ (CC) 10/29, 10/30, MMP 11/04, 11/12, 11/14, 11/19).
- apply: dropped 4 (was), kept 0, retagged 5 (2 `DJ (CC)`).
- Internal notes removed: maintainer quotes, "Same programs." of the removed entries.
- Open: none.

## PIETOOLS_demos/DEMO5_Hinf_optimal_estimator.m
- Class condense. TAG `YP, DJ, MMP, 11/21/2026` (prep also proposed DB, from the 2025 entry
  `DB, 12/29/2024` whose wording changed by a final `;`; that entry stays, DB made no change since 2025).
- Note:
  ```
  % YP, DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: the LPI is declared
  %   by lpiprogram(PIE.vars,PIE.dom), P = poscopvar(prog,T,4,popts) with
  %   independent lower and upper kernels (popts.sep = 0), Z =
  %   lpivar_cdopvar(prog,C2',4) and lpi_ineq_sop, in place of poslpivar
  %   (sep = 1), lpivar and lpi_ineq; lpigetsol_sop reads the solution and
  %   getObserver_sop computes L by quadratic programming; the estimator is
  %   a pie_struct; the predefined alternative is lpiscript_sop; PIESIM takes
  %   uinput.ic as one vector, no ndiff, and returns solution.timedep.primary
  %   and regulated as cells; help text revised.
  ```
- Header entries removed: 7 (YP 02/17, DJ (CC) 10/29, 10/30, MMP 11/04, 11/19 twice, 11/20).
- apply: dropped 5 (was), kept 0, retagged 13 (11 `DJ (CC)`).
- By hand: step comment of `popts.sep = 0;` changed from "set the lower and upper kernels equal" to
  "lower and upper kernels independent" (poscopvar: sep true replaces the lower and upper basis
  operators by the full integral, i.e. equal kernels; the code sets sep = 0).
- Open (maintainer): `popts.sep` went from 1 (2025 `opts.sep = 1`, and DJ 10/30) to 0 in bb367fe1
  (11/18, settings commit) without a marker or note; the note and comment now describe sep = 0 – confirm
  that sep = 0 is intended.

## PIETOOLS_demos/DEMO6_Hinf_optimal_control.m
- Class condense. TAG `YP, DJ, MMP, 11/21/2026` (prep's DB: the reworded 2025 entry, as in DEMO5).
- Note:
  ```
  % YP, DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: the LPI is declared
  %   by lpiprogram(PIE.vars,PIE.dom), P = poscopvar(prog,T',4) +
  %   1e-3*eye_copvar_sop(P), Z = lpivar_cdopvar(prog,B2',2) and
  %   lpi_ineq_sop, in place of poslpivar, lpivar and lpi_ineq; lpigetsol_sop
  %   reads the solution and getController_sop computes K; the closed loop is
  %   a pie_struct; the predefined alternative is lpiscript_sop; PIESIM takes
  %   uinput.ic as an array, no ndiff, and returns solution.timedep.primary
  %   and regulated as cells; help text revised.
  ```
- Header entries removed: 8 (YP 02/17, DJ (CC) 10/29, 10/30, MMP 11/04, 11/12, 11/14 twice, 11/19).
- apply: dropped 7 (was), kept 0, retagged 14 (10 `DJ (CC)`).
- Internal notes removed: maintainer quotes, "Same program." of the removed entries.
- Open: none.

## PIETOOLS_demos/DEMO7_observer_based_control.m
- Class condense. TAG `YP, DJ, MMP, 11/21/2026` (YP has no entry; the PIESIM changes of 01/31 and
  02/18/2026 are unmarked).
- Note:
  ```
  % YP, DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: K and L are computed
  %   by lpiscript_sop with the settings lpisettings_sop('heavy') and
  %   converted to copvar objects by op2copvar_sop; the estimator output zhat
  %   is the first output of C1, selected by a copvar matrix; PIESIM takes
  %   uinput.ic as an array, no ndiff, and returns solution.timedep.primary
  %   and regulated as cells; help text revised.
  ```
- Header entries removed: 4 (DJ (CC) 10/29, 10/30, MMP 11/08, 11/19).
- apply: dropped 0 (was), kept 4 (was) (2025 code: `lpisettings('heavy')`, two `lpiscript`, `piess`
  with `C1(1,:)`), retagged 6 (5 `DJ (CC)`).
- Open: none.

## PIETOOLS_demos/DEMO8_H2norm.m
- Class condense. TAG `DJ, MMP, 11/21/2026`.
- Note:
  ```
  % DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: the LPI is declared by
  %   lpiprogram(PIE.vars,PIE.dom), W = poscopvar(prog,T',3) and lpi_ineq_sop,
  %   with the trace of C1*W*C1' by trace_rn_sop, in place of poslpivar,
  %   lpi_ineq and trace(Aux.P); lpigetsol_sop reads the solution; help text
  %   revised.
  ```
- Header entries removed: 5 (DJ (CC) 10/30, MMP 11/04, 11/12, 11/14, 11/19).
- apply: dropped 3 (was), kept 0, retagged 7 (5 `DJ (CC)`).
- Internal notes removed: maintainer quotes, "Same program." of the removed entries.
- Open: none. (The 2025 comment "W>=0 with default polynomial degrees up to 3" stays although the
  degree 3 is now passed explicitly.)

## PIETOOLS_demos/DEMO9_L2_gain_analysis.m
- Class condense. TAG `YP, DJ, MMP, 11/21/2026` (DJ also: psatz 2 -> 3 on 08/19/2026, superseded).
- Note:
  ```
  % YP, DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: PIESIM takes
  %   uinput.ic = 0 and returns the 2D state in solution.timedep.primary{4},
  %   the output in regulated{1}; the LPI is declared by
  %   lpiprogram(PIE.vars,PIE.dom), P = poscopvar(prog,T,1) and lpi_ineq_sop
  %   with its default options, in place of poslpivar and lpi_ineq with
  %   opts.psatz = 2; lpigetsol_sop reads the solution; the predefined
  %   executive is lpiscript_sop; help text revised.
  ```
- Header entries removed: 7 (YP 02/17, DJ (CC) 10/29, 10/30, MMP 11/04, 11/12, 11/14, 11/19).
- apply: dropped 3 (was), kept 0, retagged 5 (3 `DJ (CC)`).
- Internal notes removed: maintainer quotes, "Same program." of the removed entries.
- Open: none.

## PIETOOLS_demos/snippets_from_manual/GettingStarted.m
- Class condense. TAG `DJ, MMP, 11/21/2026`. The 2025 file had no log; the note follows the
  `Initial coding` line (added on 11/20 with the license).
- Note:
  ```
  % DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: the stability test runs
  %   lpiscript_sop; help text, license and authors added.
  ```
- Header entries removed: 1 (MMP 11/19 release comment pass).
- apply: dropped 0, kept 1 (was) (2025 `lpiscript` call), retagged 1 (`DJ (CC)`).
- Open: none.

## PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch11_Demos.m
- Class condense. TAG `YP, MMP, 11/21/2026`.
- Note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: the grids returned by
  %   PIESIM (grid, grids) are plotted directly, the physical grid being the
  %   output itself rather than its field phys; help text revised.
  ```
- Header entries removed: 1 (MMP 11/19). apply: 0 / 0 / 0 (no markers; YP's 10 changed lines
  `grid.phys` -> `grid` are unmarked).
- Open: none.

## PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch2_Introduction.m
- Class condense. TAG `YP, DJ, DB, MMP, 11/21/2026` (DB: `grid` -> `grids`, 11/06/2026).
- Note:
  ```
  % YP, DJ, DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: x is declared
  %   before phi; PIESIM takes the initial values of x, phi1, phi2 as one
  %   vector uinput.ic, no ndiff, and returns the physical grid itself and the
  %   states and outputs in solution.timedep.primary{2} and regulated{1};
  %   lpiscript_sop runs the executives, with lpisettings_sop('heavy') for
  %   stability and the L2-gain; help text revised.
  ```
- Header entries removed: 4 (YP 02/23, DJ (CC) 10/29, MMP 11/08, 11/19).
- apply: dropped 0, kept 4 (was) (2025 `lpisettings('heavy')` and three `lpiscript` calls), retagged 4.
- Open: none.

## PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch3_PI_Operators.m
- Class asnew (23% of the 2025 code left). Edited by hand (Edit tool), no apply.
- Removed: log entries DJ 11/08/2026 (2 lines) and MMP 11/19/2026; 28 (was) lines (the 2025 opvar /
  opvar2d declarations of A, B, C, E, F); 18 `% DB, 11/06/2026` markers (3 on comment lines appended
  without a blank, 1 after the comment "dummy variable s_dum").
- Kept: license, help header, support line and `% DJ, 12/30/2024: Initial coding;` (the only dated
  line left), the step comments on pi4_copvar_sop arguments and dummy-variable names.
- codeident OK, `w/crlf`, numstat equal, no trailing blanks added.
- Open: none.

## PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch4_Declaring_PDEs.m
- Class condense. TAG `DJ, MMP, 11/21/2026`.
- Note:
  ```
  % DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: the observed output of
  %   4.1.1 is y = X(1), in place of X(0), which the boundary condition X(0) = 0
  %   made identically zero; help text and step comments added.
  ```
- Header entries removed: 2 (DJ 11/08, MMP 11/19); `DJ, 12/11/2025` is a 2025 entry and stays.
- apply: 0 / 0 / 0 (the changed line `y==subs(X,s,1)];` is unmarked).
- Open: none.

## PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch5_PIEs.m
- Class condense. TAG `DB, MMP, 11/21/2026`.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the operators of the PIEs
  %   of 5.5 are 'copvar' objects built by pi4_copvar_sop on s, s_dum and
  %   dom = [0,1], the constant and polynomial blocks included (B1 of PIE2;
  %   A, Tu, B1, D12, D22 of PIE3), in place of opvar objects, scalars and
  %   polynomials; help text revised.
  ```
- Header entries removed: 1 (MMP 11/19).
- apply: dropped 0, kept 25 (was) (the 2025 opvar declarations), retagged 18.
- By hand: 2 markers `% DB, 11/06/2026` appended without a blank (a comment line on pi4_copvar_sop,
  and the `A = pi4_copvar_sop(...)` line of 5.5.4) retagged with a blank before them; neither cpass
  prep nor check sees a marker without a preceding blank (xcheck did).
- Open: none.

## PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch6_PIESIM.m
- Class condense. TAG `YP, VJ, DJ, MMP, 11/21/2026`.
- Note:
  ```
  % YP, VJ, DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: Section 6.4.2
  %   added (radial heat equation, 1D PIESIM library example 39), later
  %   sections renumbered; 6.4.3 uses 2D library example 10; PIESIM takes
  %   opts.ifexact, initial values as arrays uinput.ic, boundary inputs as
  %   cells uinput.w{k}, no ndiff, and returns the physical grid itself
  %   (grid{1}, grid{2} in 2D), the states in solution.timedep.primary and
  %   final.primary, the outputs in regulated{1}; the DDE and its closed loop
  %   start from x(0) = [1,1]; lpiscript_sop computes the controller; help
  %   text revised.
  ```
- Header entries removed: 4 (VJ 02/14, YP 02/17, DJ (CC) 10/29, MMP 11/19).
- apply: dropped 0, kept 1 (was) (2025 `lpiscript` call), retagged 1.
- Candidate (code, not changed): `ax.FontSize = 14` in 6.4.5 (added since 2025) has no semicolon and
  prints the axes object.
- Open: none.

## PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch7_LPI_Programming.m
- Class condense. TAG `DB, MMP, 11/21/2026`.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the LPI is declared with
  %   'copvar' operators: the spaces of P, Z and the identities are read from
  %   T, C2, B1, C1 by copvar_space_list and declared on dom by poscopvar
  %   (multiplier term only), lpivar_cdopvar and eye_copvar_sop; D21 is put on
  %   the registry of Z; lpi_ineq_sop with like = P, dD = 2 and lpi_eq_sop
  %   replace lpi_ineq (psatz = 1) and lpi_eq; getsol_lpivar_sop,
  %   lpigetsol_sop and getObserver_sop read the solution and the gain; help
  %   text revised.
  ```
- Header entries removed: 1 (MMP 11/19).
- apply: dropped 0, kept 31 (was) (2025 code), retagged 49.
- By hand: removed 4 commented-out lines with non-(was) markers `% DB, 11/04/2026` (`% PIEo =
  pie2opvar(PIE);` and the three `% T = PIEo.T; ...` lines): a development step of 11/04 replaced on
  11/06, never in 2025 (apply would have retagged them).
- Note: the kept (was) lines `% T = PIE.T; ...` (2025) repeat the current code lines `T = PIE.T; ...`,
  which DB re-marked on 11/06 (the PIE fields are now copvar objects); kept as the rule says.
- Open: none.

## PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch8_Declaring_PDEs.m
- Class condense (comments only since 2025). TAG `MMP, 11/21/2026`.
- Note: `% MMP, 11/21/2026: Changes since PIETOOLS 2025: help and comments revised.`
- Header entries removed: 1 (MMP 11/19). apply: 0 / 0 / 0.
- Open: none.

## PIETOOLS_demos/DEMO2_volterra_operator_norm.m
- Class deferred: uncommitted edits of another session (` M` at the start of this pass). Not touched.
