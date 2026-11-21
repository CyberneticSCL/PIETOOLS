# G5_examples_b – ledger of the comment pass of 11/21/2026

Group: PIETOOLS_examples/Examples_Library/Stability_Examples and command_line_input_Examples, 40 files
(`G5_examples_b.tsv`). Procedure and tools: `BRIEF.md`. Every file was clean (`git status --porcelain`)
and `i/lf w/crlf` before its edit. Counts "apply" = dropped (was) lines / kept (was) lines / retagged
markers, as printed by `cpass.pl apply`.

## Files

### Stability_Examples/PIETOOLS_PDE_Ex_Diagne.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 0/1/9 (prep: 9 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the BCs read
  %   [x1(1); x2(0); x3(0)] = K [x1(0); x2(1); x3(1)], the incoming traces
  %   in terms of the outgoing ones, in both formats (were [x1(0); x2(1);
  %   x3(1)] = K [x1(1); x2(0); x3(0)]); Lm uses g in place of the constant
  %   10; default gains k1 = -0.5, k2 = -4.5 (were 1, 2); new output st, the
  %   recommended settings of the stability test; stability_form = 'P' set
  %   in the base workspace; help and comments revised.
  ```
- Removed: header entries of 11/19 (settings search, "maintainer: no default settings") and the release-pass entry. Kept: the (was) line of the 2025 signature; the unmarked commented-out 2025 code of DB (k1 = 1, k2 = 2; the 2025 K partition; the 2025 ny1/on0 lines).
- Observed: the 2025 step comments `% BCs: 0 = [1,-K01;0,-K11] * x(0)` no longer show the permutation Ppm in C0_raw, C1_raw (2025 text, left).

### Stability_Examples/PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 2/0/17 (prep: 17 RETAG, 2 DROP-WAS – the (was) signature with INDEX and the (was) error line with indx are intermediate code of DB, not 2025 code); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: new first input INDEX:
  %   13 or 13.1 gives the states u_{t}, u_{ss} (the 2025 system), 13.2 the
  %   state u, 4th order in s, declared with pde_var, any other value an
  %   error; new output st, the recommended settings of the stability tests
  %   by case; stability_form = 'P' (13, 13.1) or 'Q' (13.2) set in the base
  %   workspace; help and comments revised.
  ```
- Removed: the 11/19 settings entry, the release-pass entry, the 11/20 entry on the error message (a fix of intermediate code, not of 2025 code). No (was) line of the 2025 signature `(GUI,params)` exists to keep.

### Stability_Examples/PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m
- class condense; TAG `MMP, 11/21/2026`; apply 4/3/15 (prep: 15 RETAG, 3 KEEP-WAS, 4 DROP-WAS – the dropped (was) lines hold 2025 comment text, not code); check PASS, codeident OK.
- Note:
  ```
  % MMP, 11/21/2026: Changes since PIETOOLS 2025: the ODE takes the input
  %   -Bxr*x_{s}(a) in place of +Bxr*x_{s}(a), in both formats and in the GUI
  %   save file (the sign of the paper is taken as a typo), so that the
  %   system is stable for lam < 9.7946; new output st, the recommended
  %   settings of the stability tests; stability_form = 'Q' set in the base
  %   workspace; help and comments revised.
  ```
- Removed: three 11/19 entries (maintainer quote, Chebyshev node counts, settings search), the release-pass entry, the 11/21 entry on the GUI .mat file (folded into the note; the .mat changed in a8bd220e).
- Kept: the help lines on the stability bound and the sign of the paper (lam < 9.7946; +Bxr stable only for lam < -5.7799) – a property of the system and a literature attribution, not a test record.

### Stability_Examples/PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Tang.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 0/1/9 (prep: 9 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the term-based BC at
  %   s = 0 reads x(0) = -xo, as the batch format (it read x(0) = xo), and
  %   the help states the BCs x(0) = -xo, x(1) = k*xo; new output st, the
  %   recommended settings of the stability test; stability_form = 'Q' set
  %   in the base workspace; help and comments revised.
  ```
- Removed: the 11/19 settings entry and the release-pass entry. Batch convention checked: B*[x(a); x(b); ...] = Bx*xo, so Bx = [-1; k] gives x(0) = -xo, x(1) = k*xo.

### Stability_Examples/PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_in_BC.m
- class condense; TAG `MMP, 11/21/2026`; apply 0/1/7 (prep: 7 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % MMP, 11/21/2026: Changes since PIETOOLS 2025: new output st, the
  %   recommended settings of the stability tests; stability_form = 'P' set
  %   in the base workspace; help and comments revised.
  ```
- Removed: the 11/19 settings entry and the release-pass entry.

### Stability_Examples/PIETOOLS_PDE_Ex_Lamare.m
- class condense; TAG `MMP, 11/21/2026`; apply 0/1/12 (prep: 12 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % MMP, 11/21/2026: Changes since PIETOOLS 2025: new output st, the
  %   recommended settings of the stability test by case; stability_form =
  %   'P' set in the base workspace; help and comments revised.
  ```
- Removed: the 11/19 settings entry and the release-pass entry.

### Stability_Examples/PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 0/1/7 (prep: 7 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: default lam = 4.65
  %   (was 4.6); new output st, the recommended settings of the stability
  %   tests; stability_form = 'Q' set in the base workspace; help and
  %   comments revised.
  ```
- Removed: the 11/19 settings entry and the release-pass entry.

### Stability_Examples/PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 0/1/8 (prep: 8 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the batch format
  %   imposes x(1) = 0 in place of x_{s}(1) = 0, as the term-based format
  %   and the help; new output st, the recommended settings of the stability
  %   tests; stability_form = 'P' set in the base workspace; help and
  %   comments revised.
  ```
- Removed: the 11/19 settings entry, the release-pass entry, and the unmarked commented-out line `% PDE_t.BC{2}.term{1}.D = 1;` (not 2025 code: active in 86ad0cc1, commented out in fa8c1ef5).
- Observed, not fixed (code): `if logvals` after `logval = app.loadData(data);` (2025 code) – an undefined variable whenever GUI = 1.

### Stability_Examples/PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi_2.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 0/1/7 (prep: 7 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: new output st, the
  %   recommended settings of the stability test; stability_form = 'P' set
  %   in the base workspace; default R written as 30 (was 21+9, the same
  %   value); help and comments revised.
  ```
- Removed: the 11/19 settings entry and the release-pass entry.

### Stability_Examples/PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq.m
- class condense; TAG `MMP, 11/21/2026`; apply 0/1/8 (prep: 8 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % MMP, 11/21/2026: Changes since PIETOOLS 2025: new output st, the
  %   recommended settings of the stability test; stability_form = 'P' set
  %   in the base workspace; help and comments revised.
  ```
- Removed: the 11/19 settings entry and the release-pass entry.

### Stability_Examples/PIETOOLS_PDE_Ex_Reaction_Diffusion_Gahlawat.m
- class condense; TAG `MMP, 11/21/2026`; apply 0/1/9 (prep: 9 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % MMP, 11/21/2026: Changes since PIETOOLS 2025: new output st, the
  %   recommended settings of the stability tests; stability_form = 'P' set
  %   in the base workspace; help and comments revised.
  ```
- Removed: the 11/19 settings entry and the release-pass entry.

### Stability_Examples/PIETOOLS_PDE_Ex_Reaction_Diffusion_Robin.m
- class condense; TAG `MMP, 11/21/2026`; apply 0/1/8 (prep: 8 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % MMP, 11/21/2026: Changes since PIETOOLS 2025: new output st, the
  %   recommended settings of the stability tests; stability_form = 'P' set
  %   in the base workspace; help and comments revised.
  ```
- Removed: the 11/19 settings entry and the release-pass entry.

### Stability_Examples/PIETOOLS_PDE_Ex_Saba.m
- class condense; TAG `MMP, 11/21/2026`; apply 0/2/22 (prep: 22 RETAG, 2 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % MMP, 11/21/2026: Changes since PIETOOLS 2025: case 4.5 takes sig2 =
  %   1.048, the value of its source (Peet, Automatica 2021, Sec. 10, Ex. 5),
  %   in place of 2.048; new output st, the recommended settings of the
  %   stability test by case; stability_form = 'P' set in the base
  %   workspace; help and comments revised.
  ```
- Removed: the 11/19 settings entry, the 11/19 sig2 entry (eigenvalue measurement, maintainer quote), the release-pass entry. Literature citation kept in the note.

### Stability_Examples/PIETOOLS_PDE_Ex_Timoshenko_Beam_1.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 0/1/10 (prep: 10 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: new parameter c
  %   (default 0), the damping term -c*x3 = -c*phi_{t} in the dynamics of x3,
  %   in both formats; new output st, the recommended settings of the
  %   stability test; stability_form = 'P' set in the base workspace; help
  %   and comments revised.
  ```
- Removed: the 11/19 settings entry and the release-pass entry.

### Stability_Examples/PIETOOLS_PDE_Ex_Timoshenko_Beam_2.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 3/2/15 (prep: 15 RETAG, 2 KEEP-WAS, 3 DROP-WAS – the dropped (was) lines hold intermediate help/comment text of DB and the 2025 "Unstable." help line, not code); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the bending term
  %   x4_{ss} (E*II*phi_{ss}) enters the dynamics of x3 in both formats (the
  %   batch A2 was zero and the term-based format lacked it); new parameter
  %   c (default 0), the damping term -c*x3 in the same equation; new output
  %   st, the recommended settings of the dual stability test;
  %   stability_form = 'Q' set in the base workspace; help and comments
  %   revised.
  ```
- Removed: the 11/19 entry on the bending term (commit hash 72de7c6b, maintainer quote, "grew linearly"), the 11/19 settings entry, the release-pass entry. Kept: the help lines on stability with literature citations.

### Stability_Examples/PIETOOLS_PDE_Ex_Transport_Eq.m
- class condense; TAG `MMP, 11/21/2026`; apply 0/1/8 (prep: 8 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % MMP, 11/21/2026: Changes since PIETOOLS 2025: new output st, the
  %   recommended settings of the stability tests; stability_form = 'Q' set
  %   in the base workspace; help and comments revised.
  ```
- Removed: the 11/19 settings entry and the release-pass entry.

### Stability_Examples/PIETOOLS_PDE_Ex_Wave_Eq_Boundary_Damped.m
- class condense; TAG `MMP, 11/21/2026`; apply 0/1/10 (prep: 10 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % MMP, 11/21/2026: Changes since PIETOOLS 2025: new output st, the
  %   recommended settings of the stability test; stability_form = 'P' set
  %   in the base workspace; help and comments revised.
  ```
- Removed: the 11/19 settings entry and the release-pass entry.

### Stability_Examples/PIETOOLS_PDE_Ex_Wave_Eq_Datko_Boundary_Damped_1.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 0/1/9 (prep: 9 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the term-based format
  %   does not set the field PDE_t.PDE (an assignment of an older input
  %   format, commented out); new output st, the recommended settings of the
  %   dual stability test; stability_form = 'P' set in the base workspace;
  %   help and comments revised.
  ```
- Removed: the 11/19 settings entry and the release-pass entry. Kept: the unmarked line `%PDE_t.PDE.A{4}.coeff = 1; ...` (2025 code commented out by DB in 8f3de467).

### Stability_Examples/PIETOOLS_PDE_Ex_Wave_Eq_Datko_Boundary_Damped_2.m
- class condense; TAG `MMP, 11/21/2026`; apply 0/1/10 (prep: 10 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % MMP, 11/21/2026: Changes since PIETOOLS 2025: new output st, the
  %   recommended settings of the stability test; stability_form = 'P' set
  %   in the base workspace; help and comments revised.
  ```
- Removed: the 11/19 settings entry and the release-pass entry.

### Stability_Examples/PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m
- class condense (comments only since 2025); TAG `MMP, 11/21/2026`; apply 0/0/0 (prep: no marked line); check PASS, codeident OK.
- Note: `% MMP, 11/21/2026: Changes since PIETOOLS 2025: help and comments revised.`
- Removed: the release-pass entry.
- Observed: this file is referenced by no .m file and keeps +Bxr, while PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m now has -Bxr (a duplicate of the 2025 system; not edited).

### command_line_input_Examples/example_2D.m
- class condense; TAG `DJ, MMP, 11/21/2026` (code change by DJ in 791da3fb without a note; the only new note was the release-pass entry); apply 0/0/0 (prep: no marked line); check PASS, codeident OK.
- Note:
  ```
  % DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: the delay of 0.5 acts
  %   on the ODE input, u1(t-0.5), in place of the boundary input u2, which
  %   enters the BC at s1 = 0 undelayed; help and comments revised.
  ```
- Removed: the release-pass entry.

### command_line_input_Examples/l2gain_1.m, l2gain_3.m, l2gain_4.m
- class condense (comments only since 2025); TAG `MMP, 11/21/2026`; apply 0/0/0 each (prep: no marked line); check PASS, codeident OK (all three).
- Note (each): `% MMP, 11/21/2026: Changes since PIETOOLS 2025: help and comments revised.`
- Removed: the release-pass entry.

### command_line_input_Examples/l2gain_2.m
- class condense; TAG `MMP, 11/21/2026`; apply 0/1/2 (prep: 2 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % MMP, 11/21/2026: Changes since PIETOOLS 2025: the PDE reads x_t = x_ss +
  %   s*w, as its comments state (it declared -x_ss, a backward heat
  %   equation); the stated gain is 0.2083 (was 0.3333); help and comments
  %   revised.
  ```
- Removed: the release-pass entry and the 11/20 entry (manual section reference, "at w = 0").

### command_line_input_Examples/stability_1.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 0/0/0 (prep: no marked line); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the script does not
  %   clear the workspace (was clc; clear); the PDE reads x_t = v*x_s with
  %   v = -1 by default (was x_t = x_s); the entries of a workspace cell
  %   params of assignment strings are evaluated after the defaults, and a
  %   true workspace variable GUI loads the GUI save file; the results are
  %   named PDE and PIE (were pde and pie); help and comments revised.
  ```
- Removed: the release-pass entry. Kept: `%clc; clear; ` (2025 code commented out by DB, unmarked) and the commented-out sys()/state block of 2025 (DB edited its sign).
- Observed: the 2025 comment `% PDE: x_{t} = x_{s}` above the declaration no longer matches the code x_t = v*x_s (2025 text, left).

### command_line_input_Examples/stability_2.m, stability_5.m, stability_9.m, stability_10.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 0/1/1 each (prep: 1 RETAG, 1 KEEP-WAS each); check PASS, codeident OK (all four).
- Note (each):
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the script clears the
  %   workspace except a cell params of assignment strings, whose entries
  %   are evaluated after the default parameters; the results are named PDE
  %   and PIE (were pde and pie); help and comments revised.
  ```
- Removed: the release-pass entry and the 11/20 clearvars entry (stability_10: with its help remark on x(s=0) = -xo, covered by "help and comments revised").

### command_line_input_Examples/stability_3.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 0/1/1 (prep: 1 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the script clears the
  %   workspace except a cell params of assignment strings, whose entries
  %   are evaluated after the default parameters; default gains k1 = -0.5,
  %   k2 = -4.5 (were 1, 2); the BCs read [x1(1); x2(0); x3(0)] = K [x1(0);
  %   x2(1); x3(1)] (were [x1(0); x2(1); x3(1)] = K [x1(1); x2(1); x3(0)]);
  %   the results are named PDE and PIE (were pde and pie); help and
  %   comments revised.
  ```
- Removed: the release-pass entry and the 11/20 clearvars entry.
- Observed: Lm is built with the constant 10 in place of g (2025 code), so a params entry for g does not reach Lm; the library function PIETOOLS_PDE_Ex_Diagne.m uses g (equal at the default g = 10).

### command_line_input_Examples/stability_4.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 0/1/1 (prep: 1 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the script clears the
  %   workspace except a cell params of assignment strings, whose entries
  %   are evaluated after the default parameters; the BCs read x1(0) =
  %   qb*x2(0), x2(1) = pb*x1(1) (were x1(0) = qb*x2(1), x2(1) = pb*x1(0));
  %   the results are named PDE and PIE (were pde and pie); help and
  %   comments revised.
  ```
- Removed: the release-pass entry and the 11/20 clearvars entry.

### command_line_input_Examples/stability_6.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 0/1/1 (prep: 1 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the script clears the
  %   workspace except a cell params of assignment strings, whose entries
  %   are evaluated after the default parameters; the results are named PDE
  %   and PIE2 (were pde and pie); help and comments revised.
  ```
- Removed: the release-pass entry and the 11/20 clearvars entry.

### command_line_input_Examples/stability_7.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 1/1/3 (prep: 3 RETAG, 1 KEEP-WAS, 1 DROP-WAS – the dropped (was) line is intermediate code of DB with lam = 4.65 and c on one line); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the script clears the
  %   workspace except a cell params of assignment strings, whose entries
  %   are evaluated after the default parameters, and c is computed from lam
  %   after them; default lam = 4.65 (was 4.66); the results are named PDE
  %   and PIE (were pde and pie); help and comments revised.
  ```
- Removed: the release-pass entry and the 11/20 entry.

### command_line_input_Examples/stability_8.m
- class condense; TAG `DB, DJ, MMP, 11/21/2026` (DJ changed the ODE input in 791da3fb without a note); apply 0/1/1 (prep: 1 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % DB, DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: the script clears
  %   the workspace except a cell params of assignment strings, whose
  %   entries are evaluated after the default parameters; the BCs act at
  %   s = a, b and the ODE input at s = a, with a = 0, b = 1 (a was 1 and
  %   unused); the ODE input reads Bxr*x_{s}(a), the derivative taken before
  %   the substitution (it differentiated x(0) in s); the results are named
  %   PDE and PIE (were pde and pie); help and comments revised.
  ```
- Removed: the release-pass entry and the 11/20 entry. Kept: help line 1 with the stability range of the coded (+Bxr) system.
- Observed: this script keeps +Bxr, while PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m has -Bxr since 11/19 (the help states the script is unstable at its default lam).

### command_line_input_Examples/stability_11.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 1/1/2 (prep: 2 RETAG, 1 KEEP-WAS, 1 DROP-WAS – the dropped (was) line `%R = 2,9;...` is intermediate code of DB); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the script clears the
  %   workspace except a cell params of assignment strings, whose entries
  %   are evaluated after the default parameters; default R = 2.9 (was
  %   2.7); the BC at s = 1 reads x(1) = 0 (was x_{s}(1) = 0); the results
  %   are named PDE and PIE (were pde and pie); help and comments revised.
  ```
- Removed: the release-pass entry and the 11/20 entry.

### command_line_input_Examples/stability_12.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 0/1/1 (prep: 1 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the script clears the
  %   workspace except a cell params of assignment strings, whose entries
  %   are evaluated after the default parameters; defaults Cm = [0 0 0;
  %   s 0 0; s -s^2 0], R = 30 (were [0 0 0; s 0 0; -s^2 0 0], 21); the
  %   results are named PDE and PIE (were pde and pie); help and comments
  %   revised.
  ```
- Removed: the release-pass entry and the 11/20 clearvars entry. Kept: the unmarked `%Cm = [0,0,0;s,0,0;-s^2,0,0];     R = 21;` (the 2025 assignments, reordered and commented out by DB).

### command_line_input_Examples/stability_14.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 0/1/1 (prep: 1 RETAG, 1 KEEP-WAS); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the script clears the
  %   workspace except a cell params of assignment strings, whose entries
  %   are evaluated after the default parameters; new parameter c (default
  %   0), the damping term -c*x3 in the dynamics of x3; the results are
  %   named PDE and PIE (were pde and pie); help and comments revised.
  ```
- Removed: the release-pass entry and the 11/20 clearvars entry.

### command_line_input_Examples/stability_15.m
- class condense; TAG `DB, MMP, 11/21/2026`; apply 3/2/5 (prep: 5 RETAG, 2 KEEP-WAS, 3 DROP-WAS – the dropped (was) lines hold comment text, not code); check PASS, codeident OK.
- Note:
  ```
  % DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the script clears the
  %   workspace except a cell params of assignment strings, whose entries
  %   are evaluated after the default parameters; the bending term x4_{ss}
  %   enters the dynamics of x3; the results are named PDE and PIE (were pde
  %   and pie); help and comments revised.
  ```
- Removed: the 11/19 bending-term entry (file reference, literature sections; the stability statement stays in the help), the release-pass entry, the 11/20 clearvars entry.

### command_line_input_Examples/stability_13.m
- class asnew (18% of the 2025 code left); no note, no TAG; codeident OK; no dated marker or log entry left except `% Initial coding SS - 11/17/2022`.
- Removed (by hand, Edit tool): the release-pass entry, the 11/20 entry (clearvars, eq_PDE ';'), the two (was) lines (`%clc; clear; ...` of 2025, `%  diff(x2,t,1)==...]` of DB), the markers of the clearvars line and of the eq_PDE line. Kept: license, help header, contact and Initial coding lines; the commented-out sys()/state block of 2025 at the end (alternative syntax, no developer note).

## Summary

- Classes: 39 condense (all edited), 1 asnew (stability_13.m, edited), 0 review, 0 deferred. Skipped: none.
- TAGs: `MMP, 11/21/2026` (15 files), `DB, MMP, 11/21/2026` (22), `DJ, MMP, 11/21/2026` (example_2D.m),
  `DB, DJ, MMP, 11/21/2026` (stability_8.m).
- apply totals over the 39 condense files: 14 (was) lines dropped, 37 kept, 222 markers retagged; every count equals the
  prep verdicts; no REVIEW-BEGIN/END line in the group.
- Verification of all 40 files: check PASS (39), codeident OK (40), `i/lf w/crlf` unchanged (40),
  `git diff --numstat` = `git diff -w --numstat` (40), no 11/19 or 11/20 marker left; MATLAB parse check
  (`syntax_check`, lock runner): 40 files checked, 0 with parse problems.
- Internal notes removed: every new header entry (settings searches, "maintainer: no default settings", maintainer
  quotes, Chebyshev node counts, eigenvalue measurements, a commit hash, manual section references, the release-pass
  entries); one unmarked commented-out line not in 2025 (`% PDE_t.BC{2}.term{1}.D = 1;`, Ahmadi). The new help text and
  step comments of the release pass were kept: none is a development or test record.

## Open questions for the maintainer

1. The code changes of DB (the 23 files whose TAG names DB) carry no per-line markers; the pass retags existing markers
   only, so these lines stay unmarked (the note names the changes). The unmarked commented-out 2025 code of DB (Diagne, Datko_1, stability_1,
   stability_12) is kept as 2025 code, unmarked.
2. Help text keeps two computed stability ranges: Heat_Eq_with_ODE_Das (lam < 9.7946; the +Bxr of the paper stable only
   for lam < -5.7799) and stability_8 line 1 (lam < -5.7799). Kept as system properties; drop if read as test numerics.
3. Not fixed (code, observed only): `if logvals` in PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi.m (2025 typo, errors at
   GUI = 1); stability_8.m and the unreferenced PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m keep +Bxr while
   PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m has -Bxr; stability_3.m builds Lm with 10 in place of g.

## Rerun with the fixed cpass.pl (coordinator request, 11/21/2026)

The first preps ran on a private copy of cpass.pl with private temp files (no shared-diff race). After the fix
("DJ (CC)" initials, markers without a blank, per-process temp files): prep rerun on all 40 files – code diffs identical
to the first preps, proposed TAGs equal to the TAGs used, every marked line carries its TAG, no DROP-WAS, REVIEW,
"(CC)" or glued marker, no header entry other than the note; apply rerun on the 39 condense files – no file changed;
check PASS (39), codeident OK (40). Every note re-read against the fresh code diff: no correction needed.
