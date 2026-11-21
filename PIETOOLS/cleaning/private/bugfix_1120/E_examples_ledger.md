# E_examples ledger – bug-fix pass 11/20/2026

Agent E_examples. Items: C74, C75, C76, C77, C79, C84 (Euler-Bernoulli beam file only), C85, C86, C87,
C89, C90, C92 (stability_8.m header text only), C93, C94, C95. HEAD 17d68747.
Scratch: `/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/E_examples/`
(below: `$S`). `CL/` = `PIETOOLS_examples/Examples_Library/command_line_input_Examples/`,
`EL/` = `PIETOOLS_examples/Examples_Library/`.

Pre-edit state: every file below clean in `git status --porcelain`; all `i/lf w/crlf` (`$S/eol_before.txt`).
New-file regime (no developer notes; `new_files_classified.tsv`): `CL/stability_13_2.m`, `16_1`, `16_3`, `17`,
`19`, `36`, `39`, `PIETOOLS_demos/DEMO12_2DControl.m`, `DEMO13_Hinf_OFcontrol.m`. All other files: full protocol.

(items recorded below as they are completed)

Harness (scratch): `cl_run.m` runs a command-line script inside a function workspace, as a user runs it
from the command line, with `params` present or absent; `pde_dump.m` prints every term of the PDE
(coefficients C, D, loc, I). Scripts: `r1_before.m`, `r1b_before.m` (before any edit, outputs
`r1_before.txt`, `r1b_before.txt`, baseline dumps `r1_base.mat`), `r2_after_cl.m` (after, `r2_after_cl.txt`).

## C76 – FIXED (21 files)

Caller: no code in the repository calls the `CL/` scripts (grep of all `.m` files outside `cleaning/`);
`examples_PDE_library_PIETOOLS.m` calls the `EL/*/PIETOOLS_PDE_Ex_*` functions with `(GUI,params)` as
arguments. The `CL/` scripts are run by the user; each help documents "If a cell array params of
assignment strings exists in the workspace, its entries are evaluated after the default parameters"
(params added by DB on 07/08/2026, fa8c1ef5, "as in the old library"; DB commented out `clc; clear;`
in stability_1, 13_3, 18 only). The 21 files test no other workspace variable (`exist(` only for params).

Reproduction (measured, HEAD, `r1_before.txt`): stability_11 with params `{'R=2.5;'}`: R after the run
= 2, PDE dump identical to the default run (x_ss coefficient `[0.5 0;0 0.5]` in both); stability_5 with
`{'lam=1;'}`: lam after the run = 9.86, dumps identical (reaction coefficient 9.86 in both).

Fix: `clc; clear; clear stateNameGenerator` -> `clc; clearvars -except params; clear stateNameGenerator`
(keeps the clearing of every other variable, which also protects against workspace variables
shadowing functions; only the documented input survives).
- Existing (SS), full protocol (header log line, (was) line, col-77 markers): stability_2 (l.44-49),
  3 (46-51), 4 (42-47), 5 (41-46), 6 (38-43), 7 (41-47), 8 (39-45), 9 (41-46), 10 (40-46), 11 (41-46),
  12 (43-48), 13 (42-47), 14 (58-63), 15 (60-65).
- New (DB), code only: stability_13_2 l.52, 16_1 l.39, 16_3 l.38, 17 l.46, 19 l.50, 36 l.41, 39 l.44.

Verification (measured, `r2_after_cl.txt`): for each of the 21 files, one parameter overridden by params
and compared with an independent route, a scratch copy of the file with the default value replaced in
the source (`$S/ref76/`): value after the run = override in 21/21; PDE dump with params = dump of the
copy in 21/21; PIE with params == PIE of the copy (copvar `eq`) in 21/21; dump with params differs from
the default dump in 20/21 – stability_15 is the exception, as its help states (k, aa, II, g reach no
equation; params is evaluated, k = 2 after the run). Default runs equal the pre-edit baseline in 20/21
(stability_11 differs by C74 only). stability_15.m:65 (unused k, aa, II, g; stated in the help since
11/19): left as is – no wrong result, removal is outside the minimal fix.

## C77 – FIXED

`CL/stability_7.m` (existing, SS): c moved after the params loop (l.52-53 (was)/new without c, l.61 new
c line with comment `% c(s,lam) after params`); help l.12-13: "(c is computed from lam before that)"
(was) -> "and c is computed from lam after them". An override of c itself is now overwritten by the
computation (inferred from the code; the help says so). Masked at HEAD by C76, so no separate before run.
Verification (measured): params `{'lam=1;'}`: constant term of c 5.35 (= 0.7 + 4.65) -> 1.7 (= 0.7 + 1),
equal to the copy with `lam = 1;` in the source (dump and PIE).

## C74 – FIXED

`CL/stability_11.m:48` (existing, SS) `R = 2,9;` -> `R = 2.9;` (l.51-52, (was) + marker; header log shared
with C76). Reproduction (measured): console printed `R =` / `2`; x_ss coefficient 0.5 = 1/2. After: x_ss
coefficient 0.34482759 = 1/2.9; the console no longer prints R. Reference: help l.9 "R = 2.9" and the
library counterpart `PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi.m:82` `R = 2.9`.

## C89 – FIXED (3 files)

`CL/stability_16_1.m:86` (new, DB): `...-gam*subs(diff(x3,s,1),s,1)]` -> `...==0];`, comment `= 0` added;
`CL/stability_13.m:59` (existing, SS): `;]` -> `;];` (l.62-63, (was) + marker; header log shared with C76);
`CL/stability_17.m:63` (new, DB): `;] ` -> `;];`.
Reproduction (measured, `r1b_before.txt`): the three consoles print the equation arrays before
initialize (16_1: the 8 BCs; 13: the 2 PDEs; 17: the 2 BCs); a scratch copy with `==0` and `;` gives an
identical PDE dump (the term concatenated with equations becomes `term == 0`, as CANDIDATES states) and
no such output. After (measured): the files print nothing before "Encountered", dumps equal the pre-edit
baseline.

## C93 – FIXED (text only)

`CL/stability_10.m:11` (existing; line by DB 07/08/2026): `x(s=0) = xo` -> `x(s=0) = -xo` (l.11-12,
(was) + marker after one space, the line passing column 75; header log shared with C76). Reference: code
l.58 `subs(X,s,0)==-x` (SS 2022), help l.3, library table case 10 and
`PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Tang.m` BC 1 `0 = x(0) + xo`.

## C92 – FIXED (stability_8.m header text only); code line 62 (DJ) and the Das files HELD

Measured, without PIETOOLS (`c92_spectrum.m`, `c92_sweep.m`; outputs `c92_spectrum.txt`, `c92_sweep.txt`):
system as coded, xo_t = A xo + Bxr X_s(0), X_t = lam X + X_ss + Bpv xo, X(0) = X(1) = 0.
- Chebyshev collocation, lam = pi^2-1: rightmost eigenvalue real, +2.9565924541 at N = 16, 24, 32
  (N = 8 to 64 agree to 3e-9). Second route, sine Galerkin (K = 25 ... 800): +2.9408 ... +2.9561,
  converging to the same value at rate 1/K.
- Stability limit (+Bxr): spectral abscissa crosses 0 once on [-200, 30] (342 grid points), at
  lam = -5.779902 (bisection, N = 48; N = 64 gives 3e-10 there); negative for every grid lam <= -5.8.
- Reference only (the -Bxr sign of the 11/19 Das library file): crossing at lam = 9.794615, below
  pi^2 = 9.869604 (at lam = pi^2-0.01 the abscissa is +0.054); so "stable for lam < pi^2" holds for
  neither sign. For the maintainer (Das files are DJ's: not edited).
Fix: `CL/stability_8.m:1` (existing; line by DB 07/08/2026) "(stable for lam < pi^2 = 9.8696)" (was) ->
"(as coded: stable for lam < -5.7799, unstable at lam = pi^2-1)" (l.1-2, markers after one space).
Header log shared with C76. The library table entry of case 8 (`examples_PDE_library_PIETOOLS.m`, Das,
DJ) carries the same claim: not in this item, not edited.

## C75 – FIXED

Intent established: the inline comment (l.44) and the block comment (l.38) say `x_{t} = x_{ss} + s*w(t)`;
the user manual (2025, Sec. 4.2, GUI) shows the same system (z = int x ds, disturbance s w(t)) with
+x_ss as the library file `PIETOOLS_PDE_Ex_Heat_Eq_with_Distributed_Disturbance` (library case 22,
`A2 = 1`); the 2025 release and `git log -L` show `-diff(x,s,2)` unchanged since SS 11/17/2022
(cdc476f6), copied into the legacy commented block; -x_ss is a backward heat equation (ill-posed).
Reproduction (measured, `r1b_before.txt`): x_ss coefficient -1; against case 22: T, Bw, Cz equal, A
equal to -A(case 22). Fix: `CL/l2gain_2.m:44` (existing, SS) `==-diff(x,s,2)` -> `==diff(x,s,2)`
(l.46-47 (was) + marker, header log l.33-34). Verification (measured): coefficient +1; PIE == PIE of
case 22 (T, A, Bw, Cz equal). Not changed (outside the item): the block comment l.38 states
"(gamma = 0.3333)"; the same system in the library file states 0.2083 = 5/24 (corrected 11/05/2026).

## C84 – FIXED (Euler-Bernoulli beam file); the Saba file HELD (DJ)

`EL/Stability_Examples/PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq.m:191` (existing; line by DB 06/19/2026,
840864ec): `num2str(indx)` -> `num2str(index)`, and the missing space before "does not exist" added
(same line). l.193-194 (was) + markers after one space (the line passes column 75); header log l.76-77.
Reproduction (measured, `r1b_before.txt`): `(14,0,{})` raised `MATLAB:UndefinedFunction`
"Unrecognized function or variable 'indx'." After (measured, `r3_after_misc.txt`): "Example 14 does not
exist, please choose from 13.1 and 13.2"; INDEX 13, 13.1, 13.2 still return PDE and st.

## C79 – FIXED

`EL/2D/PIETOOLS_PDDE_Ex.m` (existing, DB): `tau = [0.25;0.3];` (with its comment) moved before the params
loop (l.85-87 new, markers); old l.90-93 (comment, unused `pvar s t`, comment, tau) kept as (was) lines
(l.96-99); help l.10-11 (was) -> l.12 "The delays tau are set before PARAMS is evaluated, so PARAMS can
change them."; header log l.70-71. `pvar s t` was unused (the code uses `'t'` as a char and s1, s2).
Reproduction (measured): params `{'tau=[0.5;0.6];'}` gave the default dump (transport coefficients 4,
3.3333). After (measured): 2 and 1.6666667, equal to a copy with `tau = [0.5;0.6]` in the source
(`$S/ref79/`, resolved first on the path); default dump unchanged; st unchanged (light, dD = 0).
Run with the file settings: SKIPPED – the 11/19 library record (PHASE_C_TABLE.md, PDE 50) gives the H2
norm program at 503,633 decision variables and its dual at 1.66e6, expected over 5 minutes; the PDE was
built and converted to a PIE (measured).

## C85 – FIXED (removed the unused st line)

Evidence (measured, `c85_demo3.m`, `c85_demo3.txt`): DEMO3 as is (lpiprogram default tier, light):
sqrt(gam) = 0.3183149800, 469 decision variables, ratio to 1/pi 1.000016 (header: the exact value 1/pi;
manual 2025, Sec. 11.3: 0.4271 with the 2025 code, minimum 1/pi = 0.3183). Scratch copy passing st
(`lpiprogram(H2op.vars,H2op.dom,st)`, heavy): the same 0.3183149800 and 469 variables. lpiprogram stores
light by default and heavy when given (`r3_after_misc.txt`); the two tiers differ, but no lpi.ineq field
differs between them, which accounts for the identical program. So the line has no effect at either
tier; passing st would also have edited l.103, a line of DJ (1bf27bbe, 11/06). Fix: l.99-100 (MMP
11/18 and 11/20) kept as (was) lines (l.101-102); header log l.68-69. Not touched: l.79-82 (commented
copy of the BVP; the protocol keeps old code). DEMO3 open in the maintainer editor: edited on disk with
the Edit tool. After (measured): 0.3183149800, 469 variables, no st in the workspace, 2.2 s.

## C86 – FIXED

`PIETOOLS_demos/DEMO5_Hinf_optimal_estimator.m:160` (existing; MMP 11/18): `opts.method='qp'` ->
`opts.method='qp';` (l.161-162 (was) + marker; header log l.90). Semantics: a statement without a
terminating semicolon displays its value; 'qp' is the default method of getObserver_sop (l.92-93 there)
and no file under `PIESIM/` reads `.method` (grep), so the line affects the console only.
After (measured, `c86_demo5.txt`): DEMO5 runs to the end (18.9 s), gam_val = 1.0004121, no struct display
of opts in the console (only the echo of the line, ECHO ON).

## C87 – FIXED (2 files)

`EL/PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion.m:104`, `..._Control.m:99` (existing; MMP 11/19):
`st.sop.deg_Q{1} = [6 2 2];` kept as (was) line; the clause "Q of degree 6 in the multiplier and 2 in
the integral kernels" removed from the settings comment ((was) lines + new line); header logs.
Read: st.sop.deg_Q is read only by `executives_sopvar/PIETOOLS_stability_sop.m:173` (grep of the
executives, LPI and settings folders; exec_sop_settings only defaults it); these files set H2_estimator,
H2_control only. Reproduction and verification (measured, library path via `libcheck.m` copied to
scratch, `lib_before.csv`, `lib_after.csv`): case 51 H2_estimator optimal 1.2695689, 1293 decision
variables, m = 103, before and after; case 52 H2_control optimal 1.0907997, 2471, m = 142, before and
after – identical programs, so the field had no effect.

## C90 – FIXED

`EL/PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance2.m:68` (existing; DB): `pvar s1 s2` kept as (was)
line (l.69); header log l.61. Semantics: s1, s2 are not referenced after l.68 (the PDE uses `pvar s t`,
l.84). Run with the file settings (veryheavy, measured): H2_norm 0.10157583, H2_norm_dual 0.10157588,
11,847 decision variables, identical before and after.

## C94 – FIXED (new file, code and its comment only)

`PIETOOLS_demos/DEMO12_2DControl.m:112`: prints `sqrt(norm(full(norm_err.C{1}.params{1})))`; comment
l.108-109 now "Residual err = Z - K*P, a map L2 -> R^ne: ||err||^2 is the 2-norm of the ne x ne matrix
err*err'". Reproduction (measured, `c94_demo12_before.txt`): norm_err is one block with one 1 x 1
parameter; the printed value 5.626e-23 is ||Z - KP||^2, labelled as the norm. Semantics: err maps
L2[s1,s2] -> R^ne, err*err' is the ne x ne matrix whose 2-norm is ||err||^2. After (measured,
`c94_demo12.txt`): prints 7.501e-12 = sqrt(5.626e-23). Second route for the expression, on the nonzero
operator Zsol of the same run: sqrt(norm(params{1} of Zsol*Zsol')) = 1.056720592518e-02, Gauss-Legendre
quadrature (20 x 20) of the squared kernel of Zsol' = 1.056720592518e-02 (relative difference 1.6e-16).
For err itself both routes sit at rounding level (QP objective -4.1e-20, its cancellation floor).
DEMO12 runs in 10 s.

## C95 – FIXED (new file, text only)

`PIETOOLS_demos/DEMO13_Hinf_OFcontrol.m:220` (line by DB): label 'Hinf_OFcontrol_demo12' ->
'Hinf_OFcontrol_demo13'. After (measured, `c95_demo13.txt`): DEMO13 runs to the end (11.9 s) and prints
"Hinf_OFcontrol_demo13: gam = 5.6252653 (certified: 1, dw = 0)".

## Declan hold – origin of every removed line (coordinator request, 11/20 evening)

`git blame -w -C -C -M` at HEAD over every removed line of the 31 files; lines attributed to DJ:
- `clc; clear; clear stateNameGenerator` in stability_2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15 (plain blame
  names DJ in 2, 5, 6, 7, 8, 10, 11 and DB's fa8c1ef5, which moved the line under the new help, in 4, 9,
  12, 14, 15): content from f7ae0e1f (DJ, 01/05/2025, "Update command line input examples"), the
  conversion of the 19 files of the folder from sys/state to pde_var. Its edit of this line is the same
  in all 19 files: SS's `clc; clear;` (86c8c3a0, 11/16/2022: stability_2, 5, 6; 2bdfad9e, 11/17/2022:
  stability_7, 8, 10, 11) gains `clear stateNameGenerator`, the reset of the ID counter used by pde_var
  (`parser/internal/stateNameGenerator.m`). A format-conversion sweep; the fix changes only SS's `clear;`
  and keeps DJ's `clear stateNameGenerator` verbatim. EDIT STANDS. (Reversal, if wanted: the (was)
  line above each changed line holds the original.)
- `EL/2D/PIETOOLS_PDDE_Ex.m` l.90-91 (`% Declare independent variables (time and space)`, `pvar s t`):
  `git log -L` and plain blame give 403b1bbf (DB, 01/20/2025, "Added PDDE example"); copy detection traces
  the two lines to `PIETOOLS_demos/DEMO8_H2norm.m` of 65fce72c (DJ, 12/28/2024), from which DB copied
  them into his example. DB logic in this file. EDIT STANDS.
All other removed lines: DB (fa8c1ef5, d676cac3, 598dfeb3, 840864ec, 403b1bbf, 9efb8fff), SS (cdc476f6),
Talitsky (ed16b798, DEMO12 l.112) or MMP.

## Checks over all 31 edited files

`syntax_check` (list `$S/edited_files.txt`, output `$S/syntax_check.txt`): 31 files, 0 parse problems.
`git ls-files --eol`: all `w/crlf`, as before. Markers at column 77, or one space after code reaching
column 75 (8 lines: stability_8 l.1-2, stability_10 l.11-12, stability_13 l.63, beam l.193-194); none on
a `...` line; no markers, log lines or (was) lines in the 9 new-regime files.
Note: the Edit tool drops trailing blanks in the lines it writes; it removed them from two unchanged lines
of stability_11 (l.48-49, inside an edited block). Restored byte-exactly with a single perl substitution
on those two lines (CRLF kept; afterwards `git diff --ignore-space-at-eol` and `git diff` give the same
6/2 lines). Every other edited file: `git diff --numstat` equals `--ignore-space-at-eol`.
