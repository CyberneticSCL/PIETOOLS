# dE ledger – report for Declan, examples and demos (11/21/2026)

Agent dE. Items: C78, C80, C81, C82, C83, C84 (Saba file), C88, C92 (line 62 of CL/stability_8.m and the
Das files), C96, C97, C98, C99, C100, C115. HEAD 905b0d26 (one commit after fa0deda4; it changes three
executives only, no file of these items). No repository file edited.

`EL/` = `PIETOOLS_examples/Examples_Library/`, `CL/` = `EL/command_line_input_Examples/`.
Scratch `$S` = `/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/dE/`.
Each `$S/<item>_repro.m` writes `$S/<item>_repro.txt`; run with
`bash .../scratchpad/mlrun.sh "run('<windows path of $S>\<item>_repro.m')"`. Helpers in `$S`:
`dE_opeval.m` (opvar-form kernels of a 1-D copvar at 4 x 4 points, via getfield_opvar1d, so that two
operators compare whatever the names of their variables), `dE_cmp_pie.m` (max kernel difference per PIE
operator), `dE_das_alpha.m` (spectral abscissa of the Das system by Chebyshev collocation or finite
differences; no PIETOOLS function), `dE_terms_str.m`, `dE_run_cl.m` (runs a CL script in a function
workspace, as at the prompt), `stub_gui/PIETOOLS_PDE_GUI.m` (stub for GUI = 1; the scripts cd to `$S`
first, since the PIETOOLS root as current folder shadows any path entry), `ref/` (scratch copies of example
files with one change, written by the scripts).

| Item | Status | Consequence |
|---|---|---|
| C78 | CONFIRMED (run) | a: silent wrong result (lam override ignored); b: silent wrong result (labelled defaults ignored), dead code (first loop) |
| C80 | CONFIRMED (run) | error on a reachable path (params b ≠ 1); not silent as CANDIDATES states |
| C81 | NOT DJ (DB, 9efb8fff) | measured: error on GUI = 1 (load fails), not another PDE |
| C82 | CONFIRMED (run) | silent wrong result (snippet PIE3; B on u, no w; T kernels swapped – new) |
| C83 | CONFIRMED (run) | silent wrong result (term format: coupling vanishes); error (GUI = 1) |
| C84 | CONFIRMED (run) | error on a rare path |
| C88 | CONFIRMED (reading) | dead code |
| C92 | CONFIRMED (run) | misleading text; two signs among the copies |
| C96 | CONFIRMED (reading) | misleading text |
| C97 | CONFIRMED (run) | misleading example (formulation) |
| C98 | NOT DJ (DB, 7e551368, both files) | – |
| C99 | CONFIRMED (reading) | misleading name and table text |
| C100 | CONFIRMED (reading) | misleading text |
| C115 | CONFIRMED (run) | misleading text (source claim); GUI file differs from the code (new) |

## C78 – CONFIRMED (run)

HEAD: (a) `EL/Stability_Examples/PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat.m:90` (c from lam) before the loop
92-98; help 36-37 (MMP, 978570e9, 11/20) states it. (b) `EL/Stability_Examples/PIETOOLS_PDE_Ex_Diagne.m`:
first loop 81-87, defaults 91-92 ("% Default Parameters", line 90), derived Lm, Mm, K 95-125, second loop
130-136; help 33 names Mm, Lm, K00-K11 as settable.

Author. (a) Plain blame: Sachin ba9d8807 (folder rename). `-C -C -M`: 12b6036d (DJ, 2022-10-19, cfun -> c),
e21afbef (DJ, 2022-08-20, file created with cfun before the loop). The order exists in the library at
95e49eaf/90ce5fa5 (SS uploads of the 2021 library, 2021-09-29) and adfc1a00 (2021b); the 2020a library
(c587fd94) has no params; the 2021 library header: "Adjusted to function format. DJ - 05-31-2021". DJ.
(b) 6a5d184e (DJ, 2022-08-18, "Partially decompose PDE examples library") created the function with both
loops; the pre-decomposition library (6a5d184e^) had one loop, after the derived coefficients. DJ.

Reproduction: `$S/C78_repro.m` -> `C78_repro.txt` (measured):
- (a) c(0), c(0.5): default [5.35 4.8625]; params {'lam=1;'} [5.35 4.8625]; copy with lam = 1 in the
  source [1.7 1.2125]; params rebuilding c from lam = 1 [1.7 1.2125]; closed form c(0) = 0.7 + lam.
- (b) PDE (A0, A1, B) equal to the default for {'Vs=2;'}, {'Hs=2;'}, {'g=9.81;'}, {'Cf=2;'}, {'Bs=1;'},
  {'k1=1;','k2=2;'}: 1 each; for {'K11=0.5;'}, {'Mm=zeros(3);'}: 0 (take effect). Copy with Vs = 2 in the
  source: differs (max|A0 diff| 0.944); copy with k1 = 1, k2 = 2: differs (max|B diff| 34.75).
- Read: the first loop has no effect – lines 89-127 assign every variable they read, and the second loop
  re-evaluates the same strings. testing_codes lb_meta.m (case 7) rebuilds c to work around (a).

Consequence: silent wrong result for a user who overrides lam (library case 7 lists lam; library help:
"add argument 'lam = 10'") or the labelled defaults of Diagne; dead code (first Diagne loop).

Report entries:

**C78a** `EL/Stability_Examples/PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat.m:90` (DJ, e21afbef, 2022-08-20; order
of the 2021 library function format, DJ 05/31/2021) – c(s) is built from lam on line 90, before the params
loop (92-98). Measured: params {'lam=1;'} gives c(0) = 5.35 = 0.7 + 4.65, the default; a copy with lam = 1
in the source gives 1.7 = 0.7 + 1. A lam override, the parameter of library case 7, is ignored without a
message (the help states it since 11/20/2026; testing_codes lb_meta.m rebuilds c to work around it).
Fix: build c after the params loop, as CL/stability_7.m since 11/20/2026.
Reproduction:
```matlab
[~,P] = PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat(0,{'lam=1;'});
pvar s;  double(subs(P.A0,s,0))      % 5.35; expected 0.7 + 1 = 1.7
```

**C78b** `EL/Stability_Examples/PIETOOLS_PDE_Ex_Diagne.m:81-87, 130-136` (DJ, 6a5d184e, 2022-08-18) – the
first params loop runs before the defaults (91-92), the second after Lm, Mm and K are derived (95-125).
Measured: params {'Vs=2;'}, {'Hs=2;'}, {'g=9.81;'}, {'Cf=2;'}, {'Bs=1;'}, {'k1=1;','k2=2;'} each return the
default PDE (A0, A1, B equal), while copies with Vs = 2 or k1 = 1, k2 = 2 in the source change A0 (by up
to 0.944) or B (by up to 34.75); Mm and K11, named in the help, take effect. Overrides of the values the
code labels "Default Parameters" are ignored without a message; the first loop has no effect.
Fix: move the first loop below line 92.
Reproduction:
```matlab
[~,P0] = PIETOOLS_PDE_Ex_Diagne(0,{});   [~,P1] = PIETOOLS_PDE_Ex_Diagne(0,{'Vs=2;'});
isequal(P0.A0,P1.A0) && isequal(P0.A1,P1.A1) && isequal(P0.B,P1.B)     % true
```

## C80 – CONFIRMED (run); class corrected to error on a reachable path

HEAD: `EL/2D/PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D.m:126-127` `subs(x2,s1,1)`, `subs(x3,s1,1)`; lines
116-117, 120, 127-129 use [a,b] or b; the commented term-based block uses loc = b (172, 177). Help 49-50
(MMP, 11/19-20): "a and b ... set the domain of s1, but the BCs take x2, x3 at s1 = 1, so b = 1 is assumed".

Author: b84cbab2 (DJ, 2025-01-05, "Update 2D PDE stability examples"; conversion to pde_var), plain and
`-C -C -M` blame agree. DJ.

Reproduction: `$S/C80_repro.m` -> `C80_repro.txt` (measured; copy `ref/dE_Wave2D_fixb.m` with b for 1 in
the two lines):
- (a,b) = (0,1): file and copy initialize and convert; PIE T, A equal (copvar eq).
- {'b=2;'}: file: initialize error "The spatial position "BC{3}.term{2}.loc" is not appropriate; for any
  variable s in [a,b], PIETOOLS can only evaluate s=a, s=b, or s=s."; copy: initialize (loc = 2) and
  convert succeed.
- {'a=1;','b=2;'}: file: initialize accepts loc = 1 = a for x2 (BC 3) and x3 (BC 4), where BCs 1-2 set
  x2 = x3 = 0; convert then stops: convert_PIETOOLS_BVP:sourceRankDeficient "Boundary conditions with a
  state of lower dimension along a rank-deficient (e.g. periodic) direction are not supported."; copy:
  initialize and convert succeed.

Consequence: error on a reachable path (any params with b ≠ 1); no silent wrong result was found (the
(1,2) case builds the wrong BCs but stops in convert, with an unrelated message).

**C80** `EL/2D/PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D.m:126-127` (DJ, b84cbab2, 2025-01-05) – the two
feedback BCs take x2, x3 at s1 = 1, while s1 ranges over [a,b] (params; the other lines use b). Measured:
params {'b=2;'}: initialize stops ("The spatial position "BC{3}.term{2}.loc" is not appropriate; ...");
{'a=1;','b=2;'}: initialize places the two terms at s1 = 1 = a, where BCs 1-2 set x2 = x3 = 0, and convert
stops (convert_PIETOOLS_BVP:sourceRankDeficient); a copy with b for 1 in the two lines initializes and
converts in both cases, and equals the file at (a,b) = (0,1) (PIE T, A). Error on a reachable path: any
params with b ≠ 1, though the help names a and b as parameters (it states "b = 1 is assumed" since
11/19/2026). Fix: subs(x2,s1,b), subs(x3,s1,b) in lines 126-127, as the commented term-based block.
Reproduction:
```matlab
PDE_t = PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D(0,{'b=2;'});
PDE = initialize(PDE_t);     % error: spatial position "BC{3}.term{2}.loc" is not appropriate
```

## C81 – NOT DJ (DB, 9efb8fff, 2025-01-06)

HEAD: `EL/PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion.m:89` (was 87):
`load(fullfile(root,'PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq_GUI.mat'))`.
Author: `git log --diff-filter=A` -> 9efb8fff (DanBraghini, 2025-01-06, "Add H2 estimator example") adds
the file with this line (copied from the beam example of DJ, e21afbef, which `-C -C` reports). The header
line "Initial coding DJ - 08/18/2022", the basis of the DJ credit in CANDIDATES, was added on 11/20/2026 by
978570e9 (MMP, release comment pass; `git log -S`): it misattributes the file. Not in the report.
Measured for the maintainer (`$S/C81_repro.m` -> `C81_repro.txt`, GUI app stubbed): GUI = 1 raises
MATLAB:load:couldNotReadFile – the beam GUI file lies in `EL/Stability_Examples/`, not in `EL/` – so the
consequence is an error on the GUI path, not another PDE in the GUI as CANDIDATES states. The sibling
`_Control` file raises "No GUI file associated with this example".

## C82 – CONFIRMED (run); additional defect in T

HEAD: `PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch5_PIEs.m:147`
`PIE3 = piess({T,0,Tu},A,{0,B1},{C1;C2},{0,D12;0,D22});`; T at line 136 (and 169):
`pi4_copvar_sop([],[],[],{[],-s,-s_dum},...)`; comment 146 (MMP, 11/20) describes the code.
Manual 2025, Sec. 5.5.3, p. 57 (PDF p. 58): formula (Tv)(s) = -int_0^s theta v(theta) dtheta - int_s^1
s v(theta) dtheta; code T.R.R1 = -s, T.R.R2 = -theta; call `piess({T,0,Tu},A,{B1,0},...)`; "(B1 w)(s) =
s(2-s)".
Author: line 147 and the T lines: 3b9962b5 (DJ, 2024-12-30, "Create PIETOOLS_Code_Illustrations_Ch5_PIEs.m",
`git log -S"{0,B1}"`); line 136 ported to pi4_copvar_sop by DB (f2f855d8, 2026-11-06) with the same
kernels. DJ.

Reproduction: `$S/C82_repro.m` -> `C82_repro.txt` (measured). Reference: convert of the Sec. 5.4 PDE
(snippet lines 88-98), a route independent of piess; operators compared by their opvar-form kernels.
- Max kernel difference against the PDE PIE, snippet {0,B1} / manual {B1,0}: T 0.76 / 0.76; Tu, A, Cz, Cy,
  Dzu, Dyu 0 / 0; Bu 0.988 / 0; Tw, Bw, Dzw, Dyw: different dimensions / 0.
- dim of Bw: snippet [0 0;1 0] (no w), manual and PDE [0 1;1 0]. Bu kernel Q2: snippet 2*s - s^2, PDE 0.
- T: snippet R1 = -s, R2 = -s_dum; PDE R1 = -s_dum, R2 = -s, as the manual formula and the hand
  derivation x(s) = u - int_0^s t v(t) dt - s int_s^1 v(t) dt (x(0) = u, x_s(1) = 0, v = x_ss).

Consequence: silent wrong result – PIE3 and the closed loops of lines 154-174 are another system than the
PIE of the Sec. 5.4 PDE that the manual describes; the snippet runs without error.

**C82** `PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch5_PIEs.m:147`, with T of lines
136, 169 (DJ, 3b9962b5, 2024-12-30; line 136 ported by DB, f2f855d8, 2026-11-06, kernels unchanged) –
`piess({T,0,Tu},A,{0,B1},...)` puts B1 = s(2-s) on u and leaves PIE3 without a disturbance (Tw, Bw, Dzw,
Dyw have no columns); the 2025 manual (Sec. 5.5.3, p. 57) passes {B1,0}. Measured against convert of the
Sec. 5.4 PDE (lines 88-98): Tu, A, Cz, Cy, Dzu, Dyu agree; Bu differs (Q2 = 2s - s^2 against 0); T differs
(R1 = -s, R2 = -s_dum against R1 = -s_dum, R2 = -s, which the manual formula and a hand derivation also
give; the manual code declares the T of the snippet). With {B1,0} every operator but T agrees. PIE3 and the
closed loops of lines 154-174 are another system than the PIE the manual describes; the snippet runs.
Fix: {B1,0}; T = pi4_copvar_sop([],[],[],{[],-s_dum,-s},[0 0;1 1],s,s_dum,dom) (and the manual code).
Reproduction (after snippet lines 88-98 and 107-147; PIE is the converted Sec. 5.4 PDE):
```matlab
getfield_opvar1d(PIE3.Bw,'dim')                                   % [0 0;1 0]: no w
getfield_opvar1d(PIE3.Bu,'Q2'), getfield_opvar1d(PIE.Bu,'Q2')     % 2*s - s^2, 0
getfield_opvar1d(PIE3.T,'R.R1'), getfield_opvar1d(PIE.T,'R.R1')   % -s, -s_dum
```

## C83 – CONFIRMED (run)

HEAD: `EL/Stability_Examples/PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m:120-122` (term coupling: x = 2,
loc = 0, no D, C = +Bxr), 150 (GUI file name); batch E0 line 95 selects x_s(a).
Unreferenced: `grep -rn "Example_Heat_Eq_with_ODE_Das" --include=*.m` outside cleaning/ finds only the file
itself (lines 1, 36, 150); the library calls the Ex_ file (`examples_PDE_library_PIETOOLS.m:287`).
Author: 120-122: 3e087c84 (DJ, 2022-07-08, "Implement new terms format", in the library; the Tang case of
the same commit has D = 1); the file: 6e55cb23 (DJ, 2022-08-29, "Fix PDE example file.", a copy of the Ex_
file of that date that differs only in the function name and the GUI file name); line 150: 6e55cb23. SS added D = 1 to the Ex_
file only (9bd434ab, 2022-12-04). DJ.

Reproduction: `$S/C83_repro.m` -> `C83_repro.txt` (measured):
- (a) term vs batch format: max kernel difference T 0, A 1.337.
- (b) term format vs the same with Bxr = 0: T 0, A 0; copvar eq: A == A0 true, T == T0 true; Q1 of A (ODE
  row, PDE columns) 0.
- (c) control: the Ex_ file (D = 1), term vs batch: T 0, A 0.
- (d) GUI = 1 (stub): MATLAB:load:couldNotReadFile, `..._Example_Heat_Eq_with_ODE_Das_GUI.mat` does not exist.

Consequence: silent wrong result (term format: coupling has no effect) and error (GUI = 1), for a user who
calls the file directly.

**C83** `EL/Stability_Examples/PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m:120-122, 150` (DJ: 3e087c84,
2022-07-08, lines 120-122; 6e55cb23, 2022-08-29, the file and line 150) – unreferenced copy of
PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m (no file outside cleaning/ calls it). Measured: the term format
couples x(0) (loc = 0, no D), zero by the BC – its PIE equals the PIE with Bxr = 0 (T, A equal) and differs
from the PIE of the batch format of the same file, which couples x_s(0) (A differs by up to 1.34); in the
Ex_ file, where SS added D = 1 (9bd434ab, 2022-12-04), the two formats give the same PIE. GUI = 1: load
fails (MATLAB:load:couldNotReadFile; PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das_GUI.mat does not exist).
Silent wrong result (term format) and error (GUI = 1) for a user who calls the file. Fix: remove the copy
(the library calls the Ex_ file), or D = 1 and the GUI file name of the Ex_ file.
Reproduction:
```matlab
[Pt,Pb] = PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das(0,{});
Pt0 = Pt;  Pt0.x{1}.term{2}.C = zeros(4,2);
P1 = convert_PIETOOLS_PDE(Pt);  P0 = convert_PIETOOLS_PDE(Pt0);
P1.A == P0.A            % true: the coupling has no effect
```

## C84 (Saba file) – CONFIRMED (run)

HEAD: `EL/Stability_Examples/PIETOOLS_PDE_Ex_Saba.m:110`
`error(['Example ',num2str(indx),'does not exist, ...'])`; the input is index (line 1).
Author: plain blame Sachin ba9d8807 (rename); `-C -C -M`: adfc1a00 (2021b); `git log -S"num2str(indx)"`:
first at 90ce5fa5 (SS upload of the 2021 library, 2021-09-29), whose header credits the function format to
DJ ("Adjusted to function format. DJ - 05-31-2021"); moved into the Saba function by DJ (6a5d184e,
2022-08-18; e21afbef). DJ (function format); the beam copy of the line was DB (fixed 11/20).
Reproduction: `$S/C84_repro.m` -> `C84_repro.txt` (measured): INDEX 4.1 returns PDE_t, PDE_b, st; INDEX 5
and 4.6: MATLAB:UndefinedFunction "Unrecognized function or variable 'indx'."; the intended message reads
"Example 5does not exist, ..." (no space).
Consequence: error on a rare path (a wrong INDEX); the undefined-variable error replaces the message.

**C84** `EL/Stability_Examples/PIETOOLS_PDE_Ex_Saba.m:110` (DJ: 2021 library function format, DJ
05/31/2021; moved into this file by 6a5d184e, 2022-08-18) – the error branch calls num2str(indx); the
input is index. Measured: INDEX 5 and 4.6 raise MATLAB:UndefinedFunction "Unrecognized function or
variable 'indx'." in place of the intended message, which also lacks a space before "does not exist".
Error on a rare path (an INDEX outside 4, 4.1-4.5). Fix: index for indx, and the space (as
PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq.m since 11/20/2026).
Reproduction:
```matlab
PIETOOLS_PDE_Ex_Saba(5,0,{})       % Unrecognized function or variable 'indx'.
```

## C88 – CONFIRMED (reading)

HEAD: `EL/2D/Hinf_Gain_Examples/PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DDDD.m:152-161`, `_DDDN.m:152-160`,
`_DNDN.m:153-161` (r_lim, L2_gain); outputs [PDE_t,st] (line 1 of each).
Author: 827b9ca1 (DJ, 2024-08-02, "Add 2D examples"), every line of the three blocks. DJ.
Reproduction: by reading – MATLAB discards the local variables of a function on return; every statement
of the blocks ends in ';', so nothing is displayed. No run.
Consequence: dead code.

**C88** `EL/2D/Hinf_Gain_Examples/PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DDDD.m:152-161`, `_DDDN.m:152-160`,
`_DNDN.m:153-161` (DJ, 827b9ca1, 2024-08-02) – r_lim and L2_gain ("For comparison, compute the L2-gain")
are local variables of a function with outputs [PDE_t,st], and every line of the block ends in ';'. Dead
code: the comparison value never reaches the user. Fix: remove the block or print L2_gain.
Reproduction: by reading: a function discards its local variables on return; no line of the block displays.

## C92 (line 62 of CL/stability_8.m, the Das files) – CONFIRMED (run)

HEAD: `CL/stability_8.m:65` (62 at 78641839) `eq_ODE=diff(x,t)==A*x+Bxr*subs(diff(X,s),s,a);`, Bxr at 63;
header line 1 corrected on 11/20 (E_examples) to "stable only for lam < -5.7799, so unstable at lam =
pi^2-1". Older copy `EL/Stability_Examples/PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m`: +Bxr in the header
(25), "Stable for lam<pi^2" (32), batch E0 (95), term C (122, at x(0), see C83). Library file
`PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m`: -Bxr in both formats since 99f7dff5 (MMP, 11/19).
Author: stability_8:65: 791da3fb (DJ, 2026-11-08, "Update converter."), which replaced
`Bxr*diff(subs(X,s,a),s)` (f7ae0e1f, DJ, 2025-01-05, from SS 2bdfad9e, 2022-11-17) by the x_s(a) coupling
and kept the sign +. Older copy: header, claim and E0 from adfc1a00 (2021b), in the file since 6e55cb23
(DJ). The sign + is that of the source (Das et al. 2018, Ex. 2, per the header of the library file).
Reproduction: `$S/C92_repro.m` -> `C92_repro.txt` (measured):
- stability_8 (run as at the prompt): ODE term 2 after initialize: x = 2, D = 1, loc = 0, C = +Bxr.
- The line before 791da3fb, `Bxr*diff(subs(X,s,0),s)`, is rejected by the current parser ("The equation
  cannot be set: the left-hand side and right-hand side depend on different spatial variables.").
- Older copy: term: x = 2, no D, loc = 0, C = +Bxr; batch: E0 = +Bxr*[0 0 I 0]. Library file: term C = -Bxr
  with D = 1, batch E0 = -Bxr*[0 0 I 0].
- Spectral abscissa at lam = pi^2-1 (outside PIETOOLS, Chebyshev N = 48): +Bxr 2.9565924529, -Bxr
  -0.7122829559, no coupling -0.6803879062. Limits in C115.
Consequence: misleading text (older copy, line 32: true only for its term format, where the coupling
vanishes) and two signs among the copies; stability_8 is unstable at its default lam (its header says so
since 11/20).

**C92** `CL/stability_8.m:65` (DJ, 791da3fb, 2026-11-08; line 62 at 78641839) and
`EL/Stability_Examples/PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m:25, 32, 95` (2021b, in this file since
6e55cb23, DJ, 2022-08-29) – these copies use the ODE input +Bxr*x_s(0) of Das et al. 2018 (Ex. 2); the
library file PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m uses -Bxr since 11/19/2026. Measured: the initialized
PDE of stability_8 has the coupling term C = +Bxr, D = 1, loc = 0; at the default lam = pi^2-1 the +Bxr
system has spectral abscissa +2.9566 (Chebyshev collocation, outside PIETOOLS), so stability_8 and the
batch format of the older copy are unstable at their defaults (limits: +Bxr lam < -5.7799, -Bxr lam <
9.7946, Bxr = 0 lam < pi^2; C115). The older copy keeps "Stable for lam<pi^2" (line 32), true only for its
term format, where the coupling vanishes (C83). Misleading text and two signs among the copies; the header
of stability_8 states the coded range since 11/20/2026. Fix: none (sign decision, C115).
Reproduction:
```matlab
stability_8;  PDE.x{1}.term{2}      % x = 2, D = 1, loc = 0, C = Bxr of line 63
% spectral abscissa: the lines of C115 give sa(pi^2-1,1) = 2.9566
```

## C96 – CONFIRMED (reading)

HEAD: `PIETOOLS_demos/DEMO9_L2_gain_analysis.m:32` "% PIETOOLS - DEMO8". Author: 65fce72c (DJ, 2024-12-28,
squashed commit that adds DEMO9_L2_gain_analysis.m; header "DJ, 11/22/2024: Initial coding"). DJ.
Reproduction: by reading. Consequence: misleading text (DEMO8 is DEMO8_H2norm.m).

**C96** `PIETOOLS_demos/DEMO9_L2_gain_analysis.m:32` (DJ, 65fce72c, 2024-12-28) – the license title reads
"PIETOOLS - DEMO8" in DEMO9 (DEMO8 is DEMO8_H2norm.m). Misleading text. Fix: "PIETOOLS - DEMO9".
Reproduction: by reading: line 32.

## C97 – CONFIRMED (run)

HEAD: `EL/Hinf_Optimal_Controller_Examples/PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_with_Control.m:96`
`z==int((1-s)*v2,s,[0,1])`; u distributed (line 92).
Author: 0e6e7bd4 (DJ, 2025-01-15, "Update Bernoulli beam control example."): u made distributed and u
removed from z; 1aefcaa8 (DJ, 2025-01-13) had `z==[u; int(0.5*(1-s)^2*v2,s,[0,1])]` with a scalar u; the
header kept z = [u; ...] until 11/20 (F2 ledger). DJ.
Reproduction: `$S/C97_repro.m` -> `C97_repro.txt` (measured): PIE: Dzu == 0 true, Dzw == 0 true, Bu == 0
false, Cz == 0 false. PIETOOLS_Hinf_control_sop: 'stripped' gam = 0.000281056 (status unknown, not
certified, 4.8 s); 'light' gam = 0.00025643 (optimal, certified, 1.9 s). Record (not rerun):
`cleaning/private/example_settings_1119/PHASE_C_TABLE.md`, file settings (veryheavy): 0.000169037 at 8,369
variables. Inferred: with D12 = 0 the problem is singular (no penalty on u); the bound decreases as the
settings grow.
Consequence: misleading example (formulation).

**C97** `EL/Hinf_Optimal_Controller_Examples/PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_with_Control.m:96` (DJ,
0e6e7bd4, 2025-01-15) – z = int (1-s) v2 ds has no control term (the first version, 1aefcaa8, 2025-01-13,
had z = [u; int 0.5(1-s)^2 v2 ds] with a scalar u). Measured: in the PIE Dzu (D12) = 0 and Dzw = 0;
PIETOOLS_Hinf_control_sop returns gam = 2.81e-4 at 'stripped' (status unknown) and 2.56e-4 at 'light'
(optimal); the 11/19/2026 record gives 1.69e-4 at the file settings (veryheavy). Inferred: with D12 = 0 the
synthesis problem is singular (no penalty on u), and the bound decreases as the settings grow. Misleading
example (formulation). Fix: none (a formulation decision).
Reproduction:
```matlab
PDE = PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_with_Control(0,{});  PIE = convert(PDE,'pie');
PIE.Dzu == 0                                                          % true
[~,~,gam] = PIETOOLS_Hinf_control_sop(PIE,lpisettings_sop('light'))   % 2.56e-4
```

## C98 – NOT DJ (both files: DB, 7e551368, 2022-11-29)

HEAD: `EL/Hinf_Gain_Examples/PIETOOLS_PDE_Ex_Parabolic_Eq_with_Disturbance.m:149` and
`EL/Hinf_Gain_Examples/PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Disturbance.m:149`
`disp('Warning: With sedumi solver, this example requires at least settings veryheavy.')`; st (line 153)
'stripped' / 'light' (MMP, 11/19).
Author: plain and `-w -C -C -M` blame of line 149 in both files: 7e551368 (DanBraghini, 2022-11-29, "Update
DEMO and examples"), which adds the line to both files (`git log -S"requires at least settings
veryheavy"`). Split as asked: C98a (Parabolic) NOT DJ, C98b (Reaction_Diffusion) NOT DJ. Not run (not in
the report).

## C99 – CONFIRMED (reading)

HEAD: `EL/2D/Hinf_Optimal_Observer_Examples/PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DNDN_EdgeObserve.m:112-113`
`y1==int(x,s2,[c,d])+w2; y2==int(x,s1,[a,b])+w3;`; commented term block 138-144 senses x(b,s2), x(s1,d);
library table case 48 (`examples_PDE_library_PIETOOLS.m:953-954`): "y1 = x(s2=d) + w1; y2 = x(s1=b) + w2".
Author: 006d40eb (DJ, 2025-01-05, "Update 2D PDE examples."): its diff replaces the term-based edge sensing
(827b9ca1, DJ, 2024-08-02) by the integrals; the table lines are 827b9ca1 (DJ). DJ.
Reproduction: by reading. Consequence: misleading name and table text.

**C99** `EL/2D/Hinf_Optimal_Observer_Examples/PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DNDN_EdgeObserve.m:112-113`
(DJ, 006d40eb, 2025-01-05) – the code senses y1 = int x ds2 + w2 and y2 = int x ds1 + w3; the file name, the
commented term-based block (lines 138-144) and the table of library case 48 (examples_PDE_library_PIETOOLS.m:
953-954; DJ, 827b9ca1, 2024-08-02: y1 = x(s2=d) + w1, y2 = x(s1=b) + w2) describe edge sensing, which
006d40eb replaced by integrals in the conversion to pde_var. Misleading name and table text: the estimator
is synthesized for integral sensing. Fix: sense the edges as the commented block, or rename the file and
correct the table.
Reproduction: by reading: lines 112-113.

## C100 – CONFIRMED (reading)

HEAD: `PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch6_PIESIM.m:272` title
`...,s_{3})$`; x3 = pde_var(s2,[c,d]) in `PIESIM/examples_pde_library_PIESIM_2D.m` (case 10, line 669);
snippet comment 232 "x3(t=tf,s2)"; panel x-label 270 '$s_{2}$'.
Author: 3e01f297 (DJ, 2025-01-01, "Create PIETOOLS_Code_Illustrations_Ch6_PIESIM.m"). DJ.
Reproduction: by reading. Consequence: misleading text (plot title).

**C100** `PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch6_PIESIM.m:272` (DJ, 3e01f297,
2025-01-01) – the title of the x3 panel reads x_3(t, s_3); x3 is a function of s2
(examples_pde_library_PIESIM_2D.m, case 10: x3 = pde_var(s2,[c,d]); x-label of the panel, line 270: s_2).
Misleading text (plot title). Fix: s_{2}.
Reproduction: by reading: line 272.

## C115 – CONFIRMED (run)

HEAD: `EL/Stability_Examples/PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m`: -Bxr in term (147, D = 1, loc = 0) and
batch (115) formats; help 37-39 and log 68-72 (MMP, 11/19): "Stable for lam < 9.7946", +Bxr "stable only
for lam < -5.7799". GUI file `PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das_GUI.mat` loaded at 187.
Author: claim "Stable for lam<pi^2" and +Bxr: adfc1a00 (2021b), in this file since e21afbef (DJ,
2022-08-20); D = 1: 9bd434ab (SS, 2022-12-04); -Bxr and the new claim: 99f7dff5 (MMP, 11/19/2026); GUI file:
blob 1f232087 unchanged since adfc1a00 (2021b), renamed by DJ in e21afbef. Reported to DJ as CANDIDATES
asks (outcome section).
Reproduction: `$S/C115_repro.m` -> `C115_repro.txt` (measured; `dE_das_alpha.m` uses no PIETOOLS function):
- Coded: term C == -Bxr (D = 1, loc = 0), batch E0 == -Bxr*[0 0 I 0].
- eig(A) = -0.680388 ± 5.86573i, -1.27893, -1.82009.
- -Bxr: alpha(pi^2-1) = -0.71228 (cheb N = 16-48), -0.71225 (fd N = 400); one crossing on [-100,30]:
  9.7946151983 (cheb 48), 9.7946152019 (cheb 64), fd N = 100/200/400: 9.793785/9.794408/9.794563,
  Richardson 9.7946152009. Claim of the file 9.7946: agrees.
- +Bxr: alpha(pi^2-1) = +2.9565924541; crossing -5.7799019945 (cheb 48), -5.7799019808 (cheb 64), fd
  Richardson -5.7799001640. Claim -5.7799: agrees.
- Bxr = 0: crossing 9.8696044009 = pi^2 (fd Richardson 9.8696043997); the claim lam < pi^2 of the source
  is the limit of the uncoupled system.
- GUI file, data.eqns{5-8}: "-1.5368*\partial_sX1(0)", "+0.8871*\partial_sX2(0)", "+1.0656*\partial_sX1(0)",
  "+1.1882*\partial_sX1(0)" – the entries of +Bxr (Bxr = [-1.5368 0;0 0.8871;1.0656 0;1.1882 0]).
- The five report lines below, run verbatim: sa(pi^2-1,+1) = 2.956592, sa(pi^2-1,-1) = -0.712283, limits
  -5.779902 (+Bxr), 9.794615 (-Bxr).
Consequence: misleading text in the source and in the file before 11/19 (lam < pi^2 holds for neither
sign); GUI = 1 loads the +Bxr system, unstable at the default lam, while PDE_t and PDE_b carry -Bxr. The
numbers stated in the file at HEAD are right.

**C115** `EL/Stability_Examples/PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m` (claim "Stable for lam<pi^2" and
+Bxr: 2021b, in this file since e21afbef, DJ, 2022-08-20; -Bxr since 99f7dff5, MMP, 2026-11-19) – measured
outside PIETOOLS (Chebyshev collocation, N = 48 and 64, agreeing to 4e-9; second-order finite differences,
N = 100-400 with Richardson extrapolation, agreeing to 2e-6): the coupled system loses stability at lam =
9.794615 with -Bxr (the file since 11/19/2026: "lam < 9.7946"), at lam = -5.779902 with +Bxr (Das et al.
2018 and the file before 11/19/2026), and at lam = pi^2 with Bxr = 0 (eig(A): largest real part -0.680).
The claim lam < pi^2 of the source holds for neither sign; it is the limit of the uncoupled system, which
the term format of 2022 built (coupling at x(0), zero by the BC, until SS added D = 1 in 9bd434ab,
2022-12-04). The GUI file PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das_GUI.mat (2021b) stores +Bxr (data.eqns{5}:
"... -1.5368*\partial_sX1(0)"), so GUI = 1 shows the system that is unstable at the default lam, while
PDE_t and PDE_b carry -Bxr. Misleading text (the source claim) and a GUI file that differs from the code;
the sign of Bxr needs a decision. Fix: none (sign decision; the GUI file and the copies of C92 follow it).
Reproduction (matrices from the file; Bxr = -E0(:,5:6) at HEAD; s = (x+1)/2):
```matlab
[~,Pb] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das(0,{});  A = Pb.A;  Bpv = Pb.E;  Bxr = -Pb.E0(:,5:6);
N = 48; x = cos(pi*(0:N)'/N); c = [2;ones(N-1,1);2].*(-1).^(0:N)'; D = (c*(1./c)')./(x-x'+eye(N+1));
D = 2*(D-diag(sum(D,2)));  i = 2:N;  m = N-1;  L = D(i,:)*D(:,i);  d0 = D(end,i);
sa = @(lam,sig) max(real(eig([A, sig*kron(Bxr,d0); kron(Bpv,ones(m,1)), kron(eye(2),lam*eye(m)+L)])));
fzero(@(l) sa(l,-1),[9.5 10]), fzero(@(l) sa(l,1),[-6 -5.5]), sa(pi^2-1,1)   % 9.7946, -5.7799, 2.9566
```
