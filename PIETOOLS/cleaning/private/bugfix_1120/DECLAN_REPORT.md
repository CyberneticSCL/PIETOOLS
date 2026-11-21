# Defects in code by Declan Jagt, and in PIESOS – ndopvar, 11/21/2026

The release comment pass of 11/19–11/20/2026 noticed these defects by reading; the bug-fix pass of 11/20/2026 left them unpatched, because the lines are by Declan Jagt or lie in PIESOS. For each item below, the location is given on ndopvar at fa0deda4 (line numbers hold at 905b0d26), the author commit was established with `git blame -w -C -C -M` and `git log -S`, and the behaviour was measured in MATLAB R2025b unless the entry says "by reading". Each entry gives the author commit, what happens, the consequence and how a user reaches it, a fix where the code, its help or the user manual determines one, and a minimal reproduction. Item numbers are those of `cleaning/private/release_comment_pass/CANDIDATES.md`; C117 is new.

| Item | Location | Consequence |
|---|---|---|
| C82 | `PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch5_PIEs.m:136, 147` | silent wrong result |
| C78 | `Stability_Examples/PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat.m:90`, `PIETOOLS_PDE_Ex_Diagne.m:81` | silent wrong result |
| C83 | `Stability_Examples/PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m:120, 150` | silent wrong result; error |
| C62 | `PIESOS/polyopvar/@intop/isvalid.m:291` | silent wrong result, rare path |
| C24 | `PIETOOLS_DDE.m:282` | silent change of later 2-D settings |
| C37 | `converters/convert_PIETOOLS_DDF.m:101, 119, 131` | error, reachable path |
| C80 | `2D/PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D.m:126-127` | error, reachable path |
| C53(a) | `PIESIM/PIESIM_routines/PIESIM_Setup/PIESIM_input_check.m:152` | error, reachable path |
| C18, C19 | `pietools_update.m:126, 170, 234` | error, rare path |
| C20 | `parser/internal/display_PDE.m:553` | error, rare path |
| C58 | `PIESIM/.../Postprocessing/PIESIM_plot_solution.m:148, 229, 314` | error, rare path; misleading legend |
| C65 | `PIESOS/internal/otimes_params.m:656` | error, rare path |
| C84 | `Stability_Examples/PIETOOLS_PDE_Ex_Saba.m:110` | error, rare path |
| C41, C42, C43, C45 | `converters/legacy`, the 2-D conversion route | errors, and a shadowed function, on a route without callers |
| C26, C27, C28, C29, C88, C117 | parser, 2-D examples, tensopvar | dead code, no effect |
| C30, C33, C34, C47, C72, C96, C97, C99, C100, C101 | parser, converters, demos, examples, copvar | misleading message or text |
| C92 | `CL/stability_8.m:65`, the older Das copy | +Bxr: unstable at the default lam |
| C63, C64, C66–C71, C73, C111 | PIESOS, lines of other authors | Section 8 |

Paths abbreviated `Stability_Examples/`, `2D/` and `CL/` lie in `PIETOOLS_examples/Examples_Library/` (`CL/` = `command_line_input_Examples/`).

## 1. Silent wrong results

**C82** `PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch5_PIEs.m:147`, with T at lines 136 and 169 (3b9962b5, 2024-12-30; line 136 ported to pi4_copvar_sop by DB in f2f855d8, 2026-11-06, kernels unchanged) – `PIE3 = piess({T,0,Tu},A,{0,B1},...)` puts B1 = s(2-s) on the input u and gives PIE3 no disturbance (Tw, Bw, Dzw, Dyw have no columns); the 2025 manual (Sec. 5.5.3, p. 57) passes `{B1,0}`. Compared with the PIE that convert returns for the PDE of Sec. 5.4 (snippet lines 88-98), Tu, A, Cz, Cy, Dzu and Dyu agree, Bu differs (Q2 = 2s - s^2 against 0), and so does T: the snippet has R1 = -s, R2 = -s_dum, the converted PIE R1 = -s_dum, R2 = -s, as do the manual formula (Tv)(s) = -int_0^s theta v(theta) dtheta - int_s^1 s v(theta) dtheta and a derivation by hand. The manual code declares the T of the snippet. With `{B1,0}` every operator except T agrees. PIE3 and the closed loops of lines 154-174 describe a system other than the one of the manual, and the snippet runs without error. Fix: `{B1,0}`, and `T = pi4_copvar_sop([],[],[],{[],-s_dum,-s},[0 0;1 1],s,s_dum,dom)`, in the snippet and in the manual code.
```matlab
% after snippet lines 88-98 and 107-147; PIE = the converted PDE of Sec. 5.4
getfield_opvar1d(PIE3.Bw,'dim')                                   % [0 0;1 0]: no w
getfield_opvar1d(PIE3.Bu,'Q2'), getfield_opvar1d(PIE.Bu,'Q2')     % 2*s - s^2, 0
getfield_opvar1d(PIE3.T,'R.R1'), getfield_opvar1d(PIE.T,'R.R1')   % -s, -s_dum
```

**C78** (a) `Stability_Examples/PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat.m:90` (e21afbef, 2022-08-20; the order dates from the function format of the 2021 library, DJ 05/31/2021) – c(s) is computed from lam at line 90, before the params loop (lines 92-98). params {'lam=1;'} gives c(0) = 5.35, the value for the default lam = 4.65; a copy with lam = 1 in the source gives c(0) = 1.7 = 0.7 + lam. An override of lam, the parameter that examples_PDE_library_PIETOOLS lists for this file (case 7), is ignored without a message. Fix: compute c after the params loop, as `CL/stability_7.m` does since 11/20/2026.
```matlab
[~,P] = PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat(0,{'lam=1;'});
pvar s;  double(subs(P.A0,s,0))      % 5.35; expected 0.7 + 1 = 1.7
```
(b) `Stability_Examples/PIETOOLS_PDE_Ex_Diagne.m:81-87, 130-136` (6a5d184e, 2022-08-18) – the first params loop runs before the values labelled "Default Parameters" (lines 91-92), the second after Lm, Mm and K are derived from them (lines 95-125). params {'Vs=2;'}, {'Hs=2;'}, {'g=9.81;'}, {'Cf=2;'}, {'Bs=1;'} and {'k1=1;','k2=2;'} each return the default PDE (A0, A1, B equal), while copies with Vs = 2, or k1 = 1 and k2 = 2, in the source change A0 by up to 0.944 or B by up to 34.75; Mm and K11, which the help names, take effect. Overrides of the labelled defaults are ignored without a message, and the first loop has no effect. Fix: move the first loop below line 92.
```matlab
[~,P0] = PIETOOLS_PDE_Ex_Diagne(0,{});   [~,P1] = PIETOOLS_PDE_Ex_Diagne(0,{'Vs=2;'});
isequal(P0.A0,P1.A0) && isequal(P0.A1,P1.A1) && isequal(P0.B,P1.B)     % true
```

**C83** `Stability_Examples/PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m:120-122, 150` (3e087c84, 2022-07-08, lines 120-122; 6e55cb23, 2022-08-29, the file and line 150) – an unreferenced copy of `PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m` (no file outside cleaning/ calls it). Its term format couples x(0) (loc = 0, no derivative), which the boundary condition sets to zero: its PIE equals the PIE with Bxr = 0 (T and A equal) and differs from the PIE of the batch format of the same file, which couples x_s(0) (A differs by up to 1.34). In the Ex_ file, where SS added the derivative (9bd434ab, 2022-12-04), the two formats give the same PIE. GUI = 1 fails to load `PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das_GUI.mat`, which does not exist. Silent wrong result in the term format and an error for GUI = 1, for a user who calls the file. Fix: delete the copy (the example library calls the Ex_ file), or set D = 1 and the GUI file name of the Ex_ file.
```matlab
[Pt,Pb] = PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das(0,{});
Pt0 = Pt;  Pt0.x{1}.term{2}.C = zeros(4,2);
P1 = convert_PIETOOLS_PDE(Pt);  P0 = convert_PIETOOLS_PDE(Pt0);
P1.A == P0.A            % true: the coupling has no effect
```

**C62** `PIESOS/polyopvar/@intop/isvalid.m:291`, local is_consistent (6fd03c15, 2026-10-15) – from the second pass on, `safe_vars = [safe_vars,find(safe_idcs)]` records positions in the reduced list rem_vars, not variable numbers, so a resolved variable is missed and another is marked resolved. For the degree-9 row below, whose limits form the cycle s2 -> s3 -> s5 -> s6 -> s7 -> s9 -> s2, isvalid returns true and the constructor accepts the row, while the 2-cycle [0,4, 0,3, 0,2] is rejected ("The intervals of integration are incompatible."). Silent wrong result on a rare path: an intmat whose limits depend on themselves passes the check. Measured: no wrong verdict for any row of degree up to 3, nor for 20000 random rows of each degree 4 to 9; inferred from the loop: a missed cycle needs at least 6 variables, and an acyclic row is never rejected. Fix: `rem_vars(safe_idcs)` for `find(safe_idcs)` (rem_vars is reduced at line 292); the test helper `testing_codes/private/PIESOS/Testfolder/intop/remove_badints.m:110` repeats the line.
```matlab
row = [0,10, 0,3, 0,5, 0,1, 0,6, 0,7, 0,9, 0,1, 0,2];
K = intop({polynomial(1)},{row},{'s'},[0,1]);    % accepted; expected: error
isvalid(K)                                        % 1
K2 = intop({polynomial(1)},{[0,4, 0,3, 0,2]},{'s'},[0,1]);   % error, as intended
```

**C24** `PIETOOLS_DDE.m:282` (9abf156e, 2025-01-06) – 9abf156e turned `settings.sosineq_on=0;` into the script variable `sosineq_on = 0;`. No executive reads it: the settings come from lpisettings_sop (line 284), which has no such field, and the container executives reject the field (`executives_sopvar/utility_functions/exec_sop_settings.m:133-140`). The variable stays in the base workspace, where settings_2d_stock reads it (`settings_sopvar/settings_2d_stock.m:144-145`): once PIETOOLS_DDE.m has run, `settings_2d_stock('light').use_sosineq` is 0 instead of 1 for the rest of the session. No effect on the DDE run; a silent change of every 2-D stock settings struct built later in the session. Fix: delete the line.
```matlab
s0 = settings_2d_stock('light');   evalin('base','sosineq_on=0;');    % the statement of line 282
s1 = settings_2d_stock('light');   [s0.use_sosineq, s1.use_sosineq]   % [1 0]
```

## 2. Errors on reachable paths

**C37** `converters/convert_PIETOOLS_DDF.m:101, 119, 131` and `converters/legacy/convert_PIETOOLS_DDE2PIE_legacy.m:108, 145-159` (8e178278, 2024-12-07) – the distributed-delay kernels are substituted and integrated in s1, while the manual (Chapter 4, Table 4.1: Adi(s), B1di(s), ...) and the help of `PIETOOLS_DDE.m` (lines 159-178) and `PIETOOLS_DDF.m` (lines 144, 421-434) declare them as functions of pvar s. A kernel in s is not integrated: convert_PIETOOLS_DDE stops with "Input contains independent variables." at convert_PIETOOLS_DDF:119, and the legacy route stops with the same message at pi4_copvar_sop:99. Error on a reachable path: every DDE, NDS or DDF with a kernel in s, as documented; constant kernels and kernels in s1 convert. Fix: rename the kernel variable s to s1 before the integration, or document s1 in the manual and the help.
```matlab
pvar s s1
DDE = struct('A0',-1,'tau',1);   DDE.Adi{1} = s;
PIE = convert_PIETOOLS_DDE(DDE,'pie');                     % error at convert_PIETOOLS_DDF:119
DDE.Adi{1} = s1;   PIE = convert_PIETOOLS_DDE(DDE,'pie');  % runs
```

**C80** `2D/PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D.m:126-127` (b84cbab2, 2025-01-05) – the two feedback BCs take x2 and x3 at s1 = 1, while s1 ranges over [a,b], parameters of the file, and the other lines use b. params {'b=2;'}: initialize stops ("The spatial position "BC{3}.term{2}.loc" is not appropriate; ..."); {'a=1;','b=2;'}: initialize places the two terms at s1 = 1 = a, where BCs 1-2 set x2 = x3 = 0, and convert stops with convert_PIETOOLS_BVP:sourceRankDeficient, a message about periodic directions. A copy with b for 1 in the two lines initializes and converts in both cases, and equals the file at (a,b) = (0,1). Error on a reachable path: any params with b ~= 1 (the help states "b = 1 is assumed" since 11/19/2026). Fix: subs(x2,s1,b) and subs(x3,s1,b), as in the commented term-based block.
```matlab
PDE_t = PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D(0,{'b=2;'});
PDE = initialize(PDE_t);     % error: spatial position "BC{3}.term{2}.loc" is not appropriate
```

**C53(a)** `PIESIM/PIESIM_routines/PIESIM_Setup/PIESIM_input_check.m:152` (d6d2d88e, 2022-10-10) – for a pde_struct without spatial variables, the block of lines 140-155 adds a temporary 1-D state to PDE and then calls `initialize_PIETOOLS_PDE(psize,true)`, psize being struct() with the field dim; PIESIM stops with "The input PDE is not appropriately specified. ..." (initialize_PIETOOLS_PDE:75). With PDE passed instead (scratch copy), the run stops later, at PIESIM_discretize_icf:104 ("Unrecognized function or variable 'acheb_glob'."), in 2021b code that assigns acheb_glob only in the loop over PDE states. Error on a reachable path: PIESIM on an ODE-only pde_struct, the case the block exists for. Fix: PDE for psize at line 152 – necessary, not sufficient.
```matlab
X = pde_var(2);   PDE = initialize(diff(X,'t')==[-1 0;0 -2]*X);
opts = struct('N',8,'tf',0.5,'dt',1e-3,'intScheme',1,'plot','no');
sol = PIESIM(PDE,opts,struct('ic',[1 1]));   % error [initialize_PIETOOLS_PDE:75 < PIESIM_input_check:152]
```

## 3. Errors on rare paths

**C18** `pietools_update.m:169-170` (7ced3647, 2022-12-15; the line itself is from the 2021b install script) – the error for a failed mkdir formats with c, which the script assigns only at line 179; 7ced3647 replaced the earlier definition `c = uigetdir(pwd)` by download_dir. On a copy with stubbed downloads and a failing mkdir, the script stops with "Unrecognized function or variable 'c'." at line 170; with a variable c in the base workspace (the file is a script), the message names that value instead of the folder d. Error on a rare path: mkdir of PIETOOLS/tbxmanager fails while tbxmanager is not installed. Fix: d for c.
Reproduction: a copy of pietools_update.m with every download, file and path call stubbed; `exist('tbxmanager')` returns 0 and mkdir returns false.

**C19** `pietools_update.m:126, 234-235` (7ced3647, 2022-12-15) – default_dir is assigned only in the fresh-install branch (line 126; 7ced3647 moved it there from the top of the script), and the savepath fallback reads it in both branches. On the stubbed copy, with savepath failing and the folder dialog cancelled, the update branch stops at line 234 with "Unrecognized function or variable 'default_dir'.", so the path is not saved (by reading, line 246 is skipped as well, so the current folder stays in the new PIETOOLS folder); the fresh-install branch saves to `<start folder>\pathdef.m`. Error on a rare path: an update by a user who cannot write the default pathdef.m and cancels the dialog. Fix: current_dir (line 54, the start folder in both branches) for default_dir at lines 234-235.
Reproduction: the stubbed copy of C18; savepath() returns 1, uigetdir returns 0.

**C20** `parser/internal/display_PDE.m:553-643`, local construct_term (d295c40f, 2022-07-08; use_Cij 74790c13, 2022-11-29; caller line 275, 423134af, 2026-02-17) – use_Cij is assigned only inside `if isfield(PDE_term,'C') && ~isempty(PDE_term.C)`, and line 275 requests it. A term without C in a PDE marked initialized ends in "Output argument "use_Cij" (and possibly others) not assigned a value in the execution with "display_PDE>construct_term" function." Error on a rare path: initialize gives every term a C (no factor without C in a linear and a nonlinear test PDE), so only a PDE edited by hand and marked initialized reaches it. Fix: `use_Cij = false;` before the if, as in display_term.m:306.
```matlab
pvar s t; x = pde_var(s,[0,1]);
PDE = initialize([diff(x,t)==diff(x,s,2)+s*x; subs(x,s,0)==0; subs(x,s,1)==0],true);
PDE.x{1}.term{1} = rmfield(PDE.x{1}.term{1},'C');   PDE.is_initialized = true;
display_PDE(PDE)      % Output argument "use_Cij" ... not assigned
```

**C58** `PIESIM/PIESIM_routines/Postprocessing/PIESIM_plot_solution.m:148-150, 229, 314, 316, 359` and `PIESIM/PIESIM_routines/2D/Postprocessing/PIESIM_plot_solution_2D.m:220` (d00ccc83, 2024-12-26; 770742f5, 2024-12-28; fe5c062a, 2025-01-01) – the observed outputs are coloured `colors_y(psize.no0-ii+1,:)` from 7 rows: with 8 outputs and opts.plot = 'yes', PIESIM stops with "Index in position 1 exceeds array bounds. Index must not exceed 7." (line 229 via PIESIM.m:285; line 220 of the 2-D file); with fewer, the legend of the figure "Observed Output Evolution" labels them $z_1(t)$, $z_2(t)$, the symbol of the regulated outputs. The 1-D file returns at line 149, before `figs = {}` (line 166), when opts.plot is not 'yes', so a direct call with an output argument errors (PIESIM requests none). Also: `colors_PDE(n+psize.n0,:)` (lines 314, 316) has 7 rows and no guard, so 8 PDE states with opts.plot = 'yes' stop at line 314 (the 2-D file guards this index); line 359 adds psize.n0 to the character code of num2str(ii), which gives wrong titles from ii + n0 = 10 on (masked by line 314). Error on a rare path (opts.plot = 'yes'; the default is 'no'), and a misleading legend. Fix: `colors_y(ii,:)` or a modulo index; `$y_` in the legend; `figs = {}` before the return; the guard of the 2-D file at lines 314, 316; `num2str(ii+psize.n0)` at line 359.
```matlab
pvar t s;   syms sx
x = pde_var('state',1,s,[0,1]);   y = pde_var('sense',8);
PDE = initialize([diff(x,t)==diff(x,s,2); y==(1:8)'*int(x,s,[0,1]); subs(x,s,0)==0; subs(x,s,1)==0]);
opts = struct('N',6,'tf',0.01,'dt',1e-3,'intScheme',1,'plot','yes');
sol = PIESIM(PDE,opts,struct('ic',sin(pi*sx)));   % error [PIESIM_plot_solution:229 < PIESIM:285]
```

**C65** `PIESOS/internal/otimes_params.m:656`, local otimes_quadPolyCell (6fd03c15, 2026-10-15) – the check of B tests `~isa(A,'struct')`. otimes_params(2,S), with S the struct form of a quadPoly (help, lines 10-14), is rejected with "Both parameters must be specified as objects of type 'quadPoly'." because A was converted to a quadPoly at line 650, while otimes_params(S,2) returns 2*S; an invalid B with a struct A passes the check and fails at line 661 ("Dot indexing is not supported for variables of this type." for a char, "Property ns is not gettable" for a polynomial). Error on a rare path: direct calls only (`@intop/mtimes_tensopvar.m:791` passes two structs). Fix: `isa(B,'struct')`.
```matlab
S = struct('C',{{sparse(3)}},'Zs',{{}},'Zt',{{}},'ns',{{}},'nt',{{}},'dim',[1,1]);
C = otimes_params(S,2);  full(C.C{1})    % 6
otimes_params(2,S)                       % error: Both parameters must be ... 'quadPoly'.
otimes_params(S,'abc')                   % error at line 661: Dot indexing ...
```

**C84** `Stability_Examples/PIETOOLS_PDE_Ex_Saba.m:110` (function format of the 2021 library, DJ 05/31/2021; moved into this file by 6a5d184e, 2022-08-18) – the error branch calls num2str(indx), while the input is index. INDEX 5 and INDEX 4.6 raise "Unrecognized function or variable 'indx'." in place of the intended message, which also lacks a space before "does not exist". Error on a rare path: an INDEX other than 4 and 4.1-4.5. Fix: index for indx, and the space, as in `PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq.m` since 11/20/2026.
```matlab
PIETOOLS_PDE_Ex_Saba(5,0,{})       % Unrecognized function or variable 'indx'.
```

## 4. The legacy 2-D conversion route (C41, C42, C43, C45)

convert_PIETOOLS_PDE_2D and convert_PIETOOLS_PDE_nonlinear_2D in `converters/legacy` have no caller other than each other, and convert_PIETOOLS_PDE converts the same PDEs. Since the opvar classes were retired (11/08/2026), both stop at their first opvar construction, so C42, C43 and C45 cannot be reached at HEAD. Retiring the route removes all four items; the entries record what each line does if the route is kept.

**C41** `converters/legacy/convert_PIETOOLS_PDE_2D.m:583-616`, `convert_PIETOOLS_PDE_nonlinear_2D.m:438-471, 844`, `converters/legacy/2D/Compute_Tmap_PDE_2D_Separable.m:545, 1230`, `converters/legacy/2D/PIETOOLS_FTC_Expansion_2D.m:229-247` (675e5b8e, 2022-08-10; d8918672, 2024-09-07; 0214a3fd, 2024-10-16; 77a4d086, 2026-02-21) – the files construct objects of the classes opvar and opvar2d and call clean_opvar and mat2opvar, none of which exists any more. A pde_struct passed to either converter stops with "Undefined function 'opvar'" ('opvar2d' in 2-D) at PIETOOLS_FTC_Expansion_2D:229 (246) – measured for 1-D and 2-D heat equations and a nonlinear 1-D PDE; a 2-D PDE with a spatially varying boundary input stops instead with "The expansion of the PDE state in terms of the PIE state fails; perhaps the PDE is ill-posed, or not currently supported.", since Compute_Tmap_PDE_2D_Separable:522-524 catches the missing-class error and blames the PDE. The opvar lines of DB (convert_PIETOOLS_PDE_2D.m:927, the two periodic-map files) lie on the same route.
```matlab
pvar t s
x = pde_var('state',1,s,[0,1]);
PDE = initialize([diff(x,t)==diff(x,s,2); subs(x,s,0)==0; subs(x,s,1)==0]);
PIE = convert_PIETOOLS_PDE_2D(PDE,'silent');   % Undefined function 'opvar' [PIETOOLS_FTC_Expansion_2D:229]
PIE = convert_PIETOOLS_PDE(PDE);               % runs
```

**C42** `converters/legacy/2D/Compute_Tmap_PDE_2D_Separable.m:961, 1016, 1409, 1444, 1453` (d8918672, 2024-09-07; line 1409 d7fd1c0e, 2026-02-15) – compute_input_diffs indexes Dop_x and Dop_y, one entry per scalar state variable after the repelem of lines 892-893, with tmp_log, a 4 x 1 flag per row block, so it reads the derivative orders of state variables 1-4, not of the states that w(k) enters; print_reorder_diff_summary reads varname_LHS_t (line 1409; varnames_LHS_t elsewhere), keeps only the last variable in Ldiff (line 1444 assigns instead of appending), and computes ceil(log(dval(jj)/log(10))) (line 1453). Unreachable at HEAD (the file stops at line 214 first). Fix: varnames_LHS_t at line 1409; `Ldiff = [Ldiff,'(',...]` at line 1444; log(dval(jj))/log(10) at line 1453.
Reproduction: by reading – line 1409 reads a variable that no line of the function assigns; a logical mask selects by position, so a 4 x 1 mask on a longer vector selects entries 1-4.

**C43** `converters/legacy/2D/expand_PIETOOLS_PDE_BCs.m:545` (675e5b8e, 2022-08-10) – the integral-removal loop looks up a state term in `PDE.([Robj,'_tab'])`, Robj being left by the term loop of lines 206-220 at the object of the last term of the BC. In a direct call, for the BC int_0^s1 x2(th,0) dth + w = 0 (x2 the second state, w the only input), line 545 runs with Robj = 'w', finds no row, and line 547 stops with "Index exceeds array bounds."; the same BC written w + int(...) = 0 expands. Unreachable at HEAD through the converters (their callers stop first). Fix: PDE.x_tab at line 545.
```matlab
pvar t s1 s2
x1 = pde_var('state',1,[s1;s2],[0,1;0,1]);  x2 = pde_var('state',1,[s1;s2],[0,1;0,1],[0;2]);  w = pde_var('input',1);
PDE = initialize([diff(x1,t)==diff(x1,s1,2)+diff(x1,s2,2); subs(x1,s1,0)==0; subs(x1,s1,1)==0; subs(x1,s2,0)==0; subs(x1,s2,1)==0;
                  diff(x2,t)==diff(x2,s2,2); int(subs(x2,s2,0),s1,0,s1)+w==0; subs(x2,s2,1)==0]);
PDEx = expand_PIETOOLS_PDE_BCs(PDE);   % Index exceeds array bounds. [expand_PIETOOLS_PDE_BCs:547]
```

**C45** `converters/legacy/convert_PIETOOLS_PDE_2D.m:1942`, `converters/legacy/convert_PIETOOLS_PDE_nonlinear_2D.m:1649` (a7ad0239, 2025-06-18) – each file keeps a local print_convert_summary(PIE,obj,tdiff_list), while `converters/print_convert_summary.m` (moved out of convert_PIETOOLS_PDE.m by DB on 11/18/2026) defines the same name with the arguments (PIE,obj,aux_idcs,F_aux). No effect: a local function takes precedence for the calls in its file (lines 905, 676), and both calls are unreachable at HEAD.
Reproduction: by reading – MATLAB resolves a call inside a file to its local function first.

## 5. Dead code and lines without effect

**C26** `parser/@pde_struct/expand_delays.m:547`, local add_delay_state (00a638da, 2022-10-16) – `PDE.BC_tab(2+nvars+tau_var_idx) = 1` is a single, column-major index, so it sets a cell of row mod(k-1,m)+1 – for m >= k the ID in column 1 of row k – not column k of the row appended at lines 545-546, except for a one-row table. No effect: expand_delays re-initializes the PDE (line 168), and initialize rebuilds BC_tab from zeros (initialize.m:530-533, 641-642). Fix: `PDE.BC_tab(end,2+nvars+tau_var_idx) = 1;`, or delete the line.
Reproduction: by reading – a single subscript indexes in column-major order; initialize rebuilds the table.

**C27** `parser/pde_var.m:293-346`, local getPDEvarID (02f99523, 2024-06-23) – the only call (line 258) has been commented out since the commit that added the function, with stateNameGenerator in its place, and no handle to it is created. Dead code. Fix: delete getPDEvarID.
Reproduction: by reading – line 258 is a comment; a local function is visible only within its file.

**C28** `parser/@pde_struct/times.m:56-59`, `parser/@pde_struct/mtimes.m:65-67` (a85b45da, 2026-02-17) and `mtimes.m:84-185` (4e29b50f, 2024-06-29) – the tests isa(...,'opvar') name a class that is not on the path (`which('opvar','-all')` empty, `exist('opvar','class')` = 0), so they are false for every input: the branch mtimes.m:84-185 is never taken, and times.m:58 and mtimes.m:67 act as else. Dead code. Fix: delete the opvar tests and the branch mtimes.m:84-185.

**C29** `parser/@pde_struct/times_pdes.m:180`, local times_term (b8fad485, 2026-02-13) – obj2_size reads `PDE.(obj1){obj2_num}.size`; obj1 and obj2 are both 'x' (lines 164-168 and 174-178 raise an error for any other object), so the value is that of `PDE.(obj2){obj2_num}`, and `if ~strcmp(obj1,obj2)` (line 185) is never true. No effect while only states can be multiplied. Fix: obj2 for obj1 at line 180.
Reproduction: by reading – lines 164-178 admit only 'x' for both factors.

**C88** `2D/Hinf_Gain_Examples/PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DDDD.m:152-161`, `_DDDN.m:152-160`, `_DNDN.m:153-161` (827b9ca1, 2024-08-02) – r_lim and L2_gain ("For comparison, compute the L2-gain") are local variables of a function with outputs [PDE_t,st], and every line of the block ends in ';'. Dead code: the comparison value never reaches the user. Fix: delete the block, or print L2_gain.
Reproduction: by reading – a function discards its local variables on return, and no line of the block displays.

**C117** (new) `PIESOS/polyopvar/@tensopvar/initialize.m:62-64` (6fd03c15, 2026-10-15) – converts factors of the classes opvar and opvar2d, which are not in the repository; isa(...) is false for every object, so the branch never runs. Dead code. Fix: delete the branch.
Reproduction: by reading; `exist('opvar','class')` = 0.

## 6. Misleading messages and text

**C30** `parser/@pde_struct/combine_vars.m:194` (ed59ca7c, 2022-10-17) – `fprintf('\n','No variables have been merged.\n')` prints one newline: the format has no conversion, so the text argument is not printed. The message is lost in every call without a merge and without silent_merge, also from expand_delays (line 190). Fix: `fprintf('\nNo variables have been merged.\n')`.
```matlab
pvar s t;   x = pde_var(s,[0,1]);
PDE = initialize([diff(x,t)==diff(x,s,2); subs(x,s,0)==0; subs(x,s,1)==0],true);
PDE2 = combine_vars(PDE);    % no "No variables have been merged."
```

**C33** `parser/@pde_struct/display_PDE_crude.m:177, 179, 199` (423134af, 2026-02-17) – only line 197 requests use_Cjj, for the first factor of terms 2 to n; the first term and every later factor are not counted. Uninitialized PDEs x_t = s*x_ss (shown with C11(s)) and x_t = x_ss + x*(s*x) (C122(s)) print no hint 'Call "PDE.C{i,j}"', while display_PDE prints it for both once they are initialized; for a product term (x_t = x_ss + (s*x)*x, C121(s)) the hint names PDE.C{i,j}, where display_PDE names PDE.C{i,j,k}. Missing or misleading hint when an unfinished PDE is displayed. Fix: request the third output at lines 177, 179 and 199 and OR it into use_C, as at lines 197 and 208.
```matlab
pvar s t;   x = pde_var(s,[0,1]);
E = diff(x,t)==s*diff(x,s,2)     % shows C11(s), no 'Call "PDE.C{i,j}"' line
```

**C34** `parser/@pde_struct/mpower.m:46-47` (e5b96e56, 2026-02-18) – the check accepts deg = 0 and its message asks for a "real nonnegative integer"; power rejects 0 (power.m:50-51) after a first check whose message also says "nonnegative" (power.m:47-48). x^(-1) gives "Degrees must be specified as real nonnegative integer.", x^0 gives "Degree-0 powers of PDE terms are currently not supported.". The stated admissible degrees include 0, which is refused. Fix: test deg < 1 in mpower, and ask for a positive integer in both messages.
```matlab
pvar s;   x = pde_var(s,[0,1]);
x^(-1)    % "... real nonnegative integer."
x^0       % "Degree-0 powers of PDE terms are currently not supported."
```

**C47** `converters/operator_converters/dpvars2sdvars.m:64` (6d39aee9, 2026-10-05) – the licence title line reads "PIETOOLS - quad_expand", copied from `PIESOS/internal/quad_expand.m`; help prints it. Fix: "PIETOOLS - dpvars2sdvars".

**C72** `PIESOS/polyopvar/@tensopvar/isvalid.m:127-128`, `initialize.m:65-66` (6fd03c15, 2026-10-15) – both tests accept 'sopvar', 'tensopvar' and 'intop' factors, and both messages name only 'sopvar': tensopvar({T}) with a tensopvar T is accepted, tensopvar({5}) errors with "Factors in tensor-PI operator must be of type 'sopvar'." (initialize:66), and isvalid with ops{1,1} = 5 returns "... must be specified as objects of type 'sopvar'.". Reached by tensopvar(cell) with an invalid factor, by assignment to ops, and through `@tensopmat/isvalid.m:160`. Fix: name the three classes in both messages.
```matlab
T = tensopvar(mat2sopvar(1,struct('in',{{'s'}},'out',{{'s'}}),struct('in',[0,1],'out',[0,1])));
Tn = tensopvar({T});     % accepted
Tb = tensopvar({5});     % error: Factors in tensor-PI operator must be of type 'sopvar'.
```

**C96** `PIETOOLS_demos/DEMO9_L2_gain_analysis.m:32` (65fce72c, 2024-12-28) – the licence title reads "PIETOOLS - DEMO8" in DEMO9 (DEMO8 is DEMO8_H2norm.m). Fix: "PIETOOLS - DEMO9".

**C97** `Hinf_Optimal_Controller_Examples/PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_with_Control.m:96` (0e6e7bd4, 2025-01-15) – z = int (1-s) v2 ds has no control term (the first version, 1aefcaa8, 2025-01-13, had z = [u; int 0.5(1-s)^2 v2 ds] with a scalar u). In the PIE, Dzu (D12) = 0 and Dzw = 0; PIETOOLS_Hinf_control_sop returns gam = 2.81e-4 at 'stripped' (status unknown) and 2.56e-4 at 'light' (optimal), and a record of 11/19/2026 gives 1.69e-4 at the settings of the file (veryheavy). Inferred: with D12 = 0 the synthesis problem does not penalize u, and the bound decreases as the settings grow. A misleading example of H-infinity optimal control. Fix: none determined (a formulation decision, e.g. the z of 1aefcaa8).
```matlab
PDE = PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_with_Control(0,{});  PIE = convert(PDE,'pie');
PIE.Dzu == 0                                                          % true
[~,~,gam] = PIETOOLS_Hinf_control_sop(PIE,lpisettings_sop('light'))   % 2.56e-4
```

**C99** `2D/Hinf_Optimal_Observer_Examples/PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DNDN_EdgeObserve.m:112-113` (006d40eb, 2025-01-05) – the code senses y1 = int x ds2 + w2 and y2 = int x ds1 + w3, while the file name, the commented term-based block (lines 138-144) and the table entry of examples_PDE_library_PIETOOLS for this file (case 48, lines 953-954; 827b9ca1: y1 = x(s2=d) + w1, y2 = x(s1=b) + w2) describe edge sensing, which 006d40eb replaced by integrals in the conversion to pde_var. The estimator is synthesized for integral sensing. Fix: sense the edges, as the commented block does, or rename the file and correct the table.
Reproduction: by reading – lines 112-113.

**C100** `PIETOOLS_demos/snippets_from_manual/PIETOOLS_Code_Illustrations_Ch6_PIESIM.m:272` (3e01f297, 2025-01-01) – the title of the x3 panel reads x_3(t, s_3), while x3 is a function of s2 (`PIESIM/examples_pde_library_PIESIM_2D.m`, case 10: x3 = pde_var(s2,[c,d]); the x-label of the panel, line 270, reads s_2). Fix: s_{2}.

**C101** `sopvar/@copvar/subsref.m:104-118`, `subsasgn.m:113-120, 151-152`, `getfield_opvar1d.m:79, 132`, `getfield_opvar2d.m:89, 124`, `setfield_opvar1d.m:72, 98, 104`, `setfield_opvar2d.m:96, 177`, and the same files of @cdopvar (3b92e68b, 2026-11-15; 2f2af6a1, 2026-11-16; the file headers name MMP) – the errors and the warning of the field interface describe the operator through the retired classes: "'P' is a field of an 'opvar', but this copvar is not opvar-like: ...", "The operator does not have the structure of an 'opvar2d', ...", "An 'opvar' has one spatial variable: 'I' must have one row.", the warning copvar:opvarField "... as for an 'opvar' or 'opvar2d'.", and 'dopvar', 'dopvar2d' for a cdopvar. The classes no longer exist; a user meets the subsref and subsasgn texts by indexing (P.Q1, P.R00, P.I = ...), the getfield and setfield texts by direct calls. The warning is given once per session (at the first read and the first assignment). Fix: name the form instead, e.g. "'%s' is a field of the 1-D form R^m x L2^n[a,b] -> R^p x L2^q[a,b] (fields P, Q1, Q2, R), but this copvar does not have that form: ...", identifiers unchanged.
```matlab
pvar t s;   x = pde_var(1);   z = pde_var('state',1,s,[0,1]);
PIE = convert_PIETOOLS_PDE(initialize([diff(x,t)==-x; diff(z,t)==diff(z,s,2); subs(z,s,0)==0; subs(z,s,1)==x]));
T = PIE.T;   v = T.P;   v = T.Q1;     % copvar:opvarField "... as for an 'opvar' or 'opvar2d'", at T.P only
T2 = T([2 1],[2 1]);   v = T2.P;      % copvar:notOpvar1d "'P' is a field of an 'opvar', ..."
Td2 = cdopvar(T2);   v = Td2.R00;     % cdopvar:notOpvar2d "'R00' is a field of a 'dopvar2d', ..."
```

## 7. The sign of Bxr in the copies of the Das example (C92)

**C92** `CL/stability_8.m:65` (791da3fb, 2026-11-08) and `Stability_Examples/PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m:25, 32, 95` (2021b text, in this file since 6e55cb23, 2022-08-29) – these copies couple the ODE through +Bxr x_s(0), the sign printed in the source (Das et al. 2018, Ex. 2). On 11/19/2026 the maintainer took that sign for a typo of the paper, and the library file `Stability_Examples/PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m` uses -Bxr since then (99f7dff5; its GUI file since 11/21/2026). The spectral abscissa, computed by Chebyshev collocation and by finite differences outside PIETOOLS, crosses zero at lam = -5.7799 with +Bxr, at lam = 9.7946 with -Bxr, and at lam = pi^2 without coupling; at the default lam = pi^2 - 1 it is +2.9566 with +Bxr. stability_8 and the batch format of the older copy are therefore unstable at their default lam, and the claim "Stable for lam<pi^2" (older copy, line 32, from the source) holds for neither sign: it is the limit of the uncoupled system. Fix: -Bxr in stability_8 and in the older copy (or delete the copy, C83), as in the library file; the default lam = pi^2 - 1 then lies inside the range lam < 9.7946.
```matlab
[~,Pb] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das(0,{});  A = Pb.A;  Bpv = Pb.E;  Bxr = -Pb.E0(:,5:6);
N = 48; x = cos(pi*(0:N)'/N); c = [2;ones(N-1,1);2].*(-1).^(0:N)'; D = (c*(1./c)')./(x-x'+eye(N+1));
D = 2*(D-diag(sum(D,2)));  i = 2:N;  m = N-1;  L = D(i,:)*D(:,i);  d0 = D(end,i);
sa = @(lam,sig) max(real(eig([A, sig*kron(Bxr,d0); kron(Bpv,ones(m,1)), kron(eye(2),lam*eye(m)+L)])));
fzero(@(l) sa(l,-1),[9.5 10]), fzero(@(l) sa(l,1),[-6 -5.5]), sa(pi^2-1,1)   % 9.7946, -5.7799, 2.9566
```

## 8. PIESOS: lines of other authors

Declan maintains PIESOS; the lines below are by Sachin Shivakumar (SS), AT and CR.

**C63** `PIESOS/polyopvar/@quadPoly/horzcat.m:51`, `vertcat.m:49` (SS, 23edd836, 2026-01-29) – `elseif isnumeric(X)`: X is never assigned, so [P,3], [3,P], [P;3] and [3;P] error with "Unrecognized function or variable 'X'." Error on a reachable path: concatenation with a numeric array, which the help supports (horzcat.m:16). Fix: `isnumeric(varargin{k})`.
```matlab
P = quadPoly([1;2],{[0;1]},{},[1,1],{'s1'},{});
[P,3]     % Unrecognized function or variable 'X'. (horzcat:51)
```

**C64** `PIESOS/custom/filter_Pcell.m:115` (CR, df0172cc, 2026-10-07) – the loop for the second output reads entry.coefficient, so a numeric Pcell{1,1}, which the help supports and the first loop filters (lines 89-91), errors with "Dot indexing is not supported for variables of this type.", with one output as well (the loop runs regardless of nargout). Error on a reachable path. Fix: for a numeric block, test `entry ~= 0`.
```matlab
F = filter_Pcell({[1,1e-9;1e-9,2]},1e-6)   % error at filter_Pcell:115
```

**C66** `PIESOS/polyopvar/@quadPoly/subs.m:73` (AT, 17bb1f7d, 2026-02-25) – `if size(vars) ~= size(vals)` is true only if both dimensions differ. For P = 1 + 2*s2 + 3*s1 + 4*s1*s2: subs(P,{'s1','s2'},[1 2 3]) returns 16 (the third value ignored), subs(P,{'s1','s2'},1) errors at line 305 (index exceeds), and subs(P,{'s1','s2'},[1;2]), with equal counts, errors with "number of vars and vals should be the same". Error on a rare path; a count mismatch in the same orientation passes silently. Fix: `numel(vars) ~= numel(vals)`.
```matlab
P = quadPoly([1;2;3;4],{[0;1],[0;1]},{},[1,1],{'s1','s2'},{});
R = subs(P,{'s1','s2'},[1 2 3]);  full(R.C)   % 16, no error
subs(P,{'s1','s2'},[1;2])                     % error: number of vars and vals should be the same
```

**C67** `PIESOS/polyopvar/@quadPoly/reshape.m:1, 57` (AT, 914b3721, 2026-02-18) – only RESHAPE(A,M,N) is accepted: reshape(P,[4,1]) errors at line 57 with "Not enough input arguments.", while reshape(P,4,1) works. Error on a rare path.
```matlab
P = quadPoly(sparse([1 3;2 4]),{},{},[2,2],{},{});
reshape(P,[4,1])     % Not enough input arguments. (reshape:57)
```

**C68** `PIESOS/polyopvar/@quadPoly/plus.m:51-88` (SS, lines 54-88, 5bba3db0 to e67d118e, 2026-01-16 to 2026-02-15; DJ, lines 38-53, 9c8ac083, 2026-06-07: a double operand is converted to a 1 x 1 quadPoly, and the comment of line 51 asks for a dimension check) – there is no dimension check; the sum takes F.dim, and liftIndex lifts both operands with F.dim. For F = [s1 0; 0 s1], F + 5 returns [5+s1 0; 0 s1] (expected [s1+5 5; 5 s1+5]), F + t1 adds t1 to entry (1,1) only, 5 + F errors at line 80 ("Arrays have incompatible sizes for this operation."), a 3 x 3 G errors in liftIndex, and quadPoly([1;2]) + quadPoly([3 4]) returns dim [2 1] with a 2 x 2 coefficient matrix. Silent wrong result: P + c for a non-scalar quadPoly P with variables and a nonzero double c. Fix: the dimension check that line 51 asks for, with a 1 x 1 operand expanded first.
```matlab
C = sparse(4,2);  C(2,1) = 1;  C(4,2) = 1;
F = quadPoly(C,{[0;1]},{},[2,2],{'s1'},{});      % [s1 0; 0 s1]
quadPoly.quadPoly2polynomial(F + 5)             % [5+s1 0; 0 s1]; expected [s1+5 5; 5 s1+5]
5 + F                                           % error: Arrays have incompatible sizes
```

**C69** `PIESOS/polyopvar/@quadPoly/subsasgn_1.m:1` (SS: subsasgn.m in 0992dfd0, 2026-02-11, renamed in dc7378f0, 2026-02-18) – the file defines `function P = subsasgn(P,S)`, and MATLAB names a method after its file, so the class has a method subsasgn_1 and no subsasgn: `P.dim = [5,5]` succeeds and leaves dim [5 5] with a 2 x 1 coefficient matrix. Dead code. `sopvar/constructors/int2sopvar.m:163` and `sopvar/@sopvar/inv_1D.m:498-501` assign quadPoly properties from outside the class, so an overload named subsasgn would break them. Fix: delete the file.
```matlab
P = quadPoly([1;2],{[0;1]},{},[1,1],{'s1'},{});
any(strcmp(methods('quadPoly'),'subsasgn'))   % 0 (subsasgn_1: 1)
P.dim = [5,5];                                % no error
```

**C70** `PIESOS/polyopvar/@quadPoly/quadPoly.m:98-101` (SS, dc7378f0, 2026-02-18; lines 99-100 e67d118e, 2026-02-15) – `if ~presorted && 0` never runs, so the input presorted has no effect (the help says so): an unsorted basis Zs = {[2;0;1]} is stored as given for presorted 0 and 1, although the property comments (lines 40-43) require sorted exponent vectors. The branch would call clean(obj); quadPoly has no clean method, and `clean` resolves to extras/clean.m of YALMIP when YALMIP is on the path, which errors for a quadPoly. Dead code. Fix: delete lines 97-101.
Reproduction: by reading – `~presorted && 0` is false for every presorted; methods('quadPoly') lists no clean.

**C71** `PIESOS/polyopvar/@quadPoly/var_swap.m:125-126, 232` (SS, 0992dfd0, 2026-02-11) – sPosMap and tPosMap are assigned and never read; rebuildBasisFromDegrees returns posMap = cell(1,k), never filled. Dead code. Fix: drop the second outputs.
Reproduction: by reading – PosMap and posMap occur only at lines 125, 126, 219, 229 and 232.

**C73** `PIESOS/custom/kernel_degree.m:28` (CR, ca05d251, 2026-09-07) – the licence block says "version 3 of the License"; the other 1065 .m files outside cleaning/ say version 2, as does the LICENSE of the repository. Fix: "version 2".

**C111** `PIESOS/custom/DMB.m:110` (CR, ca05d251, 2026-09-07) – the message "Def. 4 product has an inconsistent output component count." cites a definition of a source that neither the file nor `PIESOS/documentation/polyopvar.pdf` names, and the identifiers 'polyopvar_product:*' carry a former file name. The guard fires only if otimes returns a wrong row count. Fix: name the product as the help does (the tensor product of A and B).
