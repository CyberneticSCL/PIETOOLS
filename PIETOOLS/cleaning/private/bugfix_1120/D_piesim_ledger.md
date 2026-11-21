# D_piesim ledger – bug-fix pass 11/20/2026

Agent D_piesim. Items C14, C31 (PIETOOLS_PDE.m), C48–C61 (PIESIM). HEAD 17d68747.
Scratch: `C:\Users\mpeet\AppData\Local\Temp\claude\363c499c6332a503\874c85cb-c18d-4876-bea3-39cdd053e21f\scratchpad\D_piesim\`
– scripts: dp_A1.m (C48, probes), dp_A2.m (C49–C59), dp_A3.m / dp_A4.m (follow-up checks),
dp_regress.m + dp_compare.m (stock examples before/after), mk_c14.pl (copies of PIETOOLS_PDE.m);
outputs: A1_pre.txt, A2_pre.txt, C14_pre.log (before the edits), A1_post.txt, A2_post.txt,
A3_post.txt, A4_post.txt, C14_2d_post.txt, C14_1d_post.txt, CMP_*.txt (after),
syntax_check_out.txt.

Authorship. YP re-uploads the PIESIM files through the GitHub web interface, so `git blame`, also
with `-w -C -C -M`, names her for nearly every line. Each changed line was traced with `git log -S`
over all paths (`':(top,glob)**/PIESIM*.m'`, and the 2-D branch d6e76870 squashed into 73b8813a).
Coordinator clarification (DJ hold = logic DJ wrote, not lines he imported, moved or copied)
applied: C53(a) and C58 trace to DJ feature/bug-fix commits (held); no line changed in PIESIM traces
to DJ; in PIETOOLS_PDE.m, blame names DJ for lines 12-13 (1d9b0634, the 2021b import; content from
c587fd94, 2020a release) and line 252 (73b8813a, a squash merge; content YP 665ab0a8, not changed).

Regime: new-file (no developer notes) for PIESIM_PI_discretize.m, PIESIM_disturbance_check.m,
PIESIM_1Dto2D2Mat_cheb_2D.m, PIESIM_2Dto1D2Mat_cheb_2D.m (class Rename in
`release_new_files_1120/new_files_classified.tsv`); full protocol for the other eight files.
`git ls-files --eol` unchanged for every PIESIM file and PIETOOLS_PDE.m (eol_before.txt =
eol_after.txt); `git diff --numstat` equals `-w --numstat` for each file (trailing blanks dropped by
the Edit tool were restored); parse check: 12 files, 0 problems.

Probes (measured, A1_pre.txt): `subs({sx^2},[0.5 0.25])` returns sym [0.25 0.0625] – the Symbolic
toolbox function `toolbox/symbolic/symbolic/subs.m` (not a class method) calls `subs(sym(f),...)`,
so a 1 x 1 cell is accepted; `which nu`, `which n`, `which ns` are empty.

Regression (measured): PIESIM 1-D library 1–41 (N 10, tf 0.02, dt 1e-3, BDF2), 2-D library 1–32
(N 6, tf 0.01), DDE library default example, PIE input of 1-D examples 1, 7, 39, 41: solutions and
grids isequaln before and after all edits (41/41, 32/32, 5/5). Printed output differs only by the
removed false warning of C56(e) (1-D examples 15, 30) and a warning line number shifted by the
header log (2-D examples 13, 14, 16, 17, 18). 1-D example 32 stops in PIESIM_int_bdf:289 before
and after (known, unrelated).

| Item | Status |
|---|---|
| C14 | FIXED |
| C31 | FIXED (text) |
| C48 | FIXED |
| C49 | FIXED |
| C50 | FIXED |
| C51 | FIXED |
| C52 | FIXED |
| C53 | (a) HELD – DJ logic; (b) FIXED |
| C54 | NOT REPRODUCED |
| C55 | FIXED |
| C56 | FIXED |
| C57 | FIXED |
| C58 | HELD – DJ logic |
| C59 | FIXED |
| C60 | NEEDS DECISION |
| C61 | FIXED |

## Items

**C48** FIXED. `PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_PI_discretize.m:61`
(YP, content 7688cdfc 2024-07-10; new-file regime). Before (measured, A1_pre): unit call
`PIESIM_PI_discretize(eta,1)` equals the result for `s1_dum` (sampled along theta), likewise `theta`;
heat equation x_t = x_ss, x(0)=0, x_s(1)=0 (kernel R2 = -s, one variable) declared in `eta` or
`theta`: PIESIM stops with "Unrecognized function or variable 'ns'" (T sampled as a rank-one
operator, solution NaN, then the C57 branch); in `s1`, `s`, `xi`: max|x - exact| = 2.180e-06 at
t = 0.2. Change: a one-variable kernel is a function of theta iff its name ends in '_dum'
(`if ~endsWith(var,'_dum')`), the name pie2piesim gives the dummy variable; help text updated.
Two-variable branch unchanged. After (measured, A1_post): eta, theta give the s1 result in the
unit call; DN in eta, theta: max|x - exact| = 2.180e-06, solution identical to s1 (difference 0);
DD (two-variable kernels, varname {eta,eta_dum}): 1.247e-06 for eta and s1, identical before and
after. Reference: exact solutions sin(pi s/2)exp(-pi^2 t/4), sin(pi s)exp(-pi^2 t); s1/s1_dum
stock cases isequaln before/after (regression).

**C49** FIXED. `PIESIM/PIESIM_routines/PIESIM_Setup/PIESIM_input_check.m:608` (YP, 36d1d9d4).
Before (measured): heat equation with a distributed control u(s,t) = sin(pi s)(1+pi^2+pi^2 t), no
disturbance, PIE input: stops in PIESIM_int_bdf:231 ("indices on the left side are not
compatible"; nux = 0, the control treated as finite-dimensional). Change: `~isempty(PIE.u_tab)`.
After: PIE input max|x - exact| = 1.33e-09, identical to the PDE-input run (difference 0).
Reference: exact solution sin(pi s)(1+t); the PDE route (u counted from PDE.u_tab).

**C50** FIXED. `PIESIM/PIESIM_routines/PIESIM_Setup/PIESIM_initial_setup.m:187-190, 252-255` (YP,
4ffb9738). Before (measured): 2-D example 2 with w1 = exp(-t) given as samples [t; w]: stops in
spline ("The first input must contain unique values"). The corner inputs with derivative order 1
(w_tab rows [1 1 0 0 0 1], [3 1 0 0 0 1]) get corner value 0 for a sampled (spatially constant)
input, and `length(w==1)` is nonzero for any input, so the whole 2 x m array, time row included,
was scaled. Change: `length(uinput.w{i})==1`, and row 2 (the values) scaled; u likewise. After:
samples against the same input as sym exp(-st): max difference 1.7e-14 (max|x| 1.67).

**C51** FIXED. `PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_discretize_icf_2D.m:186`
(YP, c7ded704 2025-06-01). Before (measured): 2-D example 2 with its four boundary inputs declared
as controls: "Unrecognized function or variable 'nu'". Change: `nu=psize.nu;` (one added line).
After: control solution identical to the disturbance solution (difference 0); against the exact
solution 1.13e-01 at N 6 and 2.63e-05 at N 10, equal for both. Not changed (inferred from code,
consistent with the runs): `uinput.uspace(k)`, a 1 x 1 cell, is accepted by subs (probe); the
non-separable fallback of the w branch is not needed, since PIESIM_int_bdf (lines 232-243) and
PIESIM_transform_to_solution_2D (lines 233-245) recompute those coefficients before use.

**C52** FIXED. `PIESIM/PIESIM_routines/Time_Integrators/PIESIM_int_symbolic.m:71, 100, 108, 127, 135`
(YP, 2021b; content 7688cdfc). Before (measured): 1-D example 1, opts.intScheme = 2:
'Unrecognized field name "no"'. Change: line 71 commented out (no is unused; psize.no, the ODE count
of 2025, is now n0); the inputs, handles of t since YP 12/31/2025, evaluated at st
(`uinput.w{i}(st)`, wdot, u, udot). The rerun then stopped in PIESIM_transform_to_solution:173
('Unrecognized field name "w"'): the solution reconstruction reads solcoeff.w, .u (2026), which
only PIESIM_int_bdf returned; two lines added, `solcoeff.w=coeff.w; solcoeff.u=coeff.u;` as in
PIESIM_int_bdf:352-353, and the help lists them. After (A3, A4): 1-D example 1 (N 8, tf 0.05):
max|x - exact| = 3.19e-08, BDF2 dt 1e-5 3.23e-08, difference 5.3e-10; example 7: equal to BDF2 within
1.5e-08; 2-D example 1: equal to BDF2 within 7.7e-05, error to exact 2.2e-01 (N 6) and 6.1e-03 (N 8)
for both. Limitation kept (read): sampled inputs (spline handles) and 2-D non-separable inputs are
not supported by the symbolic integrator.

**C53** (a) HELD. `PIESIM_input_check.m:150` `initialize_PIETOOLS_PDE(psize,true)` was written by DJ
in d6e76870 ("Update PIESIM new terms format.", 2022-10-10). Measured: an ODE-only pde_struct stops
in initialize_PIETOOLS_PDE:75. For Declan: pass PDE (the augmented system). Not edited.
(b) FIXED. `PIESIM_input_check.m:317` (YP, a8d9e617). Before (measured): 2-D example 1 without
uinput.ic: "Unrecognized function or variable 'n'". Change: `n0+ns` (ns = nx+ny+n2, line 233).
After: runs with the default-zero warning; solution isequal to the run with uinput.ic = 0, max|x| = 0.

**C54** NOT REPRODUCED. `PIESIM_discretize_icf.m:148, 166`: subs accepts the 1 x 1 cell (probe,
measured); the 1-D library examples with distributed inputs (e.g. 41) run. No edit.

**C55** FIXED. `PIESIM/PIESIM_routines/Time_Integrators/PIESIM_time_integrate.m:107` (YP, 2021b;
PIESIM_int_bdf took gridall from YP 2/12/2026). Before (measured): double integrator x' = [0 1;0 0]x
plus z_t = z_ss, opts.intScheme = 2: "Not enough input arguments" (PIESIM_int_bdf:116). Change: pass
gridall. After: the fallback message prints; result isequal to BDF2 with the same step (intScheme 1,
dt 1e-3); x(tf) exact ([1.1 2]); max|z - exact| = 3.0e-04 (N 6), 3.5e-05 (N 10).

**C56** FIXED. `PIESIM_input_check.m:530, 1085` (YP, a8d9e617). Before (measured): (c) DDE library
example with opts.ifexact and 3 exact entries: "Unrecognized function or variable 'ns'"; (e)
legacy-terms example 15 (n_pde = [1,0,0,0,1]): the correct count (2 entries) printed the mismatch
warning and turned ifexact off, a wrong count (1 entry) printed none. Change: n0+nx (DDE; ns exists
only in the PDE branch), n0+ns (legacy; ns = sum(psize.n), line 1078). After: (c) runs, prints the
mismatch warning (3 entries against n0+nx); (e) 2 entries: no warning, 1 entry: warning. Regression:
the false warning is gone from 1-D examples 15 and 30.

**C57** FIXED. `PIESIM/PIESIM_routines/Postprocessing/PIESIM_transform_to_solution.m:197, 328` (YP,
3acf78b7 2023; ns renamed nx in 2026). Before (measured): ODE initial value NaN (C55 system, BDF):
"Unrecognized function or variable 'ns'" (line 197); C48 reached the same line. Change: `nx>0`;
`solution.timedep.primary{2}=[]` (the 2025 code had the typo `solution.timepde.pde` for
`solution.timedep.pde`). After: runs; final.primary{2} and timedep.primary{2} empty; no field
timepde. Reference: the final-time branch (`solution.final.primary{2}=[]`).

**C58** HELD. `PIESIM_plot_solution.m:229, 148-149`, `PIESIM_plot_solution_2D.m:220`: the colour
index `colors_y(psize.noo-ii+1,:)`, the label `$z_` and the early return were written by DJ in
fe5c062a ("Update PIESIM plot functions.", 2025-01-01; YP a8d9e617 only renamed noo to no0);
`function figs = ...` and `figs = {}` come from DJ commits squashed in 65fce72c ("Change PIESIM
plot_solution", "Update PIESIM_plot_solution.m – Bugfix."). Measured: 8 observed outputs, BDF, plot
'yes': "Index in position 1 exceeds array bounds. Index must not exceed 7." (line 229);
`figs = PIESIM_plot_solution(...)` with plot 'no': 'Output argument "figs" ... not assigned'. For
Declan: `colors_y(ii,:)` or a modulo index, `$y_`, `figs = {}` before the return. Not edited.

**C59** FIXED. `PIESIM_1Dto2D2Mat_cheb_2D.m`, `PIESIM_2Dto1D2Mat_cheb_2D.m` (new-file regime),
`PIESIM_9PI2Mat_cheb_2D.m` (YP: guards 7688cdfc, sums 9c0fea6a; the DJ 12/16/2024 changes of the call
lines untouched). Before (measured, direct calls with an empty kernel): "Unrecognized function or
variable 'A_opmult'" (1Dto2D:108, 2Dto1D:107), "'A_opint_block_x1'" (9PI:142). Change: before the
guards each term is set to a zero matrix of its size (prod(N+1) x (Ncol+1), (Nrow+1) x prod(N+1),
prod(N+1) x prod(N+1)), one deal statement per file. After: the empty-kernel calls run and equal
the zero-kernel calls (A and A_2PDEstate isequal); all 32 2-D library examples isequaln.
Reachability (read): pie2piesim fills zero kernels with zero matrices, so the stock path never had
an empty kernel.

**C60** NEEDS DECISION. `PIESIM/PIESIM_routines/PIESIM_Setup/transl_mod.m`. The 2025 release ended
with `T.I = I;`; YP removed it in 8d717205 ("Fixed rescaling", 2026-05-02), together with `T.I=I` in
transl_2D (6f9b9659; added by her on 4/7/2026) and `PIE.dom=I` in rescalePIE and rescalePIE_2D
(84d6226a, 45066069): a deliberate change of the author, made while PIESIM still used opvar
objects. The help (11/19 pass) says the output is "the operator on I". Measured: with `T.I = I`
added (scratch shadow `c60_shadow/transl_mod.m`), the 41 1-D library examples and the 5 misc cases
give isequaln solutions and identical output. Options: (a) restore `T.I = I;` (2025 behaviour,
matches the help, no measured effect); (b) keep the 5/2/2026 removal and say in the help that T.I
keeps the original interval. No edit.

**C61** FIXED. `PIESIM/PIESIM_routines/PIESIM_Setup/PIESIM_disturbance_check.m:57` (YP, ffef0e70;
new-file regime). Semantics: the assignment creates a variable `nargin` that shadows the function
with the same value; nothing reads it. Change: line deleted. Regression: every run passes this
function; isequaln before/after.

**C14** FIXED. `PIETOOLS_PDE.m:246, 254, 266` (YP, 98d89c7c). Before (measured, scratch copy with a 2-D
PDE – one state in s1, one in (s1,s2) – and plot 'no'): "Index in position 1 exceeds array bounds.
Index must not exceed 1." at `solution.timedep.pde{1,2,3}=solution.timedep.primary{2,3,4}`; the 1-D
path ran. Intended layout (PIESIM.m help): primary is 1 x 4, {2} states of s1 ((N1+1) x nx x
Nsteps), {4} 2-D states ((N1+1) x (N2+1) x n2 x Nsteps); grid is the cell {grid_s1; grid_s2}, so
`length(grid)` = 2 and `grid(:,1)` is a cell. Change: `pde = primary([2,4])` (the two kinds the
comment names; line 252, `length(pde)==2`, then holds), k from grid{2}, one added line
`sgrid` (grid{1} in 2-D, grid(:,1) in 1-D), surf on sgrid. X (line 239) left unused. After
(C14_2d_post.txt): Step 4b completes; 2 surfaces; ZData equal to primary{2}(:,1,:) and to
primary{4}(:,4,1,:) (s2 node 4 of 7, s2 = 0.5), YData equal to grid{1}; 1-D copy unchanged (ZData
equal to primary{2}, YData equal to grid).

**C31** FIXED (text). `PIETOOLS_PDE.m:12-13` (SS, content c587fd94 2020a; blame names DJ 1d9b0634,
the 2021b import). Read: Step 3 runs no executive by default (every call commented out; the scratch
copies ran Step 3 without one). Change: description "declares a linear coupled ODE-PDE (Step 1),
converts it to a PIE (Step 2), runs the executive uncommented in Step 3 (none by default), and
simulates it with PIESIM (Step 4)"; old lines kept as (was).

## Files edited

PIETOOLS_PDE.m; PIESIM/PIESIM_routines/: Discretization_Of_Operators/PIESIM_PI_discretize.m,
PIESIM_Setup/PIESIM_input_check.m, PIESIM_Setup/PIESIM_initial_setup.m,
PIESIM_Setup/PIESIM_disturbance_check.m, 2D/Discretization_Of_Operators/PIESIM_discretize_icf_2D.m,
2D/Discretization_Of_Operators/PIESIM_1Dto2D2Mat_cheb_2D.m,
2D/Discretization_Of_Operators/PIESIM_2Dto1D2Mat_cheb_2D.m,
2D/Discretization_Of_Operators/PIESIM_9PI2Mat_cheb_2D.m, Time_Integrators/PIESIM_int_symbolic.m,
Time_Integrators/PIESIM_time_integrate.m, Postprocessing/PIESIM_transform_to_solution.m.
