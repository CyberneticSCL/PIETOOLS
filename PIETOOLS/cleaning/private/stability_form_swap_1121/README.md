# Silent form substitution in the executives – fix of 11/21/2026

Maintainer, 11/21/2026, on "the executive quietly swaps one form for another": "this sounds like a bug to fix".

## The defect

`executives_sopvar/PIETOOLS_stability_sop.m` (HEAD 17d68747) replaced a requested form by another without a message:

- for a PIE in two or more variables, 'Qd' by 'Pd' and 'FQ' by 'FP';
- for a PIE with the constraint K x = 0, 'Pd' by 'P' and 'Qd' by 'Q'.

The label printed in the log was set before the substitution, so it named the requested form. Any request for 'Qd'
or 'FQ' on a PIE in two or more variables – a delay adds a variable – ran a different LPI. In the harness, every
'stability_dual' row maps to 'Qd', so every such row on a PIE in two or more variables ran 'Pd'.

No active example is affected: every 2-D and 3-D example in PIETOOLS_examples/Examples_Library has its dual test
commented out (`%evalin('base','stability_dual = 1;')`), the active dual tests are all one-variable
(Examples_Library/Stability_Examples), and the candidate statements of problem_cases_1120 for case 44
(PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay) and X3 (PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D) name
stability_form 'P', so their dual is 'Pd' with or without the change.

## The change (committed 905b0d26 on ndopvar, 11/21 02:05; not pushed)

Maintainer, 11/21: "there should NEVER be a swap"; "these swaps completely change the nature of the result,
certifying different types of stability"; "these forms of the executives can run on 2D problems"; "A nonzero K
form should error if a dual form of the executive is called." An intermediate version that warned on each swap
(run 3, run3_fa0deda4_warning/) is superseded.

- PIETOOLS_stability_sop: 'Qd' and 'FQ' in two or more variables are posed in the direct form of 'Q' (Q free,
  the storage T'*Q; the self-adjointness equality; each slack sized from its own support):
  - 'Qd' is the LPI of 'Q' on the pair (T',A'): (T')'*Q = T*Q, (A')'*Q + Q'*A' = A*Q + Q'*A';
  - 'FQ' with P >= 0 (no margin) and the decay term 2*epneg*T'*T.
  In one variable both keep the declared storage P (T'*Q = P, T*Q = P): programs unchanged.
- PIETOOLS_Hinf_gain_sop: 'Qd' is posed in every dimension (its code reads the spaces of T; was replaced by 'Pd'
  in two or more variables).
- Dual forms with K x = 0 are errors, independent of dimension: stability 'Pd', 'Qd'; gain 'Pd', 'Qd'; H2 'c',
  'cco' (were replaced by 'P', 'Q', 'o', 'oco'). The dual LPIs have no term for K; on 11/07 the dual stability LPI
  certified 7 of 8 unstable PIEs with K (memory periodic-bcs-nd-converter).
- PIETOOLS_Hinf_gain_sop 'Q' with Tw nonzero is an error (was replaced by 'P', label only): the 'Q' KYP operator
  has no Tw terms (only 'P' carries Tw'*R*B1 and Tw'*R*A), so posing it would drop the boundary disturbance.

## Tests

pc_swapfix.m (this folder), run in a clean worktree at 17d68747 (C:/wt1121s) with the patched executive copied in;
results_SWAP_*.csv, logs/.

### Run 1 – 17d68747 plus the patched executive (results_SWAP_b0/b1/k/a/c.csv, logs/SWAP_*.log)

(results_SWAP_b1.csv, results_SWAP_c.csv and results_SWAP_k.csv in this folder were overwritten by their run-3
copies on fa0deda4 – identical in every column but wall time; the run-1 logs keep the run-1 times.)

One variable, before (b0) and after (b1) – identical in every column but time:

| Case | Form | Variables / rows | Status |
|---|---|---|---|
| 5, PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq | 'Qd' light | 780 / 126 | certified |
| 5 | 'FQ' light, epneg 1 | 780 / 126 | unknown (measure 1.342) |
| 8, PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das | 'Qd' light | 2,180 / 631 | certified |

K (part k): periodic heat x_t = x_ss − x on [0,1], x(0) = x(1), x_s(0) = x_s(1), one K row. 'Pd': warning
"PIE.K is nonempty: form 'Pd' is posed as 'P'.", posed 'P', optimal. 'Qd': the same warning for 'Q', posed 'Q',
optimal.

Two and three variables (part a; case 40 at K = −0.5 in part c):

| Case | Form | Variables / rows | Status |
|---|---|---|---|
| 40, PIETOOLS_PDE_Ex_2D_Heat_Eq_with_ODE (example st: light, dD 0, deg_Q{2} 0) | 'Q' | 8,518 / 844 | certified |
| 40 | 'Qd' (new) | 8,518 / 844 | certified |
| 40 | 'Pd' (what 'Qd' ran before) | 6,057 / 530 | unknown (measure 0.21) |
| 40 | 'FQ' (new), epneg 1 | 13,366 / 908 | certified |
| 40 | 'FP' (what 'FQ' ran before), epneg 1 | 38,601 / 1,163 | unknown (measure 0.65) |
| 40 at K = −0.5 (A + BK = 0.5 I, unstable) | 'Qd', 'FQ', 'Q' | 49,846; 90,118; 49,846 (last pass) | infeasible with certificate |
| 43, PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay (two variables; example st: light, dD 0, cap_inc 'auto', deg_Q{2} [1 2 2]) | 'Qd' (new) | 1,008,600 / 10,390 (62 s) | certified |
| 43 at tau = 1.1 tau* (unstable) | 'Qd' | 1,008,600 / 10,390 | infeasible (measure −1) |
| 43 | 'FQ' (new), epneg 1 | 132,609 / 5,086 | unknown (measure −0.86) |
| 56, PIETOOLS_PDE_Ex_2D_Delayed_Heat_Eq_Fridman (three variables), stripped, dD 0 | 'FQ' (new), epneg 1 | 840,156 / 33,352 (153 s) | infeasible with certificate |
| 55, PIETOOLS_PDE_Ex_3D_Coupled_Heat_Eq (three variables), stripped, dD 0 | 'FQ' (new), epneg 1 | 29,822,880 / 283,303 | above the size gate (2e6), not solved; build 222 s |

Cost. One variable: programs identical. Two or more variables: a 'Qd' or 'FQ' request now builds the direct
program in place of the 'Pd' or 'FP' program – case 40: 8,518 against 6,057 ('Qd' against 'Pd') and 13,366 against
38,601 ('FQ' against 'FP'); case 43 'Qd': 1,008,600, 7.6 times the 132,609 of the example 'Q'; case 56 'FQ': 840,156
against 258,409 for 'FP' (11/20); case 55 'FQ': 29.8M against 114,704 for 'FP' (11/20). Memory was not measured.
The added logic is a few scalar tests per call; other forms and dimensions build the same programs.

### Run 2 – b7ba5299 plus the patched executive (results_SWAP2_*.csv, logs/SWAP2_*.log)

The bug-fix pass committed five changes at 01:17 on 11/21 (4f44a021 to b7ba5299), among them
sopvar/@sopvar/ctranspose.m, lpi_classify_sop.m, exec_sop_settings.m, poscopvar.m and copquadvar.m. Every row of
run 1 except case 55 was repeated on b7ba5299 (parts b0, b1, k, c, a2): identical in every column but wall time,
and b0 (before the change) equals b1 (after). The executive file tested is the one in the main checkout.

## PIETOOLS_Hinf_gain_sop and PIETOOLS_H2_norm_sop (maintainer 11/21: "Yes, fix both")

Final design as listed under "The change": gain 'Qd' posed in every dimension; dual forms with K and gain 'Q'
with Tw nonzero are errors (identifiers PIETOOLS_Hinf_gain_sop:formK, :formTw, PIETOOLS_H2_norm_sop:formK,
PIETOOLS_stability_sop:formK). The programs of every admitted form are unchanged.

Tests: pc_swapfix2.m (this folder), parts g0 (before), g1 (the superseded warning version) and g2 (final), on
fa0deda4 in a clean worktree (C:/wt1121u),
with the stability rows of pc_swapfix.m repeated there (b0, b1, k, c, a2). Systems built in the script with
pde_var: periodic heat x_t = x_ss − x + w, z = int x (PIE.K); x_t = x_ss − x, x(0) = w, x(1) = 0 (Tw nonzero);
Dirichlet heat with w, z (no substitution); 2-D heat x_t = x_s1s1 + x_s2s2 + w on [0,1]^2, Dirichlet, z = int x.

(results below)

### Run 4 – fa0deda4, final design (results_SWAP4_*.csv, results_SWAPSURVEY.csv, logs/SWAP4_g2k_survey.log)

| System | Executive, form | Before (g0, fa0deda4) | After (g2) |
|---|---|---|---|
| periodic heat with w, z (K) | gain 'Qd' | ran 'Q': certified, gam 1 | error :formK |
| periodic heat (K) | gain 'Pd' | ran 'P': certified, gam 1 | error :formK |
| periodic heat (K) | H2 'c' | ran 'o': 0.70710678 | error :formK |
| periodic heat (K) | H2 'cco' | ran 'oco': 0.70710678 | error :formK |
| periodic heat x_t = x_ss − x (K) | stability 'Pd', 'Qd' | ran 'P', 'Q' | error :formK |
| x(0) = w (Tw nonzero) | gain 'Q' | ran 'P': unknown, 1.2090 | error :formTw |
| 2-D heat with w, z | gain 'Qd' | ran 'Pd': unknown, 4,750.7 (759,814 variables) | 'Qd': unknown, 71,504 (339,607 variables) |
| Dirichlet heat with w, z | gain 'Q', 'Qd'; H2 'c' | 0.075765897, 0.075903677, 0.18087301 | identical, same sizes |

Survey (pc_swapsurvey.m): 140 catalog rows of the PDE, DDE and DDF libraries, each with the form its test requests
(named in the example, else the default of PIETOOLS_auto_execute_sop), against the PIE it builds. No row requests a
form its PIE does not admit. Two PIEs carry K – case 61 PIETOOLS_PDE_Ex_Heat_Laplace_Beltrami_Edge (stability 'Q')
and case 65 PIETOOLS_PDE_Ex_2D_Linear_Turing_Pair (stability 'P'), both primal; none has Tw nonzero; no gain 'Qd'
test is on a PIE in two or more variables. Two rows did not build: the catalog names PIETOOLS_PDE_Ex_Timoshenko_Beam_3
(case 16) and PIETOOLS_PDE_Example_Diagne (X2), which do not exist; both are primal stability tests.

## Boundary inputs across the executives (11/21, after 905b0d26)

Maintainer, on gain 'Q' with Tw nonzero: "block this too. Do other executives have this issue with Tw?" Blocked
since 905b0d26 (PIETOOLS_Hinf_gain_sop:formTw). Read from the code of every executive in executives_sopvar: each
one either carries the boundary terms in its LPI or errors before building it. None drops them.

| Executive | Tw (boundary disturbance) | Tu (boundary input) |
|---|---|---|
| PIETOOLS_Hinf_gain_sop | 'P': terms Tw'*R*B1 + B1'*R*Tw and Tw'*R*A (lines 297–301); 'Q': error :formTw; 'Pd', 'Qd': error :Tw | no u (u = 0) |
| PIETOOLS_Hinf_estimator_sop | terms Tw'*P*B1 + B1'*P*Tw and Tw'*P*A (lines 219–222) | no u |
| PIETOOLS_H2_norm_sop | error :Tw (line 129) | no u |
| PIETOOLS_H2_estimator_sop | error :Tw (line 113) | no u |
| PIETOOLS_Hinf_control_sop | error :Tw (line 111) | error, same check |
| PIETOOLS_Hinf_OFcontrol_sop | error :Tw (line 133) | error, same check |
| PIETOOLS_H2_control_sop | error :Tw (line 106) | error :Tu (line 108) |
| PIETOOLS_stability_sop, PIETOOLS_well_posedness_sop | no input w | no input u |

## Found on the way (reported, not changed)

- PIETOOLS_PDE_Ex_2D_Heat_Eq_with_ODE (case 40) computes K = k*eye(ne,nx) before it evaluates params, so a params
  string 'k=...' has no effect; the parameter to set is K. The first unstable control (part c with 'k=%g;') solved
  the default system and certified; rerun with 'K=%g;'.
- The same silent substitutions in two other executives: PIETOOLS_Hinf_gain_sop – 'Qd' -> 'Pd' in two or more
  variables, and with K 'Qd' -> 'Q', 'Pd' -> 'P' (no message); 'Q' -> 'P' when Tw ~= 0 (label only, no warning).
  PIETOOLS_H2_norm_sop – with K, 'c' -> 'o' and 'cco' -> 'oco' (no message).
