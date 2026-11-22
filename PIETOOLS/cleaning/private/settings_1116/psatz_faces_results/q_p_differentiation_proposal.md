# Differentiating the direct form 'Q' from the storage form 'P' in every tier (proposal, 11/18/2026 01:45)

Maintainer's request: "More broadly, propose differentiations between Q and P form in all the tiers."

## What the two forms are (maintainer's correction of 11/18)
'P': V(x) = <Tx, P Tx>, P >= eppos I from the storage row (Gram, positivity exact), one slack
for T'PA + A'PT <= 0, sized like P plus dw_feas. 'Q': V(x) = <x, T'Q x>, T'Q = Q'T, Q free at
deg_Q, the storage T'Q; two slacks sized from their own support (P >= eppos2 T'T and
A'Q + Q'A <= 0). 'P' is the case Q = P T, so 'Q' contains 'P' of the same tier only when, per
direction and block, deg Q >= deg P + deg T ("the rule" below; sufficient, not necessary:
simpler certificates exist at lower degree).

## Measured so far (DECISIONS D91-D94, ab_sep.m, ab_dde_fix.m)
- 1-D, light: Q at [2 2 2] loses five DDE thresholds that [4 2 2] restores (the P form of the
  old tree certified them after slack raises at 883-1,852 variables; the Q form at [4 2 2]
  certifies without a raise at 1,549-6,196). No slack raise (dD, dw, the ladder) recovers them.
- 2-D: the Q form at light with 'faces1' is 3.6x cheaper than the old storage form on PDE
  33-37 heavy (210,361 against 751,849 variables) and certifies 40, 43, 38 heavy, 13.2 heavy,
  10 light, 31 heavy, 48 heavy that the storage form did not; it loses PDE 39 light to
  conditioning (MOSEK OPTIMAL, eta 2.4e-9, OPT_STATUS 7.14 at 648,900) where the storage form
  is clean at 933,696; Q at the 2-D cap 2 costs 3.7M variables there.
- 3-D: the Q form at light certifies the 3-D heat at 0.9 kappa* (4.07M variables, Sol).
- The slack raises of the Q form flip verdicts on some rows (PDE 43 light, PDE 36 light) and
  never on the DDE rows; a failing bisection step at veryheavy costs four solves.
- Not measured: the library under the form 'P' on the new tree (the harness maps the stability
  executive to 'Q'; the old tree's "light" rows were the former Q form with a storage, not 'P').

## Kernel degrees of the storage rows (deg P) and the rule for deg_Q
| tier | 1-D storage: multiplier, integral kernels | rule deg_Q{1} with deg T = 1 (ODE-PDE, DDE: deg T_int 0) | 2-D storage per direction | rule 2-D cap |
|---|---|---|---|---|
| extreme | 2, none | [3 0 0] | 2 (sos row 1) | 3 |
| stripped, light, custom | 2, 2 | [3 3 3] ([3 2 2] on DDEs) | 2 | 3 |
| heavy | 4, 2 | [5 3 3] ([5 2 2] on DDEs) | 2 (sos 2 2) | 3 |
| veryheavy | 6, 6 | [7 7 7] | 3 (len_com 3) | 4 |

## Proposed differentiation (the P-form knobs stay as they are)
| tier | deg_Q (1-D, 2-D, N-D) | max_raise_Q (1-D, 2-D, N-D) | raise_Q | default form | evidence or status |
|---|---|---|---|---|---|
| extreme | [3 2 0], 1, 1 | [1 1 0] | 'deg' | 'Q' | rule for the multiplier; d2 kept 2 for the DDE kernels; untested |
| stripped | [4 2 2], 1, 1 | [0 1 0] | 'deg' | 'Q' | 1-D: as light (measured on 10 rows); 2-D: one Q step (cap 2) |
| light | [4 2 2], 1, 1 (set 01:28) | [0 1 0] or [1 1 0] | 'deg' | 'Q', 2-D fallback to 'P' | 1-D base certifies the 10 rows without a raise (D94); set 0 if the rerun shows every certificate at dw 0 |
| custom | as light | as light | 'deg' | as light | custom = light again |
| heavy | [5 3 3], 2, 1 | [1 1 0] | 'both' | 'Q', 2-D fallback to 'P' | rule; cost of [5 3 3] against [4 2 2] not measured |
| veryheavy | [7 7 7], 2, 1 | [1 1 0] | 'both' | 'Q' | rule; today's [3 3 0] slack raises cost 3/4 of a failing bisection step |

Differentiations proposed beyond the degree entries:
1. A form ladder in two or more variables: 'Q' first (cheap, more certificates), then the
   storage form 'P' on UNKNOWN with |OPT_STATUS - 1| > 0.3 (PDE 39 light). Implementation: a
   fourth value of raise_Q, 'form', or a field sop.form_ladder = {'Q','P'}; one more branch in
   lpi_solve_loop_sop's caller. Not implemented.
2. A raise cap after two consecutive infeasibility certificates (cost heuristic, DDE17 at
   veryheavy): saves the two largest solves per failing step, changes no measured verdict.
3. The Q form keeps eppos2 T'T as its margin (the least L2 margin) and the P form eppos I: no
   change proposed.
4. The 2-D term entry 'faces1' stays for both forms (4x cheaper than 'auto' at the same verdicts
   on every 2-D case measured).

## What to measure before deciding (in this order, workstation, about 2 h)
1. The library DDE and NDSDDF rows at light [4 2 2] with max_raise_Q 0 against the rerun now
   running (raise count and step codes tell whether any certificate needed a raise).
2. The library 1-D rows under the form 'P' at light and heavy (a 'P' mapping in lb_exec, one
   line) to size the P form against the Q form on the same tree: which is cheaper in 1-D.
3. heavy at [5 3 3] and veryheavy at [7 7 7] on the 1-D rows: cost against today's values.
4. The 2-D fallback on PDE 36, 38, 39 light: 'Q' then 'P' (settings_2d_stock(light) or the
   row) and the cost of the second solve.

- Correction 11/21/2026 (cleaning/private/stability_form_swap_1121/AUDIT.md): on C:/wt1117 (33d30877) PIETOOLS_stability_sop solved form 'P' for a 'Q' request on a PIE in two or more variables, so the HEAD side of D91 and D93 and "the storage form of HEAD" for PDE 39 (933,696 variables) are 'P' results; the comparisons are direct 'Q' against 'P'. History above left as written.
