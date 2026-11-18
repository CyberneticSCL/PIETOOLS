# Gap study: smallest closing configurations (2026-11-17)

Scratch folder `gapstudy\`: scripts `gs_*.m` (one per case and form), runner `gs_run.sh` (copy of
`litbuild\lb_run.sh`: watchdog 1,800 s, MATLAB slot protocol via `matlab_slots\gapstudy.lock`), logs
`logs\<script>.log` (+ `.log.meta`), results `results\ex_<tag>.mat` (executive) and `results\q_<tag>.mat`
(direct test). No repository file was edited. MATLAB R2025b, MOSEK, size cap 3,000,000 decision
variables (`litbuild\shadow\lpisolve.m`).

## Harness and conventions

- **'P'** = `PIETOOLS_stability_sop(PIE,st,'P')`, **'Q exec'** = the same executive, form 'Q' (its
  direct form), both through `gs_exec.m`: a copy of `pf_exec_row` / `wp4_row` that routes the overrides
  where the code reads them. Reason: `pf_exec_row` writes a field not present in `st.sop` or
  `st.lpi.ineq` into `st` itself, and `st.lpi.ineq` of the tiers holds only dD, dw, wmin, so
  `extra.psatz`, `extra.psatz_rows`, `extra.wcross` would have been dropped silently. `gs_exec` puts the
  term fields in `st.lpi.ineq` (read by the scratch `lpi_ineq_sop` from the program settings), the
  st.sop fields (dw_feas, max_raise_feas, cap_delta) in `st.sop`, and wcross in appdata `GS_WCROSS`,
  read by a gapstudy copy of the shadow `get_lift_degs` (the executive passes no wcross to
  `lpi_ineq_sop`). The executive ignores `st.lpi.ineq.dw` (it passes dw = dw_feas + raise
  explicitly), so the slack weight raise is set as `st.sop.dw_feas` with `max_raise_feas` 0.
- **direct** = `lb_qdirect` (free Q with the shape of T at the Q degrees, T'Q = Q'T, T'Q >= 1e-4 T'T,
  -(A'Q + Q'A) >= 0, light tier, MOSEK) through `gs_q.m`. "in" = optimal with eta < 1e-8, "out" =
  infeasible, "U" = unknown / numerical problems. Q2 = Q degrees {mult 2, int [2 2], out 2, in 2};
  default Q degrees {mult 1, int [1 1], out 1, in 1}.
- Terms (`gs_t.m`): auto (the reader's rule), none, facesS (faces, scalar offset rule), faces0 (offset
  0), faces (per-direction offset 1), faces2 (per-direction offset 2), box, corners, product2;
  suffix _x0 = wcross 0.
- In the 'P' tables "dw" is the number of raises the executive took (info.dw); the slack weight
  raise base is the tier's dw_feas unless the row says otherwise.

## G3 – AGE (item 16), b = 0.9 b* = 2.0724165 (`gs_g3q`, `gs_g3p`, `gs_g3p2`)

PIE: one L2 space in (s1, s2) that owns both registry directions, so wcross 0 removes nothing (the
_x0 rows below have the sizes of their controls, as measured). Known before (LITBUILD): light 'P'
infeasible after one raise (187,984), heavy 'P' certified (751,849), direct facesS in (56,548),
faces0 in (210,340).

| form | storage / Q degrees | terms | offset | dw | wcross | status | ndv | m or eta | t (s) |
|---|---|---|---|---|---|---|---|---|---|
| direct (control) | Q length 2 | auto | rule (= faces0 size) | 0 | rule | in | 210,340 | eta 2.0e-10 | 9 |
| direct | Q length 2 | none | – | 0 | rule | out | 42,148 | eta 1.0 | 3 |
| direct | Q length 2 | box | 1 per dir | 0 | rule | out | 67,060 | eta 1.0 | 5 |
| direct | Q length 2 | product2 | 1 per dir | 0 | rule | out | 70,660 | eta 1.0 | 5 |
| direct | Q length 2 | faces | 1 per dir | 0 | rule | in | 91,972 | eta 2.5e-10 | 5 |
| direct | Q length 2 | facesS | scalar 1 | 0 | rule | in | 56,548 | eta 2.4e-10 | 3 |
| direct | Q length 2 | none / box / facesS | as above | 0 | 0 | out / out / in | 42,148 / 67,060 / 56,548 | – | 2 – 5 |
| P (control) | stock | auto | rule | 1 raise | rule | infeasible | 187,984 | m 2,968 | 21 |
| P | stock | faces0 | 0 | 1 raise | rule | infeasible | 524,944 | m 3,892 | 44 |
| P | stock | facesS | scalar 1 | 1 raise | rule | infeasible | 187,984 | m 2,968 | 18 |
| P | stock | box | 1 per dir | 1 raise | rule | infeasible | 198,352 | m 2,968 | 19 |
| P | stock | corners | 1 per dir | 1 raise | rule | infeasible | 374,608 | m 2,968 | 25 |
| P | stock | auto | rule | dw_feas 3, no raise | rule | infeasible | 751,744 | m 5,104 | 90 |
| P | stock | auto | rule | 1 raise | 0 | infeasible | 187,984 | m 2,968 | 19 |
| P | Rheavy_u | auto | rule | 0 raises | rule | **optimal, certified** | 188,320 | m 2,996 | 9 |
| P | Rheavy_s | auto | rule | 0 raises | rule | **optimal, certified** | 189,216 | m 3,022 | 8 |
| P | Rlight2 | auto | rule | 0 raises | rule | optimal, certified | 608,561 | m 6,477 | 60 |
| P | Rheavy_u | auto | rule | **dw_feas 0 (2-D), no raise** | rule | **optimal, certified** | **104,080** | m 2,084 | 7 |
| P | Rheavy_u | none | – | 1 raise | rule | infeasible | 332,176 | m 5,132 | 66 |

**Smallest closing configuration.** Direct: facesS (faces, scalar offset 1) at dw 0, 56,548 variables
(the known result; the control 'auto' takes faces0 at w = 1, 210,340, 3.7x larger). 'P': light tier
with the storage row Rheavy_u (2-D entry d1 0, d2 2, d3 [2 2;2 2], uniform) and the slack at
dw_feas 0 without raise: certified at 104,080 variables, 7.2x fewer than the known heavy closure
(751,849) and fewer than the failing light control (187,984). Rheavy_s (the same lengths, separable)
also certifies (189,216 at the light slack). The storage length is the knob: with the stock light
storage no slack change certifies (four term sets, dw 3 at 751,744, wcross 0); Rheavy_s differs from
Rlight0 (infeasible, LITBUILD) only in the first row of len_sos{2}, [1 1] -> [2 2].

**Which 'auto' choice mattered.** In the direct form the offset switch: at w = 1 'auto' takes the faces
at offset 0 where the faces at the scalar offset 1 (facesS) or per-direction offset 1 (faces) certify
at 27 % and 44 % of its size; the box family (box, product2) does not certify, so the face family is
needed. In 'P' no slack choice mattered.

## G2 – S18 (suite, library PDE 43 with a delay state), tau = 0.9 tau* = 0.9312534 (`gs_g2q`, `gs_g2q2`, `gs_g2p`, `gs_g2p2`)

Known before: direct at Q length 2, dw 0: every term set out or U (T2, T3, T5); dw 1 facesS in at
440,023 (T4); 'P' certified by Rlight2 at 1,512,873 (light slack; report Table tab:bl-mixed),
Rlight0 infeasible at 0.5 tau*.

| form | storage / Q degrees | terms | offset | dw | wcross | status | ndv | m or eta | t (s) |
|---|---|---|---|---|---|---|---|---|---|
| direct (control) | Q length 2 | auto | rule (= faces0 size) | 0 | rule | U (numerr 2) | 227,510 | eta 5.2e-10, relres 6.7e-4 | 9 |
| direct | **Q2** | auto | rule (= faces0 size) | 0 | rule | in | 940,188 | eta 2.9e-10 | 72 |
| direct | **Q2** | **facesS** | scalar 1 | **0** | rule | **in** | **235,488** | eta 1.0e-10 | 13 |
| direct | Q2 | faces0 | 0 | 0 | rule | in | 940,188 | eta 2.9e-10 | 71 |
| direct | Q length 2 | auto | rule | 1 | rule | in | 665,815 | eta 7.0e-11 | 29 |
| direct | Q length 2 | none | – | 1 | rule | out | 258,135 | eta 1.0 | 16 |
| direct | Q length 2 | faces | 1 per dir | 1 | rule | in | 692,079 | eta 8.8e-11 | 30 |
| direct | Q length 2 | faces2 | 2 per dir | 1 | rule | in | 390,435 | eta 8.6e-10 | 23 |
| direct | Q length 2 | **facesS** | scalar 1 | **1** | **0** | **in** | **374,827** | eta 1.5e-10 | 18 |
| direct | Q length 2 | box | 1 per dir | 1 | 0 | in | 399,604 | eta 8.9e-10 | 18 |
| P (control) | stock | auto | rule | 1 raise | rule | infeasible | 508,228 | m 4,757 | 51 |
| P | Rlight0 | faces0 | 0 | 1 raise | rule | infeasible | 1,272,105 | m 6,112 | 164 |
| P | Rlight0 | facesS | scalar 1 | 1 raise | rule | infeasible | 510,057 | m 4,841 | 59 |
| P | Rlight0 | box | 1 per dir | 1 raise | rule | infeasible | 511,821 | m 4,841 | 51 |
| P | Rlight0 | auto | rule | dw_feas 3, no raise | rule | infeasible | 1,723,689 | m 7,719 | 297 |
| P | stock | faces0 | 0 | 1 raise | rule | infeasible | 1,270,276 | m 6,028 | 158 |
| P | Rheavy_s | auto | rule | 1 raise | rule | unknown | 1,725,385 | m 7,773 | 272 |
| P | **Rlight2** | auto | rule | **dw_feas 0 (2-D), no raise** | rule | **optimal, certified** | **953,001** | m 7,416 | 108 |
| P | Rlight2 | auto | rule | 0 raises (light slack) | 0 | optimal, certified | 1,325,737 | m 9,375 | 279 |
| P (known closure, control) | Rlight2 | auto | rule | 0 raises (light slack) | rule | optimal, certified | 1,512,873 | m 9,661 | 293 |
| direct | Q2 | none | – | 0 | rule | out | 188,508 | eta 1.0 | 23 |
| direct | Q2 | box | 1 per dir | 0 | rule | in | 282,468 | eta 8.8e-10 | 21 |
| direct | Q2 | faces | 1 per dir | 0 | rule | in | 376,428 | eta 3.4e-11 | 20 |
| direct | **Q2** | **facesS** | scalar 1 | **0** | **0** | **in** | **215,188** | eta 1.2e-10 | 22 |
| direct | Q2 | box | 1 per dir | 0 | 0 | in | 257,093 | eta 6.9e-10 | 21 |
| direct | Q length 2 | **faces2** | 2 per dir | **1** | **0** | **in** | **336,623** | eta 5.1e-10, relres 2.3e-4 | 18 |

**Smallest closing configuration.** Direct, dw 0: yes – Q degrees 2 (the raise of the free operator)
with facesS and wcross 0, 215,188 variables (facesS with the rule 235,488); at Q length 2 nothing
certifies at dw 0. Direct, dw 1 (Q length 2): faces at per-direction offset 2 with wcross 0, 336,623
(facesS with wcross 0 374,827; the known facesS 440,023; 'auto' 665,815). 'P': Rlight2 with the slack
at dw_feas 0 (2-D) and no raise, 953,001 variables (63 % of the known 1,512,873); Rlight2 with wcross 0
at the light slack 1,325,737. Below Rlight2 nothing closed: on Rlight0 four slack changes (faces0,
facesS, box, dw 3) and Rheavy_s are infeasible or unknown. The combination Rlight2 + dw_feas 0 +
wcross 0 was not run.

**Which 'auto' choice mattered.** The offset switch and the Q degree. At w = 1 'auto' takes the faces at
offset 0: at Q length 2 that is UNKNOWN (numerr 2), at Q degrees 2 it certifies at 940,188, 4x the
facesS size at the same Q degrees. At dw 1 'auto' certifies at 665,815 (51 % above facesS at 440,023, so the rule did not resolve to the scalar offset 1 there; which offset it chose was not read out) against per-direction
offset 2 with wcross 0 at 336,623. The cross weight (wcross 0) cut 9 – 15 % everywhere it was tried and
lost no certificate here.

## G1 – S17 (suite, transport edge state driving a 2-D diffusion), r2 = 0.9 pi^2 (`gs_g1`, `gs_g1b`)

Known before: 'P' uncertifiable under every storage (proven); direct at Q length 2: 0.5 certified by
facesS (108,950); at 0.9 every set out or U at dw 0 (T2, T3, T5) and U (numerr 2) at dw 1 (T4).

| form | storage / Q degrees | terms | offset | dw | wcross | status | ndv | m or eta | t (s) |
|---|---|---|---|---|---|---|---|---|---|
| P (control) | stock | auto | rule | 1 raise | rule | unknown | 663,748 | m 5,070 | 114 |
| direct (control) | Q length 2 | auto | rule (= faces0 size) | 0 | rule | U (numerr 2) | 435,350 | eta 8.6e-9, relres 4.0e-4 | 18 |
| direct | Q2 | auto | rule (= faces0 size) | 0 | rule | in | 1,160,588 | eta 2.4e-9 | 84 |
| direct | **Q2** | **facesS** | scalar 1 | **0** | rule | **in** | **290,588** | eta 1.8e-9, relres 5.5e-5 | 19 |
| direct | Q2 | faces0 | 0 | 0 | rule | in | 1,160,588 | eta 2.4e-9 | 82 |
| direct | Q length 2 | faces2 | 2 per dir | 1 | rule | in | 636,630 | eta 6.5e-10 | 76 |
| direct | Q length 2 | facesS | scalar 1 | 1 | 0 | U (numerr 2) | 676,410 | eta 1.5e-8 | 27 |
| direct | Q length 2 | faces0 | 0 | 1 | 0 | U (numerr 2) | 1,851,450 | eta 7.9e-9 | 145 |
| direct | Q length 2 | facesS | scalar 1 | 2 | rule | refused by the cap | 3,155,350 | – | 10 |
| direct | Q length 2 | facesS | scalar 1 | 2 | 0 | U (numerr 2) | 2,624,950 | eta 5.3e-9 | 307 |
| direct | Q2 | facesS | scalar 1 | 1 | rule | in | 2,103,088 | eta 2.2e-9 | 260 |
| direct | Q2 | none | – | 0 | rule | U (numerr 2, relres 10) | 232,588 | eta 7.0e-10 | 22 |
| direct | Q2 | box | 1 per dir | 0 | rule | U (numerr 2) | 348,588 | eta 4.1e-9 | 25 |
| direct | Q2 | faces | 1 per dir | 0 | rule | in | 464,588 | eta 2.0e-9 | 32 |
| direct | Q2 | product2 | 1 per dir | 0 | rule | in | 363,088 | eta 5.4e-9 | 27 |
| direct | **Q2** | **facesS** | scalar 1 | **0** | **0** | **in** | **267,968** | eta 1.1e-9, relres 5.7e-5 | 25 |
| direct | Q2 | box | 1 per dir | 0 | 0 | in | 320,313 | eta 2.1e-9 | 25 |

**Smallest closing configuration.** 'P': not closed (control unknown after one raise, 663,748; proven
uncertifiable, not laddered). Direct: closed at 0.9 by the raise of the free operator Q to degrees 2
(mult 2, int [2 2], out 2, in 2) with facesS and wcross 0 at dw 0, 267,968 variables (facesS with the
rule 290,588). At Q length 2 no configuration certifies at 0.9: the slack raise dw 1 and dw 2 returns
UNKNOWN (numerr 2) for facesS and faces0, the one Q-length-2 certificate is faces at offset 2 at dw 1
(636,630, 2.4x the Q2 closure).

**Which 'auto' choice mattered.** The degree of the free operator Q (not an 'auto' choice of the slack)
and the offset switch: at Q2, 'auto' resolves to the faces0 size (1,160,588, certified, 4.3x the
facesS_x0 closure). The box quadratics need wcross 0 here (box at Q2: U with the rule, in with wcross 0).

## G4 – FO2r and FO2s (item 9, damped wave with a delayed displacement), tau* = 2.881401 (`gs_g4a`, `gs_g4rq`, `gs_g4rp`, `gs_g4sq`, `gs_g4sp`)

Known before (LITBUILD): FO2r infeasible under every test at every fraction, including stabP at
0.1 tau*; FO2s 'P' infeasible at 0.5 (stock, Rlight0); FO2s direct faces0 reach 0.916, facesS U at 0.9.

**Diagnosis: the undelayed wave in the Riemann form of FO2r** (`gs_sys('FO2r0')`: e_t = v_s,
v_t = 0.1 e_s - 2 v, v(0) = 0, e(pi/2) = 0; 1-D registry). Certified by every configuration tried, so
the FO2r failure is not the form of the model alone; the FO2r ladder proceeds.

| form | storage / Q degrees | terms | dw | status | ndv | m or eta | t (s) |
|---|---|---|---|---|---|---|---|
| P | stock light | auto | 0 raises | optimal, certified | 1,120 | m 114 | 1 |
| P | stock heavy | auto | 0 raises | optimal, certified | 4,084 | m 204 | 0 |
| Q exec | light | auto | 0 raises | optimal, certified | 2,984 | m 431 | 1 |
| direct | Q length 2 | auto / none / facesS / faces0 | 0 | in / in / in / in | 1,308 / 968 / 1,648 / 2,824 | eta 7e-12 – 1.4e-10 | 0 |
| direct | Q length 2 | auto | 1 | in | 3,308 | eta 1.2e-11 | 0 |

**FO2r at 0.1 tau*** (registry s1, s2; spaces L2[s1] x2 and L2[s1,s2]).

| form | storage / Q degrees | terms | offset | dw | wcross | status | ndv | m or eta | t (s) |
|---|---|---|---|---|---|---|---|---|---|
| direct (control) | Q length 2 | auto | rule (= faces0 size) | 0 | rule | out | 376,060 | eta 1.0 | 19 |
| direct | **Q2** | facesS | scalar 1 | **0** | rule | **in** | 351,677 | eta 1.6e-9 | 25 |
| direct | **Q2** | **facesS** | scalar 1 | **0** | **0** | **in** | **302,957** | eta 9.9e-10 | 25 |
| direct | Q2 | faces0 | 0 | 0 | rule | in | 1,404,377 | eta 2.6e-11 | 73 |
| direct | Q length 2 | facesS | scalar 1 | 1 | 0 | out | 563,813 | eta 0.27 | 49 |
| direct | Q length 2 | faces2 | 2 per dir | 1 | 0 | out | 510,657 | eta 0.078 | 40 |
| P (control) | stock | auto | rule | 1 raise | rule | infeasible | 840,208 | m 6,042 | 142 |
| P | Rlight2 | auto | rule | dw_feas 0, no raise | rule | optimal, certified | 1,419,841 | m 8,533 | 76 |
| P | **Rheavy_s** | auto | rule | 0 raises (light slack) | rule | **optimal, certified** | **847,204** | m 6,103 | 46 |
| P | stock | faces0 | 0 | 1 raise | rule | infeasible | 2,099,920 | m 7,617 | 255 |

**FO2s at 0.95 tau*** (the direct faces0 reach of LITBUILD was 0.916).

| form | storage / Q degrees | terms | offset | dw | wcross | status | ndv | m or eta | t (s) |
|---|---|---|---|---|---|---|---|---|---|
| direct (control) | Q length 2 | auto | rule (= faces0 size) | 0 | rule | U (numerr 2) | 299,120 | eta 3.3e-2 | 13 |
| direct | **Q2** | facesS | scalar 1 | **0** | rule | **in** | 293,997 | eta 3.4e-9, relres 2.5e-4 | 24 |
| direct | **Q2** | **facesS** | scalar 1 | **0** | **0** | **in** | **248,497** | eta 1.5e-9, relres 8.0e-5 | 25 |
| direct | Q2 | faces0 | 0 | 0 | rule | in | 1,173,657 | eta 1.3e-9 | 93 |
| direct | Q length 2 | facesS | scalar 1 | 1 | 0 | U (numerr 2) | 433,725 | eta 7.2e-2 | 31 |
| direct | Q length 2 | faces2 | 2 per dir | 1 | 0 | U (numerr 2) | 394,961 | eta 3.5e-3 | 33 |
| P (control) | stock | auto | rule | 1 raise | rule | infeasible | 691,888 | m 5,736 | 80 |
| P | Rlight2 | auto | rule | dw_feas 0, no raise | rule | unknown | 1,206,901 | m 8,634 | 128 |
| P | Rheavy_s | auto | rule | 1 raise | rule | unknown | 2,350,772 | m 9,250 | 393 |
| P | stock | faces0 | 0 | 1 raise | rule | infeasible | 1,729,120 | m 7,224 | 217 |

FO2r follow-up (`gs_g4r2`), direct test, Q2, dw 0:

| fraction | terms | wcross | status | ndv | eta | t (s) |
|---|---|---|---|---|---|---|
| 0.5 tau* | facesS | 0 | in | 302,957 | 4.1e-10 | 27 |
| **0.9 tau*** | **facesS** | **0** | **in** | **302,957** | 5.5e-11 | 27 |
| 0.1 tau* | none | – | in | 281,497 | 7.9e-11 | 37 |
| 0.1 tau* | box | 0 | in | 360,957 | 5.1e-10 | 43 |

and 'P' at 0.1 tau*: Rheavy_s with dw_feas 0 (2-D), no raise: optimal, certified, 532,276 variables
(m 4,528, 22 s).

**Smallest closing configuration.**
- FO2r: the undelayed Riemann wave certifies under every test, so the model form is not the obstacle.
  FO2r is closed at 0.9 tau* in the direct form by Q degrees 2 with facesS and wcross 0 at dw 0,
  302,957 variables (also at 0.5 and 0.1 tau*, same size); at 0.1 tau* Q2 with no weighted term
  ('none') certifies at 281,497 – not tested at 0.9. At Q length 2 nothing certifies (dw 0 'auto' out;
  dw 1 out). 'P': closed only at 0.1 tau* (partial: below 0.5): Rheavy_s with the slack at dw_feas 0,
  532,276 (Rheavy_s at the light slack 847,204; Rlight2 dw_feas 0 1,419,841); 'P' at 0.5 and 0.9 tau*
  not run under these storages.
- FO2s at 0.95 tau* (above the earlier reach 0.916): closed in the direct form by Q degrees 2, facesS,
  wcross 0, dw 0, 248,497 variables. 'P' not closed: stock (control), stock + faces0, Rlight2 at
  dw_feas 0 (unknown, 1,206,901) and Rheavy_s (unknown after one raise, 2,350,772).
- The controls above the margin (1.1 tau*) for the Q2 closures are in the section "Controls" below.

**Which 'auto' choice mattered.** The degree of Q (length 2 -> degrees 2), then the offset switch: at
Q2 'faces0' (what 'auto' resolves to at w = 1) costs 4.6x (FO2r) and 4.7x (FO2s) the facesS_x0 closure.
In 'P' (FO2r at 0.1) the storage length is the knob (Rheavy_s, the first row of len_sos{2} 1 -> 2),
as in G3.

## G5 – 3-D heat, DD x DN x DN on the unit cube, kappa = 0.9 kappa*, kappa* = 1.5 pi^2 (`gs_g5q1`, `gs_g5q2`, `gs_g5p1`, `gs_g5p2`)

Note on the case: the brief names kappa* = 3 pi^2 (Dirichlet on every face); the known results it
cites (light 'P' infeasible, 907,729; direct faces0 2.71M certified in 759 s) are for the C-series
benchmark heatNd_pie(3,0,{'DD','DN','DN'}), kappa* = 1.5 pi^2, and the 2.71M certificate is the
Cor. 35 bench LPI (heatNd_lpi at 0.985 kappa*, `pf_run_c5b`), not the executive. The runs here use the
C-series PDE built with pde_var as `pf_run_p2_fin2` built it (`gs_heat3.m`). One L2 space owns all three
directions, so wcross 0 removes nothing (sizes equal, measured).

| form | storage / Q degrees | terms | offset | dw | wcross | status | ndv | m | t (s) |
|---|---|---|---|---|---|---|---|---|---|
| P (known control, `pf_run_p2_fin2`) | stock light | auto | rule | dw_feas 1, 0 raises | rule | infeasible | 907,729 | 23,547 | 230 |
| Q exec | light (lpi.var.deg) | auto | rule | dw_feas 1 | rule | refused by the cap (built in 70 s) | 4,069,288 | – | 70 |
| Q exec | light | auto | rule | dw_feas 1 | 0 | refused by the cap | 4,069,288 | – | 69 |
| P | stock light | faces0 | 0 | dw_feas 1, no raise, cap_delta 0 | rule | infeasible | 907,729 | 23,547 | 250 |
| P | stock light | faces0 | 0 | dw_feas 2, no raise, cap_delta 0 | rule | infeasible | 1,306,897 | 24,747 | 453 |

The faces0 slack at dw 1 has the size of the 'auto' control (seven Gram blocks of 360 in both), so at
w = 1 in 3-D 'auto' already takes the plain term and the six faces at this size (inferred from equal
block sizes and counts; the offsets were not read out). In the Q exec runs the negativity slack alone
is 3,161,088 variables (seven blocks of 672), so the executive direct form is above the cap at light.

**Smallest closing configuration.** Not closed. Largest configuration tried: 'P', stock light storage,
faces at offset 0, dw_feas 2, cap_delta 0 (1,306,897, infeasible); the executive 'Q' form at light
(4,069,288) was refused by the cap. **Which 'auto' choice mattered:** none measured; the slack weight
degree (dw 1 -> 2) did not close the 'P' gap with the stock storage.

## G6 – D2 (item 8, 2-D delayed heat; 3-D registry), tau = 0.5 tau* = 0.517363 (`gs_g6a`, `gs_g6b`, `gs_g6c`, `gs_g6d`)

Known before: 'P' infeasible at dw 0, no raise taken (stock 1,668,608; Rlight0 1,680,004); direct faces0
47,433,948 (refused).

| form | storage / Q degrees | terms | offset | dw | wcross | status | ndv | m | t (s) |
|---|---|---|---|---|---|---|---|---|---|
| direct | Q length 2 | facesS | scalar 1 | 0 | 0 | refused by the cap (built in 57 s) | 6,813,852 | – | 57 |
| direct | Q length 2 | none | – | 0 | 0 | refused by the cap (built in 43 s) | 5,892,828 | – | 43 |
| P | stock light | auto | rule | dw_feas 1, then 1 raise (max_raise_feas 1), cap_delta 0 | rule | killed by the watchdog (1,805 s) during the raised solve | 1,668,608, then 2,514,256 | – | 1,807 |

## Controls of the Q-degree-2 closures (`gs_chk`)

The same direct configurations above the margin, where the system is unstable (exact references of
LITBUILD and the suite), and no weighted term below it. Every unstable case returns no certificate.

| case | fraction | Q degrees | terms | wcross | status | ndv | eta | t (s) |
|---|---|---|---|---|---|---|---|---|
| FO2r | 1.1 tau* | Q2 | facesS | 0 | out | 302,957 | 1.0 | 46 |
| FO2r | 1.1 tau* | Q2 | none | – | out | 281,497 | 1.0 | 34 |
| FO2s | 1.1 tau* | Q2 | facesS | 0 | U (numerr 2) | 248,497 | 1.2e-2 | 35 |
| S17 | 1.1 r* | Q2 | facesS | 0 | out | 267,968 | 3.9e-2 | 29 |
| S18 | 1.1 tau* | Q2 | facesS | 0 | out | 215,188 | 1.0 | 22 |
| FO2r | 0.9 tau* | Q2 | **none** | – | **in** | **281,497** | 6.1e-11 | 50 |
| FO2s | 0.95 tau* | Q2 | **none** | – | **in** | **235,353** | 2.1e-9 | 34 |

So on FO2r (0.9) and FO2s (0.95) the raise of Q alone closes the gap without any weighted slack term:
none at Q2 is the smallest closure (281,497 and 235,353), 7 % and 5 % below facesS_x0 at Q2. On S17 and
S18, none at Q2 does not certify (S17: U, numerr 2, 232,588; S18: out, 188,508).

## G5 addendum – the executive direct form with the faces at the scalar offset (`gs_g5q3`, `gs_g5q4`; coordinator request)

'Q exec', light, stock, slack terms facesS (six 3-D faces, scalar offset 1), wcross rule, dw base 0
(dw_feas 0 in 3-D, no raise). Realized terms (info.terms of the executive, both runs): plain term
(code 0, offset 0) and the factor rows [1 0 0; 2 0 0; 0 1 0; 0 2 0; 0 0 1; 0 0 2] at offsets
[1 1 1 1 1 1]; weight degree w = [1 1 1]. Q2 = st.lpi.var.deg{3} = 2.

| form | Q degrees | terms | dw | status | ndv | m | slack Gram blocks (P >= eppos T'T / negativity) | t (s) |
|---|---|---|---|---|---|---|---|---|
| Q exec | lpi.var.deg{3} (light) | facesS | 0 | infeasible | 658,536 | 25,741 | [360 64 64 48 48 48 48] / [672 112 112 92 92 92 92] | 229 |
| Q exec | deg{3} = 2 | facesS | 0 | infeasible | 666,797 | 107,857 | the same | 243 |

The program fits at 16 % of the faces0 size (4,069,288), but neither run certifies. With Q at degree 2
the slack blocks keep their sizes (the slack is sized from the operator and w stays [1 1 1]); only the
equality rows grow (m 25,741 -> 107,857). **G5 stays not closed**; the largest configuration solved is
'P' with faces0 at dw_feas 2 (1,306,897).
| P | stock light | auto, realized: plain term + the six faces (codes [0 3 4 5 6 7 8]) at offset 1 | 1 | dw_feas 2 in 3-D, no raise, cap_delta 0 (the program `gs_g6c` raised to) | rule (w = [2 2 2]) | infeasible (MOSEK PRIMAL_INFEASIBLE_CER) | 2,514,256 | 22,907 | 1,439 |

**Smallest closing configuration.** Not closed at 0.5 tau*. The direct test at Q length 2 is above the cap
even without weighted terms and with wcross 0 (5,892,828), so no direct configuration fits under 3M at
these Q degrees. Largest configuration solved: 'P', stock light storage, slack at dw_feas 2 with 'auto'
(faces at offset 1, w = [2 2 2]), 2,514,256 variables, infeasible with a certificate. **Which 'auto'
choice mattered:** none closed; at w = 2 the rule takes the faces at offset 1 (read from info.terms).

## Summary

MATLAB wall time of the study: about 190 min (sum of the `.log.meta` walls), one gapstudy process at a
time. "direct" = the free-Q test of `lb_qdirect`; "Q2" = Q degrees {mult 2, int [2 2], out 2, in 2}.
Sizes are decision variables.

| case | gap at the tiers (light, 'auto') | closed in 'P' by | closed in the direct form by | sizes (closure / failing control) |
|---|---|---|---|---|
| G3 AGE, 0.9 b* | 'P' infeasible after one raise; direct certifies ('auto' = faces0 size) | storage row Rheavy_u (2-D len_sos first row [2 2]) + dw_feas(2) = 0, no raise | facesS, dw 0 (Q length 2) | P 104,080 / 187,984; direct 56,548 / 210,340 |
| G2 S18, 0.9 tau* | 'P' infeasible (stock); direct U (numerr 2) | Rlight2 + dw_feas(2) = 0, no raise (known Rlight2 closure 1,512,873) | dw 0: Q2 + facesS + wcross 0; dw 1 at Q length 2: faces at offset 2 per direction + wcross 0 | P 953,001 / 508,228; direct 215,188 (dw 0), 336,623 (dw 1) / 227,510 |
| G1 S17, 0.9 r* | 'P' unknown (proven uncertifiable); direct U | not closed (not laddered: proven) | Q2 + facesS + wcross 0, dw 0 | direct 267,968 / 435,350 |
| G4 FO2r, 0.9 tau* | 'P' infeasible at 0.1 tau* already; direct out | only at 0.1 tau* (partial): Rheavy_s + dw_feas(2) = 0 | Q2, no weighted term ('none'), dw 0 (also Q2 + facesS + wcross 0) | P 532,276 at 0.1 / 840,208; direct 281,497 (302,957) / 376,060 |
| G4 FO2s, 0.95 tau* | 'P' infeasible; direct U | not closed (Rlight2 dw 0, Rheavy_s: unknown) | Q2, 'none', dw 0 (also Q2 + facesS + wcross 0) | direct 235,353 (248,497) / 299,120 |
| G5 3-D heat, 0.9 kappa* | 'P' infeasible; 'Q exec' 4,069,288 > cap | not closed (faces0 dw 1, dw 2) | not closed: 'Q exec' facesS dw 0 infeasible, also with Q degree 2 | largest solved: P 1,306,897; Q exec 666,797 |
| G6 D2, 0.5 tau* | 'P' infeasible at dw 0 | not closed (dw 2 infeasible, 2,514,256) | not testable: direct >= 5,892,828 > cap | – |

Every Q2 closure was checked above the margin (1.1 of tau* or r*): FO2r, S17, S18 out, FO2s U; no
certificate of an unstable system.

### Which knob closed what (measured)

1. **The degree of the free operator Q in the direct form** (Q length 2 -> Q2) closed G1, G2 at dw 0,
   FO2r and FO2s, where no slack change at Q length 2 did (S17, FO2r, FO2s: dw 1 and dw 2 U or out).
   It is the one knob that closed four gaps. On FO2r and FO2s it closes them without any weighted
   slack term. It did not help the 3-D heat in the executive 'Q' form, where the slack keeps its size
   when Q is raised (`gs_g5q4`).
2. **The 'auto' offset switch.** At w = 1 'auto' takes the faces at offset 0. Wherever a smaller set
   certified, it was the faces at the scalar offset 1 (facesS): 3.7x (G3), 4.0x (G2 at Q2), 4.3x
   (G1 at Q2), 4.6x – 4.7x (FO2r, FO2s at Q2) fewer variables than the 'auto' choice at the same
   Q degrees. At Q length 2 the reverse held on S18 at 0.5 (T2: only faces0 certifies), so the offset 0
   is not dominated in general.
3. **The cross weight** (wcross 0) lost no certificate in this study and cut 8 – 17 % (G1, G2, FO2r,
   FO2s); it changes nothing on PIEs with one space owning every direction (AGE, 3-D heat). In T5 it
   turned one certificate (S17 at 0.5, faces0) into U.
4. **The storage length in 'P'** (2-D len_sos first row 1 -> 2: Rheavy_u, Rheavy_s; or Rlight2) closed
   AGE and FO2r at 0.1 and, with Rlight2, S18; slack changes (terms, dw up to 3, wcross) never closed a
   'P' gap with the stock or Rlight0 storage. With the longer storage the light slack raise was not
   needed: dw_feas(2) = 0 certified at 55 % (AGE), 63 % (S18) and 63 % (FO2r) of the size at the light
   slack.

## Recommendations

**Application-specific overrides** (fields of the settings structure; measured to close the case named):

- Mixed-domain PIEs with a transport or delay state (S17, S18, FO2r, FO2s), direct PIE-to-PDE form:
  Q at degrees 2 (`lpivar_cdopvar` degrees {mult 2, int [2 2], out 2, in 2}; in the executive
  st.lpi.var.deg{2} = 2 – the executive 'Q' route with this field was measured in 3-D only, `gs_g5q4`, and
  did not close G5; in 2-D it was not measured), slack terms facesS (st.lpi.ineq.psatz = 0,
  psatz_rows = the 2N faces, psatz_rows_offset = 1, psatz_rows_perdir = false), wcross 0 (a scratch
  option of `get_lift_degs`, not a repository settings field yet), dw 0.
- 'P' form on 2-D PIEs where the stock light storage fails (AGE, FO2r at low delay): the 2-D storage row
  Rheavy_u or Rheavy_s in st.lpi.pos with st.sop.dw_feas = [0 0 1] (the 2-D entry 0) and no raise.
- S18-type delay PIEs in 'P': the row Rlight2 with st.sop.dw_feas(2) = 0 (953,001 against 1,512,873 at
  the light slack).
- The 3-D cases (G5, G6) have no closure under the cap; no override is recommended.

**Tier changes** (only where every measured case agrees):

- None is recommended outright. Two candidates need one regression pass before a tier change, because a
  counterexample or an unmeasured case exists:
  - st.sop.dw_feas(2) = 0 together with a longer 2-D storage (len_sos{2} first row 2): smaller and
    certified in the three cases measured here; with the stock light storage it is untested here and the
    stock light 2-D executive took its raise in every failing control.
  - Q at degrees 2 with facesS and wcross 0 as the 2-D default of the direct form: certified G1, G2, FO2r,
    FO2s and rejected every unstable control; not measured on S3, S8, S26 or S18 at 0.5 tau* (where at
    Q length 2 only faces0 certified, T2), and wcross 0 has one known UNKNOWN (T5, S17 at 0.5, faces0).
