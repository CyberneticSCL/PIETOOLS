# Phase C: default configuration of the direct PIE-to-PDE stability test (2026-11-17)

Scratch folder `phasec\`: `pc_q.m` (copy of `litbuild\lb_qdirect.m` that also records the realized
terms of both inequalities from the `lpi_ineq_sop` info), `pc_lib.m` (catalog PDEs as `wp4_cases`),
scripts `pc_1d`, `pc_suite1`-`5`, `pc_lit1`-`3`, `pc_heat3n`, `pc_heat3_C*`; runner `pc_run.sh` (copy of
`gs_run.sh`: watchdog 1,800 s, slot lock `matlab_slots\phasec.lock`, no forced start); logs
`logs\<script>.log`; results `results\q_<group>.mat`; tables built from the logs by `pc_summ.pl`.
No repository file was edited. MATLAB R2025b, MOSEK, size cap 3,000,000 (`litbuild\shadow\lpisolve.m`).

Direct test: free Q with the shape of T at the Q degrees, T'Q = Q'T, T'Q >= 1e-4 T'T,
-(A'Q + Q'A) >= 0, dw 0, epneg 0, light tier, MOSEK. "in" = optimal with eta < 1e-8, "out" =
infeasible, "U" = otherwise (solver status in brackets), "cap" = refused by the size cap. Cell = verdict,
ndv, wall time of build + solve.

- C0: Q {mult 1, int [1 1], out 1, in 1}, terms 'auto' (today's default).
- C1: Q degree 1, terms facesS (faces, scalar offset 1).
- C2: Q {mult 2, int [2 2], out 2, in 2}, terms facesS.
- C3: C2 with wcross 0.
- C4: Q degree 2, terms 'auto'.

Term labels (last table): codes / off = the term codes and their offsets (0 plain, 1 box product in 1-D,
3..2nv+2 faces); rows / roff = the number of factor rows (POSCOPVAR_DIRECT psatz_rows) and their
offsets; w = the slack weight degree read by get_lift_degs per registry direction. In 1-D facesS is the
two face rows at offset 1, 'auto' the box product at offset 1.

## Tables

### 1-D library controls

| case | C0 | C1 | C2 | C3 | C4 |
|---|---|---|---|---|---|
| PDE1 | in 327 (2 s) | in 412 (1 s) | in 891 (0 s) | in 891 (0 s) | in 746 (0 s) |
| PDE2 | out 1,308 (0 s) | in 1,648 (0 s) | in 3,564 (0 s) | in 3,564 (0 s) | in 2,984 (0 s) |
| PDE3 | U 2,943 (0 s) [unknown] | U 3,708 (0 s) [unknown] | in 8,019 (0 s) | in 8,019 (0 s) | in 6,714 (0 s) |
| PDE5 | out 435 (0 s) | out 520 (0 s) | in 891 (0 s) | in 891 (0 s) | in 746 (0 s) |
| PDE13 | in 1,740 (0 s) | in 2,080 (0 s) | in 3,564 (0 s) | in 3,564 (0 s) | in 2,984 (0 s) |
| PDE17 | in 1,308 (1 s) | in 1,648 (1 s) | in 3,564 (0 s) | in 3,564 (0 s) | in 2,984 (0 s) |
| PDE21 | in 1,308 (0 s) | in 1,648 (0 s) | in 3,564 (0 s) | in 3,564 (0 s) | in 2,984 (0 s) |

### Library 2-D

| case | C0 | C1 | C2 | C3 | C4 |
|---|---|---|---|---|---|
| PDE33 | in 190,900 (6 s) | in 52,660 (2 s) | in 210,361 (15 s) | in 210,361 (13 s) | in 840,121 (47 s) |

### Mixed suite

| case | C0 | C1 | C2 | C3 | C4 |
|---|---|---|---|---|---|
| S6 | out 534,300 (53 s) | out 144,540 (14 s) | out 320,857 (41 s) | out 274,473 (37 s) | out 1,281,097 (194 s) |
| S7 | out 810,044 (66 s) | out 229,100 (30 s) | out 326,633 (35 s) | out 279,977 (30 s) | out 1,295,417 (110 s) |
| S1 | in 522 (3 s) | in 635 (1 s) | in 1,792 (1 s) | in 1,792 (1 s) | in 1,403 (1 s) |
| S3 | in 383,270 (18 s) | in 95,930 (5 s) | U 262,688 (24 s) [optimal] eta 1.8e-08 | U 241,236 (26 s) [optimal] eta 4.3e-08 | in 1,048,988 (101 s) |
| S4 | in 383,270 (17 s) | in 95,930 (7 s) | in 262,688 (29 s) | U 241,236 (28 s) [optimal] eta 4.0e-08 | in 1,048,988 (111 s) |
| S8 | in 492,048 (28 s) | in 123,168 (9 s) | U 320,833 (33 s) [optimal] eta 4.4e-08 | in 274,449 (30 s) | U 1,281,073 (131 s) [optimal] eta 2.7e-08 |
| S9 | out 443,213 (23 s) | out 120,473 (7 s) | out 432,051 (38 s) | out 399,315 (31 s) | out 1,575,303 (127 s) |
| S26 | in 300,395 (14 s) | in 82,859 (7 s) | in 293,514 (20 s) | in 249,298 (16 s) | in 1,167,594 (72 s) |
| S2 | out 443,213 (21 s) | out 120,473 (8 s) | out 432,051 (34 s) | out 399,315 (35 s) | out 1,575,303 (143 s) |
| S5 | in 387,173 (11 s) | in 97,505 (4 s) | U 915,747 (111 s) [optimal] eta 1.0e-06 | U 807,251 (95 s) [optimal] eta 5.2e-08 | in 1,393,143 (126 s) |
| S18_09 | U 227,510 (12 s) [unknown] | out 62,282 (5 s) | in 235,488 (18 s) | in 215,188 (24 s) | in 940,188 (106 s) |
| S17_09 | U 435,350 (17 s) [unknown] | out 108,950 (5 s) | in 290,588 (29 s) | in 267,968 (29 s) | in 1,160,588 (105 s) |

### Literature set

| case | C0 | C1 | C2 | C3 | C4 |
|---|---|---|---|---|---|
| AGE_09 | in 210,340 (6 s) | in 56,548 (3 s) | U 210,361 (13 s) [optimal] eta 1.4e-08 | U 210,361 (13 s) [optimal] eta 1.4e-08 | in 840,121 (50 s) |
| FO2s_095 | U 299,120 (14 s) [unknown] | out 82,148 (6 s) | in 293,997 (23 s) | in 248,497 (23 s) | in 1,173,657 (97 s) |
| LB_09 | out 568,195 (37 s) | out 154,659 (17 s) | out 837,986 (144 s) | out 693,097 (106 s) | out 2,892,194 (356 s) |
| RF_09 | in 383,270 (12 s) | in 95,930 (4 s) | U 262,688 (19 s) [optimal] eta 1.6e-07 | U 241,236 (16 s) [optimal] eta 6.2e-07 | U 1,048,988 (70 s) [optimal] eta 3.1e-07 |
| TL_09 | in 1,152,400 (118 s) | U 288,400 (31 s) [optimal] eta 2.2e-08 | U 841,444 (298 s) [optimal] eta 2.0e-07 | U 841,444 (283 s) [optimal] eta 2.0e-07 | cap 3,360,484 (14 s) |
| XYL_06 | in 36,966 (2 s) | in 10,506 (1 s) | in 26,241 (1 s) | in 10,581 (1 s) | in 104,541 (2 s) |
| NVF_05 | in 18,588 (1 s) | in 5,100 (0 s) | in 11,660 (1 s) | in 4,700 (0 s) | in 46,460 (1 s) |
| TE | in 10,875 (4 s) | in 13,000 (4 s) | in 22,275 (4 s) | in 22,275 (4 s) | in 18,650 (4 s) |
| FO2r_01 | out 376,060 (23 s) | out 102,928 (9 s) | in 351,677 (29 s) | in 302,957 (35 s) | in 1,404,377 (89 s) |

### S6 and S7: the direct form behind the storage form

E1 = Q degree 1, 'auto', dw 1; E2 = Q degree 1, faces (per-direction offset 1), dw 1; E3 = Q degrees 2, faces (per-direction offset 1), dw 1; E4 = Q degree 1, 'auto', dw 1, wcross 0; E5 = Q degrees 2, faces at offset 0, dw 0 (S6 only). C0 repeated.

| case | C0 | E1 | E2 | E3 | E4 | E5 |
|---|---|---|---|---|---|---|
| S6 | out 534,300 (53 s) | out 1,485,740 (202 s) | out 1,609,420 (172 s) | cap 3,601,677 (17 s) | out 1,189,804 (133 s) | out 1,281,097 (159 s) |
| S7 | out 810,044 (66 s) | out 1,921,064 (329 s) | out 2,382,760 (411 s) | cap 3,625,253 (23 s) | out 1,535,576 (212 s) | not run |

### Summary (all cases with a C0 cell)

| candidate | certified / run | median ndv / C0 | max ndv / C0 | lost vs C0 | gained vs C0 |
|---|---|---|---|---|---|
| C0 | 17 / 29 | 1.00 | 1.00 | none | none |
| C1 | 17 / 29 | 0.27 | 1.26 | TL_09 | PDE2 |
| C2 | 18 / 29 | 0.98 | 3.43 | AGE_09, RF_09, TL_09, S3, S8, S5 | PDE2, PDE3, PDE5, FO2s_095, FO2r_01, S18_09, S17_09 |
| C3 | 18 / 29 | 0.90 | 3.43 | AGE_09, RF_09, TL_09, S3, S4, S5 | PDE2, PDE3, PDE5, FO2s_095, FO2r_01, S18_09, S17_09 |
| C4 | 21 / 29 | 2.74 | 5.09 | RF_09, TL_09, S8 | PDE2, PDE3, PDE5, FO2s_095, FO2r_01, S18_09, S17_09 |

### Realized terms (get_lift_degs through lpi_ineq_sop; V = the T'Q inequality, D = the derivative inequality)

| candidate | terms | cases |
|---|---|---|
| C0 | `V:codes[0 1] off[0 1] rows0 roffzeros(1,0) w0 ; D:codes[0 1] off[0 1] rows0 roffzeros(1,0) w1` | PDE1, PDE2, PDE3, PDE17, PDE21 |
| C0 | `V:codes[0 1] off[0 1] rows0 roffzeros(1,0) w1 ; D:codes[0 1] off[0 1] rows0 roffzeros(1,0) w1` | PDE5, PDE13, TE, S1 |
| C0 | `V:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[0 0] ; D:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[1 1]` | PDE33, XYL_06, S7 |
| C0 | `V:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[0 1] ; D:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[1 1]` | AGE_09, LB_09, S9, S2 |
| C0 | `V:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[1 0] ; D:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[1 1]` | FO2s_095, NVF_05, FO2r_01, S6, S26, S18_09 |
| C0 | `V:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[1 1] ; D:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[1 1]` | RF_09, TL_09, S3, S4, S8, S5, S17_09 |
| C1 | `V:codes0 off0 rows2 roff1 w0 ; D:codes0 off0 rows2 roff1 w1` | PDE1, PDE2, PDE3, PDE17, PDE21 |
| C1 | `V:codes0 off0 rows2 roff1 w1 ; D:codes0 off0 rows2 roff1 w1` | PDE5, PDE13, TE, S1 |
| C1 | `V:codes0 off0 rows4 roff1 w[0 0] ; D:codes0 off0 rows4 roff1 w[1 1]` | PDE33, XYL_06, S7 |
| C1 | `V:codes0 off0 rows4 roff1 w[0 1] ; D:codes0 off0 rows4 roff1 w[1 1]` | AGE_09, LB_09, S9, S2 |
| C1 | `V:codes0 off0 rows4 roff1 w[1 0] ; D:codes0 off0 rows4 roff1 w[1 1]` | FO2s_095, NVF_05, FO2r_01, S6, S26, S18_09 |
| C1 | `V:codes0 off0 rows4 roff1 w[1 1] ; D:codes0 off0 rows4 roff1 w[1 1]` | RF_09, TL_09, S3, S4, S8, S5, S17_09 |
| C2 | `V:codes0 off0 rows2 roff1 w1 ; D:codes0 off0 rows2 roff1 w1` | PDE1, PDE2, PDE3, PDE5, PDE13, PDE17, PDE21, TE |
| C2 | `V:codes0 off0 rows2 roff1 w2 ; D:codes0 off0 rows2 roff1 w1` | S1 |
| C2 | `V:codes0 off0 rows4 roff1 w[1 1] ; D:codes0 off0 rows4 roff1 w[1 1]` | PDE33, AGE_09, FO2s_095, RF_09, TL_09, XYL_06, NVF_05, FO2r_01, S6, S7, S3, S4, S8, S26, S18_09, S17_09 |
| C2 | `V:codes0 off0 rows4 roff1 w[1 2] ; D:codes0 off0 rows4 roff1 w[1 1]` | S9, S2 |
| C2 | `V:codes0 off0 rows4 roff1 w[1 2] ; D:codes0 off0 rows4 roff1 w[1 2]` | LB_09 |
| C2 | `V:codes0 off0 rows4 roff1 w[2 2] ; D:codes0 off0 rows4 roff1 w[1 1]` | S5 |
| C3 | `V:codes0 off0 rows2 roff1 w1 ; D:codes0 off0 rows2 roff1 w1` | PDE1, PDE2, PDE3, PDE5, PDE13, PDE17, PDE21, TE |
| C3 | `V:codes0 off0 rows2 roff1 w2 ; D:codes0 off0 rows2 roff1 w1` | S1 |
| C3 | `V:codes0 off0 rows4 roff1 w[1 1] ; D:codes0 off0 rows4 roff1 w[1 1]` | PDE33, AGE_09, FO2s_095, RF_09, TL_09, XYL_06, NVF_05, FO2r_01, S6, S7, S3, S4, S8, S26, S18_09, S17_09 |
| C3 | `V:codes0 off0 rows4 roff1 w[1 2] ; D:codes0 off0 rows4 roff1 w[1 1]` | S9, S2 |
| C3 | `V:codes0 off0 rows4 roff1 w[1 2] ; D:codes0 off0 rows4 roff1 w[1 2]` | LB_09 |
| C3 | `V:codes0 off0 rows4 roff1 w[2 2] ; D:codes0 off0 rows4 roff1 w[1 1]` | S5 |
| C4 | `V:codes[0 1] off[0 1] rows0 roffzeros(1,0) w1 ; D:codes[0 1] off[0 1] rows0 roffzeros(1,0) w1` | PDE1, PDE2, PDE3, PDE5, PDE13, PDE17, PDE21, TE |
| C4 | `V:codes[0 1] off[0 1] rows0 roffzeros(1,0) w2 ; D:codes[0 1] off[0 1] rows0 roffzeros(1,0) w1` | S1 |
| C4 | `V:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[1 1] ; D:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[1 1]` | PDE33, AGE_09, FO2s_095, RF_09, TL_09, XYL_06, NVF_05, FO2r_01, S6, S7, S3, S4, S8, S26, S18_09, S17_09 |
| C4 | `V:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[1 2] ; D:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[1 1]` | S9, S2 |
| C4 | `V:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[1 2] ; D:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[1 2]` | LB_09 |
| C4 | `V:codes[0 3 4 5 6] off[0 1 1 1 1] rows0 roffzeros(1,0) w[2 2] ; D:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[1 1]` | S5 |
| E1 | `V:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[1 1] ; D:codes[0 3 4 5 6] off[0 1 1 1 1] rows0 roffzeros(1,0) w[2 2]` | S7 |
| E1 | `V:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[2 1] ; D:codes[0 3 4 5 6] off[0 1 1 1 1] rows0 roffzeros(1,0) w[2 2]` | S6 |
| E2 | `V:codes0 off0 rows4 roff1 w[1 1] ; D:codes0 off0 rows4 roff1 w[2 2]` | S7 |
| E2 | `V:codes0 off0 rows4 roff1 w[1 2] ; D:codes0 off0 rows4 roff1 w[2 2]` | S9 |
| E2 | `V:codes0 off0 rows4 roff1 w[2 1] ; D:codes0 off0 rows4 roff1 w[2 2]` | S6 |
| E3 | `V:codes0 off0 rows4 roff1 w[2 2] ; D:codes0 off0 rows4 roff1 w[2 2]` | S6, S7 |
| E4 | `V:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[1 1] ; D:codes[0 3 4 5 6] off[0 1 1 1 1] rows0 roffzeros(1,0) w[2 2]` | S7 |
| E4 | `V:codes[0 3 4 5 6] off[0 0 0 0 0] rows0 roffzeros(1,0) w[2 1] ; D:codes[0 3 4 5 6] off[0 1 1 1 1] rows0 roffzeros(1,0) w[2 2]` | S6 |
| E5 | `V:codes0 off0 rows4 roff0 w[1 1] ; D:codes0 off0 rows4 roff0 w[1 1]` | S6 |

## 3-D heat (sizes only)

DD x DN x DN on the unit cube at 0.9 kappa*, kappa* = 1.5 pi^2 (`gs_heat3(0.9)`, the PDE route of
`pf_run_p2_fin2`). The brief gave kappa* = 3 pi^2, which is the all-Dirichlet cube; the C-series and
`heatNd_pie(3,0,{'DD','DN','DN'})` use DD x DN x DN, whose kappa* is 1.5 pi^2, so that system was run.
The "faces0 certifies at 2,713,563 variables in 759 s" record is the Cor. 35 heat benchmark
(`heatNd_lpi`, bench preset, R and Q at the bench caps, 0.985 kappa*), not this direct test. Counted with
`LB_NDV_CAP` 0 (`pc_heat3n`, build only, 1,041 s):

| candidate | ndv | build (s) | terms (V ; D) |
|---|---|---|---|
| C1 | 10,615,240 | 72 | 6 face rows at offset 1, w [1 1 1] ; same |
| C0 | 67,932,136 | 547 | codes [0 3:8] at offset 0, w [1 1 1] ; same |
| C2 | 53,052,461 | 394 | 6 face rows at offset 1, w [1 1 1] ; same |

All three exceed the 3,000,000 cap: not solved. In 3-D the direct test at dw 0 is 3.5x (C1) to 25x (C0)
the cap.

## Raise ladders (from the tables; s67 rows excluded)

A ladder tries its rungs in order and stops at the first "in"; sums run over every rung tried.

| ladder | certified | cases not certified | sum of ndv over all rungs tried | sum of t (s) |
|---|---|---|---|---|
| C0 -> C4 | 24 / 29 | LB_09, S6, S7, S9, S2 | 21404959 | 1836 |
| C1 -> C2 | 23 / 29 | LB_09, TL_09, S6, S7, S9, S2 | 6527070 | 870 |
| C1 -> C2 -> C0 | 24 / 29 | LB_09, S6, S7, S9, S2 | 10478435 | 1188 |
| C1 -> C3 | 23 / 29 | LB_09, TL_09, S6, S7, S9, S2 | 6086529 | 814 |
| C1 -> C3 -> C0 | 24 / 29 | LB_09, S6, S7, S9, S2 | 10037894 | 1132 |
| C0 -> C1 -> C2 | 24 / 29 | LB_09, S6, S7, S9, S2 | 12758058 | 1001 |

## Incidents

- Settings drift. Another session rewrote the repository tier files at 10:17:50
  (`settings_sopvar_*.m`; `lpi.ineq` gained psatz, psatz_offset, wcross) during this campaign. The first
  literature runs (10:21 – 10:30) returned `'psatz' should be 'auto', ...` on every 'auto' candidate and
  were discarded (`logs\stale\`). From 10:31 `pc_q` puts git-HEAD copies of `lpisettings_sop` and the six
  tier files first on the path (`phasec\shadow_settings\`, asserted). Control: PDE 33 C0 and C1 under the
  frozen tiers reproduce pc_1d exactly (190,900 in, 52,660 in); the explicit-term candidates (C1 – C3) of
  the stale AGE, FO2s, LB runs equal the rerun in ndv, status and eta. Runs before 10:17:50 (pc_1d,
  suite1 – 4) read the HEAD tiers; suite5 (started 10:15) reproduces the gap-study sizes and verdicts of
  S18 and S17 in every cell the gap study ran (227,510 U; 235,488 in; 215,188 in; 940,188 in; 435,350 U;
  290,588 in; 267,968 in; 1,160,588 in).
- `pc_s67` hit the 1,800 s watchdog during S9 E2 (1,328,888 variables, solving); the S9 and S2 rows of the
  addition (E1, E2, beyond the coordinator's request) are not run.
- MATLAB wall: 7,680 s in the recorded runs plus 522 s in the discarded runs (137 min of the 180 min budget).

## Reading

29 cases ran under C0 – C4 (7 1-D library, PDE 33, 12 suite, 9 literature; the 3-D heat exceeds the
cap under every candidate). They split in three:

1. **17 certified by C0.** C1 certifies 16 of them; over all 29 cases its size is a median 0.27x C0 (2-D and mixed: 0.25 –
   0.28x; 1-D: 1.19 – 1.26x, where 'auto' is the box product and facesS two rows); it loses TL (optimal,
   eta 2.2e-8, 288,400 against C0 1,152,400). C2 and C3 lose 5 – 6 of them (AGE, RF, TL, S3, S5 and S8 or
   S4), every loss an 'optimal' exit with eta between 1.4e-8 and 1.0e-6, none infeasible.
2. **7 not certified by C0 but certifiable** (PDE 2, 3, 5, FO2s, FO2r, S18, S17; C0 out or U). C1 closes
   PDE 2 only (S18, S17, FO2s, FO2r, PDE 5 out; PDE 3 U). Q degrees 2 close all seven, with facesS (C2,
   C3) at 0.62 – 1.04x the C0 size on the 2-D/mixed ones (2.0 – 2.7x on PDE 2, 3, 5) (S18 235,488 / 215,188; S17 290,588 / 267,968;
   FO2s 293,997 / 248,497; FO2r 351,677 / 302,957) and with 'auto' (C4, = faces0 at Q2) at 2.7 – 4.1x.
3. **5 out under every candidate** (S6, S7, S9, S2 of the suite, LB of the literature set): infeasible
   with eta 0.12 – 0.98 under C0 – C4, i.e. clean infeasibility, not the gate. The storage form 'P'
   certifies all five at the same parameters (baselines_results.md: S6 525,364, S7 528,666, S9 417,895,
   S2 417,895 under stock; LITBUILD: LB stabP stock 426,044).

**S6 and S7 (coordinator addition).** The direct form at the tiers loses both: C0 – C4 out, the slack
raise dw 1 out (E1 'auto' 1,485,740 / 1,921,064; E2 faces 1,609,420 / 2,382,760; E4 'auto' wcross 0
1,189,804 / 1,535,576), Q degrees 2 with faces0 at dw 0 out on S6 (E5, 1,281,097), Q degrees 2 with
faces at dw 1 above the cap (E3, 3,601,677 / 3,625,253). Nothing in the direct form closed them below the
cap; the storage form 'P' closes them (light stock 525,364 and 528,666, and every row storage of the
baselines). The raises leave the infeasibility clean (S6: C0 eta 0.63, E1 0.093, E2 0.20, E4 0.26, E5 0.15;
S7: C0 0.53, E1 0.47, E2 0.37, E4 0.37); these data indicate no direct-form rung that closes them.

## Recommendation

**Default: C1** – free operator Q of degree 1 (mult 1, int [1 1], out 1, in 1), slack terms facesS
(faces at the scalar offset 1), dw 0. Measured: 17 / 29 certified, the same count as C0, at a median
0.27x the C0 variables (maximum 1.26x, the 1-D rows of a few thousand variables); it trades TL (lost, U)
for PDE 2 (gained).

**First raise: raise the free operator degree** (Q to degrees 2, keeping facesS: C2, or C3 with wcross 0).
Measured: it closes 6 of the 7 cases C1 misses that any candidate certifies (PDE 3, PDE 5, FO2s, FO2r,
S18, S17); the ladder C1 -> C2 certifies 23 / 29 with 6.53M variables summed over all rungs tried, against
24 / 29 and 21.40M for C0 -> C4 (870 s against 1,836 s). The alternatives measured worse as the first
step: lowering the term offset (facesS -> faces0 at Q degree 1, i.e. C0) closes only TL among the C1
misses; raising the slack weight (dw 1 at Q degree 1) closed nothing on S6 and S7 (E1, E2, E4) and costs
2.2 – 2.8x the C0 size. A second rung, faces0 at Q degree 1 (C0), recovers TL: C1 -> C2 -> C0
certifies the same 24 / 29 as C0 -> C4 at 10.48M variables and 1,188 s (C1 -> C3 -> C0: 10.04M, 1,132 s).

Limits of the recommendation:
- Q degrees 2 must not be the default rung: as the first rung (C2 / C3) it loses 5 – 6 of the cases C0
  certifies, all at the eta gate (1.4e-8 – 1.0e-6), so it belongs after C1.
- C2 against C3 (wcross 0): equal counts (18 / 29), C3 at 0.90x against 0.98x the C0 median size; C3 loses
  S4 (eta 4.0e-8) where C2 certifies, and certifies S8 where C2 does not (eta 4.4e-8). These data do not
  separate them.
- The direct form does not replace 'P': S6, S7, S9, S2 and LB are infeasible in the direct form under
  every configuration run and certified by 'P' at the light tier. A stability default that runs the
  direct form alone loses these five; one that falls back to 'P' (or runs 'P' first) does not.
- In 3-D the direct test at dw 0 exceeds the 3,000,000 cap under C0, C1 and C2 (10.6M – 67.9M); no 3-D
  verdict was obtained.
