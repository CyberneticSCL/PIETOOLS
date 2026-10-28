# settings_sopvar: settings of the container executives (MMP, 10/28/2026)

`lpisettings_sop(type, derivative_strictness, simplify, solver)` returns the settings structure that the
executives of `executives_sopvar` read, for `type` = `stripped`, `light` (default), `heavy`, `veryheavy`
or `custom` (a template equal to light). The arguments are those of `lpisettings`.

```matlab
[prog,P,gam,info] = PIETOOLS_Hinf_gain_sop(PIE, lpisettings_sop('heavy'), 'P');
```

## What a tier sets

| tier | storage 1-D | storage 2-D | dw_feas [1-D 2-D] | dw_P | dw_Q | raises feas / other | dual, witness |
|---|---|---|---|---|---|---|---|
| stripped | light | light | [0 1] | 1 | 1 | 1 / 0 | off |
| light | light | light | [0 1] | 2 | 1 | 2 / 1 | on (order 24, 240 frequencies) |
| heavy | heavy | stripped | [1 2] | 2 | 1 | 2 / 1 | on (24, 240) |
| veryheavy | veryheavy | heavy | [1 2] | 2 | 1 | 3 / 2 | on (32, 480) |

- **Storage** `light`, `heavy`, ...: the degrees of the stock `settings_PIETOOLS_<level>` (1-D) and
  `settings_PIETOOLS_<level>_2D` (2-D), read by `poslpivar_settings_sop` and `poslpivar_settings_2d_sop`.
- **dw_feas**: the weight raise of the slack in stability and well-posedness (`PIETOOLS_stability_sop`,
  `PIETOOLS_well_posedness_sop`); **dw_P**: coercive gains and norms, synthesis; **dw_Q**: the Q
  (non-coercive) forms. The raise is added to the rule of `get_lift_degs`.
- **raises**: `max_raise_feas` / `max_raise`, the passes of the degree loop (`lpi_solve_loop_sop`), each
  rebuilding the program with the weight one higher after MOSEK returns UNKNOWN or INFEASIBLE.
- **dual, witness**: `keepdual` (the dual kernel of the negativity constraint) and `witness` (the 1-D
  numerical counterpart: frequency sweep, spectrum or H2 norm, on a Legendre Galerkin model).
- Every tier: `slack` 'like', `margin` 0, `eppos_stab` 1e-2 (stability and well-posedness, every space),
  `eppos = eppos2 = 1e-4` (1-D gains, norms, synthesis) and `settings_2d.eppos = 1e-2` (2-D), as
  `lpisettings` sets them; MOSEK when installed, else SeDuMi; `simplify` off.

## Measured cost and value (10/28/2026, MOSEK 11.0.30, this workstation)

Plants, each with an exact answer:
- `rd`: x_t = x_ss + f pi^2 x on [0,1], Dirichlet; stable exactly when f < 1.
- `wp`: well-posedness of `rd` at f = 0.5.
- `io1`: rd at f = 0.5 with x_t += w, w in L2, z = int x; gain sqrt(sum_{k odd} 8/(k pi)^2 / mu_k^2) = 0.18247908,
  mu_k = (k pi)^2 - pi^2/2.
- `fin`: rd at f = 0.5 with x_t += s(1-s) w, w in R, z = int x; gain G(0) = 0.033310617, H2 norm
  0.052360097 (modal sums).
- `rd2`: x_t = lap x + f 2 pi^2 x on [0,1]^2, Dirichlet; stable exactly when f < 1.
- `fin2`: rd2 at f = 0.5 with x_t += w, w in R, z = int int x; gain G(0) = 0.068606888.

Gap = certified value / exact value - 1. Variables = numel(prog.decvartable).

**1-D** (every solve at most 1.3 s; every stability LPI certifies at f = 0.5, 0.9 and 0.99, and
well-posedness at f = 0.5):

| tier | P-form stability, variables (f 0.5 / 0.9-0.99) | well-posedness | io1 Q | io1 coercive (variables) | fin Q | fin dual | fin H2 |
|---|---|---|---|---|---|---|---|
| stripped | 280 / 568 | 625 | 1.4e-5 | 6.2e-5 (810) | 1.6e-6 | 4.1e-6 | 1.3e-3 |
| light | 280 / 568 | 625 | 1.4e-5 | 8.8e-6 (1426) | 1.6e-6 | 4.1e-6 | 1.3e-3 |
| heavy | 1021 | 2367 | 1.4e-5 | 1.1e-6 (2259) | 1.4e-6 | 2.1e-6 | 1.3e-3 |
| veryheavy | 5396 | 11653 | -3.5e-9 | 8.9e-7 (8922) | 5.0e-6 | 3.7e-6 | 9.4e-6 |

**2-D** (variables, seconds):

| tier | stability f 0.5 | stability f 0.9 | fin2 coercive gain |
|---|---|---|---|
| stripped | 103,744, 3.7 s | 187,984, 9.4 s (one raise) | 5.0e-3 (106,645, 4.0 s) |
| light | 103,744, 2.3 s | 187,984, 9.6 s (one raise) | 3.8e-4 (191,605, 9.9 s) |
| heavy | 751,849, 35 s | 751,849, 40 s | 1.0e-4 (759,358, 52 s) |
| veryheavy | not run | 2,376,580, 475 s | 7.6e-5 (2,389,913, 456 s) |

The veryheavy 2-D rows are the programs of `lpisettings('heavy')` on the same executives (storage, weight
raise 2 and eppos 1e-2 coincide with the tier); 10,514 and 10,648 rows, Gram blocks up to 1,026.

Against the stock-derived levels the executives used before (`lpisettings`): at light, P-form stability
falls from 1,000 to 280-568 variables, well-posedness from 2,325 to 625, and 2-D stability at f = 0.5 from
187,984 variables and 6.1 s to 103,744 and 2.3 s; at f = 0.9 the extra pass costs 9.5 s against 6.5 s.
Gains and norms at light, heavy and veryheavy are the programs of the same `lpisettings` levels. The stock
stripped levels certify no gain or norm in 1-D, and the stock 2-D stripped storage is 4.0x the 2-D light
one (751,849 against 187,984 variables), which is why the tiers here do not follow the stock names.

## Where each value comes from

The values were recommended by three sessions on 10/28/2026, from their banked measurements, and checked
by the runs above:
- **Posopvar structure** (positive operators, slack sizing): the weight raise by kind of LPI (1-D
  feasibility certifies at dw 0 with 2.0-3.7x fewer variables in the same time; 2-D at dw 1 in 2.5 s
  against 6.6 s at dw 2; coercive 1-D at dw 2 reaches the 1.1e-6 floor; Q forms gain nothing above 1),
  the lift offset dD 1 (the default of `get_lift_degs`), the stock storage translation, and eppos 1e-2
  for stability (relative residuals unchanged over eppos 1e-6..1, measured 09/23).
- **New Posopvar and Duality** (degree loop, dual): raises by kind (a raise is one build and one solve;
  metric LPIs start at their measured weight because an UNKNOWN pass can take 15-95 s in 2-D), the dual
  and the witness off at stripped (read-back 1-3 ms per build, witness 0.2 s per gain sweep), witness order
  24 and 240 frequencies (32 and 480 at veryheavy, inferred), slack 'like' in every tier ('stock' slack
  reproduces the stock program and its UNKNOWN on every coercive form).
- **cuADMM** (solvers): MOSEK in every tier with SeDuMi as fallback, default MOSEK tolerances (12-16
  interior-point iterations, flat in the row count), simplify off. No GPU or first-order solver in any
  tier: cuADMM was 265x slower than MOSEK on the H-infinity class and never reached 1e-8; SDPT3 returned a
  nonzero termcode on 17 of 20 library examples; SDPNAL+ returned X = 0 or an indefinite X.

The runs above changed three of the recommendations:
1. stripped and light take the stock light 1-D storage: the stock stripped storage certified no H-infinity
   or H2 bound, and no coercive stability at f = 0.9;
2. heavy takes the stock stripped 2-D storage, and veryheavy the stock heavy 2-D storage (12.6x the
   variables and 73x the time of light for 2-D stability);
3. veryheavy keeps dw_P = 2: dw_P = 3 gave the coercive gap 1.0e-6 against 8.9e-7, with 33% more variables.

## Not settings

The functions size these, and the tiers leave them alone: the support the slack is sized from ('like',
passed by the executives), the lift and weight rule and dD (`get_lift_degs`), the term set ('auto': plain
+ product at w-1 in 1-D, plain + faces at w-1 in N-D), the zero-diagonal prune, the construction route
(`poscopvar_direct`) and its cache. `margin` changes the certified statement and stays 0; `gam_fixed` is
per call; `verbose` is printing.

## Open

- eppos 0 against 1e-4 in the gain LPIs (allowed by the maintainer; not measured, under 1 min).
- SeDuMi on the container executives: every number above is MOSEK.
- Pinned-gamma bisection as a fallback for an objective LPI that ends with numerr >= 1: not implemented in
  `executives_sopvar`.
- MOSEK's thread count has changed a 2-D verdict (09/29/2026); no tier pins it.
- veryheavy 1-D on fin: the Q-form gap 5.0e-6 is larger than light's 1.6e-6 at 7.4x the variables; not
  explained.
- 2-D storages between the stock stripped and heavy 2-D levels, and (2,3)/(3,2) faces, are over the
  20-minute daytime budget.
