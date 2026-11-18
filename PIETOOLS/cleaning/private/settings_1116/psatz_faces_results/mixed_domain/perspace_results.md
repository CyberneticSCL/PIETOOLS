# Per-space slack degrees ("perspace"): results (2026-11-17)

Scratch folder `gapstudy_b\`. Shadow B: `shadow\get_lift_degs.m` (copy of `psatz_faces\shadow\get_lift_degs.m`
with options perspace, dw_space; lines marked SCRATCH B), `shadow\lpi_ineq_sop.m` (repository file plus
the tier merge of wcross, perspace, dw_space, and a size-only return under appdata GS_ESTIMATE_ONLY).
Harness copies: `gs_exec_row.m`, `gs_qdirect.m`, `gs_lb_qdirect.m`; ladder `gs_driver.m`; runner
`gs_run.sh` (slot protocol, watchdog 1,800 s). No repository file edited. Log: `logs\gs_verify.log`
(46 s wall, 05:02 – 05:02).

**Status:** run 1 of the verification failed check (2) (below, kept as the record); the coordinator chose option 1 refined (11/17 05:10): per-space numbers act in owned directions only, D and w of a non-owned direction are the registry values, and check (2) became "owned D_{k,d} <= D_d; single-space identical; slack sizes reported". Run 2 passes; the ladder runs on the fixed shadow.

## Verification run 1 (gs_verify.m, 05:02; the rule as first specified)

Targets per PIE: TT = T'T (like none, default opts); Dlike = -(A'Q+Q'A) with like 0.5(T'Q+Q'T), dw 0;
Dself = -(A'Q+Q'A), like none, dw 0. Q = lpivar_cdopvar(prog,T,{mult 1, int [1 1], out 1, in 1}).

(1) Default output identical to the psatz_faces shadow (isequal on deg and terms): **pass on all 18
(6 PIEs x 3 targets)**.

(2) perspace: w_{k,d} <= w_d and D_{k,d} <= D_d, identity on single-space PIEs. Single-space identity
passes (AGE, heat2d, all targets). **The bound w_{k,d} <= w_d fails on 3 of 12 multi-space targets.**

(4) gs_exec_row EXTRA routing: perspace, dw_space, wcross land in st.lpi.ineq and reach get_lift_degs
through the tier (the executive call printed `@@ GLD ... perspace=1`); pass.

(3) Registry and per-space numbers (rows of w_s, D_s are spaces in the order listed; vars s1, s2) and
slack sizes from lpi_ineq_sop on the program holding Q (TT sized with the tier's dw 1, wmin 2):

| PIE (spaces) | target | registry w, D | per-space w_s | per-space D_s | slack ndv rule | slack ndv perspace | (2) |
|---|---|---|---|---|---|---|---|
| S3 ({s1} {s1,s2}) | TT | [2 2], [2 2] | [1 2; 2 2] | [2 1; 2 2] | 255,780 | 255,780 | pass |
| S3 | Dlike | [1 1], [2 2] | [1 1; 1 1] | [2 1; 2 2] | 242,000 | 242,000 | pass |
| S3 | Dself | [1 1], [2 2] | [1 1; 1 1] | [2 1; 2 2] | 242,000 | 242,000 | pass |
| S17 ({s1} {s1,s2}) | TT | [1 1], [1 2] | **[1 2]**; [1 1] | [1 1; 1 2] | 14,500 | 14,500 | **fail** |
| S17 | Dlike | [1 1], [2 2] | [1 1; 1 1] | [2 1; 2 2] | 250,880 | 250,880 | pass |
| S17 | Dself | [1 1], [2 2] | [1 1; 1 1] | [2 1; 2 2] | 250,880 | 250,880 | pass |
| S18 ({s1} {s1,s2}) | TT | [1 1], [2 1] | [1 1; 1 1] | [2 1; 2 1] | 130,500 | 130,500 | pass |
| S18 | Dlike | [1 0], [2 2] | **[1 1]**; [1 0] | [2 1; 2 2] | 48,020 | **62,720** | **fail** |
| S18 | Dself | [1 1], [2 2] | [1 0; 1 1] | [2 1; 2 2] | 192,080 | 165,620 | pass |
| AGE ({s1,s2}) | TT / Dlike / Dself | [1 1],[2 2] / [0 1],[2 2] / [1 1],[2 2] | = registry | = registry | 187,920 / 46,080 / 184,320 | same | pass |
| FO2s ({s1} {s1} {s1,s2}) | TT | [1 1], [2 1] | [1 1; 1 1; 1 1] | [2 1; 1 1; 2 1] | 139,345 | 139,345 | pass |
| FO2s | Dlike | [1 0], [2 2] | **[1 1; 1 0]**; [1 0] | [2 1; 2 1; 2 2] | 62,720 | **79,380** | **fail** |
| FO2s | Dself | [1 1], [2 2] | [1 0; 1 0; 1 1] | [2 1; 2 1; 2 2] | 250,880 | 192,080 | pass |
| heat2d ({s1,s2}) | TT / Dlike / Dself | [2 2],[2 2] / [1 1],[2 2] / [1 1],[2 2] | = registry | = registry | 187,920 / 184,320 / 184,320 | same | pass |

(Sizes measured; S3 TT and the D_s column read from the log. D_{k,d} <= D_d holds everywhere.)

## Why (2) fails (diagnosis, from the numbers above and the rule)

Every violation is in a direction the space does not own: the 1-D spaces {s1}, direction s2. The rule
w = max(ceil((Dmax - D - 1)/2), ...) decreases in D. Per space, D_{k,d} = Dmin_{k,d} + dD is lower than
the registry D_d in a non-owned direction (S18 Dlike, space {s1}, s2: D 2 -> 1), while Dmax_{k,d}
there keeps the registry value (the one-sided exponent of the coupling block between {s1} and
{s1,s2} enters both spaces), so the weight rises (0 -> 1). In a non-owned direction D_{k,d} is never
used as a lift degree (mult = 0 there), so the trade the rule makes is with a lift the space does not
have. The specification's invariant w_{k,d} <= w_d therefore does not follow from the specified rule;
the implementation follows the specification (inferred from the code path; the failing rows are exactly
the non-owned directions).

Consequence measured: under perspace the decay slack of S18 (Dlike, the executive-like sizing) grows
31% (48,020 -> 62,720) and that of FO2s 27% (62,720 -> 79,380); sized from the operator itself (Dself)
perspace saves 14% (S18) and 23% (FO2s); on S3 and S17 Dlike/Dself nothing changes.

## Options offered after run 1 (the coordinator chose 1, refined; see Status)

1. Per space, use D_{k,d} only in owned directions; in a non-owned direction take the registry D_d in
   the w rule (or omit the Dmax term there).
2. Clip: w_{k,d} = min(w_{k,d}, w_d) (restores (2), but then the support-reach argument of the rule
   is not met in that direction for that space; unverified whether that loses certificates).
3. Keep the rule as specified and drop invariant (2); the ladder would then measure a variant that is
   larger than the registry rule on S18 and FO2s in the like-sized form.

## Verification run 2 (owned directions only; log gs_verify.log)

(1) default identical to the psatz_faces shadow: pass on all 18. (2) owned D_{k,d} <= D_d: pass on all 18;
single-space identity (AGE, heat2d): pass. (4) unchanged (routing pass).

Per-space numbers and slack sizes (the per-space rows of w_s, D_s now carry the registry values in the
non-owned directions):

| PIE | target | registry w, D | w_s | D_s | slack ndv rule | slack ndv perspace |
|---|---|---|---|---|---|---|
| S3 | TT | [2 2], [2 2] | [1 2; 2 2] | [2 2; 2 2] | 255,780 | 255,780 |
| S3 | Dlike / Dself | [1 1], [2 2] | = registry | = registry | 242,000 / 242,000 | 242,000 / 242,000 |
| S17 | TT | [1 1], [1 2] | = registry | = registry | 14,500 | 14,500 |
| S17 | Dlike / Dself | [1 1], [2 2] | = registry | = registry | 250,880 / 250,880 | 250,880 / 250,880 |
| S18 | TT | [1 1], [2 1] | = registry | = registry | 130,500 | 130,500 |
| S18 | Dlike | [1 0], [2 2] | = registry | = registry | 48,020 | 48,020 |
| S18 | Dself | [1 1], [2 2] | = registry | = registry | 192,080 | 192,080 |
| AGE (1 space) | TT / Dlike / Dself | [1 1],[2 2] / [0 1],[2 2] / [1 1],[2 2] | = registry | = registry | 187,920 / 46,080 / 184,320 | same |
| FO2s | TT | [1 1], [2 1] | [1 1; 1 1; 1 1] | [2 1; 1 1; 2 1] | 139,345 | 139,345 |
| FO2s | Dlike | [1 0], [2 2] | = registry | = registry | 62,720 | 62,720 |
| FO2s | Dself | [1 1], [2 2] | = registry | = registry | 250,880 | 250,880 |
| heat2d (1 space) | TT / Dlike / Dself | [2 2],[2 2] / [1 1],[2 2] / [1 1],[2 2] | = registry | = registry | 187,920 / 184,320 / 184,320 | same |

Reading. With the owned-direction rule the per-space slack equals the registry slack on every decay
target of the five PIEs: in its owned directions each space attains the registry support numbers (the
largest-degree blocks are the diagonal blocks of the space owning the direction, or a coupling block that
enters both spaces). The only per-space differences are on T'T (S3: w 1 instead of 2 in s1 for {s1};
FO2s: D 1 instead of 2 in s1 for the second {s1} space), and there the tier floor wmin = 2 (TT sizing)
makes the slack identical. Hence configuration (b) is expected to reproduce (a) on these cases
(inferred from these sizes; the ladder measures it), and the content of the ladder is the per-space
raise dw_space.

## Ladder conventions

- **direct**: `gs_qdirect` (suite) / `gs_lb_qdirect` (literature): free Q with the shape of T, degrees
  {mult 1, int [1 1], out 1, in 1}, T'Q = Q'T, T'Q >= 1e-4 T'T, -(A'Q + Q'A) >= 0, both inequalities
  sized by the reader (dw 0), light tier, MOSEK. in = optimal with eta < 1e-8; out = infeasible; U = other.
- **P**: `gs_exec_row` task stabP (PIETOOLS_stability_sop 'P'), light tier, the executive's own slack
  (dw_feas 1 in 2-D, one raise allowed); storage 'stock' or 'Rlight0'. Status as the executive returns it.
- **Configurations**: (a) the rule, terms 'auto'; (b) perspace; (c) perspace + wcross 0; (d) perspace +
  dw_space 1 on one L2 space (named by its variables); (e) dw_space 2 on the space that gained most in
  (d), or 1 on every L2 space if none gained; (f) the best of (b)-(e) with the faces at offset 0 and at
  the scalar offset 1. In the AGE and S18 runs (f) chose among (b), (d), (e) only (driver defect: (c)
  was not a candidate; fixed before S17 and FO2s). On AGE and S18 the choice is unaffected: (c) equalled
  (b) on AGE and certified nothing on S18.
- **terms (realized)**: the codes and offsets `get_lift_degs` printed (`@@ GLD`). Every row has codes
  [0 3 4 5 6] (the plain term and the four faces); "faces o0" = offsets 0, "faces o1" = offsets
  [0 1 1 1 1]. Direct: coercive / decay slack. P: the slack of the last pass. The harness name follows.
- **w (per space)**: the weight degrees `deg{k}.int` of the decay slack (rows = spaces, columns s1, s2)
  as used, after dw_space and wcross; P: of the last pass. The lift degrees D equal the rule in every row.
- **size vs (a)**: ndv relative to (a) of the same form and storage.
- Uniform-raise comparisons cite `gapstudy\GAP_RESULTS.md` (the gap study) and the T4 record; not rerun.

## AGE (item 16), b = 0.9 b* = 2.0724165; one L2 space {s1,s2} (`gs_AGE_direct`, `gs_AGE_P`)

| form | storage | config | terms (realized) | harness | dw_space | wcross | w | status | ndv | m or eta | t (s) | size vs (a) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| direct | – | a | faces o0 / o0 | auto | – | rule | [1 1] | in | 210,340 | eta 2.0e-10 | 6 | – |
| direct | – | b | faces o0 / o0 | auto | 0 | rule | [1 1] | in | 210,340 | eta 2.0e-10 | 5 | 0 % |
| direct | – | c | faces o0 / o0 | auto | 0 | 0 | [1 1] | in | 210,340 | eta 2.0e-10 | 5 | 0 % |
| direct | – | d {s1,s2} | faces o0 / o1 | auto | 1 | rule | [2 2] | in | 567,460 | eta 1.8e-10 | 19 | +170 % |
| direct | – | e all +1 | faces o0 / o1 | auto | 1 | rule | [2 2] | in | 567,460 | eta 1.8e-10 | 20 | +170 % |
| direct | – | f (b) | faces o0 / o0 | faces0 | 0 | rule | [1 1] | in | 210,340 | eta 2.0e-10 | 5 | 0 % |
| direct | – | f (b) | faces o1 / o1 | facesS | 0 | rule | [1 1] | in | 56,548 | eta 2.4e-10 | 2 | −73 % |
| P | stock | a | faces o1 | auto | – | rule | [2 2] | infeasible (1 raise) | 187,984 | m 2,968 | 13 | – |
| P | stock | b | faces o1 | auto | 0 | rule | [2 2] | infeasible | 187,984 | m 2,968 | 16 | 0 % |
| P | stock | c | faces o1 | auto | 0 | 0 | [2 2] | infeasible | 187,984 | m 2,968 | 16 | 0 % |
| P | stock | d {s1,s2} | faces o1 | auto | 1 | rule | [3 3] | infeasible | 751,744 | m 5,104 | 83 | +300 % |
| P | stock | e all +1 | faces o1 | auto | 1 | rule | [3 3] | infeasible | 751,744 | m 5,104 | 74 | +300 % |
| P | stock | f (b) | faces o0 | faces0 | 0 | rule | [2 2] | infeasible | 524,944 | m 3,892 | 37 | +179 % |
| P | stock | f (b) | faces o1 | facesS | 0 | rule | [2 2] | infeasible | 187,984 | m 2,968 | 14 | 0 % |
| P | Rlight0 | a | faces o1 | auto | – | rule | [2 2] | infeasible (1 raise) | 188,545 | m 2,983 | 17 | – |
| P | Rlight0 | b | faces o1 | auto | 0 | rule | [2 2] | infeasible | 188,545 | m 2,983 | 17 | 0 % |
| P | Rlight0 | c | faces o1 | auto | 0 | 0 | [2 2] | infeasible | 188,545 | m 2,983 | 16 | 0 % |
| P | Rlight0 | d {s1,s2} | faces o1 | auto | 1 | rule | [3 3] | infeasible | 752,305 | m 5,119 | 77 | +299 % |
| P | Rlight0 | e all +1 | faces o1 | auto | 1 | rule | [3 3] | infeasible | 752,305 | m 5,119 | 77 | +299 % |
| P | Rlight0 | f (b) | faces o0 | faces0 | 0 | rule | [2 2] | infeasible | 525,505 | m 3,907 | 37 | +179 % |
| P | Rlight0 | f (b) | faces o1 | facesS | 0 | rule | [2 2] | infeasible | 188,545 | m 2,983 | 15 | 0 % |

AGE has one space, so the per-space rule is the registry rule (verified) and dw_space is the uniform
raise. The rows reproduce the LITBUILD and gap-study sizes (direct 210,340 and 56,548; P 187,984; the
gap study's P at dw_feas 3, 751,744, equals (d) here). **No configuration changes a verdict.** The
direct form is certified at the baseline, and the smallest certificate is the scalar-offset faces at
56,548 (−73 %), the gap-study result. The 'P' gap is not a slack gap: no slack raise certifies
(751,744, infeasible), whereas the gap study closes it with the storage (Rheavy_u, dw_feas 0, 104,080,
certified).

## S18 (suite, library PDE 43 with a delay state), tau = 0.9 tau* = 0.9312534; spaces {s1}, {s1,s2} (`gs_S18_direct`, `gs_S18_P`)

| form | storage | config | terms (realized) | harness | dw_space | wcross | w ({s1}; {s1,s2}) | status | ndv | m or eta | t (s) | size vs (a) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| direct | – | a | faces o0 / o0 | auto | – | rule | [1 1; 1 1] | U | 227,510 | eta 5.2e-10 | 13 | – |
| direct | – | b | faces o0 / o0 | auto | [0 0] | rule | [1 1; 1 1] | U | 227,510 | eta 5.2e-10 | 10 | 0 % |
| direct | – | c | faces o0 / o0 | auto | [0 0] | 0 | [1 0; 1 1] | U | 201,050 | eta 6.3e-10 | 10 | −12 % |
| direct | – | d {s1} | faces o0 / o0 | auto | [1 0] | rule | [2 2; 1 1] | **in** | 325,275 | eta 8.4e-11 | 16 | +43 % |
| direct | – | d {s1,s2} | faces o0 / o0 | auto | [0 1] | rule | [1 1; 2 2] | **in** | 1,084,250 | eta 1.3e-10 | 88 | +377 % |
| direct | – | e {s1} +2 | faces o0 / o0 | auto | [2 0] | rule | [3 3; 1 1] | in | 495,830 | eta 5.9e-10 | 20 | +118 % |
| direct | – | f (d {s1}) | faces o0 / o0 | faces0 | [1 0] | rule | [2 2; 1 1] | in | 325,275 | eta 8.4e-11 | 11 | +43 % |
| direct | – | f (d {s1}) | faces o1 / o1 | facesS | [1 0] | rule | [2 2; 1 1] | **in** | **93,991** | eta 1.8e-10 | 5 | **−59 %** |
| P | stock | a | faces o1 | auto | – | rule | [3 2; 3 2] | infeasible (1 raise) | 508,228 | m 4,757 | 45 | – |
| P | stock | b | faces o1 | auto | [0 0] | rule | [3 2; 3 2] | infeasible | 508,228 | m 4,757 | 45 | 0 % |
| P | stock | c | faces o1 | auto | [0 0] | 0 | [3 0; 3 2] | infeasible | 427,156 | m 4,555 | 55 | −16 % |
| P | stock | d {s1} | faces o1 | auto | [1 0] | rule | [4 3; 3 2] | infeasible | 636,676 | m 4,978 | 62 | +25 % |
| P | stock | d {s1,s2} | faces o1 | auto | [0 1] | rule | [3 2; 4 3] | infeasible | 1,503,556 | m 7,390 | 401 | +196 % |
| P | stock | e all +1 | faces o1 | auto | [1 1] | rule | [4 3; 4 3] | infeasible | 1,721,860 | m 7,635 | 407 | +239 % |
| P | stock | f (b) | faces o0 | faces0 | [0 0] | rule | [3 2; 3 2] | infeasible | 1,270,276 | m 6,028 | 114 | +150 % |
| P | stock | f (b) | faces o1 | facesS | [0 0] | rule | [3 2; 3 2] | infeasible | 508,228 | m 4,757 | 27 | 0 % |
| P | Rlight0 | a | faces o1 | auto | – | rule | [3 2; 3 2] | infeasible (1 raise) | 510,057 | m 4,841 | 44 | – |
| P | Rlight0 | b | faces o1 | auto | [0 0] | rule | [3 2; 3 2] | infeasible | 510,057 | m 4,841 | 44 | 0 % |
| P | Rlight0 | c | faces o1 | auto | [0 0] | 0 | [3 0; 3 2] | infeasible | 428,985 | m 4,639 | 36 | −16 % |
| P | Rlight0 | d {s1} | faces o1 | auto | [1 0] | rule | [4 3; 3 2] | infeasible | 638,505 | m 5,062 | 58 | +25 % |
| P | Rlight0 | d {s1,s2} | faces o1 | auto | [0 1] | rule | [3 2; 4 3] | infeasible | 1,505,385 | m 7,474 | 285 | +195 % |
| P | Rlight0 | e, f | – | – | – | – | – | not run: watchdog kill at 1,800 s during (e) | – | – | – | – |

The per-space rule changes no S18 verdict and no size: (b) = (a) in both forms, because the decay
slack has the registry degrees in the owned directions (verification run 2). **A per-space raise closes
the direct-form gap.** Raising the weight of the 1-D space {s1} alone by 1 certifies at 325,275
variables with the faces at offset 0 and at **93,991** with the faces at the scalar offset 1
(eta 1.8e-10). Raising the 2-D space alone also certifies, at 1,084,250 (3.3x the {s1} raise). The
uniform raise dw 1 at the same Q degrees certifies at 665,815 ('auto') and 440,023 (facesS) (gap study,
T4); the gap study's smallest dw-1 certificate is 336,623 (faces at per-direction offset 2, wcross 0).
The {s1} raise with facesS is 4.7x smaller than uniform facesS and 3.6x smaller than that best dw-1
certificate. It is also below the gap study's dw-0 closure by the Q-degree raise (Q2, facesS, wcross 0:
215,188). The raise of {s1} also raised its non-owned direction s2 (w [2 2] for {s1}); a raise of the
owned direction alone (with wcross 0) was not run. In 'P' no slack configuration certifies, up to
1,721,860 variables (stock). As in the gap study, the 'P' gap is the storage (Rlight2 closes it at
953,001); wcross 0 saves 16 % without changing the verdict.

## S17 (suite, transport edge state driving a 2-D diffusion), r2 = 0.9 pi^2; spaces {s1}, {s1,s2} (`gs_S17_direct`, `gs_S17_P`)

'P' run with the stock storage only (MATLAB budget; on AGE and S18, Rlight0 returned the verdicts of
stock at +0.3 % – +0.4 % of the size).

| form | storage | config | terms (realized) | harness | dw_space | wcross | w ({s1}; {s1,s2}) | status | ndv | m or eta | t (s) | size vs (a) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| direct | – | a | faces o0 / o0 | auto | – | rule | [1 1; 1 1] | U (numerr 2) | 435,350 | eta 8.6e-9 | 21 | – |
| direct | – | b | faces o0 / o0 | auto | [0 0] | rule | [1 1; 1 1] | U (numerr 2) | 435,350 | eta 8.6e-9 | 18 | 0 % |
| direct | – | c | faces o0 / o0 | auto | [0 0] | 0 | [1 0; 1 1] | U (numerr 2) | 382,650 | eta 1.1e-9 | 17 | −12 % |
| direct | – | d {s1} | faces o0 / o0 | auto | [1 0] | rule | [2 2; 1 1] | U (numerr 2) | 581,975 | eta 5.5e-9 | 27 | +34 % |
| direct | – | d {s1,s2} | faces o0 / o0 | auto | [0 1] | rule | [1 1; 2 2] | U (numerr 2) | 1,907,975 | eta 7.2e-9 | 222 | +338 % |
| direct | – | e all +1 | faces o1 / o1 | auto | [1 1] | rule | [2 2; 2 2] | U (numerr 2) | 788,950 | eta 2.0e-9 | 56 | +81 % |
| direct | – | f (c) | faces o0 / o0 | faces0 | [0 0] | 0 | [1 0; 1 1] | U (numerr 2) | 382,650 | eta 1.1e-9 | 15 | −12 % |
| direct | – | f (c) | faces o1 / o1 | facesS | [0 0] | 0 | [1 0; 1 1] | out | 98,410 | eta 8e-3 | 5 | −77 % |
| P | stock | a | faces o1 | auto | – | rule | [3 2; 3 2] | unknown (1 raise) | 663,748 | m 5,070 | 144 | – |
| P | stock | b | faces o1 | auto | [0 0] | rule | [3 2; 3 2] | unknown | 663,748 | m 5,070 | 126 | 0 % |
| P | stock | c | faces o1 | auto | [0 0] | 0 | [3 0; 3 2] | unknown | 570,580 | m 4,866 | 74 | −14 % |
| P | stock | d {s1} | faces o1 | auto | [1 0] | rule | [4 3; 3 2] | unknown | 809,476 | m 5,292 | 165 | +22 % |
| P | stock | d {s1,s2} | faces o1 | auto | [0 1] | rule | [3 2; 4 3] | unknown | 1,998,340 | m 7,902 | 730 | +201 % |
| P | stock | e, f | – | – | [1 1] | – | [4 3; 4 3] | not finished: watchdog kill at 1,800 s during (e) | – | – | – | – |

The per-space rule changes no S17 verdict and no size ((b) = (a)). **No per-space raise closes the
direct gap at Q length 2.** Raising {s1} or {s1,s2} alone, or both (e), returns UNKNOWN with MOSEK
numerical problems (numerr 2) up to 1,907,975 variables. The uniform raise behaves the same way in the
gap study: at Q length 2, dw 1 and dw 2 return UNKNOWN for facesS and faces0, and the only Q-length-2
certificate is the faces at per-direction offset 2 at dw 1 (636,630). (e) here, dw_space [1 1], is the
uniform raise dw 1 with 'auto' resolving to the scalar-offset faces; its size, 788,950, equals the T4
facesS dw-1 row, also UNKNOWN. The gap study closes S17 by raising the degrees of Q to 2 (facesS,
wcross 0, 267,968), a knob outside this ladder. 'P' stays UNKNOWN under every slack configuration run,
consistent with the record that S17 is uncertifiable in the P form.

## FO2s (item 9, S18 form of the delayed damped wave), tau = 0.95 tau* = 2.737331; spaces {s1} (z), {s1} (v), {s1,s2} (y) (`gs_FO2s_direct`, `gs_FO2s_P`)

The two 1-D spaces are listed in the order of the PIE (the second is v, of order 0; inferred from the
order of the pde_var declarations, not read out). The point 0.95 tau* lies above the LITBUILD direct
reach (faces0 bisection, 0.91602 tau*).

| form | storage | config | terms (realized) | harness | dw_space | wcross | w ({s1}; {s1}; {s1,s2}) | status | ndv | m or eta | t (s) | size vs (a) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| direct | – | a | faces o0 / o0 | auto | – | rule | [1 1; 1 1; 1 1] | U (numerr 2, feasratio −1) | 299,120 | eta 3.3e-2 | 22 | – |
| direct | – | b | faces o0 / o0 | auto | [0 0 0] | rule | [1 1; 1 1; 1 1] | U (numerr 2) | 299,120 | eta 3.3e-2 | 18 | 0 % |
| direct | – | c | faces o0 / o0 | auto | [0 0 0] | 0 | [1 0; 1 0; 1 1] | U (numerr 2) | 240,320 | eta 0.14 | 16 | −20 % |
| direct | – | d {s1} (1st) | faces o0 / o0 | auto | [1 0 0] | rule | [2 2; 1 1; 1 1] | U (numerr 2) | 410,045 | eta 2.1e-2 | 24 | +37 % |
| direct | – | d {s1} (2nd) | faces o0 / o0 | auto | [0 1 0] | rule | [1 1; 2 2; 1 1] | U (numerr 2) | 415,005 | eta 7.2e-2 | 28 | +39 % |
| direct | – | d {s1,s2} | faces o0 / o0 | auto | [0 0 1] | rule | [1 1; 1 1; 2 2] | U (numerr 2) | 1,234,820 | eta 0.19 | 134 | +313 % |
| direct | – | e all +1 | faces o0 / o1 | auto | [1 1 1] | rule | [2 2; 2 2; 2 2] | U (numerr 2) | 887,120 | eta 9.2e-4 | 82 | +197 % |
| direct | – | f (c) | faces o0 / o0 | faces0 | [0 0 0] | 0 | [1 0; 1 0; 1 1] | U (numerr 2) | 240,320 | eta 0.14 | 16 | −20 % |
| direct | – | f (c) | faces o1 / o1 | facesS | [0 0 0] | 0 | [1 0; 1 0; 1 1] | out | 70,388 | eta 1.0 | 8 | −76 % |
| P | stock | a | faces o1 | auto | – | rule | [3 2; 3 2; 3 2] | infeasible (1 raise) | 691,888 | m 5,736 | 105 | – |
| P | stock | b | faces o1 | auto | [0 0 0] | rule | [3 2; 3 2; 3 2] | infeasible | 691,888 | m 5,736 | 102 | 0 % |
| P | stock | c | faces o1 | auto | [0 0 0] | 0 | [3 0; 3 0; 3 2] | infeasible | 494,756 | m 5,330 | 71 | −28 % |
| P | stock | d {s1} (1st) | faces o1 | auto | [1 0 0] | rule | [4 3; 3 2; 3 2] | infeasible | 840,496 | m 5,969 | 140 | +21 % |
| P | stock | d {s1} (2nd) | faces o1 | auto | [0 1 0] | rule | [3 2; 4 3; 3 2] | infeasible | 866,720 | m 5,972 | 103 | +25 % |
| P | stock | d {s1,s2} | faces o1 | auto | [0 0 1] | rule | [3 2; 3 2; 4 3] | infeasible | 1,808,176 | m 8,579 | 578 | +161 % |
| P | stock | e, f | – | – | [1 1 1] | – | [4 3; 4 3; 4 3] | not finished: watchdog kill at 1,800 s during (e) | – | – | – | – |

'P' run with the stock storage only (MATLAB budget). The per-space rule changes no FO2s verdict and no
size ((b) = (a)). In the direct form every configuration returns UNKNOWN (MOSEK numerr 2, feasratio −1,
eta 9.2e-4 – 0.19) or, with the scalar-offset faces, infeasible. **No per-space raise closes the gap at
0.95 tau*.** The uniform raise ((e), dw_space [1 1 1], 887,120) is also UNKNOWN. This point lies above
the measured direct reach (0.91602 tau*), so the negative result is consistent with the earlier bisection;
whether a per-space raise lifts that reach was not measured (the gap list fixes the point at 0.95).
In the gap study the uniform raise dw 1 at Q length 2 is also UNKNOWN (facesS wcross 0, 433,725; faces2
wcross 0, 394,961), and FO2s at 0.95 closes only by raising Q to degrees 2 (facesS, wcross 0, 248,497). In 'P' every configuration is infeasible up to 1,808,176 variables; wcross
0 saves 28 %.

## D2 (item 8, 2-D heat with a 3-D delay state), tau = 0.5 tau* = 0.517363; spaces {s1,s2}, {s1,s2,s3} (`gs_D2_direct`, size only)

Direct form; Q and the symmetry equality: 1,500 variables. The two slacks counted by the estimate of
`lpi_ineq_sop` (scratch return under appdata GS_ESTIMATE_ONLY: nothing declared, never solved), dw 0.
Realized terms: 'auto' resolves to faces at offset 0 here (its sizes equal the faces0 rows).

| config | dw_space | wcross | terms (realized) | harness | coercive slack | decay slack | total ndv | vs (a) | vs cap 3,000,000 |
|---|---|---|---|---|---|---|---|---|---|
| a | – | rule | faces o0 | auto / faces0 | 7,112,448 | 40,320,000 | 47,433,948 | – | above |
| a | – | rule | faces o1 | facesS | 1,397,088 | 6,300,000 | 7,698,588 | −84 % | above |
| b | [0 0] | rule | faces o0 | auto / faces0 | 7,112,448 | 40,320,000 | 47,433,948 | 0 % | above |
| b | [0 0] | rule | faces o1 | facesS | 1,397,088 | 6,300,000 | 7,698,588 | −84 % | above |
| c | [0 0] | 0 | faces o0 | auto / faces0 | 7,112,448 | 34,126,848 | 41,240,796 | −13 % | above |
| c | [0 0] | 0 | faces o1 | facesS | 1,397,088 | 5,415,264 | 6,813,852 | −86 % | above |
| d {s1,s2} | [1 0] | rule | faces o0 | auto / faces0 | 16,003,008 | 76,785,408 | 92,789,916 | +96 % | above |
| d {s1,s2} | [1 0] | rule | faces o1 | facesS | 3,063,744 | 13,396,320 | 16,461,564 | −65 % | above |
| d {s1,s2,s3} | [0 1] | rule | faces o0 | auto / faces0 | 113,799,168 | 361,671,408 | 475,472,076 | +902 % | above |
| d {s1,s2,s3} | [0 1] | rule | faces o1 | facesS | 21,117,024 | 77,227,920 | 98,346,444 | +107 % | above |

The rule row reproduces the LITBUILD size (47,433,948). The per-space rule changes nothing ((b) = (a)).
The smallest program of the ladder, scalar-offset faces with wcross 0, has 6,813,852 variables, 2.3x the
cap; no configuration fits, so none was solved, and (e), (f) were not formed.

## Summary

| case, point | form | baseline (a) | per-space rule (b) | smallest certificate in the ladder | uniform raise (gap study / T4 / (e)) |
|---|---|---|---|---|---|
| AGE 0.9 b* | direct | in, 210,340 | = (a) | in, 56,548 (faces o1, dw_space 0) | one space: dw_space = uniform; 567,460 in |
| AGE 0.9 b* | P (stock, Rlight0) | infeasible, 187,984 / 188,545 | = (a) | none (to 751,744) | infeasible 751,744 (gap study); closed by storage Rheavy_u, 104,080 |
| S18 0.9 tau* | direct | U, 227,510 | = (a) | **in, 93,991** (dw_space [1 0] on {s1}, faces o1) | dw 1: 'auto' 665,815, facesS 440,023, best 336,623 (faces2, wcross 0) |
| S18 0.9 tau* | P (stock, Rlight0) | infeasible, 508,228 / 510,057 | = (a) | none (to 1,721,860) | closed by storage Rlight2, 953,001 (gap study) |
| S17 0.9 r2* | direct | U, 435,350 | = (a) | none (U to 1,907,975) | dw 1 U (788,950); closed by Q degree 2, 267,968 (gap study) |
| S17 0.9 r2* | P (stock) | unknown, 663,748 | = (a) | none (to 1,998,340; (e), (f) not finished) | proven uncertifiable in P |
| FO2s 0.95 tau* | direct | U, 299,120 | = (a) | none (U to 1,234,820; facesS out) | dw 1 U (433,725; (e) 887,120); closed by Q degree 2, 248,497 (gap study) |
| FO2s 0.95 tau* | P (stock) | infeasible, 691,888 | = (a) | none (to 1,808,176; (e), (f) not finished) | not closed in the gap study (Rlight2, Rheavy_s unknown) |
| D2 0.5 tau* | direct (size only) | 47,433,948 | = (a) | smallest program 6,813,852 (faces o1, wcross 0): above the cap | – |

MATLAB wall: verification 2 x 1 min; AGE 10 min; S18 33 min; S17 37 min; FO2s 37 min; D2 1 min
(about 2 h in all). Three P runs were stopped by the 1,800 s watchdog during (e): S18 with Rlight0, S17,
and FO2s.

**What the per-space rule is worth.** With the support numbers read per space in the owned directions
only, the rule returned the registry degrees on every decay slack of the five PIEs. It changed no
verdict and no size on any case (verification run 2 and configuration (b) everywhere), so as a sizing
rule it is inert on this suite. The per-space raise dw_space is the useful knob: on S18 in the direct
form, a weight raise of the 1-D space alone certifies at 93,991 variables with the scalar-offset faces,
4.7x below the uniform raise with the same terms (440,023) and below every closure in the gap study.
It closes nothing on S17, FO2s or AGE-'P', whose gaps lie in the Q degree, beyond the reach, or in the
storage. wcross 0 cut 12 – 28 % in every multi-space row and changed no verdict.
