# Results of the autonomous psatz-faces run (2026-11-15/16)

Plan: `sopvar/psatz_faces_test_plan.md`; decisions: `DECISIONS.md` (this folder); evidence note: `sopvar/psatz_terms_evidence.md`. Instrument: scratch shadows of `poscopvar_direct` (option `psatz_rows_perdir`), `lpi_ineq_sop` (tier knob rows) and `heatNd_lpi` (term case 'rows'); no repo file changed. Tree: ndopvar at faca4e51 with the uncommitted files of `git status`. Workstation, MATLAB R2025b, MOSEK 11 (TBB cap: 7 workers; MOSEK reports 24 threads). Every number below is read from the run logs (`c01_summary.txt` beside this file; full logs in the scratchpad `psatz_faces/logs/`).

Term sets (PF_TERMS): `none` plain term only; `box` the Putinar set {1, g_1, g_2}; `product` the 2^N set {1, g_1, g_2, g_1 g_2}; `faces` the plain term and the 2N faces; `facesS` the faces under the scalar offset rule (today's 'auto'); `corners` the faces and the four two-face products; a trailing 0 means offset 0 (the same Gram for every term, the structure of E1). Offsets o = 1 lower the weight by one in the acting directions only (per-direction rule). The lift is the reader's (D = Dmin + 1); dw raises the weight w of every term.

## Phase 0: instrument checks (passed)

- 1-D: the row `box` at o = 1 gives the program of the repo's 'auto' (Markov-Lukacs pair) with the same Gram sizes [12 6], m 27, nx 180 and the same At entries.
- 2-D Gram sizes at D = [2 2], w = [1 1]: plain 144; box factor 72 (lowered in its direction only); product g_1 g_2 36; face 72 at o = 1 per direction, 144 at o = 0, 36 under the scalar rule.
- One pinned-P 2-D solve at dw 0: m 1405, nx 41,472, 1 s.

## Phase 1: pinned P (C0, C1), rung sweep 0.10, 0.50, 0.90, 0.99, 0.999, 0.9999 of lam*

Verdicts: + certified (MOSEK optimal, rel_b <= 1e-6, blocks PSD, tight tolerances), - Farkas ray verified, 0 uncertain. ndec = decision variables of the certificate (sum N^2 over the Gram blocks), m = rows.

### C0, one variable, DD, lam* = pi^2 (control for T1)

| set | dw 0 (w 0) | dw 1 (w 1) | dw 2 (w 2) | ndec at dw 1 / dw 2 |
|---|---|---|---|---|
| none | ------ | ------ | ------ | 144 / 324 |
| box o1 (Markov-Lukacs pair) | ------ | ++---- | ++++++ | 180 / 468 |
| box0 | ------ | ++---- | ++++++ | 288 / 648 |
| faces o1 | ------ | ++---- | ++++++ | 216 / 612 |
| faces0 | ------ | ++---- | ++++++ | 432 / 972 |

The faces give the verdicts of the pair at every w and never fewer failures, at 1.2 to 1.5x the variables; the control agrees with T1 (evidence note) and with E11.

### C1, two variables, DD x DD, lam* = 2 pi^2

| set | dw 0 (w [1 1]) | dw 1 (w [2 2]) | ndec dw 0 / dw 1 | m dw 0 / dw 1 | max solve s |
|---|---|---|---|---|---|
| none | ------ | 0----- | 20,736 / 104,976 | 1405 / 2965 | 9 |
| box o1 | ++---- | ++++++ | 31,104 / 198,288 | 1405 / 2965 | 11 |
| box0 | ++---- | ++++++ | 62,208 / 314,928 | 2677 / 4813 | 26 |
| product o1 | ++---- | ++++++ | 32,400 / 219,024 | 1405 / 2965 | 14 |
| product0 | ++---- | ++++++ | 82,944 / 419,904 | 2965 / 5101 | 35 |
| faces o1 | ++---- | ++++++ | 41,472 / 291,600 | 1405 / 2965 | 21 |
| faces0 | ++---- | ++++++ | 103,680 / 524,880 | 2041 / 3889 | 37 |
| facesS (scalar rule) | 0----- | ++++++ | 25,920 / 187,920 | 1405 / 2965 | 11 |
| corners o1 | ++---- | ++++++ | 46,656 / 374,544 | 1405 / 2965 | 17 |
| corners0 | ++---- | ++++++ | 186,624 / 944,784 | 2113 / 3961 | 66 |

Reading. (1) At w = [2 2] every set with a box or face term certifies the top rung 0.9999; at w = [1 1] every such set certifies 0.50 and is certified infeasible at 0.90 by a Farkas ray, the faces included. The degree decides, not the term set. (2) The Putinar set {1, g_1, g_2} at per-direction offsets certifies everything the faces certify, at 68% of the faces' variables (198,288 against 291,600) and the same rows. The one separation of E1 (box 0.99 failing at 0.999 under MOSEK, faces 0.9999) does not reproduce on the container path with the lift sized from the target and per-direction offsets: H0 holds on this plant. (3) The scalar offset rule (today's 'auto') is the smallest certifying configuration at w = [2 2] (187,920) but certifies nothing at w = [1 1] (uncertain at 0.10, Farkas from 0.50), where the per-direction faces and the box set certify 0.50. (4) The corners and the offset-0 replicas add variables and no verdict. (5) No uncertain verdict anywhere except `none` and `facesS` at the lowest rung of w = 1.

## Interim, Phase 2 (C2, heat DD x DN, d = 0, R at (1,1), Q at (2, wQ)) and Phase 3 (C4, S3 at 0.9 r*), read while the jobs run

C2 so far (kappa* = 12.337006; lo = certified feasible kappa, hi = certified infeasible; per-direction offsets):

| cfg | wQ = 1: lo / hi, nx, s | wQ = 2: lo / hi, nx, s |
|---|---|---|
| none | no certified kappa (hi 1.54) | no certified kappa (hi 3.08) |
| box o1 | 12.2643306 / 12.2764630 (resolved), 61,449, 32 | 12.3250171 / 13.57 (uncertain points fill the bracket), 358,665, 379 |
| product o1 | 12.2643306 / 12.2764630, 64,009, 38 | 12.3250171 / 13.57, 395,785, 461 |
| faces o1 | 12.0886403 / 12.0946945 (resolved), 81,929, 52 | running |

At wQ = 1 the Putinar set certifies 0.9941 kappa* and the faces at the per-direction offset 0.9799: the faces at o = 1 are one degree below the box factor in their direction (2w - 1 against 2w), so this row is the faces at lower degree, not the fair match; faces0, facesS and corners follow. The box set and the 2^N set give the same lo at both wQ (the g_1 g_2 term adds nothing here), the box set at 7% fewer variables.

C4 so far (S3 edge cascade, r = 0.9 pi^2, light tier, dw pinned; verdict = the executive's):

| storage | dw | none | box | product | faces | faces0 | facesS | corners |
|---|---|---|---|---|---|---|---|---|
| stock settings_2d | 0 | infeasible | infeasible | infeasible | infeasible | infeasible | infeasible | infeasible |
| stock settings_2d | 1 | infeasible | unknown 107,800 | unknown 114,856 | unknown 151,900 | unknown 317,716 | unknown 91,924 | unknown 180,124 |
| N-D default (route) | 0 | infeasible | infeasible | infeasible | infeasible | running | | |
| N-D default (route) | 1 | infeasible | certified 108,088, 4 s | certified 115,144, 4 s | certified 152,188, 6 s | running | | |

With the N-D default storage the Putinar set, the 2^N set and the faces all certify S3 at 0.9 r* at dw 1; with the stock 2-D storage none of them does (H3). The Putinar set is the smallest certifying configuration (71% of the faces' variables).

## Phase 3 final (C4): S3 and S8 at 0.9 of the margin, light tier, dw pinned (supersedes the interim table)

Verdict = the executive's certified flag; ndv = decision variables of the whole program (storage and slack). Every dw 0 run is 'infeasible' (MOSEK) for every set and both storages; at dw 1:

| storage | none | box o1 | product o1 | faces o1 | faces0 | facesS | corners o1 |
|---|---|---|---|---|---|---|---|
| S3, stock settings_2d | infeasible 63,700 | unknown 107,800 | unknown 114,856 | unknown 151,900 | unknown 317,716 | unknown 91,924 | unknown 180,124 |
| S3, N-D default | infeasible 63,988 | certified 108,088 (4 s) | certified 115,144 (4 s) | certified 152,188 (6 s) | certified 318,004 (12 s) | certified 92,212 (5 s) | certified 180,412 (6 s) |
| S8, stock settings_2d | infeasible 187,024 | unknown 352,912 | unknown 389,776 | unknown 518,800 | unknown 933,520 | unknown 334,480 | unknown 666,256 |
| S8, N-D default | infeasible 187,648 | certified 353,536 (22 s) | certified 390,400 (23 s) | certified 519,424 (23 s) | certified 934,144 (53 s) | certified 335,104 (19 s) | certified 666,880 (31 s) |

Reading. On the executive path the term set never changes the verdict: with the N-D default storage every set with a box or face term certifies both cascades at dw 1, with the stock 2-D storage none does (H3 holds; E8's storage confound confirmed on both systems). The smallest certifying configurations are the scalar-rule faces (today's 'auto') and the Putinar box set, within 17% of each other; the per-direction faces cost 1.4x the box set, the corners 1.7x, the faces at offset 0 2.9x. The record's 318,004 (S3) and 335,104 (S8) are reproduced exactly by faces0 and facesS respectively, which is what the executive's 'auto' resolves to on each system.

## Phase 2 final at rtol 1e-3 (C2): 2-D heat DD x DN, kappa* = 12.337006, d = 0, ep = 0.1, R at (D,w) = (1,1), Q at (2,wQ)

lo = largest certified feasible kappa, hi = smallest Farkas-certified infeasible kappa; nx = PSD variables; Ks = Gram blocks (R's five or so first, then Q's). Per-direction offsets o = 1 lower the weight in the acting direction only.

| cfg | wQ = 1: lo (reach), hi, nx, Ks(Q) | wQ = 2: lo (reach), nx, Ks(Q) |
|---|---|---|
| none | no certified kappa, hi 1.542; 40,969 | no certified kappa, hi 3.084; 190,729 |
| box o1 | 12.2643306 (0.99411), hi 12.2764630, 61,449, [192 96 96] | 12.3250171 (0.99903), 358,665, [432 288 288] |
| product o1 | 12.2643306, hi 12.2764630, 64,009, [192 96 96 48] | 12.3250171, 395,785, [432 288 288 192] |
| faces o1 | 12.0886403 (0.97987), hi 12.0946945, 81,929, [192 96 x4] | 12.3250171, 526,601, [432 288 x4] |
| faces0 | 12.3250171 (0.99903), hi 12.3371738, 204,809, [192 x5] | 12.3250171, 953,609 (budget), [432 x5] |
| facesS | 12.0886403, hi 12.0946945, 51,209, [192 48 x4] | 12.3250171, 339,209, [432 192 x4] |
| corners o1 | 12.0886403, hi 12.0946945, 92,169, [192 96 x4 48 x4] | 12.3250171, 675,081, [432 288 x4 192 x4] |

Reading. (1) At wQ = 2 every set certifies lo = 12.3250171 and nothing above it within the budget: the bisection stops with uncertain points filling [lo, kappa* + tol], and tol = rtol kappa* = 0.0123 equals kappa* - lo, so the sets are tied at the resolution of this run, not shown equal (the 10/08 record separated faces at w-1, 12.336731, from the tensor set, 12.326178, at rtol 1e-4; rerun at rtol 1e-4 in progress). The sizes at wQ = 2 order facesS 339,209 (the record's 339,209 reproduced) < box 358,665 < product 395,785 < faces 526,601 < corners 675,081 < faces0 953,609. (2) At wQ = 1 the brackets are resolved and the reach orders by the total degree of the weighted terms in the acting direction: faces at o = 1 and the scalar rule and the corners (degree 2w - 1) 0.9799 < box and product (2w) 0.9941 < faces at o = 0 (2w + 1) 0.9990, which equals the wQ = 2 value of every set. This is the interleaving of T6 seen directly: the Putinar box set at per-direction offsets beats the faces at the same offset and loses to the faces one offset up. (3) For a certificate of reach 0.999 the smallest configuration is the faces at offset 0 with wQ = 1 (nx 204,809, m 3117), 43% below the box set at wQ = 2 (358,665) and 40% below the scalar-rule faces at wQ = 2 (339,209): the odd-degree face terms reach total degree 2w + 1 while the plain block stays at w, which the box set can only do by raising w. (4) The 2^N product term g1 g2 adds nothing to the box set at either wQ (same lo and hi, 4 to 10% more variables).

## Phase 1b (C1b): the stock product alone, and the converter-built pinned-P plants of E2

Stock product alone (plain + g1 g2, 'prodonly'; the 2-D `psatz = 1` of the record) on the DD x DD plant, rungs 0.10 to 0.9999 of lam*:

| cfg | dw 0 (w 1) | dw 1 (w 2) | dw 2 (w 3) | ndec dw 0 / 1 / 2 |
|---|---|---|---|---|
| prodonly o1 | ------ | 00---0 | +00000 | 22,032 / 125,712 / 436,752 |
| prodonly0 | 0----- | +000-0 | +00000 | 41,472 / 209,952 / 663,552 |

The product alone is Farkas-certified infeasible from 0.90 at w <= 2 and certifies only 0.10 at w = 3 (the rest uncertain): the record's "the product certifies nothing" rows are a property of that cone, which every other set of this run leaves behind at w = 2.

Converter-built plants, Dirichlet on four edges, absolute lam rungs 1, 2, 5, 10, 15, 19 (the grid of the 09/24 record), P = I:

| plant | set | dw 0 (w [1 1]) | dw 1 (w [2 2]) | ndec dw 0 / dw 1 |
|---|---|---|---|---|
| L1, reaction lam (1 - 2 s1) x | none | ------ | 00---- | 20,736 / 104,976 |
| L1 | box, product, faces, corners (all alike) | ++++++ | ++++++ | box 31,104 / 198,288; faces 41,472 / 291,600 |
| L1 | facesS | 00---- | ++++++ | 25,920 / 187,920 |
| L2, diffusion 0.05 + s1 (1 - s1) | none | ------ | ------ | |
| L2 | box, product, faces, corners (all alike) | +++--- | ++++-- | as L1 |
| L2 | facesS | ------ | ++++-- | |
| ODE-coupled (X_t = -2 X + int int x) | none | 0----- | +00--- | 21,025 / 105,625 |
| ODE | box, product, faces, corners (all alike) | ++++-- | ++++++ | box 31,683 / 199,803; faces 42,341 / 293,981 |
| ODE | facesS | 000--- | ++++++ | 26,501 / 189,725 |

Reading. On every plant the Putinar box set, the 2^N set, the per-direction faces and the corners give identical verdicts rung by rung, Farkas rays included; the record's separation (product nothing, faces 19 on L1 and the ODE plant, faces 10 on L2) was the stock product against the faces, and the box set now reaches what the faces reach at 75% of their variables. The scalar-rule faces (today's 'auto') are the weakest set at w = [1 1] on all three plants and tie at w = [2 2]. No uncertain verdict appears for any set with a per-direction box or face term.

## Phase 4 (C3, io2 gain through the executive): inconclusive as first run

With dw pinned at 0 and 1 the executive's 'like' slack for the KYP block had Gram blocks of 38 and 149 rows (the tier's own level for the gain LPIs is dw_P = 2, and the record's certified 2-D KYP needed (2,2) at about 5e5 variables); every configuration, both storages, returned an uncertified gain (MOSEK unknown; objectives 1e-6 to 1e8). No conclusion about the term sets; rerun at dw 2 and 3 queued (c3b).

## Phase 4b (C5): 3-D heat DD x DN x DN, kappa* = 1.5 pi^2 = 14.804407, bench degrees (RQ 'Tspan', dp 0), d = 0, one certified solve each at kappa = 0.985 kappa* (run beside C2c: timings contended)

| cfg | verdict | nx | m | Ks (R, then Q) | solve s | peak MB |
|---|---|---|---|---|---|---|
| box o1 (per direction) | Farkas infeasible (ray violation 3.6e-12) | 583,019 | 14,471 | [288 144 96 96], [552 276 200 200] | 105 | 10,541 |
| product o1 | Farkas infeasible (6.8e-12) | 615,387 | 14,471 | [288 144 96 48 96 48 32 16], [552 276 200 100 200 100 72 36] | 99 | 13,386 |
| faces o1 (per direction) | Farkas infeasible (8.2e-14) | 778,363 | 14,471 | [288 144 x2 96 x4], [552 276 x2 200 x4] | 73 | 13,386 |

Reading. At the per-direction offset every set is certified infeasible at 0.985 kappa*, where the record's bench certified [14.60, 14.70) = [0.986, 0.993) with the faces at the full degree (offset 0, blocks 288 x 7 and 552 x 7, 2.71M variables). The per-direction terms are a smaller cone than the record's (T6: box at (w, w - 1) lies inside the faces at w), so this run is below the record's degree; the offset-0 configurations (faces0, box0, product0) follow in c5b.

## Phase 2 at rtol 1e-4 (C2c): the tied configurations re-bisected (24 solves, 1800 s each)

| cfg | QD | lo (reach) | hi | nsolve, s | nx, nnz |
|---|---|---|---|---|---|
| box o1 | (2,2) | 12.3367312 (0.99998) | uncertain to kappa* + tol | 20, 463 | 358,665, 1,123,343 |
| facesS | (2,2) | 12.3367312 (0.99998) | uncertain | 16, 321 | 339,209, 904,851 |
| faces o1 | (2,2) | 12.3367312 (0.99998) | uncertain | 16, 464 | 526,601, 1,502,341 |
| product o1 | (2,2) | 12.3367312 (0.99998) | uncertain | 18, 435 | 395,785, 1,402,452 |
| corners o1 | (2,2) | 12.3367312 (0.99998) | uncertain | 16, 591 | 675,081, 2,171,081 |
| faces0 | (2,1) | 12.3359773 (0.99992) | 12.3367312 (Farkas) | 17, 144 | 204,809, 630,999 |

Reading. At wQ = 2 the five sets certify the same kappa 12.3367312 (gap 2.74e-4 = the tolerance; above it every solve is uncertain for every set), so the Putinar box set, the 2^N set, the faces at either offset rule and the corners are indistinguishable in reach at 1e-4 on this benchmark; the record's separation of 10/08 (faces at w - 1 12.336731 against the tensor set 12.326178) does not reproduce: the box set reaches the faces' value, at 6% more variables than the scalar-rule faces and 32% fewer than the per-direction faces. The faces at offset 0 with wQ = 1 are the only configuration separated by a certificate: 12.3359773 feasible, 12.3367312 Farkas-infeasible, 7.5e-4 below the wQ = 2 value at 60% of the smallest wQ = 2 size. H1 holds on C2: the box set at per-direction offsets certifies at the same weight degree as the faces.

## Phase 4 rerun (C3b): io2 coercive H-infinity gain through the executive, dw 2 and 3, closed-form gain 0.1710999

gamma (gap = gamma/gex - 1), ndv of the whole program, time; verdict = the executive's certified flag.

| cfg | stock, dw 2 | stock, dw 3 | N-D default, dw 2 | N-D default, dw 3 |
|---|---|---|---|---|
| none | unknown, 111,621 | unknown, 351,714 | unknown, 111,701 | unknown, 351,794 |
| box o1 | 0.1717287 (+3.68e-3), 211,079, 18 s | 0.1717117 (+3.58e-3), 747,764, 99 s | 0.1717374, 211,159, 21 s | 0.1717117, 747,844, 107 s |
| product o1 | 0.1717197 (+3.62e-3), 233,280, 19 s | 0.1717119, 859,320, 114 s | 0.1717386, 233,360, 18 s | 0.1717115, 859,400, 132 s |
| faces o1 | 0.1717221 (+3.64e-3), 310,537, 14 s | 0.1717111, 1,143,814, 98 s | 0.1717195, 310,617, 21 s | 0.1717110, 1,143,894, 142 s |
| facesS | 0.1717221 (+3.64e-3), 200,425, 13 s | 0.1717116, 797,938, 71 s | 0.1717299, 200,505, 14 s | 0.1717119, 798,018, 82 s |
| corners o1 | 0.1717218 (+3.63e-3), 399,341, 35 s | 0.1717111, 1,590,038, 182 s | 0.1717216, 399,421, 28 s | 0.1717112, 1,590,118, 164 s |

Reading. Every set with a box or face term certifies at dw 2, within 2e-5 of each other in gamma (1e-4 relative) and all 0.36% above the closed form at dw 3, where they agree to 1e-6: the +0.36% floor of the record belongs to the storage, not to the term set. The record's sizes reproduce exactly (facesS 200,425 = 'facesm1'; product 233,280 = the tensor set). Size order at dw 2: facesS < box (+5%) < product (+16%) < faces (+55%) < corners (+99%). H1 holds on the KYP slack: the Putinar box set certifies the faces' gain at the same weight degree.

## Phase 5 (C6): the 2-D library rows at light, stock storage, the tier's own raise loop

Stability 'P' on 33, 34, 35, 37, 39, 40; gain 'P' on 45, 46, 47; the estimator on 48. ndv = whole program.

| row | facesS (today's 'auto') | box o1 | faces o1 | product o1 |
|---|---|---|---|---|
| 33 | cert dw 0, 25,984 | cert dw 0, 31,168 | cert dw 0, 41,536 | cert dw 0, 32,464 |
| 34 | cert dw 1, 187,984 | cert dw 1, 198,352 | cert dw 1, 291,664 | cert dw 1, 219,088 |
| 35 | cert dw 0, 25,984 | cert dw 0, 31,168 | cert dw 0, 41,536 | cert dw 0, 32,464 |
| 37 | cert dw 1, 187,984 | cert dw 0, 31,168 | cert dw 0, 41,536 | cert dw 0, 32,464 |
| 39 | cert dw 0, 233,856, 87 s | cert dw 0, 280,512, 97 s | cert dw 0, 373,824, 112 s | cert dw 0, 292,176, 100 s |
| 40 | cert dw 0, 26,582 | cert dw 0, 31,764 | cert dw 0, 42,422 | cert dw 0, 33,133 |
| 45 gain | 0.1119816923, 200,425 | 0.1119846861, 211,079 | 0.1119811974, 310,537 | 0.1119759125, 233,280 |
| 46 gain | 0.1717221120, 200,425 | 0.1717287394, 211,079 | 0.1717220536, 310,537 | 0.1717196514, 233,280 |
| 47 gain | 0.2085352276, 200,425 | 0.2085350978, 211,079 | 0.2085348300, 310,537 | 0.2085347993, 233,280 |
| 48 estimator | 0.006179431421 cert dw 0, 195,815 | unknown at dw 0 and after the raise (dw 1, 717,207, 110 s) | 0.00464656466 cert dw 0, 301,799 | 0.006456606518 cert dw 0, 226,903 |

Reading. On the nine stability and gain rows the four sets certify alike; gains agree to 1e-5 relative; the box set is 5-20% larger than the scalar-rule faces and 25-30% smaller than the per-direction faces; on row 37 the scalar-rule faces needed one raise where the three per-direction sets certified at dw 0 (the w = 1 weakness of the scalar rule again). Row 48, the estimator with the substantive gap (record: 0.0061794 at light, best 0.0024154 at 2.41M), is the one case where the term set changes the outcome: the per-direction faces give 0.00465, 25% below today's bound, the 2^N set 0.00646, and the Putinar box set fails (MOSEK unknown at dw 0 and at dw 1); the pinned-dw study follows (c7).

## Phase 5b (C7): the PDE 48 estimator with the slack weight pinned, stock storage (the N-D default storage at dw 0, 1 certifies only faces0 at dw 1: 0.0073087, 110,215)

gain bound (certified) or status; ndv of the whole program; the tier's own level for this LPI is dw 2.

| cfg | dw 0 | dw 1 | dw 2 |
|---|---|---|---|
| facesS (today's 'auto') | unknown, 8,615 | unknown, 28,919 | 0.006179431 cert, 195,815, 17 s |
| box o1 | unknown, 5,415 | unknown, 34,071 | unknown (objective 0.00781), 204,999, 22 s |
| box0 | unknown | unknown, 66,327 | 0.006281592 cert, 323,367, 33 s |
| product o1 | unknown, 7,015 | unknown, 35,671 | 0.006456607 cert, 226,903, 21 s |
| faces o1 | unknown, 8,615 | unknown, 45,623 | 0.004646565 cert, 301,799, 23 s |
| faces0 | unknown, 8,615 | 0.004505899 cert, 110,135, 8 s | 0.005861550 cert, 538,535, 51 s |
| corners o1 | unknown, 15,015 | unknown, 52,023 | 0.006252202 cert, 389,415, 26 s |

Reading. This is the one case of the run where the term set and the offset rule change the certified bound. The faces at offset 0 certify one weight degree earlier than every other set, with the best bound of the run at the smallest size (0.00451 at 110,135 against today's 0.00618 at 195,815: 27% lower at 56% of the variables); the per-direction faces are next (0.00465 at dw 2); the Putinar box set at per-direction offsets never certifies at dw <= 2 (MOSEK unknown with the objective at 0.0078), while the same set at offset 0 certifies 0.00628 at dw 2. The bound of faces0 at dw 2 (0.00586) is above its dw 1 value although the cones are nested, so the dw 2 'optimal' is a stalled near-optimal acceptance of the kind recorded in the memory note on MOSEK near-optimal flags; the dw 2 spread between sets (0.0046 to 0.0065) is therefore partly solver stall, and the comparison that stands is the certified-at-dw-1 result of faces0. The ladder at dw 3 follows (c7b).

## Phase 2, the half step (C2d): box0 and product0 with Q at (2,1), rtol 1e-4

| cfg at Q (2,1) | lo (reach) | hi (Farkas) | nx | m | solves, s |
|---|---|---|---|---|---|
| box o1 (per direction; from C2) | 12.2643306 (0.99411) | 12.2764630 | 61,449 | 2131 | 11, 32 |
| box0 | 12.3284394 (0.99931) | 12.3291932 | 122,889 | 4103 | 16, 136 |
| product0 | 12.3284394 (0.99931) | 12.3291932 | 163,849 | 4567 | 16, 200 |
| faces0 (from C2c) | 12.3359773 (0.99992) | 12.3367312 | 204,809 | 3117 | 17, 144 |
| any set at Q (2,2) (from C2c) | 12.3367312 (0.99998) | uncertain | 339,209 and up | 3943 | |

Reading. With the plain block held at w = 1, the certified reach orders box o1 < box0 = product0 < faces0 < every set at w = 2, each step separated by a Farkas certificate. The faces at offset 0 beat the box factors at offset 0 although the latter carry the higher degree (2w + 2 against 2w + 1 in the acting direction) and cost 40% fewer variables: on this benchmark the one-sided face weights do something the box factor at any offset does not at the same plain degree, and the box set needs w = 2 (1.75x the variables of faces0) to pass it. The g1 g2 term again adds nothing (product0 = box0).

PDE 48 at dw 3 (C7b, stock storage): faces0 0.006788 certified at 1,682,615 (228 s); facesS 0.006993 at 767,351; box o1 unknown at 717,207 (objective 0.00823); box0 0.008986 at 1,009,815. Every certified bound at dw 3 is above the same set's value at dw 1 or 2 although the cones are nested, so the executive's estimator solves above dw 1 are near-optimal stalls of MOSEK and the exact optimum of these programs is not known from this run; what stands is the certified 0.004506 of faces0 at dw 1 (110,135) and the absence of any certificate for the per-direction box set at dw 0 to 3 with the stock storage.

## Phase 3b (C4b): reach on S3 and S8 with the N-D default storage, dw 1 pinned, rungs 0.95 to 0.99844 of pi^2

facesS (92,212 / 335,104 variables), box o1 (108,088 / 353,536) and faces o1 (152,188 / 519,424) certify every rung on both cascades, S3 in 5-11 s and S8 in 22-57 s per solve. No set loses anything up to the last rung, 0.16% below the margin; the executive path does not separate the sets even at the margin, and the box set costs 17% more than the scalar-rule faces and 29% less than the per-direction faces.

## Phase 2, the half step on a symmetric box (C2e): heat DD x DD, kappa* = 2 pi^2 = 19.739209, Q at (2,1), rtol 1e-4 (every run stopped by the 24-solve budget)

| cfg | lo (reach) | hi | nx | m |
|---|---|---|---|---|
| box o1 | 19.6013145 (0.99301) | 19.6109595 (Farkas) | 61,449 | 2139 |
| box0 | 19.6929414 (0.99766) | uncertain above | 122,889 | 4111 |
| faces0 | 19.7315331 (0.99961) | uncertain above | 204,809 | 3125 |
| facesS | 16.3472571 (0.82816) | 19.4322950 (Farkas) | 51,209 | 2139 |

Reading. The same order as on DD x DN at the same plain degree: scalar-rule faces < per-direction box < offset-0 box < offset-0 faces, so the advantage of the offset-0 faces over the offset-0 box factors is not a Neumann effect; it holds with Dirichlet on every edge. The scalar rule at w = 1 is again the weakest configuration by far.

## The half step where the theory is exact (C2f): 1-D heat DD, kappa* = pi^2 = 9.8696044, rtol 1e-5

| cfg | Q at | lo (reach) | hi (Farkas) | nx | m |
|---|---|---|---|---|---|
| pair (box o1) | (2,1) | 9.8047900 (0.99343) | 9.8048653 | 328 | 49 |
| box0 | (2,1) | 9.8611934 (0.99915) | uncertain above | 523 | 71 |
| faces0 | (2,1) | 9.8689500 (0.99993) | 9.8690253 | 783 | 60 |
| pair (box o1) | (2,2) | 9.8695525 (0.99999) | 9.8696278 | 720 | 63 |

Reading. The Markov-Lukacs prediction is observed: the faces at offset 0 with the plain block at w lie strictly between the pair at w and the pair at w + 1, each step separated by a Farkas certificate. In one variable the pair at w + 1 is the better buy (720 variables, 0.99999) than the faces at w (783, 0.99993); the box factor at offset 0 sits in between in reach and below both in size.

## The half step on the pinned 2-D plant (C2f): DD x DD, w = [1 1], rungs 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85 of lam*

| set | verdicts | ndec |
|---|---|---|
| box o1, box0, product0, faces o1, corners o1 | +++---- (0.65 certified, Farkas from 0.70) | 31,104 / 62,208 / 82,944 / 41,472 / 46,656 |
| faces0, corners0 | +++++++ (0.85 certified; Farkas at 0.90 from Phase 1) | 103,680 / 186,624 |
| facesS | ------- | 25,920 |

Reading. At the lowest weight the faces at offset 0 reach at least 0.85 where every box configuration, the per-direction faces and the per-direction corners stop at 0.65 with Farkas certificates, including the box factors at offset 0, whose weighted terms carry a higher degree than the faces'. The one-sided face weights at the full degree are what certifies at w = 1; the corners add nothing to them (corners0 = faces0 at 1.8x the size).

## Phase 4b complete (C5b): 3-D heat at the bench caps, weighted terms at offset 0

| cfg | kappa | verdict | nx | m | blocks | solve s | peak MB |
|---|---|---|---|---|---|---|---|
| faces0 (the record's bench) | 0.985 kappa* | certified (rel_b 6.0e-7) | 2,713,563 | 29,774 | 288 x 7, 552 x 7 | 759 | 19,123 |
| box0 | 0.985 kappa* | Farkas infeasible (1.3e-11) | 1,550,619 | 45,077 | 288 x 4, 552 x 4 | 768 | 27,980 |
| product0 | 0.985 kappa* | Farkas infeasible (1.7e-11) | 3,101,211 | 71,349 | 288 x 8, 552 x 8 | 4602 | 71,274 (above the 64 GB of the machine) |
| box0 | 0.992 kappa* | Farkas infeasible (1.6e-12) | 1,550,619 | 45,077 | as above | 680 | |

Reading. In three variables at the bench caps the faces at the full degree are the only set that certifies 0.985 kappa*; the box factors and all their products, at either offset, are certified infeasible there, and their rows outnumber the faces' (45,077 and 71,349 against 29,774: the even-degree weighted terms reach kernel monomials outside the plain support). The box set would need the plain block at the next cap (heavy, 'Tspan' dp 1), whose size is read in c5c; the product set is the most expensive configuration of the run and certifies nothing here.

## 3-D, the box set at the next cap (C5c): built, not solved

box o1 with R at the bench caps and Q at the heavy caps ('Tspan' dp [0 1]): m 77,157, nx 29,554,587, blocks [288 144 96 96] and [3888 2448 2040 2040], built in 295 s; above the 4.5M gate (the 3.1M product0 solve needed 71 GB), so not solved on the workstation. The box set's route to the degree where the offset-0 faces certify 0.985 kappa* (2.71M variables) costs at least 11x the variables; whether it certifies there is a Sol question.

## Phase 5c (C6b): the library stability rows at the stripped tier (the same 2-D storage as light), five sets

Every set certifies every row (33, 34, 35, 37, 39, 40). Sizes repeat the light-tier ones (the stripped tier's 2-D programs coincide with light on these rows): facesS 25,984 to 233,856 at dw 0, box +20%, faces +60%, box0 2.4x, faces0 4x; row 34 needs dw 1 for every set; row 37 needs dw 1 for the scalar-rule faces only (187,984) where faces0 (103,744), box (31,168), box0 and faces certify at dw 0.

## Mixed-domain gain (C3c): S24, ODE to 2-D edge, exact gain 0.5, light, dw 2

| cfg | stock storage: gamma (gap), ndv, s | N-D default storage: gamma (gap), ndv, s |
|---|---|---|
| facesS | 0.5000214 (+4.3e-5), 674,256, 41 | 0.5000199 (+4.0e-5), 674,560, 36 |
| box o1 | 0.5000272 (+5.4e-5), 674,526, 39 | 0.5000176 (+3.5e-5), 674,830, 40 |
| faces o1 | 0.5000363 (+7.3e-5), 1,013,520, 52 | 0.5000089 (+1.8e-5), 1,013,824, 57 |
| faces0 | 0.5000269 (+5.4e-5), 1,676,496, 118 | 0.5000265 (+5.3e-5), 1,676,800, 115 |

Reading. Every set certifies within 7e-5 of the exact gain; the differences between sets are below the solver's resolution. On this mixed-domain registry the box set and the scalar-rule faces have the same size to 0.04%; the per-direction faces cost 1.5x and the offset-0 faces 2.5x for the same verdict.

## The lift caveat (C1c): pinned-P DD x DD with the lift at Dmin (dD 0) instead of Dmin + 1

With D = [1 1], every set (box o1, faces o1, faces0, facesS) certifies 0.10 lam* and is Farkas-infeasible at 0.50 and above at w = [1 1], [2 2] and [3 3] alike (facesS nothing at w = [1 1]; faces0 uncertain at 0.50, w = [3 3]). The lift raise of the reader (D = Dmin + 1, held fixed under D1) is necessary for every term set, and at the lower lift no set is distinguished: the box = faces result of Phase 1 does not depend on the lift raise.

## C5c on Sol (submitted 11/15 17:10, job 1621028)

The C5c program (box o1, dp [0 1], kappa = 0.985 kappa*; m 77,157, nx 29,554,587, Ks [288 144 96 96 3888 2448 2040 2040]) runs on Sol as job 1621028 (public, 64 cores, 320 GB, 10 h; D30). The driver rebuilds the program on Sol, asserts the workstation's m and nx, converts once (heatNd_sdp), solves with heatNd_solve (tight; MOSEK log on, 64 threads, max time 8.5 h), and falls back to the independent rows only if the verdict is uncertain within 6 h. Result: `c5c_sol_summary.txt` in this folder when the job lands; the log is `/scratch/mpeet/pietools/psatz_faces/jobs/c5c_box3d_1621028.out`, the record `psatz_faces/results/c5c_sol.mat` beside it.

## Mixed-domain diagnostics, 11/16 (S17, S18): `mixed_domain/`

Maintainer's request of 11/16 ("Run the S17 and S18 diagnostics"). Verdicts in `mixed_domain/S17_S18_findings.md`, the argument in `mixed_domain/S17_obstruction.md`, every row in `mixed_domain/m1_m3_summary.md` and `mixed_domain/m2_chk_summary.md`, decisions D32–D34. In one line each: S17 has no coercive PI Lyapunov operator at any degree (proven; MOSEK Farkas-infeasible once ep I <= P <= I fixes the scale at ep = 0.5 and 0.1, S3 feasible at the same ep; the executive's SeDuMi 'certified' exits and the 11/12 D4 certificate are flags-only artefacts); S18 certifies at 0.5 and 0.9 tau* (767k–791k variables, 64–138 s) once the storage carries the mixed basis operators of the 2-D space at kernel degree 2, which no tier builds; on both systems every term set (none, box, facesS, faces0, corners) gives the same verdict.

C5c on Sol, result (11/16 02:43): COMPLETED by MOSEK's 8.5 h limit after 10 iterations (1,314 s to 4,853 s each on a dense 77k Schur complement, 2.56e9 nonzeros, 6.2e14 flops per factorization); last iterate PFEAS 4.1e-5, DFEAS 8.0e-6, MU 8.0e-6, status UNKNOWN; heatNd_solve st 0, rel_b 4.1e-3, PSD 3.9e-10. 9 h 41 min, 162 GB, about 1,394 CHE. Not settled; details and the extrapolation in `c5c_sol_summary.txt`.

Mixed-domain parameterization, 11/16 02:10-04:00 (`mixed_domain/parameterization.md`): the slack reader certifies every fixed positive target (T'T, T'E1T, T'E2T, A'A of S17, S3, S18; the S18 storage); T'T is outside every separable storage and inside the split-kernel one; the storage-based Q form fails S3 in 2-D; the direct PIE2PDE form (no storage) certifies S3, S17 (0.5 r2*, 435k variables, 22 s) and S18 (0.5 tau*, 228k, 11 s), infeasible at 1.1 of the margins and standing at eps 0.1. S17 is certifiable in the default notion; the obstruction holds for the 'P' form only.

## 3-D, the box set one weight degree up in the slack alone (C5d, Sol job 1706087)

3-D heat DD x DN x DN on [0,1]^3, kappa = 0.985 kappa* = 14.5823; box quadratics (s_k-a_k)(b_k-s_k) at per-direction offset 1, slack weight degree raised by one in R and Q (RQ 'custom', int [2 2 2], mult [1 1 1], subset bounds raised by the number of thetas); lift degrees and the free operator at the bench caps.

| configuration | verdict | nx | m | blocks | iterations | MOSEK s | threads | memory |
|---|---|---|---|---|---|---|---|---|
| C5d: box o1, weight +1 (Sol 1706087, sc034) | certified (OPTIMAL, st +1, rel_b 5.1e-7, PSD 3.8e-14) | 13,913,099 | 71,349 | R [1200 800 720 720], Q [2220 1480 1356 1356] | 20 | 18,235 (5.07 h) | 64 | 69.0 GB (MaxRSS 72.9 GB) |
| faces0 (record of 11/15, workstation) | certified (rel_b 6.0e-7) | 2,713,563 | 29,774 | 288 x 7, 552 x 7 | not recorded | 759 | 8 | 19 GB |
| box0 (record of 11/15, workstation) | Farkas infeasible (1.3e-11) | 1,550,619 | 45,077 | 288 x 4, 552 x 4 | not recorded | 768 | 8 | 28 GB |

Reading. The box quadratics with the slack weight degree raised by one, everything else at the bench caps, certify 0.985 kappa* on the 3-D heat benchmark (Sol job 1706087, 20 iterations, 18,235 s on 64 threads, 69.0 GB; build 101 s, job elapsed 5:07:30 from 2026-11-16 21:10:59, about 737 CHE). This is 5.1x the variables of faces0 at the lower weight degree and 24x its solve time, though the times are on different machines (Sol, 64 threads, against the workstation, 8 threads), while box0 at the lower degree is Farkas infeasible. Sol job 1621028 (C5c), which raised the free operator as well (29.6M variables, 2.1x those of C5d), did not finish in 9.5 h.
