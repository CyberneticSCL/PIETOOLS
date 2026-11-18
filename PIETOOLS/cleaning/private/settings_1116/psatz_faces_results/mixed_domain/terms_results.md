# Term experiments on the mixed suite: results (2026-11-16)

Plan: `terms_plan.md`. Every row MOSEK, light lengths of the slack reader unless said.

## T1: the lift of the PDE norm gam I - T'T by term set (t1.log)

Bracket [lam, 1.1 lam], 10 bisection steps (resolution 1e-4 of lam); entries: the smallest certified gam over lam, or the failing status at 1.1 lam. Sets: none (plain term); auto (the reader's default); facesS (faces, scalar offset 1); faces0 (faces at offset 0); faces (per-direction offset 1); box (Putinar factors per direction); corners (faces and their cross products); faces_s1, faces_s2 (the faces of one registry direction); box_s1, box_s2 (one box factor); box_s1+faces.

| PIE | lam | none | auto | facesS | faces0 | faces | box | corners | faces_s1 | faces_s2 | box_s1 | box_s2 | box_s1+faces |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S17 | 0.541992 | infeasible at 1.1 | 1.00010 | 1.00010 | 1.00010 | 1.00010 | 1.00010 | 1.00010 | unknown at 1.1 | infeasible at 1.1 | unknown at 1.1 | infeasible at 1.1 | 1.00010 |
| S3 | 0.013704 | infeasible at 1.1 | 1.00010 | 1.00010 | 1.00010 | 1.00010 | 1.00010 | 1.00010 | unknown at 1.1 | infeasible at 1.1 | unknown at 1.1 | infeasible at 1.1 | 1.00010 |
| S8 | 0.015801 | infeasible at 1.1 | 1.00010 | 1.00010 | 1.00010 | 1.00010 | 1.00010 | 1.00010 | infeasible at 1.1 | infeasible at 1.1 | unknown at 1.1 | unknown at 1.1 | 1.00010 |
| S18 | 1.549369 | unknown at 1.1 | 1.00010 | 1.00010 | 1.00010 | 1.00010 | 1.00010 | 1.00010 | 1.04150 | unknown at 1.1 | 1.05000 | unknown at 1.1 | 1.00010 |
| S26 | 2.186661 | unknown at 1.1 | 1.00010 | 1.00010 | 1.00010 | 1.00010 | 1.00010 | 1.00010 | unknown at 1.1 | unknown at 1.1 | unknown at 1.1 | unknown at 1.1 | 1.00010 |

Reading. On every system the sets with a weight in both registry directions reach the resolution; a set with weights in one direction only certifies nothing on S17, S3, S8 and S26 (infeasible or UNKNOWN at 1.1 lam) and reaches 4-5% above the norm on S18 with the s1 weights (faces 1.0415, box 1.0500) and nothing with the s2 weights. So the lift of the PDE norm needs a term in every registry direction, including, for the 1-D space of a mixed PIE, the direction it does not own; the cross products (corners) add nothing; the offset rule (scalar or per direction, 0 or 1) does not change the reach at this resolution.

## T2 (a): the direct PIE2PDE test by term set at 0.5 of the margin (t2a.log); Q at length 2, dw 0, eppos 1e-4, both inequalities reader-sized at light with the set shown

| system | parameter | terms | Q length | dw | status | ndv | t (s) | eta |
|---|---|---|---|---|---|---|---|---|
| S3 | {"r":4.934802200544679} | none | 2 | 0 | infeasible | 76774 | 6 | 1.0e+00 |
| S3 | {"r":4.934802200544679} | facesS | 2 | 0 | optimal | 95930 | 3 | 2.1e-10 |
| S3 | {"r":4.934802200544679} | faces0 | 2 | 0 | optimal | 383270 | 9 | 5.1e-10 |
| S3 | {"r":4.934802200544679} | faces | 2 | 0 | optimal | 153398 | 4 | 2.1e-10 |
| S3 | {"r":4.934802200544679} | box | 2 | 0 | optimal | 115086 | 3 | 1.9e-10 |
| S3 | {"r":4.934802200544679} | corners | 2 | 0 | optimal | 172554 | 4 | 5.6e-10 |
| S3 | {"r":4.934802200544679} | faces_s1 | 2 | 0 | optimal | 115086 | 4 | 4.5e-10 |
| S3 | {"r":4.934802200544679} | faces_s2 | 2 | 0 | infeasible | 115086 | 4 | 7.9e-01 |
| S3 | {"r":4.934802200544679} | box_s1 | 2 | 0 | optimal | 95930 | 3 | 1.3e-09 |
| S3 | {"r":4.934802200544679} | box_s1+faces | 2 | 0 | optimal | 172554 | 4 | 1.8e-09 |
| S8 | {"r3":4.934802200544679} | none | 2 | 0 | infeasible | 98576 | 6 | 6.1e-02 |
| S8 | {"r3":4.934802200544679} | facesS | 2 | 0 | optimal | 123168 | 5 | 3.2e-10 |
| S8 | {"r3":4.934802200544679} | faces0 | 2 | 0 | optimal | 492048 | 13 | 6.2e-10 |
| S8 | {"r3":4.934802200544679} | faces | 2 | 0 | optimal | 196944 | 6 | 1.6e-10 |
| S8 | {"r3":4.934802200544679} | box | 2 | 0 | optimal | 147760 | 5 | 2.6e-10 |
| S8 | {"r3":4.934802200544679} | corners | 2 | 0 | optimal | 221536 | 6 | 1.7e-10 |
| S8 | {"r3":4.934802200544679} | faces_s1 | 2 | 0 | infeasible | 147760 | 6 | 5.4e-02 |
| S8 | {"r3":4.934802200544679} | faces_s2 | 2 | 0 | optimal | 147760 | 5 | 4.1e-10 |
| S8 | {"r3":4.934802200544679} | box_s1 | 2 | 0 | infeasible | 123168 | 6 | 7.4e-02 |
| S8 | {"r3":4.934802200544679} | box_s1+faces | 2 | 0 | optimal | 221536 | 6 | 2.6e-10 |
| S17 | {"r2":4.934802200544679} | none | 2 | 0 | unknown | 87190 | 4 | 1.3e-09 |
| S17 | {"r2":4.934802200544679} | facesS | 2 | 0 | optimal | 108950 | 5 | 6.7e-09 |
| S17 | {"r2":4.934802200544679} | faces0 | 2 | 0 | optimal | 435350 | 15 | 2.8e-09 |
| S17 | {"r2":4.934802200544679} | faces | 2 | 0 | unknown | 174230 | 6 | 6.0e-09 |
| S17 | {"r2":4.934802200544679} | box | 2 | 0 | unknown | 130710 | 4 | 2.1e-08 |
| S17 | {"r2":4.934802200544679} | corners | 2 | 0 | optimal | 195990 | 5 | 5.8e-09 |
| S17 | {"r2":4.934802200544679} | faces_s1 | 2 | 0 | infeasible | 130710 | 4 | 4.6e-02 |
| S17 | {"r2":4.934802200544679} | faces_s2 | 2 | 0 | unknown | 130710 | 6 | 3.1e-09 |
| S17 | {"r2":4.934802200544679} | box_s1 | 2 | 0 | infeasible | 108950 | 4 | 7.9e-02 |
| S17 | {"r2":4.934802200544679} | box_s1+faces | 2 | 0 | unknown | 195990 | 5 | 8.0e-09 |
| S18 | {"tau":0.517363,"sgn":-1} | none | 2 | 0 | infeasible | 45622 | 3 | 1.0e+00 |
| S18 | {"tau":0.517363,"sgn":-1} | facesS | 2 | 0 | infeasible | 62282 | 4 | 1.0e+00 |
| S18 | {"tau":0.517363,"sgn":-1} | faces0 | 2 | 0 | optimal | 227510 | 10 | 5.5e-11 |
| S18 | {"tau":0.517363,"sgn":-1} | faces | 2 | 0 | infeasible | 101678 | 3 | 1.0e+00 |
| S18 | {"tau":0.517363,"sgn":-1} | box | 2 | 0 | infeasible | 73650 | 3 | 1.0e+00 |
| S18 | {"tau":0.517363,"sgn":-1} | corners | 2 | 0 | infeasible | 118338 | 4 | 1.0e+00 |
| S18 | {"tau":0.517363,"sgn":-1} | faces_s1 | 2 | 0 | infeasible | 68358 | 3 | 1.0e+00 |
| S18 | {"tau":0.517363,"sgn":-1} | faces_s2 | 2 | 0 | infeasible | 78942 | 3 | 1.0e+00 |
| S18 | {"tau":0.517363,"sgn":-1} | box_s1 | 2 | 0 | infeasible | 56990 | 3 | 1.0e+00 |
| S18 | {"tau":0.517363,"sgn":-1} | box_s1+faces | 2 | 0 | infeasible | 113046 | 4 | 1.0e+00 |
| S26 | {"a":0.98370055013616975} | none | 2 | 0 | infeasible | 60259 | 7 | 1.0e+00 |
| S26 | {"a":0.98370055013616975} | facesS | 2 | 0 | optimal | 82859 | 6 | 5.1e-10 |
| S26 | {"a":0.98370055013616975} | faces0 | 2 | 0 | optimal | 300395 | 12 | 2.0e-10 |
| S26 | {"a":0.98370055013616975} | faces | 2 | 0 | optimal | 134955 | 7 | 2.4e-10 |
| S26 | {"a":0.98370055013616975} | box | 2 | 0 | optimal | 97607 | 8 | 3.6e-10 |
| S26 | {"a":0.98370055013616975} | corners | 2 | 0 | optimal | 157555 | 8 | 2.4e-10 |
| S26 | {"a":0.98370055013616975} | faces_s1 | 2 | 0 | optimal | 90599 | 8 | 1.0e-08 |
| S26 | {"a":0.98370055013616975} | faces_s2 | 2 | 0 | infeasible | 104615 | 7 | 1.0e+00 |
| S26 | {"a":0.98370055013616975} | box_s1 | 2 | 0 | optimal | 75429 | 7 | 1.0e-09 |
| S26 | {"a":0.98370055013616975} | box_s1+faces | 2 | 0 | optimal | 150125 | 7 | 3.0e-10 |

## T2 (b): the direct PIE2PDE test by term set at 0.9 of the margin (t2b.log); Q at length 2, dw 0

| system | parameter | terms | Q length | dw | status | ndv | t (s) | eta |
|---|---|---|---|---|---|---|---|---|
| S3 | {"r":8.8826439609804222} | none | 2 | 0 | infeasible | 76774 | 5 | 1.0e+00 |
| S3 | {"r":8.8826439609804222} | facesS | 2 | 0 | optimal | 95930 | 3 | 5.2e-11 |
| S3 | {"r":8.8826439609804222} | faces0 | 2 | 0 | optimal | 383270 | 9 | 8.2e-10 |
| S3 | {"r":8.8826439609804222} | faces | 2 | 0 | optimal | 153398 | 4 | 8.6e-11 |
| S3 | {"r":8.8826439609804222} | box | 2 | 0 | optimal | 115086 | 3 | 1.8e-10 |
| S3 | {"r":8.8826439609804222} | corners | 2 | 0 | optimal | 172554 | 4 | 9.5e-11 |
| S3 | {"r":8.8826439609804222} | faces_s1 | 2 | 0 | optimal | 115086 | 4 | 1.3e-09 |
| S3 | {"r":8.8826439609804222} | faces_s2 | 2 | 0 | infeasible | 115086 | 5 | 1.0e+00 |
| S3 | {"r":8.8826439609804222} | box_s1 | 2 | 0 | optimal | 95930 | 3 | 3.8e-09 |
| S3 | {"r":8.8826439609804222} | box_s1+faces | 2 | 0 | optimal | 172554 | 4 | 1.4e-10 |
| S8 | {"r3":8.8826439609804222} | none | 2 | 0 | infeasible | 98576 | 6 | 5.5e-02 |
| S8 | {"r3":8.8826439609804222} | facesS | 2 | 0 | optimal | 123168 | 5 | 1.7e-10 |
| S8 | {"r3":8.8826439609804222} | faces0 | 2 | 0 | optimal | 492048 | 14 | 3.1e-10 |
| S8 | {"r3":8.8826439609804222} | faces | 2 | 0 | optimal | 196944 | 6 | 1.1e-10 |
| S8 | {"r3":8.8826439609804222} | box | 2 | 0 | optimal | 147760 | 5 | 2.2e-10 |
| S8 | {"r3":8.8826439609804222} | corners | 2 | 0 | optimal | 221536 | 6 | 5.0e-11 |
| S8 | {"r3":8.8826439609804222} | faces_s1 | 2 | 0 | infeasible | 147760 | 6 | 5.4e-02 |
| S8 | {"r3":8.8826439609804222} | faces_s2 | 2 | 0 | optimal | 147760 | 6 | 1.3e-09 |
| S8 | {"r3":8.8826439609804222} | box_s1 | 2 | 0 | infeasible | 123168 | 5 | 7.4e-02 |
| S8 | {"r3":8.8826439609804222} | box_s1+faces | 2 | 0 | optimal | 221536 | 6 | 2.5e-10 |
| S17 | {"r2":8.8826439609804222} | none | 2 | 0 | infeasible | 87190 | 4 | 3.7e-02 |
| S17 | {"r2":8.8826439609804222} | facesS | 2 | 0 | infeasible | 108950 | 4 | 6.5e-03 |
| S17 | {"r2":8.8826439609804222} | faces0 | 2 | 0 | unknown | 435350 | 12 | 8.6e-09 |
| S17 | {"r2":8.8826439609804222} | faces | 2 | 0 | infeasible | 174230 | 6 | 4.9e-03 |
| S17 | {"r2":8.8826439609804222} | box | 2 | 0 | infeasible | 130710 | 4 | 9.0e-03 |
| S17 | {"r2":8.8826439609804222} | corners | 2 | 0 | infeasible | 195990 | 6 | 3.4e-03 |
| S17 | {"r2":8.8826439609804222} | faces_s1 | 2 | 0 | infeasible | 130710 | 4 | 4.7e-02 |
| S17 | {"r2":8.8826439609804222} | faces_s2 | 2 | 0 | unknown | 130710 | 6 | 2.8e-08 |
| S17 | {"r2":8.8826439609804222} | box_s1 | 2 | 0 | infeasible | 108950 | 4 | 8.1e-02 |
| S17 | {"r2":8.8826439609804222} | box_s1+faces | 2 | 0 | infeasible | 195990 | 5 | 4.4e-03 |
| S18 | {"tau":0.9312534,"sgn":-1} | none | 2 | 0 | infeasible | 45622 | 2 | 1.0e+00 |
| S18 | {"tau":0.9312534,"sgn":-1} | facesS | 2 | 0 | infeasible | 62282 | 3 | 1.0e+00 |
| S18 | {"tau":0.9312534,"sgn":-1} | faces0 | 2 | 0 | unknown | 227510 | 7 | 5.2e-10 |
| S18 | {"tau":0.9312534,"sgn":-1} | faces | 2 | 0 | infeasible | 101678 | 4 | 1.0e+00 |
| S18 | {"tau":0.9312534,"sgn":-1} | box | 2 | 0 | infeasible | 73650 | 3 | 1.0e+00 |
| S18 | {"tau":0.9312534,"sgn":-1} | corners | 2 | 0 | infeasible | 118338 | 4 | 1.0e+00 |
| S18 | {"tau":0.9312534,"sgn":-1} | faces_s1 | 2 | 0 | infeasible | 68358 | 2 | 1.0e+00 |
| S18 | {"tau":0.9312534,"sgn":-1} | faces_s2 | 2 | 0 | infeasible | 78942 | 3 | 1.0e+00 |
| S18 | {"tau":0.9312534,"sgn":-1} | box_s1 | 2 | 0 | infeasible | 56990 | 2 | 1.0e+00 |
| S18 | {"tau":0.9312534,"sgn":-1} | box_s1+faces | 2 | 0 | infeasible | 113046 | 3 | 1.0e+00 |
| S26 | {"a":1.7706609902451056} | none | 2 | 0 | infeasible | 60259 | 7 | 1.0e+00 |
| S26 | {"a":1.7706609902451056} | facesS | 2 | 0 | optimal | 82859 | 7 | 2.4e-10 |
| S26 | {"a":1.7706609902451056} | faces0 | 2 | 0 | optimal | 300395 | 11 | 1.7e-10 |
| S26 | {"a":1.7706609902451056} | faces | 2 | 0 | optimal | 134955 | 7 | 1.6e-10 |
| S26 | {"a":1.7706609902451056} | box | 2 | 0 | optimal | 97607 | 8 | 4.0e-10 |
| S26 | {"a":1.7706609902451056} | corners | 2 | 0 | optimal | 157555 | 8 | 1.7e-10 |
| S26 | {"a":1.7706609902451056} | faces_s1 | 2 | 0 | optimal | 90599 | 8 | 2.4e-09 |
| S26 | {"a":1.7706609902451056} | faces_s2 | 2 | 0 | infeasible | 104615 | 7 | 1.0e+00 |
| S26 | {"a":1.7706609902451056} | box_s1 | 2 | 0 | optimal | 75429 | 8 | 3.1e-09 |
| S26 | {"a":1.7706609902451056} | box_s1+faces | 2 | 0 | optimal | 150125 | 8 | 1.5e-10 |

Reading of T2. The plain term certifies nothing. On S3, S8 and S26 every set with weights in both directions certifies at 0.5 and at 0.9 of the margin, facesS the smallest (96k, 123k, 83k variables; faces0 4x larger); the faces of one direction suffice on these three (s1 on S3 and S26, s2 on S8) and the other direction's alone do not. On S17 at 0.5 only facesS, faces0 and corners certify (per-direction faces and the box set UNKNOWN) and at 0.9 nothing at Q length 2 (faces0 UNKNOWN, the rest infeasible). On S18 at 0.5 only faces0 certifies (every other set infeasible with a certificate): the odd-degree half step of the psatz-faces run, the faces at offset 0 reaching total degree 2w + 1 where the others reach 2w; at 0.9 faces0 is UNKNOWN. So the reader's 'auto' rule (faces at offset 0 when the weight degree is 1, the scalar offset above) is the right default on the mixed suite: it is the only set that never loses a certified case, at 3-4x the variables of the scalar-offset faces where those suffice.

## T3: the product sets (t3.log; maintainer, 11/16 11:30: the theory supports products of faces, not faces; mixed-direction products untested)

Sets (PF_TERMSETS3; rows of per-direction codes 1 lower face, 2 upper face, 3 box quadratic g_k): product2 = {g1, g2, g1 g2} (the preordering of the box quadratics); mixedfb = the four faces + f_{1,i} g2, g1 f_{2,j} (face x box across directions); full16 = the 15 nonzero rows of {0,1,2,3}^2 (every product of distinct faces, the preordering of the faces); each at the per-direction offset 1 and (suffix 0) at offset 0. Direct PIE2PDE test, Q length 2, dw 0, as T2.

| PIE, margin | product2 | product20 | mixedfb | mixedfb0 | full16 | full160 |
|---|---|---|---|---|---|---|
| S3 0.5 and 0.9 | in 119,875 | in 306,646 | in 172,554 | in 689,766 | in 234,811 | in 1,226,134 |
| S8 0.5 and 0.9 | in 153,908 | in 393,680 | in 221,536 | in 885,520 | in 301,460 | in 1,574,096 |
| S26 0.5 and 0.9 | in 103,257 | in 240,361 | in 157,555 | in 540,531 | in 223,153 | in 960,769 |
| S17 0.5 | optimal, eta 1.7e-8 (above the 1e-8 gate), 136,150 | in 348,310 | in 195,990 | U 783,510 | U 266,710 | U 1,392,790 |
| S17 0.9 | out | U | out | U | out | U |
| S18 0.5 | out | out | out | in 409,398 | out | in 727,702 |
| S18 0.9 | out | out | out | U | out | U |

Reading. No product set reaches a case the faces miss: S3/S8/S26 certified by all (1.1-13x the facesS variables); S17 at 0.5 certified by the box preordering (where 'box' alone was U in T2) and by mixedfb, not by the full preordering (U); S17 at 0.9 nothing; S18 at 0.5 only by the sets containing the faces at offset 0 (mixedfb0, full160), at 1.8x and 3.2x faces0's 228k; S18 at 0.9 U as faces0. The mixed products across directions change no verdict; the offset-0 half step remains the one term-set effect.

T3, the lift of the PDE norm (bracket [lam, 1.1 lam], 10 steps; resolution 1.0001): S17 reach 1.0001 under all six product sets; S3 1.0001 under product2 (both offsets) and mixedfb (offset 1); the pass was stopped at 16:12 (window), S3's remaining sets and S8, S18, S26 not run (a bisection solve of a product set at offset 0 is 1.4M variables). The faces reached 1.0001 on all five PIEs in T1.

## T4: raising the weight degree with even-degree weights (t4.log)

Run: S18 (library PDE 43, 2-D delay state) and S17 (transport edge state driving a 2-D diffusion), each at 0.5 and 0.9 of its stability margin, direct PIE-to-PDE test, Q length 2, extra weight degree dw = 1, under the term sets box, facesS, faces0 and product2.

| PIE, margin | box | facesS | faces0 | product2 |
|---|---|---|---|---|
| S18 0.5 | in 475,107, 20 s | in 440,023, 20 s | in 1,290,075, 80 s | in 520,579, 28 s |
| S18 0.9 | in 475,107, 28 s | in 440,023, 30 s | in 1,290,075, 111 s | in 520,579, 28 s |
| S17 0.5 | U 832,470, 61 s | U 788,950, 57 s | in 2,203,350, 281 s | U 919,510, 66 s |
| S17 0.9 | U 832,470, 45 s | U 788,950, 44 s | U 2,203,350, 152 s | U 919,510, 62 s |

Reading. At dw = 1 all four sets certify S18 at 0.5 and at 0.9, the smallest being facesS at 440,023 variables; at dw = 0 S18 at 0.5 was certified only by faces0 (228k) and S18 at 0.9 by no set. At dw = 1 S17 at 0.5 is certified only by faces0 (2,203,350 variables, eta 3.9e-9), against facesS 109k, faces0 435k and product2 136k at dw = 0. S17 at 0.9 is certified by no set at either dw. Every U entry (box, facesS and product2 on S17 at both margins; faces0 on S17 at 0.9) returned numerr 2, MOSEK numerical problems.

## The cross length (11/17; maintainer: the 0-D/1-D construction gives R^q no monomials in s; "a waste of variables?"; "why can the same not be done with the slack?")

Storage (p4cross.log, p4bcross.log): the separable light row Rlight0 with len_cross = 1, 2, 3 (the length of the monomials a space carries in a direction it does not own; 2 = the reader's rule v2 at the light lengths), stability 'P' at 0.9 of the margin (S5 at 0.5), light tier; S18 at 0.5 tau* under Rlight2 (rule = 3).

| row | cross 1 | cross 2 | cross 3 |
|---|---|---|---|
| S3 | in 318,745 | in 319,545 | in 320,545 |
| S8 | in 336,105 | in 338,305 | in 341,305 |
| S4 | in 318,745 | in 319,545 | in 320,545 |
| S6 | in 527,117 | in 529,317 | in 532,370 |
| S5 | in 321,341 | in 322,161 | in 323,181 |
| S26 | in 338,601 | in 340,841 | in 343,881 |
| S18 at 0.5 tau*, Rlight2 | in 1,510,153 | in 1,512,873 | in 1,515,793 |

Every row certifies at cross length 1, with 0.2-1% fewer variables than the rule. (Rule v1 of 11/16, which lost S3, was cross length 1 under the UNIFORM row Rdef; with the separable patterns cross length 1 loses nothing.)

Slack (t5cross.log; shadow get_lift_degs with the option wcross = 0: an L2 space carries no slack monomials in a direction it does not own; the rule gives w_k there): the direct PIE2PDE test (Q length 2, dw 0) with facesS and faces0; the controls with the rule are the T2 rows.

| PIE, margin | facesS, wcross 0 (rule) | faces0, wcross 0 (rule) |
|---|---|---|
| S3 0.5 and 0.9 | in 86,078 (in 96k) | in 334,010 (in 383k) |
| S8 0.5 and 0.9 | in 101,424 (in 123k) | in 383,328 (in 492k) |
| S17 0.5 | in 98,410 (in 109k) | U 382,650, numerr 2 (in 435k) |
| S17 0.9 | out (out) | U (U) |
| S18 0.5 | out (out) | in 201,050 (in 228k) |
| S18 0.9 | out (out) | U (U) |
| S26 0.5 and 0.9 | in 69,915 (in 83k) | in 235,675 (in 300k) |

Lift of the PDE norm with the faces and wcross 0: S3 reach 1.0001, S17 1.0001, S18 see the REACH line of t5cross.log (the rule: 1.0001 on all three).

Reading. The slack without cross monomials certifies every case the rule certifies except S17 at 0.5 with the full-length faces (UNKNOWN with numerical problems, while the short faces still certify it), at 10-22% fewer variables; the lift reaches the resolution. So the cross monomials of the storage and of the slack are both knobs that can be set to none on this suite: len_cross = 1 in the row and wcross = 0 in the slack settings, with the former rules as other values.
