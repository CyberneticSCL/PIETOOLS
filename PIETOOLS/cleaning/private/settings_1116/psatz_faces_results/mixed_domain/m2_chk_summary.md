# M2: SeDuMi pairing (light, N-D default storage, dw 1 pinned)
| system | parameter | storage | set | dw | status | cert | ndv | m | t (s) | numerr | feasratio | eta |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | def | facesS | 1 | inaccurate | 1 | 120292 | 2716 | 59 | 1 | -0.0356 | 1.3e-09 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | def | box | 1 | inaccurate | 1 | 141028 | 2716 | 84 | 1 | 0.0482 | 1.5e-10 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | def | facesS | 1 | infeasible | 0 | 92212 | 2556 | 27 | 0 | -1.01 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | def | box | 1 | infeasible | 0 | 108088 | 2556 | 33 | 0 | -1.01 | 1.0e+00 |

The executive counts a SeDuMi 'inaccurate' exit (numerr 1) as certified; both S17 rows are such exits, with feasibility ratios -0.036 and 0.048. The CHK tables below test those operators.

# CHK (batch 1, chk.log): the executive at light, dw 1, the tier's own term set ('auto'); then the Step-2 family

| system | storage | solver | status | cert | ndv | m | t (s) |
|---|---|---|---|---|---|---|---|
| S17 | def | sedumi | inaccurate | 1 | 415204 | 3715 | 313 |
| D4 | def | sedumi | infeasible | 0 | 415204 | 3717 | 163 |
| S3 | def | sedumi | inaccurate | 1 | 318004 | 3493 | 139 |
| D4 | stock | mosek | infeasible | 0 | 414916 | 3713 | 14 |
| S3 | stock | mosek | unknown | 0 | 317716 | 3484 | 14 |
| D4 | def | mosek | infeasible | 0 | 415204 | 3717 | 16 |
| S3 | def | mosek | optimal | 1 | 318004 | 3493 | 9 |

Step-2 family (e1 = (L_n, 0), e2 = (0, L_n psi''), psi = s2 (1-s2)^k; g(c) = a + 2bc + dc^2 must be <= 0): rows with n <= 8 (at n = 12 the polynomial arithmetic loses precision: V(e1) comes out negative, so those rows are not evidence).

| system | storage | solver | n | k | a | b | d | b^2 - ad | V(e1) | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| S17 | def | sedumi | 4 | 1 | -6.046e-01 | -8.101e-15 | -4.051e-04 | -2.449e-04 | 4.305e-01 | ok |
| S17 | def | sedumi | 4 | 3 | -6.046e-01 | -8.785e-16 | -1.587e-04 | -9.593e-05 | 4.305e-01 | ok |
| S17 | def | sedumi | 8 | 1 | -1.633e-01 | -1.649e-10 | -4.051e-04 | -6.617e-05 | 1.161e-01 | ok |
| S17 | def | sedumi | 8 | 3 | -1.633e-01 | 2.688e-10 | -1.585e-04 | -2.589e-05 | 1.161e-01 | ok |
| D4 | def | sedumi | 4 | 1 | -7.751e-05 | -3.120e-07 | -7.449e-06 | -5.773e-10 | 4.491e-07 | ok |
| D4 | def | sedumi | 4 | 3 | -7.751e-05 | 2.648e-07 | -1.009e-05 | -7.822e-10 | 4.491e-07 | ok |
| D4 | def | sedumi | 8 | 1 | -2.098e-05 | -8.428e-08 | -7.450e-06 | -1.563e-10 | 2.887e-08 | ok |
| D4 | def | sedumi | 8 | 3 | -2.098e-05 | 7.154e-08 | -1.019e-05 | -2.138e-10 | 2.887e-08 | ok |
| S3 | def | sedumi | 4 | 1 | -1.990e-05 | -1.282e-07 | -1.107e-06 | -2.200e-11 | 1.279e-07 | ok |
| S3 | def | sedumi | 4 | 3 | -1.990e-05 | -5.116e-08 | -1.476e-07 | -2.933e-12 | 1.279e-07 | ok |
| S3 | def | sedumi | 8 | 1 | -5.844e-06 | -8.470e-09 | -2.754e-07 | -1.609e-12 | 8.226e-09 | ok |
| S3 | def | sedumi | 8 | 3 | -5.844e-06 | -3.353e-09 | -3.315e-08 | -1.937e-13 | 8.226e-09 | ok |
| D4 | stock | mosek | 4 | 1 | -9.684e-05 | -3.120e-07 | -7.449e-06 | -7.213e-10 | 5.605e-07 | ok |
| D4 | stock | mosek | 4 | 3 | -9.684e-05 | 2.648e-07 | -1.009e-05 | -9.774e-10 | 5.605e-07 | ok |
| D4 | stock | mosek | 8 | 1 | -2.621e-05 | -8.428e-08 | -7.450e-06 | -1.953e-10 | 3.603e-08 | ok |
| D4 | stock | mosek | 8 | 3 | -2.621e-05 | 7.155e-08 | -1.007e-05 | -2.639e-10 | 3.603e-08 | ok |
| S3 | stock | mosek | 4 | 1 | -3.561e-02 | 5.379e-07 | -4.334e-07 | -1.543e-08 | 2.289e-04 | ok |
| S3 | stock | mosek | 4 | 3 | -3.561e-02 | 3.860e-08 | -5.779e-08 | -2.058e-09 | 2.289e-04 | ok |
| S3 | stock | mosek | 8 | 1 | -1.046e-02 | 1.534e-07 | -1.079e-07 | -1.128e-09 | 1.472e-05 | ok |
| S3 | stock | mosek | 8 | 3 | -1.046e-02 | 1.814e-08 | -1.321e-08 | -1.382e-10 | 1.472e-05 | ok |
| D4 | def | mosek | 4 | 1 | -9.904e-05 | -3.120e-07 | -7.449e-06 | -7.376e-10 | 5.731e-07 | ok |
| D4 | def | mosek | 4 | 3 | -9.904e-05 | 2.648e-07 | -1.009e-05 | -9.995e-10 | 5.731e-07 | ok |
| D4 | def | mosek | 8 | 1 | -2.681e-05 | -8.428e-08 | -7.449e-06 | -1.997e-10 | 3.684e-08 | ok |
| D4 | def | mosek | 8 | 3 | -2.681e-05 | 7.155e-08 | -1.011e-05 | -2.710e-10 | 3.684e-08 | ok |
| S3 | def | mosek | 4 | 1 | -1.896e-05 | -1.264e-07 | -1.101e-06 | -2.085e-11 | 1.219e-07 | ok |
| S3 | def | mosek | 4 | 3 | -1.896e-05 | -5.057e-08 | -1.467e-07 | -2.780e-12 | 1.219e-07 | ok |
| S3 | def | mosek | 8 | 1 | -5.570e-06 | -8.134e-09 | -2.739e-07 | -1.526e-12 | 7.838e-09 | ok |
| S3 | def | mosek | 8 | 3 | -5.570e-06 | -3.253e-09 | -3.552e-08 | -1.978e-13 | 7.838e-09 | ok |

# CHK2 (next.log): both families; D4 with the corner rows (the 11/12 record certified D4 at 664k variables with corner rows on by default; today both corner configurations are infeasible)

| system | storage | solver | set | dw | status | cert | ndv | m | t (s) |
|---|---|---|---|---|---|---|---|---|---|
| S17 | def | sedumi | auto | 1 | inaccurate | 1 | 415204 | 3715 | 309 |
| D4 | stock | mosek | corners0 | 1 | infeasible | 0 | 746692 | 3803 | 29 |
| D4 | stock | mosek | corners | 1 | infeasible | 0 | 235204 | 2714 | 7 |
| S3 | def | mosek | auto | 1 | optimal | 1 | 318004 | 3493 | 9 |

Family B (Step 1: x1 = g, x2 = g (1-s2)^K, g = s1 L_m or s1 (1-s1) L_m): f = <xh, -D xh> on the returned operator; rows with K m <= 64 (the K = 16, m = 8 rows have V < 0: precision lost, not evidence). For S17 f does not change with K: the SeDuMi operator carries no x2 block beyond the margin eppos = 1e-4, so the uncancelled trace term it leaves is of size eppos K and a violating state needs K of order 1/eppos, beyond polynomial reach; the two-sided-bound program M7 is the numerical test instead.

| system | storage | solver | K | m | f | V | f/V | verdict |
|---|---|---|---|---|---|---|---|---|
| S17 | def | sedumi | 4 | 0 | 1.828e+01 | 7.039e+00 | 2.597e+00 | ok |
| S17 | def | sedumi | 4 | 2 | 1.634e+02 | 1.436e+01 | 1.138e+01 | ok |
| S17 | def | sedumi | 4 | 4 | 2.623e+02 | 1.427e+01 | 1.838e+01 | ok |
| S17 | def | sedumi | 4 | 8 | 4.722e+02 | 1.428e+01 | 3.307e+01 | ok |
| S17 | def | sedumi | 8 | 0 | 1.828e+01 | 7.039e+00 | 2.597e+00 | ok |
| S17 | def | sedumi | 8 | 2 | 1.634e+02 | 1.436e+01 | 1.138e+01 | ok |
| S17 | def | sedumi | 8 | 4 | 2.623e+02 | 1.427e+01 | 1.838e+01 | ok |
| S17 | def | sedumi | 8 | 8 | 4.722e+02 | 1.428e+01 | 3.307e+01 | ok |
| S17 | def | sedumi | 16 | 0 | 1.828e+01 | 7.039e+00 | 2.597e+00 | ok |
| S17 | def | sedumi | 16 | 2 | 1.634e+02 | 1.436e+01 | 1.138e+01 | ok |
| S17 | def | sedumi | 16 | 4 | 2.623e+02 | 1.427e+01 | 1.838e+01 | ok |
| D4 | stock | mosek | 4 | 0 | 6.960e-03 | 3.482e-04 | 1.999e+01 | ok |
| D4 | stock | mosek | 4 | 2 | 2.304e-02 | 2.480e-04 | 9.290e+01 | ok |
| D4 | stock | mosek | 4 | 4 | 5.975e-02 | 2.435e-04 | 2.454e+02 | ok |
| D4 | stock | mosek | 4 | 8 | 1.962e-01 | 2.446e-04 | 8.022e+02 | ok |
| D4 | stock | mosek | 8 | 0 | 6.950e-03 | 3.480e-04 | 1.997e+01 | ok |
| D4 | stock | mosek | 8 | 2 | 2.304e-02 | 2.479e-04 | 9.292e+01 | ok |
| D4 | stock | mosek | 8 | 4 | 5.974e-02 | 2.434e-04 | 2.455e+02 | ok |
| D4 | stock | mosek | 8 | 8 | 1.955e-01 | 2.344e-04 | 8.337e+02 | ok |
| D4 | stock | mosek | 16 | 0 | 6.925e-03 | 3.480e-04 | 1.990e+01 | ok |
| D4 | stock | mosek | 16 | 2 | 2.302e-02 | 2.478e-04 | 9.287e+01 | ok |
| D4 | stock | mosek | 16 | 4 | 5.977e-02 | 2.437e-04 | 2.453e+02 | ok |
| D4 | stock | mosek | 4 | 0 | 5.549e-03 | 2.779e-04 | 1.997e+01 | ok |
| D4 | stock | mosek | 4 | 2 | 1.837e-02 | 1.980e-04 | 9.279e+01 | ok |
| D4 | stock | mosek | 4 | 4 | 4.765e-02 | 1.943e-04 | 2.452e+02 | ok |
| D4 | stock | mosek | 4 | 8 | 1.565e-01 | 1.951e-04 | 8.019e+02 | ok |
| D4 | stock | mosek | 8 | 0 | 5.538e-03 | 2.777e-04 | 1.994e+01 | ok |
| D4 | stock | mosek | 8 | 2 | 1.836e-02 | 1.979e-04 | 9.282e+01 | ok |
| D4 | stock | mosek | 8 | 4 | 4.764e-02 | 1.942e-04 | 2.454e+02 | ok |
| D4 | stock | mosek | 8 | 8 | 1.567e-01 | 2.012e-04 | 7.788e+02 | ok |
| D4 | stock | mosek | 16 | 0 | 5.513e-03 | 2.776e-04 | 1.986e+01 | ok |
| D4 | stock | mosek | 16 | 2 | 1.835e-02 | 1.978e-04 | 9.276e+01 | ok |
| D4 | stock | mosek | 16 | 4 | 4.765e-02 | 1.945e-04 | 2.450e+02 | ok |
| S3 | def | mosek | 4 | 0 | 2.609e-04 | 6.264e-05 | 4.165e+00 | ok |
| S3 | def | mosek | 4 | 2 | 3.412e-03 | 4.472e-05 | 7.630e+01 | ok |
| S3 | def | mosek | 4 | 4 | 1.001e-02 | 4.375e-05 | 2.287e+02 | ok |
| S3 | def | mosek | 4 | 8 | 3.443e-02 | 4.387e-05 | 7.847e+02 | ok |
| S3 | def | mosek | 8 | 0 | 5.951e-04 | 6.879e-05 | 8.650e+00 | ok |
| S3 | def | mosek | 8 | 2 | 3.968e-03 | 4.912e-05 | 8.077e+01 | ok |
| S3 | def | mosek | 8 | 4 | 1.120e-02 | 4.805e-05 | 2.332e+02 | ok |
| S3 | def | mosek | 8 | 8 | 6.084e-02 | 4.002e-05 | 1.520e+03 | ok |
| S3 | def | mosek | 16 | 0 | 1.233e-03 | 7.338e-05 | 1.681e+01 | ok |
| S3 | def | mosek | 16 | 2 | 4.659e-03 | 5.243e-05 | 8.886e+01 | ok |
| S3 | def | mosek | 16 | 4 | 1.152e-02 | 8.659e-05 | 1.331e+02 | ok |

# M7: the stability LPI with ep I <= P <= I (PF_BOX), light, dw 1, facesS, N-D default storage, MOSEK

| system | ep | status | ndv | t (s) | numerr | feasratio | relres | eta |
|---|---|---|---|---|---|---|---|---|
| S17_transport_edge_diffusion | 0.5 | infeasible | 220964 | 11 | 0 | -1 | 1.0e+00 | 1.0e+00 |
| S17_transport_edge_diffusion | 0.1 | infeasible | 220964 | 11 | 0 | -1 | 1.0e+00 | 1.0e+00 |
| S17_transport_edge_diffusion | 0.01 | unknown | 220964 | 12 | 2 | -1 | 1.9e+01 | 6.0e-09 |
| D4 | 0.5 | infeasible | 220964 | 16 | 0 | -1 | 1.0e+00 | 1.0e+00 |
| D4 | 0.1 | infeasible | 220964 | 11 | 0 | -1 | 1.0e+00 | 1.0e+00 |
| D4 | 0.01 | unknown | 220964 | 9 | 2 | -0.96 | 3.2e-02 | 5.2e-09 |
| S3_edge1d_heat2d | 0.5 | unknown | 192884 | 9 | 2 | -0.915 | 1.5e-03 | 4.6e-09 |
| S3_edge1d_heat2d | 0.1 | optimal | 192884 | 4 | 0 | 1 | 1.2e-07 | 4.0e-10 |
| S3_edge1d_heat2d | 0.01 | optimal | 192884 | 4 | 0 | 1 | 9.3e-08 | 1.2e-09 |

Reading. With the scale of P fixed, S17 and D4 are infeasible with a MOSEK infeasibility certificate at ep = 0.5 and 0.1 and UNKNOWN at ep = 0.01, where a violating state needs a boundary layer of width about 0.01; S3 (control) is feasible at ep = 0.1 and 0.01 (UNKNOWN at 0.5: a bound tighter than the condition number of its certificate). The executive's own margin eppos = 1e-4 is in the regime where the infeasibility is beyond the solver's resolution, which is why its exits on S17 are UNKNOWN or flags-only certificates.

# Q form (PIE2PDE as the executive builds it: storage P, free Q sized from P, T'*Q = P as an equality; the qform shadow lifts the 1-D-only restriction), light, facesS, MOSEK (next4.log)

| system | parameter | storage | dw | status | ndv | m | t (s) | numerr | eta |
|---|---|---|---|---|---|---|---|---|---|
| S3_edge1d_heat2d | {"r":8.8826439609804222} | def | 0 | infeasible | 63220 | 8647 | 6 | 0 | 1.0e+00 |
| S3_edge1d_heat2d | {"r":8.8826439609804222} | def | 1 | unknown | 441345 | 10369 | 28 | 2 | 2.1e-08 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | def | 0 | infeasible | 27112 | 4694 | 1 | 0 | 1.0e+00 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | def | 1 | infeasible | 165096 | 5938 | 7 | 0 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | def | 0 | infeasible | 21092 | 4725 | 1 | 0 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | def | 1 | infeasible | 126736 | 5956 | 5 | 0 | 1.0e+00 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | sep11 | 0 | unknown | 1193755 | 14095 | 157 | 2 | 5.2e-09 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | sep11 | 1 | unknown | 4717255 | 18998 | 1454 | 2 | 3.0e-09 |

Reading. With the uniform storage the Q form certifies nothing, S3 included (infeasible at dw 0, UNKNOWN at dw 1), although S3 certifies in the 'P' form with the same storage: the equality T'*Q = P cannot be met, because T'*Q carries the lower and upper integrals of T and Q while the uniform storage has one polynomial kernel on the whole square. The direct form of the hierarchy paper (QT symmetric and >= eps T*T as an inequality, no storage) is run next (qd.log).

# Storage-cone membership (CONE, next4.log): X = P with P from the storage builder, MOSEK; 'optimal' = X is in that storage's cone, 'infeasible' = a certificate that it is not

Variants: def = POSCOPVAR with DEG omitted (the uniform patterns of the tiers: the pure multiplier and the pure complete moment of a 2-D space; the 1-D space separable too); sep_11 = every separable basis operator (the four patterns over {evaluation, complete moment}); sep_00 = every basis operator with the lower and upper (incomplete) moments as well; light = mult 0 (2-D) / 1 (1-D), int 1, joint 2; raised = mult 1, int 2, joint 4.

| X | def light | def raised | sep_11 light | sep_11 raised | sep_00 raised |
|---|---|---|---|---|---|
| S17_TtT | infeasible | infeasible | infeasible | infeasible | optimal |
| S3__TtT | infeasible | infeasible | infeasible | infeasible | optimal |
| S18_TtT | infeasible | infeasible | infeasible | infeasible | optimal |
| S17x2_TtT | infeasible | infeasible | infeasible | infeasible | optimal |
| S18_Pstar | infeasible | infeasible | infeasible | optimal | optimal |
| S18_negD | infeasible | infeasible | infeasible | infeasible | infeasible |
| S3_Pstar | optimal | optimal | optimal | optimal | optimal |

Reading. The PDE-norm form T'*T of every system, and the x2-alone form (identity in s1, Volterra kernel in s2), lie in no separable storage: only the storage with the incomplete moments (split kernels) at the raised degrees contains them. The certified S18 storage P* lies in the four-pattern separable storage at the raised degrees and above, not in the uniform one at any degree; the certified S18 derivative operator lies in no storage (it is a slack object: the reader certifies it, fixed_operators.md). The certified S3 storage lies in every variant (control).

# Direct PIE2PDE form (QD, qd.log): Q free with the shape of T (LPIVAR_CDOPVAR, every cap at the degree shown), T'*Q = Q'*T as an equality, T'*Q >= 1e-4 T'*T and -(A'*Q + Q'*A) >= 0 certified by the slack reader from their supports (light, dw as shown); no storage; MOSEK

| system | parameter | Q degree | dw | status | ndv | t (s) | numerr | relres | eta |
|---|---|---|---|---|---|---|---|---|---|
| S3_edge1d_heat2d | {"r":8.8826439609804222} | 1 | 0 | optimal | 383270 | 15 | 0 | 1.4e-07 | 8.2e-10 |
| S3_edge1d_heat2d | {"r":8.8826439609804222} | 1 | 1 | optimal | 694555 | 33 | 0 | 2.9e-07 | 6.5e-10 |
| S3_edge1d_heat2d | {"r":8.8826439609804222} | 2 | 0 | optimal | 1048988 | 74 | 0 | 2.2e-07 | 3.6e-10 |
| S3_edge1d_heat2d | {"r":8.8826439609804222} | 2 | 1 | optimal | 1900813 | 206 | 0 | 4.7e-07 | 8.0e-09 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | 1 | 0 | optimal | 435350 | 22 | 0 | 1.1e-04 | 2.8e-09 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | 1 | 1 | unknown | 788950 | 68 | 2 | 6.5e-05 | 6.4e-09 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | 2 | 0 | optimal | 1160588 | 82 | 0 | 1.7e-05 | 2.3e-09 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | 2 | 1 | optimal | 2103088 | 257 | 0 | 2.0e-05 | 3.0e-09 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | 1 | 0 | optimal | 227510 | 11 | 0 | 1.0e-06 | 5.5e-11 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | 1 | 1 | optimal | 665815 | 20 | 0 | 1.4e-06 | 9.7e-12 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | 2 | 0 | optimal | 940188 | 45 | 0 | 6.0e-07 | 1.7e-10 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | 2 | 1 | optimal | 1703613 | 251 | 0 | 2.0e-05 | 1.3e-09 |

# Direct PIE2PDE form, soundness (QD2, qd2.log): Q degree 1, dw 0; the margin eppos in the coercivity T'*Q >= eppos T'*T

| system | parameter | eppos | status | ndv | t (s) | numerr | eta |
|---|---|---|---|---|---|---|---|
| S17_transport_edge_diffusion | {"r2":10.856564841198294} | 1e-4 | infeasible | 435350 | 20 | 0 | 1.8e-02 |
| S17_transport_edge_diffusion | {"r2":8.8826439609804222} | 1e-4 | unknown | 435350 | 17 | 2 | 8.6e-09 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | 0.1 | optimal | 435350 | 24 | 0 | 2.9e-09 |
| S18_lib43_interior_delay | {"tau":1.1381986000000002,"sgn":-1} | 1e-4 | infeasible | 227510 | 11 | 0 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.9312534,"sgn":-1} | 1e-4 | unknown | 227510 | 10 | 2 | 5.2e-10 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | 0.1 | optimal | 227510 | 16 | 0 | 6.4e-11 |

Reading. At 1.1 of the margin (unstable) both systems are infeasible; at 0.5 of the margin both certify with eppos 0.1 as with 1e-4, which the 'P'-form artefacts of S17 did not survive (M7); at 0.9 of the margin Q of degree 1 is UNKNOWN for both (S18 reaches 0.9 tau* in the 'P' form with the raised four-pattern storage).

# The PDE-norm lift gam I - T'*T (POINC, poinc.log): lam = Rayleigh lower bound of ||T||^2 (Legendre products orthonormal on [0,1]: degree 6 on a 1-D space, 4 per variable on a 2-D space); bisection of 8 steps on gam in [lam, 4 lam] (resolution 1.2% of lam); reach = smallest certified gam over lam. The S18 rows of this pass used the unit interval on a [0,pi] x [0,tau] domain and are discarded (the script now reads the PIE's domain; the sharp-constant run below).

| PIE | lam | F = none, gam up to 4 lam | faces, light | faces, raised (dD 2, dw 2) |
|---|---|---|---|---|
| S17 | 0.541992 | unknown at 4 lam | reach 1.01172 | reach 1.01172 |
| S3 | 0.013704 | infeasible at 4 lam | reach 1.01172 | reach 1.01172 |
| S17x2 | 0.010266 | infeasible at 4 lam | reach 1.01172 | reach 1.01172 |
| D1_transport_1d | 0.405285 | infeasible at 4 lam | reach 1.01172 | reach 1.01172 |

Reading. lam equals the exact norm where it is known (1-D transport 4/pi^2 = 0.405285; the 2-D state of order (0,2) alone 1/pi^4 = 0.010266). Without a term the lift of the PDE norm is infeasible with a certificate at four times the norm on S3, on the 2-D state alone and on the 1-D transport, and UNKNOWN on S17; with the faces (the box factor in 1-D) it is certified to the bisection's resolution above lam at the light lengths, and the raised lengths change nothing at this resolution.

# The sharp constant (POINC2, poinc2.log): bracket [lam, 1.1 lam], 12 steps (resolution 2.4e-5 lam), the reader's default terms (faces; the box factor in 1-D), Legendre products on the PIE's own domain

| PIE | domain | lam | faces, light: gam_min/lam | faces, raised: gam_min/lam |
|---|---|---|---|---|
| S17 | [0 1;0 1] | 0.541992 | 1.00002 | 1.00002 |
| S3 | [0 1;0 1] | 0.013704 | 1.00002 | 1.00002 |
| S18 | [0 3.14159265358979;0 0.517363] | 1.549369 | 1.00002 | 1.00002 |
| D1_transport_1d | [0 1] | 0.405285 | 1.00002 | 1.00002 |

Reading. The certificate reaches the sharp constant to the resolution of the bisection on every PIE at the light lengths; the raised lengths change nothing. lam is a lower bound of ||T||^2 and the certified gam an upper bound, so ||T||^2 is pinned to five digits on each system.
