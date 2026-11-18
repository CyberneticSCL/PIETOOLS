# Variable-coefficient PDEs: stock light vs proposed row Rlight0, and slack term sets

Run 2026-11-16, 17:03 – 17:44, MATLAB R2025b, MOSEK via `pf_setup(false)`; one MATLAB process at a time; total MATLAB wall time about 41 min.
Executive: `PIETOOLS_stability_sop(PIE,st,'P')` through `pf_exec_row` (Lyapunov form P), tier `light` unless stated. `cert=1`: eta gate passed. `dw`: number of slack-weight raises the executive took. ndv: decision variables; m: SDP rows; t: seconds.

## Systems and reference thresholds

2-D systems on [0,1]^2, homogeneous Dirichlet on all four edges, lam a constant reaction coefficient.

| id | equation | lam* FD N=40 | lam* FD N=60 | Richardson (inferred) |
|---|---|---|---|---|
| A | x_t = (1+s1^2) x_s1s1 + (1+s2) x_s2s2 + lam x | 26.839174 | 26.846533 | 26.8526 |
| B | x_t = x_s1s1 + x_s2s2 + (lam - 3 s1 s2) x | 20.474475 | 20.479775 | 20.4841 |
| C | x_t = (1+s1) x_s1s1 + x_s2s2 + s2 x_s1 + lam x | 24.242392 | 24.248789 | 24.2541 |
| M2 (2-D block of M) | x_t = x_s1s1 + (1+s2^2) x_s2s2 + r2 x | 22.199946 | 22.205980 | 22.2110 |
| M1 (1-D block of M) | x_t = (1+s1) x_s1s1 + r1 x on [0,1], Dirichlet | 14.504004 | 14.507976 | 14.5113 |
| G7 (library PDE 7) | x_t = a x_ss + b x_s + (c0 + lam) x, x(0)=0, x_s(1)=0 | N=200: 4.653780 | N=400: eigs did not converge | – |

FD: second-order central differences (non-divergence form a(s) D2, central D1), N interior nodes per direction, h = 1/(N+1); lam* = - max Re eig of the discrete operator at lam = 0 (`eigs(L,4,'largestreal')`, `inh_fdmax.m`). G7: h = 1/N, ghost node at s = 1; a = s^3 - s^2 + 2, b = 3 s^2 - 2 s, c0 = -0.5 s^3 + 1.3 s^2 - 1.5 s + 0.7. The M1 run at N = 200 returned 57.65 with the warning "3 of the 4 requested eigenvalues converged" and is discarded; the G7 run at N = 400 returned NaN ("0 of the 4 requested eigenvalues converged"). The N = 60 value is the reference for every 2-D and M reach fraction. The FD values increase with N; the Richardson column (inferred, assumes O(h^2) error) lowers each reach by 0.02 – 0.03 %.

M (mixed cascade, variable-coefficient copy of S3): x1(s1) on [0,1], x2(s1,s2) on [0,1]^2,
x1_t = (1+s1) x1_s1s1 + r1 x1, x2_t = x2_s1s1 + (1+s2^2) x2_s2s2 + r2 x2; x1(0)=x1(1)=0, x2(0,s2)=x2(1,s2)=0, x2(s1,0)=x1(s1), x2(s1,1)=0. At fraction f: r1 = f*14.507976, r2 = f*22.205980 (f = 0.9: r1 = 13.057178, r2 = 19.985382). The threshold of the cascade in f is 1 (inferred from the triangular coupling; not computed).

## Item 1 – 2-D, stabP light at 0.9 lam*

| sys | lam | storage | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|
| A | 24.1619 | stock | optimal | 1 | 1 | 187984 | 2974 | 19 |
| A | 24.1619 | Rlight0 | optimal | 1 | 1 | 188545 | 3003 | 16 |
| B | 18.4318 | stock | optimal | 1 | 1 | 187984 | 2974 | 17 |
| B | 18.4318 | Rlight0 | optimal | 1 | 1 | 188545 | 3008 | 15 |
| C | 21.8239 | stock | optimal | 1 | 1 | 187984 | 2982 | 11 |
| C | 21.8239 | Rlight0 | optimal | 1 | 1 | 188545 | 3010 | 13 |

### Bisection on [0.5, 1] lam*, 8 steps, light

| sys | storage | certified fractions | not certified (status) | lam_cert | reach vs N=60 |
|---|---|---|---|---|---|
| A | stock | 0.75 (dw 0); 0.875, 0.9375, 0.96875, 0.98438, 0.98633 (dw 1) | 0.99219, 0.98828 (unknown, dw 1) | 26.479491 | 0.98633 |
| A | Rlight0 | same fractions, same dw | 0.99219, 0.98828 (unknown, dw 1) | 26.479491 | 0.98633 |
| B | stock | all 8 (0.75 dw 0; rest dw 1) | none | 20.439775 | 0.99805 (top of bracket) |
| B | Rlight0 | all 8, same dw | none | 20.439775 | 0.99805 (top of bracket) |
| C | stock | all 8 (0.75 dw 0; rest dw 1) | none | 24.201428 | 0.99805 (top of bracket) |
| C | Rlight0 | all 8, same dw | none | 24.201428 | 0.99805 (top of bracket) |

Program sizes: stock ndv 103744 at dw 0, 187984 at dw 1; Rlight0 ndv 104305 at dw 0, 188545 at dw 1. At fixed ndv, m varied by up to 11 rows between lam values (stock 2971 – 2982, Rlight0 3002 – 3013). t 3 – 26 s per call.

## Item 2 – mixed 1-D/2-D cascade M

stabP at f = 0.9 (inh_run_2), then the follow-up inh_run_5: stock heavy at f = 0.9 and an 8-step bisection of f on [0.5, 1] at light under both storages.

| f | storage | tier | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|
| 0.9 | stock | light | unknown | 0 | 1 | 508228 | 4757 | 42 |
| 0.9 | Rlight0 | light | optimal | 1 | 1 | 510057 | 5208 | 47 |
| 0.9 | stock | heavy | optimal | 1 | 0 | 1023444 | 6190 | 76 |

| storage, tier | certified f | not certified f (status) | reach |
|---|---|---|---|
| stock light | 0.625 (dw 0, ndv 317716); 0.6875, 0.71875, 0.73438, 0.73633 (dw 1, ndv 508228) | 0.75, 0.74219, 0.73828 (unknown, dw 1) | 0.73633 |
| Rlight0 light | 0.75, 0.875 (dw 0, ndv 319545); 0.9375, 0.96875, 0.97656 (dw 1, ndv 510057) | 0.98438, 0.98047, 0.97852 (unknown, dw 1) | 0.97656 |

Times per call: stock light 19 – 66 s, Rlight0 light 11 – 60 s.

Direct PIE-to-PDE test (`pf_qdirect_pie`: Q caps mult 1, int [1 1], out 1, in 1; eppos 1e-4, epneg 0) at f = 0.9:

| term set | dws | status | in/out/U | ndv | t | relres | eta |
|---|---|---|---|---|---|---|---|
| faces0 | 0 | optimal | in | 685770 | 23 | 2.8e-06 | 4.3e-11 |
| facesS | 0 | optimal | in | 192730 | 9 | 1.5e-04 | 1.4e-09 |
| box | 0 | optimal | in | 227011 | 10 | 5.1e-05 | 8.4e-10 |
| box | 1 | optimal | in | 1144267 | 90 | 3.8e-06 | 9.5e-11 |
| product2 | 1 | optimal | in | 1281391 | 139 | 3.9e-06 | 8.0e-11 |

## Item 3 – 1-D, library PDE 7 (Gahlawat 2017) and PDE 23 (Shivakumar 2019, Ex. 1)

PDE 7 built from the library function with params {'lam=...;', 'c = -0.5*s^3 + 1.3*s^2 - 1.5*s + 0.7 + lam;'}: the function evaluates c before it evaluates the params, so setting lam alone leaves c at lam = 4.65.

| lam | storage | tier | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|
| 4.0 | stock | light | optimal | 1 | 1 | 568 | 75 | 1 |
| 4.0 | Rlight0 | light | optimal | 1 | 1 | 568 | 75 | 0 |
| 4.5 | stock | light | optimal | 1 | 1 | 568 | 75 | 0 |
| 4.5 | Rlight0 | light | optimal | 1 | 1 | 568 | 75 | 0 |

Bisection of lam on [3, 4.66], 8 steps, stock:

| tier | certified lam | not certified | lam_cert | / 4.66 | / FD N=200 (4.653780) |
|---|---|---|---|---|---|
| light | all 8 (3.83 – 4.65352), dw 1, ndv 568, m 75 | none | 4.653516 | 0.99861 | 0.99994 |
| heavy | all 8, dw 0, ndv 1021, m 77 | none | 4.653516 | 0.99861 | 0.99994 |

PDE 23 at lam = 4.6, hinfP (container executive `PIETOOLS_Hinf_gain_sop`, P form), stock storage:

| tier | status | cert | dw | gam | ndv | m | t | library header (stock path) |
|---|---|---|---|---|---|---|---|---|
| light | optimal | 1 | 0 | 23.7056528 | 1177 | 98 | 1 | 25.249 |
| heavy | optimal | 1 | 0 | 23.7057986 | 1822 | 108 | 0 | 23.706 |

## Item 4 – direct test on A at 0.9 lam* (lam = 24.1619)

| term set | dws | status | in/out/U | ndv | t | relres | eta |
|---|---|---|---|---|---|---|---|
| faces0 | 0 | optimal | in | 518500 | 17 | 3.6e-07 | 2.3e-10 |
| facesS | 0 | optimal | in | 145828 | 6 | 5.8e-06 | 2.0e-09 |
| box | 1 | optimal | in | 864244 | 50 | 5.7e-07 | 2.0e-10 |
| product2 | 1 | optimal | in | 967924 | 65 | 3.2e-07 | 2.1e-10 |

## Reading

Measured.
- 2-D single domain (A, B, C): stock and Rlight0 returned the same verdict and the same dw at all 9 lam values per system, and the same certified threshold: 0.98633 lam* (A); 0.99805 lam* (B, C), the top of the bisection bracket, so the reach of B and C is not resolved above 0.99805. At dw 1 Rlight0 adds 561 decision variables (0.3 %) and 28 – 41 rows (same lam) to the stock program. Every certificate above 0.75 lam* took one weight raise (dw 1) under both storages.
- Mixed cascade M: the storages differ. At f = 0.9 stock light returned UNKNOWN after one raise and Rlight0 light certified. Bisection reach: stock light 0.73633, Rlight0 light 0.97656 of the FD threshold, at ndv 508228 vs 510057 (dw 1) and m 4758 vs 5223. Stock heavy certified f = 0.9 with ndv 1023444 (2.01 times the Rlight0 light program) in 76 s.
- All non-certified runs, under both storages and in every item, returned UNKNOWN; none returned INFEASIBLE.
- Direct PIE-to-PDE test at 0.9 (M and A): every term set certified – faces0, facesS and (on M) box at dws 0, box and product2 at dws 1. At this single fraction the runs do not separate the term sets. ndv relative to facesS at dws 0: faces0 3.6 times, box (dws 1) 5.9 times, product2 (dws 1) 6.6 – 6.7 times.
- 1-D PDE 7: stock and Rlight0 give identical programs (ndv 568, m 75, dw 1) and verdicts at lam 4.0 and 4.5. Light and heavy both certified every bisection point up to 4.65352; no value above 4.65352 was tested, so the certified threshold is not bracketed from above. 4.65352 lies below the FD estimate 4.65378 (N = 200) and below the library value 4.66.
- PDE 23: the container executive gives gam 23.7057 (light) and 23.7058 (heavy). The light value differs from the library header light value 25.249 by -6.1 %; the heavy value agrees with the header 23.706 to 5 digits. The stock lpiscript path was not run.

Inferred.
- On these single-domain variable-coefficient systems the separable patterns of Rlight0 do not change the certified threshold at 8-step bisection resolution (1/256 of lam* for A; the bracket ceiling for B, C). The difference appears only on the mixed domain M, where the reach moves from 0.736 to 0.977 at 0.4 % more decision variables.
- The FD references converge from below (N = 40 < N = 60); against the Richardson estimates the reach fractions are 0.02 – 0.03 % lower than stated.

Not done: no bisection of the direct-test term sets (each run at f or lam fraction 0.9 only); no Rlight0 run at the heavy tier; stock lpiscript path for PDE 23.

## Files
Scripts (inhom\): inh_fd.m, inh_fdmax.m, inh_sys.m, inh_2d_sys.m, inh_run_1A.m, inh_run_1B.m, inh_run_1C.m, inh_run_2.m, inh_run_3.m, inh_run_4.m, inh_run_5.m, pf_qdirect_pie.m, chain.sh.
Logs (inhom\logs\): inh_fd.log, inh_run_1A.log, inh_run_1B.log, inh_run_1C.log, inh_run_2.log, inh_run_3.log, inh_run_4.log, inh_run_5.log, chain.log.
Records (inhom\results\): thresholds.mat, item1_A.mat, item1_B.mat, item1_C.mat, item2_exec.mat, item3.mat, item5.mat, qdirect_*.mat.
