# PI operator inversion in PIETOOLS: dependency map and measured state

MMP (with Claude), 2026-10-09. Scope: the 1-D and 2-D operator inverse routines, the
synthesis reconstruction that depends on them (`getController`, `getObserver`,
`getObserver_2D`), and what the container path (`copvar` / `cdopvar`) currently does.
Numbers in Sec. 4 were measured on this workstation (MATLAB R2025b, MOSEK) with the probe
scripts named there; everything else is read from the code at commit 96ef81f6 plus the
uncommitted `executives_sopvar/` tree.

## 1. Routines

| file | author, date | status | method |
|---|---|---|---|
| `opvar/@opvar/inv.m` | PIETOOLS team 2024, DJ 03/25/2025 | live | wrapper: `inv_opvar(P[,tol])` |
| `opvar/inv_opvar.m` | SS 05/31/2022, last edit 08/2023 | live | dispatcher, Sec. 2 |
| `opvar/inv_opvar_old.m` | MP 07/01/2020, DJ to 04/07/2025 | live (called by `inv_opvar` when `isequal(R1,R2)`) | closed form for separable kernels: Theorem 9 of arXiv 1806.08071 (Peet, multi-delay H-infinity synthesis), generalized to non-self-adjoint by DJ 04/2025 through a linear solve `AMat/BMat` |
| `opvar/inv_opvar_new.m` | SS 07/2023 | dead: no caller | discrete variant of `inv_opvar` |
| `opvar/inv_opvar_2.m` | SS 02/24/2026, last 03/30/2026 | dead: no caller (only its own recursion) | Lemma 16 / Cor. 17 of arXiv 2208.13104 (Gohberg-Krein) by RK4 on a grid, Chebyshev least-squares fit back to polynomials, 4-block Schur complement |
| `sopvar/@sopvar/inv_1D.m` | SS 02/11/2026 | dead for the class: static method on `quadPoly` inputs, which the 08/2026 `sopvar` no longer stores; one caller, `ndopvar/@quadPoly/test/transfer_function_test.m` | same Gohberg-Krein, RK4, with `decomposeLR` (SVD) for the kernel split and a two-sided monomial least-squares fit on the triangle |
| `opvar/2D/@opvar2d/inv.m` | DJ 07/01/2024 | live | separable: `inv_opvar2d_separable`; otherwise `mrdivide(1,P1op,...)` with `P1op` UNDEFINED (input is `Pop`): see Sec. 4 |
| `opvar/2D/@opvar2d/inv_opvar2d_separable.m` | DJ 01/28/2024, 01/16/2025 | live | exact inverse of a separable 2-D operator |
| `opvar/2D/@opvar2d/mrdivide.m`, `mldivide.m` | DJ 07/11/2022, 08/03/2022, 01/14/2025 | live | least-squares `X*A = B` on monomial representations with degree escalation (`deg_fctr`, `deg_fctr_max`); `mldivide` through the adjoint |
| `opvar/2D/@opvar2d/is_separable.m` | DJ 2022 | live | the fifth check is the bare statement `is_full_int(5);` with no assignment: see Sec. 4 |

Container classes: `@copvar`, `@cdopvar`, `@sopvar` have no `inv`, `mrdivide` or
`mldivide` (listing of the class folders, 2026-10-09).

## 2. What `inv_opvar` does (the live 1-D routine)

```
inv_opvar(P, tol = 1e-8)
  P must be square (dim(:,1) == dim(:,2))
  isequal(P.R.R1, P.R.R2)             -> inv_opvar_old(P)          [separable: closed form]
  dim(2,:) == 0                       -> inv(P.P)                   [matrix]
  dim(1,:) == 0                       -> R-only branch:
        getsemisepmonomials(P)        R0^{-1} by polyfit of degree 4 at 101 points (orderapp = 4),
                                      Ra = R0^{-1} R1, Rb = R0^{-1} R2 split into F_i(s) G_i(theta)
                                      by grouping monomials on the lower-degree variable
        A = [G1 F1, G1 F2; -G2 F1, -G2 F2](var1)
        U = sum_k U_k, U_k = int_a^s A U_{k-1}, stopped when max|coef U_k| < tol or k = Nmax
            (Nmax from ||A||^N/N! > eps, capped at 19), coefficients < tol zeroed
        U22 = U22(b), Pm = [0 0; U22\U21 I], Uinv = series of int Uk A
        Pinv.R = {R0^{-1}, C U (I-Pm) Uinv B R0^{-1}(theta), -C U Pm Uinv B R0^{-1}(theta)}
  otherwise (4-block)                 -> Ainv = inv(P.P); TB = inv_opvar(R - Q2 Ainv Q1, tol)
                                         Pinv = [Ainv + Ainv Q1 TB Q2 Ainv, -Ainv Q1 TB; -TB Q2 Ainv, TB]
```

The 4-block branch is the Schur complement through the matrix block; Lemma 18 of
arXiv 2208.13104 uses the other one (invert the L2 block, then the matrix
`T = P - Q1 Rhat Q2`). Both are valid for a coercive operator.

`tol` has two roles: the series termination criterion and the coefficient truncation
threshold. `getController` passes `tol = 1e-4` by default (and then `clean_opvar(K,tol)`);
`getObserver` passes `1e-5`; the wrapper default is `1e-8`.

## 3. Callers and consumers

```
PIETOOLS_Hinf_control(PIE,st)   -> P, Z (lpigetsol) -> getController(P,Z)        K = Z P^{-1}
PIETOOLS_H2_control             -> getController
PIETOOLS_Hinf_estimator         -> getObserver(P,Z)                                 L = P^{-1} Z
PIETOOLS_H2_estimator           -> getObserver
PIETOOLS_Hinf_estimator_2D      -> getObserver_2D(P,Z): separable -> inv; else mldivide(P,Z,2*ones(4,2),tol)
lpiscript('hinf-controller' | 'h2-controller' | 'hinf-observer' | 'h2-observer') -> the above
DEMO5 (getObserver), DEMO6 (getController(P,Z,1e-3)), manual snippets Ch. 7.6
executives_sopvar/private/synth_build_sop.m  ('ctrl'/'est')  -> copvar2opvar(Psol), copvar2opvar(Zsol)
                                                                -> stock getController / getObserver / getObserver_2D
executives_sopvar/private/h2_build_sop.m     ('ctrl'/'est')  -> same
```

Downstream of the gain (all require `opvar` / `opvar2d`, none accept a container):
`closedLoopPIE(PIE,K[,'observer'])` (block assembly by hand on the parameters),
`piess(T, A+B2*K, ...)` in DEMO6, then `PIESIM`. Nothing in `PIESIM`, `converters` or
`lpi_programming` inverts an operator; the `inv(` calls there are on matrices.

Dependencies of the 1-D routines inside the toolbox: `opvar` arithmetic (`mtimes`, `plus`,
`ctranspose`, `subsref`), `polynomial` (`int`, `subs`, `poly2basis`, `polyfit` on sampled
values, `combine`), `clean_opvar`, `isvalid`, `mat2opvar`. Outside: `polyfit`, `integral`
(`inv_opvar_old`), `lu`, `svd` / `svds` (`decomposeLR`), `pagemtimes` (`inv_opvar_2`).

Container side: `copvar2opvar` (1-D grid to `opvar`; rectangular grids since 10/08/2026),
`copvar2opvar2d`, `opvar2copvar` (in `sopvar/Testfolder/converters/`), `sopvar2opvar`.
A 1-D L2 -> L2 block stores `params{1}` = R0, `params{2}` = R1 (lower), `params{3}` = R2
(upper), each as the coefficient matrix of `(I_m kron ZL(s)^T) C (I_n kron ZR(theta))` on
shared per-block monomial bases `ZL`, `ZR`. So the kernel split `R_i = F_i(s) G_i(theta)`
that Lemma 16 needs is a factorization of `C` (`decomposeLR` does an SVD of it), with no
monomial bookkeeping.

## 4. Measured state (probe_inv_1d.m, probe_inv_2d.m in the session scratchpad)

Plant: reaction-diffusion with the profile `s(1-s)` on `w`, `u` in-domain, the plant of
`executives_sopvar/tests/test_executives_sop.m`; `PIETOOLS_Hinf_control` at `lpisettings('light')`,
MOSEK, gamma 0.056851286 in 3.7 s. Solved `P` has `dim = [0 0; 1 1]`, `R0` degree 2,
`R1`, `R2` degree 5, `R1 ~= R2`, so the R-only Gohberg branch runs (no ODE block).

Residual = max over a 21-point grid (triangles for the kernels) of the parameters of
`P*Pinv - I` and `Pinv*P - I`.

| inverse of the solved P | time | max residual | deg of Pinv.R1 (terms) |
|---|---|---|---|
| `inv_opvar`, tol 1e-4 (getController default) | 0.1 s | 6.5e-7 | 21 (142) |
| `inv_opvar`, tol 1e-6 | 0.1 s | (between) | 23 (168) |
| `inv_opvar`, tol 1e-8 (wrapper default) | 0.1 s | 1.3e-8 | 25 (188) |
| `inv_opvar`, tol 1e-10 | 0.2 s | 1.3e-8 | 29 (240) |
| `inv_opvar_2` (RK4, defaults) | 0.2 s | 4.5e-12 | 12 (49) |

Controller: `K.Q1` from the executive (tol 1e-4 plus `clean_opvar` 1e-4) differs from
`Z*inv(P,1e-8)` by 1.3e-5 relative; closed-loop numerical gains (`closedLoopPIE` then
`pie_witness_sop` 'gain', N 24) are 0.032872613 / 0.032872612 / 0.032872612 for the three
gains against the certified 0.056851286. On this plant the inverse is not the limiting
factor.

Hand-built 3-PI operators on [0,1], tol 1e-8:

| case | `inv_opvar` residual | `inv_opvar_2` residual |
|---|---|---|
| c1: R0 = 1+0.5s, R1 = 0.3 s th, R2 = 0.3 (s+th) | 1.4e-2 | 1.6e-6 |
| c2: R0 = 1, R1 = 0.3 s, R2 = 0.3 th | 3.0e-1 (wrong operator) | 4.6e-10 |
| c3: R0 = 1+0.5s, R1 = 0.3 s, R2 = 0.3 th | 1.0e-3 | 4.6e-7 |
| c4: 2x2 non-separable, R0 = [2+s, 0.2s; 0.2s, 1.5] | 1.6e-3 | 5.8e-7 |

`inv_opvar` returns kernels of degree 38 to 64 with 400 to 1100 terms on these; `inv_opvar_2`
degree 12 with 49 terms (its fit degrees are options: mulDeg 6, kerDeg [6 6], N 100).

2-D: an `opvar2d` with `Rxx{2} ~= Rxx{3}` is reported non-separable and `inv` fails with
"Unrecognized function or variable 'P1op'". A separable one inverts to 1.7e-18.

Causes of the c1-c4 residuals and the c2 wrong operator: Sec. 5.

## 5. Verification pass (four independent checks, each reproduced in MATLAB)

**5.1 `inv_opvar` drops the upper kernel when the variable orders differ: CONFIRMED.**
`getsemisepmonomials`, line 250 of `inv_opvar.m`, reads
`[val,~,~] = unique(full(Rb.degmat(:,var1_loc_Ra)));` with `Ra`'s column position. When
`R0` is constant, `Rinv*R2` keeps only the variables of `R2` (`combine` drops the absent
one), the padding branch appends the missing variable AFTER the present one, and `Ra`,
`Rb` end up with different variable orders. The line then selects the theta column of
`Rb`, the `G{2}` loop finds no row, and `F{2}*G{2} = 0`: the whole `R2` part is removed
from the inverse (instrumented: `F{2} = s`, `G{2} = 0`; `Pinv.R.R2` has one term). The
stock residual equals the dropped term: 0.30 on c2, c6 (R2 = 0.3 th), c10; 0.20 on c8
(R2 = 0.2 s + 0.1 th^2, the s term dropped); c9 (R2 = 0.2 + 0.1 th) is right by accident
because the degree-0 group collects every row. With a non-constant `R0` both products
carry both variables and the line is harmless (c1, c3, c4 unchanged). One-token fix
(`var1_loc_Rb`), tested in a scratch copy: c2 0.30 -> 2.1e-3, c6 0.30 -> 7.6e-4, c8 0.18
-> 9.5e-4, c10 0.30 -> 9.4e-4. The 1e-3 floor remains (5.2).

**5.2 The 1e-3 floor of `inv_opvar` is a SIGN ERROR in the series for `U^{-1}`, not the
degree-4 fit: hypothesis REFUTED, cause established.** `U^{-1} = V` solves `V' = -V A`,
`V(a) = I`, whose Picard series alternates: `V = I - int A + int(int A)A - ...`. Lines
118-131 of `inv_opvar.m` compute `Ukinv = int(Ukinv*A)` with no sign and subtract every
term (`Uinv = Uinv - Ukinv`), giving `I - T1 - T2 - T3 - ...`. The error `2(T2 + T4 + ...)`
is second order in the kernel: `U*Uinv - I` on the grid is 5.0e-2 (c1), 1.9e-2 (c3),
1.0e-2 (c4), 8.3e-2 (c5: R0 = 1, exact fit, residual 1.0e-2 / 3.5e-2 anyway), 1.7e-4 on
the solved Lyapunov P, where `R1` and `R2` nearly coincide and the second term is 8.4e-5.
Raising `orderapp` 4 -> 8 -> 12 moves the `R0` residual 4.3e-5 -> 4.5e-9 -> 3.7e-13 and
leaves `R1`, `R2` unchanged to three digits; the series breaks at 4-8 terms with the
last term below `tol` (the cap 19 binds only from `normA` 1.3, and at `normA` 2.6 the loop
never breaks and the output is garbage, 4.2 / 12, without a warning). Two-line fix
(`Ukinv = -int(Ukinv*A,...)`, `Uinv = Uinv + Ukinv`), tested in a scratch copy at
`tol` 1e-8, `orderapp` 4 (R1 / R2 residual): c1 3.8e-3 / 1.4e-2 -> 7.1e-6 / 1.6e-5;
c3 -> 1.2e-5 / 6.7e-6; c4 -> 5.2e-6 / 6.9e-6; c5 -> 1.1e-9 / 1.1e-9; Lyapunov P
1.2e-8 -> 1.8e-9 (6e-14 at `tol` 1e-12). After the fix the residual is the maximum of
the `R0^{-1}` fit error (4.3e-5 at degree 4 for `R0 = 1 + 0.5 s`, never checked or
reported) and the `tol` truncation; `orderapp` 8 takes c1 to 2.6e-9 / 4.2e-9. Cost of
the stock representation: the kernel degree of `Pinv` is `deg C + deg U + deg Uinv +
deg B + orderapp` (64 on c1, 137 with 4760 terms at kernel scale 2), all of which flows
into `Z*Pinv` in `getController`.

**5.3 `inv_opvar_2` against Lemma 16 / Cor. 17 / Lemma 18: default path CORRECT, every
non-default option defective.** The mapping `C = -R0^{-1}[F1 F2]`, `B = [G1; -G2]`, RK4
for `U` and `V`, `U22` the trailing block, `Pm`, `M1`, `M2` with `R0(theta)^{-1}` on the
right, the page indices (i = s, j = theta) and the coefficient layouts (`polynomial`
column-major, `poly_terms_2D` reshape, `fitpoly_2D_cheb` ordering, the affine map, the
Chebyshev-to-monomial conversion) are all consistent; measured 7.4e-8 on an asymmetric
2x2 whose transposed reshape would give 0.4, 1e-6 on [-1,2], 1e-6 on three 4-block
shapes where `inv_opvar` gives 1e-2..1e-3. Defects, all reproduced: (1) line 292
`info = [infoFinite, infoInfinite]` concatenates structs with different fields, so
`[Pinv,info] = inv_opvar_2(P)` ERRORS on every 4-block operator (the Lyapunov shape);
(2) `fitKernelEverywhere = 0` fills one triangle and fits the whole square: residual
0.25 on c1; (3) `outIsPoly = 0` stores the sample arrays in the opvar, which no product
accepts; (4) the matrix block is inverted first with no check: `P.P = 0`, `Q1 = Q2 = 1`,
`R0 = 1` (invertible operator) returns NaN silently, Lemma 18 gives 2.2e-16; (5) the
matrix-first pivot inflates the kernels by `||Q||^2 / sigma_min(P.P)`: at `P.P = 1e-2`
residual 1.2e-4 against 1.8e-6 for the Lemma 18 pivot (equal at `P.P = 2`); (6)
`rcond(U22) < 1e-12` and a singular `R0(t_i)` only warn; (7) `kerDeg` must be 1x2,
undocumented, scalar errors; (8) a zero kernel autofills to a 49-term zero polynomial
and adds m states to the ODE; (10) the accuracy floor is the kernel fit degree, not `N`:
c1 R1 residual 1.1e-6 at `kerDeg [6 6]`, 3.0e-10 at `[10 10]`, unchanged by `N = 400`.
The Lemma 18 route assembled by opvar algebra returns the L2 block of the inverse at
degree 28 (two fitted degree-12 kernels composed) where the matrix-first route fits it
once at degree 12; the container inverse of Sec. 6 has the same property.

**5.4 `is_separable` reports a split `Ry2` or `R2y` as separable: CONFIRMED.** Line 79 of
`opvar/2D/@opvar2d/is_separable.m` is the bare `is_full_int(5);`. An `opvar2d` on
`L2[y] x L2[x,y]` with `Ry2{2} = 0.1`, `Ry2{3} = 0.2` (everything else separable) returns
`is_sep = 1`, `inv` takes `inv_opvar2d_separable`, which reads only the index-2 kernels,
and returns `inv` of the operator with the upper kernel overwritten by the lower:
`P*Pinv - I` = 0.10 on `Ry2{3}`, no warning; the same with `R2y`. The control case (all
equal) inverts to 1.4e-17. `getObserver_2D` gates on `is_separable`, so a solved 2-D
Lyapunov operator with such a split would give a wrong observer gain silently; whether
`poslpivar_2d` can produce one was not checked. With the one-line fix the operator falls
through to the non-separable branch of `inv.m`, which errors on `P1op` (Sec. 1), a
second one-token fix. Both files are frozen 2-D code: reported, not edited.

## 6. The container inverse (new, 10/09/2026)

Files: `sopvar/@sopvar/inv.m`, `sopvar/@copvar/inv.m`,
`lpi_programming_sopvar/getController_sop.m`, `getObserver_sop.m`; wired into
`executives_sopvar/private/synth_build_sop.m` and `h2_build_sop.m` (1-D; 2-D keeps
`getObserver_2D`). Test: `sopvar/Testfolder/Test_copvar_inv.m`.

Method. The block stores `R_i(s,t) = (I kron ZL(s)') C_i (I kron ZR(t))`, so the kernel
split Lemma 16 needs is one SVD of each `C_i` (`R_i = F_i G_i`, ranks `r_1`, `r_2`); `U`,
`V` by RK4 on 101 nodes; `Pm` from `U(b)`; the three parameters of the inverse evaluated
on the node grid (kernels on the whole square from the formula, which is their smooth
continuation) and fitted by Chebyshev least squares at one degree `d` raised 4, 6, ...,
16 until the relative RMS residual of all three is below 1e-8; converted to monomials on
`ZL = ZR = (0:d)'` with the multiplier in the ZR-degree-0 columns (canonical form). The
container inverse takes `Rh = R^{-1}` once and the finite-dimensional Schur complement
`T = Pm - Q1 Rh Q2` (Lemma 18), all products by the block composition. `R0^{-1}` is
rational, so the inverse is a polynomial approximation; `info.relrms`, `info.d` say how
good.

Measured (`Test_copvar_inv(false)`, 17 checks, 10/09/2026):

| case | time | d | ranks | fit relrms | max residual `P*Pinv - I` (container; opvar composition) | vs `inv_opvar_2` on the grid |
|---|---|---|---|---|---|---|
| c1 | 0.23 s | 10 | [1 2] | 1.2e-9 | 4.3e-10 ; 4.3e-10 | 1.6e-6 |
| c2 | 0.05 s | 6 | [1 1] | 8.5e-10 | 4.6e-10 ; 4.6e-10 | 9.2e-13 |
| c3 | 0.05 s | 8 | [1 1] | 3.9e-9 | 4.5e-9 ; 4.5e-9 | 3.7e-7 |
| c4 (2x2) | 0.05 s | 8 | [2 4] | 8.6e-9 | 6.3e-9 ; 6.3e-9 | 2.3e-7 |
| c5 separable | 0.06 s | 8 | [1 1] | 8.0e-9 | 4.5e-9 | 3.7e-7 |
| c6 lower only | 0.05 s | 10 | [1 0] | 3.0e-10 | 8.8e-11 | 4.0e-7 |
| c7 4-PI, R^2 x L2^2 | 0.30 s | 8 | [2 4] | 8.6e-9 | 6.3e-9 | (not applicable) |
| identity on R x L2^2 | | | | | 1.3e-15 | |

Executives on the test plant (light, MOSEK), gains judged by their defining equations:

| executive | gamma | container gain | stock gain on the same P, Z |
|---|---|---|---|
| `Hinf_control_sop` | 0.03286578 | `|K P - Z|/|Z|` = 2.4e-11; closed-loop numerical gain 0.03286346 | `getController`: 8.1e-5 (tol 1e-4 truncation); K differs 8.0e-5 |
| `Hinf_estimator_sop` | 1.39e-4 | `|P L - Z|/|Z|` = 7.0e-9 | `getObserver`: 7.4e-3; L differs 4.7e-3 |

`test_executives_sop(true)` passes its 19 checks with the new wrappers (H2_control_sop
0.05222238, H2_estimator_sop 1.63e-4 among them).

Not done: N-D (no kernel inverse exists beyond the separable case), the gains without
the rational `R0^{-1}` (Cor. 11 of arXiv 1806.08071: evaluate `K` on a grid, fit once),
and any edit of the stock routines (SS, DJ files).

## 7. Fixes to the stock routines, APPLIED 10/09/2026 (maintainer's instruction), each re-measured

| file, line | change | measured, before -> after |
|---|---|---|
| `opvar/inv_opvar.m` 121, 128 | `Ukinv = -int(Ukinv*A,...)`; `Uinv = Uinv + Ukinv` | with the two lines below: c1 1.4e-2 -> 4.5e-9, c2 0.30 -> 1.4e-9, c3 1.0e-3 -> 4.5e-9, c4 1.6e-3 -> 6.3e-9 (tol 1e-8; kernel degrees 72, 20, 50, 59) |
| `opvar/inv_opvar.m` 250 | `var1_loc_Ra` -> `var1_loc_Rb` | (in the row above) |
| `opvar/inv_opvar.m` 162 | `orderapp` 4 -> 8 | (in the row above; `R0` residual 4e-5 -> 5e-9) |
| `executives/utility_functions/getController.m` | `inv(P)` at the wrapper default; `tol` only truncates K | `|K P - Z|/|Z|` 8.1e-5 -> 8.6e-7 on the test plant (`getObserver`, unchanged code, 7.4e-3 -> 9.6e-5 through the inverse fixes) |
| `opvar/2D/@opvar2d/inv.m` 66 | `P1op` -> `Pop` | a non-separable `Rxx` operator inverts through `mrdivide`: residual 9e-16 (constant kernels, 1.3 s) and 5e-16 (`0.1 s1`, `0.2 s1_dum`, 0.6 s) |
| `opvar/2D/@opvar2d/is_separable.m` 79 | `is_full_int(5) = false;` | the split-`Ry2` operator of Sec. 5.4 returns `is_sep = 0`, `is_full_int(5) = false` |
| `opvar/inv_opvar_2.m` | `info = struct('R',..,'L2',..)`; fit regions follow `fitKernelEverywhere`; `outIsPoly = 0` returns the samples in `info.samples`; singular `P.P` errors; help text | `[Pinv,info]` on the 4-block operator returns (residual 5.8e-7); `fitKernelEverywhere = 0` residual 0.25 -> 8.0e-7; `outIsPoly = 0` returns a polynomial opvar plus samples; `P.P = 0` errors instead of NaN |

Not changed: the fit error of `R0^{-1}` in `inv_opvar` is still not checked; the stock
`getObserver` keeps `inv(P,1e-5)`.
