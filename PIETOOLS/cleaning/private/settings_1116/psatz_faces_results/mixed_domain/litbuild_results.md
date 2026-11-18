# Literature candidates built in PIETOOLS (litbuild, 2026-11-17)

Catalog: `PIETOOLS/sopvar/psatz_faces_results/mixed_domain/COUPLED_PDE_CANDIDATES.md` (17 items).
Everything here lives in the scratch folder `litbuild\`: systems `lb_sys.m`, references `lb_s0.m`
(`logs\lb_s0.log`, `results\ref.mat`), one script per run group `lb_i<item>*.m`, logs `logs\<script>.log`
(+ `.log.meta`: start, end, wall, kills), result structs `results\ex_*.mat` and `results\qdirect_*.mat`.
No repository file was edited. MATLAB R2025b, MOSEK, one MATLAB process at a time, 01:17 – 05:13 on 2026-11-17
(about 3 h 50 min of MATLAB wall time, serial).

## Conventions

- **stabP**: `PIETOOLS_stability_sop(PIE,st,'P')` through `pf_exec_row`: P >= eppos I,
  T'PA + A'PT + 2 epneg T'PT <= 0. epneg = 0 unless a row says otherwise. With epneg = 0 a certificate
  shows that a P-weighted norm of Tx does not increase (no rate); with epneg > 0 it certifies the decay
  rate epneg of ||Tx||. **hinfP**: `PIETOOLS_Hinf_gain_sop(PIE,st,'P')`.
- **cert** = 1: the executive flag (optimal or inaccurate exit through the eta gate). **dw**: slack
  weight raises the executive took (its own loop). **ndv**: decision variables; **m**: SDP rows; **t**:
  seconds wall for the call (build + solve).
- **storage** 'stock' = repository settings of the tier; 'Rlight0' = the proposed light row of
  `pf_rows` (1-D identical to stock; 2-D stock light lengths, separable patterns; 3-D the N-D default).
  The 'Rlight0' row fixes the storage lengths, so 'Rlight0' at 'heavy' changes only the other tier
  settings (item 11 shows the same program at light and heavy).
- **direct**: `lb_qdirect.m`, a copy of `pf_qdirect_pie.m` with the output folder changed: free Q with
  the shape of T (`lpivar_cdopvar`, degrees mult 1, int [1 1] (left, right), out 1, in 1),
  T'Q = Q'T, T'Q >= 1e-4 T'T, -(A'Q + Q'A) >= 0, light tier, dw 0, MOSEK. This is the PIE-to-PDE
  notion with no decay rate. Verdict "in" = optimal with eta < 1e-8, "out" = infeasible, "U" = unknown.
  Term sets from `pf_terms`: faces0 (2N faces, offset 0), facesS (faces, scalar offset rule).
- **Size cap**: `litbuild\shadow\lpisolve.m`, a scratch shadow with the repository body plus a check,
  refuses a program above 3,000,000 decision variables and records its size. **Watchdog**: `lb_run.sh`
  kills a script after 1,500 s.
- Bisection (`lb_bis.m`): fraction f of the reference threshold, certified at f -> lower end = f.
  "reach" = the largest certified fraction.

## Item 8 – delayed heat equation, 1-D control (S18) and 2-D extension

**1-D (control, = suite S18 with sgn -1).** x1(s1) on [0, pi], x2(s1,s2) on [0, pi] x [0, tau]:

    x1 = pde_var(s1,[0,pi]);  x2 = pde_var([s1;s2],[0,pi;0,tau],[2;1]);
    x1_t = x1_s1s1 + 1.9 x1 - x2(s1,tau);   x2_t = -x2_s2;
    x1 = x2 = 0 at s1 = 0, pi;   x2(s1,0) = x1(s1)

**2-D.** x1(s1,s2) on (0, pi)^2, x2(s1,s2,s3) on (0, pi)^2 x [0, tau]:

    x1 = pde_var([s1;s2],[0,pi;0,pi]);  x2 = pde_var([s1;s2;s3],[0,pi;0,pi;0,tau],[2;2;1]);
    x1_t = x1_s1s1 + x1_s2s2 + 2.9 x1 - x2(s1,s2,tau);   x2_t = -x2_s3;
    x1 = x2 = 0 on s1 = 0, pi and s2 = 0, pi;   x2(s1,s2,0) = x1(s1,s2)

The converter accepts the 2-D state as the datum of the 3-D state on the face s3 = 0 (3.2 s; PIE in
the registry s1, s2, s3; spaces L2[s1,s2] and L2[s1,s2,s3]).

**Reference [derived, exact for the modal reduction].** Each Dirichlet sine mode gives
y' = (r - k) y - y(t - tau), k = n^2 (1-D) or m^2 + n^2 (2-D). For a = r - k in (-1, 1) the delay
margin is arccos(a)/sqrt(1 - a^2); a <= -1 is stable for every delay. Minimum over the modes
(`lb_s0`): tau* = 1.034726 in 1-D (r = 1.9) and in 2-D (r = 2.9); the minimizing mode has a = 0.9 in
both.

| label | test | storage / terms | tier | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|---|
| S18 0.5 tau* | stabP | stock | light | infeasible | 0 | 1 | 508,228 | 4,757 | 42 |
| S18 0.5 tau* | stabP | Rlight0 | light | infeasible | 0 | 1 | 510,057 | 4,841 | 47 |
| S18 0.5 tau* | direct | faces0 | – | optimal, eta 2.2e-11 (in) | – | 0 | 227,510 | – | 10 |
| 2-D 0.5 tau* | direct | faces0 | – | not solved: 47,433,948 variables > cap (build 350 s, about 38 GB resident) | – | 0 | 47,433,948 | – | 350 |
| 2-D 0.5 tau* | stabP | stock | light | infeasible | 0 | 0 | 1,668,608 | 22,476 | 476 |
| 2-D 0.5 tau* | stabP | Rlight0 | light | infeasible | 0 | 0 | 1,680,004 | 29,107 | 653 |

The 1-D control reproduces the suite verdicts: stock light infeasible after one raise (508,228
variables, the suite size), direct test with faces0 certified (227,510 variables, the suite size).
In 2-D both storages return infeasible at dw 0 and the executive took no raise. Not run, and why:
the 0.9 tau* runs of the 2-D system (stabP under both storages, the direct test). The 0.5 tau*
runs were negative at about 8 – 11 min each; the direct-test program has the same structure at
0.9 tau* (only coefficients and the s3 interval change), hence the same 47.4M size (inferred, not
built).

## Item 4 – linear Turing pair, two 2-D states coupled in-domain

    x = pde_var(2,[s1;s2],[0,1;0,1]);           % x = (u, v)
    x_t = diag(1,d) (x_s1s1 + x_s2s2) + J x,   J = [20 -25; 30 -30]
    x_s1 = 0 at s1 = 0, 1;   x_s2 = 0 at s2 = 0, 1   (Neumann, all edges)

PIE: one vector-valued 2-D state; the converter returns T as a 4 x 4 copvar for the Neumann version
(1 x 1 for a Dirichlet variant built as a check, which converts as well).

**Reference [derived].** Cosine modes, k^2 = pi^2 (m^2 + n^2), m, n >= 0. Mode k is stable iff
tr_k = -10 - (1 + d) k^2 < 0 and det_k = (20 - k^2)(-30 - d k^2) + 750 > 0. tr_k < 0 for every k.
det_k is linear in d with slope -(20 - k^2) k^2, negative only for k^2 < 20, i.e. k^2 = pi^2 and
2 pi^2. det_k = 0 at d = (750 - 30 (20 - k^2)) / (k^2 (20 - k^2)): 4.461640 at k^2 = pi^2 (modes (1,0),
(0,1)) and 144.0 at k^2 = 2 pi^2. Hence d* = 4.461640: stable for d < d*. Check (`lb_s0`, modes
m, n <= 20): max Re eig = -0.759 at 0.9 d*, -0.070 at 0.99 d*, +0.069 at 1.01 d*. The Dirichlet
variant with 2J has the same d*.

| label | test | storage | tier | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|---|
| 0.9 d* (01:49) | stabP | stock | light | optimal | 1 | 0 | 2,102,408 | 19,546 | 740 |
| 0.9 d* (01:49, same script) | stabP | Rlight0 | light | not completed: script killed at the 1,500 s watchdog during the solve | – | – | 2,120,516 | – | > 760 |
| 0.9 d* (03:40) | stabP | Rlight0 | light | not built: harness assertion (see "Harness incident") | – | – | – | – | – |
| 0.9 d* (04:05) | stabP | Rlight0 | light | not completed: killed at the 1,500 s watchdog during the solve | – | – | 2,120,516 | – | > 1,100 |

ndv is the size printed by the run (`@@ LB_SIZE`, info.shape.ndv): 2,102,408 for stock.

Bisection of d/d* on [0.9, 1] (0.9 certified above), stock light, one solve per script invocation
(`lb_i4b_stock.m`, resumable):

| step | f = d/d* | d | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|
| 1 | 0.95 | 4.238558 | optimal | 1 | 0 | 2,102,408 | 19,542 | 909 |
| 2 | 0.975 | 4.350099 | optimal | 1 | 0 | 2,102,408 | 19,542 | 905 |
| 3 | 0.9875 | 4.405870 | optimal | 1 | 0 | 2,102,408 | 19,546 | 626 |
| 4 | 0.99375 | – | not run (time given to the gap study) | | | | | |

Largest certified fraction 0.9875 (d = 4.405870, stock light). No uncertified fraction was
measured: every completed step certified. Rlight0 bisection: not run (each Rlight0 solve exceeded
the watchdog; then time given to the gap study).

## Item 11 – exact H-infinity gain, 2-D heat with a 1-D disturbance

    x = pde_var([s1;s2],[0,1;0,1]);  w = pde_var('in',1,s1,[0,1]);  z = pde_var('out',1,[s1;s2],[0,1;0,1]);
    x_t = x_s1s1 + x_s2s2 + w(s1);   z = x;   x = 0 on all four edges

**Reference.** A = Laplacian (Dirichlet) is self-adjoint and strictly negative and C = I, so
||(i omega - A)^{-1} H w||^2 = <H w, (omega^2 + A^2)^{-1} H w> <= <H w, A^{-2} H w>: the gain is
attained at omega = 0 and equals gamma* = ||A^{-1} H|| (the Lidstrom – Rantzer – Morris formula with
B = 0). Separation in sin(n pi s1) gives gamma* = sup_n ||psi_n||_{L2(0,1)}, psi_n'' - (n pi)^2 psi_n = 1,
psi_n(0) = psi_n(1) = 0, attained at n = 1, psi_1 = -(1/pi^2) [1 - cosh(pi (s - 1/2)) / cosh(pi/2)].
Closed form (evaluated in `lb_s0`): gamma* = 0.04571818. Finite differences (5-point Laplacian,
N interior nodes per direction, sigma_max(A_N^{-1} H_N) sqrt(h)): 0.04571809 (N = 40),
0.04571814 (N = 60).

| label | test | storage | tier | status | cert | dw | gain | ndv | m | t | gain / gamma* - 1 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| LRM | hinfP | stock | light | optimal | 1 | 0 | 0.0459440156 | 209,445 | 4,165 | 17 | +0.494 % |
| LRM | hinfP | Rlight0 | light | optimal | 1 | 0 | 0.0459459178 | 210,006 | 4,178 | 12 | +0.498 % |
| LRM | hinfP | stock | heavy | optimal | 1 | 0 | 0.0459431727 | 837,690 | 7,416 | 111 | +0.492 % |
| LRM | hinfP | Rlight0 | heavy | optimal | 1 | 0 | 0.0459459178 | 210,006 | 4,178 | 12 | +0.498 % |

Every certified bound lies above gamma* (no bound below the true gain). Heavy lowers the stock bound
by 8.4e-7 (0.002 % of gamma*) at 4.0 times the variables.

## Item 6 – wave with partially delayed boundary damping (Xu, Yung, Li 2006)

Riemann form, e = w_s, v = w_t on (0, 1), delayed velocity y on s2 in [0, 1]; k = 1, tau = 2:

    e = pde_var(s1,[0,1],1);  v = pde_var(s1,[0,1],1);  y = pde_var(s2,[0,1],1);
    e_t = v_s1;  v_t = e_s1;  y_t = -(1/tau) y_s2
    v(0) = 0;   e(1) = -k mu v(1) - k (1 - mu) y(1);   y(0) = v(1)

**Reference [derived, exact].** w = sinh(lambda s) gives cosh(lambda) + k sinh(lambda)
(mu + (1 - mu) e^{-lambda tau}) = 0. With tau = 2 and X = e^{2 lambda} this is the quadratic
(1 + k mu) X^2 + (1 - k mu + k (1 - mu)) X - k (1 - mu) = 0, and the eigenvalues are
lambda = (1/2) log X_j + i pi n, n in Z. At k = 1, X = -1 is a root iff mu = 1/2, so the spectral
abscissa crosses 0 at mu* = 1/2 (the threshold of the paper). Abscissa (`lb_s0`; residual of the
characteristic function <= 4.4e-16): 0.0957 (mu 0.4), 0.0488 (0.45), 0 (0.5), -0.0514 (0.55),
-0.1060 (0.6), -0.3017 (0.75). At mu = 1/2 the eigenvalues lie on the imaginary axis: the system is
not exponentially stable there.

| label | test | storage / terms | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|
| mu 0.40 | stabP | stock | infeasible | 0 | 1 | 188,244 | 455 | 11 |
| mu 0.40 | stabP | Rlight0 | infeasible | 0 | 1 | 191,520 | 455 | 6 |
| mu 0.40 | direct | faces0 | infeasible (out) | – | 0 | 36,966 | – | 2 |
| mu 0.40 | direct | facesS | infeasible (out) | – | 0 | 10,506 | – | 1 |
| mu 0.45 | stabP | stock | infeasible | 0 | 1 | 188,244 | 455 | 6 |
| mu 0.45 | direct | faces0 | infeasible (out) | – | 0 | 36,966 | – | 3 |
| mu 0.49 | stabP | stock | infeasible | 0 | 1 | 188,244 | 455 | 5 |
| mu 0.49 | direct | faces0 | infeasible (out) | – | 0 | 36,966 | – | 1 |
| mu 0.50 | stabP | stock | optimal | 1 | 0 | 47,304 | 323 | 1 |
| mu 0.50 | direct | faces0 | optimal, eta 9.6e-10 (in) | – | 0 | 36,966 | – | 1 |
| mu 0.60 | stabP | stock | optimal | 1 | 0 | 47,304 | 323 | 1 |
| mu 0.60 | stabP | Rlight0 | optimal | 1 | 0 | 50,580 | 323 | 1 |
| mu 0.60 | direct | faces0 / facesS | optimal, eta 3.9e-9 / 3.9e-9 (in / in) | – | 0 | 36,966 / 10,506 | – | 1 / 1 |
| mu 0.75 | stabP | stock / Rlight0 | optimal / optimal | 1 / 1 | 0 / 0 | 47,304 / 50,580 | 323 | 1 / 1 |
| mu 0.75 | direct | faces0 / facesS | optimal, eta 1.1e-9 / 1.2e-9 (in / in) | – | 0 | 36,966 / 10,506 | – | 1 / 1 |

Decay-rate rows (stabP stock light, epneg > 0; `lb_i2b`): epneg certifies the rate of ||Tx||.

| mu | exact abscissa | epneg | status | cert | dw | ndv | t |
|---|---|---|---|---|---|---|---|
| 0.5 | 0 | 0.01 | infeasible | 0 | 1 | 188,244 | 6 |
| 0.6 | -0.105968 | 0.05 | optimal | 1 | 0 | 47,304 | 1 |
| 0.6 | -0.105968 | 0.10 | optimal | 1 | 0 | 47,304 | 1 |
| 0.6 | -0.105968 | 0.11 | infeasible | 0 | 1 | 188,244 | 5 |

The certified rate at mu = 0.6 lies in [0.10, 0.11), against the exact 0.105968. At mu = 0.5 the
epneg-0 certificate (above) holds and the 0.01 rate does not, as the spectrum on the imaginary axis
requires.

Bisection of mu on [0.5, 0.6] (certified side above; 8 steps; `lb_i6b`):

| test | certified mu | not certified mu | smallest certified mu |
|---|---|---|---|
| stabP stock | all 8: 0.55, 0.525, ..., 0.500391 (dw 0, 47,304 variables, 1 – 3 s) | none | 0.500391 (gap 3.9e-4 to mu*) |
| stabP Rlight0 | all 8 (dw 0, 50,580 variables) | none | 0.500391 |
| direct faces0 | 0.5875, 0.58125, 0.578125, 0.576562, 0.575781, 0.575391 (eta 3.2e-9 – 4.1e-9) | 0.55 (optimal, eta 1.2e-8), 0.575 (optimal, eta 1.1e-8) | 0.575391 |

The direct-test bisection was steered by the eta gate, not by infeasibility: all eight exits are
'optimal', and the two exits counted out have eta 1.2e-8 and 1.1e-8, just above the 1e-8 gate.
At 0.500391, 0.51 and 0.525 the direct test returns optimal with eta 2.9e-8, 2.5e-8 and 1.0e-8
(`lb_i6c`); at mu = 0.5 it returns eta 9.6e-10. The eta of the direct test therefore does not
separate this family at the 1e-8 gate; stabP is the clean discriminator here.

## Item 9 – damped wave with a delayed displacement (Fridman, Orlov 2009, Example 2)

z_tt = 0.1 z_xx - 2 z_t - z(x, t - tau) on (0, pi), Dirichlet. Two PIE forms were built:

- **FO2r** (Riemann form on (0, pi/2), z(0) = 0, z_s(pi/2) = 0: the odd sine modes of the Dirichlet
  problem, whose first mode is the critical one):

      e = pde_var(s1,[0,pi/2],1);  v = pde_var(s1,[0,pi/2],1);  y = pde_var([s1;s2],[0,pi/2;0,tau],[0;1]);
      e_t = v_s1;  v_t = 0.1 e_s1 - 2 v - int_0^s1 y(r,tau) dr;  y_t = -y_s2
      v(0) = 0;  e(pi/2) = 0;  y(s1,0) = e(s1)

- **FO2s** (the S18 form on (0, pi), Dirichlet):

      z = pde_var(s1,[0,pi]);  v = pde_var(s1,[0,pi],0);  y = pde_var([s1;s2],[0,pi;0,tau],[2;1]);
      z_t = v;  v_t = 0.1 z_s1s1 - 2 v - y(s1,tau);  y_t = -y_s2
      z = y = 0 at s1 = 0, pi;  y(s1,0) = z(s1)

**Reference [derived].** Mode k: y'' + 2 y' + 0.1 k^2 y + y(t - tau) = 0. A crossing at i omega
requires (c - omega^2)^2 + 4 omega^2 = 1, c = 0.1 k^2, which has a solution only for c < 1 (k = 1, 2, 3).
Per-mode margins: 2.881401 (k = 1), 3.497930 (k = 2), 8.758610 (k = 3); k >= 4 stable for every delay.
tau* = 2.881401 for both forms. Check: rightmost root of mode 1 by Newton from a grid: -0.1949 at
0.5 tau*, -0.0200 at 0.9 tau*, +0.0155 at 1.1 tau*. The paper reports 0.31 at decay rate 0.8 (a
different question); not reproduced here.

| label | form | test | storage / terms | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|---|
| 0.1 tau* | FO2r | stabP | stock | infeasible | 0 | 1 | 840,208 | 6,042 | 90 |
| 0.5 tau* | FO2r | stabP | stock | infeasible | 0 | 1 | 840,208 | 6,042 | 128 |
| 0.5 tau* | FO2r | stabP | Rlight0 | infeasible | 0 | 1 | 844,033 | 6,071 | 185 |
| 0.5 tau* | FO2r | direct | faces0 / facesS | infeasible / infeasible | – | 0 | 376,060 / 102,928 | – | 20 / 7 |
| 0.9 tau* | FO2r | stabP | stock / Rlight0 | infeasible / infeasible | 0 / 0 | 1 / 1 | 840,208 / 844,033 | 6,042 / 6,071 | 144 / 138 |
| 0.9 tau* | FO2r | direct | faces0 / facesS | infeasible / infeasible | – | 0 | 376,060 / 102,928 | – | 16 / 8 |
| 1.1 tau* | FO2r | stabP | stock / Rlight0 | infeasible / infeasible | 0 / 0 | 1 / 1 | 840,208 / 844,033 | 6,042 / 6,071 | 114 / 131 |
| 1.1 tau* | FO2r | direct | faces0 / facesS | infeasible / infeasible | – | 0 | 376,060 / 102,928 | – | 20 / 7 |
| 0.5 tau* | FO2s | stabP | stock | infeasible | 0 | 1 | 691,888 | 5,736 | 51 |
| 0.5 tau* | FO2s | stabP | Rlight0 | infeasible | 0 | 1 | 695,713 | 5,831 | 67 |
| 0.5 tau* | FO2s | direct | faces0 | optimal, eta 3.0e-9, relres 1.6e-5 (in) | – | 0 | 299,120 | – | 13 |
| 0.9 tau* | FO2s | direct | faces0 | optimal, eta 7.3e-9, relres 7.1e-5 (in) | – | 0 | 299,120 | – | 15 |
| 0.9 tau* | FO2s | direct | facesS | unknown (numerr 2) | – | 0 | 82,148 | – | 5 |
| 1.1 tau* | FO2s | direct | faces0 / facesS | infeasible / infeasible (out) | – | 0 | 299,120 / 82,148 | – | 11 / 4 |

Bisection of tau/tau* on [0.5, 1], FO2s, direct faces0 (`lb_i9c`): certified 0.75, 0.875, 0.90625,
0.91406, 0.91602; not certified 0.9375 (unknown, numerr 2), 0.92188 (unknown, numerr 2), 0.91797
(status 'optimal' from `lpi_classify_sop` with feasratio -1, relres 1.0, eta 1.0: counted out by the
gate). **Reach 0.91602 tau* = 2.639408.**

The FO2r form returns infeasible under every test, including stabP at 0.1 tau*, while FO2s at
0.5 tau* certifies in the direct test. stabP (form P, an L2-type storage on Tx) is infeasible on
both forms.

## Item 12 – heat with a dynamic (Laplace-Beltrami) edge

    v = pde_var(s1,[0,1]);  u = pde_var([s1;s2],[0,1;0,1]);
    v_t = v_s1s1 - u_s2(s1,1) + lam v;   u_t = u_s1s1 + u_s2s2 + lam u
    v_s1 = 0 at s1 = 0, 1;   u_s1 = 0 at s1 = 0, 1;   u(s1,0) = 0;   u(s1,1) = v(s1)     (sigma = 1)

**Reference [derived].** Cosine modes in s1. Mode n = 0: u = sin(k s2), v = sin(k), and the edge
equation gives k tan k = 1, so lam* = k0^2 = 0.740174 (k0 = 0.860334); modes n >= 1 subtract
(n pi)^2 from every eigenvalue. FD of the reduced modes n = 0..3 (second-order, one-sided flux at
the edge): 0.740247 (N = 40), 0.740207 (N = 60). The exact value is the reference.

| label | test | storage | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|
| 0.9 lam* | stabP | stock | optimal | 1 | 0 | 426,044 | 4,652 | 38 |
| 0.9 lam* | stabP | Rlight0 | optimal | 1 | 0 | 430,265 | 6,270 | 45 |
| 1.02 lam* | stabP | stock | unknown | 0 | 1 | 1,697,420 | 7,826 | 442 |
| 1.1 lam* | stabP | stock | unknown | 0 | 1 | 1,697,420 | 7,827 | 506 |

Bisection of lam/lam* on [0.5, 1], 8 steps, light (`lb_i12b_*`):

| storage | certified fractions | not certified | lam certified | reach |
|---|---|---|---|---|
| stock | all 8: 0.75, 0.875, ..., 0.99805 (dw 0, 426,044 variables, m 4,651 – 4,652, 29 – 38 s) | none | 0.738728 | 0.99805 (top of the bracket) |
| Rlight0 | all 8 (dw 0, 430,265 variables, m 6,244 – 6,271, 38 – 54 s) | none | 0.738728 | 0.99805 (top of the bracket) |

## Item 13 – road-field model on a box

    r = pde_var(s1,[0,1]);  f = pde_var([s1;s2],[0,1;0,1]);
    r_t = 2 r_s1s1 + f(s1,0) - r;   f_t = f_s1s1 + f_s2s2 + lam f
    r = 0 at s1 = 0, 1;  f = 0 at s1 = 0, 1 and s2 = 1;  f_s2(s1,0) = f(s1,0) - r(s1)
    (D = 2, d = 1, mu = nu = 1; lam the linearized growth rate f'(0))

**Reference [derived].** Sine modes in s1; lam enters the field only. lam* = the smallest lam with
max Re eig = 0 over the modes n = 1..3 (FD in s2 with a ghost node for the Robin condition; fzero in
lam): 13.919655 (N = 40), 13.919824 (N = 60); the minimizing mode is n = 1.

| label | test | storage | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|
| 0.9 lam* | stabP | stock | optimal | 1 | 0 | 317,716 | 3,495 | 12 |
| 0.9 lam* | stabP | Rlight0 | optimal | 1 | 0 | 319,545 | 4,025 | 11 |
| 0.9 lam* (rerun on the harness snapshot) | stabP | stock | optimal | 1 | 0 | 317,716 | 3,495 | 17 |
| 1.02 lam* | stabP | stock | infeasible | 0 | 1 | 508,228 | 4,765 | 104 |
| 1.1 lam* | stabP | stock | infeasible | 0 | 1 | 508,228 | 4,766 | 118 |

Bisection of lam/lam* (lam* = 13.919824, FD N = 60) on [0.5, 1], 8 steps, light (`lb_i13b_*`):

| storage | certified fractions | not certified | lam certified | reach |
|---|---|---|---|---|
| stock | all 8 (dw 0, 317,716 variables, 10 – 18 s) | none | 13.892636 | 0.99805 (top of the bracket) |
| Rlight0 | all 8 (dw 0, 319,545 variables, 13 – 18 s) | none | 13.892636 | 0.99805 (top of the bracket) |

## Item 5 – two-layer heat conduction on a rectangle pair

Both layers reflected onto the unit square, interface at s1 = 0:

    u1 = pde_var([s1;s2],[0,1;0,1]);  u2 = pde_var([s1;s2],[0,1;0,1]);
    u1_t = lap u1 + lam u1;   u2_t = 2 lap u2 + lam u2
    u1_s1(0,s2) = -(H/k1)(u2(0,s2) - u1(0,s2));   u2_s1(0,s2) = (H/k2)(u2(0,s2) - u1(0,s2))
    u1 = u2 = 0 at s1 = 1, s2 = 0, s2 = 1          (k1 = 1, k2 = 2, H = 1)

The converter accepts a boundary condition that couples the traces of two 2-D states on the same
edge. **Reference [derived].** Sine modes in s2 (mode n adds -k_i (n pi)^2); FD of the 1-D
transmission problem with ghost-node interface conditions, modes n = 1..3: lam* = 13.686807
(N = 40), 13.686938 (N = 60).

| label | test | storage | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|
| 0.9 lam* | stabP | stock | optimal | 1 | 1 | 751,936 | 11,911 | 149 |
| 0.9 lam* | stabP | Rlight0 | optimal | 1 | 1 | 754,180 | 12,455 | 158 |

Not bisected (the brief asks for one run below the threshold).

## Item 1 – wave-heat transmission on two rectangles (not exponentially stable)

Wave reflected onto the unit square, (u, v = u_t), heat w on the unit square, interface s1 = 0:

    u = pde_var([s1;s2],[0,1;0,1],[2;2]);  v = pde_var([s1;s2],[0,1;0,1],[1;2]);  w = pde_var([s1;s2],[0,1;0,1],[2;2]);
    u_t = v;  v_t = u_s1s1 + u_s2s2;  w_t = w_s1s1 + w_s2s2
    u(1,.) = 0;  u_s1(0,.) = -w_s1(0,.);  v(1,.) = 0;  w(1,.) = 0;  w(0,.) = v(0,.);  u = v = w = 0 at s2 = 0, 1

v has order 1 in s1 so that its trace at s1 = 0 exists with one condition in s1. The converter
accepts the system.

| label | test | storage | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|
| epneg 0 | stabP | stock | infeasible | 0 | 1 | 1,691,856 | 26,723 | 989 |
| epneg 0 | stabP | Rlight0 | stopped by hand after the stock row (projected past the 1,500 s cap) | – | – | 938,745 at dw 0 | – | – |

epneg 0.01 rows not run: a feasible program at epneg 0.01 is feasible at epneg 0 with the same P
(the decay term is nonnegative), so stock infeasibility at epneg 0 implies stock infeasibility at
epneg 0.01 at the same degrees (inferred). The verdict is uninformative for the question of the
catalog: the (u, u_t) form leaves the stock P form without a certificate even for the conservative
part (see the reading).

## Item 16 – age-structured population with diffusion

    q = pde_var([s1;s2],[0,1;0,1],[1;2]);             % s1 = age a, s2 = space x
    q_t = -q_s1 + 0.1 q_s2s2 - q;   q(0,s2) = b int_0^1 q(a,s2) da;   q(a,0) = q(a,1) = 0

**Reference [derived].** Sine modes in x; the first mode decays iff R_eff = b (1 - e^{-c}) / c < 1,
c = mu + d pi^2 = 1.986960: b* = 2.302685.

| label | test | storage / terms | tier | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|---|
| 0.1 b* | stabP | stock | light | optimal | 1 | 0 | 103,744 | 2,055 | 5 |
| 0.5 b* | stabP | stock | light | optimal | 1 | 0 | 103,744 | 2,045 | 6 |
| 0.9 b* | stabP | stock | light | infeasible | 0 | 1 | 187,984 | 2,967 | 15 |
| 0.9 b* | stabP | Rlight0 | light | infeasible | 0 | 1 | 188,545 | 2,983 | 12 |
| 0.9 b* | stabP | stock | heavy | optimal | 1 | 0 | 751,849 | 5,120 | 74 |
| 0.5 b* | direct | faces0 / facesS | – | optimal, eta 5.7e-10 / 2.9e-10 (in / in) | – | 0 | 210,340 / 56,548 | – | 9 / 4 |
| 0.9 b* | direct | faces0 / facesS | – | optimal, eta 6.9e-11 / 2.2e-10 (in / in) | – | 0 | 210,340 / 56,548 | – | 9 / 5 |
| 1.1 b* | stabP | stock | light | infeasible | 0 | 1 | 187,984 | 2,970 | 10 |
| 1.1 b* | stabP | Rlight0 | light | infeasible | 0 | 1 | 188,545 | 2,984 | 11 |

Not bisected.

## Item 10 – heat with a delayed Robin boundary term (Nicaise, Valein, Fridman 2009)

    u = pde_var(s1,[0,pi]);  y = pde_var(s2,[0,1],1);
    u_t = u_s1s1;   y_t = -(1/tau) y_s2
    u(0) = 0;   u_s1(pi) = -mu0 u(pi) - mu1 y(1);   y(0) = u(pi)          (mu0 = 1, tau = 1)

The paper gives sufficient conditions only; the catalog gives no parameter values; mu0 = 1,
mu1 in {0.5, 2}, tau = 1 were chosen here. **Reference** (method of lines, heat FD with a ghost node,
upwind transport): spectral abscissa -0.762312 / -0.762261 (mu1 0.5; N, M = 100, 200 / 200, 400),
-0.206656 / -0.198679 (mu1 2). Both stable; no threshold computed.

| label | test | storage | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|
| mu1 0.5 | stabP | stock / Rlight0 | optimal / optimal | 1 / 1 | 0 | 21,024 / 22,480 | 142 / 145 | 1 / 0 |
| mu1 2 | stabP | stock / Rlight0 | optimal / optimal | 1 / 1 | 0 | 21,024 / 22,480 | 142 / 145 | 0 / 0 |

## Item 15 – elastic-thermoelastic transmission, 1-D

Elastic part (e1 = u_x, p = u_t) on (0, 1/2), thermoelastic part (e2 = v_x, q = v_t, th) on (1/2, 1)
shifted to (0, 1/2); a = b = k = alpha = beta = 1:

    e1_t = p_s;  p_t = e1_s;  e2_t = q_s;  q_t = e2_s - th_s;  th_t = th_ss - q_s
    e1(0) = 0;  q(1/2) = 0;  th(0) = th(1/2) = 0;  p(1/2) = q(0);  e1(1/2) = e2(0)

Modification: free end at x = 0 (e1(0) = 0) instead of the clamped end of the paper. With both ends
clamped the Riemann form keeps the uniform-strain mode (e1, e2 constant, p = q = th = 0) as an
equilibrium, so its PIE has a zero eigenvalue. The exponential stability of the modified system is
not established here (the paper treats the clamped case; qualitative only).

| label | test | storage | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|
| TE | stabP | stock / Rlight0 | optimal / optimal | 1 / 1 | 0 | 7,000 / 7,000 | 763 | 3 / 2 |

## Item 2 – wave-heat transmission, 1-D (not exponentially stable)

Wave on (-1, 0) reflected onto (0, 1) (Dirichlet variant of Sec. 5 of the paper), heat on (0, 1):

    e = pde_var(s1,[0,1],1);  v = pde_var(s1,[0,1],1);  w = pde_var(s1,[0,1]);
    e_t = v_s1;  v_t = e_s1;  w_t = w_s1s1;   v(1) = 0;  e(0) = -w_s1(0);  w(1) = 0;  w(0) = v(0)

| label | test | storage | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|
| epneg 0 | stabP | stock / Rlight0 | optimal / optimal | 1 / 1 | 0 | 2,520 / 2,520 | 284 | 1 / 1 |
| epneg 0.01 | stabP | stock | infeasible | 0 | 2 | 9,000 | 478 | 5 |
| epneg 0.001 | stabP | stock | infeasible | 0 | 2 | 9,000 | 478 | 4 |

A non-increasing P-norm is certified; no decay rate down to 0.001 is certified, consistent with the
polynomial (not exponential) decay of the paper (no exactness claim: a negative verdict holds at
these degrees only).

## Item 3 – wave with an acoustic boundary condition

    ph = pde_var([s1;s2],[0,1;0,1],[2;2]);  ps = pde_var([s1;s2],[0,1;0,1],[2;1]);
    de = pde_var(s1,[0,1],2);  et = pde_var(s1,[0,1],2);
    ph_t = ps;  ps_t = lap ph;  de_t = et;  et_t = -et - de - ps(s1,0)        (m = d = k = rho = 1)
    ph = ps = 0 at s1 = 0, 1;  ph_s2(s1,0) = -et;  ph(s1,1) = 0;  ps(s1,1) = 0;  de = et = 0 at s1 = 0, 1

**Conversion error** with the edge states of order 0 in s1 (`de = pde_var(s1,[0,1],0)`, the model of the
catalog, which has no s1-derivative on the edge), verbatim:

    The order of the derivative of x_4 is higher than that to which the state is differentiable.

The edge state et is the Neumann datum of ph on s2 = 0 and must be as differentiable in s1 as ph.
Rebuilt with de, et of order 2 in s1 and de = et = 0 at s1 = 0, 1 (compatible with ph = 0 there; the
dynamics preserve these conditions because ps(0, 0) = ps(1, 0) = 0). That PIE converts.

| label | test | storage | status | cert | dw | ndv | m | t |
|---|---|---|---|---|---|---|---|---|
| epneg 0 | stabP | stock | infeasible | 0 | 1 | 2,032,912 | 19,034 | 737 |
| epneg 0 | stabP | Rlight0 | stopped by hand (projected past the cap) | – | – | 1,270,864 at dw 0 (stock) | – | – |

## Items 14, 17, 7 – not built

- Item 14 (transport-diffusion couplings, non-coercive Lyapunov functionals): qualitative, no
  open-loop threshold; not built (brief).
- Item 17 (3-D heat-wave, non-box geometry): not built; the budget was spent (brief: only if time
  remains).
- Item 7 (wave with delayed damping) is not in the priority list of the brief and was not built.

## Harness incident (03:30 – 03:45)

At 03:34 – 03:40 another session changed `psatz_faces\pf_rows.m` and `pf_row_storage.m` (opt-in
`_x<c>` rows) and added `psatz_faces\shadow\get_lift_degs.m` (opt-in option `wcross`). The
assertion in `pf_setup` then rejected the path. Every `pf_exec_row` call failed before building a
program: the lb_i4a Rlight0 row and four lb_i4b_stock invocations. These are recorded as CONV_ERROR
rows in `results\ex_i4.mat`. The message is the assertion text, not a conversion error, and those
rows were rerun; they are not results. From 03:42 on, every run used a frozen snapshot
`litbuild\harness\`. It holds copies of `pf_setup`, `pf_exec_row`, `pf_rows`, `pf_row_storage`,
`pf_terms`, `shadow\lpi_ineq_sop`, `shadow\poscopvar_direct`, `shadow\private\`,
`heat\heatNd_lpi` and `route_row\exec_tools_sop`, without the `get_lift_degs` shadow.
`psatz_faces` is off the path, so `get_lift_degs` resolves to the repository. Check: RF at 0.9 lam*
on the snapshot reproduces the earlier program and verdict (317,716 variables, m 3,495, optimal).
The item 13 Rlight0 bisection (03:37 – 03:39) ran on the modified `pf_rows`/`pf_row_storage` but
before the `get_lift_degs` shadow existed. Its programs have the size of the 03:00 Rlight0 run
(319,545 variables). No other run fell in 03:30 – 03:45.

## Summary

| item | built | threshold known | stock verdict | Rlight0 verdict | reach stock | reach Rlight0 | direct test |
|---|---|---|---|---|---|---|---|
| 8, 1-D control (S18) | yes | tau* 1.034726 [derived, exact modal] | 0.5 tau*: infeasible (dw 1) | infeasible (dw 1) | – | – | faces0 in (0.5 tau*) |
| 8, 2-D (3-D transport) | yes | tau* 1.034726 [derived] | 0.5 tau*: infeasible (dw 0, 1.67M) | infeasible (dw 0, 1.68M) | – | – | not solved: 47.4M > cap |
| 4 Turing (Neumann) | yes | d* 4.461640 [derived, exact] | certified 0.9 – 0.9875 d* (2.10M, 626 – 909 s) | 0.9 d*: not completed (2.12M, killed at 1,500 s twice) | >= 0.9875 (no uncertified step measured) | not run | not run |
| 11 exact gain | yes | gamma* 0.04571818 [closed form, FD-checked] | 0.0459440 light, 0.0459432 heavy (+0.49 %) | 0.0459459 (+0.50 %) | – | – | – |
| 6 XYL wave, boundary delay | yes | mu* 1/2 [exact, k 1, tau 2] | infeasible 0.40 – 0.49; certified 0.5 – 0.75; rate 0.10 certified, 0.11 not (exact 0.10597) | same | 0.500391 (top) | 0.500391 (top) | out 0.40 – 0.49; in 0.5, 0.6, 0.75; bisection steered by eta: 0.575391 |
| 9 delayed damped wave | yes, two forms | tau* 2.881401 [derived] | infeasible at every tau, both forms | infeasible, both forms | – | – | FO2s faces0: in to 0.91602 tau*, out at 1.1; FO2r out everywhere |
| 12 Laplace-Beltrami edge | yes | lam* 0.740174 [derived, exact] | certified to 0.99805; UNKNOWN at 1.02, 1.1 | certified to 0.99805 | 0.99805 (top) | 0.99805 (top) | not run |
| 13 road-field | yes | lam* 13.919824 [FD N = 60] | certified to 0.99805; infeasible at 1.02, 1.1 | certified to 0.99805 | 0.99805 (top) | 0.99805 (top) | not run |
| 5 two-layer heat | yes (two 2-D states, shared edge) | lam* 13.686938 [FD N = 60] | 0.9: certified (dw 1) | certified (dw 1) | not bisected | not bisected | not run |
| 1 wave-heat rectangles | yes | none (not exp. stable) [read] | epneg 0: infeasible (dw 1) | stopped by hand | – | – | not run |
| 16 age + diffusion | yes | b* 2.302685 [derived] | light: certified 0.1, 0.5 b*, infeasible 0.9, 1.1; heavy: certified 0.9 | light: infeasible 0.9, 1.1 | not bisected | – | in at 0.5 and 0.9 (faces0, facesS) |
| 10 delayed Robin heat | yes | none; FD abscissa -0.762 (mu1 0.5), -0.199 (mu1 2) | certified both | certified both | – | – | not run |
| 15 thermoelastic | yes (free end at x = 0) | none (qualitative) | certified (epneg 0) | certified | – | – | not run |
| 2 wave-heat 1-D | yes | none (polynomial decay) [read] | epneg 0 certified; 0.01, 0.001 infeasible | epneg 0 certified | – | – | not run |
| 3 acoustic BC | yes, edge states raised to order 2 (order 0 rejected) | none | epneg 0: infeasible (dw 1, 2.03M) | stopped by hand | – | – | not run |
| 7 | no (not in the brief's list) | – | – | – | – | – | – |
| 14 | no (qualitative; brief) | – | – | – | – | – | – |
| 17 | no (budget; brief: only if time remains) | – | – | – | – | – | – |

"top": every bisection step certified, so the reach is the largest fraction the bracket allows.

## Reading

Measured.

- The converter accepted every system built:
  - a 2-D state as the face datum of a 3-D state (item 8);
  - a vector-valued 2-D state with Neumann conditions (item 4);
  - a 1-D input in a 2-D equation with a 2-D output (item 11);
  - two 2-D states coupled through traces and normal derivatives on a shared edge (items 1, 5);
  - an integral renewal condition on a mixed-order 2-D state (item 16);
  - a Volterra integral of a 2-D transport state (item 9, FO2r);
  - 1-D states in different variables coupled at points (items 6, 10).

  One system was rejected: a 1-D edge state of order 0 as the Neumann datum of a 2-D state
  (item 3).
- Item 8: the 1-D control reproduces the suite. Stock light is infeasible at 508,228 variables,
  and the direct test with faces0 certifies at 227,510. In 2-D, stabP is infeasible at 0.5 tau*
  under both storages: 1.67M and 1.68M variables, 476 s and 653 s, no raise taken. The direct-test
  program has 47,433,948 variables (350 s and about 38 GB to build) and was not solved.
- Item 4: stock light certifies 0.9, 0.95, 0.975 and 0.9875 d* at 2,102,408 variables, 626 – 909 s
  per solve. Rlight0 did not finish one solve within 1,500 s (twice).
- Item 11: the certified gains are 0.0459432 – 0.0459459, 0.49 – 0.50 % above
  gamma* = 0.0457182. Stock heavy lowers the stock light bound by 8.4e-7 with 4.0 times the
  variables.
- Item 6:
  - stabP separates mu = 0.49 (infeasible) from mu = 0.5 (certified with epneg 0) under both
    storages.
  - At mu = 0.6 stabP certifies the decay rate 0.10 and refuses 0.11 (exact 0.105968).
  - The direct test also separates 0.49 from 0.5. On the stable side its eta is 1e-9 – 3e-8,
    so the 1e-8 gate steered its bisection to 0.575.
- Item 9:
  - The Riemann form FO2r is infeasible under every test, including stabP at 0.1 tau*.
  - On the S18-type form FO2s, the direct test with faces0 certifies up to 0.91602 tau* (2.639408)
    and is infeasible at 1.1 tau*.
  - stabP is infeasible on FO2s as well.
- Items 12 and 13: every bisection step certified, up to 0.99805 of the threshold, under both
  storages. Above the threshold, item 13 is infeasible at 1.02 and 1.1 lam*. Item 12 exits UNKNOWN
  after a raise (1,697,420 variables).
- Item 16: stabP light fails at 0.9 b* under both storages and certifies at 0.1 and 0.5 b*. Stock
  heavy (751,849 variables) and the direct test both certify at 0.9 b*.
- Items 1 and 3 (2-D waves in the (u, u_t) form): stabP stock is infeasible after a raise at
  epneg 0.
- Rlight0 against stock:
  - The verdicts agree on every row where both finished (items 2, 5, 6, 8, 9, 10, 11, 12, 13, 15,
    16).
  - Rlight0 programs have 0.1 – 7 % more variables and up to 35 % more rows (item 12: 6,270
    against 4,652).
  - Item 4 is the only case where Rlight0 did not finish: more than 1,100 s of solve, against
    740 s for stock.

Inferred, not measured.

- stabP bounds Tx in a P-weighted L2 norm. For FO2s, BPS2 and AC, Tx = (z, z_t, ...) in L2, which
  omits ||z_s||^2. A high-frequency family with z = eps z_t makes the derivative of any such
  functional positive (argument sketched, not proven). That would explain why stabP fails on every
  wave in this form while the direct test certifies FO2s: the quadratic form of the direct test
  acts on the fundamental state.
- FO2r carries the energy as its L2 norm and still fails everywhere. Its Volterra delay term is the
  remaining difference from FO2s; this was not isolated.
- The item 1 and item 3 verdicts therefore do not test the question of the catalog (absence of
  exponential stability). The stock P form does not certify these representations even at
  epneg 0.
- The item 12 UNKNOWN above the threshold is a solver exit. It is neither a certificate nor an
  infeasibility certificate.

Skipped, and why.

- Item 8, 2-D at 0.9 tau* (four runs) and the 2-D direct test: the direct-test program is 47.4M
  variables, above the 3M cap, and the 0.5 tau* runs were negative at 8 – 11 min each.
- Item 4: bisection step 4 (0.99375) and the whole Rlight0 bisection were not run; the time went
  to the gap study (maintainer, 04:30). The stock bisection had already been cut to 4 steps on
  [0.9, 1], because one solve costs 626 – 909 s.
- Items 1 and 3, Rlight0 rows: stopped by hand after the stock row, because they were projected to
  run past the cap. The epneg 0.01 rows of items 1 and 3 were not run; the epneg-0 result implies
  infeasibility at the same degrees.
- Items 5 and 16: not bisected (the brief asks for one run below the threshold).
- Direct test not run on items 1 – 5, 10 – 13 and 15.
- Items 7, 14 and 17 not built.
