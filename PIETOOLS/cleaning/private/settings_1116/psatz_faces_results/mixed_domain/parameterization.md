# Which positive operators the parameterizations miss on mixed domains (2026-11-16)

Maintainer, 11/16: treat the S17 and S18 failures as a parameterization defect or as the
need for a weaker certificate; PIE2PDE stability is the default test; first find fixed
operators whose positivity cannot be certified. Tables: `fixed_operators.md` (the slack
reader on fixed targets), `m2_chk_summary.md` (the Q form as built, the storage cones,
the direct PIE2PDE form). Scripts: `harness/pf_fixed.m`, `pf_cone.m`, `pf_qdirect.m`,
the `pf_run_*.m` batches, `qform/` (the executive with the Q forms allowed in N-D).

## 1. The slack reader certifies every fixed target tried

The PDE-norm forms T'T, T'E1T, T'E2T and A'A of S17, S3 and S18, and the certified S18
storage P*, are certified positive by `lpi_ineq_sop` sized from their own support (light
tier, with or without terms, and at raised degrees): every row 'optimal'. The slack
parameterization – every basis-operator pattern, lower and upper integrals included,
degrees read from the support – is not where the gap is.

## 2. The storage cone is: T'T lies outside every separable storage

| fixed X | uniform, light | uniform, raised | four separable patterns, light | four, raised | every pattern with split kernels, raised |
|---|---|---|---|---|---|
| T'T of S17, S3, S18; T'T of S17's x2 alone | out | out | out | out | in |
| S18 certified storage P* | out | out | out | in | in |
| S3 certified storage | in | in | in | in | in |

"out" is a MOSEK infeasibility certificate for X = P over the storage's Gram variables.
The PDE norm itself – the operator T'T that PIE2PDE stability compares P against – is not
in any storage the tiers build in two variables (uniform patterns, no incomplete moments:
d1 = 0 in the proposal's row), nor in the four-pattern family that certifies S18. It needs
the incomplete-moment patterns: lower and upper integrals, the (c) factor of the proposal.

## 3. The Q form as the executive builds it fails on the control

With a storage P and the equality T'Q = P (the 1-D construction, allowed in 2-D by the
shadow), S3 is infeasible at dw 0 and UNKNOWN at dw 1 although it certifies in the 'P'
form with the same storage; S17 and S18 are infeasible with the uniform storage and with
the four-pattern storage. The equality is the obstacle: T'Q carries the Volterra kernels of
T and of Q, and the storage cannot hold them (Section 2).

## 4. The direct PIE2PDE form certifies S3, S17 and S18

Corollary 7.1 of the hierarchy paper asks for Q with QT symmetric, QT >= eps T*T and
QA + A*Q* <= 0. Built directly (`pf_qdirect.m`): Q free with the shape of T
(LPIVAR_CDOPVAR, every cap at the degree shown), the symmetry T'Q = Q'T as an equality,
the coercivity T'Q >= eps T'T and the decay -(A'Q + Q'A) >= 0 as inequalities certified by
the slack reader from their own supports; no storage operator at all.

| system | Q degree 1, dw 0 | degree 1, dw 1 | degree 2, dw 0 | degree 2, dw 1 | 1.1 x margin | eps = 0.1 | 0.9 x margin |
|---|---|---|---|---|---|---|---|
| S3, 0.9 r* | certified, 383k var., 15 s | certified | certified | certified | – | – | – |
| S17, 0.5 r2* | certified, 435k, 22 s | UNKNOWN | certified, 1.16M, 82 s | certified, 2.10M, 257 s | infeasible | certified | UNKNOWN |
| S18, 0.5 tau* | certified, 228k, 11 s | certified | certified | certified | infeasible | certified | UNKNOWN |

Every "certified" is MOSEK optimal with eta below 1e-8; the unstable cases are
infeasible with a certificate; the certificates stand at eps = 0.1, which the 'P'-form
artefacts of S17 did not (M7). So S17 is certifiable, in the default notion of the
toolbox: the obstruction of `S17_obstruction.md` is a statement about the direct 'P'
form (V a bounded quadratic form in the PDE state), and the Q form escapes it because
T'Q, symmetric, need not factor as T'PT: the (1,1) block of A'Q + Q'A carries
-2 <x1_s1, Q11 x1_s1>, a dissipation in the transport's fundamental state that no
bounded form in x1 provides. The conjectured extension of the obstruction to the Q form
in `S17_S18_findings.md` was wrong and is withdrawn there.

## 5. What this says about the parameterization

1. The slack (every pattern, degrees from the support) certifies every fixed positive
   operator tried, and both inequalities of the direct PIE2PDE form at Q of degree 1.
2. The storage is the defect: the PDE norm T'T lies outside every separable storage, so
   a PIE2PDE form that routes QT through a storage equality cannot work in two variables
   (Section 3), and the 'P' form with the uniform storage misses S18 (needs the mixed
   patterns at degree 2) and can never certify S17 (no form of that class can).
3. For the default test the storage is not needed: pose QT >= eps T*T as an inequality on
   the free Q and let the reader size it. The programs are of the size of the 'P' form
   (S18: 228k against 767k variables, 11 s against 64 s) and reach further (S17).
4. Open: the reach at the margin (0.9 r2*, 0.9 tau* UNKNOWN at Q of degree 1; degree 2
   and dw 1 not run there), the UNKNOWN of S17 at degree 1, dw 1, and the executive's
   restriction of the Q forms to 1-D, which this construction replaces rather than lifts.
   Nothing in the repository changed; `qform/` and `pf_qdirect.m` are scratch.

## 6. Is there a positive PI operator whose positivity cannot be certified?

Maintainer, 11/16 04:20. Two answers, one proven and one measured.

Proven, in one variable (09/28/2026, DEMO3): for every gamma and every length, the Poincare
operator gamma H1'H1 - H2'H2 lies outside the certificate cone without the box weight, and
inside it with the weight from length 3 (exact rational points). So a strictly positive PI
operator with no plain-term certificate at any length exists; the weights are not a
refinement but a necessity. With the weights, containment above the sharp constant at a finite
length is conjectured (the necessity conjecture, 1-D; assumed by maintainer decision 10/07).

Measured, on the mixed domains (`m2_chk_summary.md`, POINC): the PDE-norm lift
gamma I - T'T, strictly positive for gamma > ||T||^2, is the mixed-domain Poincare operator.

| PIE | lam (Rayleigh, exact where known) | F = none, up to 4 lam | faces at light | faces at raised lengths |
|---|---|---|---|---|
| S17 | 0.5420 | UNKNOWN | reach 1.0117 (resolution) | 1.0117 |
| S3 | 0.01370 | infeasible, certificate | 1.0117 | 1.0117 |
| x2 of S17 alone, order (0,2) | 0.010266 = 1/pi^4 | infeasible, certificate | 1.0117 | 1.0117 |
| 1-D transport, x(0) = 0 | 0.405285 = 4/pi^2 | infeasible, certificate | 1.0117 | 1.0117 |

Without a term, four times the norm is not certifiable, with MOSEK certificates on three of the
four; with the faces (in 1-D the box factor) the reach is the resolution of the bisection. The
sharp-constant run (bracket [lam, 1.1 lam], 12 steps, resolution 2.4e-5) is appended below
when it ends. So on mixed domains the operators that cannot be certified at the reader's
lengths are, as in 1-D, those at the boundary of the cone (semidefinite, with the sharp
constant) and those posed without terms; no strictly positive fixed operator with terms has
failed at the reader's lengths in any test of 11/16. The operators the STORAGE cannot
represent are a different and larger class (Section 2: T'T itself without incomplete moments).
