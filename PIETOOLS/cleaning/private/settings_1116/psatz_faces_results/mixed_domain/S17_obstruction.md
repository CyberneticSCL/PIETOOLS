# S17: no coercive PI Lyapunov operator exists (argument, 2026-11-16)

Derived by MMP (session af90e0a1) on 11/16/2026 from the structure of S17; not yet
checked by a second reader. Numerical corroboration (M7, `pf_box.m`): with the scale of P
fixed by ep I <= P <= I, MOSEK returns an infeasibility certificate for S17 and for the
heat-source variant D4 at ep = 0.5 and 0.1 (UNKNOWN at 0.01), while S3 is feasible at
ep = 0.1 and 0.01. The polynomial state families of `pf_check_fam.m` cannot expose the
executive's flags-only certificates, whose violation is scaled by the margin eppos = 1e-4
(a layer of width about eppos; see `S17_S18_findings.md`).

## Setting

S17 is the cascade

    x1_t = -x1_{s1},                 x1(0,t) = 0,                 s1 in [0,1],
    x2_t = x2_{s2 s2} + r2 x2,       x2(s1,0,t) = x1(s1,t),  x2(s1,1,t) = 0,

with x2 of order 0 in s1 (no s1-derivative in its equation and no boundary condition in
s1). The state is x = (x1, x2) in X = L2(0,1) x L2((0,1)^2), the domain of the generator
is the set of x with x1 in H^1, x1(0) = 0, x2(s1, .) in H^2(0,1) for almost every s1,
x2(s1,0) = x1(s1) and x2(s1,1) = 0. The fundamental state is xh = (x1_{s1}, x2_{s2 s2})
in L2 x L2, and x = T xh.

A bounded self-adjoint operator P on X is a **PI operator** when its blocks are finite
sums of the PIETOOLS basis-operator products with polynomial kernels:

    (P11 x1)(s1)     = p(s1) x1(s1) + integrals of x1 against polynomial kernels,
    (P12 x2)(s1)     = int_0^1 q(s1,t2) x2(s1,t2) dt2 + int int k12(s1,t1,t2) x2(t1,t2) dt1 dt2,
    (P22 x2)(s1,s2)  = r(s1,s2) x2(s1,s2) + operators integrating x2 in s1, in s2, or in both,

with p, q, k12, r and every kernel polynomial. P is **coercive** when P >= eps I on X for
some eps > 0; for a PI operator this forces r(s1,s2) >= eps on [0,1]^2, since the
integral parts are compact. The Lyapunov derivative of V(x) = <x, P x> along the
cascade is dV/dt = 2 <P x, A x> on the domain, a bounded quadratic form in xh.

**Proposition.** For every r2 in R and every coercive PI operator P on X, there is a
state x in the domain of the generator of S17 with dV/dt > 0.

In words: the stability LPI of S17 is infeasible in the PIETOOLS class at every degree,
and the executive's UNKNOWN exits are the solver's reading of a program without a
feasible point, not a degree or term-set limit. The argument never uses r2, so the
obstruction is about the Dirichlet edge coupling from a state that is not smoothed in
s1, not about stability.

**Proof (direct, two state families).** Write P22 = r + K22 with K22 the integral
parts, and P12 as above with q the kernel of its part that is local in s1.

*Step 1: the trace at s2 = 0 forces q(s1,0) <= -r(s1,0).* Fix a smooth g with g(0) = 0
and a smooth psi on [0,1] with psi(0) = 1, psi(1) = 0, psi'(1) = 0, and for eta in (0,1)
take x1 = g, x2(s1,s2) = g(s1) psi(s2/eta) for s2 <= eta and 0 beyond. This x lies in the
domain. Expanding dV/dt in powers of 1/eta, the terms of order 1/eta come only from the
multiplier r of P22 and from the adjoint of the local part of P12, because every other
block integrates x2 in s2 and its contribution is O(1):

    dV/dt = (2/eta) [ |psi'(0)| int g^2 (r + q)(s1,0) ds1 - (int_0^1 psi'^2) int g^2 r(s1,0) ds1 ]
            + (2/eta) |psi'(0)| <g, K_0 g> + O(1),

where K_0 is the compact operator on L2(0,1) with kernel k12(s1,t1,0). The ratio
|psi'(0)| / int psi'^2 is unbounded over the admissible psi (a profile with slope -K on
[0, alpha/K] and the constant slope -(1-alpha)/(1-alpha/K) after it has ratio about
1/alpha). Hence dV/dt <= 0 for every eta and psi requires

    int g^2 (r + q)(s1,0) ds1 + <g, K_0 g> <= 0     for every admissible g.

Taking g a narrow bump of unit norm at s1* in (0,1] sends <g, K_0 g> to 0 and gives
r(s1*,0) + q(s1*,0) <= 0, so q(s1,0) <= -r(s1,0) <= -eps on (0,1].

*Step 2: the local cross block cannot be dominated when x2 has no s1-diffusion.* Let
phi_n be an orthonormal family in L2(0,1) with phi_n -> 0 weakly (shifted Legendre
polynomials, or sqrt(2) sin(n pi s1)), let psi be smooth with psi(0) = psi(1) = 0, and
for c in R take

    xh1 = phi_n,   x1 = T phi_n,   x2 = x1(s1)(1 - s2) + c phi_n(s1) psi(s2),

so that x2 = T2 x2_{s2 s2} + T21 xh1 and x lies in the domain; ||x1|| -> 0 as n -> inf.
Then

    dV/dt = 2 <x1_t, P12 x2> + O(c^2) + o(1) = -2c int_0^1 int_0^1 phi_n(s1)^2 q(s1,s2) psi(s2) ds2 ds1 + O(c^2) + o(1),

because x1_t = -phi_n, the integral-in-s1 part of P12 applied to phi_n psi vanishes
weakly, every term of P11 and of the adjoint block carries a factor x1 or x1_t paired
with x1 and is o(1), and every term quadratic in the oscillating part of x2 is O(c^2);
no term is of order c n^p with p > 0 because the generator does not differentiate x2 in
s1. By Step 1, q(s1,0) <= -eps, so for psi >= 0 supported near s2 = 0 the double
integral is negative and bounded away from zero uniformly in n (phi_n^2 converges weakly
to a density of unit mass). Choosing the sign of c against it and |c| small, then n
large, gives dV/dt > 0. QED

**Status (11/16 04:00).** The proposition is about the direct 'P' form, V = <x, P x> a bounded form in the PDE state. The direct PIE2PDE form (Q free, T'Q symmetric, T'Q >= eps T'T) certifies S17 at 0.5 r2* (`parameterization.md`): T'Q need not factor as T'PT, and -2 <x1_s1, Q11 x1_s1> in A'Q + Q'A supplies the dissipation in the fundamental state that Step 2 shows no bounded form can. The sentence in `S17_S18_findings.md` extending the argument to the Q form was wrong.

**Where the argument stops.** For S3 and for every cascade whose target diffuses in s1,
the term 2 <x2_t, P22 x2> carries -2 c^2 <r phi_n' psi, phi_n' psi>, of order c^2 n^2 or
larger, and the quadratic in c of Step 2 stays nonpositive; Step 1 then only requires
the local cross block, which the storages carry. The argument applies as stated to the
heat-source variant D4 (x1 heat, x2 of order 0 in s1), with x1_t = +phi_n in Step 2.

**Consequence for the executives.** A solver exit on S17 or D4 that the executive counts
as certified is a numerical artefact of a program without a feasible point: the operators
returned carry no x2 block beyond the margin eppos, and the positive value of their
derivative form sits on states with a boundary layer of width about eppos, which the eta
gate does not see (M2, CHK2).
