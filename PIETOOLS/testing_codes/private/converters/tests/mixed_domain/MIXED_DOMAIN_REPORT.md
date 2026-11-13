# Coupled PDEs on mixed spatial domains: stress test of PIETOOLS

11/12/2026. Systems whose states live on spaces of different dimension –
R, L2[s1], L2[s2], L2[s1,s2] and L2[s1,s2,s3] – tested in four layers: the
parser and converter, the `copvar` operations on the resulting PIE operators,
the container executives, and PIESIM. The opvar classes never represented
these systems, so every baseline is independent of PIETOOLS: the PDE written
by hand, exact answers by separation of variables, a Chebyshev collocation of
the PDE itself, `dde23`, and the literature where it exists.

**Tree.** Layers 1, 2 and 4 and the exact baselines ran on a frozen export of
the committed `ndopvar` HEAD 6c6f60dd (11/12/2026 01:37), MATLAB R2025b,
MOSEK 11.0. The shared checkout changed under the runs (other sessions edited
`settings_sopvar/`, `lpi_programming_sopvar/` and `executives_sopvar/` from
02:21, and `poscopvar_direct.m` failed to parse at 02:31), so every number
below comes from the export unless marked *live*.
The executive layer ran on the same export. Two further file sets add the
uncommitted psatz corner-row files of another session to the export, and the
rows that use them say so: set A, copied 03:40 from the live checkout
(SHA-1 prefixes `poscopvar_direct.m` cce8a04c0a, `get_lift_degs.m`
5fc41d4f87, `lpi_ineq_sop.m` 5fd6f1dff8), and set B, copied 04:34, which
differs from set A in `get_lift_degs.m` only (2cf69eda27, corrected corner
offsets). Both enable the rows by `st.lpi.ineq.psatz = 'corners'`.

**Files** (`testing_codes/private/converters/tests/mixed_domain/`):
`test_mixed_domain_pde.m` (the test, four layers), `mixed_domain_systems.m`
(the catalogue, with the hand-written semantics of each system),
`mixed_domain_colloc.m` (the collocation reference), this report.

## Findings

1. **Converter: no wrong conversion found.** 25 accepted systems and the
   control C1 reproduce their PDE to at most 2.5e-14, including an ODE ->
   1-D -> 2-D chain, a
   dynamic (Wentzell) edge, 1-D states on two edges, a 1-D edge -> 2-D face
   -> 3-D chain, Robin edge data, vector-valued edge states, non-unit
   domains, reversed variable order, and library PDEs 42 and 43. Six
   ill-posed systems are rejected with messages that name the cause,
   including a corner condition on a derivative and an edge clash in 3-D.
2. **Converter limitation, misleading message** (D2). A higher-dimensional
   state in the BC of a lower-dimensional state, `x1(0) = int int x2`, is
   rejected as a BC that "fixes the variables s1, s2 simultaneously (e.g. at
   a corner or edge)". The PIE exists; no corner is involved.
3. **Library defect** (D1). `PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay.m`
   (library PDE 43) codes `+a1*x(t-tau)`; its header and the library table
   have `-a1` and cite Caliskan and Ozbay (2009), whose text was not read.
   As coded, the system is unstable for every tau, so the executives'
   "infeasible" on it, recorded earlier as conservatism, is the correct
   verdict.
4. **PIESIM fails on two classes of PIE whose conversion is correct.**
   (D3) Any PIE with an auxiliary state (a direction with Neumann or periodic
   conditions at both ends, so PIE.K is nonempty): 4 of 4 crash or return
   NaN, down to a single 2-D heat state with Neumann conditions in s1.
   (D4) Integral BCs in 2-D: NaN, down to a single 2-D heat state with
   `int x ds1 = 0`; this is why PIESIM 2-D library examples 8 and 9 fail.
   Both examples blew up in the Phase 0 record of 11/02 (max |x| 2.0e41 and
   5.6e19), unflagged because that harness compared no 2-D exact solution.
   On the eight mixed-domain PIEs tested without these two features
   (library 2-D examples 10 and 11, S3, S4, S5, S8, S25 with an input, PDE 43
   with the sign corrected), PIESIM agrees with the references to 1e-8, and
   to 3.4e-7 on the growing delay mode of PDE 43 as coded.
5. **copvar operations: no defect.** Composition, adjoint, sum and the
   Lyapunov term `T'A + A'T` agree with their definitions to 1.4e-14 on 25
   mixed-space PIEs and on C1.
6. **Executives: no unsound result.** In the 41 light cases of the test (37
   in the run of record, 4 on the control C1) and in 68 further executive
   calls (light, heavy, veryheavy, with and without corner rows; 16 of them
   at 0.98 or 1.02 of a margin, or a = +-0.02 for the ODE cascades), no
   unstable PIE was certified and no certified gain or H2 norm lies below
   its exact value.
7. **Stability with a lower-dimensional state as Dirichlet data on an edge
   is conservative, by an amount that depends on the source.** At light,
   8 of 12 stable cases certify at dw 0 – ODE cascades, two-way edge
   couplings (Wentzell, edge flux), in-domain couplings, PDE 42.
   These eight also certify just inside the margin (0.98 of it; a = 0.02
   for the ODE cascades), and none is certified just beyond it (1.02;
   a = -0.02). With a diffusing 1-D state as the edge data, the light
   tier certifies to about half the margin (S8 between 0.575 and 0.6 of r3*;
   S3 0.48 r*, measured by another session), and the heavy tier certifies
   both at 0.9 (S3 1.02M dv, 97 s; S8 1.34M dv, 157 s). With a transport state
   as the edge data (S17), no setting tried certifies – light down to 0.1
   r2*, heavy up to 8.75M dv – although each part of S17 certifies alone, and
   so does the same cascade with a heat state as the source.
8. **A delayed PDE is certified only where its stability does not depend on
   the delay.** Library PDE 43 with the sign corrected (S18; the delay a 2-D
   transport state) is certified in the 2 parameter settings tried where its
   first mode is stable for every delay, and in none of the 4 tried where
   that stability depends on the delay, although the scalar delay equation
   of the same mode certifies up to 0.9 tau* at dw 0. The heavy tier at
   tau = 0.5 (6.7M dv, 49 min) does not change this.
9. **H-infinity: the default form 'Q' fails on 2-D PIEs; form 'P' is tight.**
   Form 'Q', the default of `PIETOOLS_Hinf_gain_sop`, certifies none of the
   four mixed-domain gains, nor the gain of the control C1, a single 2-D
   heat state: MOSEK stops with numerical error 2, at light, and on C1 at
   heavy. Form 'P' certifies all five at dw 0, at most 0.37% above the exact
   gain, in 12 to 70 s per case against 39 to 148 s for 'Q'.
10. **H2: only the coercive form 'oco' certifies a 2-D norm, 11–15% above
    the exact value.** 'o' fails on all three norms; on C1 the default 'c'
    fails and 'cco' is infeasible. 'oco' certifies two of three at light
    (+12.3%, +15.2%; the third, +12.2%, stops at UNKNOWN), and the heavy tier
    leaves the gap at +12.5% (S24) and +10.9% (C1).
11. **The 3-D chain does not solve on the workstation.** S13 (1-D edge ->
    2-D face -> 3-D heat) at light: 12.5M dv, m 85,307; MOSEK did not return
    within 90 min at 8 threads.

**Proposal for the owner of `executives_sopvar`** (reported, not patched):
in two or more spatial variables, make 'P' the default form of
`PIETOOLS_Hinf_gain_sop` – 'Qd' is already posed as 'Pd' there – and 'oco'
the default of `PIETOOLS_H2_norm_sop`. On these tests the change turns 5 of 5
failed gains into certified bounds within 0.37%, and 3 of 3 failed H2 norms
into bounds within 15.2% (two certified).

## Catalogue

The PIE state space lists the blocks of the fundamental state in PIE order,
auxiliary blocks included. "Margin" is the exact critical value of the
parameter named; the executives run at 0.9 and 1.1 times it.

| System | PIE state space | Coupling exercised | Baseline |
|---|---|---|---|
| S1 | R x L2[s1] | ODE in a Dirichlet BC of a 1-D heat | spectrum {-a} u {-n^2 pi^2}; margin a* = 0 |
| S2 | R x L2[s2] x L2[s1,s2] | ODE as data on an edge of a 2-D heat, Neumann in s1 (auxiliary state, K) | cascade; margin a* = 0 |
| S3 | L2[s1] x L2[s1,s2] | 1-D state as Dirichlet data on an edge (cascade) | margin r* = pi^2 |
| S4 | L2[s1] x L2[s1,s2] | dynamic (Wentzell) edge: x1 = x2(s1,0), x1_t = x1'' + c x2_s2(s1,0) + r x1 | margin r* = pi^2 + k0^2, k0 tan k0 = c |
| S5 | R x L2[s1] x L2[s1,s2] | ODE -> 1-D (Neumann data) -> 2-D edge; the ODE in a 2-D BC with coefficient (1-s2) | cascade; margin a* = 0 |
| S6 | L2[s1] x L2[s1] x L2[s1,s2] | in-domain 1-D <-> 2-D (integral over s2, multiplier), Neumann in s2 (aux, K) | margin g* = pi^2 |
| S7 | R x R x L2[s1] x L2[s2] x L2[s1,s2] | ODE <-> 2-D in-domain, Neumann on every edge (3 aux, K) | margin g* = sqrt(ab) = 1 |
| S8 | L2[s1] x L2[s2] x L2[s1,s2] | 1-D states on the edges s2 = 0 and s1 = 1 | margin r3* = pi^2 |
| S9 | R x L2[s2] x L2[s1,s2] | ODE sets the edge value, the edge flux drives the ODE (two-way, K) | margin a* = c; H-inf of 1/(s - a + c sqrt(s) coth sqrt(s)) |
| S10 | R x L2[s1,s2] | ODE as Dirichlet data next to Dirichlet edges: corner clash | must be rejected |
| S11 | L2[s1] x L2[s1,s2] | corner clash in a derivative: x2_s1(0,s2) = 0 needs x1'(0) = 0 | must be rejected |
| S12 | L2[s1] x L2[s1,s2] | s1 on [0,1] for x1 and on [0,2] for x2 | must be rejected |
| S13i | L2[s1] x L2[s1,s2] x L2[s1,s2,s3] | edge clash in 3-D | must be rejected |
| S13 | L2[s1] x L2[s1,s2] x L2[s1,s2,s3] | 1-D edge -> 2-D face -> 3-D | semantics |
| S14 | L2[s1] x L2[s1,s2] | input in a 2-D BC group with a lower-dimensional source | must be rejected (documented) |
| S15 | R x L2[s1] | input and ODE in one 1-D BC | semantics |
| S16 | L2[s1] x L2[s1] x L2[s1,s2] | L2[s1,s2] input; scalar output and L2[s1] trace output | semantics |
| S17 | L2[s1] x L2[s1,s2] | transport edge state into a state of order 0 in s1 | margin r2* = pi^2 |
| S18 | L2[s1] x L2[s1,s2] | library PDE 43, delay as a 2-D transport state, both signs | tau* = 1.034726 (sign -1); unstable for every tau (sign +1) |
| S19 | L2[s1] x L2[s1,s2] | higher-dimensional state in a 1-D BC | PIE exists; converter rejects |
| S20 | L2[s1] x L2[s1,s2] | S4 with the 2-D state declared on [s2;s1] | semantics |
| S21 | L2[s1] x L2[s1,s2] | S4 on [-1,2] x [0.5,3], edge coefficient (2+s1) | semantics |
| S22 | L2^2[s1] x L2[s1,s2] | vector edge state, matrix coupling | semantics |
| S23 | L2[s1] x L2[s1] x L2[s1,s2] | H-inf with an L2[s1,s2] disturbance (S6, g = 1) | gain 0.0093385782 (series) |
| S24 | R x L2[s2] x L2[s1,s2] | H-inf, ODE -> 2-D edge (K) | gain 1/(2a) = 0.5 |
| S25 | L2[s1] x L2[s1,s2] | H-inf with an L2[s1] disturbance on the Wentzell edge | gain 0.0692045143 (series) |
| S26 | R x L2[s1] x L2[s2] x L2[s1,s2] | library PDE 42 (Kang and Fridman 2017): ODE + heat, both delayed | margin a* = pi^2/4 - 1/2 |
| S27 | R x L2[s1] x L2[s1,s2] | x1' as edge data: x2 in H^{2,2} needs x1''' | must be rejected |
| S28 | L2[s1] x L2[s1,s2] | Robin edge condition with 1-D data | semantics |
| S29 | R x L2[s1] x L2[s2] x L2[s1,s2] | PIESIM 2-D library ex 8: integral BCs | semantics; exact solution |
| S30 | L2[s2] x L2[s1,s2] | PIESIM 2-D library ex 9: integral and periodic BCs | semantics; exact solution |
| S31 | L2[s1,s2] | one 2-D state, integral BC in s1 | semantics; exact solution |
| C1 | L2[s1,s2] | control for the executives: one 2-D heat state, Dirichlet on every edge, scalar w as a uniform source, z = int int x; no lower-dimensional or auxiliary state | gain 0.0351442537, H2 norm 0.1150629240 (series over the modes) |

## Layer 1 – parser and converter

**Method.** For each accepted system, states x satisfying the BCs are built
by hand (polynomials, random coefficients, two draws). With xf_b = F_b
D^{diff_b} x_k (F_b the integral over the variables of state k that block b
lacks; F_b = I for a main block), each PIE operator is applied block by block
through `apply_sopvar`, the kernel definition, exactly on polynomials, and
compared with the PDE written by hand: T xf + Tw wf against x (and F_b x),
A xf + Bw wf against the right-hand side (and F_b of it), Cz xf + Dzw wf
against z, K xf against 0, and `pde2pie.x.C` against xf. Without auxiliary
blocks a second check needs no hand-built state: random xf, then x = T xf must
satisfy the BCs, D^{diff} x must return xf, and A xf must equal the
right-hand side of x. This second check verifies T as an inverse of the
derivative without using `pde2pie`, so the two directions are never tested
against each other.

**Result.** 25 of 25 accepted systems pass: largest error 2.5e-14 (S13,
D^{diff} T xf), every other quantity at most 2.2e-14 (S21, BC residual). The
control C1, added after the run of record, passes as well (5.5e-15). All six
systems marked "must be rejected" are rejected; the messages:

| System | Message (abridged) |
|---|---|
| S10, S11, S13i | "The boundary conditions k are not satisfied by the states coupled to states of lower dimension: these states must satisfy the compatibility (corner) conditions ..." (`eliminate_BVP_sources:275`) |
| S12 | "The same spatial variable seems to exist on different domains ..." (`initialize:397`) |
| S14 | "The boundary conditions [3 4 5 6] involve both a state of lower dimension and an input; this is not supported ..." (`solve_PIETOOLS_BVP:304`) |
| S27 | "The order of the derivative of x_1 is higher than that to which the state is differentiable." (`convert_PDE_factor:269`) |

S11 shows that the compatibility check covers derivative conditions, not only
values. S19 is the one valid system rejected (D2).

## Layer 2 – copvar operations

On T and A of 25 PIEs (every accepted system), with random polynomial blocks
f, g on the in- and out-spaces: (A*T)f against A(Tf); (T'*A)f against
T'(Af); (T+A)f against Tf + Af; (T'*A + A'*T)f against T'(Af) + A'(Tf); and
<Tf,g> - <f,T'g>, <Af,g> - <f,A'g> with the inner product integrated exactly
over each block's variables. The adjoint is checked against the inner
product, not against another class routine. Largest error 1.4e-14 (S21, A*T
and T'*A, on [-1,2] x [0.5,3]); 25 of 25 pass, and C1 (2.5e-16).

## Exact baselines

Derived by separation of variables; each margin checked against the
rightmost eigenvalue of the collocation model (N = 24) at the margin, which
must vanish.

| Quantity | Exact value | Collocation |
|---|---|---|
| S3 r* (x1 modes r - m^2 pi^2; x2 modes -(m^2+n^2) pi^2) | pi^2 | max Re -9.3e-14 |
| S4 r*: mode sin(m pi s1), x2 = sinh(mu(1-s2)), mu tanh mu = -c, mu = i k, lam = r - m^2 pi^2 - k^2 | pi^2 + 0.860333589^2 = 10.609778285 | -1.7e-11 |
| S6 g*: modes sin(m pi s1) x {1, cos(n pi s2)}; lam = -m^2 pi^2 +- g | pi^2 | 2.5e-11 |
| S7 g*: mean mode [-a g; g -b] | 1 (a = b = 1) | 1.3e-11 |
| S8 r3* | pi^2 | -1.6e-11 |
| S9 a*: s1-mean mode, lam + c sqrt(lam) coth sqrt(lam) - a = 0 | c = 1 | 2.6e-11 |
| S9 at a = 0.5: rightmost root -k^2, k^2 - k cot k + 0.5 = 0 | -0.372600855 | -0.37260086 |
| S26 a*: mode cos(pi s1/2), b' = (a - pi^2/4) b + 0.5 b(t-1) | pi^2/4 - 0.5 | not modelled |
| S18 tau* (sign -1): lam = 0.9 - exp(-lam tau) crosses at omega = sqrt(0.19) | 1.034726 | not modelled |
| h1: S9, w into the ODE, z = X: 1/(s - 0.5 + sqrt(s) coth sqrt(s)), sup at omega = 0 | 2 | 2 |
| h3: S24: (cosh mu - 1)/(mu sinh mu) /(s+1), sup at omega = 0 | 0.5 | 0.5 |
| h2: S23, series over odd m of sqrt2 (2/(m pi)) g/((s+m^2 pi^2)^2 - g^2) | 0.0093385782 | – |
| h4: S25, series over odd m of sqrt2 (2/(m pi))/(mu^2 + mu coth mu) | 0.0692045143 | – |

The collocation model discretizes the PDE (not the PIE): BCs at boundary
nodes, boundary values eliminated, `expm` in time; it shares no code with
PIETOOLS. Its self-test: derivative error 5.7e-14, second derivative 1.6e-11,
quadrature 5.6e-17 on x^3. A sign error in its first-derivative matrix was
caught by this cross-check (S4 and S9 disagreed; the pure-heat margins,
second derivatives only, agreed) and corrected before any comparison above.

## Layer 3 – container executives

**Method.** `PIETOOLS_stability_sop` (form 'P') runs at 0.9 and 1.1 times
each exact margin: a stable PIE must be certified, an unstable one must not
be. `PIETOOLS_Hinf_gain_sop` (form 'Q', the default, and form 'P') and
`PIETOOLS_H2_norm_sop` (forms 'o' and 'oco') run against the exact gains and
H2 norms: a certified bound below the exact value is unsound; a bound more
than 5% above it, or no certificate, is a FAIL. Settings
`lpisettings_sop('light')`, MOSEK 11. While the solver returns UNKNOWN or
INFEASIBLE, the executive raises the degree dw of the negativity slack by one,
up to the limit of the tier; each entry gives the status of the last solve,
its dw, the decision variables (dv) and the time of the executive call. Run
of record: `test_mixed_domain_pde('all')`, 11/12 05:44–06:40, 125 checks.

**Soundness.** No case certified an unstable PIE (13 of 13 rejected), and
every certified gain or H2 norm lies at or above its exact value – in the
tables of this section and in every further run below.

Stability, light:

| System | Coupling | Stable case | Unstable case |
|---|---|---|---|
| S1, S2, S5 | ODE into a 1-D or 2-D BC (cascade) | a = 0.5: certified at dw 0, at most 418k dv, 22 s | a = -0.5: not certified |
| S4 | Wentzell edge (two-way) | 0.9 r*: certified at dw 0, 318k dv, 8 s | 1.1 r*: not certified |
| S6 | in-domain 1-D <-> 2-D (K) | 0.9 g*: certified at dw 0, 525k dv, 29 s | not certified |
| S7 | ODE <-> 2-D in-domain (K) | 0.9 g*: certified at dw 0, 529k dv, 33 s | not certified |
| S9 | ODE <-> 2-D edge flux (K) | 0.9 a*: certified at dw 0, 418k dv, 19 s | not certified |
| S26 | library PDE 42 | 0.9 a*: certified at dw 0, 337k dv, 17 s | not certified |
| S3 | 1-D state as edge data | 0.9 r*: **not certified** (unknown at dw 1, 508k dv, 51 s) | not certified |
| S8 | 1-D states on two edges | 0.9 r3*: **not certified** (unknown at dw 1, 1.34M dv, 139 s) | not certified |
| S17 | transport state as edge data | 0.9 r2*: **not certified** (unknown at dw 1, 664k dv, 73 s) | not certified |
| S18 | PDE 43, sign corrected | 0.9 tau*: **not certified** (infeasible at dw 1, 508k dv, 33 s) | 1.1 tau*: not certified |
| S18, sign +1 | PDE 43 as coded, tau = 0.5 | – | not certified (infeasible at dw 1) |

The four stable cases left uncertified are the three in which a
lower-dimensional state is Dirichlet data on an edge of the 2-D state, and
PDE 43. Where the 1-D state and the 2-D state drive each other through the
edge (S4, S9), the light tier certifies at dw 0.

Near the margins, light, one MATLAB process per case (scratch runners; see
*Records*):

| System | Stable case | Unstable case |
|---|---|---|
| S4 | 0.98 r*: certified at dw 0, 318k dv, 18 s | 1.02 r*: not certified (infeasible at dw 1) |
| S6 | 0.98 g*: certified at dw 0, 525k dv, 50 s | 1.02 g*: not certified (unknown at dw 1) |
| S7 | 0.98 g*: certified at dw 0, 529k dv, 56 s | 1.02 g*: not certified (infeasible at dw 1) |
| S9 | 0.98 a*: certified at dw 0, 418k dv, 35 s | 1.02 a*: not certified (unknown at dw 1) |
| S26 | 0.98 a*: certified at dw 0, 337k dv, 31 s | 1.02 a*: not certified (unknown at dw 1) |
| S1, S2, S5 | a = 0.02: certified at dw 0 (S1, S2) and dw 1 (S5) | a = -0.02: not certified |

At heavy, the tier that certifies S3 at 0.9 r*, S3 at 1.02 r* is infeasible at
dw 0 and dw 1 (2.91M dv); the 40-min cap stopped the dw 2 solve.

H-infinity gain, light (forms 'Q' and 'P'):

| Gain | Exact | Form 'Q' (default) | Form 'P' |
|---|---|---|---|
| h1: S9 (a = 0.5), w into the ODE, z = X (K) | 2 | 4.4e3, not certified | 2.0002451 (+1.2e-4), dw 0, 52 s |
| h2: S23 (g = 1), w in L2[s1,s2] (K) | 0.0093385782 | 57.6, not certified | 0.0093394469 (+9.3e-5), dw 0, 70 s |
| h3: S24 (a = 1), ODE -> 2-D edge (K) | 0.5 | 198, not certified | 0.50002137 (+4.3e-5), dw 0, 38 s |
| h4: S25, w in L2[s1] on the Wentzell edge | 0.0692045143 | 526, not certified | 0.069458567 (+3.7e-3), dw 0, 39 s |
| h5: C1, control: one 2-D heat state, no lower-dimensional or auxiliary state | 0.0351442537 | 412, not certified | 0.035172999 (+8.2e-4), dw 0, 12 s |

In form 'Q' MOSEK stops with UNKNOWN (numerical error 2) at dw 0 and at
dw 1: dual norm 3e6 to 2e9 in 73 to 148 s per mixed-domain case, 2e10 and 9e9
in 39 s on C1. The control C1 fails the same way, so the failure of the
default form is a property of 2-D PIEs, not of mixed domains. At heavy, 'Q'
still fails on C1: UNKNOWN at dw 1, gam 7.71, dual norm 6e9 at dw 0, 1.35M
dv, 190 s.

H2 norm, light (forms 'o' and 'oco'):

| H2 norm | Exact | Form 'o' | Form 'oco' |
|---|---|---|---|
| n1: S9 (a = 0.5, c = 1) | 0.8692300301 | 9.97e3, not certified | 0.97499 (+12.2%), not certified: UNKNOWN at dw 1, 2.26M dv, 343 s |
| n3: S24 (a = 1) | 0.3377221541 | 293, not certified | 0.37915 (+12.3%), certified at dw 1, 2.26M dv, 300 s |
| n5: C1 (control) | 0.1150629240 | 341, not certified | 0.13256781 (+15.2%), certified at dw 0, 188k dv, 17 s |

In form 'oco' the raise from dw 0 to dw 1 lowers the solver objective by
0.4% (n1) and 0.5% (n3); on n1 the primal and dual objectives at dw 1 agree to
2e-4, with numerical error 2. On C1 form 'o' fails and form 'oco' is 15.2%
above, so neither is specific to mixed domains. At heavy, form 'oco' on
S24 certifies 0.379978 (+12.5%) at dw 0 (1.34M dv, 265 s): the higher degrees
of the heavy tier do not close the gap. On C1 at heavy, 'oco' certifies
0.127655 (+10.9%) at dw 0 (752k dv, 75 s), and 'o' still fails (UNKNOWN at
dw 1, gam 15.8, 2.14M dv, 381 s).

The default form of `PIETOOLS_H2_norm_sop` is 'c'; on S9 and S24 (nonempty K)
it is posed as 'o'. C1 has no K, so there 'c' runs as itself, at light: UNKNOWN
at dw 0 and 1 (gam 7.7e5, dual norm 2e13, 78 s), and 'cco' is primal
infeasible at dw 0 and 1 (122 s). On C1, 'oco' is the only H2 form that
certifies. In 'c' and 'cco' the adjoint of the output map, the constant
function on L2[s1,s2], enters as a distributed input; the infeasibility of
the dual forms for distributed inputs was recorded in earlier container work
(inferred to be the same cause, not traced).

C1 values: gain by quadrature of the integral of f(t)^2, f(t) = sum over odd
m of 8/(m^2 pi^2) exp(-m^2 pi^2 t) (the impulse response is f^2, every residue
positive), equal to a double series over the modes to 10 digits; H2 norm =
sqrt of the integral of f^4. C1 entered the test after the run of record.

**Beyond light.** The four stable systems left uncertified, run at other
parameters and settings on the same export, by scratch runners that call the
executives as the test does (code and logs: see *Records* below). "Corner rows"
are the psatz corner-product rows of another session, file sets A and B of
*Tree* above.

| System | Parameter | Settings | Result |
|---|---|---|---|
| S3 | 0.5 r* | light | not certified (unknown at dw 1, 508k dv, 57 s) |
| S3 | 0.5 r* | light, corner rows A | certified at dw 0, 572k dv, 25 s |
| S3 | 0.9 r* | light, corner rows A | not certified (unknown at dw 1, 2.29M dv, 223 s) |
| S3 | 0.9 r* | light, corner rows B | not certified (unknown at dw 1, 762k dv, 70 s) |
| S3 | 0.9 r* | heavy | **certified** at dw 0, 1.02M dv, 97 s |
| S3 | 0.9 r* | veryheavy | **certified** at dw 0, 5.07M dv, 1322 s |
| S8 | 0.5, 0.55, 0.575 r3* | light | certified at dw 1, 1.34M dv, 131–155 s |
| S8 | 0.6, 0.7 r3* | light | not certified (unknown at dw 1, 1.34M dv, 132–150 s) |
| S8 | 0.5 r3* | light, corner rows A | certified at dw 0, 1.68M dv, 131 s |
| S8 | 0.9 r3* | light, corner rows A | not certified (unknown at dw 1, 5.31M dv, 1539 s) |
| S8 | 0.9 r3* | light, corner rows B | not certified (unknown at dw 1, 2.08M dv, 185 s) |
| S8 | 0.9 r3* | heavy | **certified** at dw 0, 1.34M dv, 157 s |
| S17 | 0.1, 0.25, 0.5 r2* | light | not certified (unknown at dw 1, 664k dv, 58–90 s) |
| S17 | 0.5, 0.9 r2* | light, corner rows A | not certified (unknown at dw 1, 2.99M dv, 485 and 564 s) |
| S17 | 0.9 r2* | heavy | not certified (unknown at dw 0, 1 and 2; at dw 2: 8.75M dv, m 14,188, solve 4349 s; call 5208 s) |
| S18 | tau = 0.1 | light | not certified (unknown at dw 1, 508k dv, 51 s) |
| S18 | tau = 0.25 | light | not certified (infeasible at dw 1, 508k dv, 41 s) |
| S18 | tau = 0.5 | heavy | not certified (unknown at dw 2, 6.70M dv, 2963 s) |
| S18 by the delay syntax (`subs(x1,t,t-tau)`; the converter builds the transport state) | tau = 0.5 | light | not certified (infeasible at dw 1, 664k dv, 54 s); sign +1: not certified, as required |

The light reach of S8 lies in [0.575, 0.6] r3*, by bisection. The session
*New Posopvar and Duality* measured the light reach of S3 by bisection:
0.4844 r* without corner rows, 0.6719 r* with them (its runs, not repeated
here). Set B gives programs of 33% (S3) and 39% (S8) of the decision variables
of set A, and certifies neither case at 0.9. S17 and S18 are certified at no
parameter tried. At every r2, the S17 solves end with MOSEK UNKNOWN while the
primal norm grows to 1e6 and the dual objective falls toward 0 (feasratio
0.08 at 0.25 r2*); this trace suggests a program on the boundary of
feasibility (inferred from the trace, not shown).

**Localizing S17 and S18.** The parts of S17, the S18 delay mode alone, and
S18 at other coefficients, at light:

| Subsystem | Result |
|---|---|
| S17 x1 alone: x_t = -x_s1, x(0) = 0 (nilpotent) | certified at dw 0, 280 dv, 1 s |
| S17 x2 alone: order 0 in s1, heat in s2, r2 = 0 and 0.9 pi^2 | certified: dw 0, 141k dv, 3 s; dw 1, 256k dv, 10 s |
| S17 with zero edge data (the two parts side by side), r2 = 0.9 pi^2 | certified at dw 1, 664k dv, 48 s |
| S17 with a Dirichlet 1-D heat state as the edge data, r2 = 0.9 and 1.1 pi^2 | 0.9: certified at dw 1, 664k dv, 49 s; 1.1: not certified, as required |
| y' = 0.9 y - y(t - tau), the sin(s1) mode of S18 | certified at tau = 0.1, 0.5, 0.9 tau* (dw 0, 339 dv, under 1 s); 1.1 tau*: not certified, as required |
| S18 with a0 = 0, a1 = 0.5 (every mode stable for every tau), tau = 0.5 | certified at dw 0, 318k dv, 7 s |
| S18 x1 without the delay term: x_t = x_s1s1 + 0.9 x on [0, pi] (rightmost eigenvalue -0.1) | certified at dw 1, 568 dv, 2 s |
| S18 with a1 = 0 (the delay state a driven transport), a0 = 0.9, tau = 0.5 | certified at dw 0, 318k dv, 13 s |
| S18 with a1 = 1, tau = 0.5, a0 = 0.9, 1.2, 1.5 (stable: tau* = 1.679, 1.398, 1.209) | not certified: unknown (a0 = 0.9), infeasible (1.2, 1.5) at dw 1, 508k dv, 67–79 s |

S17 fails only when the edge data comes from the transport state: each part
certifies, and so does the same cascade with a heat state as the source.
S18 fails where the stability of its first mode depends on the delay. The
mode sin(s1) obeys y' = (a0 - 1) y - a1 y(t - tau): for 0 <= a1 < 1 - a0
it is stable for every tau, and for a1 > 1 - a0 its stability depends on
tau. The executive certifies S18 in the two settings tried with a1 < 1 - a0
(a0 = 0, a1 = 0.5; a0 = 0.9, a1 = 0) and in none of the four with
a1 > 1 - a0 (a1 = 1; a0 = 0.9, 1.2, 1.5 at tau = 0.5; a0 = 1.9 at tau = 0.1,
0.25, 0.5, 0.93), although the scalar equation of that mode certifies up to
0.9 tau* at dw 0 (a0 = 1.9). Neither the slow decay (rate 0.1, certified
without the delay term) nor the transport state (certified with a1 = 0) fails
on its own.

**3-D.** S13 (1-D edge -> 2-D face -> 3-D heat) at light: 12,496,976 dv in
seven Gram blocks of 1336, m = 85,307. MOSEK did not return within the 90-min
cap. The workstation (24 cores, 64 GB) ran MOSEK at 8 threads in every run:
each log reports "TBB Warning: The number of workers is currently limited to
7". Not run on the ASU Sol cluster, which needs approval.

## Layer 4 – PIESIM

PIESIM (PDE input, BDF order 4, dt 1e-5, final time 0.1), compared at its
grid with: the exact solutions (library examples and S31, N = 24), the
collocation model (N = 28, interpolated barycentrically; PIESIM N = 16), and
`dde23` (RelTol 1e-12) on the first mode of PDE 43, from a history twice
continuously differentiable into the solution (a kink in the history travels
through the transport state and limits the spectral accuracy: 3.5e-6 with a
history only once differentiable).

| Check | Error (max over the grid / max of the reference, per state) | |
|---|---|---|
| library ex 10 (ODE, L2[x], L2[y], L2[x,y]) | 7.1e-11, 7.1e-11, 3.0e-10, 1.4e-9 | pass |
| library ex 11 (4 states, inputs in 1-D BCs) | 2.2e-14, 5.2e-12, 2.3e-10, 4.6e-9 | pass |
| S3 | 6.7e-9, 7.0e-9 | pass |
| S4 (Wentzell) | 7.7e-9, 7.7e-9 | pass |
| S5 (ODE -> 1-D -> 2-D) | 6.5e-11, 4.6e-10, 2.8e-9 | pass |
| S8 (two edges) | 6.7e-9, 6.7e-9, 8.3e-9 | pass |
| S18, sign -1 (stable) | 1.6e-10 | pass |
| S18, sign +1 (library as coded, growing) | 3.4e-7 | pass |
| S25: input w(t,s1) = sin(pi s1) sin(2t) on the Wentzell edge, output z | 7.7e-9, 7.7e-9; z 7.7e-9 | pass |
| Atotal, rightmost Re eig, N = 16: Neumann at s1 = 0, Dirichlet at s1 = 1 (control) | -12.337, exact -12.337 | pass |
| library ex 8 | NaN in every state | FAIL (D4) |
| library ex 9 | NaN | FAIL (D4) |
| S31 (one state, integral BC) | NaN at N = 16; 4.2e163 at N = 24 | FAIL (D4) |
| 2-D heat, Neumann at both ends of s1 | crash: `A_cell_2PDEstate` undefined, `PIESIM_3PI2Mat_cheb:81` | FAIL (D3) |
| S2 | crash: `A_cell_2PDEstate` undefined, `PIESIM_Poly2Mat_cheb:77` | FAIL (D3) |
| S9 | crash: `A_cell` undefined, `PIESIM_PI2Mat_opint_cheb:95` | FAIL (D3) |
| S6 | NaN (at dt 1e-4: errors 9.5e-2 and 1.05) | FAIL (D3) |
| Atotal, rightmost Re eig, N = 16: S31 | +2.77e5, exact -19.74 | FAIL (D4) |
| Atotal, rightmost Re eig, N = 16: S6 (g = 1) | +2.77e5, exact -8.87 | FAIL (D3) |

Two non-defects, measured: (i) at dt 1e-4 the passing cases had errors near
7e-7, independent of N in PIESIM (16, 24) and in the collocation (28, 36),
and dt 1e-5 lowered them 100-fold, so the time step limits the accuracy at dt
1e-4. (ii) Library ex 10 has error 2.4e-4 at N = 16 and 8.3e-9 at N = 24
(2.3e-8 at N = 32): the y-domain [0,5] holds 2.5 periods of the solution, so
N = 16 under-resolves it; the Phase 0 value at N = 8 (max |x| 65.8 against an
exact 9.05) is attributed to the same cause, not rerun at N = 8.

## Defects and minimal reproductions

Each is reported, not patched: D1 is in an example file, D2 in the
converter (Declan Jagt), D3 and D4 in PIESIM.

**D1 – library PDE 43 has the wrong sign.**
`PIETOOLS_examples/Examples_Library/2D/PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay.m`:

    diff(x1,'t')==c*diff(x1,s1,2)+a0*x1+a1*subs(x2,s2,tau)

with `a1 = 1`. The header of the file and the library table ("x_{t} = ... -
a1*x(t-tau,s1)", "stable for tau<1.0347", both citing Caliskan 2009) have the
minus sign; the paper text was not read. With x2(t,s1,s2) = x1(t-s2,s1) (transport, x2(s1,0) = x1),
the mode sin(s1) of x1 obeys a' = (a0 - c) a + a1 a(t - tau); with the plus
sign, f(lam) = lam - 0.9 - exp(-lam tau) has f(0) = -1.9 < 0 and f -> +inf,
so a positive real root exists for every tau >= 0. Measured: the PIE as coded
passes Layer 1 (it is the coded PDE), the executive returns "infeasible" at
tau = 0.5, and PIESIM matches `dde23` on the growing mode to 3.4e-7. With the minus sign
the margin is tau* = atan2(sqrt(0.19), 0.9)/sqrt(0.19) = 1.034726, the value
the library header states. Fix: `-a1*subs(x2,s2,tau)` (or `a1 = -1`); the
header already states the intended system.

**D2 – a higher-dimensional state in a lower-dimensional BC is rejected with
a corner/edge message.**

    pvar s1 s2 t
    x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
    PDE = [diff(x1,t)==diff(x1,s1,2); diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2)+x1;
           subs(x1,s1,0)==int(x2,[s1;s2],[0,1;0,1]); subs(x1,s1,1)==0;
           subs(x2,s1,0)==0; subs(x2,s1,1)==0; subs(x2,s2,0)==0; subs(x2,s2,1)==0];
    PIE = convert(PDE,'pie');
    % split_BVPs:141 "The condition BC{1} fixes the variables s1, s2
    % simultaneously (e.g. at a corner or edge); ..."

Cause: `split_BVPs` marks a variable fixed when any term of the BC is
evaluated at a boundary or integrated over the full domain, so the full
integral of x2 fixes s2 as well as s1. The PIE exists: x2 does not depend on
x1 through its BCs, so x1 = T11 x1f + (1-s1) int int T22 x2f. Proposed fix:
lift such terms as sources, the converse of `lift_BVP_sources` – solve the
BVP of x1 with the term as an input, then compose with the T-map of x2 – when
the BCs of the higher-dimensional state do not involve the lower one; at
least, name the higher-dimensional term in the error.

**D3 – PIESIM fails on PIEs with auxiliary states.** The converter turns a
direction with conditions at both ends that leave the state rank deficient
(Neumann-Neumann, periodic) into an auxiliary block (the integral of the
state over that direction) and a constraint PIE.K. PIESIM reads K only in
`rescalePIE_2D`. Minimal reproduction, one state:

    pvar s1 s2 t;  syms st sx sy real
    x = pde_var(1,[s1;s2],[0,1;0,1]);
    PDE = initialize([diff(x,t)==diff(x,s1,2)+diff(x,s2,2);
          subs(diff(x,s1),s1,0)==0; subs(diff(x,s1),s1,1)==0;
          subs(x,s2,0)==0; subs(x,s2,1)==0]);
    uex = cos(pi*sx)*sin(pi*sy)*exp(-2*pi^2*st)+sin(pi*sy)*exp(-pi^2*st);
    opts = struct('N',16,'tf',0.1,'intScheme',1,'Norder',4,'dt',1e-5,'plot','no','ploteig','no');
    PIESIM(PDE,opts,struct('ic',subs(uex,st,0)));
    % Unrecognized function or variable 'A_cell_2PDEstate' (PIESIM_3PI2Mat_cheb:81)

S2 and S9 crash the same way in other routines; S6 (auxiliary block in
L2[s1]) does not crash but returns NaN, and at dt 1e-4 it returns finite
states with relative errors 0.095 and 1.05, without a warning. Proposal for
the PIESIM owner: either map the PDE initial condition through
`PIE.pde2pie.x.C` to every block and carry the auxiliary blocks, or refuse a
PIE with nonempty K with a message.

**D4 – PIESIM returns NaN on integral BCs in 2-D.** S31 (one state,
`int(x,s1,[0,1])==0`, `x_s1(1,s2) = 0`, Dirichlet in s2; exact
cos(pi s1) sin(pi s2) exp(-2 pi^2 t)) returns NaN in every grid value; so do
library ex 8 and ex 9. Their PIEs pass Layer 1 (S29, S30, S31: errors at most
2.5e-15), so the conversion is correct and the defect is in the
discretization or the time stepping. Recorded since 11/02 in
`testing_codes/private/baselines/2026-11-02_opvar_phase0/F1/pass*/index.csv`
(rows 2,8 and 2,9: max_abs_final 2.0e41 and 5.6e19, err_exact NaN).

Diagnosis: the blow-up starts in the discrete generator
`Dop.Atotal = Tcheb_inv*Acheb` of `PIESIM_discretize_ops_2D`, whose
rightmost eigenvalues are large and positive against a negative exact
spectrum – S31 +1.71e4 at N = 8 and +2.77e5 at N = 16 (exact -19.74);
library ex 9 +427 and +6.94e3; ex 8 +9.49e3 and +1.54e5 – while with ordinary
BCs the same steps return the exact rightmost eigenvalue (-12.337). Doubling
N multiplies each spurious eigenvalue by 16.2, the N^4 scaling of the
fundamental state of a state second order in both directions. The
auxiliary-state PIE S6 has the same spurious eigenvalues as S31 to five
digits, so one discretization path for integral-type terms in T appears to
underlie D4 and the non-crashing part of D3 (inferred from the spectra, not
traced in the code). cond(Tcheb) is 3e10 at N = 16, against 3e8 for the
control. Related, measured once, 1-D: `x_t = x_ss`, `int x ds = 0`,
`x_s(1) = 0` stops with "Unrecognized function or variable 'ns'"
(`PIESIM_transform_to_solution:161`).

## Literature

- Library PDE 42: Kang and Fridman, "Boundary control of delayed ODE–heat
  cascade under actuator saturation", Automatica 83 (2017) 252–261. The
  margin a* used here is derived from the mode-0 delay equation, not taken
  from the paper.
- Library PDE 43: Caliskan and Ozbay, "Stability analysis of the heat
  equation with time-delayed feedback", IFAC Proceedings Volumes 42(6)
  (2009) 220–224; the library states its result as "stable for tau <
  1.0347"; the mode-1 characteristic equation gives 1.034726. The paper
  text was not read.
- Wentzell edge (S4, S25): Favini, Goldstein, Goldstein and Romanelli, "The
  heat equation with generalized Wentzell boundary condition", J. Evol.
  Equ. 2 (2002) 1–19 – the model and its semigroup; the eigenvalue equation
  k tan k = c and the margins here are derived, not taken from it.
- The R^n x L2[x] x L2[y] x L2[x,y] PIE framework: Jagt and Peet,
  arXiv:2109.06423 (ACC 2022); delayed 1-D PDEs as coupled 1-D/2-D PIEs:
  Jagt and Peet, arXiv:2211.05326 (IEEE L-CSS 2023). Whether these papers
  report stability margins for PDEs 42 and 43 was not checked.
- No published stability margin or H-infinity gain was found for the other
  coupled mixed-domain systems (S2–S9, S17, S23–S25); their exact values
  are derived here and confirmed by the collocation model.

## How to run

    addpath(fullfile(<PIETOOLS>,'testing_codes','private','converters','tests','mixed_domain'))
    R = test_mixed_domain_pde;                 % convert, ops, exact: about 4 min
    R = test_mixed_domain_pde('piesim');       % about 6 min
    R = test_mixed_domain_pde('exec');         % about 50 min (MOSEK, 8 threads)
    R = test_mixed_domain_pde('all');          % 55 min at 6c6f60dd

Statuses: `pass`; `rejected` (an ill-posed system rejected, as required);
`KNOWN` (S19, the limitation D2); `FAIL`. At 6c6f60dd the run of record of
'all' gave 125 checks: 97 pass, 6 rejected, 1 KNOWN, 21 FAIL. The 21 FAIL
lines are D3 and D4 in 'piesim' (9 lines) and, in 'exec', the 4 stable cases
left uncertified, the 4 gains in form 'Q' and the 4 H2 norms (12 lines).
The control C1 adds 8 checks (convert, ops, 2 exact, 4 exec), run with

    R = test_mixed_domain_pde({'convert','ops','exec'},{'C1_heat2d_scalar_w'});  % 2 min
    R = test_mixed_domain_pde('exact');

5 pass; 'Q', 'o' and 'oco' fail (3 lines). With C1, 'all' gives 133 checks:
102 pass, 6 rejected, 1 KNOWN, 24 FAIL (the run of record plus these 8;
'all' was not rerun as one call). The final files, run with
`{'convert','ops','exact'}` on the export: 73 checks, 66 pass, 6 rejected,
1 KNOWN, 0 FAIL, 141 s.

## Records

`testing_codes/private/baselines/2026-11-12_mixed_domain/`: `results.csv`,
one row per executive call of the runs beyond the run of record (73 rows:
68 returned, 2 stopped at a time cap, 1 MATLAB crash and 2 runner errors,
the last three rerun); `test_runs.csv`, the check lines of the run of
record (N7) and of the C1 run of the test (N13a); `private/`, the scratch
runners and queue scripts as they ran; `README.md`, the trees, the machine
and one line per run. The MATLAB logs stay out of git.
