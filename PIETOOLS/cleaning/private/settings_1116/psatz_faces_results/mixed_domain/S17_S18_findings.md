# S17 and S18: what blocks the certificate (2026-11-16)

Read with `README.md` (the runs), `m1_m3_summary.md`, `m2_chk_summary.md` (every row) and
`S17_obstruction.md` (the argument). Every verdict is the executive's at light, dw pinned,
MOSEK unless the row says SeDuMi; "certified" is the executive's flag (optimal or inaccurate
exit through the eta gate), which the CHK and M7 runs test further.

## Verdicts

| question | S17 (r2 = 0.5 r2*) | S18 (tau = 0.5 tau*) |
|---|---|---|
| term set (M1: none, box, facesS, faces0, corners; N-D default storage) | dw 0 infeasible for every set; dw 1 UNKNOWN (numerr 2) for every set with a term | infeasible for every set at dw 0 and 1 |
| solver pairing (M2: SeDuMi, facesS and box, dw 1) | 'inaccurate' exits counted as certified, feasibility ratios -0.036 and 0.048 | infeasible |
| every basis operator, kernels split or separable per direction, light degrees (M3) | UNKNOWN for all four structures, dw 0 and 1 | infeasible for all four at dw 0; at dw 1 infeasible for three, UNKNOWN with both split |
| every basis operator at raised degrees, mult 1, int 2 per variable, joint 4 (M6, M8) | not run: the argument below excludes a certificate | **certified**: separable (M8) 767,305 variables, 64 s; s2 split (M6) 790,884, 90 s, eta 1.7e-9; dw 1 3,045,924, 740 s; both split, dw 1, 3,125,056, 789 s |
| the two-operator default storage at the raised degrees (M9) | – | infeasible (109,921 variables): the mixed operators are needed, not the degree alone |
| soundness at the margin (M8, s2 split, raised) | – | 1.1 tau* (unstable) infeasible, 146 s; 0.9 tau* certified, 138 s |
| the returned operators on the state families of the argument (CHK, CHK2) | the SeDuMi operator has no x2 block beyond eppos = 1e-4 (f independent of the layer width K to four digits); the violation it leaves is of size eppos K, beyond polynomial reach | – |
| two-sided bound ep I <= P <= I (M7) | infeasible with a MOSEK certificate at ep = 0.5 and 0.1, UNKNOWN at 0.01; S3 (control) feasible at 0.1 and 0.01 | – |

## S17: a structural obstruction, not a degree, term-set or solver limit

`S17_obstruction.md` proves that no coercive PI operator – multiplier plus polynomial
integral kernels, the PIETOOLS class at any degree – has a nonpositive Lyapunov derivative
along S17, for any r2. Two state families carry the proof: a boundary layer in s2 fed by
the edge datum shows that the trace term at s2 = 0 must be cancelled by a cross block local
in s1 (q(s1,0) <= -r(s1,0)); an oscillation in s1 of the fundamental state shows that such a
cross block cannot be dominated when x2 carries no s1-diffusion. The record matches the
argument at every point:

- every solve exits UNKNOWN or infeasible (M1, M3; the 11/12 stress test to 8.75M
  variables), never with a MOSEK certificate;
- with the scale of P fixed by ep I <= P <= I, MOSEK returns an infeasibility certificate
  at ep = 0.5 and 0.1 and loses it at 0.01, where a violating state needs a layer of width
  about 0.01; the executive's margin 1e-4 is further into that regime, which is why its exits
  are UNKNOWN. S3, whose target diffuses in s1, is feasible at the same ep (control);
- the two SeDuMi exits that the executive counts as certified (M2) carry no x2 block beyond
  the margin, so their derivative operator violates positivity only on states with a layer
  of width about eppos, which neither the eta gate nor a polynomial family resolves – a
  flags-only certificate of the kind recorded for PDE 24 'Q';
- the heat-source variant D4 (x1 heat, x2 of order 0 in s1), which the 11/12 record lists as
  certified at 664k variables, is infeasible on today's tree in every configuration tried
  (the tier's 'auto' set at 415k, the corner rows at offset 0 at 747k and at offset 1 at
  235k; stock and N-D default storages; MOSEK and SeDuMi) and Farkas-infeasible under the
  two-sided bound, as the argument requires; the 11/12 program is not reproducible (the
  psatz defaults and the lift floor changed on 11/13–11/15) and its certificate was the flag.

So the S17 row of the stress test is not conservatism of a setting: the class has no
certificate. A Lyapunov functional for a Dirichlet edge coupling from a transport state
would need a term that is not a quadratic form in the PDE state (the argument uses that
x1_t = -x1_{s1} pairs the cross block with the un-smoothed fundamental state): a different
Lyapunov form, not a larger program.

## S18: the storage – its basis operators and its degree

| storage of P (light tier otherwise) | degrees | S18 at 0.5 tau* |
|---|---|---|
| DEG-omitted default: the pure multiplier and the pure full integral of the 2-D space | light (mult 0, int 1, joint 2) | infeasible, every term set (M1, 92k–318k variables) |
| every basis operator, separable or split per direction | light | infeasible (M3, 37k–264k) |
| DEG-omitted default | raised (mult 1, int 2, joint 4) | infeasible (M9, 110k) |
| every basis operator, separable | raised | certified (M8, 767k, 64 s) |
| every basis operator, s2 split | raised | certified (M6, 791k, 90 s); 0.9 tau* certified; 1.1 tau* infeasible |
| stock 2-D storage, heavy tables (11/12 stress test) | heavy | UNKNOWN (6.7M, 2,963 s) |

The certificate needs the mixed basis operators of the 2-D space – a multiplier in one
direction and an integral in the other, which the DEG-omitted default drops – at kernel
degree 2 per variable, joint cap 4, and a multiplier of degree 1. The lower and upper
split of the delay direction, expected from the |s2 - t2| kernel of the complete-type
delay functional, is not needed up to 0.9 tau* (not tested beyond). No tier builds this
storage: POSCOPVAR with DEG omitted keeps two operators per 2-D space and makes every
direction a 2-D space owns separable (poscopvar.m, block of 11/12/2026), and the 2-D
entries of lpi.pos are int 1, joint 2 at light, heavy and veryheavy; the stock 2-D storage
carries more operator groups, and its heavy run of 11/12 at 6.7M variables returned UNKNOWN
where the 767k-variable program certifies. The term set is irrelevant to S18 (M1).

Recommendation, reported and not patched (tier files: Repo-cleanup session; the
DEG-omitted block of poscopvar: MMP, 11/12): a per-space storage entry that includes the
mixed basis operators of a 2-D space and raises its kernel degree to 2, as a tier above
light or as an option of the delay-type PIEs.

## For the psatz-faces question

On the mixed-domain systems the negatives were never the term set: S3 and S8 were the
storage (N-D default), S18 is the storage (mixed basis operators at degree 2), S17 is the
Lyapunov class. The 2N + 1 faces, the Putinar set and the corners give the same verdict on
every one of them.

## Costs (workstation, MOSEK at 8 threads unless SeDuMi)

| batch | solves | wall |
|---|---|---|
| M1 | 20 | 6 min |
| M2 | 4 (SeDuMi) | 4 min |
| M3 | 16 | 3 min |
| CHK, CHK2 | 11 executive solves + families | 20 min |
| M6 | 3 (0.8M to 3.1M variables) | 27 min |
| M7 | 9 | 2 min |
| M8 | 3 | 6 min |
| M9 | 2 | 1 min |

## Addendum, 11/16 04:00 (supersedes the two "do not help" claims above)

The direct PIE2PDE form – Q free with the shape of T, T'Q = Q'T, T'Q >= eps T'T and
-(A'Q + Q'A) >= 0 as reader-sized inequalities, no storage – certifies S17 at 0.5 r2* (Q of
degree 1, 435,350 variables, 22 s; also at degree 2 with dw 0 and 1) and S18 at 0.5 tau*
(228k variables, 11 s), both infeasible at 1.1 of the margin and both standing at eps = 0.1
(`parameterization.md`, `m2_chk_summary.md`). The obstruction holds for the 'P' form only;
the paragraph "The other forms do not help in 2-D" is withdrawn for the Q form. S17 is
certifiable in the default notion of the toolbox once the Q form is posed without a storage.
