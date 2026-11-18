# Are the single faces needed or useful? Summary of the autonomous run of 2026-11-15

Read with `RESULTS.md` (every table), `DECISIONS.md` (D1-D27) and the `c*_summary.txt` files (one per job) in this folder; the plan is `sopvar/psatz_faces_test_plan.md`, the evidence note `sopvar/psatz_terms_evidence.md`. Nothing in the repository was changed; the instrument lives in scratch shadows (session af90e0a1, `scratchpad/psatz_faces/`). Every verdict below is certified (MOSEK optimal with rel_b <= 1e-6 and PSD blocks, or a verified Farkas ray) unless it is marked as the executive's own flag.

## The hypotheses and their verdicts

- **H0 (the maintainer's): the box-set negatives of the record were degree-structure artefacts.** Confirmed on every pinned-P plant. With the lift sized from the target and the box factors lowered in their own direction only, the Putinar set {1, g_1, g_2} gives the verdicts of the faces rung by rung on the DD x DD plant, the local-coefficient plants L1 and L2 and the ODE-coupled plant (Phase 1, 1b), at 68-75% of the faces' variables; the record's separation (box 0.99 against faces 0.9999 at one MOSEK status) does not reproduce. The only set that fails is the stock product g_1 g_2 alone, with Farkas certificates.
- **H1: at per-direction offsets the box set certifies at the same weight degree as the faces.** Holds at w >= 2 everywhere in two variables: the heat benchmark DD x DN at rtol 1e-4 (every set 12.3367312), the io2 gain (every set within 1e-4 relative), the library rows, S3 and S8, and the mixed-domain gain S24 (every set within 7e-5 of the exact 0.5). Fails at w = 1 and in 3-D at the bench caps, where the faces at offset 0 certify and every box configuration is Farkas-infeasible (next item).
- **H2: the per-direction offset for the faces beats the scalar rule.** No. At w >= 2 the scalar-rule faces (today's 'auto') are the smallest certifying configuration in every executive case (5-17% below the box set, 32-55% below the per-direction faces) with the same verdicts; at w = 1 neither offset rule for the faces is right: the faces at offset 0 are.
- **H3: on the executive path the term set does not change the verdict once the storage is the N-D default.** Confirmed on S3 and S8 at 0.9 r* and up to 0.9984 r*: every set with a box or face term certifies with the N-D default storage, none with the stock 2-D storage.

## The one effect that survives: the odd-degree half step

With the plain block at weight w, the face terms at the full degree reach total degree 2w + 1; the box factors reach 2w at the per-direction offset and 2w + 2 at offset 0. The certified reach orders strictly, with a Farkas ray between each pair:

| plain block at w = 1 | per-direction box (2w) | offset-0 box (2w + 2) | offset-0 faces (2w + 1) | any set at w = 2 |
|---|---|---|---|---|
| 1-D heat DD | 0.99343 (328 var.) | 0.99915 (523) | 0.99993 (783) | 0.99999 (720) |
| 2-D heat DD x DN | 0.99411 (61,449) | 0.99931 (122,889) | 0.99992 (204,809) | 0.99998 (339,209 and up) |
| 2-D heat DD x DD | 0.99301 (61,449) | 0.99766 (122,889) | 0.99961 (204,809) | |
| 2-D pinned DD x DD (rungs) | 0.65 (31,104) | 0.65 (62,208) | >= 0.85 (103,680) | 0.9999 (187,920 and up) |
| PDE 48 estimator (executive) | never certified at dw 0-3 | 0.00628 at dw 2 (323,367) | 0.00451 at dw 1 (110,135) | 0.00618 facesS at dw 2 (195,815) |
| 3-D heat, bench caps, 0.985 kappa* | infeasible (583,019) | infeasible (1,550,619) | certified (2,713,563) | box at the next cap: 29.6M, not settled on Sol (10 MOSEK iterations in 9.5 h at 64 threads, 162 GB; c5c_sol_summary.txt) |

So the faces are not needed in the asymptotic sense (T2 of the evidence note) and not needed at w >= 2 in two variables, where every set ties; they are needed within budget at the lowest weight and in three variables at the bench caps, and they are the cheapest certificate of a given reach in 2-D whenever the half step suffices (the heat benchmark at 0.9999: 204,809 against 358,665 for the box set at w = 2; the estimator: 110,135 against 195,815). In one variable the pair at the next even degree is the better buy (720 against 783), as Markov-Lukacs predicts.

## Complexity, measured at equal verdict

| regime | smallest certifying configuration | box set relative | per-direction faces relative | corners relative |
|---|---|---|---|---|
| 2-D, w >= 2 (pinned, heat, gains, S3, S8, library) | scalar-rule faces | 0% to +17% | +32% to +55% | +70% to +100% |
| 2-D, w = 1 (heat, pinned, estimator) | offset-0 faces | does not certify, or 1.75x at w = 2 | does not certify | +80%, same verdict |
| 3-D, bench caps | offset-0 faces (2.71M) | infeasible; 29.6M at the next cap | infeasible | not run |

The corners never changed a verdict in this run (the S3 gain of the record was the storage), the 2^N product set never beat the box set (same verdicts, +5% to +16%, and 71 GB in 3-D), and the 4^N set was not run because no case asked for a term the faces or the box set lacked.

## Recommendation for the settings proposal

1. Default term set in N >= 2: the faces, at offset 0 when the weight is 1 and under the scalar rule otherwise, which is what 'auto' does today; it was the best or tied configuration in every case of this run. In one variable the Markov-Lukacs pair (unchanged).
2. The per-direction offset rule belongs to the box factors (it is their Markov-Lukacs-consistent form and ties the faces at w >= 2 in 2-D); for the face rows it is a loss at every weight. The proposal's default o_{r,k} = 1 where r_k is nonzero should therefore apply per factor type, with face rows at offset 0 at the lowest weight.
3. Keep the named sets `none`, `box`, `faces`; keep `corners` opt-in; drop `full` (never run, no case needs it) and `product` as a named set (never better than `box`).
4. The box set is sound and equivalent at w >= 2 in 2-D and 5-17% larger than the scalar-rule faces; it has N + 1 blocks against 2N + 1, which did not translate into fewer variables or faster solves here. It is not a default anywhere on this evidence.
5. Coupled and mixed domains: the term count stays 2N + 1 by the registry; what fixed S3 and S8 was the storage (the N-D default), not any term set (H3).

## What this run did not settle

- The 3-D box set at the cap where the faces certify (29.6M variables) ran on Sol on 11/15-16 (job 1621028) and did not finish: 10 interior-point iterations in 9.5 h at 64 threads and 162 GB, residuals converging, no verdict (c5c_sol_summary.txt). Any 3-D executive comparison also needs Sol (the light slack is 1728 x 7 per term).
- The PDE 48 estimator above dw 1 is solver-limited (certified bounds worsen with nested enlargements), so the true optimum of its programs is unknown; the certified 0.00451 of faces0 at dw 1 stands.
- Every pinned-P plant was a heat equation on the unit square; wave or transport plants were not run.
- The 2-D caps of the stock settings_2d (7/16, 9/16) were held fixed in the stock runs and absent in the N-D-default runs; their interaction with the term set is untested.

## Addendum, 11/16: the two mixed-domain systems no setting certified

Item 5 above holds on S17 and S18 as well (`mixed_domain/S17_S18_findings.md`): the term set changes no verdict on either. S17 (transport state as Dirichlet edge data of a target of order 0 in s1) admits no coercive PI Lyapunov operator at any degree – proven in `mixed_domain/S17_obstruction.md` and corroborated by MOSEK infeasibility certificates once the scale of P is fixed – so its row in the stress test is not conservatism; the executive's certified flags on it (SeDuMi) and on its heat-source variant D4 (11/12 record) are artefacts of a margin-scaled violation. S18 (the delay as a 2-D transport state) certifies to 0.9 tau* once the storage carries the mixed basis operators of the 2-D space at kernel degree 2 (767k variables, 64 s at 0.5 tau*), which no tier builds and the 6.7M-variable heavy run of 11/12 did not have.

Addendum, 11/16 04:00: S17 certifies in the direct PIE2PDE form (Q free, T'Q symmetric and >= eps T'T, no storage), as does S18 at a third of the 'P' form's size; the obstruction of `mixed_domain/S17_obstruction.md` is a statement about the 'P' form. The storage, not the slack, is the parameterization gap: the PDE norm T'T itself lies outside every separable storage the tiers build (`mixed_domain/parameterization.md`).
