# Test suites: what to run after a change

The full registry is ~35 cases and hours with the 2-D and scaling rows. These eight lists are
cut by **the question each one answers**, so a change is checked by the lists whose question it
can break. Lists overlap on purpose — a case belongs to every question it answers.

```matlab
bl_check('smoke')             % run one list against the banked expectations
bl_check('twod','bank')       % run it and record expectations for cases not yet banked
```

**H2 is dropped from the testing regime (2026-09-26, maintainer decision).** That's seven cases:
four H2-norm, the H2 estimator and controller, and 2-D H2. They're commented out of
`bl_cases.m` and removed from every list below. Their rows stay in `bl_expect.tsv` and
`results/2026-09-24`, so uncommenting restores them. The drop took out the suite's only
trivial-point sentinels (`h2cco_rd1`, `h2_2dc_rd`: rel_b = 1 with a clean status). The times
below were re-summed from the measured per-case times, not re-run.

Cost is set by what gets *solved*, not by how many cases there are: every 1-D case solves in
4–115 ms, so a 1-D list costs ~45 s of MATLAB startup plus 1–3 s of PDE→PIE conversion per case.

## The lists

| list | question | cases | case time | with startup |
|---|---|---|---|---|
| `smoke` | Is the pipeline alive? | 4 | **13 s** measured | ~1 min |
| `structure` | Did a core change alter any assembled program? | 16 | ~39 s | ~1.5 min |
| `executives` | Does every executive family still run? | 7 | ~19 s | ~1 min |
| `known_fail` | Did a change fix, or silently move, a known failure? | 2 | ~3 s | ~1 min |
| `nonlinear` | Is the PIESOS path intact? | 3 | ~4 min | ~5 min |
| `twod` | 2-D assembly and the psatz generators | 4 | ~4 min | ~5 min |
| `objective2d` | The expensive 2-D objective executives | 3 | ~26 min | ~27 min |
| `scaling` | Solver comparison and size questions only | 4 | ~2.7 min | ~3.5 min |

"Case time" is the sum of per-case wall times. For `smoke` it's a measured list run. For
`structure` it's the sum of that list's measured per-case times
(`results/2026-09-24/check_*.tsv`). The rest are summed from the baseline. "With startup"
adds ~45 s of MATLAB start and path setup. **`objective2d` exceeds the 20-minute limit for runs
on the shared workstation while its owner is online; schedule it for after 4:30 AM.**

**`smoke`** — one case per structural shape the SDP extraction treats differently: free
variables present (`stab_rd1`, K.f=43), none at all (`stabpde_rd1`, K.f=0), an objective with
`lpi_ineq` so `sossolve` adds blocks *during* the solve (`hinf_rd1`), and two positive operators
with an O(1) right-hand side (`wellposed_rd1`, ‖b‖=17.7, where every other feasibility case sits
at `eppos2` ≈ 1e-6).

**`structure`** — all 16 1-D executive cases, compared *exactly* on m, K.f, block list and
nnz(At). None of those depends on solver behaviour, so any difference means the assembled
program changed. This is the list that cleared the `lpi_eq` fix: 22/22 identical, including
nnz(At), when it still held the six 1-D H2 cases.

**`executives`** — one case per executive family (Q-form and direct stability, H∞ primal and
dual, estimator, controller, well-posedness). Only checks that each builds and solves.

**`known_fail`** — the 1-D cases whose *current* behaviour is a documented failure, so the
expected result *is* the failure: `hinfduco_rd1` (γ = 8251.9) and `ctrl_rd1` (rel_b ~1e-3,
‖x‖ ~ 7e7). A change that repairs one of these is news. So is one that moves it without
repairing it — and a pass/fail suite that reports only passing cases can see neither.

**`nonlinear`** — the polyopvar/PIESOS chain is a separate code path from the LPI executives. A
one-character typo in `innerprod.m` broke all of it and nothing in the linear suite could see
that. Three forms of the Fisher local-stability test at R = 4.

**`twod`** — 2-D stability with psatz off and with the four linear generators
(`eq_use_psatz = [3;4;5;6]`). With psatz off Mosek fails (rel_b 2.08) where SeDuMi certifies, so
those two rows are sentinels for the solver as much as for the formulation. With the generators
Mosek certifies at 1e-08, so those two rows are the regression test for `poslpivar_2d`
(`0cf1add1`).

**`objective2d`** — 283–574 s *each*. Only for changes to the 2-D objective code. `hinf2_rd` is
the one case in the suite with an external ground truth (closed-form gain 0.1711; the executive
returns 1.272).

**`scaling`** — answers questions about solvers and size, never about correctness; `structure`
already covers the same executive at n=1. The n=24/32 rungs are in `bl_big.m`, since the Mosek
arm can't solve them at all.

## What did you change → what to run

This follows the dependency chain in CLAUDE.md: a change low in the chain has to be tested
through everything above it.

| you changed | run |
|---|---|
| `pvar`/`polynomial`, `dpvar` | `smoke`, `structure`, `nonlinear` |
| `opvar`, `dopvar` (1-D operator classes) | `smoke`, `structure` |
| `opvar2d`, `dopvar2d` | `smoke`, `twod` |
| `lpi_programming` (`poslpivar`, `lpi_eq`, `lpivar`, `lpiprogram`) | `smoke`, `structure`, `twod` |
| `poslpivar_2d`, 2-D settings | `twod` (+ `objective2d` if the change touches objective executives) |
| one executive | `executives`, plus the class's rows in `structure` |
| `polyopvar`, PIESOS | `smoke`, `nonlinear` |
| `sossolve`, solver plumbing | `smoke`, `structure`, `known_fail` |
| anything claimed to fix a known failure | `known_fail` |
| a performance or scaling change | `scaling` (after `structure` passes) |

Before committing a change to a core data structure, the minimum is `smoke` + `structure` —
under 90 s together.

## How verdicts are decided

| compared | strictness | why |
|---|---|---|
| m, K.f, block list, nnz(At) | **exact** | properties of the assembled program; any change means the program changed |
| rel_b | within **10×** of banked, and on the same side of 1e-4 | Mosek varies ~1e-6 relative run to run, while real changes move rel_b by orders of magnitude |

| verdict | meaning |
|---|---|
| `PASS` | structure exact, rel_b in band (for a known failure: still fails as documented) |
| `FAIL` | structure changed, or a passing case stopped passing |
| `FIXED` | a known failure now passes — find out why before celebrating |
| `MOVED` | a known failure changed but still fails |
| `ERR` | did not build or solve |
| `NEW` | no banked expectation yet; run with `'bank'` to record one |

Expectations live in **`bl_expect.tsv` at the package root**, banked from the measured
baseline; `results/2026-09-24/bl_expect.tsv` is the frozen snapshot of the same file, so don't
edit that one. `'bank'` **appends to the tracked root file**, so commit or revert it
afterwards. It only adds cases that have no expectation — it never overwrites one, so a
regression cannot be silently absorbed by re-banking.

**Not yet banked**: `stab2_rd_psz` and `stab2_rd_psz90` (the psatz rows in `twod`). The first
`bl_check('twod','bank')` records them.

## Sol lists (added 2026-09-27)

Three lists in `bl_suites` with check `cuadmm`. They are run by `bl_regime` and
`sol/sol_harness.slurm` on dumps copied from the workstation, never by `bl_check`, which refuses
them: it would build `scale_stab_n24/n32` through the executive, a 45 / 142 GB Mosek solve. Sol has
**no Mosek**, so every case below is chosen for how its correct answer is known. The choice comes
from a read-only census of every registered case (09-27), with each claim then re-checked against
its cited file and line. "Mosek-F" means certified by Mosek under the harness's repaired-PSD test.
**PDE stability by theory does not make the LPI feasible**: the LPI is only sufficient, so
feasibility is always Mosek's certificate, except where noted.

### `sol_compare`: desktop vs Sol, clean cases with a known answer

| id | m | how the answer is known | why it is in |
|---|---|---|---|
| `stab_tr1` | 96 | Mosek-F | cheapest certified case (1,973 it, 1.5 s): the smoke job and instrument check |
| `stab_rd1_hv` | 135 | Mosek-F | anchor of the scale ladder: 14,092 it to F at 1e-7; clean 4,023 it at 1e-6 |
| `stab_rd1_tight` | 135 | Mosek-F, 0.95 λ\* | near the boundary: 2.7× the iterations of 0.5 λ\* at the same m |
| `stab_wave1` | 520 | Mosek-F | largest certified 1-D feasibility case. "wave" is a misnomer: it is a damped 2-state diffusion (`private/bl_b_stab1.m:17-19`) |
| `hinf_rd1_hv` | 267 | Mosek bracket [0.1824014, 0.1824816] | H∞ primal anchor; cuADMM bisection 1.0002× Mosek; also has a held-out infeasible control at γ_I |
| `hinf_rd1` | 242 | Mosek bracket [0.1825453, 0.1826255] | the acceptance case: exercises a real lean error that the endgame caught |
| `hinfdu_rd1` | 244 | Mosek bracket [0.1824704, 0.1825555] | H∞ dual; 1.0016× Mosek at 60k caps, 1.7% loose at 20k |
| `sent101_f1p010000` | 135 | infeasible by theorem (1.01 λ\*); Mosek Farkas | negative control; cuADMM correctly not-F. Registered in `bl_cases` 09-27 |
| `scale_stab_n08`, `_n16` | 8,640 / 34,560 | by construction + Mosek-F | the top of the ladder where Mosek still solves (n16: 118 s) |

Analytic L2 gain of the 1-D plant: 0.1824790804. The census recomputed it (odd-mode series); it is
not in any file yet. Both the primal and the dual bracket contain it.

### `sol_large`: size and capability

| id | m | answer | note |
|---|---|---|---|
| `scale_stab_n24`, `_n32` | 77,760 / 138,240 | by construction | beyond Mosek (45 / 142 GB Schur). Built by `bl_big` without solving |
| `scale_sent_f1p01_n16/24/32` | 34,560 … 138,240 | infeasible by theory | built by `bl_big(ns,1.01)`: the negative control at size |
| `scale_hinf_n08` | 17,025 | Mosek bracket [0.1801925, 0.1825249] | largest objective case with a bracket; single-point F at 1.2 γ_F was marginal under the 09-26 rule |
| `stab2_rd_psz` | 4,492 | Mosek-F | 2-D; cuADMM not certified at 10k it (29.5 ms/it). The capability question for 2-D |
| `Th_n3` | 80,550 | **none** | 2-D, 26.1M nnz, blocks [90 2856]: a size point only. Converted from the ladder's dump; its kept 1e-6 point is `regime_0927/keep/` |

**`scale_stab` does not test convergence at scale.** It is n decoupled copies. cuADMM takes 519
iterations to 1e-4 at every n = 1..32, and 4,023 to 1e-6 at n = 1..16, with identical residuals. It
tests memory and cost per iteration at size, and soundness where the answer is known by
construction. Whether Th_n3's three states are coupled has not been checked.

**The coupled ladder `scale_rot` (added 09-27) fills that gap.** It is built by
`bl_big(ns,frac,true)`: `x_t = x_ss + A x` with `A = Q diag(μ) Qᵀ`, where Q is dense orthogonal
(from `rng(n)`) and μ runs from `frac·π²` down to `−frac·π²`. Every state is coupled, so the SDP is
not block-separable. The answer follows from an argument, not a measurement:
- a constant orthogonal change of state variables maps PIETOOLS' Gram basis `kron(I_n, Z(s))` to
  itself, so it preserves LPI feasibility;
- the decoupled LPI is feasible exactly when each scalar one is;
- the scalar LPI is monotone in λ.

So `frac = 0.5` is feasible and `frac = 1.01` (one unstable mode) is infeasible by theory. Checked
where Mosek can: Mosek-F at n = 2, 4, 8, and block `S2` checks n = 16. Measured to 1e-6, cuADMM takes
**3,937 / 3,261 / 6,519** iterations at n = 2 / 4 / 8, against 4,023 at every n for `scale_stab`,
so convergence now depends on size. At n = 4 it has 1.56× the nonzeros of `scale_stab_n04`, with
the same m and blocks. `scale_rot_f0p50_n04` is in the Sol-0 smoke block `X0`: 10,894 iterations,
52.2 s, F on the workstation.

**S2 result (09-27):** Mosek-F at n = 16 too. cuADMM certifies n = 8 (21,984 iterations,
10.3 ms/it), but slows down sharply beyond that:
- n = 16: about 110 ms/it, unfinished in 24 min;
- n = 24: about 380 ms/it;
- n = 32: no iteration printed in an hour.

The cause is measured: the Cholesky factor of AAᵀ, which cuADMM builds once and solves against
every iteration, is about 9% dense for the coupled ladder (54M / 260M / 812M nonzeros at
n = 16 / 24 / 32), against 0.1–0.6% decoupled. Nested dissection is worse than AMD. So cuADMM's
cost is set by that fill, not by m. See `results/2026-09-27/REGIME_REPORT.md` §4.

### `problematic`: kept apart for later

| id | kind | problem | question to investigate |
|---|---|---|---|
| `hinfco_rd1` | accuracy limit | never certifies an upper end, even at 19× Mosek's γ_F; plateaus near 7.5e-6 | why does a 176-row coercive program plateau when the non-coercive sibling certifies (K.f 44 → 1)? |
| `hinfduco_rd1` | formulation suspected | certified infeasible up to 76.8, about 421× the analytic gain; Mosek numerr 2 | is `PIETOOLS_Hinf_gain_dual_coercive` sound on a plant with known gain 0.18248? |
| `ctrl_rd1` | formulation suspected | certified [128, 1266], at least 700× what K = 0 achieves (inferred); Mosek and SeDuMi upper ends differ 3× | does the `Hinf_control` LPI admit K = 0 at γ = 0.2? Also, rel_b 1e-3 against η 1e-10: which residual metric is meaningful at large γ? |
| `est_rd1` (+ `_pin`) | unknown truth | degenerate: z fully measured, so γ\* ≈ 0; Mosek-F down to 3.2e-4, never I | redesign the plant so γ\* > 0 has a certified lower end, or eliminate the pinned variable |
| `wellposed_rd1` | tolerance edge | Mosek-F, but cuADMM misses psd_tol by 18% and psd_clip (η 2.53e-7) | tolerance or the missing interior? Does tol 1e-8 move it? |
| `stabpde_rd1` | tolerance edge | direct form; F only at the 50k cap (pinf 1.24e-7 > τ); η 7.24e-8 sets the calibration edge | why is a c = 0 program gap-limited (2,284 it at 1e-4, 33,236 at 1e-6)? Near the direct-form boundary? |
| `nl_fisher_opt` | accuracy limit | plateau at m ≈ 5k; not certified at 25k it (replicated) | budget or true plateau? Records no γ (pass `opts.hi` on Sol) |
| `nl_fisher_gamfix` | accuracy limit | Mosek-OK; never run under the certified standard | one tol-1e-7 run (≈14 min, inferred) decides whether nonlinear has any Sol case |
| `stab2_rd` | unknown truth | psatz off; Mosek fails (rel_b 2.08); a SeDuMi Farkas result is claimed only in memory notes | verify the claimed Farkas certificate on the dump (A'y NSD, b'y > 0) under two row orders |
| `stab2_rd_psz90` | never run | Mosek 1e-8 only in a standalone log | build and bank it through the harness |
| `hinf2_rd` | reference failure | Mosek objective 1.272 against closed form 0.1711, numerr 2; no bracket | is it LPI conservatism or solver failure? Does psatz linear4 help here too? |
| `hinf2nc_rd` | reference failure | 23× the closed form, numerr 2; most expensive case (778 s) | worth anything before `hinf2_rd` has a bracket? |
| `hinf2du_rd` | reference failure | γ = Inf, numerr 2; its dump does not reproduce Mosek (6.5%); banked before the 09-26 fix (`0346ab43`) | re-bank, then check the γ against 0.1711 |
| `scale_hinf_n04` | transient | at 2,000 it n02 and n08 read F but n04 reads i: objective iterates are NOT n-invariant | a converged tol-1e-7 run per n; is the single γ pin row the cause? |

Also kept out: `sentLPI_f1p000000` (exactly 1.000 λ\*). Mosek returns F there with ||x|| 43× larger
than at 0.9999. Its clipped η of 6.5e-8 passes psd_clip, so only the size guard rejects it. Settle
its truth before relying on that guard. Redundant for cuADMM: `stabdual_rd1` (= `stab_rd1` to every
digit), `stabpded_rd1` (= `stabpde_rd1`), `stab2dual_rd` (= `stab2_rd`), `nl_fisher_nobnd`.

### Gaps: classes with no clean case

- **Estimator, controller, coercive H∞, nonlinear PIESOS, 2-D objective**: none certified by
  cuADMM with a known answer (see `problematic`).
- **Direct-form stability and well-posedness**: only tolerance-edge cases.
- **Dual executives on a non-self-adjoint plant**: none. `stabdual_tr1` would fill it.
- **Near-boundary ladder**: the Mosek-F rungs 0.975 … 0.999902 λ\* have no cuADMM runs.
- **Large and coupled**: `scale_rot` (09-27), with the answer from the invariance argument, checked
  with Mosek to n = 8, and n = 16 pending (`S2`). Its n = 24/32 rungs have no Mosek check.
- **Objective above m = 17k, and 2-D with a known answer above m = 4.5k**: none.
- **Certified lower ends on Sol**: cuADMM has produced no exact Farkas certificate anywhere, so
  infeasible ends come from banked Mosek γ_I, theory (`sent101`, the scale sentinels, `hinf2_rd`
  below 0.1711) or the analytic 1-D gain.
