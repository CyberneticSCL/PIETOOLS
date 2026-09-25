# Test suites: what to run after a change

The full registry is ~40 cases and hours with the 2-D and scaling rows. These eight lists are
cut by **the question each one answers**, so a change is checked by the lists whose question it
can break. Lists overlap on purpose — a case belongs to every question it answers.

```matlab
bl_check('smoke')             % run one list against the banked expectations
bl_check('twod','bank')       % run it and record expectations for cases not yet banked
```

Cost is set by what gets *solved*, not by how many cases there are: every 1-D case solves in
4–115 ms, so a 1-D list costs ~45 s of MATLAB startup plus 1–3 s of PDE→PIE conversion per case.

## The lists

| list | question | cases | case time | with startup |
|---|---|---|---|---|
| `smoke` | Is the pipeline alive? | 4 | **13 s** measured | ~1 min |
| `structure` | Did a core change alter any assembled program? | 22 | **53 s** measured | ~1.5 min |
| `executives` | Does every executive family still run? | 9 | ~22 s | ~1 min |
| `known_fail` | Did a change fix, or silently move, a known failure? | 4 | ~6 s | ~1 min |
| `nonlinear` | Is the PIESOS path intact? | 3 | ~4 min | ~5 min |
| `twod` | 2-D assembly and the psatz generators | 4 | ~4 min | ~5 min |
| `objective2d` | The expensive 2-D objective executives | 4 | ~27 min | ~28 min |
| `scaling` | Solver comparison and size questions only | 4 | ~2.7 min | ~3.5 min |

"Case time" is the sum of per-case wall times. It is measured for `smoke` and `structure`
(`results/2026-09-24/check_*.tsv`) and summed from the baseline for the rest. "With startup"
adds ~45 s of MATLAB start and path setup. **`objective2d` exceeds the 20-minute limit for runs
on the shared workstation while its owner is online; schedule it for after 4:30 AM.**

**`smoke`** — one case per structural shape the SDP extraction treats differently: free
variables present (`stab_rd1`, K.f=43), none at all (`stabpde_rd1`, K.f=0), an objective with
`lpi_ineq` so `sossolve` adds blocks *during* the solve (`hinf_rd1`), and two positive operators
with an O(1) right-hand side (`wellposed_rd1`, ‖b‖=17.7, where every other feasibility case sits
at `eppos2` ≈ 1e-6).

**`structure`** — all 22 1-D executive cases, compared *exactly* on m, K.f, block list and
nnz(At). None of those depends on solver behaviour, so any difference means the assembled
program changed. This is the list that cleared the `lpi_eq` fix: 22/22 identical, including
nnz(At).

**`executives`** — one case per executive family (Q-form and direct stability, H∞ primal and
dual, H2 in both gramians, estimator, controller, well-posedness). Only checks that each builds
and solves.

**`known_fail`** — the 1-D cases whose *current* behaviour is a documented failure, so the
expected result *is* the failure: `h2cco_rd1` (rel_b = 1.000 with numerr = 0, which the status
word does not catch), `hinfduco_rd1` (γ = 8251.9), and the two controllers (rel_b ~1e-3,
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
