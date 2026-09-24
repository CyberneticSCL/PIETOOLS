# The multi-executive interface and the benchmark suite

Written 09/23/2026 (CC).

`README.md` says what the package claims, `GUIDE.md` how to use the 2-D stability demo, and
`ALGORITHM.md` describes the mechanism and lists open questions. This file covers the
**interface added on 09/23/2026** — `pielr_solve` over several PIETOOLS executives — and the
**benchmark suite** that exists so a change to the solver can be measured rather than argued.

Everything below that is stated as a number was measured in the run that produced it. Where a
value comes from somewhere else it is marked as such and not relied on.

---

## 1. Why this exists

Every measurement the package had recorded was taken on **one program**: a 2-D
reaction-diffusion at one degree setting, one operating point, one executive. A change measured
there says nothing about its effect anywhere else, and several of the package's own open
questions are not answerable from one program at all.

So, in order:

1. `pielr_solve(PIE, lpi, opts)` — one entry point over several executives and both dimensions.
2. `pielr_bench` — a catalogue of problems, each scored against an interior-point reference
   **measured in the same run** and judged by the **same gate**.
3. `pielr_bench_diff(Ra, Rb)` — the diff between two runs. This is the point of the other two.

---

## 2. Which executives, and why those

PIETOOLS has 27 executives. The selection criteria, all checked against the source:

- the executive must build its constraint through `lpi_eq(..., X + slack, 'symmetric')`, so the
  constraint lives in a Gram block and the residual is checkable at the operator level;
- there must be a non-collapsing quantity to divide the residual by (see §4);
- a reference answer should exist that is independent of the method.

| LPI | status | why |
|---|---|---|
| `stability` | **done**, 1-D + 2-D | the existing target; direct form |
| `stability-dual` | **done**, 1-D + 2-D | differs from the primal in three lines |
| `l2gain` | **done**, 1-D | carries γ, so it scores on a continuous error, not a binary gate |
| `l2gain-dual` | **done**, 1-D | ditto |
| `l2gain` 2-D | not adapted | `PIETOOLS_Hinf_gain_2D` is a *different formulation* from the 1-D one — `Pop` directly, with `Twop`, rather than a free `lpivar` `Qop` coupled by `Top'*Qop = Rop`. It needs its own transcription, not a re-dimensioned copy. |
| `h2norm`, `h2norm-dual` | next | same shape as `l2gain`, more blocks (`pos=7, eq=3`), four coercive/non-coercive variants |
| `well-posedness` | next | feasibility, `pos=8`, no reference scalar |
| the four synthesis executives | deferred | they carry extra non-positivity decision operators and the deliverable is a recovered controller/observer, so "did it certify" is the wrong question; needs a closed-loop harness |

**`l2gain` turned out to be the most valuable arm**, for three reasons that only became clear
once it ran: the reference is a number rather than a yes/no; γ is a free decision variable, so
the program has `Kf > 0`, which no stability program does and which exposed three latent bugs
(§6); and comparing γ to the interior-point γ measures how much the low-rank restriction
actually costs.

### Two deliberate differences from `lpiscript`

1. **`'stability'` here is the DIRECT form** (`PIETOOLS_PDEstability` / `PIETOOLS_stability_2D`).
   `lpiscript('stability')` calls `PIETOOLS_PIE2PDEstability`, the Q form — a different and
   weaker LPI. Both current drivers in this package build the direct form, so the name is kept
   and documented rather than silently pointing at the other one.
2. **`lpiscript` has no 2-D dispatch at all** — every case routes to a 1-D executive. Here the
   dimension comes from the PIE and the adapter picks the matching transcription.

---

## 3. Structure

```
pielr_solve(PIE,lpi,opts)
  pielr_lpi(lpi)              adapter: build / residual / normaliser, per executive
    build_stab_1d             PIETOOLS_PDEstability(_dual), minus the solve
    build_stab_2d_st2         PIETOOLS_stability_2D(_dual), minus the solve
    build_l2gain_1d           PIETOOLS_Hinf_gain(_dual), minus the solve
  pielr_rawdata               equality rows + objective + layout
    pielr_layout              the x-vector map, read from the program
  pielr_discover              rank ladder, DR start, LM, face certification
    bm_setup/bm_proj/bm_dr/bm_lm2/bm_resid     unchanged core
    restrict_solve                             layout-aware, adapter-gated
  pielr_bisect_obj            bisection on gamma, for the LPIs that carry one
  pielr_opcheck               THE gate
  pielr_refine                shrink the accepted face
```

Everything executive-specific is in the adapter: an assembly, a residual expression, and a
normaliser. Everything else is dimension- and executive-agnostic. That is what lets one core
serve four LPIs in two dimensions; the package previously forked its whole core into `tests_1d/`
to reach 1-D, and two of those seven files have since drifted out of sync with `private/`.

---

## 4. The one behaviour change: the gate's denominator

This is the only thing about the solver that changed, and it was forced.

`tests_1d/gate1d.m` accepts on `relP = max|Res| / max|Pop|`. `opcheck_2d.m` accepts on
`max|Res| / max|Qop|`. The two halves of the package accepted on **different quantities**, and
`tests_1d/opcheck.m` states in its own comments why the 2-D choice is wrong:

> `rel` is INVALID wherever the candidate drives `Dop -> 0` as an operator: it then reads 0/0 and
> improves for a purely artifactual reason (measured on the Example 21 beam,
> `max|Dop| = 2.4e-10` manufactured a clean degree threshold that does not exist). `relP` divides
> by `max|Pop|` instead, which cannot collapse because `Pop >= eppos2*I` by construction.
> **Accept on `relP`, never on `rel` alone.**

A suite spanning both dimensions cannot be built on two different criteria. `pielr_opcheck`
takes the normaliser from the adapter, which must also state **why it cannot collapse**:

- `stability`: `Pop = poslpivar(..) + eppos*I`, so `max|Pop| >= eppos`.
- `l2gain`: `Dop` carries `-gam*Iw` and `-gam*Iz` on its diagonal, so `max|Dop| >= gam > 0`.
  `Rop` has **no** `eppos*I` in this executive and must not be used.

The old quantity is still reported as `rel_d`, and `maxRes`, `maxNrm` and `maxDop` are now
returned alongside. `opcheck_2d` computed those on every call and returned only the ratio, which
is why no existing log can say whether two residuals differ in the numerator or the denominator.

Nothing else about the solver was changed, so the existing behaviour remains the baseline.

---

## 5. The benchmark

```matlab
pielr_path
Ra = pielr_bench(struct('tier',1,'tag','baseline'));
% ... change one thing ...
Rb = pielr_bench(struct('tier',1,'tag','with-the-change'));
pielr_bench_diff(Ra,Rb)
```

**No reference value is stored anywhere.** For each case the *identical assembled program* goes
to the interior-point solver in the same run, and **both answers are judged by the same gate**.
A solver's own feasibility flag never appears. The `.note` fields in the catalogue carry
literature values as context; nothing scores against them.

Tiers, because the 2-D cases cost hours and the rest cost seconds:

| tier | adds | cost |
|---|---|---|
| 1 | 1-D PDE and DDE, 17 cases | ~8 min total |
| 2 | 2-D PDE at degree 2 and 3 | minutes each |
| 3 | 2-D PDE at degree 4 — the program every existing measurement was taken on | hours |

The catalogue varies the three axes that move the problem: the **executive** (feasibility vs an
objective), the **model type**, and the **operating point / degree**.

**The DDE family is new here and matters structurally.** The 2-D benchmark prints
`dims [n0 nx ny n2] = [0 0 0 1]`, so nothing measured before touched an operator with a
finite-dimensional part, and the square-root mechanism behind the rank claims says nothing about
one. DDE-derived PIEs have a substantial ODE block and cost seconds.

### Baseline, tier 1, `lmit = 400`, 4 threads

```
cases 17 | IPM certifies 14 | low-rank certifies 11 | agree 11 | LOW-RANK MISSED 3
gamma ratio (low-rank / interior-point, 1.000 is tight):
  gain-heat-dist       1.0010
  gain-transport       1.0819
  gain-transport-dual  1.1091
  gain-reacdiff        2.8264
```

- **Missed 3**: `rd1d-lam0.5`, `rd1d-n2`, `advdiff1d` — all three interior-point-certifiable at
  `rel` 9.5e-08 / 2.7e-07 / 9.2e-09, so they are search failures, not infeasible programs.
- **3 further cases** (`rd1d-lam0.9`, `advdiff1d-dual`, `lib-heat-ode-das`) the interior-point
  method cannot certify either. Those are not counted against the method.
- All four DDE cases certify, at ranks `[1 1 1]`, `[2 2 2]`, `[2 2 2]`, `[3 3 3]`.
- Every γ ratio is `>= 1`, as a sound upper bound must be.

### First result: the budget sweep, and what the exit reasons say

`lmit ∈ {400, 2000, 8000}`, tier 1, same machine, 4 threads. This is the experiment the
package could not previously run, because `bm_lm2` returns its exit reason and no caller
captured it.

| | lmit 400 | lmit 2000 | lmit 8000 |
|---|---|---|---|
| low-rank certifies | 11 / 17 | 11 / 17 | **12 / 17** |
| missed (IPM certifies, low-rank does not) | 3 | 3 | **2** |
| γ ratio, worst | 2.8264 | 1.0163 | **1.0006** |
| γ ratio, median | 1.0955 | 1.0104 | **1.0002** |

**γ accuracy is budget-limited; feasibility verdicts largely are not.** A 5× budget takes the
worst γ from 2.83× to 1.02×, and 20× to 1.0006×, while the verdict column barely moves.

**The reported rank was also budget-dependent**, which matters because rank is the headline
number of every low-rank claim here: `gain-transport` `[3 3 3] → [2 2 2]`, `gain-transport-dual`
`[4 4 4] → [3 3 3]`, `gain-reacdiff` `[3 6 5] → [2 3 3]`. A rank reported at `lmit = 400` is
partly an artefact of the cut.

The **LM exit census** separates the three misses into three different causes:

| case | 400 | 2000 | 8000 | reading |
|---|---|---|---|---|
| `advdiff1d` | `maxit:21` | `maxit:13` | `maxit:1` | purely **cut**; certifies at 8000 |
| `rd1d-n2` | `maxit:18 stagnant:5` | `maxit:6 stagnant:11 tol:6` | `maxit:2 stagnant:11 **tol:10**` | not budget-limited: it **exits on `tol` ten times and still fails the gate** |
| `rd1d-lam0.5` | `maxit:16` | `maxit:12` | `maxit:3` | mostly converged, still short |
| `lib-heat-ode-das` | `maxit:18 stagnant:5` | `maxit:7 stagnant:16` | `stagnant:23` | fully **flattened**; the IPM cannot certify it either |
| `rd1d-lam0.9` | `maxit:20` | `maxit:12` | `damping:2 stagnant:21` | flattened; IPM also fails |

`rd1d-n2` is the sharp one. Ten attempts exit via `tol` — the whitened residual crossing
`1e-6` — and the operator gate still rejects them. That is the units mismatch in §6 with an
operational consequence, not just a discrepancy on paper: the stopping rule is stopping on a
quantity that is not the acceptance criterion.

**One caveat on the γ ratios.** `gain-reacdiff` reads **0.9950** at `lmit = 8000` — below 1,
where a sound upper bound cannot go. The explanation is in the same row: the interior-point
point has `ipm_rel = 9.181e-07`, i.e. it only just clears the `1e-6` gate itself. The reference
is judged by the same gate as the candidate, so on a hard case it is only gate-accurate, and a
ratio within roughly the gate tolerance of 1 is not meaningful. Read ratios as "1.0 ± gate",
not as exact.

### Second result: 1e-6 is a round number, not a calibrated threshold

Measured on the **historical benchmark configuration** — `nb_rd2d(1, 0.50·λ*)`, `set2d_deg(4,[])`,
`N = [8 744]`, `m = 5920`, which is the exact program every prior measurement in this package
was taken on — using the **interior-point** solution, i.e. the full-rank answer:

| quantity | value |
|---|---|
| `maxRes` (absolute residual) | 4.7164e-11 |
| `max\|Pop\|` | 1.3486e-05 |
| `max\|Qop\|` | 6.6555e-05 |
| `rel = maxRes / max\|Pop\|` (corrected denominator) | **3.4973e-06 — FAILS 1e-6** |
| `rel_d = maxRes / max\|Qop\|` (prior denominator) | **7.0865e-07 — passes** |
| ratio | **4.94** |

Same point, same absolute residual; the verdict flips on the choice of denominator alone.
`max|Qop|` is about five times `max|Pop|` here, so the prior gate is about five times more
lenient on this program.

**What follows, and what does not.** It does *not* show the previously reported rank-[3 4]
result was wrong. It shows that a *relative* residual threshold is only meaningful together
with its denominator, and that 1e-6 was never calibrated against the corrected one. Under the
corrected denominator **not even the exact interior-point solution clears 1e-6 at this
operating point** — so a 1e-6 gate there would be asking the low-rank search to beat the full
solve, which is not an acceptance criterion anyone should want.

**The fix the suite makes possible.** `ipm_rel` is measured for every case, and across tier 1 it
spans 4.15e-10 to 9.18e-07 on the cases that certify — three and a half orders of magnitude. A
fixed absolute threshold cannot be right across that range. The acceptance test should be
relative to what the reference achieves on the same program — score `rel / ipm_rel`, or accept
at `rel <= K · ipm_rel` — and the benchmark computes both sides in the same run. That is a
change to the gate, so it is **not** made here: the point of this pass was to establish the
baseline, and re-tuning the threshold is the first change the suite should be used to evaluate.

Related, and also measured: at degree 4 the IPM solve costs ~880 s at every operating point
tried (883 / 865 / 884 s at frac 0.10 / 0.25 / 0.50), and at degree 3 a single low-rank run
cost 1025 s against the same program's 117 s IPM solve — a factor **8.7** the wrong way.

### Third result: tier 2, and how to choose the gate's `k`

The gate is now `rel <= max(abs, k·ref)` where `ref` is credible (`ref <= refmax`), and the
reported score is `rel/ref`. Two constants to pick.

**Tier 2, `lmit` 400, `maxrank` 4, 2 seeds:**

| case | lr | ipm@1e-6 | rel | ipm_rel | score | rank |
|---|---|---|---|---|---|---|
| `rd2d-deg3-f010` | yes | no | 1.114e-05 | 1.900e-06 | 5.863 | **[2 2]** |
| `rd2d-deg3-f010-dual` | yes | yes | 3.686e-06 | 8.186e-07 | 4.503 | **[2 3]** |
| `rd2d-deg2` | no | no | — | 4.935e+00 | — | negative control |
| `rd2d-deg3` | no | no | — | 6.283e-05 | — | negative control |

Two 2-D cases certify at rank [2 2] and [2 3], within 6× and 4.5× of the reference residual.
Under the old absolute gate the primal would have failed at 1.114e-05 — **and so would its own
reference**. That is the case the change was made for.

**Choosing `k`, measured closed-loop on tier 1.** Closed-loop matters: `k` changes *which point
the search returns*, since the ladder stops at the first rung that certifies and the bisection
keeps whatever certifies. Re-scoring banked points answers an easier question.

| k | certs | binds | score median | γ ratios |
|---|---|---|---|---|
| 1 | 11 | 0 | 81.29 | 1.0819 1.1091 1.0010 2.8264 |
| 10 | 11 | 2 | 81.29 | 1.0819 1.1091 1.0010 1.6856 |
| 100 | 13 | 6 | 95.39 | 1.0779 1.1083 1.0010 **0.6143** |
| 1000 | 14 | 11 | 667.6 | 1.0148 1.0410 1.0010 **0.1936** |

`k` binds only where `k·ref > abs`, i.e. `ref > 1e-6/k`, so at tier 1 `k=1` binds nowhere and is
exactly the old absolute gate.

**The γ column falsifies the large values.** A face restriction can only *lose* feasible points,
so the restricted optimum cannot fall below the full optimum — a ratio under 1 is a
contradiction, not tolerance. At `k` = 100 and 1000 the gate accepts points that do not satisfy
the LPI and the bisection then drives γ down freely. **`k ≤ 10` validated, `k ≥ 100` refuted, by
the suite's own data.** This is now a standing assertion: any γ ratio < 1 prints `*** UNSOUND`
and names the cause. The l2gain cases therefore calibrate `k` automatically, which is a reason
to keep one in any run that touches the gate.

**`refmax`, and why the first choice was wrong.** It was set to 1e-4 on tier-1 evidence —
credible references ran to 9.181e-07, non-credible ones started at 1.500, six orders of gap. The
2-D cases populated that gap. `rd2d-deg3` at frac 0.50 is a negative control whose reference
reaches only 6.283e-05; at `refmax = 1e-4` the relative branch applied, the threshold became
6.283e-04, and a point at `rel = 1.020e-04` was accepted — **with score 1.623, the best in the
suite.** That is the lesson in one row: a good ratio on a program nobody can solve is not a good
result.

The binding quantity is **`k·refmax`**, the loosest threshold the gate can ever apply, not `k`
alone. `refmax = 1e-5` is measured from both sides: admits 2-D deg 4 frac 0.50 (`ref` 3.497e-06,
the case that motivated all of this) and 2-D deg 3 frac 0.10 (1.900e-06); rejects 2-D deg 3
frac 0.50 (6.283e-05). With `k = 10` the worst threshold is 1e-4. Tier 1 is unaffected — no
tier-1 reference sits in (1e-5, 1e-4) — and was re-run to confirm.

**One reporting caveat.** The reference column is `ipm@1e-6`, not `ipm`: the two columns answer
different questions and must not be read as a head-to-head. The low-rank verdict is the relative
gate; the reference is reported against the plain 1e-6, which asks whether 1e-6 is reachable on
that program at all. Judging the reference by the relative gate would be vacuous, its own
residual being the reference. A row with `lr = yes, ipm@1e-6 = no` is expected wherever 1e-6 is
out of reach; the comparison is the score.

---

## 6. What building this found

Each of these is reproducible from the source or from a run; none is inferred.

**`raw_data`'s layout model is wrong whenever `lpivar` is used.** It assumes
`[free prefix][Gram blocks contiguously]`. A free `lpivar` operator is stored as `'poly'` entries
in `sos.var`, **interleaved between** the `'sos'` entries. Measured on the 1-D l2gain program:
`raw_data` returns `Kf=1, Ns=[10 17 8]`, accounting for 454 of 497 coordinates. Every Gram block
after the first is then read at the wrong offset.

**`raw_data` concatenates inequality rows as equalities.** `sos.expr.type{1}` is `'ineq'` on that
program — the `gam >= 0` constraint. Treating it as an equality pins γ instead of bounding it.

**`raw_data` drops the objective entirely**, so an optimisation executive silently degrades to
feasibility.

**`prog.solinfo.x` is not the decision vector.** It is SeDuMi's primal for the cone the *solver*
saw: `sossolve` turns each `'ineq'` expression into a slack cone variable, so on the 1-D l2gain
program `x` has 498 entries against 497 decision variables, in the solver's ordering. The rows of
`expr.At` — which `bm_setup`, the layout and the gate all index against — are in **decvartable**
order, i.e. `solinfo.RRx`. Using `x`, the executive's *own* solution violated its *own* equality
rows by 1.83e+01 and `x(jc)` read 2.562 for a program whose reported γ is 0.516333. The two
coincide only when the program has no free variables and no inequality — exactly the stability
case that had been the only thing tested (`||x - RRx|| = 0` over 361 entries there).

**Three `Kf > 0` bugs**, latent because every stability program has `Kf = 0`:
`wblk` starts its offset at 0 where `bm_resid` starts at `Kf`; the Douglas–Rachford start omits
the free coordinates entirely; and `padw` pads the rank-ladder warm start with **zero** columns,
whose Jacobian block `2*W*Ssym(:,rows)*kron(Y(:,c),I)` is identically zero, so LM can never move
them and the rung reproduces the previous rung's point. The package's own T4 pins that last one
as defect B2 and `fitw` is its tested repair.

**Nothing in the pipeline optimises anything.** `bm_lm2` minimises `||W(A(X)-b)||`, in which the
objective vector `c` does not appear, and `restrict_solve`'s determined solve has `np == rM`, so
a face has no degree of freedom left. Run straight on l2gain the search returned γ = 1.42095
against the interior-point 0.516333. Bisecting γ — which turns the optimisation LPI back into the
feasibility LPI the machinery is built for — gives 0.558528, and each certified γ is an upper
bound, so every bisection iterate is already a valid answer.

**`bm_lm2`'s `tol` exit is in the wrong units.** It compares the **whitened** residual `||Wr||`
against a tolerance the header calibrates against the **unwhitened** `raw_rel`. Measured: exit at
`||Wr|| = 9.99e-07 < 1e-6` with `raw_rel = 1.26e-06` still above it, a factor 1.26 on that
program. The factor is bounded only by the conditioning of `W`.

**`examples_PDE_library_PIETOOLS` cannot be driven non-interactively.** It ends in an unguarded
`input()` prompt (line 976). The example files under `Examples_Library/` take `(GUI, params)`
directly and are what a harness should call.

**MATLAB R2025b crashes** with an access violation in its JIT in two situations, both hit here
and both reproducible, and both on the *second* case of a run:

- `evalc` inside a function invoked through a function handle stored in a struct array (the
  original shape of the benchmark loaders; trace names `inEvalCmdWithLocalReturn`). A single
  `evalc` at a fixed call site is fine.
- calling `pielr_path` — which re-runs `pietools_path_update`, i.e. `addpath(genpath(...))`
  followed by a `path()` rewrite — and *then* running the suite in the same process, after
  functions have been JIT-compiled (trace names the `LXE compiler thread`).

Neither is a package defect. Run `pielr_path` in its own process, or set the path directly, and
do not put an `evalc` behind a stored handle. 45 files parse clean under `mtree` in both cases,
so these are not syntax faults.

---

## 7. What is *not* done

- 2-D `l2gain` (a separate formulation; see §2).
- `h2norm`, `well-posedness`, and the four synthesis executives.
- The γ bisection rebuilds `bm_setup` per trial, and with `pre=1` that is a **dense m×m
  eigendecomposition**. At 1-D sizes it is free; at 2-D sizes it is the dominant cost and a
  12-step bisection multiplies it by 12. Only one column of `At` changes between trials, so this
  is a rank-one update in principle. Unexploited, and the first thing to fix if the 2-D gain arm
  is wanted.
- `pielr_certify` is untouched and remains the behaviour baseline for the 2-D stability path.
  It still has the zero pad, the old gate denominator, and the discarded diagnostics. Retiring it
  should wait until the suite has measured what each of those is worth.
- `tests_1d/` still holds a forked copy of the core, and `bm_lm2.m` and `bm_resid.m` have
  diverged from `private/`. By that suite's own `t1d_corefiles` contract ("If these diverge from
  ../private/, every other result here becomes a statement about a copy"), T0 fails until they
  are re-synced.
