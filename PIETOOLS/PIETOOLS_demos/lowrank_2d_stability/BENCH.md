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

---

## 8. The objective is the wrong question to ask the solver

Raised by the maintainer, then by the parallel cuADMM session: *are you using an objective in
these test cases?* The answer turned out to be "in one arm only", and that asymmetry is the
mechanism behind the single unsoundness the suite has measured.

### The asymmetry

Nothing in the low-rank pipeline optimises. `bm_lm2` minimises `||W(A(X)-b)||`, in which the
objective vector `c` never appears, so `pielr_bisect_obj` **pins γ and tests feasibility**,
bisecting to recover the optimum. The reference, however, called `lpisolve` with
`lpisetobj(prog,gam)` still active — it *minimised*. The two arms were not solving the same kind
of problem, and once `score = rel/ipm_rel` made the reference residual set the acceptance
threshold, that stopped being cosmetic.

### Why it is one fault and not two

The natural reading is that there are two independent faults: an inflated reference (explaining a
loose *accept*) and something else entirely (explaining a γ *below* the reference, which a face
restriction cannot produce). That reading is wrong, and the bisection loop says why:

```matlab
if ~isempty(Rk) && Rk.ok, hi = g; else, lo = g; end
```

A trial is accepted **iff the gate accepts it**. So the bisection does not converge to the true
optimum — it converges to *the lowest γ the gate will still admit*. Threshold slack is not merely
permissive, it is a **search direction**: the bisection actively hunts for the most optimistic γ
inside the slack. Inflating the reference therefore does produce a sub-1 ratio, by exactly this
route, and fixing the reference removes the cause rather than hiding it.

The generalisable rule, which the cuADMM session stated after this exchange and which its own
harness already satisfies by using an absolute `rel_b < 1e-6`: **a bisection's per-trial accept
test must not depend on a quantity that varies with solver behaviour.** Absolute, or relative to
something structural like `||b||` — never relative to another solve's residual, or the tolerance
becomes an optimisation direction.

### The fix

`st.gamfix` on both objective builders (`build_l2gain_1d`, `build_poincare_1d`) **substitutes** a
numeric γ rather than pinning it with an appended row. PIETOOLS documents this route itself
("a specific gain test ... results in a feasibility test instead of an optimization problem").
Substitution rather than pinning because pinning leaves the `γ >= 0` cone block and the objective
in the data while making that block redundant, and redundant equality rows are a known
conditioning hazard in this stack. The structural signature distinguishes them, measured by the
cuADMM session on its own build: substitution takes `Kf` 44→43, `Ns` `[10 17 8 1]`→`[10 17 8]`
and `m` 150→149 — *down* by one; pinning would take `m` up by one.

`pielr_ipm_ref` takes an optional rebuild handle. Given one, it solves with the objective to get
γ\*, then re-solves at fixed γ\* as a pure feasibility program and reports **that** residual,
keeping γ\* from the optimising solve. Both are retained (`rel_obj` and `rel`) so the inflation is
a measured per-case quantity rather than a constant anyone quotes. Without the handle the
behaviour is exactly as before.

### The floor, which this puts in question

With the reference bisected, `max(abs, k*ref)` can collapse onto the absolute floor and the gate
silently stops being relative. Whether it does depends on the inflation factor, and the two
candidate numbers give opposite answers on `gain-reacdiff` tier 2 (reference 4.2845e-07, BM
accepted at rel 4.280e-06):

| inflation | bisected ref | threshold `max(1e-6, 10·ref)` | verdict on BM | what binds |
|---|---|---|---|---|
| 10× | 4.3e-08 | 1e-6 | reject | the **floor**; `k` never binds |
| 3.5× | 1.22e-07 | 1.22e-06 | reject | **`k`**; the floor never binds |

BM is rejected either way, but *what rejects it* differs, and only one of those says the floor
needs retuning. The cuADMM session first reported ~10× and then corrected itself to 2.6×–6.3×
(clustering ~3.5×) on n = 1 plant at 2 settings — too weak to settle it. This is why
`pielr_ipm_ref` records `rel_obj` and `rel` separately: **the floor is decided after the
18-case distribution, not before it.**

Independent of the inflation factor, there is a separate reason to suspect 1e-6. Achievable
fixed-γ residuals on these programs, measured by that session: SeDuMi 3.0e-10 to 1.1e-08, Mosek
through PIETOOLS 3.3e-08 to 7.8e-08; and the three 1-D stability cases reach 2.6e-08 to 4.9e-08
here. A floor of 1e-6 therefore sits one to two orders above what both arms can actually achieve,
so on the easy cases it binds and the gate reverts to the fixed threshold it was moved away from.
A floor should be anchored to achievable residual, not to a round number.

### Adjudicating BM's sub-1 γ without trusting either arm

Feasibility at fixed γ is objective-independent, so it settles whether BM's γ is achievable
without reference to any objective value. `gamfix_probe` rebuilds the identical program at a
ladder of γ spanning BM's answer and the reference and asks the interior-point solver at each
rung. Two outcomes, opposite conclusions:

- IPM **infeasible** at BM's γ → the accept is spurious; gate slack produced an invalid bound.
- IPM **feasible** at BM's γ → the reference was simply not optimal and BM is right.

The rungs are judged at the **absolute** 1e-6 (`pielr_ipm_ref` never sets `g.ref`), so the ladder
does not inherit the fault it is diagnosing — and the residuals are printed alongside the verdict,
because at 1e-6 the boolean is far coarser than the number.

### What the measurements actually showed — and it was not the reference

Two hypotheses went into this, both plausible, and the measurements refuted both.

**Reference inflation does not reproduce.** `rel_obj / rel_fixed` at the reference γ on
`gain-reacdiff`: **0.98×** at tier 2 and **0.55×** at tier 3. The fixed-γ solve is the same or
*worse*, against the 2.6–6.3× the cuADMM session measured on its own plant. That session has
since withdrawn the figure as plant-specific. The `gamfix` route stays — a like-for-like
reference is right on principle — but it is **not** the fix and is not claimed as one.

**The fixed-γ ladder locates the true boundary.** γ substituted, every rung judged at the
absolute 1e-6 so the ladder cannot inherit the fault it is diagnosing:

| γ | vs ref | rel | maxRes | feasible |
|---|---|---|---|---|
| 7.691067 | 0.9364 | 9.7757e-04 | 9.2412e-05 | no |
| 7.933943 | 0.9659 | 5.2542e-04 | 4.9274e-05 | no |
| **8.095860** | **0.9856** | **2.2808e-04** | 2.1463e-05 | **no** ← BM's answer |
| 8.154864 | 0.9928 | 1.1703e-04 | 1.0840e-05 | no |
| 8.213868 | 1.0000 | 4.2894e-07 | 1.3713e-07 | **yes** ← reference |
| 8.378145 | 1.0200 | 1.0147e-08 | 9.1608e-10 | yes |
| 9.035255 | 1.1000 | 9.4383e-09 | 9.4150e-10 | yes |

Monotone below the reference with a four-order cliff exactly at it. The reference **is** the
feasibility boundary and BM's γ is inside the infeasible region.

**The controlled experiment: the threshold is the whole cause.** Identical bisection, identical
settings, identical seeds; the only variable is the acceptance threshold.

| threshold | γ | γ / ref | sound |
|---|---|---|---|
| absolute 1e-6 | 8.216934 | **1.0004** | yes |
| 10·ref = 4.28e-06 | 8.095859 | **0.9856** | no |

**This refutes `k = 10`, not merely `k >= 100`.** §5's "k = 10 validated" is withdrawn: it was
inferred from which cases produced γ < 1 at each `k`, which is a weaker test than varying `k`
on one case with everything else held fixed. `score = rel/ipm_rel` is retained as a
**diagnostic**; the gate goes back to an absolute threshold.

### The denominator is a second, independent source of slack

One number resisted explanation: at γ = 8.095859 the low-rank point has operator `maxRes`
4.2005e-07 while the IPM's best at the same substituted γ is 2.1336e-05 — apparently beating a
convex solve 50× on an infeasible program. Two explanations were tested and refuted here:
**denominator inflation between the two points** (maxNrm BM/IPM = 1.043, not ≫1) and
**target-vs-achieved γ** (`cert.gam` is overwritten from `R.aux.gam` at `pielr_solve:251`, so the
reported γ is read back out of the solution, not the bisection's target).

The cuADMM session supplied the explanation, with both quantities on the same point from its own
panel: three genuinely different **numerators** (induced L2, Hilbert–Schmidt, coefficient-space
Frobenius) agree within **1.5×**, while changing only the **denominator** moves the answer
**23–65×**. The ratio of operator to row residual ranged 0.15× to 65× across its configs — a 430×
spread, in both directions.

The reason the maxNrm comparison could not detect this: it compared *this gate's own* denominator
between two points, which is stable by construction. The row residual's denominator is `||b||`,
**fixed by the program**, whereas `||Dop||` and `||Pop||` are **solution-dependent** — that session
measured `||Dop||` swinging 24× (3.94e-06 vs 1.62e-07) from a settings change with the physics
unchanged. Two measures whose denominators are respectively constant and solution-dependent cannot
track each other, and they diverge most where the solution is unusual — exactly the regime a gate
adjudicates.

So the "BM beats the IPM" reading is dropped; the ground truth was never in doubt, since the ladder
establishes infeasibility independently of any normaliser.

**This leaves an open design question, and it is now the important one.** Both normalisers here are
solution-dependent (`max|Dop|` on the negativity row, `max(|Top'*Qop|, |Rop|)` on the coupling row),
so the gate inherits that swing. The absolute-threshold fix addresses the *threshold*; it does not
address the *denominator*. `||b||` is fixed but lives in coefficient space, and the reason the gate
works in operator space at all is that coefficient-space residuals failed to detect operator-level
nonsense (§4). A fixed operator-space normaliser is the thing to find; there is not one yet.

### Closing the 50×: the gate inverts the ranking

`twores_u.m` reproduces the unsound accept (the relative gate attached to the adapter, which is
what `pielr_solve` does and a direct `pielr_bisect_obj` call does *not* — that asymmetry is why
the earlier direct run returned the sound 8.216934) and scores its **row** residual. Both arms at
γ = 8.095859 on the same program:

| point | row ‖A'x−b‖/‖b‖ | gate operator rel | op/row |
|---|---|---|---|
| low-rank | **3.0702e-05** | 4.2800e-06 | 7.17 |
| IPM | **9.3484e-06** | 2.2676e-04 | 0.0412 |

**BM never beat a convex solve.** In the quantity the SDP constrains, the interior-point point is
**3.3× better**. The 50× was an artifact of the operator measure.

Worse, the gate gets the ordering backwards: it **accepts** the low-rank point (4.28e-06, just
inside threshold) and **rejects** the interior-point one (2.27e-04) at the same γ, while the row
residual says the rejected point is the less-violating of the two. The operator residual
understates BM's violation 7.2× and overstates the IPM's 24× — a 174× disagreement between the
two arms on one program. This is the branch the test's own decision rule named: *the gate and the
SDP disagree, and the gate is the defect.* The threshold fix is necessary but not sufficient.

Two corollaries:

- **The program is genuinely infeasible at that γ.** Both row residuals sit at ~1e-5 against
  ~1e-8 at a feasible γ (8.378145), independently of any normaliser.
- **It is not the denominator.** maxNrm 9.81e-02 vs 9.41e-02, ratio 1.04. The whole gap is in the
  **numerator**: `max|Res|` over operator coefficients versus `‖A'x−b‖`, two measures of the same
  residual **165× apart**. That bounds the cuADMM session's decomposition (numerator ~1.5×,
  denominator 23–65×) — its three numerators were all *operator* norms of the same reconstructed
  residual and so were bound to agree. The row residual is not an operator norm, and the gap
  lives exactly there.

Scope of the claim: one case, two points. The contrast with the **sound** point matters — there
op/row was 0.133 and 0.129, consistent to 3%. So the measures track near feasibility and diverge
*and invert* at an infeasible point, which is precisely where a gate operates. That is a
hypothesis from two data points, not an established law, and it is the next thing to test across
the suite.

**Consequence for the package:** `max|Res|` on the reconstructed operator must not be treated as a
proxy for SDP feasibility. The row residual `‖A'x−b‖/‖b‖` should be computed and reported
alongside it, and acceptance should not rest on the operator residual alone. It also means an
operator residual and a `rel_b` from another session are not comparable quantities.

### Does the inversion generalise? Yes — and on exact comparisons

`inversion.m` puts a low-rank point and an interior-point point **on the same system** for every
case and scores both residuals on both. For the feasibility LPIs there is no objective, so the two
arms solve the *identical* program and the row residuals are residuals of one linear system —
an exact comparison with nothing to reconcile. Inversion is declared when `r_op = op_bm/op_ipm`
and `r_row = row_bm/row_ipm` fall on opposite sides of 1, i.e. the two measures disagree about
which of the two points is better.

| case | r_op | r_row | inversion | op/row (BM, IPM) | within-case spread |
|---|---|---|---|---|---|
| rd1d-lam0.5 | 120 | 3.56 | — | 3.67, 0.109 | 34× |
| rd1d-lam0.9 | 63.7 | 3.97 | — | 1.87, 0.117 | 16× |
| rd1d-n2 | 2.22 | **0.0734** | **INVERT** | 2.29, 0.0756 | 30× |
| advdiff1d | 4.72 | **0.606** | **INVERT** | 2.64, 0.339 | 7.8× |
| advdiff1d-dual | 1893 | 555 | — | 0.0185, 0.00542 | 3.4× |
| lib-transport | 47.1 | 3.32 | — | 0.55, 0.0387 | 14× |

**Two outright inversions in six**, on identical linear systems — so the `gain-reacdiff` result was
not a one-off, and these carry no cross-system caveat at all.

The stronger statistic is the last column. `op/row` is **never constant**: it differs by 3.4× to
34× between two points on the *same program*. If the operator and row residuals were one quantity
in different units this ratio would be fixed and inversion would be impossible. It is not, so
inversion is a structural possibility of the gate rather than an accident of one case.

**A complication, recorded because it cuts against a tidy story.** On all six of these, BM's
`op/row` is *higher* than the IPM's — the operator measure consistently makes the low-rank point
look relatively worse. On `gain-reacdiff` it was the reverse (0.139 vs 24.3). So the bias is not
even consistent in sign. That case compared two systems differing by the γ column while these are
exact, so the two are not strictly comparable; but it means "the operator residual is optimistic"
would be the wrong lesson. The right one is that it is *unrelated* to the constraint violation at
the scale that matters.

**Run note.** Both this and the tier-3 block of the absolute-gate sweep were cut short by the
R2025b JIT access violation (`evalc` around a handle stored in a struct array — the pattern
already recorded in §6). Tiers 1 and 2 completed and are reported; tier 3 is **not** a result and
is not quoted. The inversion run was re-executed in crash-isolated chunks to cover the rest.

#### The complete table (supersedes the six-case partial above)

Re-run in crash-isolated chunks after the JIT fault. `same = EXACT` means both arms solved the
*identical* program; `gam` means γ was fixed at BM's own answer, so the systems differ by the γ
column only.

| case | same | r_op | r_row | inv | op/row (BM, IPM) |
|---|---|---|---|---|---|
| rd1d-lam0.5 | EXACT | 120 | 3.56 | — | 3.67, 0.109 |
| rd1d-lam0.9 | EXACT | 63.7 | 3.97 | — | 1.87, 0.117 |
| rd1d-n2 | EXACT | 2.22 | 0.0734 | **INVERT** | 2.29, 0.0756 |
| advdiff1d | EXACT | 4.72 | 0.606 | **INVERT** | 2.64, 0.339 |
| advdiff1d-dual | EXACT | 1893 | 555 | — | 0.0185, 0.00542 |
| lib-transport | EXACT | 47.1 | 3.32 | — | 0.55, 0.0387 |
| lib-heat-ode-das | — | — | — | — | no certificate |
| lib-beam-eb | EXACT | 0.0115 | 0.00201 | — | 0.0442, 0.00772 |
| lib-wave-damped | EXACT | 14.3 | 1.09 | — | 0.548, 0.0416 |
| gain-transport | gam | 1.39e3 | 996 | — | 0.0349, 0.025 |
| gain-transport-dual | gam | 2.57e3 | 1.99e3 | — | 0.0575, 0.0447 |
| gain-heat-dist | gam | 18.2 | 9.6 | — | 0.422, 0.223 |
| gain-reacdiff | gam | 0.0189 | 3.28 | **INVERT** | 0.139, 24.3 |
| dde-scalar | EXACT | 67.6 | 7.29 | — | 0.553, 0.0596 |
| dde-2state | EXACT | 1.71e4 | 1.34e3 | — | 0.0874, 0.00684 |
| dde-2delay | EXACT | 6.28e3 | 1.29e3 | — | 0.0453, 0.00933 |
| dde-2state-dual | EXACT | 44.4 | 7.91 | — | 0.016, 0.00285 |
| poincare | gam | 104 | 111 | — | 0.714, 0.76 |

**3 inversions in 17 comparable cases**, two of them on identical systems. `op/row` over the 34
points spans **0.00285 to 24.3 — a factor of 8,526**.

**Degeneracy control** (the cuADMM session's, and it is necessary): a candidate pair spanning one
ray cannot invert under *any* degree-1 homogeneous measure, so a non-inversion there is
uninformative rather than negative. Its discriminator is `op/row` agreeing to 3+ digits between
the two points. Only **poincare** comes close (0.714 vs 0.76, 6% apart), so its non-inversion is
the least informative row; every other pair differs by ≥1.29× and its negative is real. That
session's first attempt at the same test returned `op/row` *exactly* constant and 0 of 10
disagreements, purely because all five of its candidates were solved points lying on one ray —
it nearly filed a non-reproduction. The control is not optional.

**Independent replication.** That session then reproduced the inversion on different physics with
an **analytic** infeasibility boundary (Dirichlet heat PIE, Lyapunov operator pinned to the
identity, f = λ/λ\* = 1.20 where infeasibility is Poincaré, i.e. a theorem rather than a solver
verdict): **16 rank disagreements in 66 pairs**, `op/row` spread 20.5×. So the effect is not an
artifact of these plants, of `pielr_opcheck`, or of two interior-point solvers happening to agree.

### The fixed normaliser: safe at X=0, and the coupling row is vacuous either way

`fixnorm.m`, `Nrm_fixed = |Top|·|Aop|` in this gate's own max-coefficient units, residual at the
**trivial point** per equality:

- All 13 feasibility cases: `rel_fixed` between **5.0e-03 and 5.1e-02** — four to five orders
  above the 1e-6 threshold, so **X = 0 is comfortably rejected**.
- All four objective cases: **`1.000e+00, 0.000e+00`** — the negativity row rejects at 1.0, and
  the coupling row `Top'*Qop - Rop` reads **exactly zero**.

So the concern that a fixed denominator would re-admit the trivial point is **not borne out**,
because the gate takes the *worst* ratio across equalities and the negativity row always carries
plant-data constants. It would be a real defect in a gate that took the best or the mean.

The coupling row's exact zero is not new to the fixed normaliser — it reads zero under the
current solution-dependent one too (measured earlier as `[1, 0]`). That row is **vacuous at the
trivial point under either choice** and contributes nothing to rejecting it; it should come out of
the accept test and be reported separately.

Two caveats on this table, both limiting what it can support:

- `|Top| = 1.0000` in every case. This is **not** a normalisation — `pielr_norm_pie` converts and
  repackages but does not rescale (checked). It is that `pielr_maxop` is a max absolute
  coefficient and these PIEs' `T` parameters have unit entries. On this suite `|T|·|A|` therefore
  reduces to `|A|`.
- **This instrument cannot test the cuADMM session's warning against `|T'A+A'T|`.** Its collapse
  is spectral — loss of definiteness at the stability boundary — and a max-coefficient norm cannot
  see it. The two columns here are close or equal (identical on `rd1d-lam0.9` at λ = 0.9λ\*, right
  where a collapse would show), which neither confirms nor refutes the warning. Taking that
  question further requires an induced norm this package does not have.

---

## 9. The baseline

113 cells, **one MATLAB process each** so the R2025b JIT fault costs one row rather than a run,
CSV-append with resume. 18 1-D/DDE cases across tiers 0-3, two further seed sets at tier 2, and
the five 2-D cells. 8.5 hours, **109 rows**; four cells lost, all 2-D.

### The absolute gate is necessary but NOT sufficient

§8 concluded from tiers 1 and 2 that the absolute threshold repaired soundness. **Tiers 0 and 3
refute that.** Bounds still below 1, where a face restriction can only *raise* an optimum:

| case | tier | ratio |
|---|---|---|
| gain-transport | 0 | **0.7466** |
| gain-transport-dual | 0 | **0.7878** |
| gain-reacdiff | 3 | **0.9893** |

The earlier statement was true of the tiers measured and was generalised past them. The threshold
was only an **amplifier**: because the operator residual does not track the constraint, an
absolute threshold on it still admits infeasible points. **No threshold choice repairs a gate
whose measured quantity is wrong.** The numerator is the defect, and it is unfixed.

Where the gate is sound the bounds are tight — tier 2 gives 1.00013 to 1.00186 across all five
objective cases, and Poincaré lands at 1.00013 of the analytic 1/π.

### Complexity: no power law, and the IPM wins everywhere

| fit | exponent | R² | range |
|---|---|---|---|
| `t_bm ~ m` | 1.00 | **0.173** | m 21..1834 |
| `t_bm ~ Ntot` | 0.37 | **0.045** | Ntot 145..46720 |
| `t_ref ~ m` | 0.35 | **0.313** | |
| `t_ref ~ Ntot` | 0.28 | **0.370** | |

Over a 20-fold range in m and 300-fold in `Ntot`, **no power law describes either arm**. At that
range a bad fit means the model is wrong, not that the constant is uncertain, so no exponent is
quoted. Cost on this suite is not size-determined.

`t_bm / t_ref`: **min 1.06, median 9.91, max 341.** The interior-point arm is faster on *every*
case. The low-rank arm has no wall-clock argument in 1-D, and any case for it has to rest on
something else — memory at scale, or 2-D, where it currently cannot run at all.

### 2-D is unmeasured at this budget

Four of five cells failed: three hit the 7200 s/cell timeout, one segfaulted. Only `rd2d-deg2`
produced a row. The 2-D block consumed 6¼ of the 8.5 hours and returned one data point. Any 2-D
claim in this package rests on the historical measurements, not on this baseline.

### The numerator defect, at baseline scale

**9 rank inversions in 56 comparable rows (16%)**, up from 3 in 17. `op/row` over 129 points spans
6.124e-05 to 234.4 — **a factor of 3.83e6**.

Certification: low-rank 56, reference 54, **agreeing on 71 of 73 rows** — the low-rank arm
certifies marginally *more* often than the reference.

### Seed stability is good

Across three seed sets at tier 2, the certify/fail verdict is **identical for all 18 cases** and
bounds agree to ≤7.5e-4 relative. `rel_op` itself varies by up to ~660× on one case
(`lib-beam-eb`, 1.05e-12 to 6.96e-10) without changing any verdict. So the verdicts are not
seed-flaky even though the residual is noisy — worth recording, since the repo's `Test_*` scripts
are seed-flaky and the assumption does not transfer.

---

## 10. The whitener: my objection was wrong, and here is what it is actually doing

§1 of the opening theory pass argued the BM preconditioner (`bm_setup` with `pre=1`) should go, on
two grounds: it forms `full(Ssym*Ssym')` and eigendecomposes it — O(m²) memory, O(m³) time, which
contradicts `theory.tex`'s "never touches an m×m matrix" — and it makes `bm_lm2` minimise `‖W r‖`,
which is not the quantity the gate measures. Paired A/B over the 18 1-D/DDE cases at tier 2, both
arms back to back so a machine slowdown cannot masquerade as an effect. **The prediction was
registered before the run and it failed.**

### 1. The cost objection is refuted at these sizes

`t_setup` is **0.0025 s to 0.019 s** against solve times of 1–500 s: about **0.01 %** of runtime.
The O(m³) term is real but asymptotic; at m ≤ 1834 it is free. It would only bite at the 2-D sizes
(m ~ 6000) where §9 showed the method cannot complete at all, so it is not the binding problem.

### 2. It pays for itself on time

`pre=1` is faster on **12 of 17** complete pairs, by up to **3.1×** (`rd1d-lam0.5`, 1.47 s vs
4.49 s). `pre=0` wins on three, two are ties. Whatever the setup costs, it is repaid several times
over in the solve.

### 3. On residual it is a wash — so the alignment argument is not supported

`pre=1` reaches the lower `rel_op` on **7** cases, `pre=0` on **8**. If the whitened objective were
systematically pulling the search away from the gate's quantity, `pre=0` should win consistently.
It does not. The de-alignment is real as a statement about *units* but does not show up as a worse
final residual.

### 4. Turning it off loses a certificate

`rd1d-lam0.9` certifies with `pre=1` and **fails** with `pre=0`. No case goes the other way.

### 5. What it is really doing: revealing rank deficiency

The measurement the critique missed. `P.rankA` is the numerical rank of the equality system, and
**A is rank-deficient in 14 of 17 cases**:

| case | m | rankA | dropped |
|---|---|---|---|
| gain-reacdiff | 160 | 127 | **33 (21 %)** |
| lib-heat-ode-das | 320 | 279 | 41 (13 %) |
| rd1d-n2 | 236 | 210 | 26 (11 %) |
| gain-heat-dist | 146 | 127 | 19 (13 %) |
| dde-2state-dual | 248 | 237 | 11 |
| gain-transport-dual | 131 | 121 | 10 |
| lib-beam-eb / lib-wave-damped / dde-2delay / gain-transport | 212/212/227/130 | 203/203/218/121 | 9 |
| advdiff1d-dual | 61 | 55 | 6 |
| rd1d-lam0.5 / rd1d-lam0.9 | 59 | 55 | 4 |
| dde-2state | 71 | 70 | 1 |
| lib-transport / dde-scalar / poincare | 53/61/52 | — | **0** |

So these executives build equality systems carrying **up to 21 % redundant rows**, and `pre=1` is
not merely a preconditioner — it is a **rank-revealing reduction** that solves the well-posed
reduced system while `pre=0` works with the dependent rows still in. That is the most plausible
explanation for the speed result, and it connects to a known upstream hazard: redundant equality
systems are one of the two documented triggers of `sossolve`'s "Size b mismatch".

### Verdict

**Keep the whitener; `pre=1` stays the default.** The opening objection was wrong on cost at
these sizes and unsupported on residual. `opts.pre` is retained as a knob because the O(m³) term
will matter if the 2-D arm is ever made to run.

What survives of the critique is narrower and still real: `bm_lm2`'s **exit test** compares
`‖W r‖` against a threshold calibrated in operator units, which is why `tol` fired two orders above
the gate and needed `lmtol = gate/100`. That is a defect in the *exit test*, not an argument
against the preconditioner, and it is the next thing to fix.
