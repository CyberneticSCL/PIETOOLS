# `tests_1d` — the 1-D regression suite for the low-rank certifier

    >> run_regression_1d

One call, no arguments. It prints a numbered `PASS`/`FAIL` line per test with the
measured numbers underneath, and ends with a sentinel a human or a script can read at
a glance:

    REGRESSION_1D: 10/10 PASS   [ 245.5 s ]
    REGRESSION_1D:  8/10 PASS  FAILED: T5 T6   [ 248.9 s ]

`run_regression_1d` also returns the per-test struct array if you ask for it.

**Measured runtime: 4.1 min** (245.5 s) on this box — Windows 11, R2025b,
`maxNumCompThreads(4)`, one MATLAB. Three consecutive full runs measured 355.9, 276.0
and 245.5 s; the spread is MATLAB's own JIT and file caching warming up, not the tests.
Fast is a feature: a suite nobody runs catches nothing. The most expensive test is T6
(59 s, five LPI builds and ten searches), then T7 (27 s) and T8 (22 s); the shared
fixtures cost 16 s before any test runs, and T0–T5 together cost 25 s.

---

## Why a 2-D package has a 1-D test suite

This package certifies 2-D stability, and everything in it above `private/` is
2-D-specific. **The core under test is not.** The Burer-Monteiro engine —

    bm_setup   bm_dr   bm_lm2   bm_proj   bm_resid   bm_report   raw_data

(the list lives in `t1d_corefiles.m`) — never mentions a spatial dimension. It sees
only the SeDuMi triple `(Atf, bf, Ns, Kf)` and a rank profile `rv`. Neither does the
face machinery in `restrict_solve`, apart from which operator gate it calls at the end.
**Every defect this campaign found lived in those files or in `restrict_solve`.**

So a 1-D driver over the same core tests the shipped code rather than a
re-implementation — and it does it at `m = 56` with a 0.06 s interior-point solve
instead of the 2-D sizes, where one confirmation run costs minutes to hours. That is
the whole argument for putting the suite here, and it has exactly one load-bearing
assumption: that the seven files in this directory really are the shipped ones.
**T0 asserts that byte for byte on every run.** If they ever diverge the suite says so
first and loudly, because a divergence means every other number it prints is a
statement about a copy.

The system is 1-D reaction-diffusion `phi_t = phi_ss + lam*phi` on `[0,1]` with
Dirichlet ends, stable iff `lam < lam* = pi^2 = 9.8696`. Having an analytic boundary is
what makes T8's attribution check non-vacuous: above `lam*` no certificate exists, so a
miss there *must* be charged to the relaxation.

## Provenance

The 1-D driver was not written for this suite. It was built for the lambda-reach study
in `scratchpad/reach1d/` and is copied here rather than rewritten:

| file | origin |
|---|---|
| `bm_setup` `bm_dr` `bm_lm2` `bm_proj` `bm_resid` `bm_report` `raw_data` | `../private/`, **byte-identical**, asserted by T0(a) |
| `restrict_solve_1d.m` | `../private/restrict_solve.m`, differing in the declaration and **one** line (the gate call), asserted by T0(f) |
| `gate1d` `opcheck` `op2gal` `polyeval2` `build_stab_st` `rd_pie_lam` | `scratchpad/reach1d/`, verbatim plus a provenance header |
| `t1d_set` | `scratchpad/ladder/lset.m`, verbatim, renamed so the two harnesses can share a path |
| `t1d_build` `t1d_ipm` | `scratchpad/reach1d/r1build.m`, `r1ipm.m`, renamed |
| `t1d_bm` | `scratchpad/reach1d/r1bm.m` + the `.w0` hook, `R.log` and `R.V` (see its header) |
| the rest (`run_regression_1d`, `t1d_path`, `t1d_subfun`, `t1d_codelines`, `t1d_corefiles`, `t1d_ownnames`, `t1d_whichcore`) | new |

`t1d_bm` carries a verbatim copy of `pielr_certify`'s `fitw` subfunction, because
MATLAB cannot call a subfunction from outside its file. **T9(b) asserts that copy still
matches the shipped one code-for-code**, so it cannot drift into testing its own fill
rule. T4 and T9 go further and *extract and execute* the shipped `padw` and `fitw` out
of `pielr_certify.m` (`t1d_subfun`), so the fill behaviour under test is the shipped
behaviour, not a transcription.

## What each test asserts

Each test names the defect it guards and the measured number behind it. The letters
`B1`…`B7` are the defect list this suite was commissioned from.

| # | guards | asserts | threshold |
|---|---|---|---|
| **T0** | **B6** path shadowing, and the suite's own premise | the seven core files are byte-identical to `../private`; they resolve *inside* `tests_1d`; every name this directory owns resolves uniquely; the PIETOOLS entry points resolve uniquely and inside the tree; `t1d_ownnames` covers every `.m` here; `restrict_solve_1d` differs from the shipped file in exactly 2 code lines (declaration + gate call); `R1_GATE` unset | 7/7 identical, 7/7 resolved, 18/18 and 6/6 unique, exactly 2 differing lines |
| **T1** | **B4** gate integrity | a known certificate is accepted (`relP` 1.89e-09); the **zero Gram** is rejected (2.66); a **1.02× scaled** certificate is rejected (5.96e-03); an **indefinite twin** — the certificate perturbed along a null direction of the symmetric constraint map — has an *identical* residual and is rejected by the **PSD half only** | accept `relP < 1e-7`; reject zero `> 1e-2`; reject scaled `> 1e-4`; twin `relP < 1e-6` and within 1e-3 **relative** of the certificate, `psd` false, `min(mineig/‖Q‖_F) < -1e-2` (measured −0.220) |
| **T2** | reproducibility | two identical runs agree bit for bit: `q`, `r`, `relP`, attempt/LM/DR/gate counts, and the start log | `isequal` on `q`; exact equality on all counters |
| **T3** | **B5** b-normalisation | `‖bf‖ = 1` after `bm_setup`; `b → 10b` leaves `bf`, `W` and the reported `raw_rel` unchanged and scales `nb0` by exactly 10; `b ∈ range(A)` | `‖bf‖−1 < 1e-14`; `‖W1−W2‖_F = 0` exactly; `raw_rel` relative difference `< 1e-12` (measured 1.5e-16); `bout < 1e-12` (measured 4.5e-15) |
| **T4** | **B2** zero-pad fixed point | on a settled rank-1 iterate the **shipped** `padw` gives padded columns an **identically zero** Jacobian block, live columns do not, and 150 LM steps leave the padded entries *exactly* zero; the **shipped** `fitw` fill is mobile | `max‖J(:,padded)‖ == 0` **exactly** (pinned defect); `min‖J(:,live)‖ > 1e-2` (measured 0.432); `max|w(padded)| == 0` over all 42 padded entries; `min‖J(:,surplus)‖ > 1e-6` (measured 4.09e-04) |
| **T5** | **B3** `restrict_solve` false negatives | the base face certifies and the solve on it is **determined** (`np == rM`, `clip = 0`); the wide face **provably contains** it; a well-conditioned widening still certifies; a random face does not; and the ill-conditioned widening's current failure is **pinned** | `np == rM` (3 == 3); containment `< 1e-12` and the base Gram re-gates at 4.04e-07; ≥2 of 3 single-block widenings certify; any lost widening carries `clip > 1e10`; random face `relP > 1e-2` (measured 2.56) |
| **T6** | **B1** `bm_lm2` truncation | over `lam/lam* = 0.41…0.45` the verdict is **contiguous** and **identical at twice the LM budget** | exactly 1 contiguous run (at the shipped budget of 400 this grid gives `[1 0 0 1 0]`, **2 runs** — the test bites); verdicts at 4000 and 8000 `isequal`; all 5 certify |
| **T7** | rank scaling in `n` | the search certifies at small rank for `n = 1` and `n = 2` decoupled replicas, and respects the replica bound | `max r(1) ≤ 2` (measured 1); `max r(2) ≤ 4` (measured 2); `r(2) ≤ 2·r(1)` |
| **T8** | the attribution rule | every miss on T6's grid is charged to the **search** or the **relaxation** by re-solving the same program with the IPM; and above `lam*` the charge lands on the relaxation | attributed count == miss count; at `1.05 lam*` both arms miss and `IPM relP > 1` (measured 10.36); the reference arm certifies inside the region, so the verdict discriminates |
| **T9** | **B7** the `.w0` hook | the shipped dispatch puts the `w0` branch, and its `fitw` call, strictly before the seeds branch; `t1d_bm`'s local `fitw` still matches the shipped one; `fitw`'s surplus columns are nonzero; and at run time `w0` leads at **every** rung | guard/`fitw`/seeds at lines 414/415/416 in order; `fitw` 15 code lines, equal; `min` surplus column norm `> 0` (measured 9.7e-04); first start is `w0` at every rung, ≥2 rungs observed |

## Findings from the suite's own first runs

**1. The superset invariant fails in 1-D — pinned by T5(d).** *"A face that contains a
certifying face must still certify"* does **not** hold for the shipped restricted solve
here. Take the rank-`[1 1 1]` face the search certifies at `lam/lam* = 0.10`
(`relP` 4.04e-07) and widen every block by **one arbitrary orthonormal direction**:

* containment is exact — `max‖V − Vw Vwᵀ V‖ = 7.1e-16`, and the base Gram re-gated on
  the wide face still certifies at 4.04e-07, so a certificate **provably exists** there;
* the restricted map becomes nearly singular — `cond(M) = 7.8e+07`, smallest QR pivot
  `1.9e-08` of the largest — while `np == rM` (9 == 9) still reports a *determined*
  system;
* the solve returns a point with `clip = 6.5e+11`, and its PSD projection certifies
  nothing: `relP = 2.658`, which is exactly the **zero-Gram** value.

No alternative solve rescues it. Backslash, `lsqminnorm`, truncated SVD at 1e-6 and
least squares on **all** rows were each measured on this face; the best any managed was
`clip 2.6e+03`, `relP 2.636` — still rejected. Tightening the `rM` tolerance makes it
worse (at 1e-6, `rM` drops to 8 and `clip` rises to 3.05e+14).

The mechanism: the restricted solve demands *exact* satisfaction of a selected square
subsystem, but acceptance only ever required the **gate**. The base certificate
satisfies the full equality system to `rel_eq_full = 2.1e-07`, not exactly, so the exact
solution of the widened subsystem is a *different* point — and because the added
direction is nearly unconstrained, it is an enormous one. In 2-D with `m` in the
thousands the padded directions stayed reasonably constrained, which is why the repair
verified there did not expose this.

It is **latent, not a live wrong answer**: the shipped search never builds such a face.
`t1d_bm`/`pielr_certify` set `V{i} = orth(Y_i)`, so every face column is a direction the
iterate actually uses. (Note the interaction with B2: because a zero-padded column has
`orth` drop it, the `prev` rung hands `restrict_solve` an `r`-column face while
reporting rank `r+1`.) T5(d) therefore asserts the **current** behaviour, and prints
`*** the B3 latent defect is fixed — update T5(d) ***` if it ever starts certifying.
The failure is localised: widening block 1 (the 11×11 block) alone triggers it
(`clip 1.33e+14`); widening block 2 or 3 alone does not.

**2. Three live path-shadowing mechanisms (B6), all of them found by writing T0.**

* `pietools_path_update.m` is a single `addpath(genpath(root))`, and
  `.claude/worktrees/<name>/PIETOOLS` is *inside* the tree — so **90 path entries of a
  second checkout** were on the path, including a `private/bm_lm2.m` **measured to
  differ** from the current file. It happened to be ordered behind the real one;
  `genpath` order is not an invariant worth relying on. `t1d_path` removes them.
* Out-of-tree copies: the scratchpad harnesses this suite descends from carry
  same-named ancestors of nearly every file here. `t1d_path` removes any path entry
  outside the tree that provides a name in `t1d_ownnames`.
* **The current folder**, which is searched *before* the path and cannot be removed
  from it. This bit a probe run from the scratchpad: it silently called the
  scratchpad's `bm_setup`. `run_regression_1d` therefore `cd`s here for the run, and
  T0(b) verifies the resolution afterwards regardless.

**3. `MATLAB`'s `regexp` does not implement `\b`.** `'^\s*function\b'` matches nothing;
MATLAB uses `\<` and `\>`. This silently made `t1d_subfun`'s slice run to end of file
and sweep up three neighbouring subfunctions. Fixed, and worth knowing before writing
any other source-level assertion here.

## Deliberate non-assertions

**The Galerkin check is reported, never asserted.** `op2gal` discretises an operator's
quadratic form on Gauss-Legendre nodes, which is a genuinely Gram-free second opinion —
but the `t<s` / `t>s` indicator costs about **8% relative quadrature error** even on the
*analytic* `P = I` certificate (`min eig(-Dop) = −7.3e-03` against `max = +9.5e-02` at
`N = 60`). Worse, on this fixture the **zero Gram scores better** on `eDmin`
(−2.7e-08) than a real certificate does (−3.4e-07), so it cannot discriminate in the
direction that matters. `gate1d`'s `do_gal` branch consequently rejects perfectly good
interior-point certificates, and T1 prints `ePmin`/`eDmin` as a diagnostic only. The
independence T1 needs comes from the null-direction control instead, which is exact.

**The 2-D law `r* = 2n` does not transfer.** T7 does not assert it: measured at
`heavy`, `lam/lam* = 0.10`, `lmit 4000`, the 1-D search first certifies at `r = [1 1 1]`
for `n = 1` and `r = [2 2 2]` for `n = 2`. What T7 guards is the replica upper bound
`r(n) ≤ n·r(1)`, which is what the block-diagonal construction on decoupled replicas
implies, plus that minimum rank stays `O(1)` while `Ns` doubles from `[11 20 11]` to
`[22 40 22]`.

**T6 runs at `lmit = 4000`, not the shipped 400.** At 400 this grid is ragged by
measurement, so a suite that ran there would be permanently red for a reason already
known and documented (`pielr_certify`'s `.lmit` header). 4000 is the budget at which
the reach is stable, and T6 asserts it *stays* stable — contiguous, and unchanged at
8000. Re-running the grid at 400 costs 440 s against 37 s at 4000, which is the other
reason it is not in the suite; the measured 400 verdict vector `[1 0 0 1 0]` is recorded
in `t_T6`'s header instead.

## Caveats

* **It is a 1-D suite.** It cannot see anything that only appears in `opvar2d`, in the
  36-cell `opcheck_2d` scan, or at 2-D problem sizes. It is a fast guard on the shared
  core, not a substitute for a 2-D confirmation run.
* **`restrict_solve_1d` is a copy**, not the shipped file — necessarily, since the
  operators here are `opvar`. T0(f) pins the divergence at the declaration plus one
  line; a shipped change that this copy does not get will fail T0 rather than pass
  silently.
* **T2 establishes determinism within one session only.** It says nothing about
  reproducibility across MATLAB versions or machines, and several assertions
  (`max r ≤ 2`, the `clip` bands) are behaviour pins that a different BLAS could move.
* **T5(d) and T4(a) pin defects, not desired behaviour.** Both print what they expect
  and why. A `FAIL` on either may mean the underlying bug was *fixed* — read the line
  before assuming a regression.
* `t1d_subfun` writes extracted subfunctions into `tempname` directories and leaves
  them on the path for the session. Harmless — `fitw`/`padw` are not names the suite
  owns — but a second `run_regression_1d` in the same session accumulates a few stale
  path entries.
* The `tbxmanager` `startup.m` error at MATLAB launch on this box is unrelated noise.
  Never match a wait loop on the word `Error`; match the `REGRESSION_1D:` sentinel.
