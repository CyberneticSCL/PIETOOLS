# The algorithm: what it does, how it was exercised, and what is open

Written 09/22/2026 (CC) as a handoff for an independent design critique.

**How to read this file.** `README.md` says what the package claims and `GUIDE.md` says how to use
it. This file describes the *mechanism* and the *methods used to exercise it*, and lists open
questions. It deliberately does **not** tell you which claims are true. Numbers appear only as
"this log records X", with the artifact named so you can re-run or re-read it and reach your own
conclusion. Every conclusion the author reached is collected in section 6 as a numbered
hypothesis, each with the condition that would refute it and the known weakness in its current
evidence. Treat those as open, including the ones that look settled.

---

## 1. The problem being solved

For a 2-D PIE with operators `(T, A)`, exponential stability is sought via a positive PI operator
`P` satisfying

```
T' P A + A' P T  =  -D ,     P >= 0 ,  D >= 0
```

`poslpivar_2d` parameterises each positive operator as `Z' Q Z` over a monomial basis `Z`, making
the unknowns **Gram matrices** `Q`. The LPI becomes an SDP: linear equations `A(X) = b` on
block-diagonal `X = blkdiag(X_1..X_B)` with each `X_i >= 0`.

On the benchmark used throughout (`nb_rd2d(1,lam)`, `set2d_deg(4,[])`): `B = 2`, `Ns = [8 744]`,
`m = 5,920` equations, 277,176 Gram unknowns.

**The design premise** (the author's, for you to evaluate): if a certificate exists whose Gram
blocks have rank `r << Ns(i)`, it lies in a small face of the cone, and solving on that face is
much cheaper than solving the SDP.

---

## 2. The pipeline

`pielr_certify(PIE, opts)` — 657 lines, 15 helpers in `private/` (1,835 lines total).

```
norm_pie            normalise the PIE input (T,A; vars/dom defaulted from T)
assemble            poslpivar_2d + lpi_eq_2d, direct form, psatz 0
                      -> SeDuMi data (Atf, bf) via raw_data
local_partition     recover per-block structure (B, Ns, Kf)
DISCOVERY           find a face V = {V_1..V_B} with orthonormal columns
                      route 'bm' (default) or 'mintrace'
CERTIFY             restrict_solve: X_i = V_i S_i V_i', solve the small SDP, lift
refine_blocks       shrink each V_i while the gate still passes
GATE                opcheck_2d: operator residual over all 36 opvar2d cells
                      < opts.gate (default 1e-6) AND every Gram block PSD
```

### Discovery, route `bm` (`disc_bm`, lines 391–501)

A rank ladder: `opts.rank` first (default estimate `2*n2`), then upward to `opts.maxrank`. Per
rung, per start (an optional `opts.w0` factor first, then random seeds `[11 22 33]`):

1. `bm_setup` — package the SeDuMi data, normalise `bf` by `nb0 = norm(bf)`, build a whitener `W`.
   With `pre=1` this forms a dense `m×m` matrix and eigendecomposes it; `pre=0` substitutes a
   sparse identity.
2. `bm_proj` (Frobenius projection onto `{A(X)=b}`) then `bm_dr` (Douglas–Rachford between that
   affine set and the rank-limited PSD set) to produce a starting factor.
3. `bm_lm2` — Levenberg–Marquardt on `F(w) = W (A(YY') - b)` with `X_i = Y_i Y_i'`. Selects the
   cheaper normal equations: primal `J'J` when `nv <= m`, dual `JJ'` otherwise.
4. Gate the BM point itself via `opcheck_2d`; if it passes, accept it directly.
5. Otherwise set `V_i = orth(Y_i)` and pass that face to `restrict_solve`.

### Two distinct acceptance paths

The code can accept a certificate by either of two mechanisms, which solve different problems:

| path | what it solves | availability condition |
|---|---|---|
| `restrict_solve` on the face | the equations **exactly**, restricted to the face | `np >= rM` |
| the BM point `YY'` directly | **minimises** the residual; no face involved | always |

`pielr_certify` attempts the face path first and falls back to the BM point, recording the note
`'accepted the BM point directly'` when it does. Which path fired on any given run is therefore
visible in `cert.notes`, and is worth checking before attributing a result to the face machinery.

---

## 3. Stated design rationale

Summarised from the code's own headers, as rationale on record — not as endorsement.

**The gate is at the operator level, not an equality residual or solver flag.** `pielr_certify`'s
header gives the reason: on this LPI family `norm(b) ~ 5e-6`, so `X = 0` satisfies the equality
rows to ~1e-9 relative while certifying nothing; and a least-squares point is reported there at
equality residual 1.7e-11 with a Gram block having a negative eigenvalue. `opcheck_2d` pushes the
Gram vector back through the real operator and checks all 36 `opvar2d` cells, asserting
`npar == 36`.

**Proposer / checker split.** Discovery is heuristic; certification is exact on the proposed face;
the gate arbitrates. The header states the intended consequence: reported ranks are upper bounds
exhibited by verified certificates, and a rank that fails is reported "not reached" rather than as
a proven floor.

**`b` normalisation.** `bm_setup` divides `bf` by `nb0` unconditionally, so residuals inside the BM
machinery are in normalised units while `opts.face.S` is in original units. See the UNITS section
of `pielr_certify`'s header.

---

## 4. Methods used to exercise it, and what the logs contain

Each entry states what was run and what the artifact records. Interpretation is left open.

**4.1 Rank ladder at a fixed operating point.** `lam/lam* = 0.50` on `nb_rd2d(1,·)`, ranks swept
with `maxrank` pinned so the ladder could not escalate. `r2d.log` records per-attempt rows for
rank [2 2], [3 3] and [4 4], an accepted result at rank [3 4] with `op_rel = 9.457291e-07`, and
the note that the BM point was accepted directly. `vcert.log` records a re-run of that saved
certificate against a freshly rebuilt program, reporting the same value at ratio 1.0000.

**4.2 Budget varied with rank held fixed.** `r2dbud.log`: rank pinned [2 2], seeds 11 and 22,
`lmit` ∈ {100, 200, 400, 1000}, with the 4000 point taken from `r2d.log`. Both `raw_rel` and
`op_rel` are recorded per attempt, so their relationship across budget can be read directly.
Note when using this table: the rungs are **nested** — LM is deterministic from a seed, so a
400-iteration run is a prefix of a 1000-iteration run, not an independent sample. Runs also used
differing `maxNumCompThreads` (2 here, 4 for the 4000 point, 6 for `ar2`), and the sensitivity of
`op_rel` to thread count was never measured.

**4.3 Budget at higher ranks.** `r2dattr.log`: rank pinned [4 4] at `lmit = 400`, both seeds.
`r2dr3.log`: rank pinned [3 3] at `lmit = 12000`, seed 11, which ran to the full budget.

**4.4 Anisotropy retest.** `ar2.log`: the anisotropic variant from `demo_your_system.m`, rank
pinned [2 2], seeds 11/22/33, `lmit` ∈ {400, 2000, 8000}. Two seeds return bit-identical values at
2000 and 8000; one does not. `bm_lm2` has three early exits (`tol`, `stagnant`, `damping`) and the
function returns its reason, but the runs did not print it — **which exit produced those identical
rows was not determined.**

**4.5 Frank–Wolfe as a face proposer.** `scratchpad/fwface/`: a port of the outer loop of
Talitckii & Peet (JMLR 2025) Alg. 4, swept over domain variants, trace bounds, ranks 1–11 and
iteration counts to 2e6, with a full-space control. Caution when using these: they were produced
by subagents and two load-bearing values were never independently reproduced.

**4.6 Known harness defect affecting all of the above.** `pielr_certify` returns `op_rel = NaN`
when no rung certifies. Drivers that read `cert.op_rel` therefore get NaN, and at least one
(`scratchpad/aniso_rank2/ar2.m`) printed a conclusion computed from it. **Read the per-attempt
rows in the logs, not the summary lines.**

---

## 5. Open questions

Posed without preferred answers.

1. **Which acceptance path is doing the work?** Check `cert.notes` on the runs in §4.1 and
   determine whether the face path or the BM-point fallback produced each accepted certificate. If
   the fallback dominates, what does that imply for the design premise in §1?
2. **Is the descent objective aligned with the acceptance test?** `bm_lm2` minimises a Frobenius
   equality residual; the gate is an operator residual. The per-attempt rows in `r2dbud.log` carry
   both. Are they related monotonically? If not, should the objective change?
3. **What is the relationship between `np` and `rM`?** A width-`r` face carries
   `np = sum_i r_i(r_i+1)/2` free parameters; `restrict_solve` solves the equations exactly.
   `rM` was measured on the 1-D program (48) but **not** on the 2-D one. Measure it, and determine
   whether the face path can succeed at the ranks in use.
4. **Does anything guarantee the search?** Burer–Monteiro no-spurious-minima results require
   `r(r+1)/2 > m`, i.e. `r >~ sqrt(2m) ≈ 97` at `m = 5,920`; the package runs `r = 2..4`.
   `theory.tex` §"What can be said in its favor" states this. What, if anything, justifies the
   search at these ranks?
5. **Is discovery affordable at scale?** `bm_setup` with `pre=1` forms a dense `m×m` whitener and
   eigendecomposes it. The `pre=0` path exists. Neither has been run at the largest problem size.
6. **Should the rank ladder be per-block?** The accepted certificate in §4.1 is rank [3 4] over
   blocks of size [8 744], while the ladder sweeps one rank across all blocks.
7. **Do thread count and seed choice move `op_rel` enough to matter?** Never measured; several
   comparisons in §4.2 mix thread counts, and the seed index also selects the initial magnitude
   from `[1e-1 1 1e1 1e2 1e4 1e6]`.
8. **Is `bm_lm2`'s stopping rule sound?** Default was `maxit = 400` with `tol = 1e-16`, which never
   fires; a windowed stagnation test (`rtol = 1e-3` per 100 iterations) was added. Can it stop a
   descending run? Which exit fired in §4.4?

---

## 6. Hypotheses to test

These are the author's claims, restated as hypotheses. **None is established.** Each is scoped to
a specific program and operating point, states what would refute it, and names the known weakness
in the evidence that currently bears on it. Unless stated otherwise the program is
`nb_rd2d(1,lam)` with `set2d_deg(4,[])`, giving `B = 2`, `Ns = [8 744]`, `m = 5,920`; the stability
boundary is `lam* = 2*pi^2 = 19.7392`; and the gate is `opcheck_2d`: operator residual `< 1e-6`
across all 36 `opvar2d` cells **and** every Gram block PSD.

**H1 — a rank-[3 4] certificate exists at half the stability boundary.**
At `lam = 0.50*lam* = 9.8696`, there exists `X = blkdiag(X_1,X_2)` with `rank(X_1) <= 3`,
`rank(X_2) <= 4` satisfying the gate.
*Refuted by:* loading `r2d_cert050.mat`, re-gating it against a freshly assembled program, and
obtaining `op_rel >= 1e-6` or any non-PSD block.
*Evidence:* `r2d.log` records `op_rel = 9.457291e-07`; `vcert.log` records a re-run at ratio
1.0000. *Weakness:* the margin is 5.4% inside the gate, and the sensitivity of `op_rel` to thread
count, rebuild and settings has not been measured. The accepted point came from the BM-point
fallback, not the face path — see H5.

**H2 — at `lam = 0.50*lam*`, per-block rank [2 2] is limited by rank, not budget.**
With rank pinned at [2 2], `op_rel` is bounded below by roughly 1e-5 as `lmit` grows, so no
practically attainable budget reaches the gate.
*Refuted by:* any rank-[2 2] run at that operating point that gates.
*Evidence:* `r2dbud.log`, `lmit` in {100, 200, 400, 1000}, plus the 4000 point from `r2d.log`.
*Weakness, and it is substantial:* the rungs are **nested**, not independent — LM is deterministic
from a seed, so the 400-iteration run is a prefix of the 1000-iteration run. The tabulated series
is also a best-of-two-seeds envelope whose top rung switches seed. An extrapolated exponent
previously quoted from this data is not identified by it: across estimators it ranges -0.10 to
-0.35, implying anywhere from 4.5e06 to 1.2e12 iterations, and dropping a single rung collapses
the fit. Treat the *shape* of this evidence as weak even if the conclusion turns out right.

**H3 — at the same operating point, rank [4 4] is limited by budget, not rank.**
`op_rel` crosses the 1e-6 gate somewhere between `lmit = 400` and `lmit = 4000`.
*Refuted by:* a rank-[4 4] run at `lmit` far above 4000 that fails to gate, or one at
`lmit <= 400` that gates.
*Evidence:* `r2dattr.log` (lmit 400, both seeds) against the accepted run in `r2d.log`.
*Weakness:* two budget points, two seeds, and differing thread counts between the two logs.

**H4 — the descent objective and the acceptance test are not monotonically related.**
For rank [2 2] at `lam = 0.50*lam*`, as `lmit` increases, `raw_rel` — the normalised Frobenius
equality residual that `bm_lm2` actually minimises — decreases monotonically while `op_rel`, the
gate quantity, does not.
*Refuted by:* per-attempt rows showing both quantities moving monotonically together.
*Evidence:* `r2dbud.log` carries both columns for every attempt.
*Weakness:* nested rungs again, and the effect size on `op_rel` is comparable to the unmeasured
thread-count and seed variability.

**H5 — `restrict_solve`'s exact path cannot succeed at the ranks in use on the 2-D program.**
A width-`r` face carries `np` free parameters, where `np` sums `r_i(r_i+1)` halved over blocks,
and `restrict_solve` solves the equations exactly — so the path requires `np >= rM`, where
`rM = rank(M)` for the restricted system. For rank [3 4], `np = 16`. The hypothesis is that
`rM` greatly exceeds 16 on this program, hence the face path is structurally unavailable and any
acceptance must come from the BM-point fallback.
*Refuted by:* measuring `rM` on the 2-D program and finding it comparable to `np`, or exhibiting
a low-rank face on which `restrict_solve` gates.
*Evidence:* `r2d.log` records `face op rel = 1` alongside an accepted BM point. **`rM` has never
been measured on the 2-D program** — the value 48 that appears elsewhere is from the 1-D `heavy`
program, which has `B = 3`. This hypothesis currently rests on one symptom plus an unverified
dimension argument.

**H6 — Frank–Wolfe cannot supply a usable face, and neither can proposers of its type.**
On the 1-D `heavy` program, no face built from the leading eigenvectors of a Frank–Wolfe iterate
passes `gate1d` at any rank 1 through 11 or any iteration count up to 2e6; and the stated reason
— narrow faces are generically infeasible because `np < rM`, wide ones instead fail positivity —
generalises to any proposer that reads a low-rank subspace off a computed point.
*Refuted by:* any such face that gates, or any proposer of that class that succeeds.
*Evidence:* `scratchpad/fwface/`, including a full-space control that does gate.
*Weakness:* produced by subagents; two load-bearing numbers — the `dist(b, range M)` crossover and
the random-face principal-angle null — were never independently reproduced. The generalisation
beyond Frank–Wolfe is an argument, not a measurement.

**H7 — the "+1 rank per block under anisotropy" rule is a threshold effect, not a rank law.**
On the anisotropic variant `x_t = x_s1s1 + 0.5*x_s2s2 + 2x` on the unit square with Dirichlet
conditions, rank [2 2] fails the gate, but by a margin small enough that the failure is not robust.
*Refuted by:* a rank-[2 2] run that gates, making it not a rule at all; or runs showing a stable
margin well outside measurement variability, making it a genuine rank law.
*Evidence:* `ar2.log`, seeds 11/22/33 at `lmit` in {400, 2000, 8000}; the best value reaches
1.000e-06 against a 1e-6 gate.
*Weakness:* the deciding value is printed at four significant figures and sits on the threshold, so
whether it passes cannot be determined from the log. Two seeds return bit-identical values at two
budgets, and **which of `bm_lm2`'s three exits caused that was never determined.**

**H8 — the earlier reported reach of `0.15*lam*` was an artifact of the iteration cap.**
That figure was measured with `bm_lm2` at `maxit = 400` and `tol = 1e-16`, a tolerance that never
fires, so the run stopped at a fixed cap rather than at convergence.
*Refuted by:* reproducing a reach near `0.15*lam*` under a budget large enough that the stopping
rule demonstrably fires on convergence.
*Evidence:* H1 and H3 together, if both hold.
*Weakness:* depends entirely on H1 and H3 — it is an inference from them, not a separate
measurement.

An adversarial audit of several of these was started and stopped part-way. Its partial output is in
`scratchpad/audit/` and is itself unreviewed: treat it as raw material, not as findings.

## 7. Evidence map

| what | where |
|---|---|
| rank × budget at 0.50, saved certificate | `scratchpad/reach2d/r2d.log`, `r2d_cert050.mat` |
| budget ladder, rank pinned [2 2] | `scratchpad/reach2d/r2dbud.log` |
| rank-4 at lmit 400 | `scratchpad/reach2d/r2dattr.log` |
| rank-3 at lmit 12000 | `scratchpad/reach2d/r2dr3.log` |
| re-run of the saved certificate on a fresh build | `scratchpad/reach2d/vcert.log` |
| anisotropy ladder | `scratchpad/aniso_rank2/ar2.log` |
| Frank–Wolfe study | `scratchpad/fwface/` |
| partial, unreviewed audit output | `scratchpad/audit/` |
| 1-D regression suite | `tests_1d/` |

**Path resolution.** `pielr_path` strips 89 entries from a nested `.claude/worktrees/` copy of the
repo and asserts single resolution of `poslpivar_2d`, `lpi_eq_2d`, `lpiprogram` and `monomials`;
it errors if any resolves more than once. `tests_1d/t1d_path.m` does the equivalent for the 1-D
suite, and `tests_1d` T0 asserts it. Worth confirming for yourself before trusting any run.

---

## 8. Starting points

`pielr_certify.m`'s header is the contract. `opcheck_2d.m` is the acceptance test and everything
depends on it being right. Questions 1 and 3 in §5 are the ones whose answers would most change
how the rest of the package should be read.

Nothing in the working tree is committed; `git status` shows the annotated changes. A caution
carried from the measurements in §4: a reported result within a small factor of the gate, taken
under a fixed iteration cap, may be reporting the cap rather than the method — check the budget a
number was measured at before relying on it.
